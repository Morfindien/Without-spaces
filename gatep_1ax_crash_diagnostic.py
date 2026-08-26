#!/usr/bin/env python3
"""
CASE-Q013 — 1ax exact-zero crash diagnostic.

Purpose:
- hold the recovered Gate-P environment fixed;
- run N_mscf=1 at exact zero and tiny positive amplitudes;
- capture native SIGSEGV diagnostics with GDB;
- distinguish an exact-boundary defect from a broader 1ax failure;
- never interpret crashes/NaNs/infinities as physical likelihood evidence.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import re
import shlex
import subprocess
import sys
from pathlib import Path
from typing import Dict, Any

AMPLITUDES = [0.0, 1e-14, 1e-12, 1e-10, 1e-8, 1e-6]
FROZEN = {
    "maxiclass_commit": "b5a3af9818d04d2fe696e54972545e8c6805f4ea",
    "act_commit": "0e0cd2c703c62a0e980470b572602233b27750e1",
    "spt_commit": "2cec6e762a8c540484dd5acafc529f0035856350",
    "candl_like": "2.0.3",
    "cobaya": "3.6.2",
    "reference_2ax_zero": 317.4201095748994,
    "reference_R1": 317.4238643933,
    "reproduction_tolerance": 0.1,
}

BASE_YAML = r"""
theory:
  classy:
    path: {maxiclass_path}
    ignore_obsolete: true
    extra_args:
      non linear: hmcode
      N_ur: 2.0328
      N_ncdm: 1
      T_ncdm: 0.71611
      lensing: 'yes'
      P_k_max_h/Mpc: 1.0
      do_shooting: 'yes'
      do_shooting_mscf: 'yes'
      attractor_ic_scf: 'no'
      loop_over_background_for_closure_relation: 'no'
      background_Nloga: 100000
      tol_shooting_deltaF: 0.01
      tol_shooting_deltax: 0.01
      N_mscf: 1
      n_axion_mscf: '3'
      theta_ini_mscf: '2.8'
      theta_prime_ini_mscf: '0.'
likelihood:
  act_dr6_cmbonly.ACTDR6CMBonly:
    input_file: dr6_data_cmbonly.fits
    lmax_theory: 9000
    ell_cuts:
      TT: [600, 8500]
      TE: [600, 8500]
      EE: [600, 8500]
    stop_at_error: true
    params:
      A_act:
        prior: {{min: 0.5, max: 1.5}}
        ref: 1.0011501
        proposal: 0.003
      P_act:
        prior: {{min: 0.9, max: 1.1}}
        ref: 1.0016445
        proposal: 0.01
  candl_like:
    external: !!python/name:candl.interface.CandlCobayaLikelihood ''
    data_set_file: {spt_index}
    variant: lite
    clear_internal_priors: true
    additional_args: {{}}
    feedback: true
    wrapper: null
prior:
  cal_dip_prior: 'lambda A_act: stats.norm.logpdf(A_act, loc=1.0, scale=0.003)'
  gaussian_Tcal: 'lambda Tcal: stats.norm.logpdf(Tcal, loc=1.0, scale=0.0036)'
params:
  logA:
    prior: {{min: 1.61, max: 3.91}}
    ref: 3.0516831
    proposal: 0.002
    drop: true
  A_s:
    value: 'lambda logA: 1e-10*np.exp(logA)'
  n_s:
    prior: {{min: 0.8, max: 1.2}}
    ref: 0.97006392
    proposal: 0.003
  H0:
    prior: {{min: 20.0, max: 100.0}}
    ref: 67.003038
    proposal: 1.0
  omega_b:
    prior: {{min: 0.005, max: 0.1}}
    ref: 0.022462565
    proposal: 0.0001
  omega_cdm:
    prior: {{min: 0.001, max: 0.99}}
    ref: 0.12107093
    proposal: 0.001
  m_ncdm:
    value: 0.06
  tau_reio:
    prior:
      dist: norm
      loc: 0.051
      scale: 0.006
    ref: 0.056456949
    proposal: 0.003
  Tcal:
    prior: {{min: 0.8, max: 1.2}}
    ref: 0.99859358
    proposal: 0.003
  Ecal:
    prior: {{min: 0.8, max: 1.2}}
    ref: 1.0066561
    proposal: 0.003
  log10_maxion_ac:
    prior: {{min: -4.5000137334, max: -3.0004340775}}
    ref: -3.75022390545
    proposal: 0.08
  fraction_maxion_ac:
    prior: {{min: 0.0, max: 0.3}}
    ref: {amp}
    proposal: 0.015
packages_path: {packages_path}
timing: true
debug: true
output: {output}
sampler:
  evaluate: {{}}
"""

def run(cmd, cwd=None, env=None, timeout=900):
    p = subprocess.run(
        cmd, cwd=cwd, env=env, text=True,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        timeout=timeout, check=False
    )
    return p.returncode, p.stdout

def slug_amp(x: float) -> str:
    if x == 0:
        return "0"
    return f"{x:.0e}".replace("+", "").replace("-", "m")

def parse_log(text: str) -> Dict[str, Any]:
    low = text.lower()
    segv = (
        "program received signal sigsegv" in low
        or "segmentation fault" in low
        or "sigsegv" in low
    )
    abort = "sigabrt" in low or "program received signal sigabrt" in low
    nan = bool(re.search(r"(?<![A-Za-z])nan(?![A-Za-z])", low))
    inf = bool(re.search(r"(?<![A-Za-z])(?:\+|-)?inf(?:inity)?(?![A-Za-z])", low))
    classy_loaded = re.findall(r"`classy` module loaded successfully from ([^\n\r]+)", text)
    frame_lines = [
        ln for ln in text.splitlines()
        if re.match(r"^#\d+\s+", ln.strip())
    ]
    source_like = [
        ln for ln in frame_lines
        if re.search(r"\.(?:c|h|cpp|cc|py):\d+", ln)
    ]

    loglikes = []
    for m in re.finditer(r"log-likelihood\s*=\s*([-+0-9.eE]+)", text):
        try:
            loglikes.append(float(m.group(1)))
        except ValueError:
            pass

    act = None
    spt = None
    m = re.search(r"act_dr6_cmbonly[^\n]*Computed log-likelihood\s*=\s*([-+0-9.eE]+)", text, re.I)
    if m:
        act = float(m.group(1))
    m = re.search(r"candl_like[^\n]*Computed log-likelihood\s*=\s*([-+0-9.eE]+)", text, re.I)
    if m:
        spt = float(m.group(1))

    finite = (
        not segv and not abort and not nan and not inf
        and act is not None and spt is not None
        and math.isfinite(act) and math.isfinite(spt)
    )

    return {
        "sigsegv": segv,
        "sigabrt": abort,
        "nan_seen": nan,
        "inf_seen": inf,
        "classy_loaded_from": classy_loaded[-1].strip() if classy_loaded else None,
        "native_frames": frame_lines[-120:],
        "source_frames": source_like[-60:],
        "act_loglike": act,
        "spt_loglike": spt,
        "usable_finite": finite,
    }

def build_yaml(workspace: Path, outdir: Path, amp: float) -> Path:
    maxiclass_path = workspace / "external/mAxiCLASS"
    spt_index = workspace / "external/spt_candl_data/spt_candl_data/SPT3G_D1_TnE_v0/SPT3G_D1_TnE_index.yaml"
    packages_path = workspace / ".cobaya"
    yaml_path = outdir / f"1ax_f_{slug_amp(amp)}.yaml"
    text = BASE_YAML.format(
        maxiclass_path=maxiclass_path,
        spt_index=spt_index,
        packages_path=packages_path,
        output=outdir / f"1ax_f_{slug_amp(amp)}",
        amp=repr(amp),
    )
    yaml_path.write_text(text)
    return yaml_path

def gdb_cmd(yaml_path: Path) -> list[str]:
    return [
        "gdb", "--batch", "--return-child-result",
        "-ex", "set pagination off",
        "-ex", "set confirm off",
        "-ex", "set print pretty on",
        "-ex", "set print frame-arguments all",
        "-ex", "handle SIGPIPE nostop noprint pass",
        "-ex", "run",
        "-ex", "echo \\n===== THREAD APPLY ALL BT FULL =====\\n",
        "-ex", "thread apply all bt full",
        "-ex", "echo \\n===== CURRENT FRAME =====\\n",
        "-ex", "frame",
        "-ex", "info args",
        "-ex", "info locals",
        "-ex", "echo \\n===== REGISTERS =====\\n",
        "-ex", "info registers",
        "-ex", "echo \\n===== DISASSEMBLY AROUND PC =====\\n",
        "-ex", "x/24i $pc-48",
        "--args", sys.executable, "-m", "cobaya.run", str(yaml_path)
    ]

def classify(results):
    by_amp = {r["amplitude"]: r for r in results}
    z = by_amp.get(0.0)
    positives = [r for r in results if r["amplitude"] > 0]

    if not z:
        return "INCOMPLETE"

    if z["diagnostic"]["sigsegv"] and positives and all(
        r["diagnostic"]["usable_finite"] for r in positives
    ):
        return "EXACT_ZERO_BOUNDARY_DEFECT_STRONGLY_CONFIRMED"

    if z["diagnostic"]["sigsegv"] and any(
        r["diagnostic"]["usable_finite"] for r in positives
    ):
        return "EXACT_ZERO_DEFECT_SUPPORTED_WITH_TRANSITION_REGION"

    if z["diagnostic"]["sigsegv"] and all(
        r["diagnostic"]["sigsegv"] for r in positives
    ):
        return "BROADER_1AX_NATIVE_CRASH"

    if z["diagnostic"]["usable_finite"]:
        return "ZERO_CRASH_NOT_REPRODUCED"

    return "UNRESOLVED_1AX_FAILURE"

def write_report(root: Path, all_results):
    classification = classify(all_results)
    zero = next((r for r in all_results if r["amplitude"] == 0.0), None)
    source_frames = zero["diagnostic"]["source_frames"] if zero else []

    lines = [
        "# CASE-Q013 — 1ax Exact-Zero Crash Diagnostic",
        "",
        f"**Technical classification: {classification}**",
        "",
        "Physical Gate-P verdict: **NOT EVALUATED**",
        "",
        "## Frozen state",
        "",
        f"- mAxiCLASS commit: `{FROZEN['maxiclass_commit']}`",
        f"- ACT commit: `{FROZEN['act_commit']}`",
        f"- SPT commit: `{FROZEN['spt_commit']}`",
        f"- candl-like: `{FROZEN['candl_like']}`",
        f"- Cobaya: `{FROZEN['cobaya']}`",
        "",
        "## Boundary matrix",
        "",
        "| fraction_maxion_ac | exit | SIGSEGV | NaN | inf | finite ACT+SPT |",
        "|---:|---:|---:|---:|---:|---:|",
    ]
    for r in sorted(all_results, key=lambda x: x["amplitude"]):
        d = r["diagnostic"]
        lines.append(
            f"| {r['amplitude']:.1e} | {r['exit_code']} | {d['sigsegv']} | "
            f"{d['nan_seen']} | {d['inf_seen']} | {d['usable_finite']} |"
        )

    lines += [
        "",
        "## Exact-zero source frames",
        "",
        "```text",
    ]
    lines += source_frames[:80] if source_frames else ["NO_SOURCE_LEVEL_FRAME_CAPTURED"]
    lines += [
        "```",
        "",
        "## Decision rule",
        "",
        "- Exact zero SIGSEGV + all tiny positive amplitudes finite => exact-zero implementation defect strongly confirmed.",
        "- Exact zero SIGSEGV + some tiny positive amplitudes finite => boundary/transition defect supported.",
        "- Exact zero and all positive amplitudes SIGSEGV => broader 1ax implementation defect.",
        "- Exact zero finite => previous crash not reproduced in this recovered runtime.",
        "",
        "Any crash, NaN, inf, or missing likelihood component is **INVALID FOR PHYSICAL/NESTING INTERPRETATION**.",
    ]
    (root / "CRASH_DIAGNOSTIC_REPORT.md").write_text("\n".join(lines) + "\n")

    payload = {
        "case": "CASE-Q013",
        "classification": classification,
        "frozen": FROZEN,
        "results": all_results,
    }
    (root / "CRASH_DIAGNOSTIC_RESULT.json").write_text(
        json.dumps(payload, indent=2, allow_nan=False)
    )

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workspace", required=True)
    ap.add_argument("--mode", choices=["surgical", "death"], default="surgical")
    ap.add_argument("--tag", default="frozen-optimized")
    ap.add_argument("--summarize-only", action="store_true")
    args = ap.parse_args()

    ws = Path(args.workspace).resolve()
    root = ws / "gatep_1ax_crash_diagnostic"
    root.mkdir(parents=True, exist_ok=True)

    result_files = list(root.glob("runs/*/result.json"))
    if args.summarize_only:
        results = [json.loads(p.read_text()) for p in result_files]
        if results:
            write_report(root, results)
            return 0
        print("No result files found.")
        return 1

    build_dir = os.environ.get("MAXICLASS_BUILD_DIR")
    if not build_dir:
        raise SystemExit("MAXICLASS_BUILD_DIR is not set")

    env = os.environ.copy()
    required_prefix = str(Path(build_dir).resolve())
    env["PYTHONPATH"] = ":".join([
        required_prefix,
        str(ws / "external/DR6-ACT-lite"),
        str(ws / "external/spt_candl_data"),
        env.get("PYTHONPATH", ""),
    ]).rstrip(":")

    # Verify the exact Python import context used by GDB/Cobaya.
    rc, txt = run([
        sys.executable, "-c",
        "import os,classy; "
        "print(classy.__file__); "
        "assert os.path.realpath(classy.__file__).startswith(os.path.realpath(os.environ['MAXICLASS_BUILD_DIR'])+os.sep)"
    ], env=env)
    if rc != 0:
        print(txt)
        raise SystemExit("Hard failure: frozen local classy does not win import routing")

    results = []
    for amp in AMPLITUDES:
        run_dir = root / "runs" / f"{args.tag}_f_{slug_amp(amp)}"
        run_dir.mkdir(parents=True, exist_ok=True)
        yaml_path = build_yaml(ws, run_dir, amp)

        cmd = gdb_cmd(yaml_path)
        (run_dir / "command.txt").write_text(shlex.join(cmd) + "\n")

        try:
            rc, out = run(cmd, cwd=ws, env=env, timeout=1200)
        except subprocess.TimeoutExpired as e:
            out = (e.stdout or "") + "\nTIMEOUT\n"
            rc = 124

        (run_dir / "gdb.log").write_text(out)
        diag = parse_log(out)

        result = {
            "tag": args.tag,
            "mode": args.mode,
            "amplitude": amp,
            "exit_code": rc,
            "yaml": str(yaml_path.relative_to(ws)),
            "diagnostic": diag,
        }
        (run_dir / "result.json").write_text(
            json.dumps(result, indent=2, allow_nan=False)
        )
        results.append(result)

    # Include any earlier pass in the aggregate report (e.g. optimized + debug-symbols).
    all_results = [json.loads(p.read_text()) for p in root.glob("runs/*/result.json")]
    write_report(root, all_results)

    classification = classify(all_results)
    print((root / "CRASH_DIAGNOSTIC_REPORT.md").read_text())

    # Diagnostic success means the run produced an interpretable technical
    # classification. It does NOT mean the cosmology itself passed.
    return 0 if classification != "INCOMPLETE" else 2

if __name__ == "__main__":
    raise SystemExit(main())
