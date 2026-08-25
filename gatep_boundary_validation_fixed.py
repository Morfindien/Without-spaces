#!/usr/bin/env python3
from __future__ import annotations

import csv
import hashlib
import json
import math
import os
import platform
import re
import shutil
import signal
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import yaml

ROOT = Path.cwd()
OUT = ROOT / "gatep_boundary_validation"
RUNS = OUT / "runs"
PROV = OUT / "provenance"
OUT.mkdir(exist_ok=True)
RUNS.mkdir(exist_ok=True)
PROV.mkdir(exist_ok=True)

MAXI = ROOT / "external" / "mAxiCLASS"
SPT = ROOT / "external" / "spt_candl_data"
PACKAGES = ROOT / ".cobaya"

FROZEN = {
    "case_id": "CASE-Q013",
    "hypothesis_id": "HYP-004C",
    "internal_result_id": "BV-Q013-GATEP-BOUNDARY-VALIDATION-001",
    "mAxiCLASS_commit": "b5a3af9818d04d2fe696e54972545e8c6805f4ea",
    "ACT_DR6_commit": "0e0cd2c703c62a0e980470b572602233b27750e1",
    "SPT_D1_commit": "2cec6e762a8c540484dd5acafc529f0035856350",
    "candl_like": "2.0.3",
    "cobaya": "3.6.2",
    "R1_effective_chi2": 317.4238643933,
    "globality_tolerance": 0.1,
}

BASE = {
    "logA": 3.0516831,
    "n_s": 0.97006392,
    "H0": 67.003038,
    "omega_b": 0.022462565,
    "omega_cdm": 0.12107093,
    "tau_reio": 0.056456949,
    "Tcal": 0.99859358,
    "Ecal": 1.0066561,
    "A_act": 1.0011501,
    "P_act": 1.0016445,
}

AC1 = [-3.7000772228, -3.37525565015, -3.0504340775]
AC2 = [-4.4500137334, -4.1250454781, -3.8000772228]
AC_ONE = [-4.4500137334, -3.75022390545, -3.0504340775]


def finite(x):
    return isinstance(x, (int, float)) and math.isfinite(x)


def theory_base():
    return {
        "path": str(MAXI),
        "ignore_obsolete": True,
        "extra_args": {
            "non linear": "hmcode",
            "N_ur": 2.0328,
            "N_ncdm": 1,
            "T_ncdm": 0.71611,
            "lensing": "yes",
            "P_k_max_h/Mpc": 1.0,
        },
    }


def mscf_args(n):
    d = {
        "do_shooting": "yes",
        "do_shooting_mscf": "yes",
        "attractor_ic_scf": "no",
        "loop_over_background_for_closure_relation": "no",
        "background_Nloga": 100000,
        "tol_shooting_deltaF": 0.01,
        "tol_shooting_deltax": 0.01,
        "N_mscf": n,
    }
    if n == 1:
        d.update({
            "n_axion_mscf": "3",
            "theta_ini_mscf": "2.8",
            "theta_prime_ini_mscf": "0.",
        })
    elif n == 2:
        d.update({
            "n_axion_mscf": "3,3",
            "theta_ini_mscf": "2.8,2.8",
            "theta_prime_ini_mscf": "0.,0.",
        })
    return d


def likelihood_block():
    return {
        "act_dr6_cmbonly.ACTDR6CMBonly": {
            "input_file": "dr6_data_cmbonly.fits",
            "lmax_theory": 9000,
            "ell_cuts": {"TT": [600, 8500], "TE": [600, 8500], "EE": [600, 8500]},
            "stop_at_error": True,
            "params": {
                "A_act": {"prior": {"min": 0.5, "max": 1.5}, "ref": BASE["A_act"], "proposal": 0.003},
                "P_act": {"prior": {"min": 0.9, "max": 1.1}, "ref": BASE["P_act"], "proposal": 0.01},
            },
        },
        "candl_like": {
            "external": "__CANDL_EXTERNAL__",
            "data_set_file": str(SPT / "spt_candl_data/SPT3G_D1_TnE_v0/SPT3G_D1_TnE_index.yaml"),
            "variant": "lite",
            "clear_internal_priors": True,
            "additional_args": {},
            "feedback": True,
            "wrapper": None,
        },
    }


def base_params():
    return {
        "logA": {"prior": {"min": 1.61, "max": 3.91}, "ref": BASE["logA"], "proposal": 0.002, "drop": True},
        "A_s": {"value": "lambda logA: 1e-10*np.exp(logA)"},
        "n_s": {"prior": {"min": 0.8, "max": 1.2}, "ref": BASE["n_s"], "proposal": 0.003},
        "H0": {"prior": {"min": 20.0, "max": 100.0}, "ref": BASE["H0"], "proposal": 1.0},
        "omega_b": {"prior": {"min": 0.005, "max": 0.1}, "ref": BASE["omega_b"], "proposal": 0.0001},
        "omega_cdm": {"prior": {"min": 0.001, "max": 0.99}, "ref": BASE["omega_cdm"], "proposal": 0.001},
        "m_ncdm": {"value": 0.06},
        "tau_reio": {"prior": {"dist": "norm", "loc": 0.051, "scale": 0.006}, "ref": BASE["tau_reio"], "proposal": 0.003},
        "Tcal": {"prior": {"min": 0.8, "max": 1.2}, "ref": BASE["Tcal"], "proposal": 0.003},
        "Ecal": {"prior": {"min": 0.8, "max": 1.2}, "ref": BASE["Ecal"], "proposal": 0.003},
    }


def common_config(label):
    return {
        "theory": {"classy": theory_base()},
        "likelihood": likelihood_block(),
        "prior": {
            "cal_dip_prior": "lambda A_act: stats.norm.logpdf(A_act, loc=1.0, scale=0.003)",
            "gaussian_Tcal": "lambda Tcal: stats.norm.logpdf(Tcal, loc=1.0, scale=0.0036)",
        },
        "params": base_params(),
        "packages_path": str(PACKAGES),
        "timing": True,
        "debug": True,
        "output": str(RUNS / label / label),
        "sampler": {"evaluate": {}},
    }


def config_2ax(label, ac1, ac2):
    c = common_config(label)
    c["theory"]["classy"]["extra_args"].update(mscf_args(2))
    c["params"].update({
        "log10_ac1": {"prior": {"min": -3.7500772228, "max": -3.0004340775}, "ref": ac1, "proposal": 0.06, "drop": True},
        "log10_ac2": {"prior": {"min": -4.5000137334, "max": -3.7500772228}, "ref": ac2, "proposal": 0.08, "drop": True},
        "f_ax1": {"prior": {"min": 0.0, "max": 0.3}, "ref": 0.0, "proposal": 0.015, "drop": True},
        "f_ax2": {"prior": {"min": 0.0, "max": 0.3}, "ref": 0.0, "proposal": 0.015, "drop": True},
        "log10_maxion_ac": {"value": 'lambda log10_ac1, log10_ac2: "%.16g,%.16g" % (log10_ac1, log10_ac2)', "derived": False},
        "fraction_maxion_ac": {"value": 'lambda f_ax1, f_ax2: "%.16g,%.16g" % (f_ax1 if f_ax1 > 1e-12 else 1e-12, f_ax2 if f_ax2 > 1e-12 else 1e-12)', "derived": False},
    })
    return c


def config_1ax_raw(label, ac, f=0.0):
    c = common_config(label)
    c["theory"]["classy"]["extra_args"].update(mscf_args(1))
    c["params"].update({
        "log10_maxion_ac": {"prior": {"min": -4.5000137334, "max": -3.0004340775}, "ref": ac, "proposal": 0.08},
        "fraction_maxion_ac": {"prior": {"min": 0.0, "max": 0.3}, "ref": f, "proposal": 0.015},
    })
    return c


def config_1ax_safe(label, ac, f_phys=0.0):
    c = common_config(label)
    c["theory"]["classy"]["extra_args"].update(mscf_args(1))
    c["params"].update({
        "log10_ac_phys": {"prior": {"min": -4.5000137334, "max": -3.0004340775}, "ref": ac, "proposal": 0.08, "drop": True},
        "f_ax_phys": {"prior": {"min": 0.0, "max": 0.3}, "ref": f_phys, "proposal": 0.015, "drop": True},
        "log10_maxion_ac": {"value": "lambda log10_ac_phys: log10_ac_phys", "derived": False},
        "fraction_maxion_ac": {"value": "lambda f_ax_phys: f_ax_phys if f_ax_phys > 1e-12 else 1e-12", "derived": False},
    })
    return c


def config_lcdm_shooting0(label):
    c = common_config(label)
    c["theory"]["classy"]["extra_args"].update(mscf_args(0))
    return c


def config_lcdm_plain(label):
    return common_config(label)


def dump_config(cfg, path):
    text = yaml.safe_dump(cfg, sort_keys=False, width=140)
    text = text.replace("external: __CANDL_EXTERNAL__", "external: !!python/name:candl.interface.CandlCobayaLikelihood ''")
    path.write_text(text)


def parse_log(text):
    def last_float(pattern):
        vals = re.findall(pattern, text, flags=re.MULTILINE)
        if not vals:
            return None
        try:
            return float(vals[-1])
        except ValueError:
            return None
    return {
        "chi2_ACT": last_float(r"chi2_act_dr6_cmbonly\.ACTDR6CMBonly\s*=\s*([^\s]+)"),
        "chi2_SPT": last_float(r"chi2_candl_like\s*=\s*([^\s]+)"),
        "log_likelihood": last_float(r"log-likelihood\s*=\s*([^\s]+)"),
    }


def gaussian_penalty():
    tau = ((BASE["tau_reio"] - 0.051) / 0.006) ** 2
    aact = ((BASE["A_act"] - 1.0) / 0.003) ** 2
    tcal = ((BASE["Tcal"] - 1.0) / 0.0036) ** 2
    return {"tau": tau, "A_act": aact, "Tcal": tcal, "sum": tau + aact + tcal}


def run_case(label, cfg, expected_crash=False, gdb_on_crash=False):
    d = RUNS / label
    d.mkdir(parents=True, exist_ok=True)
    y = d / f"{label}.yaml"
    dump_config(cfg, y)
    log = d / f"{label}.console.log"
    env = os.environ.copy()
    env["PYTHONFAULTHANDLER"] = "1"
    env["COBAYA_PACKAGES_PATH"] = str(PACKAGES)
    cmd = [sys.executable, "-X", "faulthandler", "-m", "cobaya.run", str(y), "--force"]
    start = time.time()
    with log.open("w") as f:
        proc = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT, env=env, cwd=ROOT)
    runtime = time.time() - start
    text = log.read_text(errors="replace")
    parsed = parse_log(text)
    rc = proc.returncode
    crashed = rc < 0 or rc in (128 + signal.SIGSEGV, 139)
    usable = rc == 0 and finite(parsed["chi2_ACT"]) and finite(parsed["chi2_SPT"])
    result = {
        "label": label,
        "exit_code": rc,
        "runtime_seconds": runtime,
        "expected_crash": expected_crash,
        "crashed": crashed,
        "finite_likelihood": usable,
        **parsed,
        "yaml": str(y.relative_to(OUT)),
        "log": str(log.relative_to(OUT)),
    }
    if usable:
        pen = gaussian_penalty()
        result["gaussian_penalties"] = pen
        result["likelihood_chi2"] = parsed["chi2_ACT"] + parsed["chi2_SPT"]
        result["effective_chi2"] = result["likelihood_chi2"] + pen["sum"]
        result["delta_vs_R1"] = result["effective_chi2"] - FROZEN["R1_effective_chi2"]
    if crashed and gdb_on_crash and shutil.which("gdb"):
        gdblog = d / f"{label}.gdb.log"
        gcmd = [
            "gdb", "--batch", "-ex", "set pagination off", "-ex", "run",
            "-ex", "thread apply all bt full", "--args",
            sys.executable, "-X", "faulthandler", "-m", "cobaya.run", str(y), "--force"
        ]
        try:
            with gdblog.open("w") as f:
                subprocess.run(gcmd, stdout=f, stderr=subprocess.STDOUT, env=env, cwd=ROOT, timeout=600)
            result["gdb_log"] = str(gdblog.relative_to(OUT))
        except subprocess.TimeoutExpired:
            result["gdb_status"] = "TIMEOUT_AFTER_600S"
    return result


def git_sha(path):
    try:
        return subprocess.check_output(["git", "-C", str(path), "rev-parse", "HEAD"], text=True).strip()
    except Exception:
        return "NOT DOCUMENTED"


def json_safe(value):
    """Convert non-finite floats to explicit strings so strict JSON remains valid.

    This changes reporting only. Numerical values used by the validation logic above
    remain untouched.
    """
    if isinstance(value, float):
        if math.isnan(value):
            return "NaN"
        if math.isinf(value):
            return "Infinity" if value > 0 else "-Infinity"
        return value
    if isinstance(value, dict):
        return {k: json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(v) for v in value]
    return value


def main():
    cases = []
    for i, (a1, a2) in enumerate(zip(AC1, AC2), 1):
        label = f"2ax_zero_{i}"
        cases.append((label, config_2ax(label, a1, a2), False, False))

    cases.append(("1ax_raw_zero", config_1ax_raw("1ax_raw_zero", AC_ONE[1], 0.0), True, True))

    for i, ac in enumerate(AC_ONE, 1):
        label = f"1ax_safe_zero_{i}"
        cases.append((label, config_1ax_safe(label, ac, 0.0), False, True))

    cases.append(("1ax_safe_1e-8", config_1ax_safe("1ax_safe_1e-8", AC_ONE[1], 1e-8), False, True))
    cases.append(("lcdm_shooting_N0", config_lcdm_shooting0("lcdm_shooting_N0"), False, False))
    cases.append(("lcdm_plain", config_lcdm_plain("lcdm_plain"), False, False))

    results = []
    for label, cfg, expected_crash, gdb in cases:
        print(f"\n=== {label} ===", flush=True)
        results.append(run_case(label, cfg, expected_crash, gdb))

    by = {r["label"]: r for r in results}
    z2 = [r for r in results if r["label"].startswith("2ax_zero_") and r.get("finite_likelihood")]
    z1 = [r for r in results if r["label"].startswith("1ax_safe_zero_") and r.get("finite_likelihood")]
    raw = by["1ax_raw_zero"]
    shoot0 = by["lcdm_shooting_N0"]
    plain = by["lcdm_plain"]

    values = [r["effective_chi2"] for r in z2 + z1 if finite(r.get("effective_chi2"))]
    spread = max(values) - min(values) if len(values) >= 2 else None
    max_r1 = max(abs(v - FROZEN["R1_effective_chi2"]) for v in values) if values else None

    validated = (
        len(z2) == 3 and len(z1) == 3 and plain.get("finite_likelihood") and
        spread is not None and spread <= FROZEN["globality_tolerance"] and
        max_r1 is not None and max_r1 <= FROZEN["globality_tolerance"]
    )

    if validated:
        status = "BOUNDARY_VALIDATED"
        next_action = "Return to Gate-P globality recovery and re-minimize 1ax from the validated physical f=0 boundary before any Gate-P verdict."
    else:
        status = "NEEDS_TECHNICAL_DEBUGGING"
        next_action = "Use the per-run faulthandler/GDB logs to isolate the remaining mAxiCLASS boundary failure. Do not issue a Gate-P verdict."

    summary = {
        "case_id": FROZEN["case_id"],
        "hypothesis_id": FROZEN["hypothesis_id"],
        "internal_result_id": FROZEN["internal_result_id"],
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "technical_status": status,
        "physical_gatep_verdict": "NOT EVALUATED BY THIS PROGRAM",
        "frozen": FROZEN,
        "tests": results,
        "cross_boundary_spread_effective_chi2": spread,
        "max_abs_delta_vs_R1": max_r1,
        "raw_1ax_zero_reproduced_crash": raw.get("crashed", False),
        "lcdm_shooting_N0_finite": shoot0.get("finite_likelihood", False),
        "lcdm_plain_finite": plain.get("finite_likelihood", False),
        "next_action": next_action,
        "interpretation_guard": "A green GitHub workflow means the diagnostic program completed. It is not a physical Gate-P pass."
    }
    (OUT / "BOUNDARY_VALIDATION_RESULT.json").write_text(json.dumps(json_safe(summary), indent=2, allow_nan=False))

    with (OUT / "BOUNDARY_VALIDATION_RESULTS.csv").open("w", newline="") as f:
        cols = ["label", "exit_code", "crashed", "finite_likelihood", "chi2_ACT", "chi2_SPT", "likelihood_chi2", "effective_chi2", "delta_vs_R1", "runtime_seconds"]
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for r in results:
            w.writerow({k: r.get(k) for k in cols})

    lines = [
        "# CASE-Q013 — Gate-P Boundary Validation",
        "",
        f"**Technical status:** `{status}`",
        "",
        "**Physical Gate-P verdict:** NOT EVALUATED BY THIS PROGRAM.",
        "",
        "This workflow tests the numerical nesting boundary only.",
        "",
        "## Core diagnostic",
        "",
        f"- Raw 1ax direct-zero crash reproduced: **{raw.get('crashed', False)}** (exit {raw.get('exit_code')}).",
        f"- Standalone N_mscf=0 shooting finite: **{shoot0.get('finite_likelihood', False)}**.",
        f"- Plain frozen-mAxiCLASS LCDM finite: **{plain.get('finite_likelihood', False)}**.",
        f"- Valid 2ax zero-boundary evaluations: **{len(z2)}/3**.",
        f"- Valid 1ax physical-zero / solver-floor evaluations: **{len(z1)}/3**.",
        f"- Cross-boundary effective-chi2 spread: **{spread if spread is not None else 'NOT AVAILABLE'}**.",
        f"- Max |delta chi2| relative to R1={FROZEN['R1_effective_chi2']}: **{max_r1 if max_r1 is not None else 'NOT AVAILABLE'}**.",
        "",
        "## Individual tests",
        "",
        "| Test | Exit | Finite | chi2 ACT | chi2 SPT | effective chi2 | delta vs R1 |",
        "|---|---:|---|---:|---:|---:|---:|",
    ]
    for r in results:
        lines.append(f"| {r['label']} | {r['exit_code']} | {r.get('finite_likelihood')} | {r.get('chi2_ACT','')} | {r.get('chi2_SPT','')} | {r.get('effective_chi2','')} | {r.get('delta_vs_R1','')} |")
    lines += [
        "",
        "## Acceptance rule",
        "",
        "All 3 two-field zero-boundary evaluations, all 3 one-field physical-zero/solver-floor evaluations, and plain LCDM must be finite. The zero-boundary effective chi2 values must agree mutually and with R1 to <= 0.1.",
        "",
        "## Next action",
        "",
        next_action,
    ]
    (OUT / "BOUNDARY_VALIDATION_REPORT.md").write_text("\n".join(lines) + "\n")

    prov = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "python": sys.version,
        "platform": platform.platform(),
        "mAxiCLASS_actual_commit": git_sha(MAXI),
        "ACT_actual_commit": git_sha(ROOT / "external/DR6-ACT-lite"),
        "SPT_actual_commit": git_sha(SPT),
        "command": "python scripts/gatep_boundary_validation.py",
        "frozen_contract": FROZEN,
        "base_anchor": BASE,
        "notes": [
            "No ACT or SPT lensing likelihood is added.",
            "The 1ax_safe tests keep physical f_ax_phys=0 while regularizing only solver-facing fraction_maxion_ac to 1e-12, matching the established 2ax zero-fraction handling.",
            "lcdm_plain uses the same frozen mAxiCLASS commit without MSCF shooting solely to isolate the N_mscf=0 shooting failure.",
            "The executed A_act Gaussian prior is preserved and remains a separate provenance issue after globality recovery."
        ]
    }
    (PROV / "BOUNDARY_VALIDATION_PROVENANCE.json").write_text(json.dumps(prov, indent=2))
    with (PROV / "pip-freeze.txt").open("w") as f:
        subprocess.run([sys.executable, "-m", "pip", "freeze"], stdout=f, check=False)

    act_file = PACKAGES / "data/ACTDR6CMBonly/v1.0/dr6_data_cmbonly.fits"
    if act_file.exists():
        h = hashlib.sha256()
        with act_file.open("rb") as f:
            for chunk in iter(lambda: f.read(1024 * 1024), b""):
                h.update(chunk)
        (PROV / "ACT_DR6_DATA_SHA256.txt").write_text(f"{h.hexdigest()}  {act_file}\n")

    print(json.dumps({"technical_status": status, "artifact_dir": str(OUT), "next_action": next_action}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
