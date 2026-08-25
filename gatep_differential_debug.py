#!/usr/bin/env python3
from __future__ import annotations

import difflib
import hashlib
import json
import math
import os
import platform
import re
import shutil
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import yaml

ROOT = Path.cwd()
OUT = ROOT / "gatep_differential_debug"
RUNS = OUT / "runs"
PROV = OUT / "provenance"
OUT.mkdir(exist_ok=True)
RUNS.mkdir(exist_ok=True)
PROV.mkdir(exist_ok=True)

sys.path.insert(0, str(ROOT / "scripts"))
import gatep_globality_recovery as GR
import gatep_boundary_validation as BV

R1 = 317.4238643933
PREVIOUS_2AX_ZERO = 317.4201095748994
TOL = 0.1
AC1 = -3.37525565015
AC2 = -4.1250454781


def finite(x):
    return isinstance(x, (int, float)) and math.isfinite(x)


def sha256_file(path: Path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def git_sha(path: Path):
    try:
        return subprocess.check_output(
            ["git", "-C", str(path), "rev-parse", "HEAD"],
            text=True
        ).strip()
    except Exception:
        return "NOT_DOCUMENTED"


def cmd_output(cmd):
    try:
        return subprocess.check_output(
            cmd, text=True, stderr=subprocess.STDOUT
        ).strip()
    except Exception as exc:
        return f"ERROR: {exc!r}"


def dump_cfg(cfg, path: Path):
    text = yaml.safe_dump(cfg, sort_keys=False, width=140)
    text = text.replace(
        "external: __CANDL_EXTERNAL__",
        "external: !!python/name:candl.interface.CandlCobayaLikelihood ''"
    )
    text = text.replace(
        "external: '__CANDL_EXTERNAL__'",
        "external: !!python/name:candl.interface.CandlCobayaLikelihood ''"
    )
    path.write_text(text)


def find_first_error(text: str):
    patterns = [
        r"(?im)^.*(?:serious error|traceback|exception|segmentation fault|error in classy|classyerror|computationerror|likelihooderror|nan|infinity|\binf\b).*$",
        r"(?im)^.*(?:failed|failure|could not|cannot|invalid).*$",
    ]
    hits = []
    for pat in patterns:
        for m in re.finditer(pat, text):
            line = m.group(0).strip()
            if line and line not in hits:
                hits.append(line)
            if len(hits) >= 16:
                return hits
    return hits


def parse_console(text: str):
    def last_float(pattern):
        vals = re.findall(pattern, text, flags=re.MULTILINE)
        if not vals:
            return None
        try:
            return float(vals[-1])
        except Exception:
            return None

    return {
        "chi2_ACT_log": last_float(
            r"chi2_act_dr6_cmbonly\.ACTDR6CMBonly\s*=\s*([^\s]+)"
        ),
        "chi2_SPT_log": last_float(
            r"chi2_candl_like\s*=\s*([^\s]+)"
        ),
        "log_likelihood": last_float(
            r"log-likelihood\s*=\s*([^\s]+)"
        ),
        "first_error_lines": find_first_error(text),
    }


def result_from_table(prefix: Path):
    candidates = [
        Path(str(prefix) + ".1.txt"),
        Path(str(prefix) + ".minimum.txt"),
    ]
    for p in candidates:
        if not p.exists():
            continue
        try:
            row = GR.parse_table(p)
            eff = GR.effective(row)
            return {
                "result_file": str(p.relative_to(OUT)),
                "table": row,
                **eff,
            }
        except Exception as exc:
            return {
                "result_file": str(p.relative_to(OUT)),
                "table_parse_error": repr(exc),
            }
    return {"result_file": None}


def run_case(label: str, cfg: dict, runner: str):
    d = RUNS / label
    d.mkdir(parents=True, exist_ok=True)
    prefix = d / label

    cfg["output"] = str(prefix)
    cfg["sampler"] = {"evaluate": {}}

    yml = d / f"{label}.yaml"
    dump_cfg(cfg, yml)

    log = d / f"{label}.console.log"
    env = os.environ.copy()
    env["PYTHONFAULTHANDLER"] = "1"
    env["COBAYA_PACKAGES_PATH"] = str(ROOT / ".cobaya")

    if runner == "cli":
        cmd = ["cobaya-run", str(yml), "--force"]
    elif runner == "module":
        cmd = [
            sys.executable, "-X", "faulthandler",
            "-m", "cobaya.run", str(yml), "--force"
        ]
    else:
        raise ValueError(runner)

    t0 = time.time()
    with log.open("w") as f:
        proc = subprocess.run(
            cmd,
            stdout=f,
            stderr=subprocess.STDOUT,
            env=env,
            cwd=ROOT,
        )
    runtime = time.time() - t0

    text = log.read_text(errors="replace")
    console = parse_console(text)
    table = result_from_table(prefix)

    eff = table.get("effective_chi2")
    finite_table = finite(eff)
    finite_logs = (
        finite(console.get("chi2_ACT_log"))
        and finite(console.get("chi2_SPT_log"))
    )
    usable = proc.returncode == 0 and (finite_table or finite_logs)

    result = {
        "label": label,
        "runner": runner,
        "exit_code": proc.returncode,
        "runtime_seconds": runtime,
        "usable_finite": usable,
        "yaml": str(yml.relative_to(OUT)),
        "yaml_sha256": sha256_file(yml),
        "console_log": str(log.relative_to(OUT)),
        **console,
        **table,
    }

    if finite(eff):
        result["delta_vs_R1"] = eff - R1
        result["delta_vs_previous_2ax_zero"] = eff - PREVIOUS_2AX_ZERO

    return result


def config_globality_known_good():
    cfg = GR.base("2ax")
    refs = {
        **GR.LCDM,
        "log10_ac1": AC1,
        "log10_ac2": AC2,
        "f_ax1": 0.0,
        "f_ax2": 0.0,
    }
    GR.set_refs(cfg, refs)
    return cfg


def config_boundary_2ax():
    return BV.config_2ax("boundary_template", AC1, AC2)


def config_plain_lcdm():
    return BV.config_lcdm_plain("plain_template")


def normalized_yaml(cfg):
    text = yaml.safe_dump(cfg, sort_keys=True, width=160)
    lines = []
    for line in text.splitlines():
        s = line.strip()
        if s.startswith("output:"):
            continue
        if s.startswith("debug:"):
            continue
        if s.startswith("timing:"):
            continue
        lines.append(line)
    return "\n".join(lines) + "\n"


def write_config_diff(a, b):
    aa = normalized_yaml(a).splitlines(keepends=True)
    bb = normalized_yaml(b).splitlines(keepends=True)
    diff = "".join(
        difflib.unified_diff(
            aa,
            bb,
            fromfile="globality_known_good",
            tofile="boundary_validation",
            n=5,
        )
    )
    p = PROV / "KNOWN_GOOD_VS_BOUNDARY_CONFIG.diff"
    p.write_text(
        diff if diff
        else "NO SEMANTIC CONFIG DIFFERENCE AFTER NORMALIZATION\n"
    )
    return str(p.relative_to(OUT))


def gdb_1ax_raw():
    d = RUNS / "1ax_raw_zero_gdb"
    d.mkdir(parents=True, exist_ok=True)

    cfg = BV.config_1ax_raw(
        "1ax_raw_zero_gdb",
        -3.75022390545,
        0.0,
    )
    cfg["output"] = str(d / "1ax_raw_zero_gdb")
    cfg["sampler"] = {"evaluate": {}}

    yml = d / "1ax_raw_zero_gdb.yaml"
    dump_cfg(cfg, yml)

    log = d / "1ax_raw_zero_gdb.gdb.log"

    if not shutil.which("gdb"):
        return {"status": "GDB_NOT_AVAILABLE"}

    env = os.environ.copy()
    env["PYTHONFAULTHANDLER"] = "1"
    env["COBAYA_PACKAGES_PATH"] = str(ROOT / ".cobaya")

    cmd = [
        "gdb",
        "--batch",
        "-ex", "set pagination off",
        "-ex", "set print frame-arguments all",
        "-ex", "run",
        "-ex", "thread apply all bt full",
        "--args",
        sys.executable,
        "-X", "faulthandler",
        "-m", "cobaya.run",
        str(yml),
        "--force",
    ]

    try:
        with log.open("w") as f:
            proc = subprocess.run(
                cmd,
                stdout=f,
                stderr=subprocess.STDOUT,
                env=env,
                cwd=ROOT,
                timeout=900,
            )

        text = log.read_text(errors="replace")
        return {
            "status": "DONE",
            "gdb_exit_code": proc.returncode,
            "log": str(log.relative_to(OUT)),
            "first_error_lines": find_first_error(text),
        }

    except subprocess.TimeoutExpired:
        return {
            "status": "TIMEOUT_900S",
            "log": str(log.relative_to(OUT)),
        }


def environment_provenance():
    data = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "github_repository": os.getenv(
            "GITHUB_REPOSITORY", "NOT_DOCUMENTED"
        ),
        "github_sha": os.getenv(
            "GITHUB_SHA", "NOT_DOCUMENTED"
        ),
        "github_run_id": os.getenv(
            "GITHUB_RUN_ID", "NOT_DOCUMENTED"
        ),
        "github_run_attempt": os.getenv(
            "GITHUB_RUN_ATTEMPT", "NOT_DOCUMENTED"
        ),
        "runner_os": os.getenv(
            "RUNNER_OS", "NOT_DOCUMENTED"
        ),
        "runner_arch": os.getenv(
            "RUNNER_ARCH", "NOT_DOCUMENTED"
        ),
        "platform": platform.platform(),
        "python": sys.version,
        "pip_freeze": cmd_output(
            [sys.executable, "-m", "pip", "freeze"]
        ),
        "cobaya_run_path": shutil.which("cobaya-run") or "NOT_FOUND",
        "mAxiCLASS_commit": git_sha(
            ROOT / "external" / "mAxiCLASS"
        ),
        "ACT_commit": git_sha(
            ROOT / "external" / "DR6-ACT-lite"
        ),
        "SPT_commit": git_sha(
            ROOT / "external" / "spt_candl_data"
        ),
    }

    (PROV / "environment.json").write_text(
        json.dumps(data, indent=2)
    )
    (PROV / "pip-freeze.txt").write_text(
        data["pip_freeze"] + "\n"
    )
    return data


def classify(results):
    by = {r["label"]: r for r in results}

    cli = by["known_good_globality_cli"]
    mod = by["known_good_globality_module"]
    bnd = by["boundary_2ax_module"]
    lcdm = by["plain_lcdm_module"]

    if (
        cli["usable_finite"]
        and finite(cli.get("effective_chi2"))
        and abs(
            cli["effective_chi2"] - PREVIOUS_2AX_ZERO
        ) <= TOL
    ):
        if not mod["usable_finite"]:
            return "ROOT_CAUSE_IDENTIFIED_RUNNER_INVOCATION"
        if not bnd["usable_finite"]:
            return "ROOT_CAUSE_SCOPE_CONFIG_OR_WRAPPER"
        if lcdm["usable_finite"]:
            return "REFERENCE_RECOVERED"
        return "REFERENCE_2AX_RECOVERED_LCDM_STILL_FAILS"

    if (
        not cli["usable_finite"]
        and not mod["usable_finite"]
        and not bnd["usable_finite"]
    ):
        return "SHARED_ENVIRONMENT_OR_PIPELINE_FAILURE"

    return "MIXED_DIAGNOSTIC_STATE"


def render_report(status, results, diff_path, gdb, env):
    rows = []
    for r in results:
        rows.append(
            f"| {r['label']} | {r['runner']} | "
            f"{r['exit_code']} | {r['usable_finite']} | "
            f"{r.get('chi2_ACT_log')} | "
            f"{r.get('chi2_SPT_log')} | "
            f"{r.get('effective_chi2')} | "
            f"{r.get('delta_vs_R1')} |"
        )

    errors = []
    for r in results:
        errors.append(f"### {r['label']}")
        hits = r.get("first_error_lines") or []
        if hits:
            errors.extend(
                f"- `{x[:500]}`" for x in hits
            )
        else:
            errors.append(
                "- No error-pattern line extracted; "
                "inspect the full console log."
            )

    report = f"""# CASE-Q013 — Gate-P Differential Debug

**Technical status: {status}**

Physical Gate-P verdict: **NOT EVALUATED BY THIS PROGRAM.**

This workflow performs a same-run differential comparison between the previously successful Globality Recovery configuration path and the newer Boundary Validation path.

## Frozen references

- R1 effective chi2: `{R1}`
- Previous successful 2ax zero-boundary effective chi2: `{PREVIOUS_2AX_ZERO}`
- Required reproduction tolerance: `{TOL}`
- mAxiCLASS commit: `{env['mAxiCLASS_commit']}`
- ACT commit: `{env['ACT_commit']}`
- SPT commit: `{env['SPT_commit']}`

## A/B results

| Test | Runner | Exit | Finite usable | ACT chi2(log) | SPT chi2(log) | Effective chi2(table) | Delta vs R1 |
|---|---:|---:|---:|---:|---:|---:|---:|
{chr(10).join(rows)}

## Interpretation key

- `REFERENCE_RECOVERED`: old known-good 2ax path and plain LCDM are finite again.
- `ROOT_CAUSE_IDENTIFIED_RUNNER_INVOCATION`: same known-good configuration works through `cobaya-run` but fails through `python -m cobaya.run`.
- `ROOT_CAUSE_SCOPE_CONFIG_OR_WRAPPER`: known-good configuration works in the current runner while Boundary Validation configuration does not.
- `SHARED_ENVIRONMENT_OR_PIPELINE_FAILURE`: even the old known-good configuration fails in the same current environment.

## Configuration diff

Artifact:

`{diff_path}`

## Earliest extracted error lines

{chr(10).join(errors)}

## Raw 1ax zero GDB diagnostic

```json
{json.dumps(gdb, indent=2)}
```

## Acceptance rule

A shell exit code of zero is **not** enough.

Any NaN, inf, -inf, crash, missing likelihood component, missing result table, or non-finite objective is:

**INVALID FOR NESTING**

Do not issue a Gate-P physical verdict from this workflow unless the reference chain is first recovered.
"""

    (OUT / "DIFFERENTIAL_DEBUG_REPORT.md").write_text(report)
    return report


def main():
    env = environment_provenance()

    known_cli = config_globality_known_good()
    known_module = config_globality_known_good()
    boundary = config_boundary_2ax()
    plain = config_plain_lcdm()

    diff_path = write_config_diff(
        config_globality_known_good(),
        config_boundary_2ax(),
    )

    results = []

    print("=== known_good_globality_cli ===", flush=True)
    results.append(
        run_case(
            "known_good_globality_cli",
            known_cli,
            "cli",
        )
    )

    print("=== known_good_globality_module ===", flush=True)
    results.append(
        run_case(
            "known_good_globality_module",
            known_module,
            "module",
        )
    )

    print("=== boundary_2ax_module ===", flush=True)
    results.append(
        run_case(
            "boundary_2ax_module",
            boundary,
            "module",
        )
    )

    print("=== plain_lcdm_module ===", flush=True)
    results.append(
        run_case(
            "plain_lcdm_module",
            plain,
            "module",
        )
    )

    print("=== 1ax_raw_zero_gdb ===", flush=True)
    gdb = gdb_1ax_raw()

    status = classify(results)

    payload = {
        "case_id": "CASE-Q013",
        "hypothesis_id": "HYP-004C",
        "internal_result_id":
            "BV-Q013-GATEP-DIFFERENTIAL-DEBUG-001",
        "status": status,
        "physical_verdict": "NOT_EVALUATED",
        "R1_effective_chi2": R1,
        "previous_2ax_zero_effective_chi2":
            PREVIOUS_2AX_ZERO,
        "tolerance": TOL,
        "results": results,
        "config_diff": diff_path,
        "gdb_1ax_raw_zero": gdb,
        "environment": env,
    }

    (OUT / "DIFFERENTIAL_DEBUG_RESULT.json").write_text(
        json.dumps(payload, indent=2, allow_nan=False)
    )

    report = render_report(
        status,
        results,
        diff_path,
        gdb,
        env,
    )

    print(report)

    # Diagnostic workflows intentionally return 0 so the complete
    # artifact is always uploaded. The JSON/Markdown status carries
    # the technical verdict.
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
