#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import importlib
import json
import math
import os
from pathlib import Path
import sys
import traceback

REPO = Path(__file__).resolve().parents[1]
MAXI = REPO / "external" / "mAxiCLASS"
ACT = REPO / "external" / "DR6-ACT-lite"
SPT = REPO / "external" / "spt_candl_data"
EXPECTED_MAXI_COMMIT = "b5a3af9818d04d2fe696e54972545e8c6805f4ea"
EXPECTED_ACT_COMMIT = "0e0cd2c703c62a0e980470b572602233b27750e1"
EXPECTED_SPT_COMMIT = "2cec6e762a8c540484dd5acafc529f0035856350"
EXPECTED_BBN_SHA256 = "321b3adc1ffd711d67be1f5a56948191a1121bf1e4864aa0212c4aa556ee670e"


def die(msg: str, code: int = 2) -> None:
    print(f"FATAL: {msg}", file=sys.stderr, flush=True)
    raise SystemExit(code)


def git_commit(path: Path) -> str:
    import subprocess
    return subprocess.check_output(["git", "-C", str(path), "rev-parse", "HEAD"], text=True).strip()


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def prepare_import_routing() -> None:
    # Critical v8 fix: never give Cobaya theory.classy.path. Put the frozen repos
    # first on sys.path/PYTHONPATH so Cobaya's normal import resolves the same binary
    # as this preflight import.
    for p in (MAXI, ACT, SPT):
        if not p.exists():
            die(f"Missing frozen dependency: {p}")
    wanted = [str(MAXI), str(ACT), str(SPT)]
    sys.path[:] = wanted + [p for p in sys.path if p not in wanted]
    os.environ["PYTHONPATH"] = os.pathsep.join(wanted + [os.environ.get("PYTHONPATH", "")])


def validate_frozen_state() -> dict:
    maxi_commit = git_commit(MAXI)
    act_commit = git_commit(ACT)
    spt_commit = git_commit(SPT)
    if maxi_commit != EXPECTED_MAXI_COMMIT:
        die(f"mAxiCLASS commit drift: {maxi_commit}")
    if act_commit != EXPECTED_ACT_COMMIT:
        die(f"ACT commit drift: {act_commit}")
    if spt_commit != EXPECTED_SPT_COMMIT:
        die(f"SPT commit drift: {spt_commit}")

    bbn = MAXI / "external" / "bbn" / "sBBN_2025.dat"
    if not bbn.is_file():
        die(f"Missing frozen BBN table: {bbn}")
    bbn_sha = sha256(bbn)
    if bbn_sha != EXPECTED_BBN_SHA256:
        die(f"sBBN_2025.dat SHA256 drift: {bbn_sha}")

    classy = importlib.import_module("classy")
    classy_file = Path(classy.__file__).resolve()
    if MAXI.resolve() not in classy_file.parents:
        die(f"Frozen classy lost import routing: {classy_file}")

    import candl
    import spt_candl_data
    import act_dr6_cmbonly
    import cobaya

    return {
        "python": sys.executable,
        "classy_file_preflight": str(classy_file),
        "candl_file": str(Path(candl.__file__).resolve()),
        "candl_version": getattr(candl, "__version__", None),
        "cobaya_version": getattr(cobaya, "__version__", None),
        "spt_file": str(Path(spt_candl_data.__file__).resolve()),
        "act_file": str(Path(act_dr6_cmbonly.__file__).resolve()),
        "mAxiCLASS_commit": maxi_commit,
        "ACT_commit": act_commit,
        "SPT_commit": spt_commit,
        "sBBN_2025_sha256": bbn_sha,
    }


def build_info(fraction: float, output_prefix: str) -> dict:
    from scipy import stats
    from candl.interface import CandlCobayaLikelihood

    # IMPORTANT: no "path" key under theory.classy. That was the v1 failure.
    info = {
        "theory": {
            "classy": {
                "ignore_obsolete": True,
                "extra_args": {
                    "non linear": "hmcode",
                    "N_ur": 2.0328,
                    "N_ncdm": 1,
                    "T_ncdm": 0.71611,
                    "lensing": "yes",
                    "P_k_max_h/Mpc": 1.0,
                    "do_shooting": "yes",
                    "do_shooting_mscf": "yes",
                    "attractor_ic_scf": "no",
                    "loop_over_background_for_closure_relation": "no",
                    "background_Nloga": 100000,
                    "tol_shooting_deltaF": 0.01,
                    "tol_shooting_deltax": 0.01,
                    "N_mscf": 1,
                    "n_axion_mscf": "3",
                    "theta_ini_mscf": "2.8",
                    "theta_prime_ini_mscf": "0.",
                },
            }
        },
        "likelihood": {
            "act_dr6_cmbonly.ACTDR6CMBonly": {
                "input_file": "dr6_data_cmbonly.fits",
                "lmax_theory": 9000,
                "ell_cuts": {"TT": [600, 8500], "TE": [600, 8500], "EE": [600, 8500]},
                "stop_at_error": True,
                "params": {
                    "A_act": {"prior": {"min": 0.5, "max": 1.5}, "ref": 1.0011501, "proposal": 0.003},
                    "P_act": {"prior": {"min": 0.9, "max": 1.1}, "ref": 1.0016445, "proposal": 0.01},
                },
            },
            "candl_like": {
                "external": CandlCobayaLikelihood,
                "data_set_file": str(SPT / "spt_candl_data" / "SPT3G_D1_TnE_v0" / "SPT3G_D1_TnE_index.yaml"),
                "variant": "lite",
                "clear_internal_priors": True,
                "additional_args": {},
                "feedback": True,
                "wrapper": None,
            },
        },
        "prior": {
            "cal_dip_prior": lambda A_act: stats.norm.logpdf(A_act, loc=1.0, scale=0.003),
            "gaussian_Tcal": lambda Tcal: stats.norm.logpdf(Tcal, loc=1.0, scale=0.0036),
        },
        "params": {
            "logA": {"prior": {"min": 1.61, "max": 3.91}, "ref": 3.0516831, "proposal": 0.002, "drop": True},
            "A_s": {"value": "lambda logA: 1e-10*np.exp(logA)"},
            "n_s": {"prior": {"min": 0.8, "max": 1.2}, "ref": 0.97006392, "proposal": 0.003},
            "H0": {"prior": {"min": 20.0, "max": 100.0}, "ref": 67.003038, "proposal": 1.0},
            "omega_b": {"prior": {"min": 0.005, "max": 0.1}, "ref": 0.022462565, "proposal": 0.0001},
            "omega_cdm": {"prior": {"min": 0.001, "max": 0.99}, "ref": 0.12107093, "proposal": 0.001},
            "m_ncdm": {"value": 0.06},
            "tau_reio": {"prior": {"dist": "norm", "loc": 0.051, "scale": 0.006}, "ref": 0.056456949, "proposal": 0.003},
            "Tcal": {"prior": {"min": 0.8, "max": 1.2}, "ref": 0.99859358, "proposal": 0.003},
            "Ecal": {"prior": {"min": 0.8, "max": 1.2}, "ref": 1.0066561, "proposal": 0.003},
            "log10_maxion_ac": {"prior": {"min": -4.5000137334, "max": -3.0004340775}, "ref": -3.75022390545, "proposal": 0.08},
            "fraction_maxion_ac": {"prior": {"min": 0.0, "max": 0.3}, "ref": fraction, "proposal": 0.015},
        },
        "packages_path": str(REPO / ".cobaya"),
        "timing": True,
        "debug": True,
        "output": output_prefix,
    }
    return info


def jsonable(x):
    if isinstance(x, (str, int, bool)) or x is None:
        return x
    if isinstance(x, float):
        return x if math.isfinite(x) else str(x)
    if isinstance(x, dict):
        return {str(k): jsonable(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [jsonable(v) for v in x]
    try:
        import numpy as np
        if isinstance(x, np.ndarray):
            return [jsonable(v) for v in x.tolist()]
        if isinstance(x, np.generic):
            return jsonable(x.item())
    except Exception:
        pass
    if hasattr(x, "__dict__"):
        return {k: jsonable(v) for k, v in vars(x).items() if not k.startswith("_")}
    return repr(x)


def run_one(fraction: float, outdir: Path) -> int:
    outdir.mkdir(parents=True, exist_ok=True)
    provenance = validate_frozen_state()
    (outdir / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")

    # Cobaya is imported only after frozen classy has already been verified and
    # placed first on sys.path. This is the core routing guarantee.
    from cobaya.model import get_model

    info = build_info(fraction, str(outdir / "cobaya" / "run"))
    result = {
        "fraction_maxion_ac": fraction,
        "status": "STARTED",
        "provenance": provenance,
    }
    try:
        with get_model(info) as model:
            classy_component = model.theory["classy"]
            cobaya_classy_file = Path(classy_component.classy_module.__file__).resolve()
            result["classy_file_cobaya"] = str(cobaya_classy_file)
            if cobaya_classy_file != Path(provenance["classy_file_preflight"]):
                die(
                    "Cobaya/preflight classy mismatch: "
                    f"{cobaya_classy_file} != {provenance['classy_file_preflight']}"
                )

            point = {
                "logA": 3.0516831,
                "n_s": 0.97006392,
                "H0": 67.003038,
                "omega_b": 0.022462565,
                "omega_cdm": 0.12107093,
                "tau_reio": 0.056456949,
                "Tcal": 0.99859358,
                "Ecal": 1.0066561,
                "log10_maxion_ac": -3.75022390545,
                "fraction_maxion_ac": fraction,
                "A_act": 1.0011501,
                "P_act": 1.0016445,
            }
            post = model.logposterior(point, as_dict=True)
            result["posterior"] = jsonable(post)
            vals = []
            def collect(obj):
                if isinstance(obj, dict):
                    for v in obj.values(): collect(v)
                elif isinstance(obj, (list, tuple)):
                    for v in obj: collect(v)
                elif isinstance(obj, (int, float)) and not isinstance(obj, bool):
                    vals.append(float(obj))
            collect(result["posterior"])
            result["all_numeric_outputs_finite"] = all(math.isfinite(v) for v in vals)
            result["status"] = "FINITE" if result["all_numeric_outputs_finite"] else "NONFINITE"
    except BaseException as exc:
        result["status"] = "PYTHON_EXCEPTION"
        result["exception_type"] = type(exc).__name__
        result["exception"] = str(exc)
        result["traceback"] = traceback.format_exc()
        (outdir / "result.json").write_text(json.dumps(result, indent=2) + "\n")
        raise

    (outdir / "result.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2), flush=True)
    return 0 if result["status"] == "FINITE" else 3


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--fraction", type=float, required=True)
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()
    prepare_import_routing()
    return run_one(args.fraction, Path(args.outdir).resolve())


if __name__ == "__main__":
    raise SystemExit(main())
