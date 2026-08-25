#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import importlib
import importlib.util
import json
import os
import platform
import subprocess
import sys
from pathlib import Path
from typing import Any

ROOT = Path.cwd()
OUT = ROOT / "gatep_environment_recovery"
OUT.mkdir(parents=True, exist_ok=True)

EXPECTED = {
    "maxiclass": os.environ.get(
        "MAXICLASS_COMMIT",
        "b5a3af9818d04d2fe696e54972545e8c6805f4ea",
    ),
    "act": os.environ.get(
        "ACT_COMMIT",
        "0e0cd2c703c62a0e980470b572602233b27750e1",
    ),
    "spt": os.environ.get(
        "SPT_COMMIT",
        "2cec6e762a8c540484dd5acafc529f0035856350",
    ),
    "candl_like": "2.0.3",
    "cobaya": "3.6.2",
}

REPOS = {
    "maxiclass": ROOT / "external/mAxiCLASS",
    "act": ROOT / "external/DR6-ACT-lite",
    "spt": ROOT / "external/spt_candl_data",
}


def run(cmd: list[str]) -> str:
    return subprocess.check_output(
        cmd, text=True, stderr=subprocess.STDOUT
    ).strip()


def git_head(path: Path) -> str:
    return run(["git", "-C", str(path), "rev-parse", "HEAD"])


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def version(dist_name: str) -> str:
    from importlib.metadata import version as dist_version
    return dist_version(dist_name)


def module_info(name: str) -> dict[str, Any]:
    spec = importlib.util.find_spec(name)
    if spec is None:
        raise RuntimeError(f"Module {name!r} is not importable")
    module = importlib.import_module(name)
    file = getattr(module, "__file__", None)
    return {
        "name": name,
        "file": str(Path(file).resolve()) if file else None,
        "version": getattr(module, "__version__", None),
    }


def package_root_for_class() -> Path:
    # mAxiCLASS setup.py installs the data package with distribution package
    # name "class". We locate it without writing `import class`, which is invalid
    # Python syntax.
    spec = importlib.util.find_spec("class")
    if spec is None:
        raise RuntimeError(
            "The mAxiCLASS package-data module 'class' is missing. "
            "The frozen source was built but not installed correctly."
        )
    if spec.submodule_search_locations:
        return Path(next(iter(spec.submodule_search_locations))).resolve()
    if spec.origin:
        return Path(spec.origin).resolve().parent
    raise RuntimeError("Could not resolve installed mAxiCLASS data package root")


def assert_equal(label: str, got: str, expected: str):
    if got != expected:
        raise RuntimeError(f"{label}: expected {expected}, got {got}")


def main() -> int:
    result: dict[str, Any] = {
        "status": "PASS",
        "python": sys.version,
        "platform": platform.platform(),
        "executable": sys.executable,
        "checks": {},
    }

    try:
        # Frozen source commits.
        for key, path in REPOS.items():
            if not path.exists():
                raise RuntimeError(f"Missing frozen repository: {path}")
            head = git_head(path)
            result["checks"][f"{key}_commit"] = head
            assert_equal(f"{key} commit", head, EXPECTED[key])

        # Frozen package versions.
        cobaya_v = version("cobaya")
        candl_v = version("candl-like")
        result["checks"]["cobaya_version"] = cobaya_v
        result["checks"]["candl_like_version"] = candl_v
        assert_equal("Cobaya version", cobaya_v, EXPECTED["cobaya"])
        assert_equal("candl-like version", candl_v, EXPECTED["candl_like"])

        # Crucial imports.
        classy = module_info("classy")
        candl = module_info("candl")
        candl_data = module_info("candl_data")
        act = module_info("act_dr6_cmbonly")
        result["checks"]["classy"] = classy
        result["checks"]["candl"] = candl
        result["checks"]["candl_data"] = candl_data
        result["checks"]["act_dr6_cmbonly"] = act

        # Verify that classy is from the environment we just installed, not a
        # stale unrelated package. The exact .so can live in site-packages, but
        # the package-data tree must come from the same frozen installation.
        classy_file = Path(classy["file"]).resolve()
        if not classy_file.exists():
            raise RuntimeError(f"classy import path does not exist: {classy_file}")

        class_root = package_root_for_class()
        result["checks"]["class_data_root"] = str(class_root)

        bbn = class_root / "external" / "bbn" / "sBBN_2025.dat"
        result["checks"]["sBBN_2025_path"] = str(bbn)
        if not bbn.is_file():
            raise RuntimeError(
                "Required frozen CLASS BBN data file is missing: "
                f"{bbn}. Do not run Gate-P."
            )
        if bbn.stat().st_size <= 0:
            raise RuntimeError(f"BBN data file is empty: {bbn}")
        with bbn.open("rb") as f:
            f.read(1)
        result["checks"]["sBBN_2025_sha256"] = sha256(bbn)
        result["checks"]["sBBN_2025_size"] = bbn.stat().st_size

        # Frozen ACT data.
        packages = Path(
            os.environ.get("COBAYA_PACKAGES_PATH", ROOT / ".cobaya")
        ).resolve()
        act_file = (
            packages
            / "data"
            / "ACTDR6CMBonly"
            / "v1.0"
            / "dr6_data_cmbonly.fits"
        )
        result["checks"]["act_data_file"] = str(act_file)
        if not act_file.is_file() or act_file.stat().st_size <= 0:
            raise RuntimeError(f"Frozen ACT DR6 data missing: {act_file}")
        result["checks"]["act_data_sha256"] = sha256(act_file)

        # Verify editable SPT package really points into frozen checkout.
        candl_data_file = Path(candl_data["file"]).resolve()
        spt_root = REPOS["spt"].resolve()
        try:
            candl_data_file.relative_to(spt_root)
        except ValueError:
            raise RuntimeError(
                "candl_data is importable but is not loaded from the frozen "
                f"SPT checkout. Loaded: {candl_data_file}; expected under {spt_root}"
            )

        # Record selected package state.
        result["checks"]["pip_freeze"] = run(
            [sys.executable, "-m", "pip", "freeze"]
        ).splitlines()

    except Exception as exc:
        result["status"] = "FAIL"
        result["error"] = repr(exc)

    (OUT / "preflight.json").write_text(
        json.dumps(result, indent=2, sort_keys=True)
    )

    lines = [
        "# Gate-P Environment Recovery Preflight",
        "",
        f"Status: **{result['status']}**",
        "",
    ]
    if result["status"] == "FAIL":
        lines += [f"Error: `{result.get('error')}`", ""]
    for key, value in result.get("checks", {}).items():
        if key == "pip_freeze":
            continue
        lines.append(f"- **{key}:** `{value}`")
    (OUT / "PREFLIGHT_REPORT.md").write_text("\n".join(lines) + "\n")

    print(json.dumps(result, indent=2, sort_keys=True))
    return 0 if result["status"] == "PASS" else 2


if __name__ == "__main__":
    raise SystemExit(main())
