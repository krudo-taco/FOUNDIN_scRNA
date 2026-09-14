#!/usr/bin/env python3
"""Write a compact dependency and runtime environment report."""
import argparse
import json
import os
import platform
import subprocess
import sys
from pathlib import Path


def module_version(name):
    try:
        module = __import__(name)
        return getattr(module, "__version__", "unknown")
    except Exception as exc:
        return f"unavailable: {exc}"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    rscript = subprocess.run(["Rscript", "--version"], capture_output=True, text=True)
    r_versions = subprocess.run(
        [
            "Rscript",
            "-e",
            "local_libs <- c('r-library', '../05-claim-support/r-library'); .libPaths(c(local_libs[file.exists(local_libs)], .libPaths())); pkgs <- c('qs','stringfish','RApiSerialize','RcppParallel','data.table','Matrix','Seurat','limma'); cat(jsonlite::toJSON(setNames(lapply(pkgs, function(p) if (requireNamespace(p, quietly=TRUE)) as.character(utils::packageVersion(p)) else 'NOT_INSTALLED'), pkgs), auto_unbox=TRUE))",
        ],
        capture_output=True,
        text=True,
    )
    report = {
        "python": sys.version,
        "platform": platform.platform(),
        "executable": sys.executable,
        "threads": {
            "OPENBLAS_NUM_THREADS": os.environ.get("OPENBLAS_NUM_THREADS"),
            "OMP_NUM_THREADS": os.environ.get("OMP_NUM_THREADS"),
            "MKL_NUM_THREADS": os.environ.get("MKL_NUM_THREADS"),
        },
        "python_packages": {
            "numpy": module_version("numpy"),
            "pandas": module_version("pandas"),
            "scipy": module_version("scipy"),
            "statsmodels": module_version("statsmodels"),
            "h5py": module_version("h5py"),
        },
        "rscript_version": (rscript.stderr or rscript.stdout).strip(),
        "r_packages": json.loads(r_versions.stdout) if r_versions.returncode == 0 and r_versions.stdout.strip() else {"status": "unavailable"},
    }
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
