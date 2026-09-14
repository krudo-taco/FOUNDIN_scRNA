"""Portable staging helpers for rerun wrappers."""
import hashlib
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path


PACKAGE_ROOT = Path(__file__).resolve().parents[2]


def load_json(path):
    return json.loads(Path(path).read_text())


def resolve_path(path, base=PACKAGE_ROOT):
    path = Path(path)
    if path.is_absolute():
        return path
    return (base / path).resolve()


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def ensure_clean_dir(path):
    path = Path(path)
    if path.exists():
        shutil.rmtree(path)
    path.mkdir(parents=True)
    return path


def copy_tree(src, dest):
    src, dest = Path(src), Path(dest)
    if dest.exists():
        shutil.rmtree(dest)
    ignore = shutil.ignore_patterns("__pycache__", "*.pyc", "work", "data")
    shutil.copytree(src, dest, ignore=ignore)


def copy_file(src, dest):
    dest = Path(dest)
    dest.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dest)


def verify_hash_records(records, search_roots=None):
    search_roots = [Path(p) for p in (search_roots or [])]
    checked = []
    for record in records:
        candidates = []
        if search_roots:
            for root in search_roots:
                candidates.append(root / record["basename"])
        elif record.get("path"):
            candidates.append(Path(record["path"]))
        source = next((p for p in candidates if p.exists()), None)
        if source is None:
            raise FileNotFoundError(f"Missing required input {record['basename']}")
        observed = sha256(source)
        if observed != record["sha256"]:
            raise ValueError(f"SHA256 mismatch for {source}: {observed} != {record['sha256']}")
        checked.append({"path": str(source), "bytes": source.stat().st_size, "sha256": observed})
    return checked


def patch_text(path, replacements):
    path = Path(path)
    text = path.read_text()
    for old, new in replacements.items():
        if old not in text:
            raise ValueError(f"Patch needle not found in {path}: {old[:80]}")
        text = text.replace(old, new)
    path.write_text(text)


def run_command(args, cwd, log_path, env=None, timeout=None):
    run_env = os.environ.copy()
    run_env.update({"OPENBLAS_NUM_THREADS": "1", "OMP_NUM_THREADS": "1", "MKL_NUM_THREADS": "1"})
    if env:
        run_env.update(env)
    log_path = Path(log_path)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("w") as log:
        proc = subprocess.run(
            args,
            cwd=str(cwd),
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            env=run_env,
            timeout=timeout,
        )
        log.write(proc.stdout)
    if proc.returncode:
        sys.stderr.write(proc.stdout[-4000:])
        raise subprocess.CalledProcessError(proc.returncode, args)
    return {"args": args, "cwd": str(cwd), "log": str(log_path), "returncode": proc.returncode}
