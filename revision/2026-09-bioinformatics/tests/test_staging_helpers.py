#!/usr/bin/env python3
"""Small portability checks for staging helpers."""
import hashlib
import json
import tempfile
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src" / "repro_package"))
from staging import verify_hash_records


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def test_search_roots_override_historical_path():
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        historical = root / "historical"
        configured = root / "configured"
        historical.mkdir()
        configured.mkdir()
        (historical / "input.txt").write_text("wrong")
        wanted = configured / "input.txt"
        wanted.write_text("right")
        record = {"path": str(historical / "input.txt"), "basename": "input.txt", "sha256": digest(wanted)}
        checked = verify_hash_records([record], search_roots=[configured])
        assert checked[0]["path"] == str(wanted)


def test_missing_configured_root_file_fails_even_if_historical_exists():
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        historical = root / "historical"
        configured = root / "configured"
        historical.mkdir()
        configured.mkdir()
        old = historical / "input.txt"
        old.write_text("right")
        record = {"path": str(old), "basename": "input.txt", "sha256": digest(old)}
        try:
            verify_hash_records([record], search_roots=[configured])
        except FileNotFoundError:
            return
        raise AssertionError("configured search_roots should not fall back to historical record.path")


def main():
    tests = [test_search_roots_override_historical_path, test_missing_configured_root_file_fails_even_if_historical_exists]
    for test in tests:
        test()
    print(json.dumps({"status": "PASS", "staging_helper_tests": len(tests)}, indent=2))


if __name__ == "__main__":
    main()

