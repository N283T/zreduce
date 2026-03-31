#!/usr/bin/env python3
"""Run cctbx reduce2 on all 11 test structures.

Usage:
    cd scripts/reduce2
    pixi run python run_reduce2_batch.py
"""

import os
import subprocess
import sys
import time

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.join(SCRIPT_DIR, "..", "..")
DATA_DIR = os.path.join(PROJECT_ROOT, "examples", "data")
OUTPUT_DIR = os.path.join(PROJECT_ROOT, "examples", "result", "reduce2")

STRUCTURES = [
    "AF-C1P619-F1-model_v6.cif",
    "AF-P0A9J6-F1-model_v6.cif",
    "AF-P0DSH8-F1-model_v6.cif",
    "AF-P22523-F1-model_v6.cif",
    "AF-P76347-F1-model_v6.cif",
    "1rqf.cif.gz",
    "2cf8.cif.gz",
    "2hnt.cif.gz",
    "3rk2.cif.gz",
    "6fys.cif.gz",
    "fold_test2_model_0.cif",
]


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    runner = os.path.join(SCRIPT_DIR, "run_reduce2.py")

    results = []
    for name in STRUCTURES:
        stem = name.replace(".cif.gz", "").replace(".cif", "")
        input_path = os.path.join(DATA_DIR, name)
        output_path = os.path.join(OUTPUT_DIR, f"{stem}_reduce2.cif")

        if os.path.exists(output_path):
            print(f"  SKIP {stem} (already exists)")
            results.append((stem, "skip", 0))
            continue

        print(f"  RUN  {stem} ...", end="", flush=True)
        t0 = time.time()
        proc = subprocess.run(
            [sys.executable, runner, input_path, output_path],
            capture_output=True,
            text=True,
        )
        elapsed = time.time() - t0

        if proc.returncode == 0:
            print(f" OK ({elapsed:.1f}s) {proc.stdout.strip()}")
            results.append((stem, "ok", elapsed))
        else:
            print(f" FAIL ({elapsed:.1f}s)")
            print(f"    stderr: {proc.stderr[:200]}")
            results.append((stem, "fail", elapsed))

    print(f"\n{'=' * 50}")
    print(f"{'Structure':40s} {'Status':6s} {'Time':>6s}")
    print(f"{'-' * 50}")
    for stem, status, elapsed in results:
        print(f"{stem:40s} {status:6s} {elapsed:5.1f}s")


if __name__ == "__main__":
    main()
