#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.12"
# dependencies = []
# ///
"""Run ChimeraX addh on test structures to generate reference data.

Executes ChimeraX in nogui mode for each structure in examples/data/,
saves output to examples/result/chimerax_addh/.

Usage:
    ./scripts/chimerax_addh.py                        # all structures
    ./scripts/chimerax_addh.py --structures AF-P0A9J6  # filter by name
    ./scripts/chimerax_addh.py --template              # enable template option
    ./scripts/chimerax_addh.py --dry-run               # show commands without running

Options:
    --template      Use idealized coordinates for atom typing in non-standard
                    residues (ChimeraX addh template true). Default: false.
    --chimerax PATH Path to ChimeraX binary. Default: chimerax (from PATH).
"""

import argparse
import gzip
import os
import shutil
import subprocess
import sys
import tempfile

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.join(PROJECT_ROOT, "examples", "data")
OUTPUT_DIR = os.path.join(PROJECT_ROOT, "examples", "result", "chimerax_addh")

# Common ChimeraX install locations (macOS)
_CHIMERAX_CANDIDATES = [
    shutil.which("chimerax"),
    shutil.which("ChimeraX"),
]
# Search /Applications for versioned installs
try:
    chimerax_apps = [
        e
        for e in os.listdir("/Applications")
        if e.startswith("ChimeraX") and e.endswith(".app")
    ]

    # Sort by version number (extract digits after "ChimeraX-")
    def _version_key(name: str) -> list[int]:
        import re

        m = re.search(r"(\d+(?:\.\d+)*)", name)
        return [int(x) for x in m.group(1).split(".")] if m else [0]

    for entry in sorted(chimerax_apps, key=_version_key, reverse=True):
        _CHIMERAX_CANDIDATES.append(
            os.path.join("/Applications", entry, "Contents", "bin", "ChimeraX")
        )
except OSError:
    pass

DEFAULT_CHIMERAX = next(
    (c for c in _CHIMERAX_CANDIDATES if c and os.path.isfile(c)),
    "chimerax",
)

STRUCTURES = [
    "AF-C1P619-F1-model_v6",
    "AF-P0A9J6-F1-model_v6",
    "AF-P0DSH8-F1-model_v6",
    "AF-P22523-F1-model_v6",
    "AF-P76347-F1-model_v6",
    "1rqf",
    "2cf8",
    "2hnt",
    "3rk2",
    "6fys",
    "fold_test2_model_0",
]


def find_input(stem: str) -> str | None:
    """Find input CIF file (plain or gzipped)."""
    for ext in (".cif", ".cif.gz"):
        path = os.path.join(DATA_DIR, f"{stem}{ext}")
        if os.path.exists(path):
            return path
    return None


def run_addh(
    stem: str,
    input_path: str,
    *,
    template: bool,
    chimerax_bin: str,
    dry_run: bool,
) -> bool:
    """Run ChimeraX addh on a single structure. Returns True on success."""
    output_path = os.path.join(OUTPUT_DIR, f"{stem}_addh.cif")
    output_gz = f"{output_path}.gz"

    # ChimeraX nogui cannot open .cif.gz directly — decompress to temp file
    needs_decompress = input_path.endswith(".gz")

    try:
        if needs_decompress:
            tmp_fd, tmp_path = tempfile.mkstemp(suffix=".cif")
            os.close(tmp_fd)
            with gzip.open(input_path, "rb") as f_in:
                with open(tmp_path, "wb") as f_out:
                    shutil.copyfileobj(f_in, f_out)
            open_path = tmp_path
        else:
            tmp_path = None
            open_path = input_path

        # Build ChimeraX command sequence
        template_flag = "true" if template else "false"
        cxc_commands = [
            f"open {open_path}",
            f"addh template {template_flag}",
            f"save {output_path} format mmcif",
            "exit",
        ]
        cxc_script = " ; ".join(cxc_commands)

        cmd = [chimerax_bin, "--nogui", "--cmd", cxc_script]

        if dry_run:
            print(f"  [dry-run] {' '.join(cmd)}")
            return True

        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=600,
        )
        if result.returncode != 0:
            print(f"  ERROR (exit {result.returncode})")
            if result.stderr:
                for line in result.stderr.strip().splitlines()[-5:]:
                    print(f"    {line}")
            return False

        # Gzip the output
        if os.path.exists(output_path):
            with open(output_path, "rb") as f_in:
                with gzip.open(output_gz, "wb") as f_out:
                    shutil.copyfileobj(f_in, f_out)
            os.remove(output_path)
            print(f"  -> {output_gz}")
            return True
        else:
            print(f"  ERROR: output file not created")
            return False

    except subprocess.TimeoutExpired:
        print(f"  ERROR: timeout (600s)")
        return False
    except FileNotFoundError:
        print(f"  ERROR: ChimeraX not found at '{chimerax_bin}'")
        sys.exit(1)
    finally:
        if tmp_path and os.path.exists(tmp_path):
            os.remove(tmp_path)


def main():
    parser = argparse.ArgumentParser(description="Run ChimeraX addh on test structures")
    parser.add_argument(
        "--structures",
        nargs="+",
        help="Filter structures by substring match",
    )
    parser.add_argument(
        "--template",
        action="store_true",
        help="Enable template option for non-standard residue atom typing (default: false)",
    )
    parser.add_argument(
        "--chimerax",
        default=DEFAULT_CHIMERAX,
        help=f"Path to ChimeraX binary (default: {DEFAULT_CHIMERAX})",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Show commands without running",
    )
    args = parser.parse_args()

    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # Filter structures
    targets = STRUCTURES
    if args.structures:
        targets = [s for s in STRUCTURES if any(f in s for f in args.structures)]
        if not targets:
            print(f"No structures matched: {args.structures}")
            sys.exit(1)

    template_str = "true" if args.template else "false"
    print(f"ChimeraX addh (template={template_str})")
    print(f"  binary: {args.chimerax}")
    print(f"  output: {OUTPUT_DIR}")
    print(f"  structures: {len(targets)}")
    print()

    ok = 0
    fail = 0
    for stem in targets:
        input_path = find_input(stem)
        if input_path is None:
            print(f"  {stem}: input not found, skipping")
            fail += 1
            continue

        print(f"  {stem}...", end=" ", flush=True)
        if run_addh(
            stem,
            input_path,
            template=args.template,
            chimerax_bin=args.chimerax,
            dry_run=args.dry_run,
        ):
            ok += 1
        else:
            fail += 1

    print(f"\nDone: {ok} succeeded, {fail} failed")


if __name__ == "__main__":
    main()
