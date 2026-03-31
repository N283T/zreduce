#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "gemmi",
# ]
# ///
"""Split wwPDB components.cif.gz into individual component files for cctbx.

cctbx's chemical_components module expects files at:
  {data_dir}/{first_char_lower}/data_{COMP_ID}.cif

Usage:
    ./scripts/reduce2/split_ccd.py [components.cif.gz] [output_dir]
"""

import os
import sys

import gemmi


def main():
    ccd_path = sys.argv[1] if len(sys.argv) > 1 else os.path.expanduser(
        "~/pdb/data/monomers/components.cif.gz"
    )
    output_dir = sys.argv[2] if len(sys.argv) > 2 else os.path.expanduser(
        "~/pdb/data/monomers/components_split"
    )

    print(f"Reading {ccd_path}...")
    doc = gemmi.cif.read(ccd_path)
    print(f"Found {len(doc)} blocks")

    count = 0
    for block in doc:
        comp_id = block.name
        if not comp_id:
            continue

        # Create subdirectory: first char lowercase
        subdir = os.path.join(output_dir, comp_id[0].lower())
        os.makedirs(subdir, exist_ok=True)

        # Write individual CIF file
        out_path = os.path.join(subdir, f"data_{comp_id}.cif")
        out_doc = gemmi.cif.Document()
        out_doc.add_copied_block(block)
        out_doc.write_file(out_path)
        count += 1

        if count % 5000 == 0:
            print(f"  {count} components written...")

    print(f"Done: {count} components written to {output_dir}")


if __name__ == "__main__":
    main()
