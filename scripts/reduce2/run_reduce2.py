#!/usr/bin/env python3
"""Run cctbx reduce2 on an mmCIF file.

Requires:
  - cctbx-base installed via pixi (scripts/reduce2/pixi.toml)
  - chem_data symlinked into pixi env's site-packages
    (from ~/.pixi/envs/cctbx-base/site-packages/chem_data)

Usage:
    cd scripts/reduce2
    pixi run python run_reduce2.py ../../examples/data/input.cif output.cif
"""

import os
import sys


def _setup_chem_data():
    """Find chem_data in site-packages and configure cctbx paths."""
    import site

    for sp in site.getsitepackages():
        chem_data = os.path.join(sp, "chem_data")
        if os.path.isdir(chem_data):
            break
    else:
        raise RuntimeError(
            "chem_data not found in site-packages. "
            "Symlink it from ~/.pixi/envs/cctbx-base/site-packages/chem_data"
        )

    os.environ["MMTBX_CCP4_MONOMER_LIB"] = os.path.join(chem_data, "geostd")

    import libtbx.load_env
    from libtbx.path import absolute_path

    libtbx.env.repository_paths.insert(0, absolute_path(chem_data))

    import mmtbx.chemical_components as cc

    cc.data_dir = cc.find_data_dir()
    assert cc.data_dir is not None, f"chemical_components not found in {chem_data}"


def main():
    if len(sys.argv) < 3:
        print(f"Usage: {sys.argv[0]} input.cif output.cif [extra_args...]")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    _setup_chem_data()

    from iotbx.cli_parser import run_program
    from mmtbx.programs.reduce2 import Program

    run_program(
        program_class=Program,
        args=[
            input_file,
            f"output.filename={output_file}",
            "output.overwrite=True",
            "use_neutron_distances=True",
            "add_flip_movers=True",
        ]
        + sys.argv[3:],
        logger=open(os.devnull, "w"),
    )

    import iotbx.pdb

    pdb_out = iotbx.pdb.input(output_file)
    h_count = sum(
        1
        for a in pdb_out.construct_hierarchy().atoms()
        if a.element.strip() in ("H", "D")
    )
    print(f"Output: {output_file} ({h_count} H atoms)")


if __name__ == "__main__":
    main()
