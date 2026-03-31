# reduce2 reference environment

cctbx reduce2 (hydrogen placement) for generating reference data.

## Setup

```bash
cd scripts/reduce2

# Install cctbx-base via pixi
pixi install

# Copy chem_data from a cctbx-base environment with chem_data installed.
# If you previously had a global pixi env:
#   cp -R ~/.pixi/envs/cctbx-base/site-packages/chem_data \
#     .pixi/envs/default/lib/python3.*/site-packages/chem_data
#
# Or download from cctbx GitHub releases and extract into the same path.
```

## Usage

```bash
# Single file
pixi run python run_reduce2.py input.cif output.cif

# Batch (all 11 structures)
pixi run python run_reduce2_batch.py
```
