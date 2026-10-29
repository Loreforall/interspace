# Reproduction

Python 3.10+; install `python -m pip install -e ".[test,plots]"`. Run `python -m pytest -q`. For a fresh engineering figure and summary: `python -m examples.run_stage12_release --output release/generated`. Historical stage experiment runners remain under `examples/`. All examples use illustrative inputs; do not interpret the dimensional values as a calibrated design. `release/generated/` includes machine-readable CSV/JSON and PNG. Matplotlib output uses explicit axes, grids and numerical diagnostics.
