# INTERSPACE: Stage 10 extension and Stage 11 initial verification

## Scope delivered
- Ring-band kinetic energy history and first threshold-crossing times (a *proxy* for wave arrival; the fixed ring itself has no kinetic energy).
- Ideal downstream electrical storage with capacity, efficiency, conversion losses and overflow.
- Repeatable 3D off-axis impact sweeps across contact position, piezo coupling and time step.
- Energy-partition assertions, validation of input parameters, and time-step comparison output.

## Limitations / incomplete acceptance gates
- Impact is still **nearest-node vertical-only contact**, not full surface-to-surface 3D contact.
- Energy storage is a downstream accounting model, not a physically integrated rectifier/charging circuit; storage does **not** affect transducer forces.
- Ring is fixed; wave arrival at ring is estimated from the last *free* latitude band, not from an experimentally calibrated flux measurement.
- No calibrated material constants, mesh-refinement convergence, experimental comparison, uncertainty propagation, or matched flat-vs-dome comparison.
- Stage 9 prescribed deployment violates its illustrative actuator constraints and has **not** been repaired.
- Therefore Stage 10 is an *extended prototype*, not physically complete, and Stage 11 is a *preliminary numerical study*, not validated optimization.

## Run
```bash
pip install -e '.[test]'
pytest -q
python -m examples.run_stage11
```

## Scientific interpretation
A mechanical dome with transducers can exchange energy with a load. This does not demonstrate energy filling vacuum gaps between atoms, a free-standing field barrier, or physically achievable rapid deployment. Resistive output is an electrical-energy proxy; the storage ledger is idealized.
