# Stage 9 — Deployable dome: initial inverse-dynamics reference

The stowed network is nearly flat; a quintic deployment schedule lifts the free nodes into a hemisphere. The ring remains fixed. This implementation computes **ideal actuator forces** from prescribed acceleration and internal spring/damper forces. It measures instantaneous mechanical power, positive and negative work, and numerical energy-balance residuals. Power and force budgets are reported as pass/fail checks, **not physically enforced**: if they fail, the prescribed trajectory is infeasible under those budgets.

This is a first deployment model, not a finished mechanical design. It does not include actuator placement, buckling, self-collision, material calibration, membrane bending, motor electrical efficiency, impact during deployment, or the energy-harvesting ring. Those belong to the remaining Stage 9 refinement and Stages 10–12.

Run: `python -m examples.run_deployment` and `pytest -q`.
