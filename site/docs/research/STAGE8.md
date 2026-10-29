# Stage 8 — 3D ring-supported hemispherical reference model

## Physical interpretation

This is a *material* dome represented by a three-dimensional triangular mass-spring-damper network. Its equatorial ring is rigidly supported. A centered spherical projectile impacts the dome pole. Deformation waves carry energy toward the ring; relative motion produces heat. It is **not** a free-standing field, atomic-gap storage mechanism, deployed dome, calibrated material, or electric energy harvester.

## Scope and limitations

- SI units; geometrical 3D mesh, centered normal impact, fixed ring, explicit projectile contact.
- Energy ledger: projectile kinetic, network kinetic, spring elastic, contact elastic, heat; total numerical residual.
- No harvesting efficiency is inferred from heat. The fixed ring transfers reaction forces to its supports; its idealized stationary boundary does no work.
- Semi-implicit Euler is a starter integrator; its residual must be measured and refined. Default network spring parameters are illustrative, **not material-calibrated**.
- Next upgrades: deployable membrane geometry and activation work, off-center contact and finite projectile contact patch, mesh convergence, continuum shell calibration, ring electromechanical transducers, thermal response, experimental comparison.

## Run

```bash
python -m pip install -e '.[test,plots]'
pytest -q
python -m examples.run_hemisphere
```

Outputs in `results/stage8/`.
