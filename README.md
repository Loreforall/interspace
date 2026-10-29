# INTERSPACE — Vercel engineering workstation

The deployable site is in `site/`. It includes the Stage 12 Python archive, the recorded simulation, equations, original reports, CSV numerical data, and browsable solver source.

## Vercel

1. Push this folder to GitHub.
2. Import the repository into Vercel.
3. Framework Preset: **Other**. Root Directory: **site**. Build Command: **none**. Output Directory: **.** (site root).
4. Deploy, then add your domain in Vercel Project Settings > Domains.

Alternatively, when using the repository root as Vercel Root Directory, set Output Directory to `site` and leave Build Command empty.

## Local preview

`cd site && python -m http.server 8080`

## Interpretation

The visualizer replays Stage 12 baseline results; it is not a browser-side physics solver. Source code and numerical datasets are provided for independent review. Model consistency is not evidence that the speculative free-standing barrier is feasible.
