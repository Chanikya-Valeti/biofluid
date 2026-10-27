# Biofluid analysis test cases

This repository collects reproducible computational test cases for biofluid analysis. Each case is stored in its own self-contained directory and includes the input files, run instructions, and supporting documentation needed to reproduce the simulation.

## Available test cases

All cases use OpenFOAM.com v2412.

| Case | Biofluid problem | Solver | Reference / check |
| --- | --- | --- | --- |
| [`porousValve`](porousValve/) | Valve opening and closing in a cylindrical tube, represented by time-dependent porous resistance regions | `pimpleFoam` | — |
| [`womersleyPipe`](womersleyPipe/) | Pulsatile flow in a rigid artery (α = 5.5) | `pimpleFoam` | Exact Womersley solution |
| [`carreauBloodPipe`](carreauBloodPipe/) | Shear-thinning (Carreau) blood flow in a small artery | `simpleFoam` | Exact non-Newtonian tube-flow solution |
| [`stenosedArtery`](stenosedArtery/) | Flow through a 75 % cosine stenosis, Newtonian or Carreau blood | `simpleFoam` | Poiseuille inflow/recovery, recirculation |
| [`arterialBifurcation`](arterialBifurcation/) | Symmetric bifurcation with Murray's-law branches (`snappyHexMesh`) | `simpleFoam` | Flow split, Murray's-law wall shear stress |

The four newer cases each include a `validate.py` script that `Allrun` calls after the solver; it compares the result with an analytical or physical reference and returns a non-zero exit status if a check fails. The scripts need `python3` with `numpy`; `matplotlib` is optional and only used for the figures.

Additional biofluid-analysis test cases will be added as separate top-level directories.
