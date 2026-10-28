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
| [`pulsatileStenosis`](pulsatileStenosis/) | Pulsatile flow through the 75 % stenosis; TAWSS and OSI | `pimpleFoam` | Womersley inflow, cycle periodicity |
| [`arterialBifurcation`](arterialBifurcation/) | Symmetric bifurcation with Murray's-law branches (`snappyHexMesh`) | `simpleFoam` | Flow split, Murray's-law wall shear stress |
| [`curvedArtery`](curvedArtery/) | Dean flow in a curved artery (De = 113) | `simpleFoam` | Curved-pipe friction-factor correlations |
| [`sidewallAneurysm`](sidewallAneurysm/) | Saccular aneurysm on a straight artery (`snappyHexMesh`) | `simpleFoam` | Poiseuille parent flow, neck mass balance, low sac shear |
| [`fdaNozzle`](fdaNozzle/) | FDA benchmark medical-device nozzle, Re_throat = 500 | `simpleFoam` | Published centreline velocity at the expansion |
| [`graetzMassTransfer`](graetzMassTransfer/) | Oxygen transport to the wall of an arteriole (Pe = 250) | `simpleFoam` + `scalarTransport` | Exact Graetz solution, Sherwood number |
| [`bloodParticles`](bloodParticles/) | Tracer and debris particles in arterial flow | `icoUncoupledKinematicParcelFoam` | Exact particle trajectories, settling and deposition |

Every case except `porousValve` includes a `validate.py` script that `Allrun` calls after the solver; it compares the result with an analytical or physical reference and returns a non-zero exit status if a check fails. The scripts need `python3` with `numpy`; `matplotlib` is optional and only used for the figures.

Additional biofluid-analysis test cases will be added as separate top-level directories.
