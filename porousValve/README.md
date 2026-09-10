# Porous valve OpenFOAM case

An incompressible-flow OpenFOAM testcase that represents a valve as a porous disk in a cylindrical tube. The disk is divided into a core and three annular zones. Time-limited Darcy resistance sources close and reopen the zones in sequence; leaflet mechanics are not resolved.

## Requirements

- OpenFOAM.com v2412
- `blockMesh`, `topoSet`, `checkMesh`, and `pimpleFoam`
- MPI only for the parallel workflow

Load the OpenFOAM v2412 environment before running the case. The case was checked with the v2412 executables on macOS.

## Case layout

```text
0/          Initial velocity and pressure
constant/   Fluid, turbulence, and porous-source dictionaries
system/     Mesh, zone, solver, and decomposition dictionaries
Allrun      Reproducible serial or parallel workflow
Allclean    Removes regenerated output only
```

The mesh is generated from `system/blockMeshDict`; `system/topoSetDict` then creates the five cell zones needed by `constant/fvOptions`. No `constant/polyMesh` is stored in the repository.

## Model

The 0.25 m long cylindrical domain has a radius of 0.02 m. `INLET` and `OUTLET` use zero-gradient velocity and fixed pressure values of 1 and 0 m²/s² respectively; `WALL` is no-slip. The fluid is Newtonian and laminar with kinematic viscosity `3.3e-06 m²/s`.

`pimpleFoam` advances the case from 0 to 1 s with a fixed `0.001 s` time step. The baseline porous region has Darcy coefficient `1e6 m⁻²`; closing sources add `1e12 m⁻²`. The Forchheimer coefficient is zero.

| Zone | Radius range (m) | Closed interval (s) |
| --- | ---: | --- |
| `porousCore` | 0–0.005 | 0.50–0.60 |
| `porousRing1` | 0.005–0.010 | 0.35–0.70 |
| `porousRing2` | 0.010–0.015 | 0.25–0.80 |
| `porousRing3` | 0.015–0.020 | 0.15–0.90 |

The porous disk extends from `y=-0.01 m` to `y=0.01 m` (20 mm thickness). The specified resistance changes abruptly at the listed times.

## Run

```bash
./Allrun
./Allrun parallel
```

The parallel run reads the processor count (`40`) from `system/decomposeParDict`, decomposes the case, runs `pimpleFoam -parallel`, and reconstructs the latest time.

To remove only generated results, logs, processor directories, and regenerated mesh data:

```bash
./Allclean
```

## Post-processing

Open the case in ParaView with `paraFoam`. The `Images/` directory contains the supplied geometry figures and a demonstration video:

- `Images/Cylindrical_Tube.jpeg` — cylindrical tube geometry
- `Images/Porous_block.jpeg` — porous-region illustration
- `Images/movie.mp4` — demonstration output

The previously supplied ParaView state file is not part of this repository.

## Validation

In a clean temporary copy with OpenFOAM.com v2412, `blockMesh`, `topoSet`, and `checkMesh -allTopology -allGeometry` completed successfully. The mesh contains 50,000 hexahedral cells, three patches, and five porous cell zones. A one-time-step `pimpleFoam` smoke test completed successfully.

## Limitations

This is a numerical porous-resistance representation of a valve. It does not resolve leaflets, fluid-structure interaction, or wall motion. The resistance coefficients and timing should be validated for the intended application.

## Author, citation, and licence

Author: Chanikya Valeti.

If you use this testcase, please cite this repository and the OpenFOAM version used. No DOI is currently assigned.

No reuse licence is currently granted. Copyright © Chanikya Valeti; all rights reserved.
