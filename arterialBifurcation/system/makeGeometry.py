#!/usr/bin/env python3
"""Generate constant/triSurface/bifurcation.stl: an idealised, symmetric,
planar arterial bifurcation with Murray's-law daughter radii.

Geometry (all in the y-z plane; flow enters along +y):
    parent   radius Rp = 3 mm, from y = -6 Dp up to the branch point (origin)
    daughters radius Rd = Rp / 2^(1/3) = 2.381 mm (Murray's law),
             +/- 35 degrees from the y-axis (70 degree bifurcation angle),
             length 10 Dd measured from the branch point
A sphere of radius Rp at the branch point blends the outer junction walls.

The STL has four named solids: INLET, OUTLET1 (+z branch), OUTLET2 (-z
branch) and WALL, which become the mesh patches.

Requires: numpy, trimesh, manifold3d  (pip install numpy trimesh manifold3d)
The generated STL is stored in the repository, so this script only needs to
be run when the geometry is changed.
"""
import math
from pathlib import Path

import numpy as np
import trimesh

RP = 0.003
DP = 2 * RP
RD = RP / 2 ** (1.0 / 3.0)
DD = 2 * RD
HALF_ANGLE = math.radians(35.0)
L_PARENT = 6.0 * DP
L_DAUGHTER = 10.0 * DD
SEGMENTS = 96


def cylinder(radius, length, axis, start):
    """Closed cylinder from `start` along unit vector `axis`."""
    cyl = trimesh.creation.cylinder(radius=radius, height=length, sections=SEGMENTS)
    # trimesh cylinders are centred on the origin along z
    cyl.apply_translation([0, 0, length / 2.0])
    rot = trimesh.geometry.align_vectors([0, 0, 1], axis)
    cyl.apply_transform(rot)
    cyl.apply_translation(start)
    return cyl


def main():
    y_axis = np.array([0.0, 1.0, 0.0])
    d1 = np.array([0.0, math.cos(HALF_ANGLE), math.sin(HALF_ANGLE)])
    d2 = np.array([0.0, math.cos(HALF_ANGLE), -math.sin(HALF_ANGLE)])

    parent = cylinder(RP, L_PARENT + RP, y_axis, [0, -L_PARENT, 0])
    branch1 = cylinder(RD, L_DAUGHTER, d1, [0, 0, 0])
    branch2 = cylinder(RD, L_DAUGHTER, d2, [0, 0, 0])
    blend = trimesh.creation.icosphere(subdivisions=5, radius=RP)

    shape = trimesh.boolean.union([parent, branch1, branch2, blend], engine="manifold")
    assert shape.is_watertight, "union is not watertight"

    # Label faces: caps are flat discs normal to the vessel axis at the ends
    centres = shape.triangles_center
    normals = shape.face_normals
    labels = np.full(len(centres), "WALL", dtype=object)
    caps = {
        "INLET": (np.array([0, -L_PARENT, 0]), -y_axis),
        "OUTLET1": (L_DAUGHTER * d1, d1),
        "OUTLET2": (L_DAUGHTER * d2, d2),
    }
    for name, (origin, n) in caps.items():
        on_plane = np.abs((centres - origin) @ n) < 1e-7
        aligned = normals @ n > 0.999
        labels[on_plane & aligned] = name

    out = Path(__file__).resolve().parent.parent / "constant" / "triSurface" / "bifurcation.stl"
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w") as f:
        for name in ("INLET", "OUTLET1", "OUTLET2", "WALL"):
            idx = np.where(labels == name)[0]
            f.write(f"solid {name}\n")
            for i in idx:
                nx, ny, nz = normals[i]
                f.write(f"  facet normal {nx:.6e} {ny:.6e} {nz:.6e}\n    outer loop\n")
                for v in shape.triangles[i]:
                    f.write(f"      vertex {v[0]:.8e} {v[1]:.8e} {v[2]:.8e}\n")
                f.write("    endloop\n  endfacet\n")
            f.write(f"endsolid {name}\n")
            print(f"{name:8s}: {len(idx)} triangles")
    print(f"Wrote {out}")
    print(f"Rp = {RP*1e3:.3f} mm, Rd = {RD*1e3:.3f} mm, volume = {shape.volume*1e9:.1f} mm^3")


if __name__ == "__main__":
    main()
