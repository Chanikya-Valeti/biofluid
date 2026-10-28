#!/usr/bin/env python3
"""Generate constant/triSurface/aneurysm.stl: an idealised sidewall saccular
aneurysm on a straight parent artery.

Geometry (flow along +y; the aneurysm bulges towards +z):
    parent artery  radius R = 2 mm (D = 4 mm, internal-carotid size),
                   from y = -8 D to y = +12 D
    aneurysm sac   sphere of radius Ra = 3 mm (6 mm dome diameter) centred at
                   y = 0, z = c, with c chosen so that the neck (the sphere's
                   intersection with the artery wall in the symmetry plane)
                   is 4 mm long:  c = R + sqrt(Ra^2 - (Wn/2)^2)
    dome height (neck plane to dome top) = c + Ra - R = 5.24 mm,
    aspect ratio (height / neck width) = 1.31

The STL has the named solids INLET, OUTLET and WALL.

Requires: numpy, trimesh, manifold3d  (pip install numpy trimesh manifold3d)
The generated STL is stored in the repository, so this script only needs to
be run when the geometry is changed.
"""
import math
from pathlib import Path

import numpy as np
import trimesh

R = 0.002
D = 2 * R
RA = 0.003
NECK = 0.004
C = R + math.sqrt(RA**2 - (NECK / 2) ** 2)
Y_IN = -8 * D
Y_OUT = 12 * D
SEGMENTS = 96


def main():
    length = Y_OUT - Y_IN
    artery = trimesh.creation.cylinder(radius=R, height=length, sections=SEGMENTS)
    artery.apply_transform(trimesh.geometry.align_vectors([0, 0, 1], [0, 1, 0]))
    artery.apply_translation([0, (Y_IN + Y_OUT) / 2, 0])
    sac = trimesh.creation.icosphere(subdivisions=4, radius=RA)
    sac.apply_translation([0, 0, C])

    shape = trimesh.boolean.union([artery, sac], engine="manifold")
    assert shape.is_watertight, "union is not watertight"

    centres, normals = shape.triangles_center, shape.face_normals
    labels = np.full(len(centres), "WALL", dtype=object)
    for name, y, n in (("INLET", Y_IN, -1.0), ("OUTLET", Y_OUT, 1.0)):
        labels[(np.abs(centres[:, 1] - y) < 1e-7) & (normals[:, 1] * n > 0.999)] = name

    out = Path(__file__).resolve().parent.parent / "constant" / "triSurface" / "aneurysm.stl"
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w") as f:
        for name in ("INLET", "OUTLET", "WALL"):
            idx = np.where(labels == name)[0]
            f.write(f"solid {name}\n")
            for i in idx:
                nx, ny, nz = normals[i]
                f.write(f"  facet normal {nx:.6e} {ny:.6e} {nz:.6e}\n    outer loop\n")
                for v in shape.triangles[i]:
                    f.write(f"      vertex {v[0]:.8e} {v[1]:.8e} {v[2]:.8e}\n")
                f.write("    endloop\n  endfacet\n")
            f.write(f"endsolid {name}\n")
            print(f"{name:7s}: {len(idx)} triangles")
    print(f"Wrote {out}")
    print(f"sac centre z = {C*1e3:.3f} mm, dome height = {(C+RA-R)*1e3:.3f} mm, "
          f"aspect ratio = {(C+RA-R)/NECK:.2f}")


if __name__ == "__main__":
    main()
