#!/usr/bin/env python3
"""Write the exact Womersley inflow for the timeVaryingMappedFixedValue INLET.

The flow rate is Q(t) = Qm [1 + 0.6 sin(2 pi t)]; the fully developed
profile in a rigid tube of radius R0 is the steady Poiseuille profile for
Qm plus the Womersley profile for the oscillating part:

    u(r, t) = 2 Um (1 - r^2/R0^2)
            + Re{ G/(i w) [1 - J0(L r/R0) / J0(L)] e^{i w t} },  L = i^{3/2} alpha

with G chosen so that the oscillating part carries 0.6 Qm sin(w t).

Output: constant/boundaryData/INLET/points (a polar point cloud on the inlet
plane) and constant/boundaryData/INLET/<t>/U for 0 <= t <= T_END every DT,
each followed by the exact mean velocity for setAverage.
Only the Python standard library is needed. Allrun runs this script.
"""
import cmath
import math
from pathlib import Path

D = 0.006
R0 = D / 2
Y_IN = -7.0 * D            # inlet plane (see makeBlockMeshDict.py)
NU = 3.3e-6
QM = 1.5551e-6
Q1 = 0.6 * QM
FREQ = 1.0
T_END = 3.0
DT = 0.005
N_R, N_TH = 30, 64         # point cloud resolution

W = 2 * math.pi * FREQ
ALPHA = R0 * math.sqrt(W / NU)
LAM = (1j) ** 1.5 * ALPHA


def j(n, z, terms=80):
    total, term = 0j, (z / 2) ** n / math.factorial(n)
    for k in range(terms):
        total += term
        term *= -(z / 2) ** 2 / ((k + 1) * (k + 1 + n))
    return total


F = 1 - 2 * j(1, LAM) / (LAM * j(0, LAM))
G_HAT = (-1j * Q1) * 1j * W / (math.pi * R0**2 * F)     # Q1 sin(wt) = Re(-i Q1 e^{iwt})
UM = QM / (math.pi * R0**2)
J0L = j(0, LAM)


def main():
    case = Path(__file__).resolve().parent.parent
    base = case / "constant" / "boundaryData" / "INLET"
    base.mkdir(parents=True, exist_ok=True)

    radii = [R0 * i / N_R for i in range(N_R + 1)]
    pts = [(0.0, 0.0)]
    for r in radii[1:]:
        for k in range(N_TH):
            th = 2 * math.pi * k / N_TH
            pts.append((r, th))
    with open(base / "points", "w") as f:
        f.write(f"{len(pts)}\n(\n")
        for r, th in pts:
            f.write(f"({r*math.cos(th):.9g} {Y_IN:.9g} {r*math.sin(th):.9g})\n")
        f.write(")\n")

    # Radial shape functions (time independent)
    shape = {r: 1 - j(0, LAM * r / R0) / J0L for r in radii}
    n_t = int(round(T_END / DT))
    for it in range(n_t + 1):
        t = it * DT
        e = cmath.exp(1j * W * t)
        u_r = {r: 2 * UM * (1 - (r / R0) ** 2) + (G_HAT / (1j * W) * shape[r] * e).real
               for r in radii}
        tdir = base / f"{t:g}"
        tdir.mkdir(exist_ok=True)
        with open(tdir / "U", "w") as f:
            f.write(f"{len(pts)}\n(\n")
            for r, _ in pts:
                f.write(f"(0 {u_r[r]:.9g} 0)\n")
            f.write(")\n")
            # Area-averaged value enforced by setAverage: the exact mean
            # velocity Q(t) / (pi R0^2). This removes the small flow-rate error
            # of interpolating the profile onto the polygonal inlet faces.
            f.write(f"(0 {(QM + Q1 * math.sin(W * t)) / (math.pi * R0**2):.9g} 0)\n")
    print(f"Wrote {len(pts)} points and {n_t + 1} time levels to {base}")
    print(f"alpha = {ALPHA:.3f}, mean velocity {UM:.5f} m/s")


if __name__ == "__main__":
    main()
