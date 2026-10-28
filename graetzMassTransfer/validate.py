#!/usr/bin/env python3
"""Compare the graetzMassTransfer OpenFOAM result with the Graetz solution.

Problem: fully developed Poiseuille flow u(r) = 2U(1 - r^2/R^2) carries a
dissolved species (C = 1 at the inlet) towards a wall held at C = 0. The
classical Graetz solution neglects axial diffusion:

    u dC/dy = D_C (1/r) d/dr (r dC/dr),   C(r, 0) = 1,   C(R, y) = 0.

It is solved here by eigen-decomposition of a fine (400-point) radial
finite-difference operator. At
Pe = U D / D_C = 250, axial diffusion changes the bulk concentration by less
than 0.1 % beyond y = D.

Checks (all must pass):
  1. Bulk (flow-weighted) concentration at y/D = 1, 2, 4, 7, 10, 15 within 2 %
     of the Graetz solution (relative to the local value).
  2. Fully developed Sherwood number, from the decay of the bulk
     concentration between y = 10 D and 15 D, within 2 % of 3.657.

Only numpy is required (matplotlib is optional and used for the figure).
"""
import math
import sys
from pathlib import Path

import numpy as np

D = 1.0e-4
R = D / 2.0
U = 0.005
DC = 2.0e-9
PE = U * D / DC
SH_INF = 3.6568
STATIONS = {"y01": 1.0, "y02": 2.0, "y04": 4.0, "y07": 7.0, "y10": 10.0, "y15": 15.0}
TOLERANCE = 0.02


def graetz_bulk(y_over_d, n_r=400):
    """Bulk concentration of the parabolic Graetz problem at the given y/D.

    The radial operator is discretised with second-order finite differences
    and the resulting linear ODE system dC/dy = M C is solved exactly by
    eigen-decomposition.
    """
    r = np.linspace(0.0, R, n_r + 1)
    dr = r[1]
    u = 2.0 * U * (1.0 - (r / R) ** 2)
    n = n_r                                       # unknowns: r[0] ... r[n_r-1]
    A = np.zeros((n, n))
    A[0, 0], A[0, 1] = -4.0 / dr**2, 4.0 / dr**2
    for i in range(1, n):
        A[i, i - 1] = 1.0 / dr**2 - 1.0 / (2.0 * r[i] * dr)
        A[i, i] = -2.0 / dr**2
        if i + 1 < n:
            A[i, i + 1] = 1.0 / dr**2 + 1.0 / (2.0 * r[i] * dr)
    M = DC * A / u[:n, None]
    eigval, eigvec = np.linalg.eig(M)
    inv = np.linalg.inv(eigvec)
    c0 = np.ones(n)
    coeff = inv @ c0
    w = (u * r)[:n]
    norm = np.trapezoid(u * r, r)
    out = {}
    for t in y_over_d:
        c = np.real(eigvec @ (coeff * np.exp(eigval * t * D)))
        cf = np.append(c, 0.0)
        out[t] = np.trapezoid(u * r * cf, r) / norm
    return out


def bulk_from_profile(path):
    data = np.loadtxt(path)                      # x  C  Ux  Uy  Uz (axis to wall)
    r, c, u = data[:, 0], data[:, 1], data[:, 3]
    if r[-1] < 0.999 * R:                         # wall point not sampled
        r, c, u = np.append(r, R), np.append(c, 0.0), np.append(u, 0.0)
    return np.trapezoid(u * c * r, r) / np.trapezoid(u * r, r)


def latest_dir(path):
    best = None
    for d in path.iterdir():
        try:
            t = float(d.name)
        except ValueError:
            continue
        if best is None or t > best[0]:
            best = (t, d)
    return best[1] if best else None


def main():
    case = Path(__file__).resolve().parent
    prof = case / "postProcessing" / "profiles"
    if not prof.is_dir():
        sys.exit("No postProcessing data found: run ./Allrun first.")
    last = latest_dir(prof)

    exact = graetz_bulk(list(STATIONS.values()))
    numeric = {}
    worst = 0.0
    print(f"Peclet number Pe = U D / D_C  : {PE:.0f}")
    print(f"Iteration                     : {last.name}")
    print(" y/D   Cb OpenFOAM   Cb Graetz   error")
    for name, yd in STATIONS.items():
        cb = bulk_from_profile(last / f"{name}_C_U.xy")
        numeric[yd] = cb
        err = abs(cb - exact[yd]) / exact[yd]
        worst = max(worst, err)
        print(f" {yd:4.0f}   {cb:10.5f}   {exact[yd]:10.5f}   {100*err:5.2f} %")

    # Fully developed Sherwood number from the bulk decay between 10 D and 15 D:
    # d(ln Cb)/dy = -4 Sh / (D Pe)
    sh = PE * math.log(numeric[10.0] / numeric[15.0]) / (4.0 * 5.0)
    sh_err = abs(sh - SH_INF) / SH_INF
    print(f"Fully developed Sherwood no.  : {sh:.4f} (Graetz limit {SH_INF}, {100*sh_err:.2f} %)")

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        yy = np.linspace(0.05, 16, 200)
        ex = graetz_bulk(list(yy))
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.semilogy(yy, [ex[t] for t in yy], "k-", label="Graetz solution")
        ax.semilogy(list(numeric), list(numeric.values()), "o", label="OpenFOAM")
        ax.set_xlabel("y / D")
        ax.set_ylabel("bulk concentration C_b")
        ax.set_title(f"Mass transfer to the wall, Pe = {PE:.0f}")
        ax.legend()
        fig.tight_layout()
        fig.savefig(case / "graetz_validation.png", dpi=150)
        print("Figure written                : graetz_validation.png")
    except ImportError:
        pass

    checks = {
        f"bulk concentration (< {100*TOLERANCE:g} %)": worst < TOLERANCE,
        f"fully developed Sherwood number (< {100*TOLERANCE:g} %)": sh_err < TOLERANCE,
    }
    for name, ok in checks.items():
        print(f"  [{'PASS' if ok else 'FAIL'}] {name}")
    ok = all(checks.values())
    print("RESULT                        :", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
