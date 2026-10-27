#!/usr/bin/env python3
"""Compare the carreauBloodPipe OpenFOAM result with the exact fully developed
solution for a Carreau(-Yasuda) fluid in a straight tube.

For steady, fully developed flow driven by the kinematic pressure gradient G,
the shear stress varies linearly with radius:

    tau(r) / rho = G r / 2 = nu(gammaDot) gammaDot

nu(gammaDot) gammaDot increases monotonically with gammaDot for n > 0, so the
local shear rate is found by bisection and the velocity follows from

    u(r) = integral_r^R gammaDot(s) ds.

Only numpy is required (matplotlib is optional and used for the figure).
Exit status is 0 when every error is below the tolerance, 1 otherwise.
"""
import math
import sys
from pathlib import Path

import numpy as np

# Case parameters (keep in sync with 0/p, constant/transportProperties and
# system/makeBlockMeshDict.py)
R = 0.002                 # tube radius [m]
G = 0.02 / 0.01           # kinematic pressure gradient [m/s^2]
NU0 = 5.2830e-05          # zero-shear kinematic viscosity [m^2/s]
NUINF = 3.2547e-06        # infinite-shear kinematic viscosity [m^2/s]
K = 3.313                 # time constant lambda [s]
N = 0.3568                # power-law index
A = 2.0                   # Yasuda exponent (2 = Carreau)
TOLERANCE = 0.02          # 2 % of the peak value


def nu_carreau(sr):
    return NUINF + (NU0 - NUINF) * (1.0 + (K * sr) ** A) ** ((N - 1.0) / A)


def shear_rate(r):
    """Shear rate magnitude at radius r (vectorised bisection)."""
    tau = G * np.asarray(r, dtype=float) / 2.0
    lo = np.zeros_like(tau)
    hi = tau / NUINF + 1.0           # nu >= nuInf gives an upper bound
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        too_big = nu_carreau(mid) * mid > tau
        hi = np.where(too_big, mid, hi)
        lo = np.where(too_big, lo, mid)
    return 0.5 * (lo + hi)


def exact_profile(n_points=4001):
    r = np.linspace(0.0, R, n_points)
    sr = shear_rate(r)
    # u(r) = int_r^R sr ds, trapezoidal rule from the wall inwards
    seg = 0.5 * (sr[1:] + sr[:-1]) * np.diff(r)
    u = np.concatenate([np.cumsum(seg[::-1])[::-1], [0.0]])
    q = np.trapezoid(2.0 * math.pi * r * u, r) if hasattr(np, "trapezoid") \
        else np.trapz(2.0 * math.pi * r * u, r)
    return r, u, sr, q


def latest_dir(path):
    dirs = []
    for d in path.iterdir():
        try:
            dirs.append((float(d.name), d))
        except ValueError:
            pass
    return max(dirs)[1] if dirs else None


def main():
    case = Path(__file__).resolve().parent
    prof_dir = case / "postProcessing" / "profile"
    flow_dir = case / "postProcessing" / "flowRate"
    if not prof_dir.is_dir() or not flow_dir.is_dir():
        sys.exit("No postProcessing data found: run ./Allrun first.")

    last = latest_dir(prof_dir)
    # 'sets' raw output for fields (U nu): columns  x  nu  Ux  Uy  Uz
    data = np.loadtxt(last / "radial_nu_U.xy")
    r_num = data[:, 0]
    nu_num = data[:, 1]
    u_num = data[:, 3]                      # axial (y) component

    r_ex, u_ex, sr_ex, q_ex = exact_profile()
    u_at = np.interp(r_num, r_ex, u_ex)
    nu_at = nu_carreau(np.interp(r_num, r_ex, sr_ex))
    u_err = np.max(np.abs(u_num - u_at)) / u_ex[0]
    # The shear rate vanishes on the axis, where the Carreau viscosity varies
    # too steeply to be represented by a cell average; check nu for r > R/5.
    outer = r_num > 0.2 * R
    nu_err = np.max(np.abs(nu_num[outer] - nu_at[outer]) / nu_at[outer])

    dat = sorted(flow_dir.glob("*/surfaceFieldValue*.dat"))[-1]
    q_num = np.loadtxt(dat, comments="#")[-1, 1]
    q_err = abs(q_num - q_ex) / q_ex

    # Newtonian reference (nu = nuInf) for context
    u_newt_c = G * R**2 / (4.0 * NUINF)
    q_newt = math.pi * G * R**4 / (8.0 * NUINF)

    print(f"Profile time                : {last.name}")
    print(f"Wall shear rate (exact)     : {sr_ex[-1]:.1f} 1/s")
    print(f"Centreline velocity         : OpenFOAM {u_num[0]:.5f}  exact {u_ex[0]:.5f} m/s"
          f"  (Newtonian nuInf: {u_newt_c:.5f})")
    print(f"Flow rate                   : OpenFOAM {q_num*6e7:.3f}  exact {q_ex*6e7:.3f} mL/min"
          f"  (Newtonian nuInf: {q_newt*6e7:.3f})")
    print(f"Max velocity error          : {100*u_err:.3f} % of centreline value")
    print(f"Max viscosity error, r>R/5  : {100*nu_err:.3f} % (local)")
    print(f"Flow-rate error             : {100*q_err:.3f} %")

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))
        ax1.plot(u_ex, r_ex / R, "k-", label="Carreau, exact")
        ax1.plot(G * (R**2 - r_ex**2) / (4 * NUINF), r_ex / R, "k--", lw=1,
                 label=r"Newtonian ($\nu_\infty$), exact")
        ax1.plot(u_num, r_num / R, "o", ms=3, label="Carreau, OpenFOAM")
        ax1.set_xlabel("u [m/s]")
        ax1.set_ylabel("r / R")
        ax1.set_title("Axial velocity")
        ax1.legend(fontsize=8)
        ax2.semilogx(nu_carreau(sr_ex) * 1060 * 1e3, r_ex / R, "k-", label="exact")
        ax2.semilogx(nu_num * 1060 * 1e3, r_num / R, "o", ms=3, label="OpenFOAM")
        ax2.set_xlabel("dynamic viscosity [mPa s]")
        ax2.set_title("Carreau viscosity")
        ax2.legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(case / "carreau_validation.png", dpi=150)
        print("Figure written              : carreau_validation.png")
    except ImportError:
        pass

    ok = max(u_err, nu_err, q_err) < TOLERANCE
    print("RESULT                      :", "PASS" if ok else "FAIL",
          f"(tolerance {100*TOLERANCE:g} %)")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
