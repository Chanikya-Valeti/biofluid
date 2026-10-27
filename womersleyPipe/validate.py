#!/usr/bin/env python3
"""Compare the womersleyPipe OpenFOAM results with the exact Womersley solution.

Fully developed flow in a rigid tube driven by -dp/dx = G0 cos(omega t) has

    u(r, t) = Re{ G0/(i omega) [1 - J0(L r/R) / J0(L)] e^{i omega t} },
    Q(t)    = Re{ G0 pi R^2/(i omega) [1 - 2 J1(L) / (L J0(L))] e^{i omega t} },

with L = i^{3/2} alpha and Womersley number alpha = R sqrt(omega/nu).

Only numpy is required (matplotlib is optional and used for the figure).
Exit status is 0 when every error is below the tolerance, 1 otherwise.
"""
import math
import sys
from pathlib import Path

import numpy as np

# Case parameters (keep in sync with 0/p, constant/transportProperties and
# system/blockMeshDict)
R = 0.004            # tube radius [m]
NU = 3.3e-6          # kinematic viscosity [m^2/s]
FREQ = 1.0           # [Hz]
G0 = 0.04 / 0.02     # kinematic pressure-gradient amplitude [m/s^2]
COMPARE_FROM = 4.0   # compare only the final cycle (t > 4 s)
TOLERANCE = 0.02     # 2 % of the peak value

OMEGA = 2.0 * math.pi * FREQ
ALPHA = R * math.sqrt(OMEGA / NU)
LAMBDA = (1j) ** 1.5 * ALPHA


def bessel_j(n, z, terms=80):
    """Bessel function of the first kind J_n(z) for complex z (power series)."""
    z = np.asarray(z, dtype=complex)
    total = np.zeros_like(z)
    term = (z / 2.0) ** n / math.factorial(n)
    for k in range(terms):
        total += term
        term = term * (-(z / 2.0) ** 2) / ((k + 1) * (k + 1 + n))
    return total


def u_exact(r, t):
    shape = 1.0 - bessel_j(0, LAMBDA * r / R) / bessel_j(0, LAMBDA)
    return np.real(G0 / (1j * OMEGA) * shape * np.exp(1j * OMEGA * t))


def q_exact(t):
    factor = 1.0 - 2.0 * bessel_j(1, LAMBDA) / (LAMBDA * bessel_j(0, LAMBDA))
    return np.real(G0 * math.pi * R**2 / (1j * OMEGA) * factor * np.exp(1j * OMEGA * t))


def main():
    case = Path(__file__).resolve().parent
    profile_dir = case / "postProcessing" / "profile"
    flow_dir = case / "postProcessing" / "flowRate"
    if not profile_dir.is_dir() or not flow_dir.is_dir():
        sys.exit("No postProcessing data found: run ./Allrun first.")

    # --- Velocity profiles --------------------------------------------------
    times, profiles = [], []
    for d in profile_dir.iterdir():
        try:
            t = float(d.name)
        except ValueError:
            continue
        if t > COMPARE_FROM + 1e-9:
            data = np.loadtxt(d / "radial_U.xy")  # columns: x  Ux  Uy  Uz (axis = y)
            times.append(t)
            profiles.append(data)
    if not times:
        sys.exit("No velocity profiles from the final cycle were found.")
    order = np.argsort(times)
    times = [times[i] for i in order]
    profiles = [profiles[i] for i in order]

    u_peak = max(np.max(np.abs(u_exact(p[:, 0], t))) for t, p in zip(times, profiles))
    profile_err = 0.0
    for t, p in zip(times, profiles):
        err = np.max(np.abs(p[:, 2] - u_exact(p[:, 0], t))) / u_peak
        profile_err = max(profile_err, err)

    # --- Flow rate ----------------------------------------------------------
    dat = sorted(flow_dir.glob("*/surfaceFieldValue*.dat"))[-1]
    q = np.loadtxt(dat, comments="#")
    mask = q[:, 0] > COMPARE_FROM + 1e-9
    tq, qn = q[mask, 0], q[mask, 1]
    qe = q_exact(tq)
    flow_err = np.max(np.abs(qn - qe)) / np.max(np.abs(qe))

    print(f"Womersley number alpha     : {ALPHA:.3f}")
    print(f"Profiles compared          : {len(times)} (t = {times[0]:g} ... {times[-1]:g} s)")
    print(f"Peak exact velocity        : {u_peak:.5f} m/s")
    print(f"Max velocity error         : {100 * profile_err:.3f} % of peak")
    print(f"Peak exact flow rate       : {np.max(np.abs(qe)) * 6e7:.3f} mL/min")
    print(f"Max flow-rate error        : {100 * flow_err:.3f} % of peak")

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))
        rr = np.linspace(0, R, 200)
        for t, p in list(zip(times, profiles))[::2]:
            line, = ax1.plot(u_exact(rr, t), rr / R, lw=1)
            ax1.plot(p[:, 2], p[:, 0] / R, "o", ms=2.5, color=line.get_color())
        ax1.set_xlabel("u [m/s]")
        ax1.set_ylabel("r / R")
        ax1.set_title(f"Velocity profiles (alpha = {ALPHA:.2f})\nlines: exact, dots: OpenFOAM")
        ax2.plot(tq, qe * 6e7, "k-", label="exact")
        ax2.plot(tq[::10], qn[::10] * 6e7, "o", ms=3, label="OpenFOAM")
        ax2.set_xlabel("t [s]")
        ax2.set_ylabel("Q [mL/min]")
        ax2.set_title("Flow rate, final cycle")
        ax2.legend()
        fig.tight_layout()
        fig.savefig(case / "womersley_validation.png", dpi=150)
        print("Figure written             : womersley_validation.png")
    except ImportError:
        pass

    ok = profile_err < TOLERANCE and flow_err < TOLERANCE
    print("RESULT                     :", "PASS" if ok else "FAIL",
          f"(tolerance {100 * TOLERANCE:g} %)")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
