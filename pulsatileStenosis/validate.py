#!/usr/bin/env python3
"""Post-process and check the pulsatileStenosis OpenFOAM result.

Checks (all must pass):
  1. Mass conservation: inlet and outlet flow rates agree at every time step
     of the last cycle.
  2. Upstream flow: the exact Womersley profile imposed at the inlet is
     carried unchanged to the probe 4 D upstream of the throat (3 D into the
     domain): its centreline velocity follows the Womersley solution for
     Q(t) = Qm [1 + 0.6 sin(2 pi t)] within 3 % of its peak.
  3. Periodicity: the spatially averaged wall shear stress in cycle 3 matches
     cycle 2 (within 2 % of its peak), so the start-up transient has gone.
  4. Physical consistency of the cycle-averaged indices: 0 <= OSI <= 0.5
     everywhere, and the post-stenotic region has oscillatory shear
     (OSI > 0.1 somewhere downstream of the throat).

Reported: peak TAWSS at the throat, TAWSS upstream, peak OSI and its location.
Kinematic stresses are converted with rho = 1060 kg/m^3.
Only numpy is required (matplotlib is optional and used for the figure).
"""
import math
import re
import sys
from pathlib import Path

import numpy as np

RHO = 1060.0
NU = 3.3e-6
D = 0.006
R0 = D / 2
QM = 1.5551e-6
Q1 = 0.6 * QM
FREQ = 1.0
OMEGA = 2 * math.pi * FREQ
ALPHA = R0 * math.sqrt(OMEGA / NU)
LAM = (1j) ** 1.5 * ALPHA


def bessel_j(n, z, terms=80):
    z = complex(z)
    total, term = 0j, (z / 2) ** n / math.factorial(n)
    for k in range(terms):
        total += term
        term *= -(z / 2) ** 2 / ((k + 1) * (k + 1 + n))
    return total


def womersley_centreline(t):
    """Centreline velocity for Q(t) = QM + Q1 sin(omega t)."""
    u_mean = 2.0 * QM / (math.pi * R0**2)
    q_hat = -1j * Q1                                 # Q1 sin = Re(-i Q1 e^{iwt})
    f = 1.0 - 2.0 * bessel_j(1, LAM) / (LAM * bessel_j(0, LAM))
    g_hat = q_hat * 1j * OMEGA / (math.pi * R0**2 * f)
    u_hat = g_hat / (1j * OMEGA) * (1.0 - 1.0 / bessel_j(0, LAM))
    return u_mean + np.real(u_hat * np.exp(1j * OMEGA * np.asarray(t)))


def load_dat(path):
    return np.loadtxt(path, comments="#")


def load_probes(path):
    """OpenFOAM probes file for a vector field -> (time, array[n_t, n_probe, 3])."""
    times, vals = [], []
    for line in path.read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        nums = [float(x) for x in re.findall(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", line)]
        times.append(nums[0])
        vals.append(np.array(nums[1:]).reshape(-1, 3))
    return np.array(times), np.array(vals)


def cycle(t, f, k):
    sel = (t > k - 1 + 1e-9) & (t <= k + 1e-9)
    return t[sel] - (k - 1), f[sel]


def main():
    case = Path(__file__).resolve().parent
    pp = case / "postProcessing"
    if not (pp / "probes").is_dir():
        sys.exit("No postProcessing data found: run ./Allrun first.")

    qi = load_dat(sorted((pp / "inletFlow").glob("*/surfaceFieldValue*.dat"))[-1])
    qo = load_dat(sorted((pp / "outletFlow").glob("*/surfaceFieldValue*.dat"))[-1])
    n = min(len(qi), len(qo))
    t_q, q_in, q_out = qi[:n, 0], -qi[:n, 1], qo[:n, 1]
    if t_q[-1] < 3.0 - 1e-6:
        sys.exit(f"The run stopped at t = {t_q[-1]:.3f} s; three cycles (t = 3 s) are needed.")
    sel = t_q > 2.0
    mass_err = np.max(np.abs(q_in[sel] - q_out[sel])) / QM

    t_p, u_p = load_probes(sorted((pp / "probes").glob("*/U"))[-1])
    t3, u3 = cycle(t_p, u_p[:, 0, 1], 3)
    u_ex = womersley_centreline(t3)
    wom_err = np.max(np.abs(u3 - u_ex)) / np.max(np.abs(u_ex))
    _, u_throat = cycle(t_p, u_p[:, 1, 1], 3)

    wss = load_dat(sorted((pp / "wallMeanWSS").glob("*/surfaceFieldValue*.dat"))[-1])
    t2w, w2 = cycle(wss[:, 0], wss[:, 1], 2)
    t3w, w3 = cycle(wss[:, 0], wss[:, 1], 3)
    w2i = np.interp(t3w, t2w, w2)
    period_err = np.max(np.abs(w3 - w2i)) / np.max(np.abs(w3))

    wd = sorted((pp / "wallData").glob("*"))[-1]
    mean_vec = np.loadtxt(wd / "wallShearStressMean_wall.raw", comments="#")
    mean_mag = np.loadtxt(wd / "mag(wallShearStress)Mean_wall.raw", comments="#")
    xyz = mean_vec[:, :3]
    tawss = mean_mag[:, 3]
    osi = 0.5 * (1.0 - np.linalg.norm(mean_vec[:, 3:6], axis=1) / np.maximum(tawss, 1e-30))
    y = xyz[:, 1]
    up = (y > -5 * D) & (y < -3 * D)
    down = y > 0.0
    i_osi = int(np.argmax(np.where(down, osi, -1)))

    print(f"Womersley number alpha        : {ALPHA:.2f}")
    print(f"Mean / peak Reynolds number   : {2*QM/(math.pi*R0*NU):.0f} / "
          f"{2*(QM+Q1)/(math.pi*R0*NU):.0f}")
    print(f"Mass imbalance (cycle 3, max) : {100*mass_err:.4f} % of mean flow")
    print(f"Centreline 4 D upstream vs Womersley (cycle 3, max): {100*wom_err:.2f} % of peak")
    print(f"Cycle 3 vs cycle 2 mean WSS   : {100*period_err:.3f} % of peak")
    print(f"Peak throat centreline velocity: {u_throat.max():.4f} m/s")
    print(f"TAWSS upstream (mean)         : {tawss[up].mean()*RHO:.4f} Pa "
          f"(steady Poiseuille at Qm: {4*NU*QM/(math.pi*R0**3)*RHO:.4f} Pa)")
    print(f"TAWSS peak (throat region)    : {tawss.max()*RHO:.4f} Pa at y = {y[np.argmax(tawss)]/D:.2f} D")
    print(f"OSI peak downstream           : {osi[i_osi]:.3f} at y = {y[i_osi]/D:.2f} D")

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, axs = plt.subplots(1, 3, figsize=(14, 3.8))
        tt = np.linspace(0, 1, 200)
        axs[0].plot(tt, womersley_centreline(tt), "k-", label="Womersley (exact)")
        axs[0].plot(t3[::10], u3[::10], "o", ms=3, label="OpenFOAM, cycle 3")
        axs[0].set_xlabel("t / T"); axs[0].set_ylabel("u_y [m/s]")
        axs[0].set_title("Centreline, 4 D upstream"); axs[0].legend(fontsize=8)
        order = np.argsort(y)
        axs[1].plot(y[order] / D, tawss[order] * RHO, ".", ms=1)
        axs[1].set_xlabel("y / D"); axs[1].set_ylabel("TAWSS [Pa]"); axs[1].set_title("Time-averaged WSS")
        axs[2].plot(y[order] / D, osi[order], ".", ms=1, color="C3")
        axs[2].set_xlabel("y / D"); axs[2].set_ylabel("OSI"); axs[2].set_title("Oscillatory shear index")
        for a in axs[1:]:
            a.axvspan(-1, 1, color="0.9")
        fig.tight_layout()
        fig.savefig(case / "pulsatile_stenosis_results.png", dpi=150)
        print("Figure written                : pulsatile_stenosis_results.png")
    except ImportError:
        pass

    checks = {
        "mass conservation (< 0.1 %)": mass_err < 1e-3,
        "Womersley inflow (< 3 %)": wom_err < 0.03,
        "cycle-to-cycle periodicity (< 2 %)": period_err < 0.02,
        "0 <= OSI <= 0.5": osi.min() >= -1e-6 and osi.max() <= 0.5 + 1e-6,
        "oscillatory shear downstream (OSI > 0.1)": osi[i_osi] > 0.1,
    }
    for name, ok in checks.items():
        print(f"  [{'PASS' if ok else 'FAIL'}] {name}")
    ok = all(checks.values())
    print("RESULT                        :", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
