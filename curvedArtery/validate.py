#!/usr/bin/env python3
"""Post-process and check the curvedArtery OpenFOAM result (Dean flow).

In a curved tube the centrifugal force drives a pair of counter-rotating
secondary (Dean) vortices, shifts the axial velocity peak towards the outer
wall and raises the friction factor above the straight-tube (Poiseuille)
value. The Dean number De = Re sqrt(R/Rc) controls these effects.

Checks (all must pass):
  1. Mass conservation between INLET and OUTLET.
  2. Straight inlet leg: Poiseuille pressure gradient 8 nu U / R^2 (< 3 %).
  3. Developing bend flow settles: at this Dean number the pressure gradient
     approaches its developed value through a damped oscillation, so the
     change between successive 45-degree intervals must shrink downstream
     (120-165 -> 165-210 -> 210-255 degrees).
  4. Curved-tube friction factor ratio f_c/f_s from 165-255 degrees agrees
     with the Mishra and Gupta (1979) correlation
     f_c/f_s = 1 + 0.033 (log10 De)^4 (< 5 %; the spread between published
     correlations at this Dean number is about 3 %).
  5. The axial velocity peak in the plane of the bend lies on the outer side.

Only numpy is required (matplotlib is optional and used for the figure).
"""
import math
import sys
from pathlib import Path

import numpy as np

D = 0.004
R = D / 2
RC = 7 * R
U = 0.2475
NU = 3.3e-6
RE = U * D / NU
DE = RE * math.sqrt(R / RC)


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


def section_p(pp, name):
    f = sorted((pp / name).glob("*/surfaceFieldValue*.dat"))[-1]
    return np.atleast_2d(np.loadtxt(f, comments="#"))[-1, 1]


def main():
    case = Path(__file__).resolve().parent
    pp = case / "postProcessing"
    if not (pp / "lines").is_dir():
        sys.exit("No postProcessing data found: run ./Allrun first.")

    q_in = abs(np.loadtxt(sorted((pp / "inletFlow").glob("*/*.dat"))[-1], comments="#")[-1, 1])
    q_out = np.loadtxt(sorted((pp / "outletFlow").glob("*/*.dat"))[-1], comments="#")[-1, 1]
    mass_err = abs(q_out - q_in) / q_in

    grad_pois = 8.0 * NU * U / R**2
    last = latest_dir(pp / "lines")
    leg = np.loadtxt(last / "inletLeg_p_U.xy")
    y, p = leg[:, 0], leg[:, 1]
    sel = (y > -0.012) & (y < -0.004)
    grad_leg = -np.polyfit(y[sel], p[sel], 1)[0]
    leg_err = abs(grad_leg - grad_pois) / grad_pois

    ps = {th: section_p(pp, f"theta{th:03d}") for th in (120, 165, 210, 255)}
    ds = RC * math.radians(45.0)
    g0 = (ps[120] - ps[165]) / ds
    g1 = (ps[165] - ps[210]) / ds
    g2 = (ps[210] - ps[255]) / ds
    dev_err = abs(g1 - g2) / g2
    dev_prev = abs(g0 - g1) / g1
    g_bend = (ps[165] - ps[255]) / (2 * ds)
    ratio = g_bend / grad_pois
    ratio_mg = 1.0 + 0.033 * math.log10(DE) ** 4
    ratio_white = 1.0 / (1.0 - (1.0 - (11.6 / DE) ** 0.45) ** 2.2)
    ratio_err = abs(ratio - ratio_mg) / ratio_mg

    bp = np.loadtxt(last / "bendPlane240_p_U.xy")       # distance p Ux Uy Uz
    dist = bp[:, 0]
    th = math.radians(240.0)
    t_vec = np.array([0.0, math.cos(th), math.sin(th)])
    u_ax = bp[:, 2:5] @ t_vec
    d_peak = dist[np.argmax(u_ax)]
    shift = (R - d_peak) / R                           # > 0: towards the outer wall

    print(f"Reynolds number Re = U D / nu : {RE:.0f}")
    print(f"Dean number De = Re sqrt(R/Rc): {DE:.1f}  (R/Rc = {R/RC:.4f})")
    print(f"Mass imbalance                : {100*mass_err:.4f} %")
    print(f"Inlet-leg -dp/ds              : {grad_leg:.4f} m/s^2 (Poiseuille {grad_pois:.4f}, "
          f"{100*leg_err:.2f} %)")
    print(f"Bend -dp/ds, 120-165 / 165-210 / 210-255 deg: {g0:.4f} / {g1:.4f} / {g2:.4f} m/s^2")
    print(f"   successive changes         : {100*dev_prev:.2f} % then {100*dev_err:.2f} %")
    print(f"Friction factor ratio f_c/f_s : {ratio:.4f}")
    print(f"   Mishra & Gupta (1979)      : {ratio_mg:.4f}  ({100*ratio_err:.2f} %)")
    print(f"   White (1929)               : {ratio_white:.4f}")
    print(f"Axial velocity peak (240 deg) : {shift:+.3f} R from the centreline "
          f"({'outer' if shift > 0 else 'inner'} side), u_max/U = {u_ax.max()/U:.3f}")

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(6, 4))
        ax.plot((dist - R) / R, u_ax / U, "b-", label="curved, 240 degrees")
        rr = np.linspace(-1, 1, 101)
        ax.plot(rr, 2 * (1 - rr**2), "k--", label="Poiseuille")
        ax.set_xlabel("position in the plane of the bend, a / R  (-1 outer wall, +1 inner wall)")
        ax.set_ylabel("u_axial / U")
        ax.set_title(f"Dean flow, De = {DE:.0f}")
        ax.legend()
        fig.tight_layout()
        fig.savefig(case / "curved_artery_results.png", dpi=150)
        print("Figure written                : curved_artery_results.png")
    except ImportError:
        pass

    checks = {
        "mass conservation (< 0.1 %)": mass_err < 1e-3,
        "inlet-leg Poiseuille gradient (< 3 %)": leg_err < 0.03,
        "bend pressure gradient settling (oscillation decays)": dev_err < dev_prev,
        "friction factor ratio vs Mishra & Gupta (< 5 %)": ratio_err < 0.05,
        "velocity peak shifted to the outer wall": shift > 0.05,
    }
    for name, ok in checks.items():
        print(f"  [{'PASS' if ok else 'FAIL'}] {name}")
    ok = all(checks.values())
    print("RESULT                        :", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
