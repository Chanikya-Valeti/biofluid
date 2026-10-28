#!/usr/bin/env python3
"""Post-process and check the fdaNozzle OpenFOAM result (Re_throat = 500).

Checks (all must pass):
  1. Mass conservation between INLET and OUTLET.
  2. Fully developed inflow: inlet-pipe centreline velocity = 2 U_inlet and
     Poiseuille pressure gradient 8 nu U / R^2 (both < 3 %).
  3. Peak centreline velocity at the sudden expansion compared with the
     published value for this benchmark condition: 0.76 m/s at a throat
     velocity of 0.4143 m/s (rho = 1056 kg/m^3, mu = 3.5 mPa s; Huang et al.
     2022, lattice-Boltzmann study of the FDA nozzle), i.e. 1.834 u_t
     (within 5 %).
  4. The jet stays axisymmetric: the cross-stream centreline velocity is
     below 1 % of the throat velocity everywhere.

Reported: centreline velocity along the nozzle, pressure drop from the inlet
pipe to the throat exit, and the jet length (distance downstream of the
expansion where the centreline velocity has decayed to half its peak).
Only numpy is required (matplotlib is optional and used for the figure).
"""
import sys
from pathlib import Path

import numpy as np

NU = 3.3e-6
D_IN = 0.012
R_IN = D_IN / 2
D_T = 0.004
RE_T = 500.0
U_T = RE_T * NU / D_T
U_IN = U_T * (D_T / D_IN) ** 2
REF_PEAK = 0.76 / 0.4143          # peak centreline velocity / throat mean velocity


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
    pp = case / "postProcessing"
    if not (pp / "lines").is_dir():
        sys.exit("No postProcessing data found: run ./Allrun first.")

    q_in = abs(np.loadtxt(sorted((pp / "inletFlow").glob("*/*.dat"))[-1], comments="#")[-1, 1])
    q_out = np.loadtxt(sorted((pp / "outletFlow").glob("*/*.dat"))[-1], comments="#")[-1, 1]
    mass_err = abs(q_out - q_in) / q_in

    last = latest_dir(pp / "lines")
    cl = np.loadtxt(last / "centreline_p_U.xy")          # y p Ux Uy Uz
    y, p, uy = cl[:, 0], cl[:, 1], cl[:, 3]
    ucross = np.hypot(cl[:, 2], cl[:, 4])

    up = (y > -0.110) & (y < -0.075)                      # inlet pipe, upstream of the cone
    u_inlet_cl = uy[up].mean()
    inlet_err = abs(u_inlet_cl - 2 * U_IN) / (2 * U_IN)
    grad = -np.polyfit(y[up], p[up], 1)[0]
    grad_pois = 8 * NU * U_IN / R_IN**2
    grad_err = abs(grad - grad_pois) / grad_pois

    u_exp = float(np.interp(0.0, y, uy))
    peak_ratio = u_exp / U_T
    peak_err = abs(peak_ratio - REF_PEAK) / REF_PEAK
    cross = ucross.max() / U_T

    down = y > 0
    half = np.where(down & (uy < 0.5 * u_exp))[0]
    jet_len = y[half[0]] if len(half) else None
    dp = float(np.interp(-0.08, y, p) - np.interp(-0.001, y, p)) * 1056.0

    print(f"Throat Reynolds number        : {RE_T:.0f} (u_t = {U_T:.4f} m/s, inlet U = {U_IN:.5f} m/s)")
    print(f"Mass imbalance                : {100*mass_err:.4f} %")
    print(f"Inlet centreline velocity     : {u_inlet_cl:.5f} m/s (2 U = {2*U_IN:.5f}, {100*inlet_err:.2f} %)")
    print(f"Inlet -dp/dy                  : {grad:.5f} m/s^2 (Poiseuille {grad_pois:.5f}, {100*grad_err:.2f} %)")
    print(f"Centreline velocity at the expansion: {u_exp:.4f} m/s = {peak_ratio:.3f} u_t "
          f"(reference {REF_PEAK:.3f} u_t, {100*peak_err:.2f} %)")
    print(f"Max cross-stream centreline velocity: {100*cross:.3f} % of u_t")
    print(f"Pressure drop, y = -80 mm to throat exit: {dp:.1f} Pa (rho = 1056 kg/m^3)")
    if jet_len is not None:
        print(f"Jet half-decay length         : {jet_len*1e3:.1f} mm = {jet_len/D_T:.1f} d downstream of the expansion")
    else:
        print("Jet half-decay length         : beyond the outlet")

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, a = plt.subplots(figsize=(8, 3.5))
        a.plot(y / D_T, uy / U_T, "b-")
        a.plot([0], [REF_PEAK], "rs", label="published peak (Huang et al. 2022)")
        a.axvspan(-62.685 / 4, -10, color="0.92", label="cone")
        a.axvspan(-10, 0, color="0.85", label="throat")
        a.set_xlabel("y / d (sudden expansion at 0)")
        a.set_ylabel("centreline u / u_t")
        a.set_title("FDA nozzle, Re_t = 500")
        a.legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(case / "fda_nozzle_centreline.png", dpi=150)
        print("Figure written                : fda_nozzle_centreline.png")
    except ImportError:
        pass

    checks = {
        "mass conservation (< 0.1 %)": mass_err < 1e-3,
        "fully developed inflow (< 3 %)": inlet_err < 0.03 and grad_err < 0.03,
        "peak centreline velocity vs published value (< 5 %)": peak_err < 0.05,
        "axisymmetric jet (cross-stream < 1 % of u_t)": cross < 0.01,
    }
    for name, ok in checks.items():
        print(f"  [{'PASS' if ok else 'FAIL'}] {name}")
    ok = all(checks.values())
    print("RESULT                        :", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
