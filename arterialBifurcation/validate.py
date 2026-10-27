#!/usr/bin/env python3
"""Post-process and check the arterialBifurcation OpenFOAM result.

Checks (all must pass):
  1. Mass conservation: Q_in = Q_out1 + Q_out2.
  2. Symmetric flow split: Q_out1 = Q_out2 (identical branches and outlet
     pressures).
  3. Upstream parent vessel is Poiseuille flow: axial pressure gradient
     8 nu U / Rp^2 and wall shear stress 4 nu U / Rp.
  4. Murray's law: with Rd = Rp / 2^(1/3) and an even split, fully developed
     Poiseuille flow in each daughter has the same wall shear stress as the
     parent. The mean WSS far down each daughter is compared with the parent
     value.

Reported: peak WSS and its location (expected near the flow divider/apex).
Kinematic quantities from OpenFOAM are converted with rho = 1060 kg/m^3.
Only numpy is required (matplotlib is optional and used for the figure).
"""
import math
import sys
from pathlib import Path

import numpy as np

RHO = 1060.0
NU = 3.3e-6
RP = 0.003
DP = 2 * RP
RD = RP / 2 ** (1.0 / 3.0)
DD = 2 * RD
U_MEAN = 0.0825
HALF_ANGLE = math.radians(35.0)
D1 = np.array([0.0, math.cos(HALF_ANGLE), math.sin(HALF_ANGLE)])
D2 = np.array([0.0, math.cos(HALF_ANGLE), -math.sin(HALF_ANGLE)])


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


def last_value(pp, name):
    dat = sorted((pp / name).glob("*/surfaceFieldValue*.dat"))[-1]
    return np.loadtxt(dat, comments="#")[-1, 1]


def main():
    case = Path(__file__).resolve().parent
    pp = case / "postProcessing"
    if not (pp / "wallData").is_dir():
        sys.exit("No postProcessing data found: run ./Allrun first.")

    q_in = abs(last_value(pp, "inletFlow"))
    q1 = last_value(pp, "outlet1Flow")
    q2 = last_value(pp, "outlet2Flow")
    mass_err = abs(q1 + q2 - q_in) / q_in
    split_err = abs(q1 - q2) / (q1 + q2)

    last = latest_dir(pp / "wallData")
    w = np.loadtxt(last / "wallShearStress_wall.raw", comments="#")
    xyz, wss = w[:, :3], np.linalg.norm(w[:, 3:6], axis=1)

    # Parent vessel, away from the inlet and the junction
    y = xyz[:, 1]
    parent = (y > -3.5 * DP) & (y < -1.5 * DP)
    tau_parent = wss[parent].mean()
    tau_pois = 4.0 * NU * U_MEAN / RP
    parent_err = abs(tau_parent - tau_pois) / tau_pois

    # Daughters: faces 7.5 - 9.5 Dd along each branch axis
    daughter_tau = []
    for d, side in ((D1, xyz[:, 2] > 0), (D2, xyz[:, 2] < 0)):
        s = xyz @ d
        sel = side & (s > 7.5 * DD) & (s < 9.5 * DD)
        daughter_tau.append(wss[sel].mean())
    murray_err = max(abs(t - tau_pois) / tau_pois for t in daughter_tau)

    lines = latest_dir(pp / "lines")
    ax = np.loadtxt(lines / "parentAxis_p_U.xy")
    grad_num = (np.interp(-3.5 * DP, ax[:, 0], ax[:, 1])
                - np.interp(-1.5 * DP, ax[:, 0], ax[:, 1])) / (2.0 * DP)
    grad_pois = 8.0 * NU * U_MEAN / RP**2
    grad_err = abs(grad_num - grad_pois) / grad_pois

    i_max = int(np.argmax(wss))

    print(f"Iteration                    : {last.name}")
    print(f"Reynolds number (parent)     : {U_MEAN * DP / NU:.0f}")
    print(f"Flow rates [mL/min]          : in {q_in*6e7:.3f}, out1 {q1*6e7:.3f}, out2 {q2*6e7:.3f}")
    print(f"Mass imbalance               : {100*mass_err:.4f} %")
    print(f"Flow split asymmetry         : {100*split_err:.4f} %")
    print(f"Parent -dp/dy                : {grad_num:.5f} m/s^2 (Poiseuille {grad_pois:.5f}, "
          f"{100*grad_err:.2f} %)")
    print(f"Parent mean WSS              : {tau_parent*RHO:.4f} Pa (Poiseuille {tau_pois*RHO:.4f}, "
          f"{100*parent_err:.2f} %)")
    print(f"Daughter mean WSS (1, 2)     : {daughter_tau[0]*RHO:.4f}, {daughter_tau[1]*RHO:.4f} Pa "
          f"(Murray {tau_pois*RHO:.4f}, max dev {100*murray_err:.2f} %)")
    print(f"Peak WSS                     : {wss[i_max]*RHO:.4f} Pa "
          f"({wss[i_max]/tau_pois:.2f} x parent) at "
          f"({xyz[i_max,0]*1e3:.2f}, {xyz[i_max,1]*1e3:.2f}, {xyz[i_max,2]*1e3:.2f}) mm")

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, ax1 = plt.subplots(figsize=(6, 6))
        front = xyz[:, 0] > 0
        sc = ax1.scatter(xyz[front, 2] * 1e3, xyz[front, 1] * 1e3, c=wss[front] * RHO,
                         s=3, cmap="viridis", vmin=0, vmax=np.percentile(wss * RHO, 99.5))
        ax1.set_aspect("equal")
        ax1.set_xlabel("z [mm]")
        ax1.set_ylabel("y [mm]")
        ax1.set_title("Wall shear stress (x > 0 half)")
        fig.colorbar(sc, ax=ax1, label="WSS [Pa]", shrink=0.8)
        fig.tight_layout()
        fig.savefig(case / "bifurcation_wss.png", dpi=150)
        print("Figure written               : bifurcation_wss.png")
    except ImportError:
        pass

    checks = {
        "mass conservation (< 0.5 %)": mass_err < 0.005,
        "symmetric flow split (< 1 %)": split_err < 0.01,
        "parent Poiseuille pressure gradient (< 5 %)": grad_err < 0.05,
        "parent Poiseuille WSS (< 5 %)": parent_err < 0.05,
        "Murray's law daughter WSS (< 7 %)": murray_err < 0.07,
    }
    for name, ok in checks.items():
        print(f"  [{'PASS' if ok else 'FAIL'}] {name}")
    ok = all(checks.values())
    print("RESULT                       :", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
