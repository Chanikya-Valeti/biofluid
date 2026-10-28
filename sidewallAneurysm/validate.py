#!/usr/bin/env python3
"""Post-process and check the sidewallAneurysm OpenFOAM result.

Checks (all must pass):
  1. Mass conservation between INLET and OUTLET.
  2. Parent artery upstream of the aneurysm is Poiseuille flow: pressure
     gradient 8 nu U / R^2 and wall shear stress 4 nu U / R (< 5 %).
  3. Downstream of the aneurysm the wall shear stress has recovered to the
     Poiseuille value (< 5 %).
  4. The sac is a low-shear region: its mean wall shear stress is below the
     parent-artery value.
  5. The slow vortex in the sac has converged: the velocity at two probes in
     the sac changed by less than 1 % over the last 300 iterations.

Reported: flow exchanged through the neck (in and out; the two differ by a
few per cent because they are integrated from cell values on a plane that
cuts through cells), mean WSS in the sac and on the dome, the difference
between the two halves of the sac (the geometry is symmetric about x = 0 but
the snappyHexMesh mesh is not, and the sac shear is small enough for this to
show), and the area of the sac with WSS below 0.4 Pa (a level
often associated with aneurysm wall remodelling).
Kinematic stresses are converted with rho = 1060 kg/m^3.
Only numpy is required (matplotlib is optional and used for the figure).
"""
import math
import sys
from pathlib import Path

import numpy as np

RHO = 1060.0
NU = 3.3e-6
R = 0.002
D = 2 * R
U = 0.2475
RA = 0.003
C = R + math.sqrt(RA**2 - 0.002**2)


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


def last_row(pp, name):
    f = sorted((pp / name).glob("*/surfaceFieldValue*.dat"))[-1]
    return np.atleast_2d(np.loadtxt(f, comments="#"))[-1]


def main():
    case = Path(__file__).resolve().parent
    pp = case / "postProcessing"
    if not (pp / "wallData").is_dir():
        sys.exit("No postProcessing data found: run ./Allrun first.")

    q_in = abs(last_row(pp, "inletFlow")[1])
    q_out = last_row(pp, "outletFlow")[1]
    mass_err = abs(q_out - q_in) / q_in

    neck = last_row(pp, "neckFlow")
    neck_in, neck_out = neck[1], neck[2]

    tau_pois = 4.0 * NU * U / R
    grad_pois = 8.0 * NU * U / R**2

    last = latest_dir(pp / "wallData")
    w = np.loadtxt(last / "wallShearStress_wall.raw", comments="#")
    xyz, wss = w[:, :3], np.linalg.norm(w[:, 3:6], axis=1)
    y, z = xyz[:, 1], xyz[:, 2]
    parent = (y > -5.5 * D) & (y < -2.0 * D)
    recovered = (y > 8.0 * D) & (y < 11.0 * D)
    sac = z > R + 1e-4
    dome = z > C
    tau_parent = wss[parent].mean()
    tau_down = wss[recovered].mean()
    parent_err = abs(tau_parent - tau_pois) / tau_pois
    down_err = abs(tau_down - tau_pois) / tau_pois
    sac_mean = wss[sac].mean()
    half_a = wss[sac & (xyz[:, 0] > 0)].mean()
    half_b = wss[sac & (xyz[:, 0] < 0)].mean()
    sym_err = abs(half_a - half_b) / (0.5 * (half_a + half_b))
    low_frac = np.mean(wss[sac] * RHO < 0.4)          # faces have similar areas

    probe_t, probe_u = [], []
    for line in sorted((pp / "sacProbes").glob("*/U"))[-1].read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        nums = [float(x) for x in line.replace("(", " ").replace(")", " ").split()]
        probe_t.append(nums[0])
        probe_u.append(np.array(nums[1:]).reshape(-1, 3))
    probe_t, probe_u = np.array(probe_t), np.array(probe_u)
    i_old = np.argmin(np.abs(probe_t - (probe_t[-1] - 300)))
    mags = np.linalg.norm(probe_u, axis=2)
    sac_change = np.max(np.abs(mags[-1] - mags[i_old]) / mags[-1])

    lines = latest_dir(pp / "lines")
    ax = np.loadtxt(lines / "parentAxis_p_U.xy")
    sel = (ax[:, 0] > -5.5 * D) & (ax[:, 0] < -2.5 * D)
    grad = -np.polyfit(ax[sel, 0], ax[sel, 1], 1)[0]
    grad_err = abs(grad - grad_pois) / grad_pois

    print(f"Reynolds number               : {U*D/NU:.0f}")
    print(f"Mass imbalance                : {100*mass_err:.4f} %")
    print(f"Parent -dp/dy                 : {grad:.4f} m/s^2 (Poiseuille {grad_pois:.4f}, "
          f"{100*grad_err:.2f} %)")
    print(f"Parent WSS upstream / downstream: {tau_parent*RHO:.4f} / {tau_down*RHO:.4f} Pa "
          f"(Poiseuille {tau_pois*RHO:.4f}; {100*parent_err:.2f} % / {100*down_err:.2f} %)")
    print(f"Neck inflow / outflow         : {neck_in*6e7:.3f} / {neck_out*6e7:.3f} mL/min")
    print(f"Fraction of parent flow entering the sac: {100*0.5*(neck_in+neck_out)/q_in:.2f} %")
    print(f"Sac WSS, x > 0 / x < 0 halves : {half_a*RHO:.4f} / {half_b*RHO:.4f} Pa "
          f"({100*sym_err:.2f} % apart)")
    print(f"Mean WSS sac / dome           : {sac_mean*RHO:.4f} / {wss[dome].mean()*RHO:.4f} Pa "
          f"({sac_mean/tau_pois:.2f} / {wss[dome].mean()/tau_pois:.2f} x parent)")
    print(f"Sac wall with WSS < 0.4 Pa    : {100*low_frac:.1f} % of faces")
    print(f"Sac probe speeds (centre, dome): {mags[-1, 0]*1e3:.3f}, {mags[-1, 1]*1e3:.3f} mm/s "
          f"(change over 300 iterations {100*sac_change:.3f} %)")

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, a = plt.subplots(figsize=(8, 3.2))
        side = xyz[:, 0] > 0
        sc = a.scatter(y[side] * 1e3, z[side] * 1e3, c=wss[side] * RHO, s=2, cmap="viridis",
                       vmin=0, vmax=np.percentile(wss * RHO, 99))
        a.set_aspect("equal")
        a.set_xlim(-12, 16)
        a.set_xlabel("y [mm]"); a.set_ylabel("z [mm]")
        a.set_title("Wall shear stress (x > 0 half)")
        fig.colorbar(sc, ax=a, label="WSS [Pa]", shrink=0.8)
        fig.tight_layout()
        fig.savefig(case / "aneurysm_wss.png", dpi=150)
        print("Figure written                : aneurysm_wss.png")
    except ImportError:
        pass

    checks = {
        "mass conservation (< 0.5 %)": mass_err < 0.005,
        "parent Poiseuille pressure gradient (< 5 %)": grad_err < 0.05,
        "parent Poiseuille WSS upstream (< 5 %)": parent_err < 0.05,
        "WSS recovered downstream (< 5 %)": down_err < 0.05,
        "low-shear sac (mean WSS < parent)": sac_mean < tau_parent,
        "sac flow converged (< 1 % change)": sac_change < 0.01,
    }
    for name, ok in checks.items():
        print(f"  [{'PASS' if ok else 'FAIL'}] {name}")
    ok = all(checks.values())
    print("RESULT                        :", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
