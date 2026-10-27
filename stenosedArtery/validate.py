#!/usr/bin/env python3
"""Post-process and check the stenosedArtery OpenFOAM result.

Checks (all must pass):
  1. Mass conservation between INLET and OUTLET.
  2. Fully developed inlet flow: centreline velocity = 2 U_mean (Poiseuille).
  3. Upstream axial pressure gradient = 8 nu U_mean / R0^2 (Poiseuille).
  4. A post-stenotic recirculation zone exists and reattaches inside the domain.
  5. Far downstream the centreline velocity has recovered to 2 U_mean.

Checks 2, 3 and 5 compare with Newtonian Poiseuille flow, so for
./Allrun carreau the script reports these quantities without testing them.

Reported: throat jet velocity, reattachment point and stenosis pressure loss.
Only numpy is required (matplotlib is optional and used for the figure).
"""
import sys
from pathlib import Path

import numpy as np

D = 0.006
R0 = D / 2.0
U_MEAN = 0.055
NU = 3.3e-6
RE = U_MEAN * D / NU


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


def is_newtonian(case):
    local = case / "constant" / "transportProperties.local"
    if local.exists():                       # written by ./Allrun carreau
        text = local.read_text()
    else:
        text = (case / "constant" / "transportProperties").read_text()
    for line in text.splitlines():
        line = line.split("//")[0].split()
        if len(line) >= 2 and line[0] == "transportModel":
            return line[1].rstrip(";") == "Newtonian"
    return True


def main():
    case = Path(__file__).resolve().parent
    pp = case / "postProcessing"
    if not (pp / "lines").is_dir():
        sys.exit("No postProcessing data found: run ./Allrun first.")
    newtonian = is_newtonian(case)

    last = latest_dir(pp / "lines")
    cl = np.loadtxt(last / "centreline_p_U.xy")      # y p Ux Uy Uz
    nw = np.loadtxt(last / "nearWall_p_U.xy")
    y, p, uy = cl[:, 0], cl[:, 1], cl[:, 3]

    def at(yq, f):
        return float(np.interp(yq, y, f))

    q_in = np.loadtxt(sorted((pp / "inletFlow").glob("*/surfaceFieldValue*.dat"))[-1])[-1, 1]
    q_out = np.loadtxt(sorted((pp / "outletFlow").glob("*/surfaceFieldValue*.dat"))[-1])[-1, 1]
    mass_err = abs(abs(q_out) - abs(q_in)) / abs(q_in)

    u_in = at(-3.5 * D, uy)
    dev_in = abs(u_in - 2 * U_MEAN) / (2 * U_MEAN)
    grad_num = (at(-3.5 * D, p) - at(-2.0 * D, p)) / (1.5 * D)
    grad_exact = 8.0 * NU * U_MEAN / R0**2
    grad_err = abs(grad_num - grad_exact) / grad_exact
    u_out = at(15.5 * D, uy)
    dev_out = abs(u_out - 2 * U_MEAN) / (2 * U_MEAN)
    u_throat = at(0.0, uy)

    # Reattachment: last sign change (negative -> positive) of the near-wall
    # axial velocity downstream of the stenosis
    yw, uw = nw[:, 0], nw[:, 3]
    neg = uw < 0
    reattach = None
    if neg.any():
        i_last_neg = np.where(neg)[0][-1]
        if i_last_neg + 1 < len(yw):
            y0, y1 = yw[i_last_neg], yw[i_last_neg + 1]
            u0, u1 = uw[i_last_neg], uw[i_last_neg + 1]
            reattach = y0 - u0 * (y1 - y0) / (u1 - u0)

    # Pressure drop between -2D and +12D in excess of the fully developed drop
    # over the same length (using the measured upstream gradient, so the
    # figure is meaningful for the Carreau model too)
    dp_total = at(-2.0 * D, p) - at(12.0 * D, p)
    dp_extra = dp_total - grad_num * 14.0 * D

    print(f"Model                        : {'Newtonian' if newtonian else 'non-Newtonian'} "
          f"(Re = {RE:.0f} based on nu = {NU:g})")
    print(f"Iteration                    : {last.name}")
    print(f"Mass imbalance               : {100 * mass_err:.4f} %")
    print(f"Inlet centreline velocity    : {u_in:.5f} m/s (Poiseuille {2*U_MEAN:.5f}, "
          f"{100*dev_in:.2f} %)")
    print(f"Upstream -dp/dy              : {grad_num:.5f} m/s^2 (Poiseuille {grad_exact:.5f}, "
          f"{100*grad_err:.2f} %)")
    print(f"Throat centreline velocity   : {u_throat:.5f} m/s ({u_throat/U_MEAN:.2f} U_mean)")
    if reattach is not None:
        print(f"Reattachment point           : y = {reattach/D:.2f} D downstream of the throat")
    else:
        print("Reattachment point           : not found")
    print(f"Outlet-region centreline vel.: {u_out:.5f} m/s (Poiseuille {2*U_MEAN:.5f}, "
          f"{100*dev_out:.2f} %)")
    print(f"Stenosis pressure loss       : {dp_extra*1060:.2f} Pa above the fully developed "
          f"drop (rho = 1060 kg/m^3)")

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 6), sharex=True)
        ax1.plot(y / D, uy / U_MEAN, "b-", label="centreline")
        ax1.plot(yw / D, uw / U_MEAN, "r-", label="r = 0.97 R0")
        ax1.axhline(0, color="k", lw=0.5)
        ax1.axhline(2, color="k", lw=0.5, ls="--")
        if reattach is not None:
            ax1.axvline(reattach / D, color="r", lw=0.8, ls=":")
        ax1.set_ylabel("u_y / U_mean")
        ax1.legend()
        ax1.set_title(f"75 % cosine stenosis, Re = {RE:.0f}")
        ax2.plot(y / D, p * 1060, "k-")
        ax2.set_xlabel("y / D (throat at 0)")
        ax2.set_ylabel("centreline p [Pa]")
        ax2.axvspan(-1, 1, color="0.9")
        ax1.axvspan(-1, 1, color="0.9")
        fig.tight_layout()
        fig.savefig(case / "stenosis_results.png", dpi=150)
        print("Figure written               : stenosis_results.png")
    except ImportError:
        pass

    checks = {"mass conservation (< 0.1 %)": mass_err < 1e-3,
              "recirculation reattaches in domain": reattach is not None
              and reattach < 15.0 * D}
    if newtonian:
        checks.update({
            "fully developed inlet (< 2 %)": dev_in < 0.02,
            "upstream Poiseuille gradient (< 3 %)": grad_err < 0.03,
            "downstream recovery (< 3 %)": dev_out < 0.03,
        })
    for name, ok in checks.items():
        print(f"  [{'PASS' if ok else 'FAIL'}] {name}")
    ok = all(checks.values())
    print("RESULT                       :", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
