#!/usr/bin/env python3
"""Compare the bloodParticles OpenFOAM result with the exact particle motion.

Carrier flow: Poiseuille flow u_y(r) = 2 Um (1 - r^2/R^2) (frozen). Each
particle obeys

    dv/dt = (u(x) - v) / tau_p * (1 + Re_p^(2/3) / 6) + g (1 - rho_f / rho_p),
    tau_p = rho_p d^2 / (18 mu),   Re_p = rho_f |u - v| d / mu,

i.e. the Schiller-Naumann sphere drag and buoyancy-corrected gravity used by
OpenFOAM's sphereDrag and gravity models. These equations are integrated
here with a fine fourth-order Runge-Kutta scheme.

Checks (all must pass):
  1. Tracers (d = 1 um) leave through OUTLET at the time L / u(r0)
     predicted for their release radius (within 2 %), at the release radius.
  2. Debris particles (d = 200 um) reach their terminal settling velocity
     (within 1 %) and deposit on the lower wall at the predicted axial
     position (within 2 % of the travelled distance).

Only numpy is required (matplotlib is optional and used for the figure).
"""
import math
import re
import sys
from pathlib import Path

import numpy as np

R = 0.002
L = 0.06
UM = 0.05
RHO_F = 1060.0
MU = 3.498e-3
RHO_P = 2000.0
G = np.array([0.0, 0.0, -9.81])
Y0 = 0.001
TRACERS = [(0.0, 0.0), (0.0004, 0.0), (0.0008, 0.0), (0.0012, 0.0), (0.0016, 0.0)]  # (x, z)
DEBRIS = [(0.0, 0.0016), (0.0, 0.0008), (0.0, 0.0)]
D_TRACER = 1e-6
D_DEBRIS = 2e-4


def fluid_velocity(p):
    r2 = p[0] ** 2 + p[2] ** 2
    return np.array([0.0, 2.0 * UM * max(0.0, 1.0 - r2 / R**2), 0.0])


def accel(p, v, d):
    tau = RHO_P * d**2 / (18.0 * MU)
    rel = fluid_velocity(p) - v
    re = RHO_F * np.linalg.norm(rel) * d / MU
    return rel / tau * (1.0 + re ** (2.0 / 3.0) / 6.0) + G * (1.0 - RHO_F / RHO_P)


def track(x0, z0, d, t_end=3.0):
    """Integrate a particle from rest; return (event, time, position, velocity)."""
    tau = RHO_P * d**2 / (18.0 * MU)
    dt = min(2e-5, tau / 5.0) if tau > 1e-6 else None
    p = np.array([x0, Y0, z0])
    v = np.zeros(3)
    t = 0.0
    if dt is None:
        # Tracer limit: relaxation time ~1e-8 s, the particle moves with the
        # fluid plus its (negligible) settling velocity
        vs = G * (1.0 - RHO_F / RHO_P) * tau
        u = fluid_velocity(p)[1]
        t_out = (L - Y0) / u
        p_out = p + np.array([0.0, L - Y0, 0.0]) + vs * t_out
        return "OUTLET", t_out, p_out, fluid_velocity(p) + vs
    while t < t_end:
        k1v = accel(p, v, d); k1p = v
        k2v = accel(p + 0.5 * dt * k1p, v + 0.5 * dt * k1v, d); k2p = v + 0.5 * dt * k1v
        k3v = accel(p + 0.5 * dt * k2p, v + 0.5 * dt * k2v, d); k3p = v + 0.5 * dt * k2v
        k4v = accel(p + dt * k3p, v + dt * k3v, d); k4p = v + dt * k3v
        p_new = p + dt / 6.0 * (k1p + 2 * k2p + 2 * k3p + k4p)
        v_new = v + dt / 6.0 * (k1v + 2 * k2v + 2 * k3v + k4v)
        if p_new[0] ** 2 + p_new[2] ** 2 >= R**2:
            # linear interpolation to the wall crossing
            f0 = math.hypot(p[0], p[2]) - R
            f1 = math.hypot(p_new[0], p_new[2]) - R
            s = f0 / (f0 - f1)
            return "WALL", t + s * dt, p + s * (p_new - p), v + s * (v_new - v)
        if p_new[1] >= L:
            s = (L - p[1]) / (p_new[1] - p[1])
            return "OUTLET", t + s * dt, p + s * (p_new - p), v + s * (v_new - v)
        p, v, t = p_new, v_new, t + dt
    return "NONE", t, p, v


def terminal_velocity(d):
    """Settling speed in still fluid with Schiller-Naumann drag (bisection)."""
    tau = RHO_P * d**2 / (18.0 * MU)
    g = 9.81 * (1.0 - RHO_F / RHO_P)
    lo, hi = 0.0, g * tau
    for _ in range(200):
        w = 0.5 * (lo + hi)
        re = RHO_F * w * d / MU
        if w / tau * (1.0 + re ** (2.0 / 3.0) / 6.0) > g:
            hi = w
        else:
            lo = w
    return 0.5 * (lo + hi)


def read_hits(case):
    """All hit records from particlePostProcessing, keyed by (typeId?, origId)."""
    base = case / "postProcessing" / "lagrangian" / "kinematicCloud" / "wallAndOutlet"
    hits = []
    for f in base.glob("*/*.dat"):
        patch = f.stem
        for line in f.read_text().splitlines():
            if not line.strip() or line.startswith("#"):
                continue
            nums = [float(x) for x in re.findall(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", line)]
            # Time currentProc (4 coords) (3 position) celli tetFacei tetPti
            # facei stepFraction origProc origId active typeId nParticle d
            # dTarget (3 U) rho age ...
            t = nums[0]
            pos = np.array(nums[6:9])
            orig_id = int(nums[15])
            d = nums[19]
            vel = np.array(nums[21:24])
            hits.append({"patch": patch, "t": t, "pos": pos, "id": orig_id,
                         "d": d, "U": vel})
    # The same hit can be repeated in later write intervals: keep one per particle
    unique = {}
    for h in hits:
        unique[(round(h["d"] * 1e9), h["id"])] = h
    return list(unique.values())


def main():
    case = Path(__file__).resolve().parent
    if not (case / "postProcessing" / "lagrangian").is_dir():
        sys.exit("No postProcessing data found: run ./Allrun first.")
    hits = read_hits(case)
    tracers = sorted([h for h in hits if h["d"] < 1e-5], key=lambda h: h["id"])
    debris = sorted([h for h in hits if h["d"] > 1e-5], key=lambda h: h["id"])

    ok_all = True
    print(f"Carrier flow: Poiseuille, Um = {UM} m/s, Re = {2*R*UM*RHO_F/MU:.0f}")
    print(f"\nTracers (d = 1 um), {len(tracers)}/{len(TRACERS)} reached the outlet")
    print("  r0/R   t_out OpenFOAM   t_out exact   error   r_out/R")
    worst_t = 0.0
    for i, (x0, z0) in enumerate(TRACERS):
        match = [h for h in tracers if abs(h["pos"][0] - x0) < 1e-4 and h["patch"] == "OUTLET"]
        _, t_ex, _, _ = track(x0, z0, D_TRACER)
        if not match:
            print(f"  {x0/R:4.2f}   missing")
            ok_all = False
            continue
        h = match[0]
        err = abs(h["t"] - t_ex) / t_ex
        worst_t = max(worst_t, err)
        print(f"  {x0/R:4.2f}   {h['t']:12.4f}   {t_ex:11.4f}   {100*err:5.2f} %   "
              f"{math.hypot(h['pos'][0], h['pos'][2])/R:.3f}")

    w_t = terminal_velocity(D_DEBRIS)
    print(f"\nDebris (d = 200 um), {len(debris)}/{len(DEBRIS)} deposited on the wall")
    print(f"  Terminal settling velocity (Schiller-Naumann): {w_t*1e3:.4f} mm/s "
          f"(Stokes: {(RHO_P-RHO_F)*9.81*D_DEBRIS**2/(18*MU)*1e3:.4f} mm/s)")
    print("  z0/R   y_dep OpenFOAM   y_dep exact   error   w at wall")
    worst_y, worst_w = 0.0, 0.0
    for x0, z0 in DEBRIS:
        match = [h for h in debris if h["patch"] == "WALL"
                 and abs(h["t"] - track(x0, z0, D_DEBRIS)[1]) < 0.1]
        ev, t_ex, p_ex, v_ex = track(x0, z0, D_DEBRIS)
        if not match or ev != "WALL":
            print(f"  {z0/R:4.2f}   missing (exact: {ev})")
            ok_all = False
            continue
        h = min(match, key=lambda m: abs(m["pos"][1] - p_ex[1]))
        err = abs(h["pos"][1] - p_ex[1]) / (p_ex[1] - Y0)
        werr = abs(-h["U"][2] - w_t) / w_t
        worst_y, worst_w = max(worst_y, err), max(worst_w, werr)
        print(f"  {z0/R:4.2f}   {h['pos'][1]*1e3:9.3f} mm   {p_ex[1]*1e3:8.3f} mm   "
              f"{100*err:5.2f} %   {-h['U'][2]*1e3:.4f} mm/s")

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(8, 3.2))
        for x0, z0 in DEBRIS:
            p = np.array([x0, Y0, z0]); v = np.zeros(3); path = [p.copy()]
            dt = 2e-4
            for _ in range(20000):
                v = v + dt * accel(p, v, D_DEBRIS); p = p + dt * v; path.append(p.copy())
                if p[2] <= -R or p[1] >= L:
                    break
            path = np.array(path)
            ax.plot(path[:, 1] * 1e3, path[:, 2] * 1e3, "k-", lw=1)
        for h in debris:
            ax.plot(h["pos"][1] * 1e3, h["pos"][2] * 1e3, "ro")
        ax.axhline(-R * 1e3, color="0.5"); ax.axhline(R * 1e3, color="0.5")
        ax.set_xlim(0, L * 1e3)
        ax.set_xlabel("y [mm]"); ax.set_ylabel("z [mm]")
        ax.set_title("200 um debris settling in Poiseuille flow (lines: exact, dots: OpenFOAM deposits)")
        fig.tight_layout()
        fig.savefig(case / "particles_validation.png", dpi=150)
        print("\nFigure written: particles_validation.png")
    except ImportError:
        pass

    checks = {
        "tracer transit times (< 2 %)": len(tracers) == len(TRACERS) and worst_t < 0.02,
        "debris settling velocity (< 1 %)": len(debris) == len(DEBRIS) and worst_w < 0.01,
        "debris deposition position (< 2 %)": len(debris) == len(DEBRIS) and worst_y < 0.02,
    }
    print()
    for name, ok in checks.items():
        print(f"  [{'PASS' if ok else 'FAIL'}] {name}")
    ok = ok_all and all(checks.values())
    print("RESULT:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
