"""What happens to a sealed N2O run tank as it heats up on the pad?

Fill a tank at T_fill with some vapour ullage, seal it, then let it warm. For each initial ullage we
plot tank pressure and remaining ullage vs. temperature, and mark:
  * the liquid-full point (ullage -> 0, pressure slope jumps ~4x), and
  * the critical temperature (36.4 C: above it there is no liquid N2O at all).

Also prints the vent rate needed to hold the tank at temperature against a given heat leak.

Run:  python examples/01_desert_tank_heating.py
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from utsrocketrysim.n2o import (
    T_CRIT, isothermal_vent_rate, liquid_full_temperature, mass_for_fill, tank_state,
)

C = 273.15
V_TANK = 0.008          # m^3 (8 L) — placeholder, set to your run tank
T_FILL = 15 + C         # K, temperature at fill
ULLAGES = [0.05, 0.10, 0.15, 0.20]
COLORS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]  # validated categorical slots 1-4
T_SWEEP = np.linspace(0 + C, 45 + C, 451)

OUT = Path(__file__).resolve().parent.parent / "output"
OUT.mkdir(exist_ok=True)

fig, (ax_p, ax_u) = plt.subplots(2, 1, figsize=(8, 7.5), sharex=True, constrained_layout=True)
ink, muted, grid = "#1a1a19", "#6b6a63", "#e4e3dc"

print(f"Tank {V_TANK*1e3:.1f} L filled at {T_FILL-C:.0f} C")
for ull, col in zip(ULLAGES, COLORS):
    m = mass_for_fill(V_TANK, T_FILL, ull)
    states = [tank_state(V_TANK, m, T) for T in T_SWEEP]
    p = np.array([s.p for s in states]) / 1e5
    u = np.array([s.ullage_fraction for s in states]) * 100
    Tc = T_SWEEP - C
    label = f"{ull*100:.0f}% ullage ({m:.2f} kg)"
    ax_p.plot(Tc, p, color=col, lw=2, label=label)
    ax_u.plot(Tc, u, color=col, lw=2, label=label)

    T_full = liquid_full_temperature(V_TANK, m)
    if T_full is not None:
        st = tank_state(V_TANK, m, T_full)
        ax_p.plot(T_full - C, st.p / 1e5, "o", ms=8, color=col, mec="white", mew=2, zorder=5)
        ax_u.plot(T_full - C, 0, "o", ms=8, color=col, mec="white", mew=2, zorder=5)
        print(f"  {label:>24}: liquid-full at {T_full-C:5.1f} C, p = {st.p/1e5:5.1f} bar")
    else:
        print(f"  {label:>24}: reaches T_crit with vapour space left (goes supercritical)")

for ax in (ax_p, ax_u):
    ax.axvspan(T_CRIT - C, 45, color="#efeee8", zorder=0)
    ax.axvline(T_CRIT - C, color=muted, lw=1, ls="--")
    ax.grid(True, color=grid, lw=0.8)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.tick_params(colors=muted)

ax_u.text(T_CRIT - C + 0.5, 25, "above 36.4 °C:\nsupercritical,\nno liquid N2O", color=muted,
          fontsize=9, va="center")
ax_p.set_ylim(25, 115)
ax_p.set_ylabel("Tank pressure (bar)", color=ink)
ax_p.set_title(
    f"Sealed N2O tank warming after fill — {V_TANK*1e3:.0f} L, filled at {T_FILL-C:.0f} °C\n"
    "● = liquid-full: ullage gone, pressure rises ~4× faster",
    loc="left", color=ink, fontsize=11,
)
ax_p.legend(frameon=False, loc="upper left")
ax_u.set_ylabel("Vapour ullage (% of tank volume)", color=ink)
ax_u.set_xlabel("Tank temperature (°C)", color=ink)
ax_u.set_ylim(-1, 50)
ax_u.set_xlim(0, 45)

path = OUT / "01_desert_tank_heating.png"
fig.savefig(path, dpi=150)
print(f"\nSaved {path}")

print("\nVent rate to hold temperature against a heat leak (two-phase tank):")
for TC in (15, 25, 30, 34):
    for Q in (100, 500):
        r = isothermal_vent_rate(TC + C, Q)
        print(f"  {TC:2d} C, {Q:3d} W  ->  {r*1e3:5.2f} g/s  ({r*3600:5.2f} kg/h)")
