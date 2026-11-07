"""A sealed N2O tank warming on the pad after filling.

Run:  python examples/tank_heating.py
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from utsrocketrysim.n2o.properties import T_CRIT
from utsrocketrysim.n2o.tank import liquid_full_temperature, mass_for_fill, tank_state

C = 273.15
V = 0.008             # tank volume, m^3 (8 L): change to your tank
T_FILL = 15 + C       # temperature at fill, K
ULLAGES = [0.05, 0.10, 0.15, 0.20]
temps = np.linspace(0 + C, 45 + C, 451)

fig, (ax_p, ax_u) = plt.subplots(2, 1, figsize=(8, 7), sharex=True, constrained_layout=True)

for ullage in ULLAGES:
    m = mass_for_fill(V, T_FILL, ullage)
    states = [tank_state(V, m, T) for T in temps]
    label = f"{ullage:.0%} ullage at fill ({m:.2f} kg)"
    line, = ax_p.plot(temps - C, [s.p / 1e5 for s in states], lw=2, label=label)
    ax_u.plot(temps - C, [s.ullage * 100 for s in states], lw=2, color=line.get_color())

    T_full = liquid_full_temperature(V, m)
    if T_full is not None:
        p_full = tank_state(V, m, T_full).p
        ax_p.plot(T_full - C, p_full / 1e5, "o", color=line.get_color(), ms=8)
        ax_u.plot(T_full - C, 0, "o", color=line.get_color(), ms=8)
        print(f"{label}: liquid-full at {T_full - C:.1f} C ({p_full / 1e5:.1f} bar)")
    else:
        print(f"{label}: goes supercritical at {T_CRIT - C:.1f} C")

for ax in (ax_p, ax_u):
    ax.axvspan(T_CRIT - C, 45, color="0.92")
    ax.axvline(T_CRIT - C, color="0.5", ls="--", lw=1)
    ax.grid(alpha=0.3)

ax_p.set_ylim(25, 115)
ax_p.set_ylabel("Tank pressure (bar)")
ax_p.set_title(f"Sealed {V * 1e3:.0f} L N2O tank filled at {T_FILL - C:.0f} °C  "
               "(● = liquid-full; shaded = supercritical)", loc="left")
ax_p.legend(loc="upper left", frameon=False)
ax_u.set_ylabel("Vapour ullage (% of volume)")
ax_u.set_xlabel("Tank temperature (°C)")
ax_u.set_xlim(0, 45)
ax_u.set_ylim(-1, 30)

out = Path("output")
out.mkdir(exist_ok=True)
plt.savefig(out / "tank_heating.png", dpi=150)
plt.show()