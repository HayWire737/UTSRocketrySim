"""Thermodynamic state of a closed (or venting) N2O tank.

Phase 1 model: *homogeneous thermal equilibrium* — liquid and vapour in the tank share one
temperature and sit on the saturation line whenever both phases exist. This is the standard
"equilibrium tank" model (Whitmore & Chandler 2010, Zimmerman et al. 2013). Its known weakness
is that it under-predicts the pressure drop during fast blowdown (real tanks go out of
equilibrium — the liquid superheats). A non-equilibrium model is planned for phase 2.

Given the tank volume V and loaded mass m, the specific volume v = V/m is fixed. Temperature then
decides everything:

    quality  x = (v - v_l) / (v_v - v_l)      (vapour mass fraction)

* 0 < x < 1  -> two-phase, p = p_sat(T)
* x <= 0     -> the liquid has expanded to fill the whole tank ("liquid-full" / hydraulic lock).
               Pressure is then set by compressing liquid: dp/dT jumps ~4x (about 1.2 -> 4.6
               bar/K near 20 C). Your burst disc / relief must be sized for this case.
* T >= T_c   -> supercritical: no liquid/vapour distinction, no liquid to feed the injector.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .properties import N2O, T_CRIT, CoolPropN2O

PHASE_TWO_PHASE = "two-phase"
PHASE_LIQUID_FULL = "liquid-full"
PHASE_VAPOUR_ONLY = "vapour-only"
PHASE_SUPERCRITICAL = "supercritical"


@dataclass(frozen=True)
class TankState:
    T: float  # K
    p: float  # Pa
    m: float  # kg, total N2O in tank
    V: float  # m^3
    phase: str
    m_l: float  # kg liquid
    m_v: float  # kg vapour
    rho_l: float  # kg/m^3 (saturated liquid; nan if not two-phase)
    rho_v: float  # kg/m^3 (saturated vapour; nan if not two-phase)

    @property
    def quality(self) -> float:
        """Vapour mass fraction x = m_v / m."""
        return self.m_v / self.m

    @property
    def liquid_volume_fraction(self) -> float:
        """Fraction of the tank volume occupied by liquid (1 - ullage)."""
        if self.phase == PHASE_TWO_PHASE:
            return (self.m_l / self.rho_l) / self.V
        if self.phase == PHASE_LIQUID_FULL:
            return 1.0
        if self.phase == PHASE_VAPOUR_ONLY:
            return 0.0
        return float("nan")  # supercritical: liquid/vapour split is undefined

    @property
    def ullage_fraction(self) -> float:
        return 1.0 - self.liquid_volume_fraction


def mass_for_fill(V: float, T_fill: float, ullage: float, props: CoolPropN2O = N2O) -> float:
    """Mass of N2O in a tank of volume V filled at T_fill with a given vapour ullage fraction."""
    sat = props.saturation(T_fill)
    return V * ((1.0 - ullage) * sat.rho_l + ullage * sat.rho_v)


def tank_state(V: float, m: float, T: float, props: CoolPropN2O = N2O) -> TankState:
    """Equilibrium tank state for fixed volume V, mass m, at temperature T."""
    rho = m / V
    if T >= T_CRIT:
        p = props.state("T", T, "D", rho, "P")
        return TankState(T, p, m, V, PHASE_SUPERCRITICAL, 0.0, m, np.nan, np.nan)

    sat = props.saturation(T)
    v, v_l, v_v = 1.0 / rho, 1.0 / sat.rho_l, 1.0 / sat.rho_v
    x = (v - v_l) / (v_v - v_l)

    if x <= 0.0:
        # Compressed liquid filling the whole tank. EOS gives the (very high) pressure.
        p = props.state("T", T, "D", rho, "P")
        return TankState(T, p, m, V, PHASE_LIQUID_FULL, m, 0.0, sat.rho_l, sat.rho_v)
    if x >= 1.0:
        p = props.state("T", T, "D", rho, "P")
        return TankState(T, p, m, V, PHASE_VAPOUR_ONLY, 0.0, m, sat.rho_l, sat.rho_v)

    return TankState(T, sat.p, m, V, PHASE_TWO_PHASE, (1 - x) * m, x * m, sat.rho_l, sat.rho_v)


def liquid_full_temperature(V: float, m: float, props: CoolPropN2O = N2O) -> float | None:
    """Temperature at which liquid expansion fills the tank (ullage -> 0), or None if the tank
    reaches the critical temperature first (low fill density)."""
    from scipy.optimize import brentq

    rho = m / V

    def f(T):
        return props.saturation(T).rho_l - rho

    T_lo, T_hi = 200.0, T_CRIT - 1e-3
    if f(T_hi) > 0:  # liquid at T_c is still denser than the mean fill density
        return None
    if f(T_lo) < 0:
        return T_lo
    return brentq(f, T_lo, T_hi, xtol=1e-6)


def isothermal_vent_rate(T: float, Q_dot: float, props: CoolPropN2O = N2O) -> float:
    """Mass of vapour that must be vented (kg/s) to hold a two-phase tank at constant T (and so
    constant p) while heat Q_dot (W) leaks in.

    Derivation (fixed volume, vent vapour only, both phases stay saturated at T):
        volume:  dm_l/rho_l + dm_v/rho_v = 0
        energy:  dU = Q - h_v * dm_vent,   dm_vent = -(dm_l + dm_v)
    Substituting and using h = u + p/rho gives  dm_l = -Q / h_fg, so

        m_dot_vent = Q_dot * (1 - rho_v/rho_l) / h_fg

    The (1 - rho_v/rho_l) factor is because the evaporated liquid's old volume is refilled with
    vapour that stays behind. See docs/physics/n2o.md for the full working.
    """
    sat = props.saturation(T)
    return Q_dot * (1.0 - sat.rho_v / sat.rho_l) / sat.h_fg
