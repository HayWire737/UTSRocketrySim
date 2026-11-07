"""N2O run tank: the state of a fixed mass of nitrous in a fixed volume.

Model: thermal equilibrium. Liquid and vapour share one temperature, and whenever both exist they
sit on the saturation line. Standard starting point for N2O tanks (Whitmore & Chandler 2010;
Zimmerman et al. 2013).

Key idea: volume V and mass m are fixed, so the average density m/V is fixed. Temperature then
decides how that mass splits between liquid and vapour, and that split sets the pressure.
"""

from dataclasses import dataclass

from scipy.optimize import brentq

from utsrocketrysim.n2o.properties import T_CRIT, T_TRIPLE, saturation, state

TWO_PHASE = "two-phase"
LIQUID_FULL = "liquid-full"
VAPOUR_ONLY = "vapour-only"
SUPERCRITICAL = "supercritical"


@dataclass(frozen=True)
class TankState:
    """N2O in a tank at one instant. SI units."""

    T: float        # temperature, K
    p: float        # pressure, Pa
    V: float        # tank volume, m^3
    m: float        # total N2O mass, kg
    m_l: float      # liquid mass, kg
    m_v: float      # vapour mass, kg
    phase: str      # one of the four phase names above

    @property
    def quality(self) -> float:
        """Vapour mass fraction, x = m_v / m."""
        return self.m_v / self.m

    @property
    def ullage(self) -> float:
        """Fraction of the tank VOLUME that is vapour. nan when supercritical (no liquid/vapour)."""
        if self.phase == TWO_PHASE:
            V_liquid = self.m_l / saturation(self.T).rho_l
            return 1.0 - V_liquid / self.V
        if self.phase == LIQUID_FULL:
            return 0.0
        if self.phase == VAPOUR_ONLY:
            return 1.0
        return float("nan")


def tank_state(V: float, m: float, T: float) -> TankState:
    """State of m kg of N2O in a V m^3 tank at temperature T (K)."""
    rho = m / V  # average density, fixed by the tank

    if T >= T_CRIT:
        p = state("P", "T", T, "D", rho)
        return TankState(T, p, V, m, 0.0, 0.0, SUPERCRITICAL)

    s = saturation(T)
    # Quality from specific volumes (v = 1/rho):  v = (1 - x) v_l + x v_v
    v, v_l, v_v = 1 / rho, 1 / s.rho_l, 1 / s.rho_v
    x = (v - v_l) / (v_v - v_l)

    if x <= 0.0:
        # The liquid has expanded to fill the whole tank. Now it's compressed liquid, and the
        # pressure comes from the equation of state, not the saturation line.
        p = state("P", "T", T, "D", rho)
        return TankState(T, p, V, m, m, 0.0, LIQUID_FULL)
    if x >= 1.0:
        p = state("P", "T", T, "D", rho)
        return TankState(T, p, V, m, 0.0, m, VAPOUR_ONLY)

    return TankState(T, s.p, V, m, (1 - x) * m, x * m, TWO_PHASE)


def mass_for_fill(V: float, T_fill: float, ullage: float) -> float:
    """Mass of N2O (kg) when a V m^3 tank is filled at T_fill (K), leaving `ullage` of the
    volume as vapour (e.g. 0.15 = 15%)."""
    s = saturation(T_fill)
    return V * ((1 - ullage) * s.rho_l + ullage * s.rho_v)


def liquid_full_temperature(V: float, m: float) -> float | None:
    """Temperature (K) at which the expanding liquid fills the tank completely.

    Returns None if the tank reaches the critical point first (low fill), meaning it goes
    supercritical instead of liquid-full.
    """
    rho = m / V

    def gap(T):
        return saturation(T).rho_l - rho  # positive while there is still room for vapour

    T_low, T_high = T_TRIPLE + 1.0, T_CRIT - 1e-3
    if gap(T_high) > 0:
        return None
    if gap(T_low) < 0:
        return T_low
    return brentq(gap, T_low, T_high, xtol=1e-6)


def isothermal_vent_rate(T: float, Q_dot: float) -> float:
    """Vapour vent rate (kg/s) needed to hold a two-phase tank at constant T while heat Q_dot (W)
    leaks in.

        m_dot = Q_dot * (1 - rho_v / rho_l) / h_fg

    Derivation in docs/physics/n2o.md.
    """
    s = saturation(T)
    return Q_dot * (1 - s.rho_v / s.rho_l) / s.h_fg