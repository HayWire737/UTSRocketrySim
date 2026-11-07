"""Cross-checks for the N2O property and tank models.

Philosophy: every model is checked against something independent of it —
another correlation, a published table, or a conservation law evaluated numerically.
"""

import numpy as np
import pytest
import CoolProp.CoolProp as CP

from utsrocketrysim.n2o import (
    CoolPropN2O,
    ESDUN2O,
    T_CRIT,
    tank_state,
    mass_for_fill,
    liquid_full_temperature,
    isothermal_vent_rate,
)
from utsrocketrysim.n2o.tank import PHASE_LIQUID_FULL, PHASE_SUPERCRITICAL, PHASE_TWO_PHASE

C = 273.15
cool, esdu = CoolPropN2O(), ESDUN2O()


# --- Saturation properties ------------------------------------------------------------------

# NIST Webbook saturation points (as tabulated in HalfCatSim v1.3.11, converted from F / kPa).
# NOTE: same underlying EOS as CoolProp, so this checks our wrapper/units, not the physics.
NIST_POINTS = [
    # T [F], p [kPa], rho_l, rho_v [kg/m^3]
    (-4, 1797.1, 998.07, 46.727),
    (32, 3122.1, 907.07, 85.171),
    (68, 5052.5, 785.10, 157.99),
    (86, 6307.5, 688.11, 232.70),
]


@pytest.mark.parametrize("T_F,p_kPa,rho_l,rho_v", NIST_POINTS)
def test_coolprop_matches_nist_table(T_F, p_kPa, rho_l, rho_v):
    T = (T_F - 32) * 5 / 9 + C
    s = cool.saturation(T)
    assert s.p == pytest.approx(p_kPa * 1e3, rel=1e-3)
    assert s.rho_l == pytest.approx(rho_l, rel=1e-3)
    assert s.rho_v == pytest.approx(rho_v, rel=1e-3)


@pytest.mark.parametrize("T_C", [-40, -20, 0, 10, 20, 25])
def test_esdu_agrees_with_eos_away_from_critical(T_C):
    """Independent correlation vs EOS. Tolerances are a few x the observed scatter."""
    a, b = cool.saturation(T_C + C), esdu.saturation(T_C + C)
    assert b.p == pytest.approx(a.p, rel=0.01)
    assert b.rho_l == pytest.approx(a.rho_l, rel=0.01)
    assert b.rho_v == pytest.approx(a.rho_v, rel=0.02)
    assert b.h_fg == pytest.approx(a.h_fg, rel=0.02)


def test_clausius_clapeyron_consistency():
    """dp/dT along the saturation line must equal h_fg / (T (v_v - v_l)) — exact thermodynamics.
    Checks the EOS backend's p, rho and h are mutually consistent."""
    T, dT = 290.0, 0.01
    s = cool.saturation(T)
    dpdT_numeric = (cool.saturation(T + dT).p - cool.saturation(T - dT).p) / (2 * dT)
    dpdT_cc = s.h_fg / (T * (1 / s.rho_v - 1 / s.rho_l))
    assert dpdT_numeric == pytest.approx(dpdT_cc, rel=1e-4)


def test_no_saturation_above_critical():
    with pytest.raises(ValueError):
        cool.saturation(T_CRIT + 1)


# --- Tank state ------------------------------------------------------------------------------

def test_tank_mass_and_volume_close():
    V = 0.010
    m = mass_for_fill(V, 15 + C, ullage=0.15)
    st = tank_state(V, m, 15 + C)
    assert st.phase == PHASE_TWO_PHASE
    assert st.m_l + st.m_v == pytest.approx(m)
    assert st.m_l / st.rho_l + st.m_v / st.rho_v == pytest.approx(V, rel=1e-9)
    assert st.ullage_fraction == pytest.approx(0.15, rel=1e-6)


def test_heating_shrinks_ullage_until_liquid_full():
    V = 0.010
    m = mass_for_fill(V, 10 + C, ullage=0.10)
    T_full = liquid_full_temperature(V, m)
    assert T_full is not None and 10 + C < T_full < T_CRIT
    below = tank_state(V, m, T_full - 0.5)
    above = tank_state(V, m, T_full + 0.5)
    assert below.phase == PHASE_TWO_PHASE and below.ullage_fraction < 0.02
    assert above.phase == PHASE_LIQUID_FULL
    # Once liquid-full, dp/dT jumps from the saturation-curve slope (~1.2 bar/K near 20 C) to the
    # liquid's (dp/dT)_rho = beta/kappa (~4.6 bar/K). Near-critical N2O liquid is fairly
    # compressible, so it's ~4x steeper — not the near-instant spike you'd see with water.
    dpdT_sat = (cool.saturation(T_full + 2).p - cool.saturation(T_full).p) / 2
    dpdT_full = (tank_state(V, m, T_full + 3).p - tank_state(V, m, T_full + 1).p) / 2
    assert dpdT_full > 3 * dpdT_sat


def test_supercritical_flagged():
    V = 0.010
    m = mass_for_fill(V, 20 + C, ullage=0.6)  # low fill density -> reaches T_c with ullage left
    assert liquid_full_temperature(V, m) is None
    assert tank_state(V, m, 40 + C).phase == PHASE_SUPERCRITICAL


# --- Venting -----------------------------------------------------------------------------------

def test_isothermal_vent_rate_holds_temperature():
    """Integrate one small step with the analytic vent rate using a direct (rho, u) EOS flash.
    Temperature must stay put; without venting it must rise."""
    V, T0, Q, dt = 0.010, 25 + C, 200.0, 1.0
    m0 = mass_for_fill(V, T0, ullage=0.3)
    st = tank_state(V, m0, T0)
    s = cool.saturation(T0)
    U0 = st.m_l * s.u_l + st.m_v * s.u_v

    def T_after(dm_vent):
        m1 = m0 - dm_vent
        U1 = U0 + Q * dt - s.h_v * dm_vent
        return CP.PropsSI("T", "Dmass", m1 / V, "Umass", U1 / m1, "NitrousOxide")

    dT_vent = T_after(isothermal_vent_rate(T0, Q) * dt) - T0
    dT_sealed = T_after(0.0) - T0
    assert dT_sealed > 1e-3
    assert abs(dT_vent) < 1e-3 * dT_sealed
