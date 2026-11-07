"""Tests for the equilibrium tank model."""

import CoolProp.CoolProp as CP
import pytest

from utsrocketrysim.n2o.properties import T_CRIT, saturation
from utsrocketrysim.n2o.tank import (
    LIQUID_FULL, SUPERCRITICAL, TWO_PHASE,
    isothermal_vent_rate, liquid_full_temperature, mass_for_fill, tank_state,
)

C = 273.15
V = 0.010  # 10 L test tank


def test_fill_gives_requested_ullage():
    m = mass_for_fill(V, 15 + C, ullage=0.15)
    t = tank_state(V, m, 15 + C)
    assert t.phase == TWO_PHASE
    assert t.ullage == pytest.approx(0.15, rel=1e-6)


def test_mass_and_volume_add_up():
    """Conservation: liquid + vapour masses = total, liquid + vapour volumes = tank."""
    m = mass_for_fill(V, 15 + C, ullage=0.3)
    t = tank_state(V, m, 25 + C)
    s = saturation(25 + C)
    assert t.m_l + t.m_v == pytest.approx(m)
    assert t.m_l / s.rho_l + t.m_v / s.rho_v == pytest.approx(V, rel=1e-9)


def test_heating_shrinks_ullage():
    m = mass_for_fill(V, 10 + C, ullage=0.2)
    ullages = [tank_state(V, m, T + C).ullage for T in (10, 15, 20, 25)]
    assert ullages == sorted(ullages, reverse=True)


def test_liquid_full_point():
    m = mass_for_fill(V, 10 + C, ullage=0.10)
    T_full = liquid_full_temperature(V, m)
    assert tank_state(V, m, T_full - 0.5).phase == TWO_PHASE
    assert tank_state(V, m, T_full + 0.5).phase == LIQUID_FULL
    # Past the liquid-full point, pressure rises several times faster per degree
    slope_sat = (saturation(T_full + 2).p - saturation(T_full).p) / 2
    slope_full = (tank_state(V, m, T_full + 3).p - tank_state(V, m, T_full + 1).p) / 2
    assert slope_full > 3 * slope_sat


def test_low_fill_goes_supercritical():
    m = mass_for_fill(V, 20 + C, ullage=0.6)
    assert liquid_full_temperature(V, m) is None
    assert tank_state(V, m, 40 + C).phase == SUPERCRITICAL


def test_vent_rate_holds_temperature():
    """Energy balance check: take one 1 s step with the vent rate and find the new temperature
    straight from the equation of state. It must not move. Without venting it must rise."""
    T0, Q, dt = 25 + C, 200.0, 1.0
    m0 = mass_for_fill(V, T0, ullage=0.3)
    t, s = tank_state(V, m0, T0), saturation(T0)
    U0 = t.m_l * s.u_l + t.m_v * s.u_v  # total internal energy, J

    def T_after(m_vented):
        m1 = m0 - m_vented
        U1 = U0 + Q * dt - s.h_v * m_vented  # heat in, minus energy carried out by the vapour
        return CP.PropsSI("T", "Dmass", m1 / V, "Umass", U1 / m1, "NitrousOxide")

    dT_vented = T_after(isothermal_vent_rate(T0, Q) * dt) - T0
    dT_sealed = T_after(0.0) - T0
    assert dT_sealed > 1e-3
    assert abs(dT_vented) < 1e-3 * dT_sealed