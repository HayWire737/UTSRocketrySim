"""Tests for N2O properties.

Run with:  pytest
Each test checks our code against something independent of it.
"""

import pytest

from utsrocketrysim.n2o.esdu import esdu_saturation
from utsrocketrysim.n2o.properties import (
    T_CRIT,
    saturation,
    saturation_from_pressure,
    state,
)

C = 273.15

# NIST Webbook saturation data (also tabulated in HalfCatSim v1.3.11), converted to SI.
# NIST uses the same equation of state as CoolProp, so these check our wrapper and units,
# not the physics.
NIST = [
    # T [C],  p [Pa],   rho_l,  rho_v [kg/m^3]
    (-20.0, 1.7971e6, 998.07, 46.727),
    (0.0, 3.1221e6, 907.07, 85.171),
    (20.0, 5.0525e6, 785.10, 157.99),
    (30.0, 6.3075e6, 688.11, 232.70),
]


@pytest.mark.parametrize("T_c, p, rho_l, rho_v", NIST)
def test_matches_nist(T_c, p, rho_l, rho_v):
    s = saturation(T_c + C)
    assert s.p == pytest.approx(p, rel=1e-3)
    assert s.rho_l == pytest.approx(rho_l, rel=1e-3)
    assert s.rho_v == pytest.approx(rho_v, rel=1e-3)


@pytest.mark.parametrize("T_c", [-40, -20, 0, 10, 20, 25])
def test_esdu_agrees_with_coolprop(T_c):
    """The real physics check: an independent correlation agrees with the EOS."""
    e, c = esdu_saturation(T_c + C), saturation(T_c + C)
    assert e.p == pytest.approx(c.p, rel=0.01)
    assert e.rho_l == pytest.approx(c.rho_l, rel=0.01)
    assert e.rho_v == pytest.approx(c.rho_v, rel=0.02)
    assert e.h_fg == pytest.approx(c.h_fg, rel=0.02)


def test_clausius_clapeyron():
    """Exact thermodynamics: along the saturation line, dp/dT = h_fg / (T * (v_v - v_l)).
    If pressure, density and enthalpy are consistent with each other, both sides agree."""
    T, dT = 290.0, 0.01
    s = saturation(T)
    dpdT_numerical = (saturation(T + dT).p - saturation(T - dT).p) / (2 * dT)
    dpdT_theory = s.h_fg / (T * (1 / s.rho_v - 1 / s.rho_l))
    assert dpdT_numerical == pytest.approx(dpdT_theory, rel=1e-4)


def test_pressure_round_trip():
    s = saturation(20 + C)
    assert saturation_from_pressure(s.p).T == pytest.approx(20 + C, abs=1e-6)


def test_liquid_full_pressure_exceeds_saturation():
    """Liquid squeezed into a full tank sits well above the saturation pressure."""
    p_full = state("P", "T", 30 + C, "D", 750)
    assert p_full > 1.2 * saturation(30 + C).p


def test_celsius_by_mistake_gives_helpful_error():
    with pytest.raises(ValueError, match="Celsius"):
        saturation(20)


def test_above_critical_temperature_rejected():
    with pytest.raises(ValueError, match="critical point"):
        saturation(T_CRIT + 1)


def test_above_critical_pressure_rejected():
    with pytest.raises(ValueError, match="critical pressure"):
        saturation_from_pressure(80e5)