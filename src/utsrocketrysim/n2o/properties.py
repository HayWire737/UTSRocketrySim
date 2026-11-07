"""Nitrous oxide thermophysical properties.

Two independent backends with the same interface, so every number can be cross-checked:

* ``CoolPropN2O`` — full Helmholtz-energy equation of state (Lemmon & Span, 2006) via CoolProp.
  Valid for liquid, vapour, two-phase and supercritical states. This is the same EOS behind
  NIST REFPROP / the NIST Webbook (and therefore HalfCatSim's N2O table). Use this one for
  simulation.

* ``ESDUN2O`` — the ESDU 91022 saturation-line correlations (as popularised by Aspire Space's
  "Modelling the nitrous run tank emptying"). Saturation only, valid roughly -90 C to the
  critical point. Use it as an independent sanity check of the EOS backend and of our code.

Units: SI throughout (K, Pa, kg/m^3, J/kg, J/kg/K).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import CoolProp.CoolProp as CP

FLUID = "NitrousOxide"

# Critical point from the Lemmon & Span EOS (CoolProp). ESDU uses slightly different values.
T_CRIT = CP.PropsSI("Tcrit", FLUID)  # ~309.52 K (36.37 C)
P_CRIT = CP.PropsSI("pcrit", FLUID)  # ~7.245 MPa
RHO_CRIT = CP.PropsSI("rhocrit", FLUID)  # ~452 kg/m^3
T_TRIPLE = CP.PropsSI("Ttriple", FLUID)  # ~182.33 K
MOLAR_MASS = CP.PropsSI("molarmass", FLUID)  # 0.0440128 kg/mol


@dataclass(frozen=True)
class SaturationState:
    """Properties of saturated liquid (_l) and saturated vapour (_v) at temperature T.

    Fields a backend cannot provide are ``nan``.
    """

    T: float
    p: float
    rho_l: float
    rho_v: float
    h_l: float
    h_v: float
    cp_l: float = float("nan")
    cp_v: float = float("nan")
    u_l: float = float("nan")
    u_v: float = float("nan")
    s_l: float = float("nan")
    s_v: float = float("nan")

    @property
    def h_fg(self) -> float:
        """Latent heat of vaporisation, J/kg (independent of enthalpy reference state)."""
        return self.h_v - self.h_l


class CoolPropN2O:
    """Helmholtz-EOS backend. The reference backend for simulation."""

    name = "CoolProp (Lemmon & Span 2006 EOS)"

    def saturation(self, T: float) -> SaturationState:
        if not (T_TRIPLE <= T < T_CRIT):
            raise ValueError(
                f"T = {T:.2f} K is outside the two-phase region "
                f"({T_TRIPLE:.2f} K to {T_CRIT:.2f} K). Above T_crit there is no liquid."
            )

        def liq(k):
            return CP.PropsSI(k, "T", T, "Q", 0, FLUID)

        def vap(k):
            return CP.PropsSI(k, "T", T, "Q", 1, FLUID)

        return SaturationState(
            T=T,
            p=liq("P"),
            rho_l=liq("D"),
            rho_v=vap("D"),
            h_l=liq("H"),
            h_v=vap("H"),
            cp_l=liq("C"),
            cp_v=vap("C"),
            u_l=liq("U"),
            u_v=vap("U"),
            s_l=liq("S"),
            s_v=vap("S"),
        )

    def saturation_from_pressure(self, p: float) -> SaturationState:
        T = CP.PropsSI("T", "P", p, "Q", 0, FLUID)
        return self.saturation(T)

    @staticmethod
    def state(input1: str, value1: float, input2: str, value2: float, output: str) -> float:
        """General single-call EOS lookup, e.g. ``state("T", 300, "D", 600, "P")``.

        Thin wrapper over CoolProp.PropsSI so the rest of the code never imports CoolProp directly
        (makes swapping in REFPROP or a tabulated backend later a one-file change).
        """
        return CP.PropsSI(output, input1, value1, input2, value2, FLUID)


class ESDUN2O:
    """ESDU 91022 saturation correlations. Independent cross-check backend.

    Provides p, rho_l, rho_v, h_l, h_v, cp_l. Enthalpies use the ESDU reference state, so only
    differences (e.g. h_fg) are comparable with CoolProp.
    """

    name = "ESDU 91022 correlations"

    # Critical constants used by the correlations themselves.
    Tc = 309.57  # K
    pc = 7.251e6  # Pa
    rhoc = 452.0  # kg/m^3

    def saturation(self, T: float) -> SaturationState:
        if not (183.0 <= T < self.Tc):
            raise ValueError(f"T = {T:.2f} K outside ESDU range (183 K to {self.Tc} K)")
        Tr = T / self.Tc
        a = 1.0 - Tr  # (1 - Tr)
        bt = 1.0 / Tr - 1.0  # (1/Tr - 1)

        # Vapour pressure
        p = self.pc * np.exp(
            (1.0 / Tr) * (-6.71893 * a + 1.35966 * a**1.5 - 1.3779 * a**2.5 - 4.051 * a**5)
        )
        # Saturated liquid density
        rho_l = self.rhoc * np.exp(
            1.72328 * a ** (1 / 3) - 0.8395 * a ** (2 / 3) + 0.5106 * a - 0.10412 * a ** (4 / 3)
        )
        # Saturated vapour density
        rho_v = self.rhoc * np.exp(
            -1.009 * bt ** (1 / 3)
            - 6.28792 * bt ** (2 / 3)
            + 7.50332 * bt
            - 7.90463 * bt ** (4 / 3)
            + 0.629427 * bt ** (5 / 3)
        )
        # Specific enthalpies (kJ/kg in the source -> J/kg)
        h_l = 1e3 * (
            -200.0
            + 116.043 * a ** (1 / 3)
            - 917.225 * a ** (2 / 3)
            + 794.779 * a
            - 589.587 * a ** (4 / 3)
        )
        h_v = 1e3 * (
            -200.0
            + 440.055 * a ** (1 / 3)
            - 459.701 * a ** (2 / 3)
            + 434.081 * a
            - 485.338 * a ** (4 / 3)
        )
        # Saturated liquid isobaric heat capacity (kJ/kg/K -> J/kg/K)
        cp_l = 1e3 * 2.49973 * (1.0 + 0.023454 / a - 3.80136 * a + 13.0945 * a**2 - 14.518 * a**3)

        return SaturationState(
            T=T, p=float(p), rho_l=float(rho_l), rho_v=float(rho_v),
            h_l=float(h_l), h_v=float(h_v), cp_l=float(cp_l),
        )


# Default backend used across the package.
N2O = CoolPropN2O()
