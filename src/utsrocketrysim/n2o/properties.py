"""Nitrous oxide (N2O) thermophysical properties.

Backend: CoolProp, using the Lemmon & Span (2006) equation of state, the same model NIST uses.
All inputs and outputs are SI: K, Pa, kg/m^3, J/kg, J/(kg K), N/m.
"""

from dataclasses import dataclass

import CoolProp.CoolProp as CP

FLUID = "NitrousOxide"

# Fixed points of N2O, read from the equation of state itself
T_CRIT = CP.PropsSI("Tcrit", FLUID) # K (~309.52 K = 36.37 C)
P_CRIT = CP.PropsSI("pcrit", FLUID) # Pa (~72.45 bar)
RHO_CRIT = CP.PropsSI("rhocrit", FLUID) # kg/m^3 (~452)
T_TRIPLE = CP.PropsSI("Ttriple", FLUID) # K (~182.3 K = -90.8 C, below this it freezes)
MOLAR_MASS = CP.PropsSI("molarmass", FLUID) # kg/mol (0.04401)


@dataclass(frozen=True)
class SaturationState:
    """Everything about saturated N2O at one temperature.

    _l = saturated liquid, _v = saturated vapour.
    """

    T: float        # temperature, K
    p: float        # vapour pressure, Pa
    rho_l: float    # density, kg/m^3
    rho_v: float
    h_l: float      # specific enthalpy, J/kg
    h_v: float
    u_l: float      # specific internal energy, J/kg
    u_v: float
    s_l: float      # specific entropy, J/(kg K)
    s_v: float
    cp_l: float     # specific heat at constant pressure, J/(kg K)
    cp_v: float
    sigma: float    # surface tension, N/m

    @property
    def h_fg(self) -> float:
        """Latent heat of vaporisation, J/kg: energy to boil 1 kg of liquid at this T."""
        return self.h_v - self.h_l


def saturation(T: float) -> SaturationState:
    """Saturated liquid and vapour properties of N2O at temperature T (kelvin)."""
    if not (T_TRIPLE <= T < T_CRIT):
        raise ValueError(
            f"T = {T:.2f} K is outside the liquid+vapour range "
            f"({T_TRIPLE:.2f} K to {T_CRIT:.2f} K). Above T_crit there is no liquid."
        )

    def liquid(prop):
        return CP.PropsSI(prop, "T", T, "Q", 0, FLUID)

    def vapour(prop):
        return CP.PropsSI(prop, "T", T, "Q", 1, FLUID)

    return SaturationState(
        T=T,
        p=liquid("P"),
        rho_l=liquid("D"), rho_v=vapour("D"),
        h_l=liquid("H"),   h_v=vapour("H"),
        u_l=liquid("U"),   u_v=vapour("U"),
        s_l=liquid("S"),   s_v=vapour("S"),
        cp_l=liquid("C"),  cp_v=vapour("C"),
        sigma=liquid("I"),
    )


def saturation_from_pressure(p: float) -> SaturationState:
    """Same as saturation(), but you give the pressure (Pa) instead of the temperature."""
    T = CP.PropsSI("T", "P", p, "Q", 0, FLUID)
    return saturation(T)


def state(output: str, input1: str, value1: float, input2: str, value2: float) -> float:
    """Any single property, any phase, e.g. state("P", "T", 300, "D", 600).

    Use this for non-saturated N2O: compressed liquid, superheated vapour, supercritical.
    """
    return CP.PropsSI(output, input1, value1, input2, value2, FLUID)


if __name__ == "__main__":
    C = 273.15
    print(f"Critical point: {T_CRIT - C:.2f} C, {P_CRIT / 1e5:.2f} bar\n")
    print("  T [C]   p [bar]  rho_l   rho_v   h_fg [kJ/kg]")
    for T_c in [-20, 0, 10, 20, 25, 30, 35]:
        s = saturation(T_c + C)
        print(f"{T_c:7.0f} {s.p / 1e5:9.1f} {s.rho_l:7.1f} {s.rho_v:7.1f} {s.h_fg / 1e3:10.1f}")