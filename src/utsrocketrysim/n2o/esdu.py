"""ESDU 91022 saturation correlations for nitrous oxide.

An independent cross-check of the CoolProp equation of state in properties.py. These are simple
curve fits of the saturation line, published by ESDU in 1991 and popularised for rocketry by
Aspire Space ("Modelling the nitrous run tank emptying").

Valid from about -90 C up to the critical point. Saturation only: no compressed liquid,
superheated vapour or supercritical states.
"""

import math
from dataclasses import dataclass

# ESDU's own critical constants. Slightly different from CoolProp's, and that's fine:
# they're part of the fit.
TC = 309.57    # K
PC = 7.251e6   # Pa
RHOC = 452.0   # kg/m^3


@dataclass(frozen=True)
class ESDUSaturation:
    """Saturated N2O properties from the ESDU correlations. Units are SI.

    Enthalpies use ESDU's own reference point, so compare h_fg with CoolProp, not h_l or h_v.
    """

    T: float
    p: float
    rho_l: float
    rho_v: float
    h_l: float
    h_v: float
    cp_l: float

    @property
    def h_fg(self) -> float:
        return self.h_v - self.h_l


def esdu_saturation(T: float) -> ESDUSaturation:
    """Saturated N2O properties at temperature T (kelvin), from ESDU 91022."""
    if not (183.0 <= T < TC):
        raise ValueError(f"T = {T:.2f} K is outside the ESDU range (183 K to {TC} K)")

    Tr = T / TC        # reduced temperature: 1.0 at the critical point
    a = 1.0 - Tr       # how far below the critical point we are
    b = 1.0 / Tr - 1.0

    p = PC * math.exp(
        (1.0 / Tr) * (-6.71893 * a + 1.35966 * a**1.5 - 1.3779 * a**2.5 - 4.051 * a**5)
    )
    rho_l = RHOC * math.exp(
        1.72328 * a**(1 / 3) - 0.8395 * a**(2 / 3) + 0.5106 * a - 0.10412 * a**(4 / 3)
    )
    rho_v = RHOC * math.exp(
        -1.009 * b**(1 / 3) - 6.28792 * b**(2 / 3) + 7.50332 * b
        - 7.90463 * b**(4 / 3) + 0.629427 * b**(5 / 3)
    )
    # ESDU gives these in kJ/kg and kJ/(kg K); multiply by 1000 for SI
    h_l = 1e3 * (-200.0 + 116.043 * a**(1 / 3) - 917.225 * a**(2 / 3)
                 + 794.779 * a - 589.587 * a**(4 / 3))
    h_v = 1e3 * (-200.0 + 440.055 * a**(1 / 3) - 459.701 * a**(2 / 3)
                 + 434.081 * a - 485.338 * a**(4 / 3))
    cp_l = 1e3 * 2.49973 * (1.0 + 0.023454 / a - 3.80136 * a + 13.0945 * a**2 - 14.518 * a**3)

    return ESDUSaturation(T=T, p=p, rho_l=rho_l, rho_v=rho_v, h_l=h_l, h_v=h_v, cp_l=cp_l)


if __name__ == "__main__":
    from utsrocketrysim.n2o.properties import saturation

    C = 273.15
    print("Percentage difference, ESDU vs CoolProp")
    print("  T [C]     p     rho_l   rho_v   h_fg    cp_l")
    for T_c in [-40, -20, 0, 10, 20, 25, 30, 35]:
        e, c = esdu_saturation(T_c + C), saturation(T_c + C)
        diff = [100 * (getattr(e, k) - getattr(c, k)) / getattr(c, k)
                for k in ("p", "rho_l", "rho_v", "h_fg", "cp_l")]
        print(f"{T_c:7.0f}" + "".join(f"{d:8.2f}" for d in diff))