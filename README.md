# UTSRocketrySim

**An open-source simulator for nitrous oxide hybrid rocket engines.**
Built by the UTS Rocketry Team for Australian university rocketry.

![Python](https://img.shields.io/badge/python-3.10%2B-blue)
![Status](https://img.shields.io/badge/status-early%20development-orange)
![Licence](https://img.shields.io/badge/licence-MIT-green)

---

## Why this exists

Nitrous oxide hybrids are the go-to engine for university teams. N2O is self-pressurising and
far safer to handle than cryogenic propellants. But it's also the hardest part to model:

- In the tank it's a **liquid and a vapour at the same time**, and the balance shifts with
  temperature. On a 35 °C desert pad, a tank can go **completely full of liquid** or
  **supercritical** (no liquid at all) before you launch.
- It **flashes to vapour inside the injector**, so simple liquid flow equations can be wrong
  by a factor of two or more.
- The fuel grain **burns back unevenly** along its length, so thrust and O/F ratio change
  throughout the burn.

Spreadsheets can't capture this well. UTSRocketrySim models it from first principles, in code
that anyone can read, check and improve.

## Features

| | Feature | Status |
|---|---|---|
| **N2O properties**: saturation, any-phase state, critical and triple-point limits | 
| Independent cross-check (ESDU 91022) and automated tests | 
| **Tank model**: pressure, ullage, liquid-full and supercritical warnings, boil-off | 
| **Tank blowdown** during the burn | 
| **Injector flow**: SPI, HEM and Dyer two-phase models, compared | 
| **Pad thermal model**: desert heating, venting | 
| **Fuel grains**: HTPB, paraffin, ABS and blends; regression rate; fit (a, n) from your test fires | 
| **3D burn-back**: import your grain from CAD (STL), watch the surface regress in 3D | 
| **Coupled engine**: chamber pressure, thrust curve, Isp, O/F shift, from NASA CEA tables | 

> **Not reinventing CEA.** Combustion chemistry comes from
> [NASA CEA](https://cearun.grc.nasa.gov/): you run it, we import the tables.

## Installation

Requires Python 3.10 or newer.

```bash
git clone https://github.com/HayWire737/UTSRocketrySim.git
cd UTSRocketrySim
python3 -m venv .venv
source .venv/bin/activate        # Windows: .venv\Scripts\activate
pip install -e .
```

Check it works:

```bash
python -m utsrocketrysim.n2o.properties
```

## Quick start

```python
from utsrocketrysim.n2o.properties import saturation, state, T_CRIT

C = 273.15  # all inputs are SI, so temperatures are in kelvin

# Saturated N2O at 20 °C
s = saturation(20 + C)
print(f"Vapour pressure: {s.p / 1e5:.1f} bar")           # 50.5 bar
print(f"Liquid density:  {s.rho_l:.1f} kg/m^3")          # 785.1
print(f"Vapour density:  {s.rho_v:.1f} kg/m^3")          # 158.0
print(f"Latent heat:     {s.h_fg / 1e3:.1f} kJ/kg")      # 169.9

# A tank that has gone liquid-full: 30 °C at an average density of 750 kg/m^3
p = state("P", "T", 30 + C, "D", 750)
print(f"Liquid-full pressure: {p / 1e5:.1f} bar")       # 79.7 bar, vs 63.1 bar saturated

print(f"Critical temperature: {T_CRIT - C:.2f} °C")      # 36.37 °C, no liquid above this
```

## Project structure

```
UTSRocketrySim/
├── src/utsrocketrysim/      the installable package
│   ├── n2o/                 nitrous oxide: properties, tank, injector
│   ├── grain/               fuel grains, regression, 3D burn-back
│   ├── combustion/          CEA import, chamber, nozzle
│   └── engine.py            ties it all together
├── tests/                   automated checks
├── examples/                scripts showing how to use it
└── docs/                    physics derivations and notes
```

## How we check the physics

Every model in this project must be tested against something **independent of itself**:

- a second, independent source (e.g. CoolProp's equation of state vs the ESDU 91022 correlations),
- published experimental data, or
- a conservation law (mass, energy) evaluated numerically.

"It runs" doesn't count as a test. If you can't check it, flag it.

## Contributing

Contributions from any university team are welcome.

1. Create a branch: `git checkout -b feature/your-feature`
2. Make your changes. Every physics function needs a docstring naming its source (paper and
   equation), plus a test.
3. Keep all units **SI** inside the code. Convert only for input and display.
4. Open a pull request into `main` for review.

## References

- Lemmon, E.W. & Span, R. (2006). Short fundamental equations of state for 20 industrial
  fluids. *J. Chem. Eng. Data* 51(3). The N2O equation of state used via CoolProp.
- ESDU 91022 (1991). *Thermophysical properties of nitrous oxide.*
- Dyer, J. et al. (2007). Modeling feed system flow physics for self-pressurizing propellants.
  AIAA 2007-5702.
- Zimmerman, J.E. et al. (2013). Review and evaluation of models for self-pressurizing
  propellant tank dynamics. AIAA 2013-4045.

## Safety

This is a design and learning tool, **not a guarantee of safety**. Validate against your own
test data, use appropriate safety factors, and follow your range safety officer and
AURC/Tripoli rules. Nitrous oxide can decompose violently. Keep it clean, compatible and
away from heat.

## Licence

MIT. Free to use, modify and share, including by other universities. See [LICENSE](LICENSE).