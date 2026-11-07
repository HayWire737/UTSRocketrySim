# UTSRocketrySim

Open-source simulator for **nitrous oxide hybrid rocket engines**, built by the UTS Rocketry Team
for Australian university rocketry.

- **Two-phase N2O**: real-fluid properties, tank pressure/ullage vs temperature, liquid-full and
  supercritical limits, venting/boil-off, blowdown, two-phase injector flow
- **Solid fuel**: regression-rate laws, fuel blends (HTPB, paraffin, ABS, additives), thrust and
  O/F curves, and **3D grain burn-back** from your CAD (STL)
- **Your own CEA**: import NASA CEA tables. This project does not reimplement chemical equilibrium.

> ⚠️ Simulation is not a guarantee of safety. Validate against your own test data and follow
> range safety practice.

## Status
| Phase | Scope | State |
|---|---|---|
| 1 | N2O properties (2 independent sources), static tank, venting | ✅ |
| 2 | Blowdown, injector models (SPI/HEM/Dyer), feed line, pad thermal | ⏳ |
| 3 | Fuel blends, regression, 2D/3D burn-back, CAD import, 3D viewer | ⏳ |
| 4 | Coupled engine, CEA import, thrust curves, validation | ⏳ |

See [`docs/ARCHITECTURE.md`](docs/ARCHITECTURE.md) for the roadmap and
[`docs/physics/`](docs/physics) for derivations.

## Install
```bash
git clone https://github.com/<owner>/UTSRocketrySim.git
cd UTSRocketrySim
python3 -m venv .venv && source .venv/bin/activate
pip install -e ".[dev]"          # add ,grain3d for 3D burn-back (phase 3)
pytest                           # all tests should pass
```

## Quick start
```python
from utsrocketrysim.n2o import N2O, tank_state, mass_for_fill, liquid_full_temperature

C = 273.15
sat = N2O.saturation(20 + C)
print(sat.p / 1e5, sat.rho_l, sat.rho_v)        # 50.5 bar, 785 kg/m3, 158 kg/m3

V = 0.008                                         # 8 L run tank
m = mass_for_fill(V, T_fill=15 + C, ullage=0.15)  # 5.74 kg
print(liquid_full_temperature(V, m) - C)          # 27.5 C -> no ullage left above this
```

```bash
python examples/01_desert_tank_heating.py         # writes output/01_desert_tank_heating.png
```

## Contributing
Branch → pull request → review. Every new physics model needs a docstring naming its source and a
test against something independent of it (another correlation, published data, or a
conservation law). Keep units SI internally.

## Licence
MIT. See [LICENSE](LICENSE).
