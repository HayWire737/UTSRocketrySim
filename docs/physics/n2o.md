# N2O physics — derivations to check

Each result here is implemented in `src/utsrocketrysim/n2o/` and tested in `tests/test_n2o.py`.

## 1. Property sources
| Backend | Basis | Range | Role |
|---|---|---|---|
| `CoolPropN2O` | Helmholtz EOS, Lemmon & Span (2006) | all phases | simulation |
| `ESDUN2O` | ESDU 91022 saturation correlations | −90 °C to T_c | independent check |

The NIST Webbook (and HalfCatSim's table) use the same EOS as CoolProp, so they match it exactly.
That checks units and wrapping, not physics. ESDU is the independent check.

Agreement (ESDU vs EOS): p_sat < 0.3%, ρ_l < 0.6%, ρ_v < 0.7%, h_fg < 0.4% from −40 to 20 °C.
Near T_c, ρ_v and h_fg drift 2–3%, and c_p,l diverges (−32% at 35 °C) because c_p → ∞ at the
critical point. **Trust the EOS above about 25 °C.**

Critical point (EOS): T_c = 309.52 K (36.37 °C), p_c = 72.45 bar, ρ_c ≈ 452 kg/m³.

Reference check: saturated liquid at 20 °C is **785 kg/m³** (vapour 158 kg/m³, p = 50.5 bar).
743 kg/m³ is the value at 25 °C.

## 2. Equilibrium tank at fixed mass
For tank volume V and mass m, the specific volume v = V/m is fixed. At temperature T:

  x = (v − v_l)/(v_v − v_l),    ullage = x·v_v / v

- 0 < x < 1: two-phase, p = p_sat(T).
- x ≤ 0: **liquid-full**. Found by solving ρ_l,sat(T) = m/V. Past that point, p comes from the
  compressed-liquid EOS, and dp/dT jumps from ~1.2 bar/K (saturation curve) to ~4.6 bar/K
  ((∂p/∂T)_ρ = β/κ_T) near 20 °C.
- T ≥ T_c: supercritical. No liquid, so there's nothing for the injector to feed as liquid.

Example (`examples/01_desert_tank_heating.py`): an 8 L tank filled at 15 °C with 20% ullage goes
liquid-full at 30.4 °C. With 10% ullage, it goes liquid-full at 24.0 °C.

## 3. Venting to hold temperature
Two-phase, fixed V, vapour vented at h_v, both phases stay saturated at T:

  volume:  dm_l/ρ_l + dm_v/ρ_v = 0          ⇒ dm_v = −(ρ_v/ρ_l) dm_l
  mass:    dm_vent = −(dm_l + dm_v) = −dm_l (1 − ρ_v/ρ_l)
  energy:  d(m_l u_l + m_v u_v) = δQ − h_v dm_vent

Substitute, and use h = u + p/ρ:
  dm_l [u_l − u_v ρ_v/ρ_l − h_v (1 − ρ_v/ρ_l)] = δQ
  bracket = (u_l + p/ρ_l) − h_v = h_l − h_v = −h_fg

  ⇒ **ṁ_vent = Q̇ (1 − ρ_v/ρ_l) / h_fg**

At 25 °C this is 0.51 g/s per 100 W of heat leak (1.8 kg/h).
Numerical check: the test does one step with a direct (ρ, u) EOS flash and confirms
|ΔT_vented| < 10⁻³ |ΔT_sealed|.
