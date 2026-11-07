# References

## Nitrous oxide properties and tank / feed modelling
- Lemmon, E.W. & Span, R. (2006). Short fundamental equations of state for 20 industrial fluids.
  *J. Chem. Eng. Data* 51(3). — the N2O EOS used by CoolProp / NIST.
- ESDU 91022 (1991). Thermophysical properties of nitrous oxide. — saturation correlations.
- Dyer, J., Doran, E., Dunn, Z., Lohner, K., Bayart, C., Sadhwani, A., Zilliac, G., Cantwell, B.,
  Karabeyoglu, A. (2007). Modeling feed system flow physics for self-pressurizing propellants.
  AIAA 2007-5702. — NHNE injector model.
- Whitmore, S.A. & Chandler, S.N. (2010). Engineering model for self-pressurizing
  saturated-N2O-propellant feed systems. *J. Propulsion & Power* 26(4).
- Zimmerman, J.E., Waxman, B.S., Cantwell, B.J., Zilliac, G. (2013). Review and evaluation of
  models for self-pressurizing propellant tank dynamics. AIAA 2013-4045.

## Hybrid combustion and regression
- Marxman, G.A., Wooldridge, C.E., Muzzy, R.J. (1964). Fundamental mechanisms of combustion in
  hybrid rocket engines. — boundary-layer diffusion-limited regression, r_dot ∝ G^0.8.
- Karabeyoglu, M.A., Altman, D., Cantwell, B.J. (2002). Combustion of liquefying hybrid
  propellants, Part 1. *J. Propulsion & Power* 18(3). — paraffin entrainment.

## Project library (in the claude.ai project "UTS Rocketry Sim")
- `InstabilityPaperJPPNov2005.pdf`, `Transient_Modeling_of_Hybrid_Rocket_Low_Frequency_.pdf` —
  hybrid low-frequency instability / transient chamber models (phase 4).
- `Surmacz_Hybridrocketpropulsiondevelopmentandapplication.pdf`
- `ArcachonHybridv10b2012peractaastronautica.pdf`
- `A_systematic_review_on_computational_modelling_and.pdf`
- `Numerical_modeling_methods_of_in-chamber_processes.pdf`
- `1-s2.0-S0016236125027619-main.pdf`, `19680010362.pdf`, `AL11-13.pdf`
- `HalfCatSim_v1.3.11.xlsx` (Half Cat Rocketry, GPL-3) — liquid biprop blowdown spreadsheet.
  Uses NIST saturation tables, SPI injector flow, and an empirical tank-pressure decay. We use it
  only as a comparison point and do not copy its code.
