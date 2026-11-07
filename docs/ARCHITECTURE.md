# UTSRocketrySim — Architecture & Roadmap

A hybrid rocket engine simulator built for checking from first principles. Every model states its
assumptions, cites its source, and has at least one test against something *independent* of
itself: a second correlation, published test data, or a conservation law evaluated numerically.

**Not in scope:** chemical equilibrium. We do **not** rebuild NASA CEA. You run CEA and we import
its tables of c*, γ, molecular weight and T_c against (P_c, O/F).

## Package layout

```
src/utsrocketrysim/
  n2o/            Nitrous oxide
    properties.py   EOS (CoolProp) + ESDU 91022 cross-check backend        [done]
    tank.py         Equilibrium tank state, liquid-full, venting           [done]
    blowdown.py     Transient tank emptying: equilibrium + non-equilibrium [phase 2]
    injector.py     SPI / HEM / Dyer NHNE two-phase injector models        [phase 2]
    feedline.py     Line losses, flashing, priming, adiabatic-compression  [phase 2]
    thermal.py      Pad heating: solar + convection into tank, vent/relief [phase 2]
  combustion/
    cea.py          Import NASA CEA output -> interpolation tables         [phase 4]
    chamber.py      Chamber filling ODE, c* efficiency, nozzle, thrust     [phase 4]
  grain/
    materials.py    Fuel blends: density, regression law (a, n), sources   [phase 3]
    regression.py   r_dot = a * G_ox^n, axial marching of mass flux        [phase 3]
    levelset.py     3D voxel signed-distance burn-back (fast marching)     [phase 3]
    cad.py          STL/OBJ import (Onshape export) -> voxel SDF           [phase 3]
    viz3d.py        PyVista surface mesh per timestep, slider, VTK export  [phase 3]
  engine.py       Couples everything: tank -> injector -> chamber <- grain [phase 4]
  io/             YAML engine definition files, CSV/plot output            [phase 4]
```

## Phases

### Phase 1 — N2O properties and static tank (this commit)
- Saturation properties from two independent sources, agreeing within 0.4% below 20 °C.
- Fixed-mass tank vs temperature: pressure, ullage, liquid-full point, supercritical flag.
- Constant-temperature vent rate under heat leak (closed-form, verified numerically).

### Phase 2 — N2O transient: blowdown, injector, feed system
- **Tank blowdown** (m, U, V) → state using a CoolProp (ρ, u) flash. Start with the equilibrium
  model (Whitmore & Chandler; Zimmerman 2013), then add a two-node non-equilibrium model
  (separate liquid and vapour control volumes with interfacial heat/mass transfer), because
  equilibrium models are known to under-predict pressure drop in fast blowdown.
- **Injector**: implement all three and compare them, because the choice changes mass flow by
  2–3×:
  - SPI (single-phase incompressible): `ṁ = Cd·A·√(2ρΔP)`. Over-predicts flow for saturated
    N2O. HalfCatSim uses this.
  - HEM (homogeneous equilibrium, isentropic flash). Under-predicts for short orifices.
  - **Dyer NHNE** (Dyer et al. 2007): blends the two by `κ = √((P1−Pc)/(Pv−P2))`. This one is
    the recommended default.
- **Feed line** priming, vapour-slug at valve opening, and adiabatic-compression temperature rise
  (a known N2O decomposition ignition hazard).
- **Pad thermal model**: solar + ambient convection into tank, so we can predict when the tank
  goes liquid-full or supercritical through a desert day, and how much it vents to stay put.

### Phase 3 — Fuel grain
- **Regression law** `r_dot = a·G_ox^n` (optionally `a·G^n·x^m`) with literature (a, n) per fuel.
  Blend density and enthalpy by mass-fraction mixing. Regression constants for a *new* blend
  cannot be predicted from composition, so we also ship a **ballistic reconstruction tool** that
  fits (a, n) from your own static fires (initial/final port geometry, burn time, ṁ_ox trace).
- **2D**: axisymmetric and arbitrary cross-sections (like openMotor's fast-marching maps).
- **3D**: grain = voxel signed-distance field φ (φ<0 = port). Each step:
  1. slice along the axis, measure port area A(z) and wetted perimeter P(z);
  2. march G(z) = (ṁ_ox + ṁ_fuel added upstream)/A(z) down the port;
  3. speed field F(z) = r_dot(G(z)); update φ ← φ − F·Δt; re-initialise φ with fast marching;
  4. extract the burning surface by marching cubes → area, port volume, fuel mass flow.
  This captures axial regression variation and lets finocyl fins burn out realistically.
- **CAD import**: export STL from Onshape → `trimesh` voxelise → φ. No hand-specified geometry.
- **3D viewer**: PyVista surface mesh with a time slider, coloured by local regression rate;
  export `.vtk` series for ParaView.

### Phase 4 — Coupled engine
- Chamber pressure ODE (filling/emptying): `dPc/dt = (RT/V_c)(ṁ_ox + ṁ_f − ṁ_nozzle)` with
  c*, γ, T from imported CEA tables at the current (Pc, O/F). c* efficiency as a user input.
- Thrust, Isp, O/F shift over the burn, total impulse, motor class.
- Engine definition files (YAML) so a design is one file in the repo.
- Validation against published static-fire data (see the papers in `docs/references.md`).

## Conventions
- SI internally, always. Convert only at input/output.
- Every physics function's docstring names its source (paper, equation number).
- Every model gets a test against an independent check, not just "it runs".
- Branch → PR → review. `main` always passes `pytest`.

## Open questions (to resolve together)
1. **"Decomposition for startup."** In a hybrid, the igniter heats the fuel surface and N2O
   reacts with the pyrolysis gases. Ignition doesn't depend on N2O decomposing in the tank.
   N2O does decompose exothermically in the hot chamber (2N2O → 2N2 + O2, ~82 kJ/mol), which
   CEA already accounts for. In the tank and lines, decomposition is a **hazard** to rule out,
   triggered by adiabatic compression, contamination or flashback. Which do we model: chamber
   ignition transient, feed-line priming, or a decomposition-hazard check (or all three)?
2. **FireHawk feed system.** Blowdown or pump-fed? Pumping saturated N2O has near-zero NPSH, so
   the tank model and the cavitation check differ a lot between the two.
3. **"Aluminium oxide" in fuel.** Al2O3 is inert (already oxidised): it adds mass and acts as an
   opacifier, but releases no energy. Aluminium *powder* is the energetic additive. Which do we
   mean, or both?
4. Licence: MIT (permissive) is the default here. HalfCatSim is GPL-3, so we re-derive from the
   papers rather than copying its formulas.
