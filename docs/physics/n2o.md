# N2O physics: derivations

## Venting to hold a tank at constant temperature

Fixed tank volume. Vapour is vented at enthalpy h_v. Both phases stay saturated at T.

Volume stays constant:  dm_l / rho_l + dm_v / rho_v = 0   →   dm_v = -(rho_v / rho_l) dm_l

Mass vented:  dm_vent = -(dm_l + dm_v) = -dm_l (1 - rho_v / rho_l)

Energy:  d(m_l u_l + m_v u_v) = dQ - h_v dm_vent

Substitute both, and use h = u + p/rho:

    dm_l [ u_l - u_v rho_v/rho_l - h_v (1 - rho_v/rho_l) ] = dQ

The bracket simplifies to h_l - h_v = -h_fg, so dm_l = -dQ / h_fg, giving

    m_dot_vent = Q_dot (1 - rho_v / rho_l) / h_fg

Checked numerically in tests/test_tank.py::test_vent_rate_holds_temperature.
At 25 °C: about 0.5 g/s per 100 W of heat leaking in (1.8 kg/h).