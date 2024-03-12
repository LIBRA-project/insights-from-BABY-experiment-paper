from simple_tritium_transport_model import ureg, Model
import numpy as np
from helpers import substract_background_from_measurements, cumulative_activity

background_1 = 0.265 * ureg.Bq

raw_measurements = {
    1: {
        1: 0.279 * ureg.Bq,
        2: 0.284 * ureg.Bq,
        3: 6.352 * ureg.Bq,
        4: 0.617 * ureg.Bq,
        "background": background_1,
    },
    2: {
        1: 0.272 * ureg.Bq,
        2: 0.278 * ureg.Bq,
        3: 6.621 * ureg.Bq,
        4: 0.581 * ureg.Bq,
        "background": background_1,
    },
    3: {
        1: 0.252 * ureg.Bq,
        2: 0.272 * ureg.Bq,
        3: 2.347 * ureg.Bq,
        4: 0.706 * ureg.Bq,
        "background": background_1,
    },
}

measurements_after_background_sub = substract_background_from_measurements(
    raw_measurements
)

replacement_times = [
    1 * ureg.day,
    2 * ureg.day,
    4 * ureg.day,
]

cumulative_release = cumulative_activity(measurements_after_background_sub)

baby_diameter = 1.77 * ureg.inches - 2 * 0.06 * ureg.inches  # from CAD drawings
baby_radius = 0.5 * baby_diameter
baby_volume = 0.125 * ureg.L
baby_cross_section = np.pi * baby_radius**2
baby_height = baby_volume / baby_cross_section
baby_model = Model(
    radius=baby_radius,
    height=baby_height,
    TBR=3.3e-4 * ureg.particle * ureg.neutron**-1,  # stefano 10/24/2023
)

fitting_param = 0.84

mass_transport_coeff_factor = 3

baby_model.k_top *= mass_transport_coeff_factor
baby_model.k_wall *= mass_transport_coeff_factor

exposure_time = 12 * ureg.hour
baby_model.irradiations = [
    [0 * ureg.hour, 0 + exposure_time],
    [24 * ureg.hour, 24 * ureg.hour + exposure_time],
]

baby_model.neutron_rate = fitting_param * (1.2e8 + 3.96e8) * ureg.neutron * ureg.s**-1
