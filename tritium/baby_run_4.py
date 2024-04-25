from simple_tritium_transport_model import ureg, Model
import numpy as np
from helpers import substract_background_from_measurements, cumulative_activity

background_2 = 0.28 * ureg.Bq
raw_measurements = {
    2: {
        1: 0.299 * ureg.Bq,
        2: 0.277 * ureg.Bq,
        3: 0.283 * ureg.Bq,
        4: 0.274 * ureg.Bq,
        "background": background_2,
    },
    3: {
        1: 0.279 * ureg.Bq,
        2: 0.276 * ureg.Bq,
        3: 0.564 * ureg.Bq,
        4: 0.290 * ureg.Bq,
        "background": background_2,
    },
    4: {
        1: 0.331 * ureg.Bq,
        2: 0.296 * ureg.Bq,
        3: 2.100 * ureg.Bq,
        4: 0.310 * ureg.Bq,
        "background": background_2,
    },
    5: {
        1: 0.269 * ureg.Bq,
        2: 0.284 * ureg.Bq,
        3: 4.050 * ureg.Bq,
        4: 0.369 * ureg.Bq,
        "background": background_2,
    },
    6: {
        1: 0.247 * ureg.Bq,
        2: 0.308 * ureg.Bq,
        3: 8.469 * ureg.Bq,
        4: 0.754 * ureg.Bq,
        "background": background_2,
    },
    7: {
        1: 0.280 * ureg.Bq,
        2: 0.292 * ureg.Bq,
        3: 3.439 * ureg.Bq,
        4: 0.553 * ureg.Bq,
        "background": background_2,
    },
    8: {
        1: 0.260 * ureg.Bq,
        2: 0.275 * ureg.Bq,
        3: 1.243 * ureg.Bq,
        4: 0.480 * ureg.Bq,
        "background": background_2,
    },
}

measurements_after_background_sub = substract_background_from_measurements(
    raw_measurements
)

replacement_times = [
    12 * ureg.hour,
    24 * ureg.hour,
    36 * ureg.hour,
    48 * ureg.hour,
    3 * ureg.day,
    4 * ureg.day,
    6 * ureg.day,
]

replacement_times = sorted(replacement_times)

# # Cumulative values
cumulative_release = cumulative_activity(measurements_after_background_sub)

# Model

baby_diameter = 1.77 * ureg.inches - 2 * 0.06 * ureg.inches  # from CAD drawings
baby_radius = 0.5 * baby_diameter
baby_volume = 0.100 * ureg.L
baby_cross_section = np.pi * baby_radius**2
baby_height = baby_volume / baby_cross_section
baby_model = Model(
    radius=baby_radius,
    height=baby_height,
    TBR=3.3e-4 * ureg.particle * ureg.neutron**-1,  # stefano 10/24/2023
)

fitting_param = 1.1

mass_transport_coeff_factor = 3

baby_model.k_top *= mass_transport_coeff_factor
baby_model.k_wall *= mass_transport_coeff_factor

exposure_time = 12 * ureg.hour

baby_model.irradiations = [
    [0 * ureg.hour, 0 + exposure_time],
    [24 * ureg.hour, 24 * ureg.hour + exposure_time],
]

baby_model.neutron_rate = fitting_param * (1.2e8 + 3.96e8) * ureg.neutron * ureg.s**-1
