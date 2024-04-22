from simple_tritium_transport_model import ureg, Model
import numpy as np
from helpers import (
    substract_background_from_measurements,
    cumulative_activity,
    background_sub,
)


background = 0.314 * ureg.Bq
background_2_2_2024 = 0.310 * ureg.Bq
raw_measurements = {
    1: {
        1: 0.288 * ureg.Bq,
        2: 0.310 * ureg.Bq,
        3: 4.455 * ureg.Bq,
        4: 0.396 * ureg.Bq,
        "background": background,
    },
    2: {
        1: 0.306 * ureg.Bq,
        2: 0.300 * ureg.Bq,
        3: 3.127 * ureg.Bq,
        4: 0.409 * ureg.Bq,
        "background": background,
    },
    3: {
        1: 0.358 * ureg.Bq,
        2: 0.286 * ureg.Bq,
        3: 5.735 * ureg.Bq,
        4: 0.442 * ureg.Bq,
        "background": background,
    },
    4: {
        1: 0.272 * ureg.Bq,
        2: 0.299 * ureg.Bq,
        3: 3.915 * ureg.Bq,
        4: 0.398 * ureg.Bq,
        "background": background,
    },
    5: {
        1: 0.301 * ureg.Bq,
        2: 0.305 * ureg.Bq,
        3: 2.930 * ureg.Bq,
        4: 0 * ureg.Bq,  # missing water!!
        "background": background,
    },
}

measurements_after_background_sub = substract_background_from_measurements(
    raw_measurements
)

# TODO find a way to replace
measurements_after_background_sub[6] = {
    1: background_sub(0.277 * ureg.Bq, background_2_2_2024),
    2: background_sub(0.310 * ureg.Bq, background_2_2_2024),
    3: background_sub(1.191 * ureg.Bq, background),
    4: background_sub(0.361 * ureg.Bq, background),
}
measurements_after_background_sub[7] = {
    1: background_sub(0.276 * ureg.Bq, background_2_2_2024),
    2: background_sub(0.283 * ureg.Bq, background_2_2_2024),
    3: background_sub(0.991 * ureg.Bq, background_2_2_2024),
    4: background_sub(0.635 * ureg.Bq, background_2_2_2024),
}

# time starts at 01/25 9:36 AM
# 01/26 9:36 AM = 24 hours = 1 * ureg.day + 0 * ureg.hour + 0 * ureg.minute
# 01/27 9:36 AM = 48 hours = 2 * ureg.day + 0 * ureg.hour + 0 * ureg.minute
replacement_times = [
    # 01/25 21:44
    0 * ureg.day + 15 * ureg.hour + 8 * ureg.minute,
    # 01/26 09:22
    1 * ureg.day
    + (9 * ureg.hour + 22 * ureg.minute)
    - (9 * ureg.hour + 36 * ureg.minute),
    # 01/26 21:46
    1 * ureg.day
    + (21 * ureg.hour + 46 * ureg.minute)
    - (9 * ureg.hour + 36 * ureg.minute),
    # 01/27 10:48
    2 * ureg.day
    + (10 * ureg.hour + 48 * ureg.minute)
    - (9 * ureg.hour + 36 * ureg.minute),
    # 01/28 11:19
    3 * ureg.day
    + (11 * ureg.hour + 19 * ureg.minute)
    - (9 * ureg.hour + 36 * ureg.minute),
    # 01/29 08:02
    4 * ureg.day
    + (8 * ureg.hour + 2 * ureg.minute)
    - (9 * ureg.hour + 36 * ureg.minute),
    # 01/31 17:00
    6 * ureg.day
    + (17 * ureg.hour + 0 * ureg.minute)
    - (9 * ureg.hour + 36 * ureg.minute),
]

replacement_times = sorted(replacement_times)

# # Cumulative values

cumulative_release = cumulative_activity(measurements_after_background_sub)

# Model

baby_diameter = 1.77 * ureg.inches - 2 * 0.06 * ureg.inches  # from CAD drawings
baby_radius = 0.5 * baby_diameter
baby_volume = 100 * ureg.mL
baby_cross_section = np.pi * baby_radius**2
baby_height = baby_volume / baby_cross_section
calculated_TBR = 4.71e-4 * ureg.particle * ureg.neutron**-1
baby_model = Model(
    radius=baby_radius,
    height=baby_height,
    TBR=calculated_TBR,
)


mass_transport_coeff_factor = 3

baby_model.k_top *= mass_transport_coeff_factor * 0.7
optimised_ratio = 3e-2
baby_model.k_wall = baby_model.k_top * optimised_ratio

exposure_time = 12 * ureg.hour

baby_model.irradiations = [
    [0 * ureg.hour, 0 + exposure_time],
    [24 * ureg.hour, 24 * ureg.hour + exposure_time],
]

# calculated from Kevin's activation foil analysis
P383_neutron_rate = 4.95e8 * ureg.neutron * ureg.s**-1
A325_neutron_rate = 2.13e8 * ureg.neutron * ureg.s**-1

neutron_rate_relative_uncertainty = 0.089
baby_model.neutron_rate = (
    P383_neutron_rate + A325_neutron_rate
) / 2  # the neutron rate is divided by two to acount for the double counting (two detectors)
