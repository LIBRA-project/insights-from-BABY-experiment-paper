from simple_tritium_transport_model import ureg, Model
import numpy as np
from helpers import substract_background_from_measurements, cumulative_activity


background = 0.29 * ureg.Bq
raw_measurements = {
    1: {
        1: 0.254 * ureg.Bq,
        2: 0.289 * ureg.Bq,
        3: 3.885 * ureg.Bq,
        4: 0.344 * ureg.Bq,
        "background": background,
    },
    2: {
        1: 0.262 * ureg.Bq,
        2: 0.283 * ureg.Bq,
        3: 3.998 * ureg.Bq,
        4: 0.408 * ureg.Bq,
        "background": background,
    },
    3: {
        1: 0.266 * ureg.Bq,
        2: 0.293 * ureg.Bq,
        3: 4.644 * ureg.Bq,
        4: 0.396 * ureg.Bq,
        "background": background,
    },
    4: {
        1: 0.275 * ureg.Bq,
        2: 0.278 * ureg.Bq,
        3: 5.098 * ureg.Bq,
        4: 0.505 * ureg.Bq,
        "background": background,
    },
    5: {
        1: 0.276 * ureg.Bq,
        2: 0.286 * ureg.Bq,
        3: 2.823 * ureg.Bq,
        4: 0.468 * ureg.Bq,
        "background": background,
    },
    6: {
        1: 0.265 * ureg.Bq,
        2: 0.281 * ureg.Bq,
        3: 1.252 * ureg.Bq,
        4: 0.333 * ureg.Bq,
        "background": background,
    },
    7: {
        1: 0.262 * ureg.Bq,
        2: 0.286 * ureg.Bq,
        3: 0.874 * ureg.Bq,
        4: 0.387 * ureg.Bq,
        "background": background,
    },
}

measurements_after_background_sub = substract_background_from_measurements(
    raw_measurements
)

# time starts at 12/05 9:30 AM
# 12/06 9:30 AM = 24 hours
# 12/07 9:30 AM = 48 hours
replacement_times = [
    # 12/05 22:58
    0 * ureg.day + 13 * ureg.hour + 28 * ureg.minute,
    # 12/06 09:09
    0 * ureg.day + 23 * ureg.hour + 39 * ureg.minute,
    # 12/06 21:58
    1 * ureg.day + 12 * ureg.hour + 28 * ureg.minute,
    # 12/07 11:19
    1 * ureg.day + 25 * ureg.hour + 49 * ureg.minute,
    # 12/08 09:53
    2 * ureg.day + 24 * ureg.hour + 23 * ureg.minute,
    # 12/09 08:28
    3 * ureg.day + 22 * ureg.hour + 58 * ureg.minute,
    # 12/11 12:48
    6 * ureg.day + 3 * ureg.hour + 18 * ureg.minute,
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

baby_model.k_top *= mass_transport_coeff_factor * 0.62
optimised_ratio = 4e-2
baby_model.k_wall = baby_model.k_top * optimised_ratio

exposure_time = 12 * ureg.hour

baby_model.irradiations = [
    [0 * ureg.hour, 0 + exposure_time],
    [24 * ureg.hour, 24 * ureg.hour + exposure_time],
]

# calculated from Kevin's activation foil analysis
P383_neutron_rate = 5.42e8 * ureg.neutron * ureg.s**-1
A325_neutron_rate = 2.33e8 * ureg.neutron * ureg.s**-1

neutron_rate_relative_uncertainty = 0.089
baby_model.neutron_rate = (
    P383_neutron_rate + A325_neutron_rate
) / 2  # the neutron rate is divided by two to acount for the double counting (two detectors)
