from settings import ureg

foil_data = {
    "20240128_BABY6_Nb-A_A-325_1054_4hr": {
        "foil_name": "A",
        "foil_mass": 0.2733 * ureg.g,
        "photon_counts": 76100.8 * ureg.count,
        "photon_counts_uncertainty": 416.398 * ureg.count,
        "time_generator_off": "1/26/2024 21:38:00",
        "start_time_counting": "1/28/2024 10:54:00",
        "counting_time": 4 * ureg.h,
        "distance_from_center_of_target_plane": 5.08 * ureg.cm,
    },
    "20240128_BABY6_Nb-A_A-325_1638_4hr": {
        "foil_name": "A",
        "foil_mass": 0.2733 * ureg.g,
        "photon_counts": 75449.6 * ureg.count,
        "photon_counts_uncertainty": 408.274 * ureg.count,
        "time_generator_off": "1/26/2024 21:38:00",
        "start_time_counting": "1/28/2024 16:38:00",
        "counting_time": 4 * ureg.h,
        "distance_from_center_of_target_plane": 5.08 * ureg.cm,
    },
    "20240127_BABY6_Nb-C_P383_1313_4hr": {
        "foil_name": "C",
        "foil_mass": 0.2869 * ureg.g,
        "photon_counts": 199604 * ureg.count,
        "photon_counts_uncertainty": 535.974 * ureg.count,
        "time_generator_off": "1/26/2024 21:38:00",
        "start_time_counting": "1/27/2024 13:13:00",
        "counting_time": 4 * ureg.h,
        "distance_from_center_of_target_plane": 5.08 * ureg.cm,
    },
    "20240127_BABY6_Nb-C_P383_1738_4hr": {
        "foil_name": "C",
        "foil_mass": 0.2869 * ureg.g,
        "photon_counts": 194856 * ureg.count,
        "photon_counts_uncertainty": 531.229 * ureg.count,
        "time_generator_off": "1/26/2024 21:38:00",
        "start_time_counting": "1/27/2024 17:38:00",
        "counting_time": 4 * ureg.h,
        "distance_from_center_of_target_plane": 5.08 * ureg.cm,
    },
}
