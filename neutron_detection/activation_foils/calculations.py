import numpy as np
from settings import *
import datetime
from data import foil_data


def n93_number(foil_mass: float):
    n93_molar_mass = 92.90637 * ureg.g / ureg.mol
    return (foil_mass / n93_molar_mass).to(ureg.particle)


def delay_time(start_time: str, end_time: str):
    # convert string to time with datetime
    start_time = datetime.datetime.strptime(start_time, "%m/%d/%Y %H:%M:%S")
    end_time = datetime.datetime.strptime(end_time, "%m/%d/%Y %H:%M:%S")
    return (end_time - start_time).total_seconds() * ureg.s


def get_neutron_flux(experiment: dict):
    overall_efficiency = (
        (geometric_efficiency * nal_gamma_efficiency * branching_ratio)
        * ureg.count
        * ureg.particle**-1
    )
    flux = experiment["photon_counts"] / (
        overall_efficiency
        * n93_number(experiment["foil_mass"])
        * Nb93_n_2n_Nb92m_cross_section_at_14Mev
    )
    print(Nb93_n_2n_Nb92m_cross_section_at_14Mev)
    print(Nb92m_decay_constant)

    X = np.exp(-Nb92m_decay_constant * 12 * ureg.h)
    flux *= (1 - (1 - (1 - X) * X) * X) ** -1

    delta_t = delay_time(
        experiment["time_generator_off"], experiment["start_time_counting"]
    )
    flux *= (
        (
            np.exp(-Nb92m_decay_constant * (delta_t + 4 * ureg.h))
            - np.exp(-Nb92m_decay_constant * delta_t)
        )
    ) ** -1

    flux *= -Nb92m_decay_constant
    return flux.to(ureg["1/(cm^2 hr)"])


for name, data in foil_data.items():
    flux = get_neutron_flux(data)
    print(f"Neutron flux for {name}, foil {data['foil_name']} is {flux:.2e}")
