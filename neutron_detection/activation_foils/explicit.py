from data import foil_data
from settings import *
from calculations import n93_number, delay_time
import numpy as np


def get_chain(irradiations):
    """
    Returns the value of
    (1 - exp(-\lambda * \Delta t_1)) * (1 - exp(-\lambda * \Delta t_2)) * ... * (1 - exp(-\lambda * \Delta t_n))
    where \Delta t_i is the time between the end of the i-th period (rest or irradiation) and the start of the next one

    Args:
        irradiations (list): list of dictionaries with keys "t_on" and "t_off" for irradiations

    Returns:
        float or pint.Quantity: the value of the chain
    """
    result = 1
    periods = [{"start": irradiations[0]["t_on"], "end": irradiations[0]["t_off"]}]
    for irr in irradiations[1:]:
        periods.append({"start": periods[-1]["end"], "end": irr["t_on"]})
        periods.append({"start": irr["t_on"], "end": irr["t_off"]})

    for period in periods:
        delta_t = period["end"] - period["start"]
        result = 1 - result * np.exp(-Nb92m_decay_constant * delta_t)
    return result


def get_neutron_flux(experiment: dict, irradiations: list):
    """_summary_

    Args:
        experiment (dict): _description_
        irradiations (list): _description_

    Returns:
        _type_: _description_
    """
    counting_time = experiment["counting_time"]
    overall_efficiency = (
        (geometric_efficiency * nal_gamma_efficiency * branching_ratio)
        * ureg.count
        / ureg.particle
    )
    number_of_Nb92m_decays_measured = experiment["photon_counts"] / overall_efficiency
    flux = (
        number_of_Nb92m_decays_measured
        / n93_number(experiment["foil_mass"])
        / Nb93_n_2n_Nb92m_cross_section_at_14Mev
    )

    flux *= get_chain(irradiations) ** -1
    time_between_generator_off_and_start_of_counting = delay_time(
        experiment["time_generator_off"], experiment["start_time_counting"]
    )
    flux *= (
        -1
        / Nb92m_decay_constant
        * (
            np.exp(
                -Nb92m_decay_constant
                * (time_between_generator_off_and_start_of_counting + counting_time)
            )
            - np.exp(
                -Nb92m_decay_constant * time_between_generator_off_and_start_of_counting
            )
        )
    ) ** -1

    return flux


if __name__ == "__main__":
    irradiations = [
        {"t_on": 0, "t_off": 12 * ureg.h},
        {"t_on": 24 * ureg.h, "t_off": 36 * ureg.h},
    ]
    for key, value in foil_data.items():
        print(key)
        print(
            f"{get_neutron_flux(value, irradiations).to(ureg.cm**-2 * ureg.h**-1):.2e}"
        )
        print()
