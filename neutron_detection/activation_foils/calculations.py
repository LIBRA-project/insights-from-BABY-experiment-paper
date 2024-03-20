import numpy as np
import sympy as sp
from settings import *
import datetime
import matplotlib.pyplot as plt

t = sp.Symbol("t")
neutron_flux = sp.Symbol("\Gamma_n")
number_of_Nb93 = sp.Symbol(r"N_{Nb93}")
decay_constant = sp.Symbol("\lambda", positive=True)
cross_section = sp.Symbol("\sigma")
overall_efficiency = sp.Symbol("\eta")


def n93_number(foil_mass: float):
    n93_molar_mass = 92.90637 * ureg.g / ureg.mol
    return (foil_mass / n93_molar_mass).to(ureg.particle)


def delay_time(start_time: str, end_time: str):
    # convert string to time with datetime
    start_time = datetime.datetime.strptime(start_time, "%m/%d/%Y %H:%M:%S")
    end_time = datetime.datetime.strptime(end_time, "%m/%d/%Y %H:%M:%S")
    return (end_time - start_time).total_seconds() * ureg.s


def N_during_irradiation(
    N0, t, decay_constant, neutron_flux, cross_section, nb_Nb93, t_0, mod=sp
):
    A = nb_Nb93 * neutron_flux * cross_section / decay_constant
    return A + (N0 - A) * mod.exp(-decay_constant * (t - t_0))


def N_during_rest(N0, t, decay_constant, t_0, mod=sp):
    return N0 * mod.exp(-decay_constant * (t - t_0))


def get_number_of_Nb92m_after_last_irradiations(irradiations: list):
    """
    Returns the number of Nb92m at the end of the last irradiation
    as a symbolic expression with sympy


    Args:
        irradiations (list): list of dictionaries with keys "t_on" and "t_off"

    Returns:
        sp.Expression: symbolic expression for the number of Nb92m at the end of the last irradiation
    """

    # there is a series of irradiation and rest, get the number of Nb92m at the end of the last irradiation
    t_on = irradiations[0]["t_on"]
    t_off = irradiations[0]["t_off"]
    N_after_irradiation = N_during_irradiation(
        0, t_off, decay_constant, neutron_flux, cross_section, number_of_Nb93, t_on
    )
    for irr in irradiations[1:]:
        N_at_end_of_rest = N_during_rest(
            N_after_irradiation, irr["t_on"], decay_constant, t_0=t_off
        )
        N_after_irradiation = N_during_irradiation(
            N0=N_at_end_of_rest,
            t=irr["t_off"],
            decay_constant=decay_constant,
            neutron_flux=neutron_flux,
            cross_section=cross_section,
            nb_Nb93=number_of_Nb93,
            t_0=irr["t_on"],
        )
        t_off = irr["t_off"]

    return N_after_irradiation


def get_neutron_flux_sympy(experiment: dict, irradiations: list):
    """
    Returns the neutron flux as a symbolic expression with sympy

    Args:
        experiment (dict): dictionary containing the experiment data
        irradiations (list): list of dictionaries with keys "t_on" and "t_off"

    Returns:
        sp.Expression: symbolic expression for the neutron flux
    """
    delta_t = sp.Symbol("delay_time")

    counting_time = experiment["counting_time"]

    number_of_Nb92m_decays_measured = experiment["photon_counts"] / overall_efficiency

    number_of_Nb92m_at_the_end_of_irradiation = (
        get_number_of_Nb92m_after_last_irradiations(irradiations).simplify()
    )

    number_of_Nb92m_after_irradiation = N_during_rest(
        N0=number_of_Nb92m_at_the_end_of_irradiation,
        t=t,
        decay_constant=decay_constant,
        t_0=irradiations[-1]["t_off"],
    ).simplify()

    number_of_Nb92m_decays_analytical = sp.integrate(
        number_of_Nb92m_after_irradiation * decay_constant,
        (
            t,
            irradiations[-1]["t_off"] + delta_t,
            irradiations[-1]["t_off"] + delta_t + counting_time,
        ),
    ).simplify()

    neutron_flux_eq = sp.solve(
        sp.Eq(number_of_Nb92m_decays_measured, number_of_Nb92m_decays_analytical),
        neutron_flux,
    )[0]

    return neutron_flux_eq.simplify()


def get_number_ofNb92m_numpy(
    irradiations: list, times: np.ndarray, neutron_flux: float, nb_Nb93: float
):
    """
    Returns the number of Nb92m at different times as a numpy array

    Args:
        irradiations (list): list of dictionaries with keys "t_on" and "t_off"
        times (np.ndarray): array of times
        neutron_flux (float): neutron flux
        nb_Nb93 (float): number of Nb93 atoms in the foil

    Returns:
        np.ndarray: array of the number of Nb92m at different times
    """
    N = np.zeros_like(times)
    N[0] = 0
    decay_constant = Nb92m_decay_constant.to(1 / ureg.s).magnitude
    cross_section = Nb93_n_2n_Nb92m_cross_section_at_14Mev.to(ureg.cm**2).magnitude
    t_on = irradiations[0]["t_on"]
    t_off = irradiations[0]["t_off"]
    previous_irradiation_off = t_off
    idx = np.where(
        np.logical_and(
            (times >= irradiations[0]["t_on"]), (times <= irradiations[0]["t_off"])
        )
    )
    N[idx] = N_during_irradiation(
        0,
        times[idx],
        decay_constant,
        neutron_flux,
        cross_section,
        nb_Nb93,
        t_on,
        mod=np,
    )
    current_N = N[idx][-1]

    for irr in irradiations[1:]:

        idx_rest = np.where((times > previous_irradiation_off) & (times <= irr["t_on"]))
        N[idx_rest] = N_during_rest(
            current_N,
            times[idx_rest],
            decay_constant,
            t_0=previous_irradiation_off,
            mod=np,
        )

        idx_irr = np.where(
            np.logical_and((times > irr["t_on"]), (times <= irr["t_off"]))
        )
        N[idx_irr] = N_during_irradiation(
            N[idx_rest][-1],
            times[idx_irr],
            decay_constant,
            neutron_flux,
            cross_section,
            nb_Nb93,
            t_0=irr["t_on"],
            mod=np,
        )
        previous_irradiation_off = irr["t_off"]
        current_N = N[idx_irr][-1]
    return N


if __name__ == "__main__":

    irradiations = [
        {"t_on": 0, "t_off": 12 * 3600},
        {"t_on": 24 * 3600, "t_off": 36 * 3600},
    ]

    times = np.linspace(0, 36 * 3600, 1000)
    neutron_flux_val = 3.51e9 * ureg["1/(cm^2 h)"]  # 1/(cm^2 h)
    nb_Nb93 = 1.78e21
    N = get_number_ofNb92m_numpy(
        irradiations, times, neutron_flux_val.to(ureg["1/(cm^2 s)"]).magnitude, nb_Nb93
    )
    plt.plot(times / 3600, N)

    N_analytical = get_number_of_Nb92m_after_last_irradiations(irradiations)

    N_analytical = (
        N_analytical.subs(
            cross_section,
            Nb93_n_2n_Nb92m_cross_section_at_14Mev.to(ureg.cm**2).magnitude,
        )
        .subs(decay_constant, Nb92m_decay_constant.to(1 / ureg.s).magnitude)
        .subs(neutron_flux, neutron_flux_val.to(ureg["1/(cm^2 s)"]).magnitude)
        .subs(number_of_Nb93, nb_Nb93)
    )

    plt.scatter([times[-1] / 3600], [N_analytical], color="red")
    plt.xlabel("Time (h)")
    plt.ylabel("Number of Nb92m")

    for irr in irradiations:
        plt.gca().axvspan(
            irr["t_on"] / 3600, irr["t_off"] / 3600, alpha=0.2, color="red"
        )
    plt.show()
