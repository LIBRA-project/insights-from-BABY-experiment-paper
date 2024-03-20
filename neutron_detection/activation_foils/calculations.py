import numpy as np
import sympy as sp
from settings import *
import datetime
from data import foil_data
import matplotlib.pyplot as plt


def n93_number(foil_mass: float):
    n93_molar_mass = 92.90637 * ureg.g / ureg.mol
    return (foil_mass / n93_molar_mass).to(ureg.particle)


def delay_time(start_time: str, end_time: str):
    # convert string to time with datetime
    start_time = datetime.datetime.strptime(start_time, "%m/%d/%Y %H:%M:%S")
    end_time = datetime.datetime.strptime(end_time, "%m/%d/%Y %H:%M:%S")
    return (end_time - start_time).total_seconds() * ureg.s


def get_neutron_flux(experiment: dict):
    """This assumes 12h irradiation + 12h rest + 12h irradiation

    Args:
        experiment (dict): _description_

    Returns:
        _type_: _description_
    """
    overall_efficiency = (
        (geometric_efficiency * nal_gamma_efficiency * branching_ratio)
        * ureg.count
        * ureg.particle**-1
    )
    number_of_Nb92m_measured = experiment["photon_counts"] / (overall_efficiency)
    flux = number_of_Nb92m_measured / (
        n93_number(experiment["foil_mass"]) * Nb93_n_2n_Nb92m_cross_section_at_14Mev
    )

    X = np.exp(-Nb92m_decay_constant * 12 * ureg.h)
    flux *= (1 - (1 - (1 - X) * X) * X) ** -1

    delta_t = delay_time(
        experiment["time_generator_off"], experiment["start_time_counting"]
    )
    counting_time = 4 * ureg.h
    flux *= (
        (
            np.exp(-Nb92m_decay_constant * (delta_t + counting_time))
            - np.exp(-Nb92m_decay_constant * delta_t)
        )
    ) ** -1

    flux *= -Nb92m_decay_constant
    return flux.to(ureg["1/(cm^2 hr)"])


def get_neutron_flux_generic(experiment: dict):
    """

    Args:
        experiment (dict): _description_

    Returns:
        _type_: _description_
    """
    irradiations = [
        {"t_on": 0, "t_off": 12 * 3600},
        {"t_on": 24 * 3600, "t_off": 36 * 3600},
    ]
    delta_t = delay_time(
        experiment["time_generator_off"], experiment["start_time_counting"]
    )
    counting_time = 4 * ureg.h
    overall_efficiency = (
        (geometric_efficiency * nal_gamma_efficiency * branching_ratio)
        * ureg.count
        * ureg.particle**-1
    )
    activity_measured = (
        experiment["photon_counts"] / overall_efficiency
    ) / counting_time
    number_of_Nb92m_measured = activity_measured / Nb92m_decay_constant

    number_of_Nb92m_after_irradiation = get_number_of_Nb92m(irradiations)

    number_of_Nb92m_after_counting = N_during_rest(
        N0=number_of_Nb92m_after_irradiation,
        t=irradiations[-1]["t_off"]
        + delta_t.to(ureg.s).magnitude
        + counting_time.to(ureg.s).magnitude,
        decay_constant=Nb92m_decay_constant.to(1 / ureg.s).magnitude,
        t_0=irradiations[-1]["t_off"],
    )

    neutron_flux = sp.solve(
        sp.Eq(sp.Symbol("N"), number_of_Nb92m_after_counting),
        sp.Symbol("\Gamma_n"),
    )[0]
    neutron_flux = (
        neutron_flux.subs(
            sp.Symbol("\lambda"), Nb92m_decay_constant.to(1 / ureg.s).magnitude
        )
        .subs(
            sp.Symbol("N_{Nb93}"),
            n93_number(experiment["foil_mass"]).to(ureg.particle).magnitude,
        )
        .subs(
            sp.Symbol("\sigma"),
            Nb93_n_2n_Nb92m_cross_section_at_14Mev.to(ureg.cm**2).magnitude,
        )
        .subs(sp.Symbol("N"), number_of_Nb92m_measured.to(ureg.particle).magnitude)
    )
    return neutron_flux * ureg["1/(cm^2 s)"]


def N_during_irradiation(
    N0, t, decay_constant, neutron_flux, cross_section, nb_Nb93, t_0, mod=sp
):
    A = nb_Nb93 * neutron_flux * cross_section / decay_constant
    return A + (N0 - A) * mod.exp(-decay_constant * (t - t_0))


def N_during_rest(N0, t, decay_constant, t_0, mod=sp):
    return N0 * mod.exp(-decay_constant * (t - t_0))


def get_number_of_Nb92m(irradiations: list):
    """
    Returns the number of Nb92m at the end of the last irradiation
    as a symbolic expression with sympy


    Args:
        irradiations (list): list of dictionaries with keys "t_on" and "t_off"

    Returns:
        sp.Expression: symbolic expression for the number of Nb92m at the end of the last irradiation
    """

    neutron_flux = sp.Symbol("\Gamma_n")
    number_of_Nb93 = sp.Symbol(r"N_{Nb93}")
    decay_constant = sp.Symbol("\lambda")
    cross_section = sp.Symbol("\sigma")

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


def get_number_ofNb92m_numpy(
    irradiations: list, times: np.ndarray, neutron_flux: float, nb_Nb93: float
):
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
    for name, data in foil_data.items():
        flux = get_neutron_flux(data)
        print(f"Neutron flux for {name}, foil {data['foil_name']} is {flux:.2e}")
        flux_analytic = get_neutron_flux_generic(data).to(ureg["1/(cm^2 hr)"])
        print(
            f"Neutron flux for {name}, foil {data['foil_name']} is {flux_analytic:.2e}"
        )

    irradiations = [
        # {"t_on": 0, "t_off": sp.Symbol("t_off1")},
        {"t_on": 0, "t_off": sp.Symbol("t_off1")},
        {"t_on": sp.Symbol("t_on2"), "t_off": sp.Symbol("t_off2")},
        # {"t_on": sp.Symbol("t_on3"), "t_off": sp.Symbol("t_off3")},
        # {"t_on": sp.Symbol("t_on4"), "t_off": sp.Symbol("t_off4")},
    ]
    analytical = get_number_of_Nb92m(irradiations)
    # print(get_number_of_Nb92m(irradiations).simplify())
    # print(
    #     sp.solve(
    #         sp.Eq(sp.Symbol("N"), get_number_of_Nb92m(irradiations)),
    #         sp.Symbol("\Gamma_n"),
    #     )
    # )

    irradiations = [
        {"t_on": 0, "t_off": 12 * 3600},
        {"t_on": 24 * 3600, "t_off": 36 * 3600},
    ]
    times = np.linspace(0, 36 * 3600, 1000)
    neutron_flux = 3.51e9 * ureg["1/(cm^2 h)"]  # 1/(cm^2 h)
    nb_Nb93 = 1.78e21
    N = get_number_ofNb92m_numpy(
        irradiations, times, neutron_flux.to(ureg["1/(cm^2 s)"]).magnitude, nb_Nb93
    )
    plt.plot(times / 3600, N)

    N_analytical = (
        analytical.subs(
            sp.Symbol("\sigma"),
            Nb93_n_2n_Nb92m_cross_section_at_14Mev.to(ureg.cm**2).magnitude,
        )
        .subs(sp.Symbol("\lambda"), Nb92m_decay_constant.to(1 / ureg.s).magnitude)
        .subs(sp.Symbol("\Gamma_n"), neutron_flux.to(ureg["1/(cm^2 s)"]).magnitude)
        .subs(sp.Symbol("N_{Nb93}"), nb_Nb93)
        .subs(sp.Symbol("t_on1"), 0)
        .subs(sp.Symbol("t_off1"), 12 * 3600)
        .subs(sp.Symbol("t_on2"), 24 * 3600)
        .subs(sp.Symbol("t_off2"), 36 * 3600)
    )

    plt.scatter([times[-1] / 3600], [N_analytical], color="red")
    plt.xlabel("Time (h)")
    plt.ylabel("Number of Nb92m")

    for irr in irradiations:
        plt.gca().axvspan(
            irr["t_on"] / 3600, irr["t_off"] / 3600, alpha=0.2, color="red"
        )
    plt.show()
