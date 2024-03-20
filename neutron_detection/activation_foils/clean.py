import sympy as sp
from data import foil_data
from settings import *
from calculations import n93_number, delay_time

t = sp.Symbol("t")
neutron_flux = sp.Symbol("\Gamma_n")
number_of_Nb93 = sp.Symbol(r"N_{Nb93}")
decay_constant = sp.Symbol("\lambda", positive=True)
cross_section = sp.Symbol("\sigma")
overall_efficiency = sp.Symbol("\eta")


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


# def get_neutron_flux(experiment, irradiations):
#     delay_time_val = delay_time(experiment, irradiations)

#     number_of_Nb92m_decays_measured = experiment["photon_counts"] / overall_efficiency

#     number_of_Nb92m_at_the_end_of_irradiation = get_number_of_Nb92m(irradiations)

#     flux_analytical = get_neutron_flux_sympy(experiment, irradiations)

#     flux = (
#         flux_analytical.subs(
#             number_of_Nb93,
#             n93_number(experiment["foil_mass"]).to(ureg.particle).magnitude,
#         )
#         .subs(decay_constant, Nb92m_decay_constant.to(ureg.s**-1).magnitude)
#         .subs(
#             cross_section,
#             Nb93_n_2n_Nb92m_cross_section_at_14Mev.to(ureg.cm**2).magnitude,
#         )
#         .subs(
#             overall_efficiency,
#             geometric_efficiency * nal_gamma_efficiency * branching_ratio,
#         )
#         .subs(delay_time, delay_time_val.to(ureg.s).magnitude)
#         .subs(counting_time, experiment["counting_time"].to(ureg.s).magnitude)
#         .subs(
#             sp.Symbol("counts", positive=True),
#             experiment["photon_counts"].magnitude,
#         )
#     )


def get_neutron_flux_sympy(experiment: dict, irradiations: list):
    """

    Args:
        experiment (dict): _description_

    Returns:
        _type_: _description_
    """
    delta_t = experiment["start_time_counting"] - experiment["time_generator_off"]

    counting_time = experiment["counting_time"]

    number_of_Nb92m_decays_measured = experiment["photon_counts"] / overall_efficiency

    number_of_Nb92m_at_the_end_of_irradiation = get_number_of_Nb92m(
        irradiations
    ).simplify()

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


if __name__ == "__main__":
    delta_t = sp.Symbol("\delta t")
    # irradiations = [
    #     {"t_on": 0, "t_off": delta_t},
    #     {"t_on": 2 * delta_t, "t_off": 3 * delta_t},
    # ]
    irradiations = [
        {"t_on": 0, "t_off": 12 * 3600},
        {"t_on": 24 * 3600, "t_off": 36 * 3600},
    ]
    experiment = {
        "photon_counts": sp.Symbol("counts", positive=True),
        "time_generator_off": 3 * delta_t,
        "start_time_counting": 3 * delta_t + sp.Symbol("delay_time"),
        "counting_time": sp.Symbol("t_{counting}"),
    }
    flux_analytic = get_neutron_flux_sympy(experiment, irradiations)

    flux_analytic = (
        (
            flux_analytic.subs(
                number_of_Nb93,
                n93_number(foil_data["20240128_BABY6_Nb-A_A-325_1054_4hr"]["foil_mass"])
                .to(ureg.particle)
                .magnitude,
            )
            .subs(decay_constant, Nb92m_decay_constant.to(ureg.s**-1).magnitude)
            .subs(
                cross_section,
                Nb93_n_2n_Nb92m_cross_section_at_14Mev.to(ureg.cm**2).magnitude,
            )
            .subs(
                overall_efficiency,
                geometric_efficiency * nal_gamma_efficiency * branching_ratio,
            )
            .subs(delta_t, (12 * ureg.h).to(ureg.s).magnitude)
            .subs(sp.Symbol("t_{counting}"), (4 * ureg.h).to(ureg.s).magnitude)
            .subs(
                sp.Symbol("delay_time"),
                (37.27 * ureg.h).to(ureg.s).magnitude,
            )
            .subs(
                sp.Symbol("counts", positive=True),
                foil_data["20240128_BABY6_Nb-A_A-325_1054_4hr"][
                    "photon_counts"
                ].magnitude,
            )
        )
        * ureg.cm**-2
        * ureg.s**-1
    )
    print(f"{flux_analytic.to(ureg.cm**-2 * ureg.h**-1):.2e}")
