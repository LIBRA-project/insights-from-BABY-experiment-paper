import pint
import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

ureg = pint.UnitRegistry()
ureg.setup_matplotlib()
ureg.define("neutron = 1 * particle = n")

SPECIFIC_ACT = 3.57e14 * ureg.Bq * ureg.g**-1
MOLAR_MASS = 6.032 / 2 * ureg.g * ureg.mol**-1


class Model:
    def __init__(
        self,
        radius: pint.Quantity = None,
        height: pint.Quantity = None,
        TBR: pint.Quantity = None,
    ) -> None:
        self.radius = radius
        self.height = height

        self.L_wall = 0.06 * ureg.inches

        self.neutron_rate = 3e8 * ureg.neutron * ureg.s**-1

        self.TBR = TBR
        self.irradiations = []

        self.k_wall = 1.9e-8 * ureg.m * ureg.s**-1  # from Kumagai
        self.k_top = 4.9e-7 * ureg.m * ureg.s**-1  # from Kumagai

        self.concentrations = []
        self.times = []

    @property
    def volume(self):
        return self.A_top * self.height

    @property
    def A_top(self):
        return np.pi * self.radius**2

    @property
    def A_wall(self):
        perimeter_wall = 2 * np.pi * (self.radius + self.L_wall)
        return perimeter_wall * self.height + self.A_top

    def source(self, t):
        for irradiation in self.irradiations:
            irradiation_start = irradiation[0]
            irradiation_stop = irradiation[1]
            if irradiation_start < t < irradiation_stop:
                return self.TBR * self.neutron_rate
        return 0 * self.TBR * self.neutron_rate

    def Q_wall(self, c_salt):
        return self.A_wall * self.k_wall * c_salt

    def Q_top(self, c_salt):
        return self.A_top * self.k_top * c_salt

    def rhs(self, t, c):
        t *= ureg.s
        c *= ureg.particle * ureg.m**-3

        return self.volume.to(ureg.m**3) ** -1 * (
            self.source(t).to(ureg.particle * ureg.s**-1)
            - self.Q_wall(c).to(ureg.particle * ureg.s**-1)
            - self.Q_top(c).to(ureg.particle * ureg.s**-1)
        )

    def run(self, t_final):
        concentration_units = ureg.particle * ureg.m**-3
        time_units = ureg.s
        initial_concentration = 0
        res = solve_ivp(
            fun=self.rhs,
            t_span=(0, t_final.to(time_units).magnitude),
            y0=[initial_concentration],
            t_eval=np.sort(
                np.concatenate(
                    [
                        np.linspace(0, t_final.to(time_units).magnitude, 1000),
                        [irr[1].to(time_units).magnitude for irr in self.irradiations],
                    ]
                )
            ),
            # method="RK45",  # RK45 doesn't catch the end of irradiations properly... unless constraining the max_step
            # max_step=(0.5 * ureg.h).to(time_units).magnitude,
            method="Radau",
        )
        self.times = res.t * time_units
        self.concentrations = res.y[0] * concentration_units

    def reset(self):
        self.concentrations = []
        self.times = []

    def integrated_release_top(self):
        top_release = self.Q_top(self.concentrations)
        integrated_top = cumulative_trapezoid(
            top_release.to(ureg.particle * ureg.h**-1).magnitude,
            self.times.to(ureg.h).magnitude,
            initial=0,
        )
        integrated_top *= ureg.particle  # attach units
        return integrated_top

    def integrated_release_wall(self):
        wall_release = self.Q_wall(self.concentrations)
        integrated_wall = cumulative_trapezoid(
            wall_release.to(ureg.particle * ureg.h**-1).magnitude,
            self.times.to(ureg.h).magnitude,
            initial=0,
        )
        integrated_wall *= ureg.particle  # attach units
        return integrated_wall


def quantity_to_activity(Q):
    return Q * SPECIFIC_ACT * MOLAR_MASS


def activity_to_quantity(A):
    return A / (SPECIFIC_ACT * MOLAR_MASS)


def plot_bars(measurements, index=None, bar_width=0.35, stacked=True):
    vial_1_vals = (
        np.array([sample[1].magnitude for sample in measurements.values()]) * ureg.Bq
    )
    vial_2_vals = (
        np.array([sample[2].magnitude for sample in measurements.values()]) * ureg.Bq
    )
    vial_3_vals = (
        np.array([sample[3].magnitude for sample in measurements.values()]) * ureg.Bq
    )
    vial_4_vals = (
        np.array([sample[4].magnitude for sample in measurements.values()]) * ureg.Bq
    )

    if index is None:
        if stacked:
            index = np.arange(len(measurements))
        else:
            group_spacing = 1  # Adjust this value to control spacing between groups
            index = (
                np.arange(len(measurements)) * (group_spacing / 2 + 1) * bar_width * 4
            )

    if stacked:
        vial_3_bar = plt.bar(
            index,
            vial_3_vals,
            bar_width,
            label="Vial 3",
            color="#FB8500",
        )
        vial_4_bar = plt.bar(
            index,
            vial_4_vals,
            bar_width,
            label="Vial 4",
            color="#FFB703",
            bottom=vial_3_vals,
        )
        vial_1_bar = plt.bar(
            index,
            vial_1_vals,
            bar_width,
            label="Vial 1",
            color="#219EBC",
            bottom=vial_3_vals + vial_4_vals,
        )
        vial_2_bar = plt.bar(
            index,
            vial_2_vals,
            bar_width,
            label="Vial 2",
            color="#8ECAE6",
            bottom=vial_3_vals + vial_4_vals + vial_1_vals,
        )
    else:
        vial_1_bar = plt.bar(
            index - 1.5 * bar_width,
            vial_1_vals,
            bar_width,
            linewidth=2,
            edgecolor="white",
            label="Vial 1",
            color="#219EBC",
        )
        vial_2_bar = plt.bar(
            index - 0.5 * bar_width,
            vial_2_vals,
            bar_width,
            linewidth=2,
            edgecolor="white",
            label="Vial 2",
            color="#8ECAE6",
        )
        vial_3_bar = plt.bar(
            index + 0.5 * bar_width,
            vial_3_vals,
            bar_width,
            linewidth=2,
            edgecolor="white",
            label="Vial 3",
            color="#FB8500",
        )
        vial_4_bar = plt.bar(
            index + 1.5 * bar_width,
            vial_4_vals,
            bar_width,
            linewidth=2,
            edgecolor="white",
            label="Vial 4",
            color="#FFB703",
        )

    return index


def replace_water(sample_activity, time, replacement_times):
    sample_activity_changed = np.copy(sample_activity)
    times_changed = np.copy(time)

    for replacement_time in sorted(replacement_times):
        indices = np.where(times_changed > replacement_time)
        # at each replacement, make the sample activity drop to zero
        sample_activity_changed[indices] -= sample_activity_changed[indices][0]

        # insert nan value to induce a line break in plots
        if indices[0].size > 0:
            first_index = indices[0][0]
            sample_activity_changed = np.insert(
                sample_activity_changed, first_index, np.nan * ureg.Bq
            )
            times_changed = np.insert(times_changed, first_index, np.nan * ureg.day)

    return sample_activity_changed, times_changed


COLLECTION_VOLUME = 10 * ureg.ml
LSC_SAMPLE_VOLUME = 10 * ureg.ml


def plot_sample_activity_top(
    model: Model,
    replacement_times,
    collection_vol=COLLECTION_VOLUME,
    lsc_sample_vol=LSC_SAMPLE_VOLUME,
    **kwargs
):
    integrated_top = quantity_to_activity(model.integrated_release_top()).to(ureg.Bq)
    sample_activity_top = integrated_top / collection_vol * lsc_sample_vol
    times = model.times
    sample_activity_top, times = replace_water(
        sample_activity_top, times, replacement_times
    )
    l = plt.plot(times.to(ureg.day), sample_activity_top, **kwargs)
    return l


def plot_sample_activity_wall(
    model: Model,
    replacement_times,
    collection_vol=COLLECTION_VOLUME,
    lsc_sample_vol=LSC_SAMPLE_VOLUME,
    **kwargs
):
    integrated_wall = quantity_to_activity(model.integrated_release_wall()).to(ureg.Bq)
    sample_activity_wall = integrated_wall / collection_vol * lsc_sample_vol
    times = model.times
    sample_activity_wall, times = replace_water(
        sample_activity_wall, times, replacement_times
    )
    l = plt.plot(times.to(ureg.day), sample_activity_wall, **kwargs)
    return l


def plot_salt_inventory(model: Model, **kwargs):
    salt_inventory = (quantity_to_activity(model.concentrations) * model.volume).to(
        ureg.Bq
    )
    l = plt.plot(model.times.to(ureg.day), salt_inventory, **kwargs)
    return l


def plot_top_release(model: Model, **kwargs):
    top_release = model.Q_top(model.concentrations)
    top_release = quantity_to_activity(top_release).to(ureg.Bq * ureg.h**-1)
    l = plt.plot(model.times.to(ureg.day), top_release, **kwargs)
    return l


def plot_wall_release(model: Model, **kwargs):
    wall_release = model.Q_wall(model.concentrations)
    wall_release = quantity_to_activity(wall_release).to(ureg.Bq * ureg.h**-1)
    l = plt.plot(model.times.to(ureg.day), wall_release, **kwargs)
    return l


def plot_integrated_top_release(model: Model, **kwargs):
    integrated_top = quantity_to_activity(model.integrated_release_top()).to(ureg.Bq)
    sample_activity_top = integrated_top / COLLECTION_VOLUME * LSC_SAMPLE_VOLUME
    l = plt.plot(model.times.to(ureg.day), sample_activity_top, **kwargs)
    return l


def plot_integrated_wall_release(model: Model, **kwargs):
    integrated_wall = quantity_to_activity(model.integrated_release_wall()).to(ureg.Bq)
    sample_activity_wall = integrated_wall / COLLECTION_VOLUME * LSC_SAMPLE_VOLUME
    l = plt.plot(model.times.to(ureg.day), sample_activity_wall, **kwargs)
    return l


def plot_irradiation(model: Model, **kwargs):
    pols = []
    for irr in model.irradiations:
        pol = plt.axvspan(irr[0].to(ureg.day), irr[1].to(ureg.day), **kwargs)
        pols.append(pol)
    return pols
