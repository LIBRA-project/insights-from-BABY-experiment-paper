import numpy as np
import matplotlib.pyplot as plt

from helpers import main, get_avg_neutron_rate


# Define the directory where your CSV files are located
data_directory = "raw_data"

res = main(data_directory, bin_time=100, energy_peak_min=1180, energy_peak_max=1300)
energy_values = res["energy_values"]
time_values = res["time_values"]
peak_count_rates = res["peak_count_rates"]
peak_count_rate_bins = res["peak_count_rate_bins"]
all_count_rates = res["all_count_rates"]
all_count_rate_bins = res["all_count_rate_bins"]

### Calculate average neutron rate for each generator
data = {
    "Day 1 Steady": {"window": [25810, 42970]},
    "Day 1 Jump": {"window": [43062, 50200]},
    "Day 2 Steady": {"window": [106870, 133300]},
}

for section_data in data.values():
    new_dict = get_avg_neutron_rate(
        res["time_values"], res["peak_time_values"], section_data["window"]
    )
    for key, value in new_dict.items():
        section_data[key] = value

fig1, ax1 = plt.subplots()
ax1.hist(energy_values, bins=300, histtype="step")
# ax1.set_title("Combined Energy Spectrum")
ax1.set_xlabel("Energy Channel")
ax1.set_ylabel("Counts")
print(time_values)
print(time_values.max())


fig2, ax2 = plt.subplots()
ax2.stairs(peak_count_rates, edges=peak_count_rate_bins)
ax2.set_xlabel("Time (s)")
ax2.set_ylabel(r"(n,$\alpha$) Count Rate (CPS)")


fig3, ax3 = plt.subplots()
ax3.stairs(all_count_rates, edges=all_count_rate_bins)
ax3.set_xlabel("Time (s)")
ax3.set_ylabel("Total Count Rate (CPS)")

## Plot relevant count rates in regions of interest
for section_data in data.values():
    ax3.hlines(
        section_data["tot count rate"][0],
        section_data["window"][0],
        section_data["window"][1],
        colors="k",
        linestyles="--",
    )
    t = ax3.text(
        np.mean(section_data["window"]),
        section_data["tot count rate"] + 7,
        "{:.0f} cps".format(section_data["tot count rate"][0]),
        ha="center",
        c="white",
    )
    t.set_bbox(dict(facecolor="black", alpha=1.0, edgecolor=None))

    ax2.hlines(
        section_data["peak count rate"][0],
        section_data["window"][0],
        section_data["window"][1],
        colors="k",
        linestyles="--",
    )
    t2 = ax2.text(
        np.mean(section_data["window"]),
        section_data["peak count rate"] + 0.2,
        "{:.1f} cps".format(section_data["peak count rate"][0]),
        ha="center",
        c="white",
    )
    t2.set_bbox(dict(facecolor="black", alpha=1.0, edgecolor=None))

# add 14 MeV peak annotation
# ax1.annotate(
#     "14 MeV peak",
#     # xy=(np.mean(energy_peak_min, energy_peak_max), 0.2),
#     xy=(1300, 0.2),
# )

# remove top and right axes for all figures
for ax in [ax1, ax2, ax3]:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(alpha=0.3)

for ext in ["pdf", "png", "svg"]:
    fig1.savefig("combined_energy_spectrum.{}".format(ext))
    fig2.savefig("n_alpha_peak_count_rate.{}".format(ext))
    fig3.savefig("overall_count_rate.{}".format(ext))

plt.show()
