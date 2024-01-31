import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from helpers import main, get_avg_neutron_rate


# Define the directory where your CSV files are located
data_runs = {}
for run, directory in zip([5], ["run_5"]):
    res = main(
        directory=directory, bin_time=1000, energy_peak_min=1180, energy_peak_max=1300
    )
    data_runs[run] = res

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
for run, res in data_runs.items():
    ax1.hist(
        res["energy_values"],
        bins=300,
        histtype="stepfilled",
        label="Run {}".format(run),
    )
# ax1.set_title("Combined Energy Spectrum")
if len(data_runs) > 1:
    ax1.legend()
ax1.set_xlabel("Energy Channel")
ax1.set_ylabel("Counts")
plt.gca().yaxis.set_major_formatter(
    mticker.ScalarFormatter(useOffset=False, useMathText=True)
)

fig2, ax2 = plt.subplots()
for run, res in data_runs.items():
    t_0_idx = np.where(res["peak_count_rates"] > 1)[0][0]
    t_0 = res["peak_count_rate_bins"][t_0_idx]
    bins_rescaled = res["peak_count_rate_bins"] - t_0
    ax2.stairs(
        res["peak_count_rates"],
        edges=bins_rescaled,
        label="Run {}".format(run),
    )
if len(data_runs) > 1:
    ax2.legend()
ax2.set_xlabel("Time (s)")
ax2.set_ylabel(r"(n,$\alpha$) Count Rate (CPS)")


fig3, ax3 = plt.subplots()
for run, res in data_runs.items():
    t_0_idx = np.where(res["all_count_rates"] > 1)[0][0]
    t_0 = res["all_count_rate_bins"][t_0_idx]
    bins_rescaled = res["all_count_rate_bins"] - t_0
    ax3.stairs(
        res["all_count_rates"],
        edges=bins_rescaled,
        label="Run {}".format(run),
    )
if len(data_runs) > 1:
    ax3.legend()
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
ax1.annotate(
    r"(n,$\alpha$) peak",
    xy=(1260, 7e4),
    xytext=(1300, 0.2e6),
    arrowprops=dict(arrowstyle="->"),
)

# remove top and right axes for all figures
for ax in [ax1, ax2, ax3]:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(alpha=0.3)
    ax.set_axisbelow(True)


for ext in ["pdf", "png", "svg"]:
    fig1.savefig("combined_energy_spectrum.{}".format(ext))
    fig2.savefig("n_alpha_peak_count_rate.{}".format(ext))
    fig3.savefig("overall_count_rate.{}".format(ext))

plt.show()
