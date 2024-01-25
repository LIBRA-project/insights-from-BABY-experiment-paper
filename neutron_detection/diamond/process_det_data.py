import numpy as np
import matplotlib.pyplot as plt
import os

# Define the directory where your CSV files are located
data_directory = "raw_data"
bin_time = 100

# Initialize an empty list to store energy values
energy_values = []
time_values = []


# Iterate through all CSV files in the directory
for filename in os.listdir(data_directory):
    if filename.endswith(".CSV"):
        # Load data from the CSV file
        csv_file_path = os.path.join(data_directory, filename)
        data = np.genfromtxt(csv_file_path, delimiter=";")
        # print(data.shape)
        # print(data)
        # Assuming the energy values are in a specific column (e.g., column 1)
        time_column = 2
        energy_column = 3

        # Extract and append the energy data to the list
        time_values.extend(data[:, time_column])
        energy_values.extend(data[:, energy_column])

# Create a histogram to represent the combined energy spectrum
# print(energy_values)
inds = np.argsort(time_values)
time_values = np.array(time_values)[inds] / 1e12
time_bins = np.arange(0, time_values[-2], bin_time)
energy_values = np.array(energy_values)[inds]
print(time_values[0])
print(time_values[-1])

peak_mask = (energy_values > 1180) & (energy_values < 1300)

peak_time_values = time_values[peak_mask]
print(peak_time_values.shape)
peak_energy_values = energy_values[peak_mask]


peak_count_rates, peak_count_rate_bins = np.histogram(peak_time_values, bins=time_bins)
peak_count_rates = peak_count_rates / bin_time


### Calculate average neutron rate for each generator
data = {
    "A325": {
        "Startup": {"window": [2550, 2710]},
        "Shutdown 1": {"window": [16600, 17300]},
        "Shutdown 2": {"window": [23850, 24500]},
    },
    "total": {
        "Day 1 Steady": {"window": [25810, 42970]},
        "Day 1 Jump": {"window": [43062, 50200]},
        "Day 2 Steady": {"window": [106870, 133300]},
    },
}

for gen, gen_data in data.items():
    print(gen)
    for section, section_data in gen_data.items():
        # Create mask to only count pulses of any energy in section time window
        tot_time_mask = np.logical_and(
            time_values > section_data["window"][0],
            time_values < section_data["window"][1],
        )

        # Create mask to only count pulses from (n,alpha) peak in section time window
        peak_time_mask = np.logical_and(
            peak_time_values > section_data["window"][0],
            peak_time_values < section_data["window"][1],
        )

        section_data["tot counts"] = len(time_values[tot_time_mask])
        section_data["tot err"] = np.sqrt(len(time_values[tot_time_mask]))
        section_data["tot count rate"] = section_data["tot counts"] / np.diff(
            section_data["window"]
        )
        section_data["tot count rate err"] = section_data["tot err"] / np.diff(
            section_data["window"]
        )

        section_data["peak counts"] = len(peak_time_values[peak_time_mask])
        section_data["peak err"] = np.sqrt(len(peak_time_values[peak_time_mask]))
        section_data["peak count rate"] = section_data["peak counts"] / np.diff(
            section_data["window"]
        )
        section_data["peak count rate err"] = section_data["peak err"] / np.diff(
            section_data["window"]
        )

        print("\tSection: {}".format(section))
        print(
            "\t\tTotal Counts: {} +/- {}".format(
                section_data["tot counts"], section_data["tot err"]
            )
        )
        print(
            "\t\tTotal Count Rate: {} +/- {}".format(
                section_data["tot count rate"],
                section_data["tot count rate err"],
            )
        )
        print(
            "\t\tPeak Counts: {} +/- {}".format(
                section_data["peak counts"], section_data["peak err"]
            )
        )
        print(
            "\t\tPeak Count Rate: {} +/- {}".format(
                section_data["peak count rate"],
                section_data["peak count rate err"],
            )
        )


print(data["total"]["Day 1 Steady"]["window"])
print(data["total"]["Day 1 Steady"]["window"][0])
print(data["total"]["Day 1 Steady"]["tot count rate"])


data["A325"]["ave"] = {
    "peak count rate": (
        data["A325"]["Shutdown 1"]["peak count rate"]
        + data["A325"]["Shutdown 2"]["peak count rate"]
    )
    / 2,
    "peak count rate err": np.sqrt(
        data["A325"]["Shutdown 1"]["peak count rate err"] ** 2
        + data["A325"]["Shutdown 2"]["peak count rate err"] ** 2
    )
    / 2,
}

print(
    "A325 Average (n,alpha) Count Rate: {} +/- {}".format(
        data["A325"]["ave"]["peak count rate"],
        data["A325"]["ave"]["peak count rate err"],
    )
)

# Get the Day 1 Steady average count rate by assuming that (total - A325) = P383
data["P383"] = {
    "ave": {
        "peak count rate": data["total"]["Day 1 Steady"]["peak count rate"]
        - data["A325"]["ave"]["peak count rate"],
        "peak count rate err": np.sqrt(
            data["A325"]["ave"]["peak count rate err"] ** 2
            + data["total"]["Day 1 Steady"]["peak count rate err"] ** 2
        ),
    }
}

print(
    "P383 Average (n,alpha) Count Rate: {} +/- {}".format(
        data["P383"]["ave"]["peak count rate"],
        data["P383"]["ave"]["peak count rate err"],
    )
)


fig1, ax1 = plt.subplots()
ax1.hist(energy_values, bins=300, histtype="step")
ax1.set_title("Combined Energy Spectrum")
ax1.set_xlabel("Energy Channel")
ax1.set_ylabel("Counts")
print(time_values)
print(time_values.max())

all_count_rates, all_count_rate_bins = np.histogram(time_values, bins=time_bins)
all_count_rates = all_count_rates / bin_time

fig2, ax2 = plt.subplots()
ax2.stairs(peak_count_rates, edges=peak_count_rate_bins)
ax2.set_title(r"(n,$\alpha$) Peak Count Rate")
ax2.set_xlabel("Time [s]")
ax2.set_ylabel("Count Rate [cps]")


fig3, ax3 = plt.subplots()
ax3.stairs(all_count_rates, edges=all_count_rate_bins)
ax3.set_title("Overall Count Rate")
ax3.set_xlabel("Time [s]")
ax3.set_ylabel("Count Rate [cps]")

## Plot relevant count rates in regions of interest
for section_data in data["total"].values():
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

# remove top and right axes for all figures
for ax in [ax1, ax2, ax3]:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(alpha=0.3)

for ext in ["pdf", "png", "svg"]:
    fig1.savefig("combined_energy_spectrum.{}".format(ext))
    fig2.savefig("peak_count_rate.{}".format(ext))
    fig3.savefig("overall_count_rate.{}".format(ext))

plt.show()
