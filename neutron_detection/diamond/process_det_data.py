import numpy as np
import matplotlib.pyplot as plt
import os

# Define the directory where your CSV files are located
run_num = 5
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
energy_values = np.array(energy_values)[inds]
print(time_values[0])
print(time_values[-1])

fig1, ax1 = plt.subplots(figsize=[10, 8])
ax1.hist(energy_values, bins=300, color="blue", alpha=1.0, histtype="step")
ax1.set_title("Combined Energy Spectrum")
ax1.set_xlabel("Energy Channel")
ax1.set_ylabel("Counts")
print(time_values)
print(time_values.max())

time_bins = np.arange(0, time_values[-2], bin_time)
all_count_rates, all_count_rate_bins = np.histogram(time_values, bins=time_bins)
all_count_rates = all_count_rates / bin_time

fig3, ax3 = plt.subplots(figsize=[10, 8])
ax3.stairs(all_count_rates, edges=all_count_rate_bins, color="blue", alpha=1.0)
ax3.set_title("Overall Count Rate")
ax3.set_xlabel("Time [s]")
ax3.set_ylabel("Count Rate [cps]")

peak_mask = (energy_values > 1180) & (energy_values < 1300)

peak_time_values = time_values[peak_mask]
print(peak_time_values.shape)
peak_energy_values = energy_values[peak_mask]

peak_count_rates, peak_count_rate_bins = np.histogram(peak_time_values, bins=time_bins)
peak_count_rates = peak_count_rates / bin_time

fig2, ax2 = plt.subplots(figsize=[10, 8])
ax2.stairs(peak_count_rates, edges=peak_count_rate_bins, color="blue", alpha=1.0)
ax2.set_title("(n,alpha) Peak Count Rate")
ax2.set_xlabel("Time [s]")
ax2.set_ylabel("Count Rate [cps]")


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

for gen in data.keys():
    print(gen)
    for section in data[gen].keys():
        # Create mask to only count pulses of any energy in section time window
        tot_time_mask = np.logical_and(
            time_values > data[gen][section]["window"][0],
            time_values < data[gen][section]["window"][1],
        )

        # Create mask to only count pulses from (n,alpha) peak in section time window
        peak_time_mask = np.logical_and(
            peak_time_values > data[gen][section]["window"][0],
            peak_time_values < data[gen][section]["window"][1],
        )

        data[gen][section]["tot counts"] = len(time_values[tot_time_mask])
        data[gen][section]["tot err"] = np.sqrt(len(time_values[tot_time_mask]))
        data[gen][section]["tot count rate"] = data[gen][section][
            "tot counts"
        ] / np.diff(data[gen][section]["window"])
        data[gen][section]["tot count rate err"] = data[gen][section][
            "tot err"
        ] / np.diff(data[gen][section]["window"])

        data[gen][section]["peak counts"] = len(peak_time_values[peak_time_mask])
        data[gen][section]["peak err"] = np.sqrt(len(peak_time_values[peak_time_mask]))
        data[gen][section]["peak count rate"] = data[gen][section][
            "peak counts"
        ] / np.diff(data[gen][section]["window"])
        data[gen][section]["peak count rate err"] = data[gen][section][
            "peak err"
        ] / np.diff(data[gen][section]["window"])

        print("\tSection: {}".format(section))
        print(
            "\t\tTotal Counts: {} +/- {}".format(
                data[gen][section]["tot counts"], data[gen][section]["tot err"]
            )
        )
        print(
            "\t\tTotal Count Rate: {} +/- {}".format(
                data[gen][section]["tot count rate"],
                data[gen][section]["tot count rate err"],
            )
        )
        print(
            "\t\tPeak Counts: {} +/- {}".format(
                data[gen][section]["peak counts"], data[gen][section]["peak err"]
            )
        )
        print(
            "\t\tPeak Count Rate: {} +/- {}".format(
                data[gen][section]["peak count rate"],
                data[gen][section]["peak count rate err"],
            )
        )


print(data["total"]["Day 1 Steady"]["window"])
print(data["total"]["Day 1 Steady"]["window"][0])
print(data["total"]["Day 1 Steady"]["tot count rate"])
## Plot relevant count rates in regions of interest
for section in data["total"].keys():
    ax3.plot(
        data["total"][section]["window"],
        [data["total"][section]["tot count rate"][0]] * 2,
        "--k",
    )
    t = ax3.text(
        np.mean(data["total"][section]["window"]),
        data["total"][section]["tot count rate"] + 7,
        "{:.0f} cps".format(data["total"][section]["tot count rate"][0]),
        ha="center",
        c="white",
    )
    t.set_bbox(dict(facecolor="black", alpha=1.0, edgecolor=None))

    ax2.plot(
        data["total"][section]["window"],
        [data["total"][section]["peak count rate"][0]] * 2,
        "--k",
    )
    t2 = ax2.text(
        np.mean(data["total"][section]["window"]),
        data["total"][section]["peak count rate"] + 0.2,
        "{:.1f} cps".format(data["total"][section]["peak count rate"][0]),
        ha="center",
        c="white",
    )
    t2.set_bbox(dict(facecolor="black", alpha=1.0, edgecolor=None))

data["A325"]["ave"] = {}
data["A325"]["ave"]["peak count rate"] = (
    data["A325"]["Shutdown 1"]["peak count rate"]
    + data["A325"]["Shutdown 2"]["peak count rate"]
) / 2
data["A325"]["ave"]["peak count rate err"] = (
    np.sqrt(
        data["A325"]["Shutdown 1"]["peak count rate err"] ** 2
        + data["A325"]["Shutdown 2"]["peak count rate err"] ** 2
    )
    / 2
)

print(
    "A325 Average (n,alpha) Count Rate: {} +/- {}".format(
        data["A325"]["ave"]["peak count rate"],
        data["A325"]["ave"]["peak count rate err"],
    )
)

# Get the Day 1 Steady average count rate by assuming that (total - A325) = P383
data["P383"] = {}
data["P383"]["ave"] = {}
data["P383"]["ave"]["peak count rate"] = (
    data["total"]["Day 1 Steady"]["peak count rate"]
    - data["A325"]["ave"]["peak count rate"]
)
data["P383"]["ave"]["peak count rate err"] = np.sqrt(
    data["A325"]["ave"]["peak count rate err"] ** 2
    + data["total"]["Day 1 Steady"]["peak count rate err"] ** 2
)

print(
    "P383 Average (n,alpha) Count Rate: {} +/- {}".format(
        data["P383"]["ave"]["peak count rate"],
        data["P383"]["ave"]["peak count rate err"],
    )
)

plt.show()
