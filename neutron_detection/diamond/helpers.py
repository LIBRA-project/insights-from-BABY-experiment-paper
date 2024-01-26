import os
import numpy as np


def get_time_energy_values(directory):
    energy_values = []
    time_values = []

    # Iterate through all CSV files in the directory
    for filename in os.listdir(directory):
        if filename.endswith(".CSV"):
            # Load data from the CSV file
            csv_file_path = os.path.join(directory, filename)
            data = np.genfromtxt(csv_file_path, delimiter=";")
            # print(data.shape)
            # print(data)
            # Assuming the energy values are in a specific column (e.g., column 1)
            time_column = 2
            energy_column = 3

            # Extract and append the energy data to the list
            time_values.extend(data[:, time_column])
            energy_values.extend(data[:, energy_column])
    return time_values, energy_values


def main(directory, bin_time, energy_peak_min, energy_peak_max):
    time_values, energy_values = get_time_energy_values(directory)

    # Create a histogram to represent the combined energy spectrum

    # sort time and energy values
    inds = np.argsort(time_values)
    time_values = np.array(time_values)[inds]
    time_values *= 1 / 1e12  # ps to s

    bin_time = 100

    time_bins = np.arange(0, time_values[-2], bin_time)
    energy_values = np.array(energy_values)[inds]

    # Define the energy range of the (n,alpha) peak
    # energy_peak_min, energy_peak_max = 1180, 1300
    peak_mask = (energy_values > energy_peak_min) & (energy_values < energy_peak_max)
    peak_time_values = time_values[peak_mask]

    peak_count_rates, peak_count_rate_bins = np.histogram(
        peak_time_values, bins=time_bins
    )
    peak_count_rates = peak_count_rates / bin_time

    all_count_rates, all_count_rate_bins = np.histogram(time_values, bins=time_bins)
    all_count_rates = all_count_rates / bin_time

    res = {
        "all_count_rates": all_count_rates,
        "all_count_rate_bins": all_count_rate_bins,
        "time_values": time_values,
        "energy_values": energy_values,
        "peak_count_rates": peak_count_rates,
        "peak_count_rate_bins": peak_count_rate_bins,
        "peak_time_values": peak_time_values,
    }
    return res


def get_avg_neutron_rate(res, window):
    section_data = {}
    # Create mask to only count pulses of any energy in section time window
    time_values = res["time_values"]
    peak_time_values = res["peak_time_values"]

    tot_time_mask = np.logical_and(
        time_values > window[0],
        time_values < window[1],
    )

    # Create mask to only count pulses from (n,alpha) peak in section time window
    peak_time_mask = np.logical_and(
        peak_time_values > window[0],
        peak_time_values < window[1],
    )

    section_data["tot counts"] = len(time_values[tot_time_mask])
    section_data["tot err"] = np.sqrt(len(time_values[tot_time_mask]))
    window_range = np.diff(window)
    section_data["tot count rate"] = section_data["tot counts"] / window_range
    section_data["tot count rate err"] = section_data["tot err"] / window_range

    section_data["peak counts"] = len(peak_time_values[peak_time_mask])
    section_data["peak err"] = np.sqrt(len(peak_time_values[peak_time_mask]))
    section_data["peak count rate"] = section_data["peak counts"] / window_range
    section_data["peak count rate err"] = section_data["peak err"] / window_range

    print("\tSection: {}".format(section_data))
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
    return section_data


if __name__ == "__main__":
    pass
