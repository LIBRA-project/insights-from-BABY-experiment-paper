from simple_tritium_transport_model import ureg
import matplotlib.ticker as mticker
import math


def background_sub(measured, background):
    """Substracts the background of a measured activity.
    Returns zero if the background is greater than measurement.

    Args:
        measured (pint.Quantity): The measured activity
        background (pint.Quantity): the background acitivity

    Returns:
        pint.Quantity: activity with substracted background
    """
    if measured > background:
        return measured - background
    else:
        return 0 * ureg.Bq


def substract_background_from_measurements(raw_measurements):
    """Substracts the background from a set of measurements.

    Args:
        raw_measurements (dict): The raw measurements. The keys are the measurement
        numbers, and the values are dictionaries with the keys being the vial
        numbers and the values being the measured activities. The dictionary should
        also contain a key "background" with the background activity.

    Returns:
        dict: The measurements with the background substracted. The keys are the
        measurement numbers, and the values are dictionaries with the keys being the
        vial numbers and the values being the measured activities with the
        background substracted.
    """
    measurements_after_background_sub = {
        i: {
            j: background_sub(act, raw_measurements[i]["background"])
            for j, act in raw_measurements[i].items()
            if j != "background"
        }
        for i in raw_measurements
    }
    return measurements_after_background_sub


def cumulative_activity(measurements):
    """Calculates the cumulative activity of a set of measurements.

    Args:
        measurements (dict): The measurements. The keys are the measurement numbers,
        and the values are dictionaries with the keys being the vial numbers and the
        values being the measured activities.

    Returns:
        list: The cumulative activity of the measurements
    """
    for measurement in measurements.values():
        if "background" in measurement:
            raise ValueError("Background should be substracted first")

    total_samples = []
    for measurement in measurements.values():
        total_sample = sum(list(measurement.values()))
        total_samples.append(total_sample)

    cumulative_values = []
    for total_sample in total_samples:
        if cumulative_values:  # if list is not empty
            cumulative_values.append(cumulative_values[-1] + total_sample)
        else:
            cumulative_values.append(total_sample)

    return cumulative_values


class ScalarFormatterClass(mticker.ScalarFormatter):
    def format_data(self, value):
        # docstring inherited
        e = math.floor(math.log10(abs(value)))
        s = round(value / 10**e, 10)
        significand = self._format_maybe_minus_and_locale(
            "%d" if s % 1 == 0 else "%1.3g",
            s,  # here is the change! the format %1.3g is used to have 2 significant digits
        )
        if e == 0:
            return significand
        exponent = self._format_maybe_minus_and_locale("%d", e)
        if self._useMathText or self._usetex:
            exponent = "10^{%s}" % exponent
            return (
                exponent
                if s == 1  # reformat 1x10^y as 10^y
                else rf"{significand} \times {exponent}"
            )
        else:
            return f"{significand}e{exponent}"
