from abc import ABC, abstractmethod
import copy
import csv
import json
import numpy as np
import datetime
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator, MultipleLocator
import pandas as pd
import peakutils
from scipy import stats
from scipy.signal import savgol_filter
from scipy.integrate import simps as sp
from scipy.integrate import trapz as tp

from SFG.spectrum.exceptions import InvalidSpectrumError, IntegrationError


class MetaSpectrum(type):

    def __call__(cls, *args, **kwargs):
        temp = super().__call__(*args, **kwargs)

        for attr in ("x", "y", "x_unit", "y_unit", "name"):
            if not hasattr(temp, attr):
                raise InvalidSpectrumError(f'Tried to instantiate spectrum object without suitable attribute {attr}!')

        return temp


class BaseSpectrum(metaclass=MetaSpectrum):

    def __init__(self, name=None, x=None, y=None, x_unit=None, y_unit=None, timestamp=None):
        self.name = name
        self.x = x
        self.y = y
        self.x_unit = x_unit
        self.y_unit = y_unit
        self.timestamp = timestamp

    # magic methods

    def __repr__(self):
        return f'{type(self).__name__} Object with name "{self.name}"'

    # working on the data
    def yield_spectral_range(self):
        """returns a list containing maximum and minimum wavenumer and the number of data points"""
        return [min(self.x), max(self.x), len(self.x)]

    def get_xrange_indices(self, lower, upper):
        """Takes a high (upper) and a low (lower) target x value as argument. Returns
        the indices of the wavenumber array of the spectrum that are the borders of this interval."""
        lower_index = np.argmax(self.x >= lower)
        upper_index = np.argmax(self.x >= upper)
        return int(lower_index), int(upper_index)

    def get_xrange(self, lower, upper):
        # todo: ensure this functions work as well for y_values
        """Returns the slice of the x values in the borders of lower to upper"""
        lower, upper = self.get_xrange_indices()
        return self.x[lower, upper + 1]

    def normalize(self, external=None):
        """Normalize the spectrum's y data either to the maximum of the y values or an
        external factor"""
        if external is None:
            return np.max(self.y)
        else:
            return self.y/external

    def integrate_slice(self, x_array, y_array):
        """Integrates the y_array which has a spacing given by the x_array. First it tries to apply
        simpson rule integration rule, but if this fails the function invokes integration via
        trapeziodal rule"""
        try:
            area = sp(y_array, x_array)
            if np.isnan(area):
                print(f"""Integration failed in spectrum {self.name} using Simpson's rule. 
                Falling back to trapezoidal rule.""")
                area = tp(y_array, x_array)
            return area
        except:
            raise IntegrationError(f'Integration not possible for {self.name}')

    # export functions
    def properties_to_dict(self):
        temp = {
                "name": self.name,
                "x_unit": self.x_unit,
                "y_unit": self.y_unit,
                "x": self.x,
                "y": self.y,
                "timestamp": self.timestamp
                }
        return temp

    def to_pandas_dataframe(self) -> pd.DataFrame:
       pd.DataFrame(data=self.properties_to_dict())

    def to_csv(self):
        self.to_pandas_dataframe().to_csv(self.name + ".csv", index=False, sep=";")

    def to_json(self) -> str:
        temp = self.properties_to_dict()
        return json.dumps(temp)


class SfgSpectrum(BaseSpectrum):
    """The SFG spectrum class is the foundation of all analysis and plotting tools. It contains a class
    SystematicName (or a derived class) which carries most of the metainformation. Besides holding the
    experimental data, it gives access to a variety of functions like normalization, peak picking etc."""

    # magic methods
    def __init__(self, wavenumbers, intensity, ir_intensity, vis_intensity, meta):

        self.x = wavenumbers[::-1]
        self.raw_intensity = intensity[::-1]
        self.vis_intensity = vis_intensity[::-1]
        self.ir_intensity = ir_intensity[::-1]
        self.meta = meta
        self.y = self.raw_intensity / (self.vis_intensity * self.ir_intensity)
        self.baseline_corrected = None

        self.x_unit = "wavenumber/ cm$^{-1}$"
        self.y_unit = "SFG intensity/ arb. u."
        self.name = self.meta["name"]

        self.regions = None
        self.set_regions()

    def __lt__(self, SFG2) -> bool:
        """Returns true if the current spectrum was measured before SFG2"""
        if self.meta["creation_time"] < SFG2.name.meta["creation_time"]:
            return True
        else:
            return False

    # spectral data processing and analysis tools

    def normalize_to_highest(self, intensity="default", external_norm="none") -> np.ndarray:
        """normalize an given array to its maximum, typically the normalized or raw intensity"""
        # todo: this function is somehow strange
        if intensity == "default":
            intensity = self.y
        if external_norm == "none":
            norm_factor = np.max(intensity)
        else:
            norm_factor = external_norm

        return intensity / norm_factor

    def integrate_peak(self, x_array: np.ndarray, y_array: np.ndarray) -> float:
        """
        Numpy integration routine for numerical peak integration with the trapezoidal rule.

        :param x_array: the x values (spacing) of the curve to integrate
        :type x_array: array-like
        :param y_array: the y values of the curve to integrate
        :type y_array: array-like
        :return: the area under the x/y curve
        :rtype: float
        """
        try:
            area = sp(y_array, x_array)
            if np.isnan(area):
                area = tp(y_array, x_array)
            return area

        except:
            raise ValueError(f'Integration not possible for {self.name}')

    def root(self) -> np.ndarray:
        """
        :return: the spectrum's normalized intensity
        :rtype: np.ndarray
        """
        return np.sqrt(self.y)

    def yield_maximum(self) -> float:
        """
        :return: maximum intensity value of the spectrum
        :rtype: float
        """
        return np.max(self.y)

    def yield_peaklist(self, mode="norm"):

        out = []
        tup = self.detailed_analysis(threshold=1.5, intensity=mode)
        for peak in tup:
            out.append(peak[0])
        return out

    def yield_spectral_range(self):
        """returns a list containing maximum and minimum wavenumer and the number of data points"""
        return [min(self.x), max(self.x), len(self.x)]

    def yield_increment(self):
        """Calculates stepsize and wavenumbers where the stepsize is changed"""
        borders = []
        stepsize = []
        current = self.x[0]
        currentstep = abs(current - self.x[1])
        borders.append(current)

        for wavenumber in self.x[1:]:
            s = abs(wavenumber - current)
            if s != currentstep:
                stepsize.append(currentstep)
                borders.append(current)
                currentstep = s
                current = wavenumber
            else:
                current = wavenumber
        borders.append(current)
        stepsize.append(s)
        increment = [borders, stepsize]
        return increment

    def yield_wn_length(self):
        return np.max(self.x) - np.min(self.x)

    # info functions

    def drop_ascii(self) -> None:
        """Create an ascii file with the wavenumbers and normalized intensities"""
        with open(self._name + ".csv", "w") as outfile:
            writer = csv.writer(outfile, delimiter=";")
            for i in zip(self.x, self.y):
                writer.writerow((i[0], i[1]))

    def convert_to_export_dataframe(self) -> pd.DataFrame:
        """
        This function returns a Pandas dataframe, suitable for data export to
        origin and similar other programs

        :return: a pandas dataframe with the spectral data
        :rtype: pd.DataFrame
        """
        data = {
            "wavenumbers": self.x,
            "normalized intensity": self.y,
            "raw intensity": self.raw_intensity,
            "IR intensity": self.ir_intensity,
            "VIS intensity": self.vis_intensity
        }
        return pd.DataFrame(data=data)

    # CH baseline correction and integration

    # todo das hier ist ganz großer Mist und richtig error-prone
    def make_ch_baseline(self, debug=False):

        # todo: interchange high and low at the slice borders function
        if np.min(self.x) > 2800:
            left = self.slice_by_borders(np.min(self.x), 2815)
        else:
            left = self.slice_by_borders(np.min(self.x), 2800)

        if np.max(self.x) >= 3030:
            right = self.slice_by_borders(3000, 3030)
        else:
            right = self.slice_by_borders(self.x[-4], self.x[-1])

        left_x = self.x[left[0]:left[1] + 1]
        left_y = self.y[left[0]:left[1] + 1]

        right_x = self.x[right[0]:right[1] + 1]
        right_y = self.y[right[0]:right[1] + 1]

        slope = (np.average(right_y) - np.average(left_y)) / \
                (np.average(right_x) - np.average(left_x))

        intercept = np.average(left_y) - slope * np.average(left_x)

        if debug:
            print(f'intercept: {intercept}, slope: {slope}, left:{left}, right: {right}')

        def baseline(x):
            return x*slope+intercept

        return baseline

    def correct_baseline(self):

        if np.max(self.x) >= 3000:
            borders = (2750, 3000)
        else:
            borders = (2750, np.max(self.x))

        func = self.make_ch_baseline()

        if self.baseline_corrected is None:
            temp = copy.deepcopy(self.y)

        else:
            # if the baseline correction already was performed, return immediately
            return

        xvals = self.x.copy()
        corr = func(xvals)


        # ensure that only the region of the spec defined in the borders is used
        # xvals is a vector being 0 everywhere except in the to-correct area where it
        # is 1 so that xvals*corr yields nonzero in the defined regions only
        np.putmask(xvals, xvals < borders[0], 0)
        np.putmask(xvals, xvals > borders[1], 0)
        np.putmask(xvals, xvals != 0, 1)

        corr *= xvals

        # apply the correction
        temp -= corr

        self.baseline_corrected = temp

    def calculate_ch_integral(self):
        self.correct_baseline()
        if max(self.x) >= 3000:
            borders = self.slice_by_borders(np.min(self.x), 3000)
        else:
            borders = self.slice_by_borders(np.min(self.x), max(self.x))

        x_array = self.x[borders[0]:borders[1] + 1]
        y_array = self.baseline_corrected[borders[0]:borders[1] + 1]
        integral = self.integrate_peak(x_array, y_array)
        return integral

    def calc_region_integral(self, region):
        borders = self.regions[region]
        borders = self.slice_by_borders(borders[0], borders[1])
        x_array = self.x[borders[0]:borders[1] + 1]
        y_array = self.y[borders[0]:borders[1] + 1]
        try:
            integral = self.integrate_peak(x_array, y_array)
            if np.isnan(integral):
                print(f'x: {x_array}, y: {y_array}')
            return integral

        except ValueError:
            print(f'Integration not possible in {self.name} in region{region}')
            return np.nan

    # auxiliary function

    def create_pointlist(self, y_array):

        output = []
        for i, (a, b) in enumerate(zip(self.x[::-1], y_array)):
            output.append((a, b, i))

        return output

    def slice_by_borders(self, lower, upper):
        """Takes a high (upper) and a low (lower) reciprocal centimeter value as argument. Returns
        the indices of the wavenumber array of the spectrum that are the borders of this interval."""
        lower_index = np.argmax(self.x >= lower)
        upper_index = np.argmax(self.x >= upper)
        return int(lower_index), int(upper_index)

    def set_regions(self):
        self.regions = {"CH": (int(np.min(self.x)), 3000),
                        "dangling": (3670, 3760),
                        "OH": (3005, 3350), "OH2": (3350, 3670)}

    # new peak picking and baseline correction
    def full_baseline_correction(self):
        lower = self.get_xrange_indices(np.min(self.x), 3030)
        upper = self.get_xrange_indices(3030, np.max(self.x))

        left = self.y[lower[0]:lower[1] + 1]
        right = (self.y[upper[0] + 1:upper[1] + 1])

        base = peakutils.baseline(self.y[lower[0]:lower[1] + 1], deg=1)
        base2 = peakutils.baseline(self.y[upper[0] + 1:upper[1] + 1], deg=1)

        return np.concatenate([left - base, right - base2])


class LtIsotherm(BaseSpectrum):
    """A class to represent experimental Langmuir trough isotherms, handling time, area, area per molecule and
    surface pressure"""

    def __init__(self, name, measured_time, time, area, apm, pressure, lift_off=None, correct=True):
        self.name = name
        self.measured_time = measured_time
        self.time = time
        self.area = area
        self.apm = apm
        self.pressure = pressure
        self.compression_factor = self.calc_compression_factor()

        # todo: remove this hack which ensures the spectrum complies with the metaclass
        self.x = None
        self.y = None
        self.x_unit = "area/ cm$^{2}$"
        self.y_unit = "surface pressure/ mNm$^{-1}$"

        self.speed = None
        self.measurement_number = None

        self.lift_off = lift_off

        if correct is True:
            p = np.min(self.pressure)
            if p < 0:
                self.pressure += np.abs(p)

        self.setup_spec()

    def __lt__(self, other):
        if self.get_maximum_pressure() < other.get_maximum_pressure():
            return True

    def setup_spec(self) -> None:
        self._x = self.area
        self._y = self.pressure

    def drop_ascii(self) -> None:
        """Drops an ascii file with semikolon-separated data in the form time;area;surface pressure. Intention
        is easy interfacing with external software like Excel or Origin"""

        with open(self.name + ".out", "w") as outfile:
            for a, b, c in zip(self.time, self.area, self.pressure):
                outfile.write(str(a) + ";" + str(b) + ";" + str(c) + "\n")

    def convert_to_export_dataframe(self):
        """This function returns a Pandas dataframe, suitable for data export to
        origin and similar other programs"""
        data = {
            "time": self.time,
            "area": self.area,
            "area per molecule": self.apm,
            "surface pressure": self.pressure
        }
        return pd.DataFrame(data=data)

    def get_maximum_pressure(self, shrinked=None):
        """Returns the maximum measured surface pressure. Note: This property is uesd for the less-then
        operator implementation of this class!"""
        if shrinked == None:
            return np.max(self.pressure)
        else:
            try:
                return np.max(shrinked)
            except:
                # todo specify the type of error numpy will throw
                raise TypeError("Can not calc maximum for this operand")

    def calc_compression_factor(self):
        max = np.max(self.area)
        return self.area / max

    def derive_pressure(self):
        """Calculates the difference quotient of the surface pressure with respect to the area.
        Useful for calculation of surface elasticity"""
        return np.diff(self.pressure) / np.diff(self.area)

    def create_pointlist(self, x_array):
        """Returns a list containing the index, the x_array (usually area or time) and the surface pressure.
        This is used for example by the GUI functions to find the closest datapoint to a mouse-defined position
        in the plot"""

        output = []
        for i, (a, b) in enumerate(zip(x_array, self.pressure)):
            output.append((a, b, i))
        return output

    def get_slice(self, x_array, lower, upper, smooth=False):
        """Returns a slice of the x_array (usually time or area) defined by the lower and upper integer
        index."""

        x_out = x_array[lower:upper + 1]
        if smooth == False:
            y_out = self.pressure[lower:upper + 1]
        else:
            y_out = self.smooth()[lower:upper + 1]

        return x_out, y_out

    def calculate_elasticity(self) -> np.ndarray:
        """Returns the surface elasticity of the isotherm"""

        xdata = self.area[::-1]
        ydata = self.pressure[::-1]
        out = []

        for i in range(len(self.pressure) - 1):
            p = abs(ydata[i + 1] - ydata[i])
            a = abs((xdata[i + 1] - xdata[i]))
            A = (xdata[i + 1] - xdata[i]) / 2

            out.append(p / a * A)

        return np.array(out[::-1])

    def cut_away_decay(self, x_array):
        max_area = np.argmax(self.area < 31.2)
        return self.get_slice(x_array, 0, max_area)

    @staticmethod
    def get_closest_index(array_datapoints, check) -> int:
        # todo: replace with numpy argmax function
        d = 1000000000000
        index = None
        for point in array_datapoints:

            d_temp = np.sqrt((check[0] - point[0]) ** 2 + (check[1] - point[1]) ** 2)
            if d_temp < d:
                d = d_temp
                index = point[2]

        return index




