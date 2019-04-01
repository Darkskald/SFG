from abc import ABC, abstractmethod
import csv
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from scipy.signal import savgol_filter
from scipy.integrate import simps as sp


class AbstractSpectrum(ABC):
    """The baseclass which ensures that all necessary properties for plotting
    are available as well as some basic functionality all spectra have in common."""

    @abstractmethod
    def setup_spec(self, x, y, x_unit=None, y_unit=None):
        """This method has to make sure that x and y values and the corresponding units are set properly"""
        pass
    @property
    def x(self):
        return self._x
    @x.setter
    def set_x(self, value):
        self._x = value

    @property
    def x_unit(self):
        return self._x_unit

    @x.setter
    def set_x_unit(self, value):
        self._x_unit = value

    @property
    def y(self):
        return self._y
    @y.setter
    def set_y(self, value):
        self._x = value

    @property
    def name(self):
        return self._name

    @x.setter
    def set_name(self, value):
        self._name = value

    @property
    def y_unit(self):
        return self._y_unit

    @x.setter
    def set_x_unit(self, value):
        self._y_unit = value

    def __repr__(self):
        return self._name

    def __str__(self):
        return self._name

    def smooth(self, a=9, b=3):
        """Performs a smooth operation of the measured pressure involving a Savitzky-Golay-filter"""
        return savgol_filter(self._y, a, b)

    # have to be overriden by the inheriting classes
    def drop_ascii(self):
        pass


# todo: ensure correct creation time injection and meta dictionary creation
class SfgSpectrum(AbstractSpectrum):
    """The SFG spectrum class is the foundation of all analysis and plotting tools. It contains a class
    SystematicName (or a derived class) which carries most of the metainformation. Besides holding the
    experimental data, it gives access to a variety of functions like normalization, peak picking etc."""

    # magic methods
    def __init__(self, wavenumbers, intensity, ir_intensity, vis_intensity, meta):
        self.wavenumbers = wavenumbers
        self.raw_intensity = intensity
        self.vis_intensity = vis_intensity
        self.ir_intensity = ir_intensity
        self.meta = meta
        self.normalized_intensity = self.raw_intensity / (self.vis_intensity * self.ir_intensity)
        self.baseline_corrected = None

        self.setup_spec()

    def setup_spec(self):
        self._name = self.meta["name"]
        self._x = self.wavenumbers
        self._y = self.normalized_intensity
        self._x_unit ="wavenumber/ cm$^{-1}$"
        self._y_unit = "SHG intensity/ arb. u."

    def __lt__(self, SFG2):
        """Returns true if the current spectrum was measured before SFG2"""
        if self.meta["creation_time"] < SFG2.name.meta["creation_time"]:
            return True
        else:
            return False

    # spectral data processing and analysis tools

    def normalize_to_highest(self, intensity="default", external_norm="none"):
        """normalize an given array to its maximum, typically the normalized or raw intensity"""
        if intensity == "default":
            intensity = self.normalized_intensity
        if external_norm == "none":
            norm_factor = np.max(intensity)
        else:
            norm_factor = external_norm

        return (intensity / norm_factor)

    def integrate_peak(self, x_array, y_array):
        """Numpy integration routine for numerical peak integration with the trapezoidal rule"""
        try:
            area = sp(y_array, x_array)
            return area
        except:
            return "Area  could not be calculated"

    def root(self):
        return np.sqrt(self.normalized_intensity)

    def yield_maximum(self):
        return np.max(self.normalized_intensity)

    def yield_peaklist(self, mode="norm"):

        out = []
        tup = self.detailed_analysis(threshold=1.5, intensity=mode)
        for peak in tup:
            out.append(peak[0])
        return out

    def yield_spectral_range(self):
        """returns a list containing maximum and minimum wavenumer and the number of data points"""
        return [min(self.wavenumbers), max(self.wavenumbers), len(self.wavenumbers)]

    def yield_increment(self):
        """Calculates stepsize and wavenumbers where the stepsize is changed"""
        borders = []
        stepsize = []
        current = self.wavenumbers[0]
        currentstep = abs(current - self.wavenumbers[1])
        borders.append(current)

        for wavenumber in self.wavenumbers[1:]:
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

    # info functions

    def drop_ascii(self):
        """Create an ascii file with the wavenumbers and normalized intensities"""
        with open(self._name + ".csv", "w") as outfile:
            writer = csv.writer(outfile, delimiter=";")
            for i in zip(self.wavenumbers, self.normalized_intensity):
                writer.writerow((i[0], i[1]))

    # CH baseline correction and integration

    def make_ch_baseline(self, average="min"):

        if np.min(self.wavenumbers > 2760):
            l_interval = self.slice_by_borders(2805, np.min(self.wavenumbers))

        else:
            l_interval = self.slice_by_borders(2800, 2750)

        l_interval_wl = self.wavenumbers[l_interval[0]:l_interval[1] + 1]
        l_interval = self.normalized_intensity[l_interval[0]:l_interval[1] + 1]

        interval = self.slice_by_borders(2960, 2895)
        interval_wl = self.wavenumbers[interval[0]:interval[1] + 1]
        interval = self.normalized_intensity[interval[0]:interval[1] + 1]

        min_index = np.argmin(interval)
        l_min_index = np.argmin(l_interval)

        if average == "min":
            slope = (interval[min_index] - l_interval[l_min_index]) / (
                        interval_wl[min_index] - l_interval_wl[l_min_index])
            intercept = l_interval[l_min_index] - slope * l_interval_wl[l_min_index]

        elif average == "min_reg":
            y2 = []
            x2 = []

            q = np.sort(l_interval)

            for i in range(3):
                y2.append(q[i])
                index = int((np.where(l_interval == q[i]))[0])
                x2.append(l_interval_wl[index])

            q = np.sort(interval)

            for i in range(3):
                y2.append(q[i])
                index = int((np.where(interval == q[i]))[0])
                x2.append(interval_wl[index])

            slope, intercept, r_value, p_value, std_err = stats.linregress(x2, y2)

        elif average == "gernot":

            if np.min(self.wavenumbers) > 2760:
                left = self.slice_by_borders(2760, 2750)
            else:
                left = self.slice_by_borders(2805, 2800)

            right = self.slice_by_borders(3010, 2990)

            left_x = self.wavenumbers[left[0]:left[1] + 1]
            left_y = self.normalized_intensity[left[0]:left[1] + 1]

            right_x = self.wavenumbers[right[0]:right[1] + 1]
            right_y = self.normalized_intensity[right[0]:right[1] + 1]

            slope = (np.average(right_y) - np.average(left_y)) / \
                    (np.average(right_x) - np.average(left_x))

            intercept = np.average(left_y) - slope * np.average(left_x)

        baseline = lambda x: slope * x + intercept
        return baseline

    def correct_baseline(self, average="min"):

        func = self.make_ch_baseline(average=average)
        temp = copy.deepcopy(self.normalized_intensity)

        for i in range(2750, 3000):
            index = np.where(self.wavenumbers == i)
            correction = func(self.wavenumbers[index])
            temp[index] = temp[index] - correction

        self.baseline_corrected = temp

    def calculate_ch_integral(self, average="min"):

        self.correct_baseline(average=average)
        borders = self.slice_by_borders(3000, np.min(self.wavenumbers))
        x_array = self.wavenumbers[borders[0]:borders[1] + 1]
        y_array = self.baseline_corrected[borders[0]:borders[1] + 1]
        integral = self.integrate_peak(x_array[::-1], y_array[::-1])
        return integral

    # auxiliary functions
    def create_pointlist(self, y_array):

        output = []
        for i, (a, b) in enumerate(zip(self.wavenumbers[::-1], y_array)):
            output.append((a, b, i))

        return output

    def slice_by_borders(self, upper, lower):
        """Takes a high (upper) and a low (lower) reciprocal centimeter value as argument. Returns
        the indices o the wavenumber array of the spectrum that are the borders of this intervall."""

        diff = 10000000
        upper_index = 0
        lower_index = -1

        for index, spectrum in enumerate(self.wavenumbers):

            temp_diff = abs(upper - self.wavenumbers[index])
            if temp_diff < diff:
                diff = temp_diff
                upper_index = index

        diff = 10000000

        for index, spectrum in enumerate(self.wavenumbers):

            temp_diff = abs(lower - self.wavenumbers[index])
            if temp_diff < diff:
                diff = temp_diff
                lower_index = index

        return upper_index, lower_index


class LtIsotherm(AbstractSpectrum):
    """A class to represent experimental Langmuir trough isotherms, handling time, area, area per molecule and
    surface pressure"""

    def __init__(self, name, measured_time, time, area, apm, pressure, lift_off=None, correct=True):
        self._name = name
        self.measured_time = measured_time
        self.time = time
        self.area = area
        self.apm = apm
        self.pressure = pressure
        self.compression_factor = self.calc_compression_factor()

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

    def setup_spec(self):
        self._x = self.area
        self._y = self.pressure
        self._x_unit = "area/ cm$^{2}$"
        self._y_unit = "surface pressure/ mNm$^{-1}$"

    def drop_ascii(self):
        """Drops an ascii file with semikolon-separated data in the form time;area;surface pressure. Intention
        is easy interfacing with external software like Excel or Origin"""

        with open(self.name + ".out", "w") as outfile:
            for a, b, c in zip(self.time, self.area, self.pressure):
                outfile.write(str(a) + ";" + str(b) + ";" + str(c) + "\n")

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
        return (self.area/max)

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

    def calculate_elasticity(self):
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

        max = self.get_maximum_pressure()
        index = int(np.where(self.pressure==max)[0][0])
        return self.get_slice(x_array, 0, index)

    @staticmethod
    def get_closest_index(array_datapoints, check):

        d = 1000000000000
        index = None
        for point in array_datapoints:

            d_temp = np.sqrt((check[0] - point[0]) ** 2 + (check[1] - point[1]) ** 2)
            if d_temp < d:
                d = d_temp
                index = point[2]

        return index


class DummyPlotter:
    """A test class to monitor the interaction of the subclasses of AbstractSpectrum with plotting routines."""
    def __init__(self, speclist):
        self.speclist = speclist

    def plot_all(self):
        for spectrum in self.speclist:
            plt.plot(spectrum.x, spectrum.y, label=spectrum.name)

        plt.xlabel(spectrum.x_unit)
        plt.ylabel(spectrum.y_unit)
        plt.legend()
        plt.show()


