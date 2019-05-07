from abc import ABC, abstractmethod
import copy
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


class SfgSpectrum(AbstractSpectrum):
    """The SFG spectrum class is the foundation of all analysis and plotting tools. It contains a class
    SystematicName (or a derived class) which carries most of the metainformation. Besides holding the
    experimental data, it gives access to a variety of functions like normalization, peak picking etc."""

    # magic methods
    def __init__(self, wavenumbers, intensity, ir_intensity, vis_intensity, meta):

        # todo: invert the wavenunmber/ intensity scale to make it ascending

        self.wavenumbers = wavenumbers
        self.raw_intensity = intensity
        self.vis_intensity = vis_intensity
        self.ir_intensity = ir_intensity
        self.meta = meta
        self.normalized_intensity = self.raw_intensity / (self.vis_intensity * self.ir_intensity)
        self.baseline_corrected = None
        self.setup_spec()

        self.regions = None
        self.set_regions()

    def setup_spec(self):
        self._name = self.meta["name"]
        self._x = self.wavenumbers
        self._y = self.normalized_intensity
        self._x_unit = "wavenumber/ cm$^{-1}$"
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

    def yield_density(self):
        return np.max(self.wavenumbers) - np.min(self.wavenumbers)

    # info functions

    def drop_ascii(self):
        """Create an ascii file with the wavenumbers and normalized intensities"""
        with open(self._name + ".csv", "w") as outfile:
            writer = csv.writer(outfile, delimiter=";")
            for i in zip(self.wavenumbers, self.normalized_intensity):
                writer.writerow((i[0], i[1]))

    # CH baseline correction and integration

    def make_ch_baseline(self):

        # todo: interchange high and low at the slice borders function
        if np.min(self.wavenumbers) > 2800:
            left = self.slice_by_borders(2810, np.min(self.wavenumbers))
        else:
            left = self.slice_by_borders(2800, np.min(self.wavenumbers))

        right = self.slice_by_borders(3030, 3000)


        left_x = self.wavenumbers[left[0]:left[1] + 1]
        left_y = self.normalized_intensity[left[0]:left[1] + 1]

        right_x = self.wavenumbers[right[0]:right[1] + 1]
        right_y = self.normalized_intensity[right[0]:right[1] + 1]

        slope = (np.average(right_y) - np.average(left_y)) / \
                (np.average(right_x) - np.average(left_x))

        intercept = np.average(left_y) - slope * np.average(left_x)

        # x = list(left_x) + list(right_x)
        # y = list(left_y) + list(right_y)

        # slope, intercept, r, p, std = stats.linregress(x, y)

        baseline = lambda x: slope * x + intercept
        return baseline

    def correct_baseline(self):

        borders = (2750, 3000)

        func = self.make_ch_baseline()

        if self.baseline_corrected is None:
            temp = copy.deepcopy(self.normalized_intensity)

        else:
            # if the baseline correction already was performed, return immediately
            return

        xvals = self.wavenumbers.copy()
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

        # todo: this must be changed to be proper numpy functions!

        self.correct_baseline()
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
        the indices of the wavenumber array of the spectrum that are the borders of this interval."""

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

    def set_regions(self):
        self.regions = {"CH": (int(np.min(self.wavenumbers)), 3000),
                        "dangling": (3670, 3760),
                        "OH": (3005, 3350), "OH2": (3350, 3670)}


class AverageSpectrum(SfgSpectrum):

    def __init__(self, wavenumbers, intensities, meta):
        # todo: invert the wavenunmber/ intensity scale to make it ascending
        self.wavenumbers = wavenumbers
        self.normalized_intensity = intensities
        self.meta = meta
        self.baseline_corrected = None
        self.regions = None
        super().setup_spec()
        super().set_regions()


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
        return (self.area / max)

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
        index = int(np.where(self.pressure == max)[0][0])
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

    def __init__(self, speclist, save=False, savedir="", savename="Default", special=None):
        self.speclist = speclist
        self.special = special
        self.save = save
        self.savedir = savedir
        self.savename = savename

    def plot_all(self):
        for spectrum in self.speclist:
            if self.special is None:
                plt.plot(spectrum.x, spectrum.y, label=spectrum.name, marker="^", linestyle="-")
            else:
                if self.special not in spectrum.name:
                    plt.plot(spectrum.x, spectrum.y, label=spectrum.name, marker="^", alpha=0.3)
                else:
                    plt.plot(spectrum.x, spectrum.y, label=spectrum.name, marker="^", linestyle="-", color="r")

        plt.xlabel(spectrum.x_unit)
        plt.ylabel(spectrum.y_unit)
        plt.minorticks_on()
        plt.legend()

        if self.save is False:
            plt.show()

        else:
            path = self.savedir + "/" + self.savename + ".png"
            plt.savefig(path)
            plt.close()


class SfgAverager:

    def __init__(self, spectra, references=None):
        self.failure_count = 0
        self.log = ""
        self.log += "Log file for averaging spectra\n"

        self.spectra = spectra
        self.references = references

        self.day_counter = {}
        self.average_spectrum = self.average_spectra()

        self.integral = self.average_spectrum.calculate_ch_integral()
        self.coverage = self.calc_coverage()

        self.benchmark()

    def average_spectra(self):
        to_average = []

        # sort spectra by density (lambda)
        self.spectra.sort(key=lambda x: x.yield_density(), reverse=True)

        root_x_scale = self.spectra[0].x[::-1]

        # get y values by interpolation and collect the y values in a list
        # collect the dates of measurement for DPPC referencing
        for item in self.spectra:
            date = item.meta["time"].date()

            if date not in self.day_counter:
                self.day_counter[date] = 1
            else:
                self.day_counter[date] += 1

            new_intensity = np.interp(root_x_scale, item.x[::-1], item.y[::-1])
            to_average.append(new_intensity)

        to_average = np.array(to_average)
        average = np.mean(to_average, axis=0)
        std = np.std(to_average, axis=0)

        # prepare meta data for average spectrum
        newname = self.spectra[0].name + "baseAV"
        in_new = [n.name for n in self.spectra]
        s_meta = {"name": newname, "made_from": in_new, "std": std}
        s = AverageSpectrum(root_x_scale[::-1], average[::-1], s_meta)

        return s

    def calc_reference_part(self):

        spec_number = len(self.spectra)
        total = 0
        self.log += f'Start of the reference calculating section: \n'

        for date in self.day_counter:
            # divide by total number of spectra in the average
            self.log += f'date {date} divided by the number of spectra {spec_number}\n'
            self.day_counter[date] /= spec_number

            # multiply the weighting factor by the integral of the day
            try:
                self.log += f"""Now multiplying the factor{self.day_counter[date]} 
                by the reference integral {self.references[date]}\n"""

                self.day_counter[date] *= self.references[date]
                total += self.day_counter[date]
            except KeyError:
                self.failure_count += 1
                self.log += f'Error: no suitable DPPC reference found four date {date}\n'

        self.log += f'Finalizing calculation. The total factor now is {total}.\n'

        return total

    def calc_coverage(self):

        if self.references is not None:
            dppc_factor = self.calc_reference_part()
            coverage = np.sqrt(self.integral / dppc_factor)
            return coverage

        else:
            print("Coverage not available for reference samples!")

    def benchmark(self):
        self.create_log()
        l = [i for i in self.spectra]
        l.append(self.average_spectrum)
        p = DummyPlotter(l, save=True, savedir="benchmark", savename=self.spectra[0].name, special="AV")
        p.plot_all()

    def create_log(self):

        name = "benchmark/" + self.spectra[0].name + ".log"

        s = f'This average contains {len(self.spectra)} SFG spectra:\n'

        self.log += 80 * "-" + "\n"
        self.log += s
        for i in self.spectra:
            self.log += (i.name + "\n")

        self.log += 80 * "-" + "\n"
        s = f'integral: {self.integral}\ncoverage: {self.coverage}\n'

        with open(name, "w") as outfile:
            outfile.write(self.log)


"""TEST"""
# sorting by density - does it work?
