import copy
import csv
import datetime
from typing import List

import numpy as np
import pandas as pd
import peakutils
from matplotlib import pyplot as plt
from scipy.integrate import simps as sp, trapz as tp

from SFG.spectrum.exceptions import CoverageCalculationImpossibleError
from SFG.spectrum.base_spectrum import BaseSpectrum


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
            return x * slope + intercept

        return baseline

    def correct_baseline(self):

        if self.baseline_corrected is None:
            temp = copy.deepcopy(self.y)

        else:
            # if the baseline correction already was performed, return immediately
            return

        if np.max(self.x) >= 3000:
            borders = (2750, 3000)
        else:
            borders = (2750, np.max(self.x))

        func = self.make_ch_baseline()

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

    def calculate_ch_integral(self, old_baseline=False) -> float:
        """Calculate the integral of the spectrum in the range of 2750-3000 wavenumbers"""

        # choose between the old and new style of baseline correction
        if old_baseline:
            self.correct_baseline()
        else:
            try:
                self.baseline_corrected = self.full_baseline_correction()
            except ValueError:
                print(f'{self} does not yield acceptable baseline by peakutils!')
                self.correct_baseline()

        # check upper border
        if max(self.x) >= 3000:
            upper = 3000
        else:
            upper = max(self.x)

        # check lower border
        if min(self.x) <= 2750:
            lower = 2750
        else:
            lower = min(self.x)

        borders = self.slice_by_borders(lower, upper)
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


class SfgAverager:
    # todo: throw an error and plot the spectra if the integral is NAN or zero!
    # todo: the benchmark function MUST display the integral value and the baseline
    # todo: replace print comments by professional logging
    """This class takes a list of SFG spectra and generates an average spectrum of them by interpolation and
    averaging. It is possible to pass a dictionary of date:dppc_integral key-value-pairs in order to calculate
    the coverage."""

    def __init__(self, spectra: List[SfgSpectrum], references=None, enforce_scale=False, name="default", debug=False,
                 baseline=False):
        self.failure_count = 0
        self.log = ""
        self.log += "Log file for averaging spectra\n"
        self.spectra = spectra
        self.references = references
        self.enforce_scale = enforce_scale
        self.name = name

        if len(self.spectra) == 0:
            print("Warning: zero spectra to average in SfgAverager!")
            self.average_spectrum = None
            self.integral = None
            self.coverage = None

        else:
            self.day_counter = {}

            self.average_spectrum = self.average_spectra(baseline=baseline)
            self.integral = self.average_spectrum.calculate_ch_integral()
            self.coverage = self.calc_coverage()

            if debug:
                if self.integral < 0:
                    self.benchmark()
                    print("Warning: negative integral value in SfgAverager!")
                    self.integral = 0
                    self.coverage = 0

    # todo: mit Gernot abklären ob von allen Baseline oder nur vom average
    def average_spectra(self, baseline=True):
        """Function performing the averaging: it ensures that all spectra are interpolated to have the same shape,
        then they are averaged. A AverageSpectrum  object is constructed and returned."""
        to_average = []

        # sort spectra by length of the wavenumber array (lambda)
        if self.enforce_scale is False:
            self.spectra.sort(key=lambda x: x.yield_wn_length(), reverse=True)
            root_x_scale = self.spectra[0].x
        else:
            root_x_scale = SfgAverager.enforce_base()

        # get y values by interpolation and collect the y values in a list
        # collect the dates of measurement for DPPC referencing
        for item in self.spectra:

            if 0 <= item.meta["time"].hour < 8:
                item.meta["time"] -= datetime.timedelta(days=1)

            date = item.meta["time"].date()

            if date not in self.day_counter:
                self.day_counter[date] = 1
            else:
                self.day_counter[date] += 1

            new_intensity = np.interp(root_x_scale, item.x, item.y)
            mask = (root_x_scale > np.max(item.x)) | (root_x_scale < np.min(item.x))
            new_intensity[mask] = np.nan
            to_average.append(new_intensity)

        to_average = np.array(to_average)
        average = np.nanmean(to_average, axis=0)
        std = np.nanstd(to_average, axis=0)

        # prepare meta data for average spectrum
        if self.name == "default":
            newname = self.spectra[0].name + "baseAV"
        else:
            newname = self.name
        in_new = [n.name for n in self.spectra]
        s_meta = {"name": newname, "made_from": in_new, "std": std}

        # with open("blabla.txt", "a") as outfile:
        # outfile.write(f'name: {newname} x: {root_x_scale}, y: {average}\n')

        s = AverageSpectrum(root_x_scale, average, s_meta)

        if baseline:
            try:
                s.y = s.full_baseline_correction()
            except (ValueError, ZeroDivisionError):
                if s.baseline_corrected:
                    s.y = s.baseline_corrected
                else:
                    s.correct_baseline()

        return s

    def calc_reference_part(self):
        """Calculate the participation of each DPPC references. This is important if the spectra to average are
        measured on different sampling days. If, for example,  5 samples are to average and 3 of them are measured
        on one day, 2 on another, the final coverage is calculated by dividing the AveragedSpectrum integral by the
        weighted sum of the DPPC integrals of the corresponding days, in our example (2/5 * DPPC_1) + (3/5 * DPPC_2)"""

        spec_number = len(self.spectra)
        total = 0
        self.log += f'Start of the reference calculating section: \n'

        for date in self.day_counter:
            # divide by total number of spectra in the average
            self.log += f'date {date} divided by the number of spectra {spec_number}\n'
            self.day_counter[date] /= spec_number

            # multiply the weighting factor by the integral of the day
            try:
                self.log += f"""Now multiplying the factor {self.day_counter[date]} 
                by the reference integral {self.references[date]}\n"""

                self.day_counter[date] *= self.references[date]
                total += self.day_counter[date]
            except KeyError:
                self.failure_count += 1
                self.log += f'Error: no suitable DPPC reference found four date {date}\n'

        self.log += f'Finalizing calculation. The total factor now is {total}.\n'

        return total

    def calc_coverage(self):
        """A convenience function  to calculate the surface coverage"""

        if self.references is not None:
            dppc_factor = self.calc_reference_part()
            coverage = np.sqrt(self.integral / dppc_factor)
            return coverage

        else:
            raise CoverageCalculationImpossibleError(
                f'Coverage not available for reference samples, integral is {self.integral}!')

    @DeprecationWarning
    def benchmark(self):
        self.create_log()
        l = [i for i in self.spectra]
        l.append(self.average_spectrum)
        p = DummyPlotter(l, save=True, savename=self.spectra[0].name, special="AV")
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

        # with open(name, "w") as outfile:
        # outfile.write(self.log)

    @staticmethod
    def enforce_base():
        reg1 = np.arange(2750, 3055, 5)
        reg2 = np.arange(3050, 3670, 20)
        reg3 = np.arange(3650, 3845, 5)
        new = np.concatenate((reg1, reg2, reg3), axis=None)
        return new


class AverageSpectrum(SfgSpectrum):

    def __init__(self, wavenumbers, intensities, meta):
        self.x = wavenumbers
        self.y = intensities
        self.x_unit = "wavenumber/ cm$^{-1}$"
        self.y_unit = "SFG intensity/ arb. u."
        self.name = meta["name"]

        self.meta = meta
        self.baseline_corrected = None
        self.regions = None

        # ensure nan-values in intensity and their corresponding wavenumbers are removed
        mask = np.isfinite(self.y)
        self.y = self.y[mask]
        self.x = self.x[mask]
        super().set_regions()


@DeprecationWarning
class DummyPlotter:
    """A test class to monitor the interaction of the subclasses of AbstractSpectrum with plotting routines."""

    # todo: remove from module
    def __init__(self, speclist, save=False, savedir="", savename="Default", special=None):
        self.speclist = speclist
        self.special = special
        self.save = save
        self.savedir = savedir
        self.savename = savename

    def plot_all(self, base=False, marker=True):
        for spectrum in self.speclist:

            if base is True:
                if isinstance(spectrum, AverageSpectrum):
                    func = spectrum.make_ch_baseline()
                    testx = np.linspace(2750, 3000, 1000)
                    testy = func(testx)
                    plt.plot(testx, testy, color="black")
                    integral = spectrum.calculate_ch_integral()
                    plt.title(str(round(integral, 6)))

            if self.special is None:
                if marker:
                    plt.plot(spectrum.x, spectrum.y, label=spectrum.name, marker="^", linestyle="-")
                else:
                    plt.plot(spectrum.x, spectrum.y, label=spectrum.name)

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
