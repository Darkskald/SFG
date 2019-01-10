from sfg import SystematicName
from langmuir import LtIsotherm

# standard utilities
import os
import shutil
import csv
import time
import datetime
import copy
import traceback
import logging

# scientific libraries
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.lines import Line2D

import sqlite3
from scipy.signal import savgol_filter
from scipy.integrate import simps as sp
from scipy import stats

class SystematicGasExName(SystematicName):
    """Special modification of the systematic name in order to fit the requirements for the GasEx
    cruise samples from June and September 2018"""

    def __init__(self, namestring, creation_time="unknown"):

        self.full_name = namestring
        self.name = self.full_name[:-4]
        self.creation_time = creation_time

        temp = self.name.split("_")

        if len(temp) > 3:
            self.measurement_date = temp[0]
            self.station = temp[1] + "_" + temp[2]
            self.type = temp[3]

            try:
                self.number = temp[4]

            except IndexError:
                self.number = 1


        else:
            self.name.date = temp[0]
            self.station = temp[1] + "_" + temp[2]
            self.type = temp[3][0]

        self.station_hash = temp[1] + "_" + temp[2][1]
        self.sample_hash_string = "_".join(temp[1:])

    def __str__(self):

        return self.name

    def __repr__(self):
        return self.__str__()


class GasExLtIsotherm(LtIsotherm):

    def __init__(self, *args):
        super().__init__(*args)
        self.sample_hash = None
        self.process_name()

    def process_name(self):
        """Function to extract metainformation from the filename"""
        temp = self.name.split("_")

        hashstring = temp[1] + "_" + temp[2] + "_" + temp[3]

        if "c" not in temp[3][0]:
            hashstring += ("_" + temp[4])

        self.sample_hash = SampleHash(hashstring)

        if len(temp) >= 5:
            self.measurement_number = temp[-1]
            self.speed = temp[-2]

        else:
            self.measurement_number = 1
        if type(self.sample_hash) != SampleHash:
            raise ValueError(f'{self} was not created correctly!')

    def __str__(self):
        return self.name + " GasExLtIsotherm Object"

    def __repr__(self):
        return self.name + " GasExLtIsotherm Object"





class SampleHash:

    def __init__(self, namestring):
        self.namestring = namestring
        self.process_list = namestring.split("_")

        if len(self.process_list) not in (3, 4):
            print(namestring)
            #raise IndexError("Invalid process list length!")

        self.station_type = None
        self.station_number = None

        self.sample_type = None
        self.sample_number = None

        self.date = None
        self.station_hash = None
        self.station_construct = None

        self.date_from_name()
        self.set_station_hash()
        self.get_type()


    def __eq__(self, other):

        if self.namestring == other.namestring:
            return True
        else:
            return False

    def date_from_name(self, year=2018):
        month = int(self.process_list[0][0:2])
        day = int(self.process_list[0][2:])
        self.date = datetime.date(year, month, day)

    def set_station_hash(self):
        try:
            self.station_hash = self.process_list[0] + self.process_list[1][1]
        except IndexError:
            print(f'Invalid sample name {self.namestring}')

    def get_type(self):
        self.station_type = self.process_list[1][0].lower()
        self.station_number = self.process_list[1][1]
        self.sample_type = self.process_list[2]

        if len(self.process_list) == 3:
            self.sample_number = 1
        elif len(self.process_list) == 4:
            self.sample_number = self.process_list[-1]

        if "s" in self.sample_type:
            self.sample_type = "s"

        elif self.sample_type in ("deep", "low"):
            self.sample_type = "c"

    def get_doy(self):

        return self.date.timetuple().tm_yday


class Sample:

    def __init__(self, sample_hash):
        self.sample_hash = sample_hash
        self.station_hash = self.sample_hash.station_hash  # todo:
        self.lt_isotherms = []  # list of corresponding Isotherms
        self.sfg_spectra = []  # list of sfg spectra

        self.ch_integral = None
        self.max_pressure = None
        self.surface_tension = None

    def __repr__(self):
        pass

    def __str__(self):
        pass


class Station:
    """A class carrying all information and data for a given cruise station, especially SFG and isotherms"""

    def __init__(self, name, type=None, parent=None):

        self.station_hash = name
        self.date = datetime.date(2018, int(name[0:2]), int(name[2:4]))
        self.station_number = name[4]
        self.type = type

        self.sfg_spectra = []
        self.lt_isotherms = []
        self.tensions = []

        self.samples = []

        self.lt_joined = {}

        self.parent = parent
        self.stats = {

            "tension_average": None,
            "tension_deep": None,
            "tension_sml": None,
            "tension_plate": None,
            "tension_screen": None,
            "pressure_average": None,
            "pressure_deep": None,
            "pressure_sml": None,
            "pressure_plate": None,
            "pressure_screen": None,
            "ch_average": None,
            "ch_deep": None,
            "ch_sml": None,
            "ch_plate": None,
            "ch_screen": None,
            "coverage_average": None,
            "coverage_deep": None,
            "coverage_sml": None,
            "coverage_plate": None,
            "coverage_screen": None
        }

    def __rpr__(self):
        return f'Station {self.station_number} on date {self.date}'

    def __str__(self):
        return self.__rpr__()

    def __lt__(self, other):

        if self.date < other.date:
            return True

        elif self.date == other.date:

            if int(self.station_number) < int(other.station_number):
                return True

    def print_stats(self):
        """Formatted output of the stations stats, calculated from the LtIsotherms belonging to the
        station."""
        s = ''
        for item in self.stats:
            s += f'{item} : {self.stats[item]}\n'
        return s

    def print_data(self):
        s = ""
        for sp in self.sfg_spectra:
            s += f'{sp.name.full_name}\n'

        for l in self.lt_isotherms:
            s += f'{str(l)}\n'

        for t in self.tensions:
            s += str(t) + "\n"

        return s

    def get_doy(self):
        doy = self.date.timetuple().tm_yday
        factor = (1 - int(self.station_number)) * 0.25
        return doy + factor

    def get_value_by_type(self, type_, value):

        dic = {"a": ("s", "p", "c"), "deep": ("c",), "sml": ("p", "s"), "s": ("s",), "p": ("p",)}
        types = dic[type_]
        out = []

        if value == "tension":

            if len(self.tensions) > 0:

                for tension in self.tensions:
                    h = SampleHash(tension[0])
                    if h.sample_type in types:
                        out.append(tension[1])

        elif value == "ch":

            for spec in self.sfg_spectra:

                if spec.get_sample_hash().sample_type in types:
                    if np.min(spec.wavenumbers) <= 2750:
                        out.append(spec.calculate_ch_integral(average="gernot"))

        elif value == "dppc":

            for spec in self.sfg_spectra:

                if spec.get_sample_hash().sample_type in types:

                    if np.min(spec.wavenumbers) <= 2750:
                        integral = spec.calculate_ch_integral(average="gernot")
                        dppc_integral = self.parent.dppc_ints[spec.name.creation_time.date()]
                        temp = integral / dppc_integral

                        if temp < 0:
                            temp = 0

                        out.append(np.sqrt(temp))


        elif value == "max_pres":

            if len(self.lt_isotherms) > 0:

                isotherms = {}

                for isotherm in self.lt_isotherms:

                    if isotherm.sample_hash.namestring not in isotherms:
                        isotherms[isotherm.sample_hash.namestring] = isotherm
                    else:
                        if isotherm.measured_time < isotherms[isotherm.sample_hash.namestring].measured_time:
                            isotherms[isotherm.sample_hash.namestring] = isotherm

                for isotherm in isotherms.values():

                    if isotherm.sample_hash.sample_type in types:
                        pres = isotherm.get_maximum_pressure()
                        if pres < 72:
                            out.append(pres)

        if len(out) > 0:
            average = np.average(out)
            std = np.std(out)
            count = len(out)

            return average, std, count

        else:
            return None

    def analyze_station_data(self):

        # tension
        if len(self.tensions) > 0:
            self.stats["tension_average"] = self.get_value_by_type("a", "tension")
            self.stats["tension_deep"] = self.get_value_by_type("deep", "tension")
            self.stats["tension_sml"] = self.get_value_by_type("sml", "tension")
            self.stats["tension_plate"] = self.get_value_by_type("p", "tension")
            self.stats["tension_screen"] = self.get_value_by_type("s", "tension")

        # surface pressure
        if len(self.lt_isotherms) > 0:
            self.stats["pressure_average"] = self.get_value_by_type("a", "max_pres")
            self.stats["pressure_deep"] = self.get_value_by_type("deep", "max_pres")
            self.stats["pressure_sml"] = self.get_value_by_type("sml", "max_pres")
            self.stats["pressure_plate"] = self.get_value_by_type("p", "max_pres")
            self.stats["pressure_screen"] = self.get_value_by_type("s", "max_pres")

        # SFG data
        if len(self.sfg_spectra) > 0:
            self.stats["ch_average"] = self.get_value_by_type("a", "ch")
            self.stats["ch_deep"] = self.get_value_by_type("deep", "ch")
            self.stats["ch_sml"] = self.get_value_by_type("sml", "ch")
            self.stats["ch_plate"] = self.get_value_by_type("p", "ch")
            self.stats["ch_screen"] = self.get_value_by_type("s", "ch")

            self.stats["coverage_average"] = self.get_value_by_type("a", "dppc")
            self.stats["coverage_deep"] = self.get_value_by_type("deep", "dppc")
            self.stats["coverage_sml"] = self.get_value_by_type("sml", "dppc")
            self.stats["coverage_plate"] = self.get_value_by_type("p", "dppc")
            self.stats["coverage_screen"] = self.get_value_by_type("s", "dppc")

    def arange_to_sample(self):
        """Matches the analytical data of the station to the corresponding sampless
        in order to keep them together for further analysis"""
        samples = []

        for s in self.sfg_spectra:
            if s.get_sample_hash() not in samples:
                samples.append(s.get_sample_hash())

        for l in self.lt_isotherms:  # type: LtIsotherm
            if l.sample_hash not in samples:
                samples.append(l.sample_hash)

        for t in self.tensions:
            if SampleHash(t[0]) not in samples:
                samples.append(SampleHash(t[0]))

        for sample in samples:
            self.samples.append(Sample(sample))

        for sample in self.samples:

            for spec in self.sfg_spectra:
                if spec.get_sample_hash() == sample.sample_hash:
                    sample.sfg_spectra.append(spec)

            for l in self.lt_isotherms:
                if l.sample_hash == sample.sample_hash:
                    sample.lt_isotherms.append(l)

            for t in self.tensions:
                if SampleHash(t[0]) == sample.sample_hash:
                    sample.surface_tension = t


class LtManager:
    """A class to perform operations on a set of LtIsotherm objects. It is an extension to the SessionControllManager
    to extend his features with isotherm handling. It relies on sqlite databases as well."""

    # todo: obsolete. The small remaining amount of code does not justify a class
    def __init__(self, database, table="lt_gasex"):

        self.database = database
        self.cursor = database.cursor()
        self.table = table
        self.isotherms = []
        self.days = None
        self.ordered_days = {}

        self.get_all_isotherms()
        #todo : the functions below have to be replaced by the sample/station hash system
        #self.join_days_isotherms()
        #self.order_by_sample()
        #self.join_same_measurement()

    def get_all_isotherms(self):
        """Fetches all isotherm data from the database and transforms them to LtIsotherm objects."""
        command = "SELECT * from " + self.table
        self.cursor.execute(command)
        result = self.cursor.fetchall()
        for i in result:
            lt = GasExLtIsotherm(i[1], i[2], i[4], i[5], i[6], i[7])
            self.isotherms.append(lt)


class SfgSpectrum:
    """The SFG spectrum class is the foundation of all analysis and plotting tools. It contains a class
    SystematicName (or a derived class) which carries most of the metainformation. Besides holding the
    experimental data, it gives access to a variety of functions like normalization, peak picking etc."""

    # magic methods
    def __init__(self, wavenumbers, intensity, ir_intensity, vis_intensity, SystematicName):
        self.wavenumbers = wavenumbers
        self.raw_intensity = intensity
        self.vis_intensity = vis_intensity
        self.ir_intensity = ir_intensity
        self.name = SystematicName
        self.normalized_intensity = self.raw_intensity / (self.vis_intensity * self.ir_intensity)
        self.peaks = self.detailed_analysis()
        self.baseline_corrected = None

    def __repr__(self):
        return self.name.full_name[:-4] + "   " + str(self.yield_spectral_range())

    def __str__(self):
        """Printable representation of the SFG object"""
        date = self.name.date
        surf = self.name.surfactant
        srange = str(self.yield_spectral_range())
        full = self.name.full_name[:-4]

        return date + "\t" + surf + "\t" + srange + "\t" + full

    def __add__(self, SFG2):
        """Definition of an addition method for SFG spectra. Returning an Added_SFG object"""

        if self.wavenumbers.all() != SFG2.wavenumbers.all():
            print("Warning! SFG data do not have the same shape!")

        wavenumbers = self.wavenumbers
        intensity = self.normalized_intensity

        values = [(a, b) for a, b in zip(wavenumbers, intensity)]
        for i in range(len(SFG2.wavenumbers)):
            values.append((SFG2.wavenumbers[i], SFG2.normalized_intensity[i]))

        values.sort()
        values = values[::-1]
        new_wavenumbers = []
        new_intensities = []
        last_wl = 0

        for i in range(len(values)):
            tup = values[i]
            if last_wl != tup[0]:
                new_wavenumbers.append(tup[0])
                new_intensities.append(tup[1])
                last_wl = tup[0]
            else:
                x = len(new_intensities) - 1
                new_intensities[x] = (new_intensities[x] + tup[1]) * 0.5

        new_wavenumbers = np.array(new_wavenumbers)
        new_intensities = np.array(new_intensities)

        if isinstance(SFG2, SfgSpectrum):

            names = [self.name.full_name[:-4], SFG2.name.full_name[:-4]]
            name = AddedName(names)
            Added = AddedSpectrum((new_wavenumbers, new_intensities), name)
            return Added

        elif isinstance(SFG2, AddedSpectrum):
            names = self.name.full_name + "_" + SFG2.name.full_name

            sensitizers = self.name.sensitizer
            for i in SFG2.name.sensitizer:
                if i not in sensitizers:
                    sensitizers.append(i)

            surfactants = self.name.surfactants
            for i in SFG2.name.sensitizer:
                if i not in surfactants:
                    surfactants.append(i)

            name = AddedName(names, sensitizers, surfactants)
            Added = AddedSpectrum((new_wavenumbers, new_intensities), name)
            Added.speccounter += self.speccounter
            Added.speccounter += SFG2.speccounter
            return Added

    def __lt__(self, SFG2):
        """Returns true if the current spectrum was measured before SFG2"""
        if self.name.creation_time < SFG2.name.creation_time:
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

    def smooth(self, points=9, order=5):
        """Apply a smooth routine to the normalized_intensity"""
        y = savgol_filter(self.normalized_intensity, points, order)
        return y

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

    def detailed_analysis(self, threshold=2, intensity="norm"):
        """Function returning peak information (central wavenumber, flanks, integral). The threshold value characterizes
        the factor that a peak must be greater than the average intensity"""
        x_array = self.wavenumbers[::-1]
        if intensity == "norm":
            y_array = self.normalized_intensity[::-1]
        elif intensity == "raw":
            y_array = self.raw_intensity[::-1]

        slopes = [(y_array[i + 1] - y_array[i]) / (x_array[i + 1] - x_array[i]) for i in range((len(x_array) - 1))]

        possible_peaks = []
        peak_tuples = []

        for i in range(1, len(y_array) - 1):
            if slopes[i - 1] > 0 and slopes[i] < 0:
                possible_peaks.append(i)

        average_intensity = np.average(y_array)

        confirmed_peaks = [i for i in possible_peaks if (y_array[i] > average_intensity * threshold)]

        for i in confirmed_peaks:

            left = 0
            right = 0
            center = i
            k = i - 2
            border = np.average(slopes) * 1 / 6

            # check for left border
            while slopes[k] > border and k >= 0:
                k -= 1
            left = k

            # check for right border
            k = i + 1

            if k >= len(slopes):
                k = len(slopes) - 1
            else:
                # if traversing to the right does not find a proper peak ending
                try:
                    while (slopes[k] < border):
                        k += 1
                except IndexError:
                    k = k - 1
                    break
            right = k

            peak_tuples.append((center, left, right))

        data_out = []
        for i in peak_tuples:
            indices = (i[0], i[1], i[2])
            center = x_array[i[0]]
            center_intensity = y_array[i[0]]
            left = x_array[i[1]]
            right = x_array[i[2]]
            peak_slice_x = x_array[i[1]:i[2] + 1]
            peak_slice_y = y_array[i[1]:i[2] + 1]
            area = self.integrate_peak(peak_slice_x, peak_slice_y)
            datapoints = len(peak_slice_x)

            data_out.append(
                (center, left, right, center_intensity, peak_slice_x, peak_slice_y, datapoints, area, indices))

        # sort peaks by peak intensity
        data_out = sorted(data_out, key=(lambda x: x[3]), reverse=True)
        return data_out

    # info functions

    def drop_ascii(self):
        """Create an ascii file with the wavenumbers and normalized intensities"""
        with open(self.name.full_name[:-4] + ".csv", "w") as outfile:
            writer = csv.writer(outfile, delimiter=";")
            for i in zip(self.wavenumbers, self.normalized_intensity):
                writer.writerow((i[0], i[1]))

    def drop_tex_peaktable(self, threshold=1.5):

        peaks = self.detailed_analysis(threshold=threshold)
        tablestring = ""

        counter = 1
        for peak in peaks:
            tablestring += "\subcaption*{" + "Peak " + str(counter) + "}\n"
            tablestring += "\\begin{tabular}{|c|c|c|c|c|c|}\n\hline\n"
            tablestring += "Wavenumber & normalized intensity & area & left border & right border" + "\\\\" + "\hline\n"

            tablestring += str(peak[0]) + " & " + str(peak[3]) + " & " + str(peak[7]) + " & " + str(
                peak[1]) + " & " + str(peak[2]) + "\\\\"
            tablestring += "\n\hline\n\end{tabular}\n"
            counter += 1

        return tablestring

    def get_sample_hash(self):
        if isinstance(self.name, SystematicGasExName) is True:
            return SampleHash(self.name.sample_hash_string)

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
    def calc_dish_area(self, diameter=5.1):
        """A auxialiary function to calculate the area of a teflon dish in square angstroms. Diameter given in cm."""
        radius = diameter * 0.5
        area = np.pi * radius ** 2
        area = area * 10 ** 16  # conversion to square angstroms
        return area

    def calc_area_per_molecule(self):
        """The function calculates the area per molecule. The area should be given in square angstroms, the
        concentration in milimole per liter and the volume in microliter"""
        if self.name.sensitizer == "-" or self.name.sensitizer in ["Benzophenone", "Humic Acid"]:
            concentration = self.name.surf_stock_concentration
            if self.name.surf_stock_concentration == "unknown":
                concentration = input("Enter surf stock concentration of spectrum " + self.name.full_name + ": \n")

            volume = float(self.name.surfactant_spread_volume)

            concentration = float(concentration) * 10 ** -3  # conversion in mol per liter
            volume = volume * 10 ** -6  # conversion in liter
            amount = volume * concentration

        else:
            concentration = self.name.surf_stock_concentration
            if self.name.surf_stock_concentration == "unknown":
                concentration = input(
                    "Enter surf stock concentration of spectrum " + self.name.full_name + ": \n")

            volume = float(self.name.surfactant_spread_volume)
            concentration = float(concentration) * 10 ** -3  # conversion in mol per liter
            volume = volume * 10 ** -6  # conversion in liter
            amount_su = volume * concentration

            sens_stock_conc = input("Enter sens stock concentration of spectrum " + self.name.full_name + ": \n")
            concentration = float(sens_stock_conc) * 10 ** -3  # conversion in mol per liter
            volume = float(self.name.sensitizer_spread_volume) * 10 ** -6  # conversion in liter
            amount_se = volume * concentration

            amount = amount_se + amount_su

        molecules = (6.022 * 10 ** 23) * amount  # number of molecules
        area_per_molecule = self.calc_dish_area() / molecules

        return area_per_molecule

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

