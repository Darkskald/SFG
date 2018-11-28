# module-internal imports
from new_gui import run_app, run_lt_app

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

# for session controll manager
# from new_gui import run_app

rcParams['mathtext.default'] = 'regular'


# noinspection PySimplifyBooleanCheck,PySimplifyBooleanCheck 1234


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
        return (self.name.creation_time < SFG2.name.creation_time)

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

            right = self.slice_by_borders(3000, 2950)

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


class AddedSpectrum(SfgSpectrum):
    """The result of the addition of two SFG spectra constructed by averaging over datapoints which are shared
    by the spectra. Note: If only one of the two spectra has the datapoint, it will be in the result spectrum with
    the same intensity value like in the original"""

    def __init__(self, wn_intenstup, name):
        self.wavenumbers = wn_intenstup[0]
        self.normalized_intensity = wn_intenstup[1]
        self.speccounter = 2
        self.name = name

    def __str__(self):
        s = "sum of " + str(self.speccounter) + " SFG spectra with systematic names" + str(self.name.full_name)
        return s

    def __add__(self, SFG2):

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
            Added.speccounter = self.speccounter + 1
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


class SystematicName:
    """A class extracting metainformation from the string filename and storing it in rather self-explanatory
    variables.(eg date of measurement, type of surfactant...)"""

    def __init__(self, namestring, creation_time="unknown"):

        self.refpath = "name_info/"

        try:
            with open(self.refpath + "Surfactants.txt") as outfile:
                pass
        except FileNotFoundError:
            self.refpath = "../name_info/"

        self.creation_time = creation_time
        # load the allowed surfactans and sensitizers from files
        self.Surfactants = {}
        self.Sensitizers = {}

        with open(self.refpath + "Surfactants.txt", "r") as infile:
            for line in infile:
                collect = line.split(":")
                self.Surfactants[collect[0]] = collect[1].strip()

        with open(self.refpath + "Sensitizers.txt", "r") as infile:
            for line in infile:
                collect = line.split(":")
                self.Sensitizers[collect[0]] = collect[1].strip()

        # now the actual processing begins
        self.full_name = namestring
        processing_list = self.full_name[:-4]  # removes .sfg
        self.processing_list = processing_list.split("_")
        self.date = self.processing_list[0]

        # setting defaults
        self.surfactant = "unknown"
        self.surfactant_spread_volume = "unknown"
        self.sensitizer = "-"
        self.sensitizer_spread_volume = 0
        self.photolysis = "none"
        self.sample_number = 1
        self.measurement = "#1"  # das ist nicht optimal, sample_number als int und measurement als string zu haben
        self.comment = "none"
        self.surf_stock_concentration = "unknown"
        if len(self.processing_list) == 2:
            self.check_boknis()
        # traversing the processing list
        for i in self.processing_list[1:]:

            if i in self.Surfactants:
                self.surfactant = self.Surfactants[i]

            elif i in self.Sensitizers:
                self.sensitizer = self.Sensitizers[i]

            # the following line could result in a problem with the date (processing_list[0])
            elif self.is_number(i) == True:
                if self.surfactant_spread_volume == "unknown":
                    self.surfactant_spread_volume = i
                else:
                    self.sensitizer_spread_volume = i

            # pH value info will be comment#todo fix that this is bullshit
            elif "p" in i and "H" in i:
                self.comment = i

            # photolysis handling
            elif "p" in i:
                if self.photolysis == "none":
                    time = i.strip("p")
                    if self.is_number(time) == True:
                        self.photolysis = time + " min."
                    elif "h" in time:
                        if self.is_number(time.strip("h")) == True:
                            self.photolysis = str(int(time.strip("h")) * 60) + " min."
                        else:
                            self.comment = i
                else:
                    self.comment = i

            # check for sample number
            elif "x" in i:
                if self.is_number(i.strip("x")) == True:
                    self.sample_number = i.strip("x")
                else:
                    self.comment = i

            # check for measurement number
            elif "#" in i:
                self.measurement = i  # hier auch ein String! ggf. aendern

            elif self.is_number(i.strip("mM")) == True:
                self.surf_stock_concentration = float(i.strip("mM"))

            elif "mM" in i:
                index = i.find("mM")
                if self.is_number(i[index - 1]):
                    self.surf_stock_concentration = float(i[index - 1])

            else:
                if self.comment == "none":
                    self.comment = i
                else:
                    self.comment += i

        # postprocessing

        if self.surfactant == "unknown" and self.sensitizer != "-":
            self.surfactant = self.sensitizer
            self.sensitizer = "-"

        if self.sensitizer == "-":
            self.sensitizer_spread_volume = 0

        if self.sensitizer == "-" and self.surfactant == "DPPC":
            self.surf_stock_concentration = 1

    def is_number(self, s):
        """Auxiliary function returning e boolean, depending on test variable can be converted in a float"""
        try:
            float(s)
            return True
        except ValueError:
            return False

    def yield_info(self):
        return (self.full_name, self.date, self.surfactant, self.surfactant_spread_volume,
                self.sensitizer, self.sensitizer_spread_volume, self.photolysis,
                self.sample_number, self.measurement, self.comment)

    def check_boknis(self):
        """Function determining if a spectrum belongs to the BE time series. Those follow a certain name convention"""

        if self.is_number(self.processing_list[0]) == True and self.is_number(self.processing_list[1]) == True:
            self.comment = "BoknisEckSample"
            self.sample_number = self.processing_list[1]
            self.surfactant = "Nat. surface sample"

    def date_split(self):
        """Splitting the date which is stored as an integer (YYYYMMDD), extracting year, month and day information
        separately"""

        year = int(self.date[0:4])
        month = int(self.date[4:6])
        day = int(self.date[6:])

        return (year, month, day)


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


class AddedName(SystematicName):
    """This class is derived from the SfgSpectrum and represents the result of the addition of Sfg intensities."""

    def __init__(self, names):
        self.full_name = ("_").join(names)


class SessionControlManager:
    """Successor of the IpyInterpreter. A class to access all experimental data, search data by certain match criteria
    and produce plots efficiently. Designed to work with a sqllite database containing the spectral raw data. Especially
    useful from interactive python environments (eg IPy)."""

    def __init__(self, database, id):
        # todo: handling of UV/Ir/Raman-Data has to be included
        self.db = sqlite3.connect(database, detect_types=sqlite3.PARSE_DECLTYPES | sqlite3.PARSE_COLNAMES)
        self.cur = self.db.cursor()

        self.session_id = id
        self.Surfactants = {}
        self.Sensitizers = {}
        self.Makros = {}
        self.stations = None
        self.get_senssurf_names()

        self.lt_manager = None
        self.spec_manager = None
        self.tensions = None

        # former IpyInterpreter functionality, tracking the primary key in parallel
        self.subset_ids = []
        self.subset = []

        self.recover_ids = []
        self.recover = []

        self.dppc_ints = None

        self.set_spec_manager()


    # database methods

    def fetch_single(self, number, database, condition3=None, condition4=None):
        """Fetches the data of a single SFG spectrum from the database. Returns an SfgSpectrum object. The
        conditional parameters allow a selection of the spectrum according to specified properties. The
        default_data kwarg controlls which database the spectrum should be extracted from."""

        # todo: A separate table for BoknisEck data is desirable

        command = "SELECT * FROM " + database + " WHERE ROWID=" + str(number)

        if condition3 is not None:
            command += " AND " + str(condition3) + "=" + str(condition4)

        self.cur.execute(command)

        result = self.cur.fetchall()[0]
        spec = self.construct_sfg(result, database)

        return spec

    def construct_sfg(self, query_result, database):
        """A function to create an SFG spectrum from the result of a query in the SQL database"""

        result = query_result
        creationtime = result[2]

        if database == "sfg_database":

            sysname = SystematicName(result[1], creationtime)
            wavenumber = np.array(result[-7].split(";")).astype(np.float)
            sfg = np.array(result[-6].split(";")).astype(np.float)
            ir = np.array(result[-5].split(";")).astype(np.float)
            vis = np.array(result[-4].split(";")).astype(np.float)

        else:

            if "DPPC" not in result[1]:
                sysname = SystematicGasExName(result[1], creationtime)
                wavenumber = np.array(result[-4].split(";")).astype(np.float)
                sfg = np.array(result[-3].split(";")).astype(np.float)
                ir = np.array(result[-2].split(";")).astype(np.float)
                vis = np.array(result[-1].split(";")).astype(np.float)

            else:
                sysname = SystematicName(result[1], creationtime)
                wavenumber = np.array(result[-4].split(";")).astype(np.float)
                sfg = np.array(result[-3].split(";")).astype(np.float)
                ir = np.array(result[-2].split(";")).astype(np.float)
                vis = np.array(result[-1].split(";")).astype(np.float)

        spec = SfgSpectrum(wavenumber, sfg, ir, vis, sysname)

        return spec

    def general_fetch(self, condition_1=None, condition_2=None, database="sfg_database"):
        """A function to fetch spectra from the database according to specified match criteria.
        The default_data kwarg controlls the datatable that is used."""

        command = "SELECT * from " + database

        if condition_1 is not None and condition_2 is not None:
            command += (" WHERE " + condition_1 + "=" + condition_2)

        self.cur.execute(command)

        for item in self.cur.fetchall():
            id = item[0]
            self.subset_ids.append(id)

            try:
                self.subset.append(self.construct_sfg(item, database))
            except:
                logging.error(traceback.format_exc())

    def general_refine(self, condition1, condition2, database):
        """Refinement of the actual subset by applying further match criteria. This is the
        actual implementation of the get method abstracting away the database access routine from the
        user."""
        temp = []
        temp_id = []

        for id in self.subset_ids:

            try:
                s = self.fetch_single(id, database, condition1, condition2)
                temp.append(s)
                temp_id.append(id)
            except IndexError:
                pass

        self.recover = self.subset
        self.recover_ids = self.subset_ids

        self.subset = temp
        self.subset_ids = temp_id

    # user interaction functions

    def plot(self):
        """Runs an external PyQt application plotting all current spectra in the subset. The external
        app gives access to further manual analysis functionality"""
        run_app(self.subset, self.session_id)

    def show(self):
        """Prints a list of the spectra of the current subset including their subset index by which they
        can be accessed easily"""
        for i, spectrum in enumerate(self.subset):
            print(str(i) + " " + str(spectrum))

    def clear(self):
        """Deletes all spectra from the current subset"""
        self.subset = []
        self.subset_ids = []

    def get(self, flagstring, database="sfg_database", ref=False):
        """Fetches spectra according to the desired properties and adds them to the subset"""
        t = self.flagstring_split(flagstring)

        if t[0] == "su" or t[0] == "se":
            condition1 = self.retranslate_name(t[0])
            condition2 = "\"" + self.retranslate_name(t[1]) + "\""

            if ref is False:
                try:
                    self.general_fetch(condition1, condition2, database=database)

                except:
                    pass

            if ref is True:
                try:
                    self.general_refine(condition1, condition2, database)

                except:
                    pass

        elif t[0] == "name":
            self.general_fetch(t[0], "\"" + t[1] + ".sfg" + "\"")

    def ref(self, flagstring):
        """Refines the current subset by applying further match criteria"""
        self.get(flagstring, database="sfg_database", ref=True)

    def by_time(self, time1, time2, database, refine=False):
        """Fetch or  refine the spectral data by time of measurement. The Number has to be given as a string, embraced
        in '' quotation marks to pass it to the SQL query. Note this function is the only user function that directly
        accesses the SQL database"""

        if refine == False:

            command = "SELECT * from " + database + " WHERE measured_time between " + time1 + " and " + time2

            self.cur.execute(command)
            keys = []
            for item in self.cur.fetchall():
                id = item[0]
                self.subset_ids.append(id)
                self.subset.append(self.fetch_single(id))
        else:
            temp = []
            temp_id = []

            for id in self.subset_ids:

                try:
                    command = "SELECT * FROM " + database + " WHERE ROWID=" + str(
                        id) + " AND measured_time between " + time1 + " and " + time2
                    self.cur.execute(command)
                    result = self.cur.fetchall()[0]
                    s = self.construct_sfg(result, database)
                    temp.append(s)
                    temp_id.append(id)
                except IndexError:
                    pass

            self.recover = self.subset
            self.recover_ids = self.subset_ids

            self.subset = temp
            self.subset_ids = temp_id

    def rec(self):
        """A recover function that resets the subset to the state before the last refinement or keep call"""
        self.subset = self.recover
        self.subset_ids = self.recover_ids

    def remove(self, numbers):
        """Removes items by (a list of) index(indices from the subset"""

        options = numbers.split(",")
        to_remove = [self.subset[int(i)] for i in options]
        newlist = [i for i in self.subset if i not in to_remove]
        new_indexlist = [self.subset_ids[int(i)] for i in options]
        self.recover = self.subset
        self.recover_ids = self.subset_ids

        self.subset = newlist
        self.subset_ids = new_indexlist

    def keep(self, flagstring):
        """Removes everything but the specified items. Specification is done by a list of indices"""
        f = flagstring
        options = f.split(",")
        new_list = [self.subset[int(i)] for i in options]
        new_indexlist = [self.subset_ids[int(i)] for i in options]
        recover = [i for i in self.subset if i not in new_list]

        self.recover = recover
        self.recover_ids = self.subset_ids

        self.subset = new_list
        self.subset_ids = new_indexlist

    # specific functionality for GasEx

    def set_lt_manager(self):
        """Creates a LtManager object to give access to analysis tools for Langmuir trough isotherms besides
        the built-in SFG functionality. It can be seen as an add-in to the SessionControlManager to extend
        his capabilities."""

        self.lt_manager = LtManager(self.db)

    def match_same_measurements(self):
        """Returns a dictionary containing the sample hash as keys and a list of all the SFG measurements that were
        performed with this sample"""

        out = {}
        for spectrum in self.subset:

            if spectrum.get_sample_hash() not in out:
                out[spectrum.get_sample_hash()] = []
                out[spectrum.get_sample_hash()].append(spectrum)

            else:
                out[spectrum.get_sample_hash()].append(spectrum)

        return out

    def collect_stations(self):
        """Creates Station objects from the station information of the currently available isotherms and
        spectra. This helps to keep track of the set of samples taken within one GasEx sampling station"""
        station_hashes = []
        types = []

        for spec in self.subset: #type: SfgSpectrum

            if isinstance(spec.name, SystematicGasExName) is True:
                h = spec.get_sample_hash()
                _hash = h.station_hash


                if _hash not in station_hashes:
                    station_hashes.append(_hash)

                    if h.station_type == "a":
                        types.append("small")

                    elif h.station_type in ("r", "c"):
                        types.append("big")

                    else:
                        print("Unknwonwn station type!")

        for isotherm in self.lt_manager.isotherms:#type: LtIsotherm
            _hash = isotherm.sample_hash.station_hash
            if _hash not in station_hashes:
                station_hashes.append(_hash)
                if isotherm.sample_hash.station_type == "a":
                    types.append("small")

                elif isotherm.sample_hash.station_type in ("r", "c"):
                    types.append("big")

                else:
                    print("Unknwonwn station type!")

        for tension in self.tensions:
            try:
                h = SampleHash(tension)
                _hash = h.station_hash

                if _hash not in station_hashes:
                    station_hashes.append(_hash)
                    if h.station_type == "a":
                        types.append("small")

                    elif h.station_type in ("r", "c"):
                        types.append("big")

                    else:
                        print("Unknwonwn station type!")

            except IndexError:
                print(tension)

        self.stations = {}
        for s in zip(station_hashes, types):
            self.stations[s[0]] = Station(s[0], s[1], parent=self)

    def map_to_stations(self):

        for spectrum in self.subset:
            if isinstance(spectrum.name, SystematicGasExName) is True:
                self.stations[spectrum.get_sample_hash().station_hash].sfg_spectra.append(spectrum)

        for tension in self.tensions:
            try:
                self.stations[SampleHash(tension).station_hash].tensions.append([tension, self.tensions[tension]])
            except IndexError:
                print(tension)


        for lt_isotherm in self.lt_manager.isotherms:
            self.stations[lt_isotherm.sample_hash.station_hash].lt_isotherms.append(lt_isotherm)

    def fetch_gasex_sfg(self):
        """Fetches all SfgSpectra from the GasEx cruise and puts them in the subset attribute"""
        self.general_fetch(database="sfg_gasex")

    def setup_for_gasex(self):
        """Convenience function, calls all necessary functions to setup the SCM for gasex data
        processing"""
        self.set_lt_manager()
        self.fetch_gasex_sfg()
        self.tensions = self.fetch_tension_data()
        self.collect_stations()
        self.map_to_stations()
        self.dppc_ints = self.get_dppc_average()

        for station in self.stations.values():
                station.analyze_station_data()


    def fetch_tension_data(self):
        out = {}
        command = f'SELECT * from gasex_surftens'
        self.cur.execute(command)
        for item in self.cur.fetchall():
            out[item[1]] = float(item[2])
        return out


    # handling DPPC normalization
    def get_dppc_average(self):
        """Takes a list of SfgSpectra as arguments. Returns a dictionary of dates as keys and the average SFG CH
        integral as values."""

        temp = {}
        out = {}
        #todo: erase the error which arises if the spectrum is not measured till 3000 cm
        for spec in self.subset:  # type: SfgSpectrum

            if not (isinstance(spec.name, SystematicGasExName)):

                if spec.name.surfactant == "DPPC":

                    if np.max(spec.wavenumbers >= 3000):

                        try:
                            temp[spec.name.creation_time.date()].append(spec)
                        except KeyError:
                            temp[spec.name.creation_time.date()] = [spec]



        for key in temp:

            intens = 0
            counter = 0

            for spec in temp[key]:
                intens += spec.calculate_ch_integral()
                counter += 1

            intens /= counter
            out[key] = intens

        return out

    # auxiliary functions

    def flagstring_split(self, flagstring):
        """This function processes the given flagstring and returns a list of SFG objects
        which are passed through the Finder methods utilizing the flags and options"""
        f = flagstring.split(" ")
        return f

    def get_senssurf_names(self):
        """Loads allowed names of surfactants and sensitizers from the control file"""

        with open("name_info/Surfactants.txt", "r") as infile:
            for line in infile:
                collect = line.split(":")
                self.Surfactants[collect[0]] = collect[1].strip()

        with open("name_info/Sensitizers.txt", "r") as infile:
            for line in infile:
                collect = line.split(":")
                self.Sensitizers[collect[0]] = collect[1].strip()

        with open("name_info/makros.txt", "r") as infile:
            for line in infile:
                collect = line.split(":")
                self.Makros[collect[0]] = collect[1].strip()

    def retranslate_name(self, stri):
        """Extracts name information from the string, returning the full name of the substance abbreviation"""

        if stri in self.Surfactants:
            return self.Surfactants[stri]
        elif stri in self.Sensitizers:
            return self.Sensitizers[stri]
        elif stri in self.Makros:
            return self.Makros[stri]
        else:
            print("Retranslation failed. Unknown expression.")

    def set_spec_manager(self):
        self.spec_manager = SpectraManager(self.db)



class LtIsotherm:
    """A class to represent experimental Langmuir trough isotherms, handling time, area, area per molecule and
    surface pressure"""

    def __init__(self, *args):

        self.name = args[0]
        self.measured_time = args[1]
        self.time = np.array(args[2].split(";")).astype(np.float)
        self.area = np.array(args[3].split(";")).astype(np.float)
        self.apm = np.array(args[4].split(";")).astype(np.float)
        self.pressure = np.array(args[5].split(";")).astype(np.float)
        self.compression_factor = self.calc_compression_factor()


        self.sample_hash = None
        self.speed = None
        self.measurement_number = None

        self.partners = [] #todo obsolete

        self.process_name()

    def __str__(self):
        return self.name + " LtIsotherm Object"

    def __repr__(self):
        return self.name + " LtIsotherm Object"

    def __lt__(self, other):
        if self.get_maximum_pressure() < other.get_maximum_pressure():
            return True

    def process_name(self): #todo odjust to create hash object
        """Function to extract metainformation from the filename"""
        temp = self.name.split("_")

        hashstring = temp[1] + "_" + temp[2] + "_" + temp[3]

        if "c" not in temp[3][0]:
            hashstring += ("_"+temp[4])


        self.sample_hash = SampleHash(hashstring)

        if len(temp) >= 5:
            self.measurement_number = temp[-1]
            self.speed = temp[-2]

        else:
            self.measurement_number = 1

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

    def same_sample(self, other):
        """Checks wether two isotherms belong to the same sample. This is the case if the same sample
        is measured several times in a row. Returns a bool."""
        if self.sample_hash == other.sample_hash:
            return True
        else:
            return False

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

    def smooth(self):
        """Performs a smooth operation of the measured pressure involving a Savitzky-Golay-filter"""
        return savgol_filter(self.pressure, 9, 5)

    def cut_away_decay(self, x_array):

        max = self.get_maximum_pressure()
        index = int(np.where(self.pressure==max)[0][0])
        return self.get_slice(x_array, 0, index)

    @staticmethod
    def process_lift_off(x, y):

        start = 0
        diff = 0
        increment = int(0.15 * len(x))
        line1 = None
        line2 = None
        x_out = None

        while start < (len(x) - increment * 2):

            left_x = x[start:start + increment]
            left_y = y[start:start + increment]

            right_x = x[start + increment:start + (2 * increment)]
            right_y = y[start + increment:start + (2 * increment)]

            slope_l, intercept_l, r_value, p_value, std_err = \
                stats.linregress(left_x, left_y)

            slope_r, intercept_r, r_value, p_value, std_err = \
                stats.linregress(right_x, right_y)

            temp = np.abs(slope_l - slope_r)

            if temp > diff:
                diff = temp
                line1 = [slope_l, intercept_l]
                line2 = [slope_r, intercept_r]
                x_out = x[start:start + (2 * increment)]

            start += 1

        return diff, line1, line2, x_out

    def find_lift_off(self, x_array):
        x, y = self.cut_away_decay(x_array)
        tup = LtIsotherm.process_lift_off(x, y)
        lift_off = LtIsotherm.calculate_intersection(tup[1], tup[2])
        return lift_off

    @staticmethod
    def calculate_intersection(line1, line2):

        x = (line2[1] - line1[1]) / (line1[0] - line2[0])

        y = line1[0] * x + line1[1]

        return x, y





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
            lt = LtIsotherm(i[1], i[2], i[4], i[5], i[6], i[7])
            self.isotherms.append(lt)



class Station:
    """A class carrying all information and data for a given cruise station, especially SFG and isotherms"""

    def __init__(self, name, type=None, parent=None):

        self.station_hash = name
        self.date = datetime.date(2018,int(name[0:2]),int(name[2:4]))
        self.station_number = name[4]
        self.type = type

        self.sfg_spectra = []
        self.lt_isotherms = []
        self.tensions = []

        self.lt_joined = {}
        
        self.parent = parent
        self.stats = {
                
                "tension_average": None,
                "tension_deep":None,
                "tension_sml":None,
                "tension_plate":None,
                "tension_screen":None,
                "pressure_average":None,
                "pressure_deep":None,
                "pressure_sml":None,
                "pressure_plate":None,
                "pressure_screen":None,
                "ch_average":None,
                "ch_deep":None,
                "ch_sml":None,
                "ch_plate":None,
                "ch_screen":None,
                "coverage_average":None,
                "coverage_deep":None,
                "coverage_sml":None,
                "coverage_plate":None,
                "coverage_screen":None
                }
        

    def __rpr__(self):
        return f'Station {self.station_number} on date {self.date}'

    def __str__(self):
        return self.__rpr__()
    
    def __lt__(self, other):
        
        if self.date < other.date:
            return True
        
        elif self.date == other.date:
            
            if int(self.station_number) < int( other.station_number):
                return True


    def join_samples(self):
        """Joins Langmuir trough measurements of the same sample. Much more comprehensive than
        the approach in the LtManager."""

        for isotherm in self.lt_isotherms:

            if isotherm.create_sample_hash() not in self.lt_joined:
                self.lt_joined[isotherm.create_sample_hash()] = [isotherm]
            else:
                self.lt_joined[isotherm.create_sample_hash()].append(isotherm)

        for isolist in self.lt_joined.values():
            isolist.sort(reverse=True)

        self.isotherm_count = len(self.lt_joined)

    def print_stats(self):
        """Formatted output of the stations stats, calculated from the LtIsotherms belonging to the
        station."""

        for item in self.stats:
            s = f'{item} : {self.stats[item]}\n'
            print(s)

    def make_average_sfg(self, dppc=False):

        to_av = []
        for spec in self.sfg_spectra:
            if isinstance(spec.name, SystematicGasExName):
                if "deep" not in spec.name.type and "low" not in spec.name.type:
                    to_av.append(spec)

        if dppc == False:

            if len(to_av) != 0:
                out = to_av[0]
                for spec in to_av[1:]:
                    out += spec
                out.name.full_name = self.name
                out = out.calculate_ch_integral(average="gernot")

            else:
                out = None

        else:
            if len(to_av) != 0:
                dates = dppc
                temp= []

                for spec in to_av: # type: SfgSpectrum
                    dppc_integral = dates[spec.name.creation_time.date()]
                    ch_integral = spec.calculate_ch_integral(average="gernot")
                    temp.append(np.sqrt((ch_integral/dppc_integral)))

                out = np.average(temp), np.std(temp)
            else:
                out = None

        return out

    def get_overview(self, sfg=True, lt=True):
        """A function returning a formatted string of information about the station. Can be used to
        check wether everything was imported correctly"""
        separator = "-" * 80
        outstring = f'Station number {self.station_number} with name {self.name} '
        outstring += f'on {self.date[0:2]}.{self.date[2:]}.2018\n'
        outstring += f'The station contains {len(self.sfg_spectra)} SFG measurements '
        outstring += f' and  {len(self.lt_isotherms)} compression isotherm numbers.\n'

        if sfg == True:
            outstring += f'{separator}\n List of the SFG measurements:\n'
            for i, spectrum in enumerate(self.sfg_spectra):
                temp = f'{i+1}: {spectrum.name.full_name}\n'
                outstring += temp

        if lt == True:
            outstring += f'{separator}\n List of the LT measurements:\n'
            for i, spectrum in enumerate(self.lt_isotherms):
                temp = f'{i+1}: {spectrum.name}\n'
                outstring += temp
        outstring += separator + "\n"

        return outstring

    def get_doy(self):
        doy = self.date.timetuple().tm_yday
        factor = (1-int(self.station_number))*0.25
        return doy+factor

    def get_value_by_type(self, type, value):

        dic = { "a": ("s", "p", "c"), "deep":("c",), "sml":("p", "s"), "s":("s", ) , "p":("p",)}
        types = dic[type]
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
                        temp = np.sqrt(integral/dppc_integral)
                        out.append(temp)
                        
  
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

        #tension
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


class Spectrum:

    def __init__(self, name, x_data, y_data):

        self.name = name
        self.x_data = x_data
        self.y_data = y_data

        self.printstring = f'Spectrum object of name {self.name}'

    def __repr__(self):
        return self.printstring

    def __print__(self):
        return self.__repr__()

    def normalize(self, external=False):
        """Normalization function, either norming the y_data to its maximum or an external factor"""

        if external is False:
            factor = np.max(self.y_data)
            return self.y_data / factor

        else:
            return self.y_data / external

    def yield_spectral_range(self):
        """Returns the minimum and maximum and the length of the x_data of the spectrum"""

        return min(self.x_data), max(self.x_data), len(self.x_data)

    def smooth(self, points=9, order=5):
        """Apply a smooth routine to the y_data"""
        y = savgol_filter(self.y_data, points, order)
        return y

    def drop_ascii(self, delimiter=";"):
        """Create an ascii file with the wavenumbers and normalized intensities"""
        with open(self.name[:-4] + ".csv", "w") as outfile:
            writer = csv.writer(outfile, delimiter=delimiter)
            for i in zip(self.x_data, self.y_data):
                writer.writerow((i[0], i[1]))

    def slice_by_borders(self, upper, lower):
        diff = 10000000
        upper_index = 0
        lower_index = -1

        for index, spectrum in enumerate(self.x_data):

            temp_diff = abs(upper - self.x_data[index])
            if temp_diff < diff:
                diff = temp_diff
                upper_index = index

        diff = 10000000

        for index, spectrum in enumerate(self.x_data):

            temp_diff = abs(lower - self.x_data[index])
            if temp_diff < diff:
                diff = temp_diff
                lower_index = index

        return upper_index, lower_index


class SpectraManager:
    """A class to handle UV, IR, and Raman specta"""

    def __init__(self, database):
        self.raman = []
        self.ir = []
        self.uv = []

    def construct_spectrum(self, item):
        """Constructs a spectrum object from the SQL query results"""
        pass

    def query(self, table):
        """Extract the data from the table"""
        pass


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
        self.sample_hash = SampleHash(sample_hash)
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


# plotting functions

def scatter_maxpressure_day(isothermlist):
    """Create a scatter plot (day of cruise vs. maximum surface pressure) of the LtIsotherms
    provided in the isothermlist."""
    for isotherm in isothermlist:

        if isotherm.type == "p":
            color = "red"
        else:
            color = "blue"
        if "1" in isotherm.station:
            marker = "o"
        elif "2" in isotherm.station:
            marker = "x"
        elif "3" in isotherm.station:
            marker = "*"
        else:
            marker = "8"
        plt.scatter(isotherm.day, isotherm.get_maximum_pressure(), color=color, s=40, marker=marker)
        plt.text(isotherm.day + 0.1, isotherm.get_maximum_pressure() + 0.1,
                 isotherm.station + isotherm.type + isotherm.number,
                 fontsize=7)
        # plt.text(isotherm.day+0.2, isotherm.get_maximum_pressure()+0.1, isotherm.name, fontsize=7)

    plt.xlabel("days of cruise")
    plt.ylabel("maximum pressure/ mNm$^{-1}$")
    plt.grid()
    # legend_elements = [Line2D([0], [0], marker='o', color='r', label='glass plate', markerfacecolor='g', markersize=40)]
    # plt.legend(handles=legend_elements)
    plt.show()


def plot_vs_time(isothermlist):
    """Plots surface pressure vs time for each of the LtIsotherms provided in isothermlist."""
    for i in isothermlist:
        plt.plot(i.time, i.pressure)

    plt.xlabel("time")
    plt.ylabel("maximum pressure/ mNm$^{-1}$")
    plt.grid()
    plt.show()


def plot_per_sample(isothermlist):
    """Matches the items provided in the list of LtIsotherms according to their samples and exports all
     plots of a sample as pdf."""

    hashes = []
    for isotherm in isothermlist:

        if isotherm.create_sample_hash() not in hashes:
            hashes.append(isotherm.create_sample_hash())
            leg = isotherm.name
            plt.plot(isotherm.time, isotherm.pressure, label=leg)

            for partner in isotherm.partners:
                leg = partner.name
                plt.plot(partner.time, partner.pressure, label=leg)

            plt.xlabel("time/ s")
            plt.ylabel("maximum pressure/ mNm$^{-1}$")
            plt.title(str(isotherm.day) + " " + isotherm.station + " " + isotherm.type + " " + str(isotherm.number))
            plt.grid()
            plt.legend()
            plt.savefig(isotherm.create_sample_hash() + ".png")
            plt.cla()

        else:
            print("Already processed!")


def sfg_pub_plot(speclist, title="default", normalized="false"):
    """Produces a pre-formatted SFG plot from a list of SFG spectrum objects"""

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.set_xlabel("Wavenumber/ cm$^{-1}$", fontsize=10)
    ax.set_ylabel("Norm. SFG intensity/ arb. u.", fontsize=10)

    inc = 0.25 / len(speclist)
    counter = 0
    for spectrum in speclist:
        eff_alpha = 0.75 + inc * counter
        if normalized == "false":
            ax.plot(spectrum.wavenumbers, spectrum.normalized_intensity, linewidth=1.5, marker="o", markersize=3,
                    alpha=eff_alpha, label=spectrum.name.full_name)
        elif normalized == "true":
            ax.plot(spectrum.wavenumbers, spectrum.normalize_to_highest(), linewidth=1.5, marker="o", markersize=3,
                    alpha=eff_alpha, label=spectrum.name.full_name)
        counter += 1

    size = fig.get_size_inches()
    ratio = size[0] / size[1]
    fig.set_size_inches(3.2 * ratio, 3.2)
    fig.tight_layout()
    return fig


def sfg_stack_plot(speclist):
    """Stacks the provided SfgSpectrum objects after normalizing them to 1 internally."""
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.set_xlabel("Wavenumber/ cm$^{-1}$", fontsize=10)
    ax.set_ylabel("Norm. SFG intensity/ arb. u.", fontsize=10)
    ax.set_yticks([])
    ax.set_yticklabels([])

    inc = 0.25 / len(speclist)
    counter = 0
    offset = 0
    for spectrum in speclist:
        eff_alpha = 0.75 + inc * counter
        ax.plot(spectrum.wavenumbers, spectrum.normalize_to_highest() + offset, linewidth=1.5, marker="o", markersize=3,
                alpha=eff_alpha, label=spectrum.name.full_name)
        counter += 1
        offset += 0.2

    size = fig.get_size_inches()
    ratio = size[0] / size[1]
    fig.set_size_inches(3.2 * ratio, 3.2)
    fig.tight_layout()
    return fig


def sfg_doublestack_plot(speclist1, speclist2):
    """Two stacked-by-offset plots, one on the bottom, on on the top of the figure."""

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax1 = fig.add_subplot(211)
    ax2 = fig.add_subplot(212)
    ax2.set_xlabel("Wavenumber/ cm$^{-1}$", fontsize=10)
    ax.set_ylabel("Norm. SFG intensity/ arb. u.", fontsize=10)

    ax.spines['top'].set_color('none')
    ax.spines['bottom'].set_color('none')
    ax.spines['left'].set_color('none')
    ax.spines['right'].set_color('none')
    ax.tick_params(labelcolor='w', top=False, bottom=False, left=False, right=False)

    ax2.set_yticks([])
    ax2.set_yticklabels([])
    ax1.set_yticks([])
    ax1.set_yticklabels([])

    inc = 0.25 / len(speclist1)
    counter = 0
    offset = 0
    for spectrum in speclist1:
        eff_alpha = 0.75 + inc * counter
        ax1.plot(spectrum.wavenumbers, spectrum.normalize_to_highest() + offset, linewidth=1.5, marker="o",
                 markersize=3,
                 alpha=eff_alpha, label=spectrum.name.full_name)
        counter += 1
        offset += 0.2

    inc = 0.25 / len(speclist1)
    counter = 0
    offset = 0
    for spectrum in speclist2:
        eff_alpha = 0.75 + inc * counter
        ax2.plot(spectrum.wavenumbers, spectrum.normalize_to_highest() + offset, linewidth=1.5, marker="o",
                 markersize=3,
                 alpha=eff_alpha, label=spectrum.name.full_name)
        counter += 1
        offset += 0.2
    ax1.legend()
    ax2.legend()
    # size = fig.get_size_inches()
    # ratio = size[0] / size[1]
    # fig.set_size_inches(3.2 * ratio, 3.2)
    # fig.tight_layout()
    return fig


def finalize_figure(fig, title="test2"):
    """Makes figures publication-ready and exports them as pdf. Takes the title for the output file as argument"""
    size = fig.get_size_inches()
    ratio = size[0] / size[1]
    fig.set_size_inches(3.2 * ratio, 3.2)
    fig.tight_layout()
    fig.savefig(title + ".pdf", dpi=600)


def lt_plot_stats_new(S):
    """Takes Session controll manager as argument. Performs the bar/max surface pressure plot plot for the LT
    isotherms"""
    S.setup_for_gasex()

    for s in S.stations.values():
        s.join_samples()
        s.count_per_type()

    fig, (ax1, u_ax) = plt.subplots(nrows=2, ncols=1)
    ax1.set_xlabel("station number")
    ax1.set_ylabel("average surface pressure/ mN/m")
    ax2 = ax1.twinx()

    ax2.set_xlabel("day of cruise")
    ax2.set_ylabel("positive samples/ percent")

    ax1.set_title("Surfactant occurence in GasEx 1 (June '18), \nmeasured by Langmuir Trough\n\n\n")
    ax1.grid(True)

    days = [i for i in range(1, 15)]
    ax3 = u_ax.twiny()

    ax3.xaxis.set_ticks_position("bottom")
    ax3.xaxis.set_label_position("bottom")
    ax3.spines["bottom"].set_position(("axes", -0.35))
    ax3.set_xlabel("day of cruise")

    u_ax.set_xlabel("station number")
    u_ax.set_ylabel("Norm. SFG intensity/ arb.u.")
    u_ax.grid(True)
    u_ax.set_ylim(-0.0001, 0.0006)

    p_percentages = []
    s_percentages = []
    t_percentages = []

    stations = []
    plates = []
    totals = []
    screens = []

    p_std = []
    s_std = []
    t_std = []

    for s in S.stations.values():

        if len(s.lt_isotherms) > 0:
            station = s.station_number
            stations.append(station)
            screens.append(s.stats["screen_av"])
            plates.append(s.stats["plate_av"])
            totals.append(s.stats["total_av"])

            t_std.append(s.stats["std_total"])
            s_std.append(s.stats["std_screen"])
            p_std.append(s.stats["std_plate"])

            ax3.scatter(s.cruise_day, s.stats["screen_av"], s=0)
            s_percentages.append([station, s.stats["percent_screen"]])
            p_percentages.append([station, s.stats["percent_plate"]])
            t_percentages.append([station, s.stats["total_percent"]])

        max_ins = []

        for spectrum in s.sfg_spectra:  # type: SfgSpectrum

            station = s.station_number
            b = spectrum.slice_by_borders(3000, 2800)
            max = np.max(spectrum.normalized_intensity[b[0]:b[1] + 1])
            max_ins.append(max)

        max_ins = np.array(max_ins)
        av = np.average(max_ins)
        std = np.std(max_ins)
        u_ax.errorbar(station, av, std, alpha=0.9, fmt="o", color="b", label="average", barsabove="true", capsize=5,
                      capthick=2)

    ax1.errorbar(stations, totals, yerr=t_std, fmt="o", color="r", barsabove="true", capsize=5, capthick=2)
    rects2 = ax2.bar([a[0] + 0.2 for a in t_percentages], [a[1] for a in t_percentages], alpha=0.45, label="total")

    ax2.legend()
    # u_ax.legend()
    plt.tight_layout()
    plt.show()


def sfg_with_lt(S):
    """Takes Session controll manager as argument. Performs the bar/max surface pressure plot plot for the LT
        isotherms"""
    S.setup_for_gasex()

    for s in S.stations.values():
        s.join_samples()
        s.count_per_type()

    fig, ax1 = plt.subplots()
    ax1.set_xlabel("station number")
    ax1.set_ylabel("norm. SFG intensity/ arb. u.")
    ax2 = ax1.twinx()

    days = [i for i in range(1, 15)]
    ax3 = ax1.twiny()
    ax3.set_xlabel("day of cruise")
    ax3.xaxis.set_ticks_position("bottom")
    ax3.xaxis.set_label_position("bottom")
    ax3.spines["bottom"].set_position(("axes", -0.25))

    ax2.set_ylabel("positive samples/ percent")
    ax1.set_title("Surfactant occurence in GasEx 1 (June '18), \nmeasured by Langmuir Trough and SFG\n\n")
    ax1.grid(True)

    p_percentages = []
    s_percentages = []
    t_percentages = []

    for s in S.stations.values():

        station = s.station_number

        for spectrum in s.sfg_spectra:  # type: SfgSpectrum

            b = spectrum.slice_by_borders(3000, 2800)
            max = np.max(spectrum.normalized_intensity[b[0]:b[1] + 1])
            mapping = {"p": ["g", "plate"],
                       "s1": ["b", "screen"],
                       "s2": ["b", "screen"],
                       "s3": ["b", "screen"],
                       "deep": "black",
                       "low": "black"}
            try:
                ax1.scatter(station, max, alpha=0.9, color=mapping[spectrum.name.type][0],
                            label=mapping[spectrum.name.type][1],
                            s=60)
            except KeyError:
                ax1.scatter(station, max, alpha=0.9, color="black", label="CTD", s=60)

        ax3.scatter(s.cruise_day, max, s=0)
        ax1.set_ylim(0, 0.0008)
        # ax1.legend()

        s_percentages.append([station, s.stats["percent_screen"]])
        p_percentages.append([station, s.stats["percent_plate"]])
        t_percentages.append([station, s.stats["total_percent"]])

    rects1 = ax2.bar([a[0] - 0.2 for a in p_percentages], [a[1] for a in p_percentages], alpha=0.35, width=0.2,
                     color="g", label="plate, trough (positive)")
    rects2 = ax2.bar([a[0] for a in s_percentages], [a[1] for a in s_percentages], alpha=0.35, width=0.2, color="b",
                     label="screen, trough (positive)")
    rects2 = ax2.bar([a[0] + 0.2 for a in t_percentages], [a[1] for a in t_percentages], alpha=0.35, width=0.2,
                     color="r", label="total, trough (positive)")

    # ax2.legend()

    h1 = [Line2D([0], [0], marker='o', color='w', markerfacecolor='g', markersize=10),
          Line2D([0], [0], marker='o', color='w', markerfacecolor='b', markersize=10),
          Line2D([0], [0], marker='o', color='w', markerfacecolor='black', markersize=10)]

    l1 = ["plate, SFG", "screen, SFG", "CTD, SFG"]
    h2, l2 = ax2.get_legend_handles_labels()
    ax1.legend(h1 + h2, l1 + l2, loc="upper center", ncol=2).draggable()
    plt.sca(ax1)
    plt.tight_layout()
    plt.show()


def lt_sfg_integral(S):
    """Takes Session controll manager as argument. Performs the bar/max surface pressure plot plot for the LT
    isotherms"""

    S.setup_for_gasex()

    for s in S.stations.values():
        s.join_samples()
        s.count_per_type()

    fig, (ax1, u_ax) = plt.subplots(nrows=2, ncols=1)
    ax1.set_xlabel("station number")
    ax1.set_ylabel("average surface pressure/ mN/m")
    ax2 = ax1.twinx()

    ax2.set_xlabel("day of cruise")
    ax2.set_ylabel("positive samples/ percent")

    ax1.set_title("Surfactant occurrence in GasEx 1 (June '18)\n\n")
    ax1.grid(True)

    days = [i for i in range(1, 15)]
    ax3 = u_ax.twiny()

    ax3.xaxis.set_ticks_position("bottom")
    ax3.xaxis.set_label_position("bottom")
    ax3.spines["bottom"].set_position(("axes", -0.35))
    ax3.set_xlabel("day of cruise")

    u_ax.set_xlabel("station number")
    u_ax.set_ylabel("Integrated SFG CH intensity/ arb.u.")
    u_ax.grid(True)
    u_ax.set_ylim(-0.00025, 0.014)

    p_percentages = []
    s_percentages = []
    t_percentages = []

    stations = []
    plates = []
    totals = []
    screens = []

    p_std = []
    s_std = []
    t_std = []

    for s in S.stations.values():

        if len(s.lt_isotherms) > 0:
            station = s.station_number
            stations.append(station)
            screens.append(s.stats["screen_av"])
            plates.append(s.stats["plate_av"])
            totals.append(s.stats["total_av"])

            t_std.append(s.stats["std_total"])
            s_std.append(s.stats["std_screen"])
            p_std.append(s.stats["std_plate"])

            ax3.scatter(s.cruise_day, s.stats["screen_av"], s=0)
            s_percentages.append([station, s.stats["percent_screen"]])
            p_percentages.append([station, s.stats["percent_plate"]])
            t_percentages.append([station, s.stats["total_percent"]])

        max_ins = []

        if len(s.sfg_spectra) > 0:

            temp = []

            for spec in s.sfg_spectra:
                temp.append(spec.calculate_ch_integral())

            average = np.average(temp)
            std = np.std(temp)
            u_ax.errorbar(s.station_number, average, yerr=std, fmt="o", color="b", barsabove="true", capsize=5,
                          capthick=2)

    ax1.errorbar(stations, totals, yerr=t_std, fmt="o", color="r", barsabove="true", capsize=5, capthick=2)
    rects2 = ax2.bar([a[0] + 0.2 for a in t_percentages], [a[1] for a in t_percentages], alpha=0.45, label="positive")

    ax2.legend()
    # u_ax.legend()
    plt.tight_layout()
    plt.show()


def lt_sfg_integral_dppc(S):
    """Takes Session controll manager as argument. Performs the bar/max surface pressure plot plot for the LT
    isotherms"""

    S.setup_for_gasex()
    figure, plots = plt.subplots(3, 1)
    t_percentages = []



    ltplot = plots[0]
   # ltplot.set_title("Surfactant occurrence in GasEx 1 (June '18)\n\n", fontweight='bold')
    ltplot.set_ylabel("Average surface pressure/\n mN/m")
    ltplot.xaxis.set_ticklabels([])

    bars = ltplot.twinx()
    bars.set_ylabel("Positive samples/\n percent")

    sfg = plots[1]
    sfg.set_ylim(-0.00025, 0.014)
    sfg.set_ylabel("Integrated SFG \nCH intensity/ arb.u.")
    sfg.xaxis.set_ticklabels([])

    dppc = plots[2]
    dppc.set_ylabel("Surface coverage/\n %")
    dppc.set_xlabel("\nStation number")

    base = ltplot.twiny()
    base.xaxis.set_ticks_position("top")
    base.xaxis.set_label_position("top")
    #base.spines["top"].set_position(("axes", +0.7))
    base.set_xlabel("Day of cruise\n")

    for station in S.stations.values(): # type: Station

        if station.cruise_day < 13:
            base.scatter(station.cruise_day, station.stats["total_av"], s=0)


        if len(station.lt_isotherms) > 0:
            ltplot.scatter(station.station_number, station.stats["total_av"], color="green")
            ltplot.text(station.station_number+0.15,station.stats["total_av"]+0.15, str(station.isotherm_count))
            t_percentages.append([station.station_number, station.stats["total_percent"]])



        if len(station.sfg_spectra) > 0:
            if station.make_average_sfg() is not None:
                sfg.scatter(station.station_number, station.make_average_sfg(), color="blue")

            if station.make_average_sfg(dppc=S.dppc_ints) is not None:
                dppc.scatter(station.station_number, station.make_average_sfg(dppc=S.dppc_ints)[0]*100, color="red")

    bars.bar([a[0] + 0.2 for a in t_percentages], [a[1] for a in t_percentages], alpha=0.45, label="positive")
    plt.show()


def sfg_plot_broken_axis(speclist, lower, upper, title="default", normalized="false"):
    """Produces a pre-formatted SFG plot from a list of SFG spectrum objects"""

    fig, (ax, ax2) = plt.subplots(1, 2, sharey=True)
    ax.set_xlabel("Wavenumber/ cm$^{-1}$", fontsize=10)
    ax.set_ylabel("Norm. SFG intensity/ arb. u.", fontsize=10)

    ax.set_xlim(lower[0], lower[1])
    ax2.set_xlim(upper[0], upper[1])

    ax.spines['right'].set_visible(False)
    ax2.spines['left'].set_visible(False)
    ax.yaxis.tick_left()
    ax.tick_params(labeltop='off')  # don't put tick labels at the top
    ax2.yaxis.tick_right()

    d = .015  # how big to make the diagonal lines in axes coordinates
    # arguments to pass plot, just so we don't keep repeating them
    kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
    ax.plot((1 - d, 1 + d), (-d, +d), **kwargs)
    ax.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)

    kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
    ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)
    ax2.plot((-d, +d), (-d, +d), **kwargs)

    inc = 0.25 / len(speclist)
    counter = 0
    for spectrum in speclist:
        eff_alpha = 0.75 + inc * counter
        if normalized == "false":
            ax.plot(spectrum.wavenumbers, spectrum.normalized_intensity, linewidth=1.5, marker="o", markersize=3,
                    alpha=eff_alpha, label=spectrum.name.full_name)
            ax2.plot(spectrum.wavenumbers, spectrum.normalized_intensity, linewidth=1.5, marker="o", markersize=3,
                     alpha=eff_alpha, label=spectrum.name.full_name)
        elif normalized == "true":
            ax.plot(spectrum.wavenumbers, spectrum.normalize_to_highest(), linewidth=1.5, marker="o", markersize=3,
                    alpha=eff_alpha, label=spectrum.name.full_name)
        counter += 1

    size = fig.get_size_inches()
    ratio = size[0] / size[1]
    fig.set_size_inches(3.2 * ratio, 3.2)
    fig.tight_layout()
    return fig


def lt_integral_average(S):
    """Takes Session controll manager as argument. Performs the bar/max surface pressure plot plot for the LT
    isotherms"""
    S.setup_for_gasex()

    for s in S.stations.values():
        s.join_samples()
        s.count_per_type()

    fig, (ax1, u_ax) = plt.subplots(nrows=2, ncols=1)
    ax1.set_xlabel("station number")
    ax1.set_ylabel("average surface pressure/ mN/m")
    ax2 = ax1.twinx()

    ax2.set_xlabel("day of cruise")
    ax2.set_ylabel("positive samples/ percent")

    ax1.set_title("Surfactant occurrence in GasEx 1 (June '18)\n\n")
    ax1.grid(True)

    days = [i for i in range(1, 15)]
    ax3 = u_ax.twiny()

    ax3.xaxis.set_ticks_position("bottom")
    ax3.xaxis.set_label_position("bottom")
    ax3.spines["bottom"].set_position(("axes", -0.35))
    ax3.set_xlabel("day of cruise")

    u_ax.set_xlabel("station number")
    u_ax.set_ylabel("Integrated SFG CH intensity/ arb.u.")
    u_ax.grid(True)
    u_ax.set_ylim(-0.00025, 0.014)

    p_percentages = []
    s_percentages = []
    t_percentages = []

    stations = []
    plates = []
    totals = []
    screens = []

    p_std = []
    s_std = []
    t_std = []

    for s in S.stations.values():  # type: Station

        if len(s.lt_isotherms) > 0:
            station = s.station_number
            stations.append(station)
            screens.append(s.stats["screen_av"])
            plates.append(s.stats["plate_av"])
            totals.append(s.stats["total_av"])

            t_std.append(s.stats["std_total"])
            s_std.append(s.stats["std_screen"])
            p_std.append(s.stats["std_plate"])

            ax3.scatter(s.cruise_day, s.stats["screen_av"], s=0)
            s_percentages.append([station, s.stats["percent_screen"]])
            p_percentages.append([station, s.stats["percent_plate"]])
            t_percentages.append([station, s.stats["total_percent"]])

        max_ins = []

        if len(s.sfg_spectra) > 0:

            av_spec = s.make_average_sfg()
            if isinstance(av_spec, AddedSpectrum):
                integral = av_spec.calculate_ch_integral()
                u_ax.scatter(s.station_number, integral)

    ax1.errorbar(stations, totals, yerr=t_std, fmt="o", color="r", barsabove="true", capsize=5, capthick=2)
    rects2 = ax2.bar([a[0] + 0.2 for a in t_percentages], [a[1] for a in t_percentages], alpha=0.45, label="positive")

    ax2.legend()
    # u_ax.legend()
    plt.tight_layout()
    plt.show()


def baseline_demo(spectrum, name="default"):
    spectrum.correct_baseline(average="gernot")

    test = np.linspace(2750, 3050, 10000)
    func = spectrum.make_ch_baseline(average="gernot")
    borders = spectrum.slice_by_borders(3000, np.min(spectrum.wavenumbers))

    f, axarr = plt.subplots(2, sharex=True)
    axarr[0].plot(spectrum.wavenumbers, spectrum.normalized_intensity, label=spectrum.name.full_name, linewidth=1.5,
                  marker="o", markersize=3)
    axarr[0].plot(test, func(test), color="r", label="Baseline")
    axarr[0].set_xlabel("Wavenumber/ cm$^{-1}$")
    axarr[0].set_ylabel("Norm. SFG intensity/ arb. u.")
    #axarr[0].set_title("Demonstration of the automatic baseline subtraction and integration")
    axarr[0].legend()

    axarr[1].plot(spectrum.wavenumbers, spectrum.baseline_corrected, label=spectrum.name.full_name, linewidth=1.5,
                  marker="o", markersize=3)
    axarr[1].fill_between(spectrum.wavenumbers[borders[0]:borders[1] + 1],
                          spectrum.baseline_corrected[borders[0]:borders[1] + 1])
    axarr[1].set_xlabel("Wavenumber/ cm$^{-1}$")
    axarr[1].set_ylabel("Norm. SFG intensity/ arb. u.")
    axarr[1].legend()
    name = spec.name.full_name
    plt.savefig(name + ".png")


def benchmark_baseline(speclist):
    for spectrum in speclist:  # type: SfgSpectrum

        standard = spectrum.calculate_ch_integral()
        regress = spectrum.calculate_ch_integral(average="min_reg")
        gernot = spectrum.calculate_ch_integral(average="gernot")
        return standard, regress, gernot


def plot_lt_isotherm(isotherm): # type: LtIsotherm
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.set_xlabel("Total area (cm$^{-2}$)\n", fontsize=20)
    ax.set_ylabel("Surface pressure (mN/m)\n", fontsize=20)
    ax.set_ylim(-0.1, 2)
    ax.plot(isotherm.area, isotherm.pressure, linewidth=3)
    plt.show()

def broken_axis(x, y, lim):

    fig, (ax, ax2) = plt.subplots(1, 2, sharey=True)

    ax.set_ylabel("Surface tension/ $mN \cdot m^{-1}$")
    fig.text(0.5, 0.015, s="Day of the year", ha="center", va="center", size=14)

    ax.set_xlim(lim[0], lim[1])
    ax2.set_xlim(lim[2], lim[3])

    ax.spines['right'].set_visible(False)
    ax2.spines['left'].set_visible(False)
    ax.yaxis.tick_left()
    ax.tick_params(labeltop='off')  # don't put tick labels at the top
    ax2.yaxis.tick_right()

    d = .012  # how big to make the diagonal lines in axes coordinates
    # arguments to pass plot, just so we don't keep repeating them
    kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
    ax.plot((1 - d, 1 + d), (-d, +d), **kwargs)
    ax.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)

    kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
    ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)
    ax2.plot((-d, +d), (-d, +d), **kwargs)

    ax.scatter(x, y, marker="o")
    ax2.scatter(x, y, marker="o")



def broken_axis_errorbar(lim):
    fig, (ax, ax2) = plt.subplots(1, 2, sharey=True)

    ax.set_ylabel("Surface tension/ $mN \cdot m^{-1}$")
    fig.text(0.5, 0.015, s="Day of the year", ha="center", va="center", size=14)

    ax.set_xlim(lim[0], lim[1])
    ax2.set_xlim(lim[2], lim[3])

    ax.spines['right'].set_visible(False)
    ax2.spines['left'].set_visible(False)
    ax.yaxis.tick_left()
    ax.tick_params(labeltop='off')  # don't put tick labels at the top
    ax2.yaxis.tick_right()

    d = .012  # how big to make the diagonal lines in axes coordinates
    # arguments to pass plot, just so we don't keep repeating them
    kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
    ax.plot((1 - d, 1 + d), (-d, +d), **kwargs)
    ax.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)

    kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
    ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)
    ax2.plot((-d, +d), (-d, +d), **kwargs)

    return ax, ax2


def plot_stats(stations, stats, ylabel="Surface tension/ $mN \cdot m^{-1}$"):

    axes = broken_axis_errorbar([153, 166, 254, 266])
    counter = 0
    colormap = { 0: "ro", 1: "bo", 2:"go", 3:"purple"}
    axes[0].set_ylabel(ylabel)

    for s in stats:

        doy = []
        out = []
        err = []

        for station in stations:

            if station.stats[s] is not None:
                doy.append(station.get_doy())
                out.append(station.stats[s][0])
                err.append(station.stats[s][1])


        axes[0].errorbar(doy, out, yerr=err, fmt=colormap[counter], color="r", barsabove="true", capsize=5,
                         capthick=1, ecolor="black", elinewidth=1.0,
                     markeredgecolor="black", markeredgewidth=0.4, antialiased=True)
        axes[1].errorbar(doy, out, yerr=err, fmt=colormap[counter], color="r", barsabove="true", capsize=5, capthick=1, ecolor="black", elinewidth=1.0,
                      markeredgecolor="black",  markeredgewidth=0.4, antialiased=True, label=s)
        counter += 1

    axes[1].legend()
    plt.show()


def tension_average(station):
        average = []
        x = station.get_doy()

        for tension in station.tensions:
            average.append(tension[1])

        average = np.array(average)
        av_out = np.average(average)
        std = np.std(average)



        return x, av_out, std


plt.style.use('seaborn-talk')
import matplotlib as mpl
mpl.rcParams['axes.linewidth']= 2


S = SessionControlManager("sfg.db", "test")
S.setup_for_gasex()
limits = [153, 166, 254, 266]

l1 = "SFG CH integral/ arb u."
l2 = "surface coverage"
l3 = "surface tension/ $mN \cdot m^{-1}$"
l4 = "surface pressure/ $mN \cdot m^{-1}$"


plot_stats(S.stations.values(), ["pressure_deep","pressure_sml"], ylabel=l4)

# axes[0].errorbar(doys1, av_deeps, yerr=av_deep_stds, fmt="ro", color="r", barsabove="true", capsize=5, capthick=1, ecolor="black", elinewidth=1.0,
#             markeredgecolor="black", markeredgewidth=0.4, antialiased=True)
# axes[1].errorbar(doys1, av_deeps, yerr=av_deep_stds, fmt="ro", color="r", barsabove="true", capsize=5, capthick=1, ecolor="black", elinewidth=1.0,
#              markeredgecolor="black",  markeredgewidth=0.4, antialiased=True, label="deep water")
#
# axes[0].errorbar(doys2, av_smls, yerr=av_sml_stds, fmt="bo", color="r", barsabove="true", capsize=5, capthick=1, ecolor="black", elinewidth=1.0,
#             markeredgecolor="black", markeredgewidth=0.4, antialiased=True)
# axes[1].errorbar(doys2, av_smls, yerr=av_sml_stds, fmt="bo", color="r", barsabove="true", capsize=5, capthick=1, ecolor="black", elinewidth=1.0,
#              markeredgecolor="black",  markeredgewidth=0.4, antialiased=True, label="SML")







