#module-internal imports
from new_gui import run_app, run_lt_app

#standard utilities
import os
import shutil
import csv
import time
import datetime

#scientific libraries
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.lines import Line2D
import sqlite3
from scipy.signal import savgol_filter
from scipy.integrate import simps as sp

#for session controll manager
#from new_gui import run_app

rcParams['mathtext.default'] = 'regular'


# noinspection PySimplifyBooleanCheck,PySimplifyBooleanCheck


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

    def normalize_to_highest(self, intensity="default", external_norm="none"):
        """normalize an given array to its maximum, typically the normalized or raw intensity"""
        if intensity == "default":
            intensity = self.normalized_intensity
        if external_norm == "none":
            norm_factor = np.max(intensity)
        else:
            norm_factor = external_norm

        return (intensity / norm_factor)

    def yield_peaklist(self, mode="norm"):

        out = []
        tup = self.detailed_analysis(threshold=1.5, intensity= mode)
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

    def smooth(self, points=9, order=5):
        """Apply a smooth routine to the normalized_intensity"""
        y = savgol_filter(self.normalized_intensity, points, order)
        return y

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
            border = np.average(slopes)*1/6

            # check for left border
            while slopes[k] > border and k >= 0:
                k -= 1
            left = k

            # check for right border
            k = i + 1

            if k >= len(slopes):
                k = len(slopes)-1
            else:
                #if traversing to the right does not find a proper peak ending
                try:
                    while (slopes[k] < border):
                        k += 1
                except IndexError:
                        k = k-1
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

            data_out.append((center, left, right, center_intensity, peak_slice_x, peak_slice_y, datapoints, area, indices))

        #sort peaks by peak intensity
        data_out = sorted(data_out, key=(lambda x: x[3]), reverse=True)
        return data_out

    def integrate_peak(self, x_array, y_array):
        """Numpy integration routine for numerical peak integration with the trapezoidal rule"""
        try:
            area = sp(y_array, x_array)
            return area
        except:
            return "Area  could not be calculated"

    def drop_ascii(self):
        """Create an ascii file with the wavenumbers and normalized intensities"""
        with open(self.name.full_name[:-4]+".csv" , "w") as outfile:
            writer = csv.writer(outfile, delimiter=";")
            for i in zip(self.wavenumbers,self.normalized_intensity):
                writer.writerow((i[0],i[1]))

    def drop_tex_peaktable(self, threshold=1.5):

        peaks = self.detailed_analysis(threshold=threshold)
        tablestring = ""

        counter = 1
        for peak in peaks:
            tablestring += "\subcaption*{" + "Peak " + str(counter) + "}\n"
            tablestring += "\\begin{tabular}{|c|c|c|c|c|c|}\n\hline\n"
            tablestring += "Wavenumber & normalized intensity & area & left border & right border"+"\\\\"+"\hline\n"

            tablestring += str(peak[0]) + " & " + str(peak[3]) + " & " + str(peak[7]) + " & " + str(
                peak[1]) + " & " + str(peak[2]) + "\\\\"
            tablestring += "\n\hline\n\end{tabular}\n"
            counter += 1

        return tablestring

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
                concentration = input("Enter surf stock concentration of spectrum "+self.name.full_name+": \n")


            volume = float(self.name.surfactant_spread_volume)

            concentration = float(concentration) * 10 ** -3  # conversion in mol per liter
            volume = volume * 10 ** -6  # conversion in liter
            amount = volume * concentration

        else:
            concentration = self.name.surf_stock_concentration
            if self.name.surf_stock_concentration == "unknown":
                concentration = input(
                    "Enter surf stock concentration of spectrum " + self.name.full_name+": \n")

            volume = float(self.name.surfactant_spread_volume)
            concentration = float(concentration) * 10 ** -3  # conversion in mol per liter
            volume = volume * 10 ** -6  # conversion in liter
            amount_su = volume * concentration

            sens_stock_conc = input("Enter sens stock concentration of spectrum " + self.name.full_name+": \n")
            concentration = float(sens_stock_conc) * 10 ** -3  # conversion in mol per liter
            volume = float(self.name.sensitizer_spread_volume) * 10 ** -6  # conversion in liter
            amount_se = volume * concentration

            amount = amount_se + amount_su

        molecules = (6.022 * 10 ** 23) * amount  # number of molecules
        area_per_molecule = self.calc_dish_area()/ molecules

        return area_per_molecule

    def yield_maximum(self):
        return np.max(self.normalized_intensity)

    def create_pointlist(self, y_array):

        output = []
        for i, (a, b) in enumerate(zip(self.wavenumbers[::-1], y_array)):
            output.append((a, b, i))

        return output

    def root(self):
        return np.sqrt(self.normalized_intensity)

    def get_sample_hash(self):
        return self.name.date+self.name.surfactant+self.name.sensitizer+self.name.surfactant_spread_volume\
               +str(self.name.sample_number)+self.name.comment

    def slice_by_borders(self, upper, lower):

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

    def make_ch_baseline(self):

        left = np.nonzero(self.wavenumbers == 2750)[0][0]

        l_interval = self.slice_by_borders(2800, 2750)
        l_interval_wl = self.wavenumbers[l_interval[0]:l_interval[1] + 1]
        l_interval = self.normalized_intensity[l_interval[0]:l_interval[1] + 1]

        interval = self.slice_by_borders(2960, 2895)
        interval_wl = self.wavenumbers[interval[0]:interval[1] + 1]
        interval = self.normalized_intensity[interval[0]:interval[1]+1]

        min_index = np.argmin(interval)
        l_min_index = np.argmin(l_interval)


        #av_intens = self.normalized_intensity[a1]+self.normalized_intensity[a2]+self.normalized_intensity[a3]
        #av_wn = self.wavenumbers[a1]+self.wavenumbers[a2]+self.wavenumbers[a3]
        #av_intens /=3
        #av_wn /=3

        slope = (interval[min_index] - l_interval[l_min_index])/(interval_wl[min_index]-l_interval_wl[l_min_index])
        intercept = l_interval[l_min_index]-slope*l_interval_wl[l_min_index]



        baseline = lambda x: slope*x+intercept
        return baseline

    def correct_baseline(self):

        func = self.make_ch_baseline()
        temp = self.normalized_intensity

        for i in range(2750,2960):
            index = np.nonzero(self.wavenumbers == i)
            correction = func(self.wavenumbers[index])
            temp[index] =temp[index]-correction

        self.baseline_corrected = temp







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
            Added.speccounter = self.speccounter+1
            return Added

        elif isinstance(SFG2, AddedSpectrum):
            names = self.name.full_name+"_"+SFG2.name.full_name

            sensitizers = self.name.sensitizer
            for i in SFG2.name.sensitizer:
                if i not in sensitizers:
                    sensitizers.append(i)

            surfactants = self.name.surfactants
            for i in SFG2.name.sensitizer:
                if i not in surfactants:
                    surfactants.append(i)

            name = AddedName(names, sensitizers, surfactants)
            Added = AddedSpectrum((new_wavenumbers, new_intensities),name)
            Added.speccounter += self.speccounter
            Added.speccounter += SFG2.speccounter
            return Added


class SystematicName:
    """A class extracting metainformation from the string filename and storing it in rather self-explanatory
    variables.(eg date of measurement, type of surfactant...)"""

    def __init__(self, namestring, creation_time="unknown"):

        self.refpath = "name_info/"

        try:
            with open(self.refpath+"Surfactants.txt") as outfile:
                pass
        except FileNotFoundError:
            self.refpath = "../name_info/"

        self.creation_time = creation_time
        # load the allowed surfactans and sensitizers from files
        self.Surfactants = {}
        self.Sensitizers = {}

        with open(self.refpath+"Surfactants.txt", "r") as infile:
            for line in infile:
                collect = line.split(":")
                self.Surfactants[collect[0]] = collect[1].strip()

        with open(self.refpath+"Sensitizers.txt", "r") as infile:
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
                if self.is_number(i[index-1]):
                    self.surf_stock_concentration = float(i[index-1])

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
            self.station = temp[1]+"_"+temp[2]
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


class AddedName(SystematicName):
    """This class is derived from the SfgSpectrum and represents the result of the addition of Sfg intensities."""
    def __init__(self,names):
        self.full_name = ("_").join(names)


class Analyzer:
    """This class takes, a list of SFG spectra as constructor argument. Its purpose
    is to perform analytical tasks with the spectral data, e.g. compare peaks, integral, datapoints
    and will be extended to handle statistics in the future"""

    def __init__(self, speclist):
        self.speclist = speclist

    def count_peak_abundance(self, number):
        peaklist = []

        for spectrum in self.speclist:
            temp = spectrum.yield_peaklist()
            if len(temp) > number:
                temp = temp[:number]
            for peak in temp:
                peaklist.append(peak)

        abundance = {}

        for peak in peaklist:

            if peak not in abundance:
                abundance[peak] = 1
            else:
                abundance[peak] += 1

        #calculate total number

        for key in abundance:
            abundance[key] /= len(self.speclist)

        return abundance

    def write_analysis(self):
        for spectrum in self.speclist:
            with open(spectrum.name.full_name+".ana", "w") as outfile:
                outfile.write("*"*80+"\n")
                outfile.write("Analytics output: "+"\n")
                tup = spectrum.detailed_analysis()
                counter = 1
                for peak in tup:
                    outfile.write("Peak "+str(counter)+":\n")
                    outfile.write("Wavenumber: " + str(peak[0]))
                    outfile.write("\nIntensity: "+str(peak[3]))
                    outfile.write("\nLeft border: "+str(peak[1]))
                    outfile.write("\nRight border: "+str(peak[2]))
                    outfile.write("\nIntegral: "+str(peak[7])+"\n")
                    outfile.write("-"*80+"\n")
                    counter += 1

    def intensity_vs_surfcoverage(self):

        out= []
        for spectrum in self.speclist:
            maxi = 0
            temp = spectrum.detailed_analysis()

            for peak in temp:
                if maxi < peak[3]:
                    maxi = peak[3]

            out.append((spectrum.calc_area_per_molecule(), maxi))

        return out

    def average_peak(self, number):

        out = []

        for i in range(number):
            wavenumber = []
            intensity = []
            counter = 0

            for spec in self.speclist:
                try:
                    peaks = spec.detailed_analysis()[i]
                except IndexError:
                    break

                wavenumber.append(peaks[0])
                intensity.append(peaks[3])

                counter += 1
                #print("found for Peak"+str(i)+":" + str(peaks[0]), str(peaks[3]))

            wavenumber = np.average(np.array(wavenumber))
            intensity = np.average(np.array(intensity))

            out.append((wavenumber, intensity))

        return out


class SessionControlManager:

    """Successor of the IpyInterpreter. A class to access all experimental data, search data by certain match criteria
    and produce plots efficiently. Designed to work with a sqllite database containing the spectral raw data. Especially
    useful from interactive python environments (eg IPy)."""

    # setup methods
    # todo: sort the functions under categories, now it is a mess!
    def __init__(self, database, id):

        self.db = sqlite3.connect(database, detect_types=sqlite3.PARSE_DECLTYPES|sqlite3.PARSE_COLNAMES)
        self.cur = self.db.cursor()
        self.table = "sfg_database"
        self.gasex_table = "sfg_gasex"
        self.session_id = id
        self.Surfactants = {}
        self.Sensitizers = {}
        self.Makros = {}
        self.stations = None
        self.get_senssurf_names()

        self.lt_manager = None

        #former IpyInterpreter functionality, tracking the primary key in parallel
        self.subset_ids = []
        self.subset = []

        self.recover_ids = []
        self.recover = []
        self.gasex_included = False

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

    # database methods

    def fetch_single(self, number, condition3=None, condition4=None, default_data=True):
        """Fetches the data of a single SFG spectrum from the database. Returns an SfgSpectrum object. The
        conditional parameters allow a selection of the spectrum according to specified properties. The
        default_data kwarg controlls which database the spectrum should be extracted from."""

        if default_data is True:
            database = self.table
        else:
            database = self.gasex_table

        command = "SELECT * FROM "+database+" WHERE ROWID="+str(number)

        if condition3 is not None:
            command += " AND "+str(condition3)+"="+str(condition4)

        self.cur.execute(command)

        result = self.cur.fetchall()[0]

        creationtime = result[2]

        if default_data is True or "DPPC" in result[1]:

            if "DPPC" not in result[1]:
                sysname = SystematicName(result[1], creationtime)
                wavenumber = np.array(result[-7].split(";")).astype(np.float)
                sfg = np.array(result[-6].split(";")).astype(np.float)
                ir = np.array(result[-5].split(";")).astype(np.float)
                vis = np.array(result[-4].split(";")).astype(np.float)
            else:
                sysname = SystematicName(result[1], creationtime)
                wavenumber = np.array(result[-4].split(";")).astype(np.float)
                sfg = np.array(result[-3].split(";")).astype(np.float)
                ir = np.array(result[-2].split(";")).astype(np.float)
                vis = np.array(result[-1].split(";")).astype(np.float)

        else:
            sysname = SystematicGasExName(result[1], creationtime)
            wavenumber = np.array(result[-4].split(";")).astype(np.float)
            sfg = np.array(result[-3].split(";")).astype(np.float)
            ir = np.array(result[-2].split(";")).astype(np.float)
            vis = np.array(result[-1].split(";")).astype(np.float)



        spec = SfgSpectrum(wavenumber, sfg, ir, vis, sysname)

        return spec

    def general_fetch(self, condition_1, condition_2, default_data=True):
        """A function to fetch spectra from the database according to specified match criteria.
        The default_data kwarg controlls the datatable that is used."""

        if default_data == True:
            database = "sfg_database"

        else:
            database ="sfg_gasex"

        command = "SELECT * from " + database + " WHERE "+condition_1+"="+condition_2
        self.cur.execute(command)
        keys = []
        for item in self.cur.fetchall():
            id = item[0]
            self.subset_ids.append(id)
            self.subset.append(self.fetch_single(id))

    def fetch_gasex_sfg(self):
        """Fetches all SfgSpectra from the GasEx cruise and puts them in the subset attribute"""
        command = "SELECT * from " + self.gasex_table
        self.cur.execute(command)
        for item in self.cur.fetchall():
            id = item[0]
            self.subset_ids.append(id)
            self.subset.append(self.fetch_single(id, default_data=False))

        self.gasex_included = True

    def general_refine(self, condition1, condition2, default_data=True):
        """Refinement of the actual subset by applying further match criteria. This is the
        actual implementation of the get method abstracting away the database access routine from the
        user."""
        temp = []
        temp_id = []

        for id in self.subset_ids:

            try:
                s = self.fetch_single(id, condition1, condition2, default_data)
                temp.append(s)
                temp_id.append(id)
            except IndexError:
                pass

        self.recover = self.subset
        self.recover_ids = self.subset_ids

        self.subset = temp
        self.subset_ids = temp_id

    def plot(self):
        """Runs an external PyQt application plotting all current spectra in the subset. The external
        app gives access to further manual analysis functionality"""
        run_app(self.subset, self.session_id)

    def show(self):
        """Prints a list of the spectra of the current subset including their subset index by which they
        can be accessed easily"""
        for i,spectrum in enumerate(self.subset):
            print(str(i)+" "+str(spectrum))

    def clear(self):
        """Deletes all spectra from the current subset"""
        self.subset = []
        self.subset_ids = []

    def get(self, flagstring, ref=False):
        """Fetches spectra according to the desired properties and adds them to the subset"""
        t = self.flagstring_split(flagstring)

        if t[0] == "su" or t[0] == "se":
            condition1 = self.retranslate_name(t[0])
            condition2 = "\""+self.retranslate_name(t[1])+"\""
            if ref is False:
                self.general_fetch(condition1, condition2)
                self.general_fetch(condition1, condition2, default_data=self.gasex_included)
            if ref is True:
                self.general_refine(condition1, condition2)
                self.general_refine(condition1, condition2, default_data=self.gasex_included)

        elif t[0] == "name":
            self.general_fetch(t[0], "\""+t[1]+".sfg"+"\"")

    def ref(self, flagstring):
        """Refines the current subset by applying further match criteria"""
        self.get(flagstring, ref=True)

    def flagstring_split(self, flagstring):
        """This function processes the given flagstring and returns a list of SFG objects
        which are passed through the Finder methods utilizing the flags and options"""
        f = flagstring.split(" ")
        return f

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
        # todo: include SFG stations here as well as in the match to station function
        # todo: include error handling (LtManager defined?)
        stations = []
        for isotherm in self.lt_manager.isotherms:

            temp = isotherm.long_station
            if temp not in stations:
                stations.append(Station(temp))
        self.stations = {s.station_hash: s for s in stations}

    def match_to_stations(self):
        """Matches the LtIsotherms in the LtManager to the generated station list. The
        stations attribute of the SessionControllManager is a dictionary with the station names
        as keys and a list with LtIsotherms matching the station as value"""

        for isotherm in self.lt_manager.isotherms:

            self.stations[isotherm.station_hash].lt_isotherms.append(isotherm)

        if self.gasex_included == True:

            for spectrum in self.subset:
                if isinstance(spectrum.name, SystematicGasExName):

                    try:
                        self.stations[spectrum.name.station_hash].sfg_spectra.append(spectrum)
                    except KeyError:
                        self.stations[spectrum.name.station_hash] = Station(spectrum.name.station)
                        self.stations[spectrum.name.station_hash].sfg_spectra.append(spectrum)


    def get_station_numbers(self):
        """Traverse the stations and assign them a number in chronological order (1-n)"""

        stations = [i for i in self.stations.values()]
        for i, s in enumerate(stations):
            print(str(i+1)+"   "+str(stations[i].name)+"\n")
        stations.sort()
        for i, s in enumerate(stations):
            s.station_number = i+1
            print(str(i+1)+"   "+str(stations[i].name)+"\n")


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

        self.day = None
        self.long_day = None
        self.type = None
        self.station = None
        self.long_station = None
        self.number = None
        self.speed = None

        self.partners = []

        self.process_name()

    def __str__(self):
        return self.name+" LtIsotherm Object"

    def __repr__(self):
        return self.name + " LtIsotherm Object"

    def __lt__(self, other):
        if self.get_maximum_pressure() < other.get_maximum_pressure():
            return True

    def get_day(self, string):

        if len(string) == 4:
            day = string[2:]
        else:
            raise ValueError("Invalid day string length at spectrum "+self.name)

        day = int(day)
        day = day-2
        return day

    def process_name(self):
        """Function to extract metainformation from the filename"""
        temp = self.name.split("_")
        self.day = self.get_day(temp[1])
        self.long_day = temp[1]
        self.type = temp[3].lower()
        self.station = temp[2]
        self.number = temp[4]
        self.long_station = self.long_day+"_"+self.station
        self.station_hash = self.long_day+"_"+self.station[1]

        try:
            self.speed = temp[5]
        except IndexError:
            print(self.name+ " has not a defined compression speed!")

    def drop_ascii(self):
        """Drops an ascii file with semikolon-separated data in the form time;area;surface pressure. Intention
        is easy interfacing with external software like Excel or Origin"""

        with open(self.name+".out", "w") as outfile:

            for a,b,c in zip(self.time, self.area, self.pressure):
                outfile.write(str(a)+";"+str(b)+";"+str(c)+"\n")

    def get_maximum_pressure(self, shrinked=None):
        """Returns the maximum measured surface pressure. Note: This property is uesd for the less-then
        operator implementation of this class!"""
        if shrinked == None:
            return np.max(self.pressure)
        else:
            try:
                return np.max(shrinked)
            except:
                #todo specify the type of error numpy will throw
                raise TypeError("Can not calc maximum for this operand")

    def derive_pressure(self):
        """Calculates the difference quotient of the surface pressure with respect to the area.
        Useful for calculation of surface elasticity"""
        return np.diff(self.pressure)/np.diff(self.area)

    def same_sample(self, other):
        """Checks wether two isotherms belong to the same sample. This is the case if the same sample
        is measured several times in a row. Returns a bool."""
        if self.create_sample_hash() == other.create_sample_hash():
            return True
        else:
            return False

    def create_sample_hash(self):
        """Creates a string which identifies the sample by the day, station and type. All isotherms
        from the same station taken by the same method (plate or screen or CTD) will therefore yield
        the same sample_hash. This is usually used to match samples together"""
        return str(self.day)+self.station+self.type+str(self.number)

    def create_pointlist(self, x_array):
        """Returns a list containing the index, the x_array (usually area or time) and the surface pressure.
        This is used for example by the GUI functions to find the closest datapoint to a mouse-defined position
        in the plot"""

        output = []
        for i, (a, b) in enumerate(zip(x_array, self.pressure)):
            output.append((a, b, i))
        return output

    def get_slice(self, x_array, lower, upper):
        """Returns a slice of the x_array (usually time or area) defined by the lower and upper integer
        index."""

        x_out = x_array[lower:upper+1]
        y_out = self.pressure[lower:upper+1]
        return x_out, y_out

    def calculate_elasticity(self):
        """Returns the surface elasticity of the isotherm"""

        xdata = self.area[::-1]
        ydata = self.pressure[::-1]
        out = []

        for i in range(len(self.pressure)-1):
            p = abs(ydata[i+1] - ydata[i])
            a = abs((xdata[i+1] - xdata[i]))
            A = (xdata[i+1] - xdata[i])/2

            out.append(p/a*A)

        return np.array(out[::-1])

    def smooth(self):
        """Performs a smooth operation of the measured pressure involving a Savitzky-Golay-filter"""
        return savgol_filter(self.pressure, 9, 5)


class LtManager:
    """A class to perform operations on a set of LtIsotherm objects. It is an extension to the SessionControllManager
    to extend his features with isotherm handling. It relies on sqlite databases as well."""

    def __init__(self, database, table="lt_gasex"):

        self.database = database
        self.cursor = database.cursor()
        self.table = table
        self.isotherms = []
        self.days = None
        self.ordered_days = {}

        self.get_all_isotherms()
        self.join_days_isotherms()
        self.order_by_sample()
        self.join_same_measurement()

    def get_all_isotherms(self):
        """Fetches all isotherm data from the database and transforms them to LtIsotherm objects."""
        command="SELECT * from "+self.table
        self.cursor.execute(command)
        result = self.cursor.fetchall()
        for i in result:
            lt = LtIsotherm(i[1], i[2], i[4], i[5], i[6], i[7])
            self.isotherms.append(lt)

    def get_days(self):
        """Extract the days on which the samples were taking from the isotherms."""
        days = []
        for isotherm in self.isotherms:
            if isotherm.day not in days:
                days.append(isotherm.day)
        return days

    def join_days_isotherms(self):
        """Matches the isotherms to the days of measurement, storing this mapping as a dictionary
        in the days attribute."""

        days = self.get_days()
        self.days = {i: [j for j in self.isotherms if j.day == i] for i in days}

    def order_by_sample(self):
        """Matches isotherms to a certain sample, ensuring that all consecutive measurements of
        the same sample are stored together."""

        for item in self.days:
            stations = []

            for isotherm in self.days[item]:
                if isotherm.station not in stations:
                    stations.append(isotherm.station)

            self.ordered_days[item] = {i: [j for j in self.days[item] if j.station == i] for i in stations}

        for day in self.ordered_days:
            for station in self.ordered_days[day]:

                types = []
                for isotherm in self.ordered_days[day][station]:
                    if isotherm.type not in types:
                        types.append(isotherm.type)
                self.ordered_days[day][station] = {i: [j for j in self.ordered_days[day][station] if j.type == i] for i in types}

    def join_same_measurement(self):
        """Stores information about the other measurements of a single sample in the isotherms partner
        attribute."""
        for i, isotherm in enumerate(self.isotherms):

            for isotherm2 in self.isotherms[i+1:]:
                if isotherm.same_sample(isotherm2) and isotherm2 not in isotherm.partners:
                    isotherm.partners.append(isotherm2)




class Station:
    """A class carrying all information and data for a given cruise station, especially SFG and isotherms"""

    def __init__(self, name):

        self.name = name
        self.station_hash = None
        self.date, self.id = self.name.split("_")
        self.cruise_day = int(self.date[2:])-2
        self.sfg_spectra = []
        self.lt_isotherms = []
        self.station_number = None
        self.lt_joined = {}
        self.stats = {
            "positive_plate": 0,
            "positive_screen": 0,
            "total_screen": 0,
            "total_plate": 0,
            "screen_av": 0,
            "plate_av": 0,
            "total" : 0,
            "std_plate": 0,
            "std_screen": 0,
            "std_total": 0,
            "percent_plate": 0,
            "percent_screen": 0,
            "total_percent" : 0,
            "total_av" : 0,
        }

        self.set_station_hash()

    def __lt__(self, other):
        """Determines which of two stations took place earlier."""

        daytag = int(self.name.split("_")[0][2:])
        monthtag = int(self.name.split("_")[0][0:2])
        stationtag = int(self.id[1])

        daytag2 = int(other.name.split("_")[0][2:])
        monthtag2 = int(other.name.split("_")[0][0:2])
        stationtag2 = int(other.id[1])

        outbool = None

        if monthtag <= monthtag2:

            if daytag == daytag2:
                if stationtag < stationtag2:
                    outbool = True

            elif daytag < daytag2:
                outbool = True

            else:
                outbool = False
        else:
            outbool = False

        return  outbool

    def __rpr__(self):
        return f'Station {self.id[1]} on date {self.date}'

    def __str__(self):
        return f'Station {self.id[1]} on date {self.date}'

    def set_station_hash(self):
        temp = self.name.split("_")
        try:
            self.station_hash = temp[0]+"_"+temp[1][1]
        except IndexError:
            print("*"*80+"\n"+str(temp))
            print(self.name)


    def count_per_type(self):
        """Counts the occurence of plate and screen samples as well as how many of those are positive
        with respect to surfactants (high surface pressure in Langmuir trough measurement). Stores all the
        statistical information in the stats dictionary."""

        # todo: handle CTD samples!
        # todo: implement a similar sample for SFG

        isos = [i[0] for i in self.lt_joined.values()]
        plate_array = []
        screen_array = []

        for isotherm in isos: # type: LtIsotherm

            p_max = isotherm.get_maximum_pressure()

            if p_max < 73:
                if isotherm.type == "p":

                    if p_max > 2:
                        self.stats["positive_plate"] += 1

                    self.stats["total_plate"] += 1
                    self.stats["plate_av"] += p_max
                    plate_array.append(p_max)

                elif isotherm.type[0] == "s":

                    if 72 > p_max > 2:
                        self.stats["positive_screen"] += 1

                    self.stats["total_screen"] += 1
                    self.stats["screen_av"] += p_max
                    screen_array.append(p_max)

        self.stats["total"] = (self.stats["total_screen"] + self.stats["total_plate"])

        try:
            self.stats["total_av"] = (self.stats["screen_av"]+self.stats["plate_av"])/self.stats["total"]
        except ZeroDivisionError:
            pass
        try:
            self.stats["screen_av"] /= self.stats["total_screen"]
        except ZeroDivisionError:
            pass
        try:
            self.stats["plate_av"] /= self.stats["total_plate"]
        except ZeroDivisionError:
            pass
        try:
            self.stats["percent_screen"] = self.stats["positive_screen"]/self.stats["total_screen"]*100
        except ZeroDivisionError:
            pass
        try:
            self.stats["percent_plate"] = self.stats["positive_plate"] / self.stats["total_plate"] * 100
        except ZeroDivisionError:
            pass
        try:
            self.stats["total_percent"] = (self.stats["positive_screen"]+self.stats["positive_plate"])/self.stats["total"]*100
        except ZeroDivisionError:
            pass
        total_array = plate_array+screen_array

        self.stats["std_plate"] = np.std(plate_array)
        self.stats["std_screen"] = np.std(screen_array)
        self.stats["std_total"] = np.std(total_array)


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

    def print_stats(self):
        """Formatted output of the stations stats, calculated from the LtIsotherms belonging to the
        station."""

        for item in self.stats:
            s = f'{item} : {self.stats[item]}\n'
            print(s)

#

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
        plt.scatter(isotherm.day, isotherm.get_maximum_pressure(), color=color, s=40, marker=marker )
        plt.text(isotherm.day+0.1, isotherm.get_maximum_pressure()+0.1, isotherm.station+isotherm.type+isotherm.number,
                 fontsize=7)
        #plt.text(isotherm.day+0.2, isotherm.get_maximum_pressure()+0.1, isotherm.name, fontsize=7)

    plt.xlabel("days of cruise")
    plt.ylabel("maximum pressure/ mNm$^{-1}$")
    plt.grid()
    #legend_elements = [Line2D([0], [0], marker='o', color='r', label='glass plate', markerfacecolor='g', markersize=40)]
    #plt.legend(handles=legend_elements)
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
            plt.title(str(isotherm.day)+" "+isotherm.station+" "+isotherm.type+" "+str(isotherm.number))
            plt.grid()
            plt.legend()
            plt.savefig(isotherm.create_sample_hash()+".png")
            plt.cla()

        else:
            print("Already processed!")


def sfg_pub_plot(speclist, title="default", normalized="false"):
    """Produces a pre-formatted SFG plot from a list of SFG spectrum objects"""

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.set_xlabel("Wavenumber/ cm$^{-1}$", fontsize=10)
    ax.set_ylabel("Norm. SFG intensity/ arb. u.", fontsize=10)

    inc = 0.25/len(speclist)
    counter = 0
    for spectrum in speclist:
        eff_alpha = 0.75+inc*counter
        if normalized == "false":
            ax.plot(spectrum.wavenumbers, spectrum.normalized_intensity, linewidth=1.5, marker="o", markersize=3,
                    alpha=eff_alpha, label=spectrum.name.full_name)
            ax.plot(test, baseline(test))
            spectrum.correct_baseline()
            ax.plot(spectrum.wavenumbers, spectrum.baseline_corrected)
        elif normalized == "true":
            ax.plot(spectrum.wavenumbers, spectrum.normalize_to_highest(), linewidth=1.5, marker="o", markersize=3,
                    alpha=eff_alpha, label=spectrum.name.full_name)
        counter +=1

    size = fig.get_size_inches()
    ratio = size[0]/size[1]
    fig.set_size_inches(3.2*ratio, 3.2)
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
        ax.plot(spectrum.wavenumbers, spectrum.normalize_to_highest()+offset, linewidth=1.5, marker="o", markersize=3,
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
    ax2 =fig.add_subplot(212)
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
        ax1.plot(spectrum.wavenumbers, spectrum.normalize_to_highest()+offset, linewidth=1.5, marker="o", markersize=3,
                alpha=eff_alpha, label=spectrum.name.full_name)
        counter += 1
        offset += 0.2

    inc = 0.25 / len(speclist1)
    counter = 0
    offset = 0
    for spectrum in speclist2:
        eff_alpha = 0.75 + inc * counter
        ax2.plot(spectrum.wavenumbers, spectrum.normalize_to_highest()+offset, linewidth=1.5, marker="o", markersize=3,
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
    S.fetch_gasex_sfg()
    S.collect_stations()
    S.match_to_stations()
    S.get_station_numbers()

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

        for spectrum in s.sfg_spectra: # type: SfgSpectrum

            station = s.station_number
            b = spectrum.slice_by_borders(3000,2800)
            max = np.max(spectrum.normalized_intensity[b[0]:b[1]+1])
            max_ins.append(max)

        max_ins = np.array(max_ins)
        av = np.average(max_ins)
        std = np.std(max_ins)
        u_ax.errorbar(station, av, std, alpha=0.9,fmt="o", color="b", label="average",barsabove="true", capsize=5, capthick=2)

    ax1.errorbar(stations, totals, yerr=t_std, fmt="o",color="r",barsabove="true", capsize=5, capthick=2)
    rects2 = ax2.bar([a[0] + 0.2 for a in t_percentages], [a[1] for a in t_percentages], alpha=0.45, label="total")

    ax2.legend()
    #u_ax.legend()
    plt.tight_layout()
    plt.show()


def sfg_with_lt(S):
    """Takes Session controll manager as argument. Performs the bar/max surface pressure plot plot for the LT
        isotherms"""
    S.fetch_gasex_sfg()
    S.collect_stations()
    S.match_to_stations()
    S.get_station_numbers()

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

        for spectrum in s.sfg_spectra: # type: SfgSpectrum

            b = spectrum.slice_by_borders(3000,2800)
            max = np.max(spectrum.normalized_intensity[b[0]:b[1]+1])
            mapping = {"p":["g","plate"],
                       "s1":["b","screen"],
                       "s2":["b","screen"],
                       "s3":["b","screen"],
                       "deep":"black",
                       "low":"black"}
            try:
                ax1.scatter(station, max,alpha=0.9, color=mapping[spectrum.name.type][0], label=mapping[spectrum.name.type][1],
                            s=60)
            except KeyError:
                ax1.scatter(station, max, alpha=0.9, color="black", label="CTD", s=60)



        ax3.scatter(s.cruise_day, max, s=0)
        ax1.set_ylim(0, 0.0008)
        #ax1.legend()

        s_percentages.append([station, s.stats["percent_screen"]])
        p_percentages.append([station, s.stats["percent_plate"]])
        t_percentages.append([station, s.stats["total_percent"]])

    rects1 = ax2.bar([a[0] - 0.2 for a in p_percentages], [a[1] for a in p_percentages], alpha=0.35, width=0.2,
                     color="g", label="plate, trough (positive)")
    rects2 = ax2.bar([a[0] for a in s_percentages], [a[1] for a in s_percentages], alpha=0.35, width=0.2, color="b",
                     label="screen, trough (positive)")
    rects2 = ax2.bar([a[0] + 0.2 for a in t_percentages], [a[1] for a in t_percentages], alpha=0.35, width=0.2,
                     color="r", label="total, trough (positive)")

    #ax2.legend()

    h1 = [Line2D([0],[0], marker='o', color='w', markerfacecolor='g',markersize=10),
          Line2D([0], [0], marker='o', color='w', markerfacecolor='b',markersize=10),
          Line2D([0], [0], marker='o', color='w', markerfacecolor='black',markersize=10)]

    l1 = ["plate, SFG", "screen, SFG", "CTD, SFG"]
    h2, l2 = ax2.get_legend_handles_labels()
    ax1.legend(h1 + h2, l1 + l2, loc="upper center", ncol=2).draggable()
    plt.sca(ax1)
    plt.tight_layout()
    plt.show()



S = SessionControlManager("sfg.db", "test")
S.set_lt_manager()
S.fetch_gasex_sfg()
S.collect_stations()
S.match_to_stations()
S.get_station_numbers()

specs = [i for i in S.subset if isinstance(i.name, SystematicGasExName)]
f=sfg_pub_plot(specs[3:6])
plt.show()



# S.get("name 20180319_PA_5_x1_#1_5mM")
# pa = S.subset[0]
# S.clear()
#
# S.get("name 20180315_SA_BX12_10_x1_#2_mixedlayer1zu2mM")
# sa_bx = S.subset[0]
# S.clear()
#
# S.get("name 20180320_PA_BX12_10_x1_#1_1.5zu1.5mM")
# pa_bx = S.subset[0]
# S.clear()
#
# S.get("name 20180524_BX12_8_x2_#2_1mM")
# bx = S.subset[0]
# S.clear()
#
# S.get("name 20160316_SA_5_x1_#1_5mM")
# sa = S.subset[0]
# S.clear()
#
# stack1 = [sa, sa_bx, bx]
# stack2 = [pa, pa_bx, bx]
#
# f = sfg_doublestack_plot(stack1, stack2)
# plt.show()














