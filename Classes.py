# import functions as fun
import os
import shutil
import csv
import time
import datetime
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
import sqlite3
from scipy.signal import savgol_filter
from scipy.integrate import simps as sp

#for session controll manager
#from new_gui import run_app

rcParams['mathtext.default'] = 'regular'


# noinspection PySimplifyBooleanCheck,PySimplifyBooleanCheck


class SfgSpectrum:
    # magic methods

    def __init__(self, wavenumbers, intensity, ir_intensity, vis_intensity, SystematicName):
        self.wavenumbers = wavenumbers
        self.raw_intensity = intensity
        self.vis_intensity = vis_intensity
        self.ir_intensity = ir_intensity
        self.name = SystematicName
        self.normalized_intensity = self.raw_intensity / (self.vis_intensity * self.ir_intensity)
        self.peaks = self.detailed_analysis()

    def __repr__(self):
        return self.name.full_name[:-4] + "              " + str(self.yield_spectral_range())

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
            sensitizers = [self.name.sensitizer, SFG2.name.sensitizer]
            surfactants = [self.name.surfactant, SFG2.name.surfactant]

            name = AddedName(names, sensitizers, surfactants)
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
        return np.sort(self.normalized_intensity)


# noinspection PyMissingConstructor
class AddedSpectrum(SfgSpectrum):
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
            sensitizers = [self.name.sensitizer, SFG2.name.sensitizer]
            surfactants = [self.name.surfactant, SFG2.name.surfactant]

            name = AddedName(names, sensitizers, surfactants)
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


# noinspection PySimplifyBooleanCheck,PySimplifyBooleanCheck,PySimplifyBooleanCheck,PySimplifyBooleanCheck
class SystematicName:

    def __init__(self, namestring, creation_time="unknown"):

        self.refpath = "name_info/"

        try:
            with open(self.refpath+"Surfactants.txt") as outfile:
                pass
        except FileNotFoundError:
            print("change dir")
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


class AddedName(SystematicName):
    """This class is derived from the SfgSpectrum and represents the result of the addition of Sfg intensities."""
    def __init__(self,names,sensitizers,surfactants):

        self.full_name = ("_").join(names)
        self.sensitizer = [s for s in sensitizers]
        self.surfactant = [s for s in surfactants]


class Analyzer:
    """This class takes, what a surprise, a list of SFG spectra as constructor argument. Its purpose
    is to perform analytical tasks with the spectral data, e.g. compare peaks, integral, datapoints
    and will be extendet to handle statistics in the future"""

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

    def __init__(self, database, id):

        self.db = sqlite3.connect(database, detect_types=sqlite3.PARSE_DECLTYPES|sqlite3.PARSE_COLNAMES)
        self.cur = self.db.cursor()
        self.table = "sfg_database"
        self.session_id = id
        self.Surfactants = {}
        self.Sensitizers = {}
        self.Makros = {}
        self.get_senssurf_names()

        self.lt_manager = None



        #former IpyInterpreter functionality, tracking the primary key in parallel
        self.subset_ids = []
        self.subset = []

        self.recover_ids = []
        self.recover = []

    def get_senssurf_names(self):

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

        if stri in self.Surfactants:
            return self.Surfactants[stri]
        elif stri in self.Sensitizers:
            return self.Sensitizers[stri]
        elif stri in self.Makros:
            return self.Makros[stri]
        else:
            print("Retranslation failed. Unknown expression.")

    def fetch_single(self, number, condition3=None, condition4=None):

        command = "SELECT * FROM "+self.table+" WHERE ROWID="+str(number)

        if condition3 is not None:
            command += " AND "+str(condition3)+"="+str(condition4)

        self.cur.execute(command)

        result = self.cur.fetchall()[0]

        creationtime = result[2]

        sysname = SystematicName(result[1], creationtime)
        wavenumber = np.array(result[-7].split(";")).astype(np.float)
        sfg = np.array(result[-6].split(";")).astype(np.float)
        ir = np.array(result[-5].split(";")).astype(np.float)
        vis = np.array(result[-4].split(";")).astype(np.float)

        spec = SfgSpectrum(wavenumber, sfg, ir, vis, sysname)

        return spec

    def clean(self):

        self.cur.close()
        self.db.close()

    def general_fetch(self, condition_1, condition_2, database="sfg_database"):
        command = "SELECT * from " + database + " WHERE "+condition_1+"="+condition_2
        self.cur.execute(command)
        keys = []
        for item in self.cur.fetchall():
            id = item[0]
            self.subset_ids.append(id)
            self.subset.append(self.fetch_single(id))

    def general_refine(self, condition1, condition2):
        temp = []
        temp_id = []

        for id in self.subset_ids:

            try:
                s = self.fetch_single(id, condition1, condition2)
                temp.append(s)
                temp_id.append(id)
            except IndexError:
                pass

        self.recover = self.subset
        self.recover_ids = self.subset_ids

        self.subset = temp
        self.subset_ids = temp_id

    def recovery(self):
        self.subset = self.recover
        self.subset_ids = self.recover_ids

    def plot(self):
        run_app(self.subset, self.session_id)

    def show(self):
        for i,spectrum in enumerate(self.subset):
            print(str(i)+" "+str(spectrum))

    def clear(self):
        self.subset = []
        self.subset_ids = []

    def get(self, flagstring, ref=False):
        t = self.flagstring_split(flagstring)
        print(t)
        if t[0] == "su" or t[0] == "se":
            condition1 = self.retranslate_name(t[0])
            condition2 = "\""+self.retranslate_name(t[1])+"\""
            if ref is False:
                self.general_fetch(condition1, condition2)
            if ref is True:
                self.general_refine(condition1, condition2)

    def ref(self, flagstring):
        self.get(flagstring, ref=True)

    def flagstring_split(self, flagstring):
        """This function processes the given flagstring and returns a list of SFG objects
        which are passed through the Finder methods utilizing the flags and options"""
        f = flagstring.split(" ")
        flag = f[0]
        options = f[1]
        return (flag, options)

    def rec(self):
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

        self.lt_manager = LtManager(self.db)


class LtIsotherm:

    def __init__(self, *args):

        self.name = args[0]
        self.measured_time = args[1]
        self.time = args[2]
        self.area = np.array([float(i) for i in args[3].strip("[]").split(",")])
        self.apm = np.array(args[4].strip("[]"))
        self.pressure = np.array(args[5].strip("[]"))

    def __str__(self):
        return self.name+" LtIsotherm Object"

    def drop_ascii(self):

        with open(self.name+".out", "w") as outfile:

            for a,b in zip(self.area, self.pressure):
                outfile.write(str(a)+";"+str(b))

    def get_maximum_pressure(self, shrinked=None):

        if shrinked == None:
            return np.max(self.pressure)
        else:
            try:
                return np.max(shrinked)
            except:
                #todo specify the type of error numpy will throw
                raise TypeError("Can not calc maximum for this operand")

class LtManager:

    def __init__(self, database, table="lt_gasex"):

        self.database = database
        self.cursor = database.cursor()
        self.table = table
        self.isotherms = []

    def get_all_isotherms(self):

        command="SELECT * from "+self.table
        self.cursor.execute(command)
        result = self.cursor.fetchall()
        for i in result:
            lt = LtIsotherm(i[1], i[2], i[4], i[5], i[6], i[7])
            self.isotherms.append(lt)


S = SessionControlManager("sfg.db", "test")
S.set_lt_manager()
S.lt_manager.get_all_isotherms()
for i in S.lt_manager.isotherms: #type: LtIsotherm
    print(i)
    #print(type(i))
    print(i.area)
    print(type(i.area))








