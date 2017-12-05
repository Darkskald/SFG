# import functions as fun
import os
import shutil
import csv
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
from scipy.signal import savgol_filter
from scipy.integrate import simps as sp

rcParams['mathtext.default'] = 'regular'


# noinspection PySimplifyBooleanCheck,PySimplifyBooleanCheck
class Importer:
    # first: make a list of day folders in the archive directory
    def __init__(self):
        self.new_folders = []
        for folder in os.listdir("archive"):
            self.new_folders.append(folder)

        # now: get the files of each folder, strip FL off, add daytag to filename
        self.new_files = []
        for folder in self.new_folders:
            daytag = (folder.split("FL"))[0].strip()
            # now creating the day_information_file
            day_info_file = daytag + ".dif"
            if self.gate_keeper("library/day_information", day_info_file) == False:
                with open("library/day_information/" + day_info_file, "w") as outfile:
                    outfile.write(daytag + "\n")
                    files = []
                    for file in os.listdir("archive/" + folder):
                        files.append(file)
                    outfile.write(str(len(files)) + "\n")
                    for file in files:
                        outfile.write(file + "\n")
                    outfile.write("#")

            # next: collect all spectra per day and copy them to the library
            for file in os.listdir("archive/" + folder):
                new_filename = daytag + "_" + file
                # the gate_keeper function prohibits double import
                if self.gate_keeper("library", new_filename) == False:
                    shutil.copy2("archive/" + folder + "/" + file, "library/" + new_filename)

    def gate_keeper(self, checkdir, filename):
        """checks if filename is already present in the directory chekdir"""
        checklist = []
        for file in os.listdir(checkdir):
            checklist.append(file)
        if filename in checklist:
            print("File " + filename + " already exists")
            return True
        else:
            return False


class LibraryManager:
    """controll and maintain the library management file"""
    def __init__(self):
        # self.entries is a list of lists extracted from the lines of the library file and splitted at the ;
        self.entries = []
        with open("library_management.txt", "r") as infile:
            for line in infile:
                buffer = line.split(";")
                self.entries.append(buffer)

    def update(self):

        files = []
        for file in os.listdir("library"):
            if file.endswith(".sfg"):
                files.append(file)

        newfiles = []
        for file in files:
            for checkfile in self.entries:
                if file == checkfile[1].strip():
                    print("File already in library!")
            else:
                newfiles.append(file)

        with open("library_management.txt", "a") as outfile:
            counter = len(self.entries)
            for file in newfiles:
                counter += 1
                sfg = FileFetcher(file).sfg
                specrange = sfg.yield_spectral_range()
                outfile.write(str(counter) + ";" + file + ";" + str(specrange[0]) + ";" + str(specrange[1]) + ";" + str(
                    specrange[2]) + "\n")


class SfgSpectrum:
    # magic methods

    def __init__(self, wavenumbers, intensity, ir_intensity, vis_intensity, SystematicName):
        self.wavenumbers = wavenumbers
        self.raw_intensity = intensity
        self.vis_intensity = vis_intensity
        self.ir_intensity = ir_intensity
        self.name = SystematicName
        self.normalized_intensity = self.raw_intensity / (self.vis_intensity * self.ir_intensity)

    def __repr__(self):
        return self.name.full_name[:-4] + "              " + str(self.yield_spectral_range())

    def __str__(self):
        """Printable represantation of the SFG object"""
        date = self.name.date
        surf = self.name.surfactant
        srange = str(self.yield_spectral_range())
        full = self.name.full_name[:-4]

        return date + "\t" + surf + "\t" + srange + "\t\t" + full

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



    def normalize_to_highest(self, intensity="default"):
        """normalize an given array to its maximum, typically the normalized or raw intensity"""
        if intensity == "default":
            intensity = self.normalized_intensity
        norm_factor = np.max(intensity)
        return (intensity / norm_factor)

    def yield_peaklist(self, intensity="default", num=6, threshold=25):
        """yield a defined list of peaks separated from each other by minimum the threshold value in wavenumber"""
        pair_get = []
        out = []
        num = num
        if intensity == "default":
            intensity = self.normalized_intensity

        for i in range(len(self.wavenumbers)):
            pair_get.append([intensity[i], self.wavenumbers[i]])

        while len(out) < (num):

            if len(out) != 0:

                for i in range(len(out)):
                    k = max(pair_get)

                    if np.abs(out[i][1] - k[1]) < threshold:
                        pair_get.remove(k)
                        break
                    if i == (len(out) - 1):
                        out.append(k)
                        pair_get.remove(k)


            else:
                k = max(pair_get)
                out.append(k)

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
        self.normalized_intensity = y

    def detailed_analysis(self, threshold=3):
        """Function returning peak information (central wavenumber, flanks, integral). The threshold value characterizes
        the factor that a peak must be greater than the average intensity"""
        x_array = self.wavenumbers[::-1]
        y_array = self.normalized_intensity[::-1]

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

            # check for left border
            while slopes[k] > 0:
                k -= 1
            left = k

            # check for right border
            k = i + 1

            while ((slopes[k] < 0) and (k + 1 <= len(slopes) - 1)):
                k += 1

            right = k

            peak_tuples.append((center, left, right))

        data_out = []
        for i in peak_tuples:
            indices = (i[0], i[1], i[2])
            center = x_array[i[0]]
            left = x_array[i[1]]
            right = x_array[i[2]]
            peak_slice_x = x_array[i[1]:i[2] + 1]
            peak_slice_y = y_array[i[1]:i[2] + 1]
            area = self.integrate_peak(peak_slice_x, peak_slice_y)
            datapoints = len(peak_slice_x)

            data_out.append((center, left, right, peak_slice_x, peak_slice_y, datapoints, area, indices))
        return data_out

    def integrate_peak(self, x_array, y_array):
        """Numpy integration routine for numerical peak integration with the trapezoidal rule"""
        area = sp(y_array, x_array)
        return area


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

    def __init__(self, namestring):

        # load the allowed surfactans and sensitizers from files
        self.Surfactants = {}
        self.Sensitizers = {}

        with open("../name_info/Surfactants.txt", "r") as infile:
            for line in infile:
                collect = line.split(":")
                self.Surfactants[collect[0]] = collect[1].strip()

        with open("../name_info/Sensitizers.txt", "r") as infile:
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
        self.sensitizer_spread_volume = "-"
        self.photolysis = "none"
        self.sample_number = 1
        self.measurement = "#1"  # das ist nicht optimal, sample_number als int und measurement als string zu haben
        self.comment = "none"
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
                    sensitizer_spread_volume = i
            # pH value info will be comment
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
            self.sensitizer_spread_volume = "-"

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

    def __init__(self,names,sensitizers,surfactants):

        self.full_name = ("_").join(names)
        self.sensitizer = [s for s in sensitizers]
        self.surfactant = [s for s in surfactants]



class FileFetcher:
    """The file fetcher class interconnects different subdirectories of the folder and simplifies
    the handling of filenames and filepaths. It changes the working directory to the file storage directory,
    usually the library, creates a DataCollector obect and returns a SFG object"""

    def __init__(self, filename, destination="library"):
        self.filename = filename
        initial_wd = os.getcwd()
        os.chdir(destination)

        self.collector = DataCollector(filename)
        self.sfg = self.collector.yield_SfgSpectrum()
        os.chdir(initial_wd)


class DataCollector:
    # noinspection PyTypeChecker
    def __init__(self, filename):
        # the filename is ONLY the filename, not containing any directory information.
        self.file = filename

        # now start extraction
        data_collect = []

        with open(self.file, "r") as infile:

            readCSV = csv.reader(infile, delimiter="\t")
            for row in readCSV:
                data_collect.append(row)

            data_package = [0] * len(data_collect[0])
            for i in range(len(data_collect[0])):
                # mind this complex list comprehension
                data_package[i] = [j[i] for j in data_collect]
                # remove useles column
            del data_package[2]

            # convert strings to float and generate numpy array
            convert = []
            for i in data_package:
                q = [float(j) for j in i]
                convert.append(q)
            convert = [np.array(i) for i in convert]
            self.data_package = convert

    def yield_SfgSpectrum(self):
        #todo this function will need exception handling. It returns an SFG spectrum object
        sysname = SystematicName(self.file)
        data = self.data_package
        sfg = SfgSpectrum(data[0], data[1], data[3], data[2], sysname)
        return sfg


# noinspection PySimplifyBooleanCheck
class Finder:
    """Powerfull class generating sfg objects of ALL files in library and traversing them for match
    criteria. Contains a list of SFG objects whitch mach the criteria. Later features to extract 
    information from libary management file and day information file will be added"""

    def __init__(self):
        self.database = []
        for file in os.listdir("library"):
            if file.endswith(".sfg"):
                sfg = FileFetcher(file).sfg
                self.database.append(sfg)

    """Basic search functions. Choose spectra from database by date, surfactant etc.."""
    def date_based(self, dateflag, subset="default"):
        if subset == "default":
            subset = self.database
        matches = []
        dateflag = dateflag.strip("d")
        for spectrum in subset:
            if spectrum.name.date == dateflag:
                matches.append(spectrum)
        return matches

    def sample_based(self, sampleflag, subset="default"):
        if subset == "default":
            subset = self.database
        matches = []
        for spectrum in subset:
            if spectrum.name.sample_number == sampleflag or str(spectrum.name.sample_number) == str(sampleflag):
                matches.append(spectrum)
        return matches

    def surfactant_based(self, surfactantflag, subset="default"):
        if subset == "default":
            subset = self.database
        matches = []
        for spectrum in subset:
            if spectrum.name.surfactant == surfactantflag:
                matches.append(spectrum)
        return matches

    def sensitizer_based(self, sensitizerflag, subset="default"):
        if subset == "default":
            subset = self.database
        matches = []
        for spectrum in subset:
            if spectrum.name.surfactant == sensitizerflag:
                matches.append(spectrum)
        return matches

    def photo_based(self, subset="default", photolyzed=True):
        if photolyzed == True:
            if subset == "default":
                subset = self.database
            matches = []
            for spectrum in subset:
                if spectrum.name.photolysis != "none":
                    matches.append(spectrum)
        else:
            if subset == "default":
                subset = self.database
            matches = []
            for spectrum in subset:
                if spectrum.name.photolysis == "none":
                    matches.append(spectrum)
        return matches

    def measurement_based(self, measurementflag, subset="default"):
        if subset == "default":
            subset = self.database
        measures = measurementflag.split(",")
        matches = [i for i in subset if i.name.measurement in measures]
        return matches

    def comment_based(self, option="BoknisEckSample", subset="default"):
        if subset == "default":
            subset = self.database
        matches = []
        for spectrum in subset:
            if spectrum.name.comment == option:
                matches.append(spectrum)
        return matches

    def sur_volume_based(self, amount, subset="default"):
        if subset == "default":
            subset = self.database
        matches = [i for i in subset if i.name.surfactant_spread_volume == amount]
        return matches

    def unclassified(self):
        subset = self.database
        matches = [i for i in subset if i.name.surfactant == "unknown"]
        return matches

    def by_year(self, yearlist, subset="default"):

        if subset == "default":
            subset = self.database

        matches = [i for i in subset if str(i.name.date_split()[0]) in yearlist]
        return matches

    def by_month(self, month, subset="default"):
        if subset == "default":
            subset = self.database
        matches = [i for i in subset if i.name.date_split()[1] == month]
        return matches

    def by_monthrange(self, begin, end, subset="default"):
        if subset == "default":
            subset = self.database
        matches = [i for i in subset if (begin <= i.name.date_split()[1] <= end)]
        return matches


# noinspection PySimplifyBooleanCheck,PySimplifyBooleanCheck
class Plotter:
    """Simple class providing an interface between the matplotlib API and the SFG objects. Depending on the specific method
    and kwargs, a variety of different methods is accessible"""
    def __init__(self, speclist, raw=False, title="Default"):
        self.speclist = speclist
        self.raw = raw
        self.title = title

    def simple_plot(self):

        fig = plt.figure()
        ax = plt.subplot(111)
        for spectrum in self.speclist:
            wl = spectrum.wavenumbers
            if self.raw == False:
                intensity = spectrum.normalized_intensity
            else:
                intensity = spectrum.raw_intensity
            ax.plot(wl, intensity, label=spectrum.name.full_name, linestyle='--', markersize=4, marker="o")

        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
        ax.grid()
        ax.set_title(self.title)
        ax.set_xlabel("Wavenumber/ $cm^{-1}$")
        ax.set_ylabel("Intensity/ a.u.")
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        plt.show()

    def raw_plot(self):
        for spectrum in self.speclist:
            wl = spectrum.wavenumbers
            intensity = spectrum.raw_intensity
            plt.plot(wl, intensity, label=spectrum.name.full_name, linestyle='--', markersize=4, marker="o")

        plt.grid(True)
        plt.title(self.title)
        plt.xlabel("Wavenumber/ $cm^{-1}$")
        plt.ylabel("Raw intensity/ a.u.")
        if len(self.speclist < 6):
            plt.legend(loc="upper right")
        else:
            plt.legend(bbox_to_anchor(1, 0.5), loc='center left', numpoints=1)
        plt.show()

    def raw_plot_plus_ir(self):
        for spectrum in self.speclist:
            wl = spectrum.wavenumbers
            intensity = spectrum.raw_intensity
            ir = spectrum.ir_intensity
            plt.plot(wl, intensity, label=spectrum.name.full_name, linestyle='dotted', markersize=4, marker="o")
            plt.plot(wl, ir, label=spectrum.name.full_name + " IR_intensity", linestyle='--')
        plt.grid(True)
        plt.title(self.title)
        plt.xlabel("Wavenumber/ $cm^{-1}$")
        plt.ylabel("Raw intensity/ a.u.")
        plt.legend(loc="upper right")
        plt.show()

    def custom_plot(self):

        fig = plt.figure()
        ax = plt.subplot(111)
        for spectrum in self.speclist:
            wl = spectrum.wavenumbers
            if self.raw == False:
                intensity = spectrum.normalized_intensity
            else:
                intensity = spectrum.raw_intensity
            photo = ""

            if spectrum.name.photolysis != "none":
                photo = spectrum.name.photolysis

            stitle = str(spectrum.name.sample_number) + " " + spectrum.name.comment + " " + str(
                spectrum.name.measurement) + " " + photo
            ax.plot(wl, intensity, label=stitle, linestyle='dotted', markersize=2, marker="o")

        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
        ax.grid()
        ax.set_title(self.title)
        ax.set_xlabel("Wavenumber/ $cm^{-1}$")
        ax.set_ylabel("Intensity/ a.u.")
        lgd = ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        plt.savefig("analysis_out/" + self.title + ".pdf", bbox_extra_artists=(lgd,), bbox_inches='tight')
        plt.close()
        with open("analysis_out/" + self.title + ".log", "w") as outfile:

            outfile.write("Logfile of Spectral analysis routine")
            outfile.write("List of included spectra\n")
            for spectrum in self.speclist:
                outfile.write(spectrum.name.full_name + "\n")


class Analyzer:
    """This class takes, what a surprise, a list of SFG spectra as constructor argument. Its purpose
    is to perform analytical tasks with the spectral data, e.g. compare peaks, integral, datapoints
    and will be extendet to handle statistics in the future"""

    def __init__(self, speclist):
        self.speclist = speclist


class Planer:
    """A class to perform task management and plan the research tasks for a period of measurements"""

    def __init__(self):

        self.tasks = []
        self.count = 0
        self.refresh()
        print("Planer initialized. The current task count is " + str(self.count))
        self.show_tasks()

    def show_tasks(self):

        index = 1
        for i in self.tasks:
            print(str(index) + " : " + i)
            index += 1

    def add_task(self, string):

        with open("tasks", "a") as infile:
            infile.write(string + "\n")
        self.refresh()

    def refresh(self):
        self.tasks = []
        with open("tasks", "r") as infile:
            for line in infile:
                self.tasks.append(line)
                self.count = len(self.tasks)

    def done(self, number):

        print(str(self.tasks[number - 1]) + "was removed!")
        del self.tasks[number - 1]
        self.update_file()

    def update_file(self):

        with open("tasks", "w") as outfile:
            for task in self.tasks:
                outfile.write(task)
