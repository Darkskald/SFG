import csv

import numpy as np


#basic classes providing functionality to deal with SFG spectra
from scipy.signal import savgol_filter


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


#auxiliary name classes carrying metainformation about the spectra
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


class AddedName(SystematicName):
    """This class is derived from the SfgSpectrum and represents the result of the addition of Sfg intensities."""

    def __init__(self, names):
        self.full_name = ("_").join(names)