import sqlite3
import os
import shutil
import datetime
import numpy as np
import csv


def debug(func):

    def debug_wrapper(*args, **kwargs):
        print(func.__name__ +" calling! Entering function scope")
        func(args, kwargs)
        print(func.__name__ + " calling! ... leaving function scope")
    return debug_wrapper


class Importer:
    # first: make a list of day folders in the archive directory
    def __init__(self, input_folder="archive", output_folder="library"):

        self.new_folders = []
        self.speclist = []

        for folder in os.listdir(input_folder):
            self.new_folders.append(folder)

        # now: get the files of each folder, strip FL off, add daytag to filename
        self.new_files = []

        for folder in self.new_folders:

            if "FL" in folder:
                daytag = (folder.split("FL"))[0].strip()
            elif "TT" in daytag:
                daytag = (folder.split("TT"))[0].strip()

            # now creating the day_information_file
            day_info_file = daytag + ".dif"
            if self.gate_keeper(output_folder + "/" + "day_information/", day_info_file) == False:

                with open(output_folder + "/" + "/day_information/" + day_info_file, "w") as outfile:
                    outfile.write(daytag + "\n")
                    files = []
                    for file in os.listdir(input_folder + "/" + folder):
                        files.append(file)
                    outfile.write(str(len(files)) + "\n")
                    for file in files:
                        outfile.write(file + "\n")
                    outfile.write("#")

            # next: collect all spectra per day and copy them to the library
            for file in os.listdir(input_folder + "/" + folder):
                new_filename = daytag + "_" + file
                # the gate_keeper function prohibits double import
                if self.gate_keeper(output_folder, new_filename) == False:
                    shutil.copy2(input_folder + "/" + folder + "/" + file, output_folder + "/" + new_filename)

        self.create_speclist(output_folder)

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

    def get_spectrum(self, filename, destination="library"):

        initial_wd = os.getcwd()
        os.chdir(destination)

        collector = DataCollector(filename)
        sfg = collector.yield_SfgSpectrum()
        sfg.creation_time = collector.creation_time
        os.chdir(initial_wd)
        return sfg

    def create_speclist(self, target):

        for file in os.listdir(target):
            if file.endswith(".sfg"):
                sfg = self.get_spectrum(file, destination=target)
                self.speclist.append(sfg)

    def generate_sql(self, flag):

        self.wizard = SqlWizard(self.speclist, flag)

    def import_Lt_data(self, filename):

        with open(filename, "r") as infile:

            collector = []

            for line in infile:

                if line[0] != "#":
                    temp = line.strip().split("\t")
                    collector.append(temp)

            return collector

    def refine_collector(self, collector):

        time = []
        area = []
        apm = []
        surface_pressure = []

        for i in collector:
            time.append(i[1])
            area.append(i[2])
            apm.append(i[3])
            surface_pressure.append(i[4])

        return np.array(time), np.array(area), np.array(apm), np.array(surface_pressure)

    def create_lt_list(self, folder):

        lt_list = []

        for file in os.listdir(folder):

            if file.endswith(".dat"):


                path = folder + "/" + file
                creation_time = datetime.datetime.fromtimestamp(os.path.getmtime(path))
                data = self.import_Lt_data(path)
                data = list(self.refine_collector(data))
                data.append(creation_time)
                data.append(file[:-4])
                lt_list.append(data)

        return lt_list

    def write_lt_to_sql(self, folder):

        lt_list = self.create_lt_list(folder)
        db = sqlite3.connect("sfg.db")
        cur = db.cursor()
        for lt in lt_list:

            time = ";".join(lt[0].astype(str))
            area = ";".join(lt[1].astype(str))
            apm = ";".join(lt[2].astype(str))
            surface_pressure = ";".join(lt[3].astype(str))

            command = \
                """
                INSERT INTO lt_gasex
                (
                name,
                measured_time,
                time,
                area,
                apm,
                surface_pressure
                )
                VALUES(?,?,?,?,?,?);
                """
            try:

                cur.execute(command, (lt[-1], lt[-2], time, area, apm, surface_pressure))
            except sqlite3.IntegrityError as e:
                print("Spectrum already in database!")
        db.commit()
        db.close()




class DataCollector:
    # noinspection PyTypeChecker
    def __init__(self, filename):
        # the filename is ONLY the filename, not containing any directory information.
        self.file = filename
        self.creation_time = self.get_creation_time()

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
        sysname = SystematicName(self.file, self.creation_time)
        data = self.data_package
        sfg = SfgSpectrum(data[0], data[1], data[3], data[2], sysname)
        return sfg

    def get_creation_time(self):
        t = datetime.datetime.fromtimestamp(os.path.getmtime(self.file))
        return t


class SqlWizard:

    def __init__(self, speclist, flag="regular"):

        self.create_database()
        self.speclist = speclist

        for spectrum in self.speclist:

            if flag == "regular":
                self.add_spectrum(spectrum)
            elif flag == "gasex":
                self.add_gasex_spectrum(spectrum)

    def create_database(self):

        #Table one for SFGs without Sensitizer
        command =\
            """
            CREATE TABLE IF NOT EXISTS sfg_database (
            id INTEGER PRIMARY KEY,
            name TEXT,
            measured_time TIMESTAMP,
            measurer TEXT,
            wavenumbers TEXT,
            sfg TEXT,
            ir TEXT,
            vis TEXT,
            surfactant TEXT,
            sensitizer TEXT,
            photolysis TEXT,
            CONSTRAINT unique_name UNIQUE(name)
            );
            """



        command2 =\
        """
        CREATE TABLE IF NOT EXISTS sfg_gasex (
            id INTEGER PRIMARY KEY,
            name TEXT,
            measured_time TIMESTAMP,
            measurer TEXT,
            wavenumbers TEXT,
            sfg TEXT,
            ir TEXT,
            vis TEXT,
            CONSTRAINT unique_name UNIQUE(name)
            );
        """

        command3 = \
            """
            CREATE TABLE IF NOT EXISTS lt_gasex (
                id INTEGER PRIMARY KEY,
                name TEXT,
                measured_time TIMESTAMP,
                measurer TEXT,
                time TEXT,
                area TEXT,
                apm TEXT,
                surface_pressure TEXT,
                CONSTRAINT unique_name UNIQUE(name)
                );
            """


        db = sqlite3.connect("sfg.db")
        cur = db.cursor()

        cur.execute(command)
        cur.execute(command2)
        cur.execute(command3)
        db.commit()
        db.close()

    def add_spectrum(self, spectrum):

        s = spectrum  # type: SfgSpectrum
        fname = s.name # type: SystematicName

        name = fname.full_name
        time = fname.creation_time
        wavenumbers = ";".join(s.wavenumbers.astype(str))
        sfg = ";".join(s.raw_intensity.astype(str))
        ir = ";".join(s.ir_intensity.astype(str))
        vis = ";".join(s.vis_intensity.astype(str))

        surfactant = fname.surfactant
        sensitizer = fname.sensitizer
        photolysis = fname.photolysis

        db = sqlite3.connect("sfg.db")
        cur = db.cursor()

        command =\
        """
        INSERT INTO sfg_database
        (
        name,
        measured_time,
        wavenumbers,
        sfg,
        ir,
        vis,
        surfactant,
        sensitizer,
        photolysis)
        VALUES(?,?,?,?,?,?,?,?,?);
        """
        try:
            cur.execute(command, (name, time, wavenumbers, sfg, ir, vis,surfactant,sensitizer,photolysis))
        except sqlite3.IntegrityError as e:
            print("Spectrum already in database!")
        db.commit()
        db.close()

    def add_gasex_spectrum(self, spectrum):

        s = spectrum  # type: SfgSpectrum
        fname = s.name # type: SystematicName

        name = fname.full_name
        time = fname.creation_time
        wavenumbers = ";".join(s.wavenumbers.astype(str))
        sfg = ";".join(s.raw_intensity.astype(str))
        ir = ";".join(s.ir_intensity.astype(str))
        vis = ";".join(s.vis_intensity.astype(str))

        db = sqlite3.connect("sfg.db")
        cur = db.cursor()

        command =\
        """
        INSERT INTO sfg_gasex
        (
        name,
        measured_time,
        wavenumbers,
        sfg,
        ir,
        vis
        )
        VALUES(?,?,?,?,?,?);
        """
        try:
            cur.execute(command, (name, time, wavenumbers, sfg, ir, vis))
        except sqlite3.IntegrityError as e:
            print("Spectrum already in database!")
        db.commit()
        db.close()


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

Importer().generate_sql("regular")
I = Importer("gasex","gasex_out")
I.generate_sql("gasex")
I.write_lt_to_sql("gasex_lt")