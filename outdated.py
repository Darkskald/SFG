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


class TexInterface:
    def __init__(self, name, sfg_list):
        self.name = name
        self.blocks = []
        self.out_dir = "tex_out/"
        self.fig_dir = "fig/"
        self.sfg_list = sfg_list

        # section control
        self.create_header()

    def add_figure(self, filename, label):
        figstring = "\\begin{figure}[h!]\n\centering\n\includegraphics[scale=0.4]"
        figstring += "{" + self.fig_dir + filename + "}\n"
        figstring += "\caption{" + label + ".}\n\end{figure}\n"
        self.blocks.append(figstring)

    def add_tabular(self, spectrum):
        corr_name = " ".join(spectrum.name.full_name.split("_"))
        corr_name = corr_name.replace("#", "no")
        tabular_string = "\\begin{table}[h!]\n"
        tabular_string += "\caption{" + "Peaklist of " + corr_name + "}" + "\n"
        tabular_string += "\\begin{center}\n"

        tabular_string += spectrum.drop_tex_peaktable()

        tabular_string += "\\end{center}\n"
        tabular_string += "\\end{table}\n"
        print("DEBUUUUUUUUUG")
        self.blocks.append(tabular_string)

    def load_template(self, name):
        pass

    def create_header(self):
        with open("texhead", "r") as infile:
            string = infile.read()
            self.blocks.append(string)

    def write_file(self):
        self.blocks.append("\n" + "\end{document}")
        outstring = "\n".join(self.blocks)
        with open(self.out_dir + self.name + ".tex", "w") as outfile:
            outfile.write(outstring)

    def compile(self):
        pass

    def add_section(self, name):
        self.blocks.append("\section{" + name + "}" + "\n")

    def traverse_sfg_list(self):
        for spectrum in self.sfg_list:
            corr_name = " ".join(spectrum.name.full_name.split("_"))
            corr_name = corr_name.replace("#", "NR. ")
            self.add_section(corr_name)
            P = Plotter([spectrum])
            s = spectrum.name.full_name[:-4].replace(".", "x")
            s = s.replace("#", "no")
            P.title = s
            P.simple_plot(mode="save")
            self.add_figure(s + ".pdf", "Simple plot with normalized intensity")
            self.add_tabular(spectrum)

    def join_speclist(self):
        P = Plotter(self.sfg_list)
        P.title = self.name
        P.simple_plot(mode="save")
        self.add_section(self.name)
        self.add_figure(self.name + ".pdf", "Joint plot of a list of spectra")


class Plotter:
    """Simple class providing an interface between the matplotlib API and the SFG objects. Depending on the specific method
    and kwargs, a variety of different methods is accessible"""
    def __init__(self, speclist, raw=False, title="Default"):
        self.speclist = speclist
        self.raw = raw
        self.title = title
        self.save_dir = "tex_out/fig"

    def simple_plot(self, mode="show"):

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
        if mode == "show":
            plt.show()
        elif mode =="save":
            plt.savefig(self.save_dir+"/"+self.title+".pdf")

    def raw_plot(self, mode="show"):
        for spectrum in self.speclist:
            wl = spectrum.wavenumbers
            intensity = spectrum.raw_intensity
            plt.plot(wl, intensity, label=spectrum.name.full_name, linestyle='--', markersize=4, marker="o")

        plt.grid(True)
        plt.title(self.title)
        plt.xlabel("Wavenumber/ $cm^{-1}$")
        plt.ylabel("Raw intensity/ a.u.")
        if len(self.speclist) < 6:
            plt.legend(loc="upper right")
        else:
            plt.legend(bbox_to_anchor=(1, 0.5), loc='center left', numpoints=1)
        if mode == "show":
            plt.show()
        elif mode =="save":
            plt.savefig(self.save_dir+"/"+self.title+".pdf")

    def raw_plot_plus_ir(self, mode="show"):
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
        if mode == "show":
            plt.show()
        elif mode =="save":
            plt.savefig(self.save_dir+"/"+self.title+".pdf")

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

    def stack_plot(self, mode="show"):
        spacer = 0
        fig = plt.figure()
        ax = plt.subplot(111)
        for spectrum in self.speclist:
            fig = plt.figure()
            ax = plt.subplot(111)
            for spectrum in self.speclist:
                wl = spectrum.wavenumbers
                intensity = spectrum.normalize_to_highest()+spacer
                ax.plot(wl, intensity, label=spectrum.name.full_name, linestyle='--', markersize=4, marker="o")
                spacer += 0.5

            box = ax.get_position()
            ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
            ax.grid()
            ax.set_title(self.title)
            ax.set_xlabel("Wavenumber/ $cm^{-1}$")
            ax.set_ylabel("Intensity/ a.u.")
            ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
            ax.axes.get_yaxis().set_ticks([])
            self.log_plot("stack")
            if mode == "show":
                plt.show()
            elif mode == "save":
                plt.savefig(self.save_dir + "/" + self.title + ".pdf")

    def log_plot(self, mode):
        # Part 1: Generate logfile
        with open(self.title+"_"+mode, "w") as outfile:
            outfile.write("Log file for plot " + self.title + "\n")
            outfile.write("The mode of the plot is: " + mode + "\n")
            outfile.write("Plot contains the following spectra: \n")
            for spectrum in self.speclist:
                outfile.write(spectrum.name.full_name + "\n")

    def bar_peaks(self, number=4):

        A = Analyzer(self.speclist)
        dic = A.count_peak_abundance(number)


        for key in dic:
            plt.bar(key, dic[key], width=3.)

        plt.xlabel("wavenumber")
        plt.ylabel("rel. peak abundance")
        plt.grid(True)
        plt.minorticks_on()
        plt.show()

    def marked_peaks(self, mode="show"):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax2 = ax.twiny()

        for spectrum in self.speclist:
            wl = spectrum.wavenumbers
            if self.raw == False:
                intensity = spectrum.normalized_intensity
            else:
                intensity = spectrum.raw_intensity
            ax.plot(wl, intensity, label=spectrum.name.full_name, linestyle='--', markersize=4, marker="o")

            peakticks = []
            for peak in spectrum.detailed_analysis(threshold=1.2):
                ax2.axvline(peak[0], color="red")
                ax2.axvline(peak[1], ls="dashed", color="blue")
                ax2.axvline(peak[2], ls="dashed", color="blue")
                peakticks.append(peak[0])

        ax.grid()
        ax.set_xlabel("Wavenumber/ $cm^{-1}$")
        ax.set_ylabel("Intensity/ a.u.")
        ax.legend()


        ax2.set_xticks(peakticks)
        ax2.set_xticklabels(peakticks, rotation=45, color='blue')
        ax2.set_xlim(ax.get_xlim())

        if mode == "show":
            plt.show()
        elif mode == "save":
            plt.savefig(self.save_dir + "/" + self.title + ".pdf")

    def relation_coverage_peaks(self):

        A = Analyzer(self.speclist)
        fig = plt.figure()
        ax = plt.subplot(111)
        for entry in A.intensity_vs_surfcoverage():
            ax.scatter(entry[0], entry[1])
        ax.set_xlabel("area per molecule/ A$^{2}$")
        ax.set_ylabel("maximum intensity")
        plt.show()

    def relation_3d(self):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        wavelengths = []
        intensities = []
        surface_concentrations = []

        for spectrum in self.speclist:
            for peak in spectrum.peaks:
                wavelengths.append(peak[0])
                intensities.append(peak[3])
                surface_concentrations.append(spectrum.calc_area_per_molecule())

        ax.scatter(wavelengths, intensities, surface_concentrations)
        ax.set_xlabel("wavelength")
        ax.set_ylabel("intensity")
        ax.set_zlabel("area per molecule")
        plt.show()


class FileFetcher:
    """The file fetcher class interconnects different subdirectories of the folder and simplifies
    the handling of filenames and filepaths. It changes the working directory to the file storage directory,
    usually the library, creates a DataCollector object and returns a SFG object"""

    def __init__(self, filename, destination="library"):
        self.filename = filename
        initial_wd = os.getcwd()
        os.chdir(destination)

        self.collector = DataCollector(filename)
        self.sfg = self.collector.yield_SfgSpectrum()
        os.chdir(initial_wd)


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



def simple_analysis():
    Substances = {}
    i = Ipy_Interpreter()

    with open("name_info/Surfactants.txt", "r") as infile:
        for line in infile:
            collect = line.split(":")
            Substances[collect[0]] = collect[1].strip()

    with open("name_info/Sensitizers.txt", "r") as infile:
        for line in infile:
            collect = line.split(":")
            Substances[collect[0]] = collect[1].strip()

    for key in Substances:
        if key != "DPPC":
            dates = []
            # erster Schritt: alle Daten des Surfactants holen
            i.get("su " + key)

            # zweiter Schritt: Datum extrahieren
            for j in i.subset:
                dates.append(j.name.date)

            # an jedem Tag alle Probennummern extrahieren
            for d in dates:
                i.get("su " + key)
                daysamples = []
                i.refine("d " + d)
                for spectrum in i.subset:
                    daysamples.append(spectrum.name.sample_number)

                    # die einzelnen Samples zusammen plotten
                for sample in daysamples:
                    plotllist = [q for q in i.subset if q.name.sample_number == sample]
                    if len(plotllist) > 1:
                        plot_title = key + " " + str(d) + " " + str(sample)
                        P = Plotter(plotllist, title=plot_title)
                        P.custom_plot()

                i.clear()

        # subset wieder säubern
        i.clear()
"""
def yield_peaklist():
    yield a defined list of peaks separated from each other by minimum the threshold value in wavenumber
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
"""
def calc_dish_area(diameter):
    """A auxialiary function to calculate the area of a teflon dish in square angstroms. Diameter given in cm."""
    radius = diameter * 0.5
    area = np.pi * radius ** 2
    area = area * 10 ** 16  # conversion to square angstroms
    return area


def calc_area_per_molecule(area, concentration, volume):
    """The function calculates the area per molecule. The area should be given in square angstroms, the
    concentration in milimole per liter and the volume in microleter"""

    concentration = concentration * 10 ** -3  # conversion in mol per liter
    volume = volume * 10 ** -6  # conversion in liter
    amount = volume * concentration
    molecules = (6.022 * 10 ** 23) * amount  # number of molecules
    area_per_molecule = area / molecules

    return area_per_molecule

"""Former UV file"""


class UV_Spectrum:
    def __init__(self, name, wavelengths, intensities):
        self.name = name
        self.wavelengths = wavelengths
        self.intensities = intensities


def get_data(filename):
    with open(filename, "r") as infile:
        next(infile)
        next(infile)
        c = csv.reader(infile, delimiter="\t")
        wavelength = []
        intensity = []
        for line in c:
            l = line[0].split(",")
            if len(l) == 2:
                print(l)
                wavelength.append(float(l[0]))
                intensity.append(float(l[1]))
    return (wavelength, intensity)


"""Former IR file"""


class IR_Spectrum:
    def __init__(self, name, wavenumbers, intensities):
        self.name = name
        self.wavenumbers = wavenumbers
        self.intensities = intensities


class Raman_Spectrum:
    def __init__(self, name, wavenumbers, intensities):
        self.name = name
        self.wavenumbers = wavenumbers
        self.intensities = intensities

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
            if spectrum.name.sensitizer == sensitizerflag:
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
class IR_Plotter:
    """This class takes a list of IR,Raman or SFG Spectra as input and provides plotting functionality"""

    def __init__(self, speclist, title="Default"):
        self.speclist = speclist
        self.title = title

    def simple_plot(self):
        """This function is made to overlap spectra of a common type and does not provide double y axis!"""

        fig = plt.figure()
        ax = plt.subplot(111)
        for spectrum in self.speclist:
            ax.plot(spectrum.wavenumbers, spectrum.intensities)

        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
        ax.grid()
        ax.set_title(self.title)
        ax.set_xlabel("Wavenumber/ $cm^{-1}$")
        ax.set_ylabel("Intensity/ a.u.")
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        plt.show()


class SqlImporter:
        # first: make a list of day folders in the archive directory
        def __init__(self, newspec_directory):
            self.speclist = []
            self.directory = newspec_directory
            for file in os.listdir(self.directory):
                if file.endswith("sfg"):
                    os.chdir("library")
                    D = DataCollector(file).yield_SfgSpectrum()
                    os.chdir("..")
                    self.speclist.append(D)

            S = SqlWizard(self.speclist)



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

colors = {
    0: "b",
    1: "g",
    2: "c",
    3: "m",
    4: "y",
    5: "k",
}


def detailed_analysis(x_array, y_array):
    x_array = x_array[::-1]
    y_array = y_array[::-1]

    slopes = [(y_array[i + 1] - y_array[i]) / (x_array[i + 1] - x_array[i]) for i in range((len(x_array) - 1))]
    possible_peaks = []
    peak_tuples = []

    for i in range(1, len(y_array) - 1):
        if slopes[i - 1] > 0 and slopes[i] < 0:
            possible_peaks.append(i)

    average_intensity = np.average(y_array)

    confirmed_peaks = [i for i in possible_peaks if (y_array[i] > average_intensity * 5)]

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
        print("center: " + str(x_array[i[0]]))
        center = x_array[i[0]]
        print("left: " + str(x_array[i[1]]))
        left = x_array[i[1]]
        print("right: " + str(x_array[i[2]]))
        right = x_array[i[2]]
        peak_slice_x = x_array[i[1]:i[2] + 1]
        peak_slice_y = y_array[i[1]:i[2] + 1]
        area = integrate_peak(peak_slice_x, peak_slice_y)
        datapoints = len(peak_slice_x)

        data_out.append((center, left, right, peak_slice_x, peak_slice_y, datapoints, area, indices))
    return data_out


def integrate_peak(x_array, y_array):
    area = 0
    for i in range(len(x_array) - 1):
        dx = abs(x_array[i + 1] - x_array[i])
        print(dx)
        square = dx * y_array[i]
        triangle = dx * abs(y_array[i + 1] - y_array[i]) * 0.5
        total = square + triangle
        area += total
    return area


# test code section
"""
i = Ipy_Interpreter()
i.get("su DPPC")
i.keep("80")
i.subset[0].smooth()
a = i.subset[0].wavenumbers
b = i.subset[0].normalized_intensity
q=detailed_analysis(a,b)
print(q)
i.plot()


i = Ipy_Interpreter()
i.get("su DPPC")

for j in range(len(i.subset)):
    i.subset[j].smooth()
    a = i.subset[j].wavenumbers
    b = i.subset[j].normalized_intensity
    try:
        q=detailed_analysis(a,b)
        print(q)
    except LookupError:
        print("Something happened with indexing")  
x = np.linspace(-0.5,3.3,500)
y = 0.8*np.sin(2*x)**2

print("Hello")
"""
i = Ipy_Interpreter()
i.get("su DPPC")
i.refine("d 20170914")
sfg = i.subset[0]

x = sfg.wavenumbers
y = sfg.normalized_intensity

analyse = detailed_analysis(x, y)

x = sfg.wavenumbers[::-1]
y = sfg.normalized_intensity[::-1]

fig = plt.figure()
ax = plt.subplot(111)
ax.plot(x, y)
ax.set_ylim(0)
j = 0
for i in analyse:
    print("Test:")
    print(i)
    center = i[0]
    left = i[1]
    right = i[2]
    area = i[6]
    indices = i[7]
    print(indices)
    centerx = x[indices[0]]
    centery = y[indices[0]]
    print(centerx, center)

    ax.fill_between(i[3], i[4])
    ax.axvline(left, color=colors[j], linewidth=2)
    ax.axvline(right, color=colors[j], linewidth=2)

    ax.annotate((float("{0:.2f}".format(centerx))), (centerx, centery)

                )
    ax.annotate("area: " + str((float("{0:.2f}".format(area)))), (centerx - 0.2, centery * 0.3))
    j += 1

plt.show()


def get_iraman_data(filename):
    with open(filename, "r") as infile:
        # function can handle IR files as well as Raman files
        c = csv.reader(infile, delimiter="\t")
        wavenumbers = []
        intensities = []
        for line in c:
            wavenumbers.append(line[0])
            intensities.append(line[1])
        return (wavenumbers, intensities)


# test code section


"""not yet implemented classes"""


class Day_meta:
    pass


class Interpreter:
    """Class providing the functionality for the command line of the plotting routine"""

    def __init__(self, command):

        commandlist = command.split(" ")
        if len(commandlist) == 1:
            if commandlist[0] == "ud":
                self.update()
            elif commandlist[0] == "x":
                sys.exit()

        elif len(command) < 3:
            print("Wrong number of parameters! Try again")
        else:
            self.type = commandlist[0]
            self.flags = commandlist[1]
            self.options = commandlist[2]

            if self.type == "plot" and self.flags == "f":
                files = self.options.split(",")
                self.plothandler(files)

            elif self.type == "plot" and self.flags == "d,su,sa":
                self.sample_surfactant_date(self.options)

            elif self.type == "plot" and self.flags == "bo":
                self.plot_boknis()

            elif self.type == "list":
                self.list_()

    def plothandler(self, filelist):

        sfg_objects = []
        for file in filelist:
            sfg = FileFetcher(file).sfg
            sfg_objects.append(sfg)
        P = Plotter(sfg_objects)
        P.simple_plot()

    def sample_surfactant_date(self, options):
        options = options.split(",")
        if len(options) != 3:
            print("sa_su_da calling. Invalid number of flags!")
        else:
            dateflag = options[0]
            surflag = self.retranslate_name(options[1])
            sampleflag = int(options[2])

            f = Finder()
            f1 = f.date_based(dateflag)
            f2 = f.surfactant_based(surflag, f1)
            f3 = f.sample_based(sampleflag, f2)

            Plotter(f3).simple_plot()

    def plot_boknis(self):
        if self.options == "none":
            f = Finder()
            sfg = f.comment_based()
            Plotter(sfg).simple_plot()

        elif "s" in self.options:
            number = int(self.options[1])
            f = Finder()
            sfg = f.comment_based()
            sfg = f.sample_based(number, sfg)
            Plotter(sfg).simple_plot()

        else:
            print("Not yet implemented")
            print(self.options, len(self.options))

    def update(self):
        Importer()

    def retranslate_name(self, stri):
        # load the allowed surfactans and sensitizers from files
        self.Surfactants = {}
        self.Sensitizers = {}

        with open("name_info/Surfactants.txt", "r") as infile:
            for line in infile:
                collect = line.split(":")
                self.Surfactants[collect[0]] = collect[1].strip()

        with open("name_info/Sensitizers.txt", "r") as infile:
            for line in infile:
                collect = line.split(":")
                self.Sensitizers[collect[0]] = collect[1].strip()
        if stri in self.Surfactants:
            return self.Surfactants[stri]
        elif stri in self.Sensitizers:
            return self.Sensitizers[stri]
        else:
            print("Retranslation failed. Unknown expression.")

    def list_(self):
        f = Finder()

        if self.flags == "su":
            surfoptions = self.options.split(",")  # if comma-separated multiplit sensitizers are given
            subsets = []
            for i in surfoptions:
                i = self.retranslate_name(i)
                f1 = f.surfactant_based(i)
                for j in f1:
                    subsets.append(j)
            self.answerset = subsets

        elif self.flags == "d":

            subsets = []
            dateoptions = self.options.split(",")
            for i in dateoptions:
                subset = f.date_based(self.options)
                for j in subset:
                    subsets.append(j)
            self.answerset = subsets


# stolen from the Analyzer Class
"""def list_peaks(self, number):
    intensities = []
    wavenumbers = []

    for i in self.speclist:
        j = i.yield_peaklist(num=number)
        for k in j:
            intensities.append(k[0])
            wavenumbers.append(k[1])
    # plt.scatter(wavenumbers,intensities)
    plt.hist(wavenumbers, rwidth=0.02, normed=True)
    plt.show()"""