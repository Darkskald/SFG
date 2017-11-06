




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


a = get_data("BX6_2.csv")
b = get_data("BX9_1.csv")
c = get_data("BX12.csv")

plt.plot(a[0], a[1], label="BX6")
plt.plot(b[0], b[1], label="BX9")
plt.plot(c[0], c[1], label="BX12")
plt.legend()
plt.show()

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
a = get_iraman_data("BX9")
b = get_iraman_data("BX9.dpt")

d = FileFetcher("20170901_BX9_5_x1_#2_5mM.sfg").sfg
c = (d.wavenumbers, d.normalized_intensity)

fig, ax1 = plt.subplots()

ax1.plot(a[0], a[1])
ax1.set_xlabel("Wavenumbers")
ax1.set_ylabel("Raman intensity")
ax1.set_xlim(1000, 4000)

ax2 = ax1.twinx()  # creating another axis object sharing x axis, twiny() also possible
ax2.plot(b[0], b[1], color="r")
ax2.set_xlim(1000, 4000)

"""
ax3 = ax1.twinx()
ax3.plot(c[0],c[1],color="black")
ax2.set_xlim(1000,4000)
"""
ax2.set_ylabel("Transmission")

# fig.tight_layout()
plt.show()

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


#stolen from the Analyzer Class
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