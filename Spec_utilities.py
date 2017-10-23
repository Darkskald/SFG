import os
import csv
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams['mathtext.default'] = 'regular'


class Spectrum:

    """Abstract class to manage the very basic properties of a spectra (wavelength/number, intensity and name
    in a convenient way. Covers IR, Raman and UV spectra at once"""
    def __init__(self, name, wavelengths, intensities):
        self.name = name
        self.wavelengths = wavelengths
        self.intensities = intensities

class UV_DataCollector:

    """Class to fetch UV spectral data from file. Handles the change to the UV directory internally"""

    def __init__(self,name):

        initial_dir = os.getcwd()
        os.chdir("UV")
        self.filename = name
        self.data = self.get_data(self.filename)
        self.spectrum = self.yield_uv()
        os.chdir(initial_dir)

    def get_data(self,filename):
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

    def yield_uv(self):

        name = self.filename[:-4]
        uv = Spectrum(name, self.data[0], self.data[1])
        return uv

class IR_DataCollector:
    """Class to fetch IR and Raman spectral data from file. Handles the change to the IR/Raman
    directory internally"""

    def __init__(self,name,type):

        initial_dir = os.getcwd()
        self.filename = name

        self.type = type
        if self.type == "IR":
            os.chdir("IR")
            self.data = self.get_iraman_data(self.filename)
            self.spectrum = self.yield_spec()
            os.chdir(initial_dir)

        elif self.type == "Raman":
            os.chdir("Raman")
            self.data = self.get_iraman_data(self.filename)
            self.spectrum = self.yield_spec()
            os.chdir(initial_dir)

    def get_iraman_data(self, filename):

        with open(filename, "r") as infile:
            # function can handle IR files as well as Raman files
            c = csv.reader(infile, delimiter="\t")
            wavenumbers = []
            intensities = []
            for line in c:
                wavenumbers.append(line[0])
                intensities.append(line[1])
            return (wavenumbers, intensities)

    def yield_spec(self):

        if self.type == "IR":
            name = self.filename[:-4]
        else:
            name = self.filename
        spec = Spectrum(name, self.data[0], self.data[1])
        return spec


class SpecDatabase:

    """A simple class for dealing with a certain set of experimental UV/Vis, IR and Raman spectra. Gives access
    to plotting, comparison and other useful functions. The type argument in the constructor determines
    the concrete procedure"""

    def __init__(self, type):

        self.database = []
        self.type = type
        if self.type == "UV":

            for file in os.listdir("UV"):
                if file.endswith(".csv"):
                    collect = UV_DataCollector(file).yield_uv()
                    self.database.append(collect)

        elif self.type == "IR":
            for file in os.listdir("IR"):
                    collect = IR_DataCollector(file, type="IR").yield_spec()
                    self.database.append(collect)

        elif self.type == "Raman":
            for file in os.listdir("Raman"):
                    collect = IR_DataCollector(file, type="Raman").yield_spec()
                    self.database.append(collect)

    def uv_plot(self):
        if self.type != "UV":
            raise TypeError("UV plot function not defined for type " + self.type)

        fig = plt.figure()
        ax = plt.subplot(111)
        for spectrum in self.database:
            ax.plot(spectrum.wavelengths, spectrum.intensities, label=spectrum.name)
        ax.set_xlabel("wavelength/ nm")
        ax.set_ylabel("absorbance")
        ax.legend()
        ax.set_xlim(250,700)
        plt.show()

    def ir_plot(self):
        if self.type != "IR" and self.type != "Raman":
            raise TypeError("IR plot function not defined for type " + self.type)
        fig = plt.figure()
        ax = plt.subplot(111)
        for spectrum in self.database:
            ax.plot(spectrum.wavelengths, spectrum.intensities, label=spectrum.name)
        ax.grid()
        ax.set_xlabel("Wavenumber/ $cm^{-1}$")
        ax.set_ylabel("Intensity/ a.u.")
        ax.legend()
        plt.show()



#test code section
s = SpecDatabase("Raman")
s.uv_plot()