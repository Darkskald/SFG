import os
import csv
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams

rcParams['mathtext.default'] = 'regular'

class UV_DataCollector:
    """Class to fetch UV spectral data from file. Handles the change to the UV directory internally"""

    def __init__(self, name):

        initial_dir = os.getcwd()
        os.chdir("UV")
        self.filename = name
        self.data = self.get_data(self.filename)
        self.spectrum = self.yield_uv()
        os.chdir(initial_dir)

    def get_data(self, filename):
        with open(filename, "r") as infile:
            next(infile)
            next(infile)
            c = csv.reader(infile, delimiter="\t")
            wavelength = []
            intensity = []
            for line in c:
                l = line[0].split(",")
                if len(l) == 2:
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

    def __init__(self, name, type):

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
        ax.set_xlim(250, 700)
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


class Spectrum:
    """Simple class to manage the very basic properties of a spectra (wavelength/number, intensity and name
       in a convenient way. Covers IR, Raman and UV spectra at once"""
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