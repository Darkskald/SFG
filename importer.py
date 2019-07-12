import datetime
import os

import pandas as pd


class Importer:

    def __init__(self):
        os.chdir("newport")
        # sfgs
        #self.regular_sfg = self.import_sfg("regular")
        #self.gasex_sfg = self.import_sfg("gasex_sfg")
        #self.boknis = self.import_sfg("boknis")

        # lts
        #self.lt = self.import_lt("lt")
        #self.gasex_lt = self.import_lt("gasex_lt")
        self.gasex_lift_off = self.import_liftoffs("liftoff_points.csv")

        # spectra
        #self.ir = None
        #self.uv = None
        #self.raman = None

        # surface tension
        self.gasex_tension = self.import_tensions("gasex_surftens.txt")

    # SFG
    def import_sfg(self, parent_dir):
        """SFG import function, calls extract_sfg_file() to collect the actual
        measurement file, extracts the file creation timestamp, measurer and name
        and packs it in a form suitable for persistence"""
        out = []

        for directory in os.listdir(parent_dir):

            date, measurer = directory.split(" ")

            for file in os.listdir(parent_dir + "/" + directory):

                if file.endswith(".sfg"):
                    creation_time = datetime.datetime.fromtimestamp(
                        os.path.getmtime(parent_dir + "/" + directory + "/" + file))
                    name = date + "_" + file[:-4]

                data = self.extract_sfg_file(parent_dir + "/" + directory + "/" + file)

                dic = {"name": name, "type": parent_dir, "measured_time": creation_time, "measurer": measurer,
                       "data": data}

                if "dppc" in dic["name"] or "DPPC" in dic["name"]:

                    if parent_dir == "boknis":
                        dic["type"] = "boknis_ref"

                    else:
                        dic["type"] = "regular"

                out.append(dic)

        return out

    def extract_sfg_file(self, file):
        """A function extracting the measurement data from a SFG spectral file """
        col_names = ['wavenumbers', 'sfg', 'ir', 'vis']
        temp = pd.read_csv(file, sep="\t", usecols=[0, 1, 3, 4], names=col_names, encoding='utf8', engine='python')
        return temp

    # LT
    def import_lt(self, parent_dir):
        """Imports LT data and collects the raw data by calling the extract_lt_function().
        It adds addtional metadata like creation timestamp, measured and name."""

        out = []

        for file in os.listdir(parent_dir):

            if file.endswith(".dat"):
                creation_time = datetime.datetime.fromtimestamp(os.path.getmtime(parent_dir + "/" + file))
                name = file[:-4]

                data = self.extract_lt_file(parent_dir + "/" + file)
                dic = {"name": name, "type": parent_dir, "measured_time": creation_time, "data": data}

                out.append(dic)

        return out

    def extract_lt_file(self, file):
        """A function extracting the measurement data from a SFG spectral file """
        col_names = ["time", "area", "apm", "surace_pressure"]
        temp = pd.read_csv(file, comment="#", sep='\t', usecols=[1, 2, 3, 4], names=col_names)
        return temp

    def import_liftoffs(self, file):
        """Imports the lift-off points determined manually from a file"""
        col_names = ["name", "lift_off"]
        temp = pd.read_csv(file, sep=';', names=col_names)
        return temp

    # other spectra
    def extract_other_spectra(self, file, sep='\t', skip=0):
        """Extracts the data from UV, IR and Raman measurements"""
        temp = pd.read_csv(file, sep=sep, skiprows=skip)
        return temp

    # other
    def import_substances(self, file):
        """Extracts the information about predefined subtances out of the corresponding file"""
        col_names = ["name", "short", "molar_mass", "sensitizing"]
        temp = pd.read_csv(file, sep=';', names=col_names)
        return temp

    def import_tensions(self, file):
        col_names = ["name", "tension"]
        temp = pd.read_csv(file, sep=';', names=col_names)
        return temp
