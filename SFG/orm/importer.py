import datetime
import json
import os

import pandas as pd


class Importer:

    def __init__(self):

        # sfgs
        self.regular_sfg = self.import_sfg("C:/Users/lange/Desktop/CharmingSFG/SFG/newport/regular")
        self.gasex_sfg = self.import_sfg("C:/Users/lange/Desktop/CharmingSFG/SFG/newport/gasex_sfg")
        self.boknis = self.import_sfg("C:/Users/lange/Desktop/CharmingSFG/SFG/newport/boknis")

        # lts
        self.lt = self.import_lt("C:/Users/lange/Desktop/CharmingSFG/SFG/newport/lt")
        self.gasex_lt = self.import_lt("C:/Users/lange/Desktop/CharmingSFG/SFG/newport/gasex_lt")
        self.gasex_lift_off = self.import_liftoffs("C:/Users/lange/Desktop/CharmingSFG/SFG/newport/liftoff_points.csv")

        # spectra
        self.ir = []
        self.uv = []
        self.raman = []
        self.import_other_spectra()

        # surface tension and salinity
        self.gasex_tension = self.import_tensions("C:/Users/lange/Desktop/CharmingSFG/SFG/newport/gasex_surftens.txt")
        self.salinity = self.import_salinity()

        #substances
        self.substances = self.import_substances("C:/Users/lange/Desktop/CharmingSFG/SFG/newport/substances.json")

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

                    dic = {"name": name, "type": parent_dir.split("/")[-1], "measured_time": creation_time, "measurer": measurer,
                           "data": data}

                    if "dppc" in dic["name"] or "DPPC" in dic["name"]:

                        if parent_dir == "C:/Users/lange/Desktop/CharmingSFG/SFG/newport/boknis":
                            dic["type"] = "boknis_ref"

                        else:
                            dic["type"] = "regular"

                    out.append(dic)
                    with open("log.txt", "a") as outfile:
                        outfile.write(dic["name"]+"\n")

        return out

    def extract_sfg_file(self, file):
        """A function extracting the measurement data from a SFG spectral file """
        col_names = ['wavenumbers', 'sfg', 'ir', 'vis']
        temp = pd.read_csv(file, sep="\t", usecols=[0, 1, 3, 4], names=col_names, encoding='utf8', engine='python')
        return temp

    # LT
    def import_lt(self, head_dir):
        """Imports LT data and collects the raw data by calling the extract_lt_function().
        It adds addtional metadata like creation timestamp, measured and name."""

        out = []

        for parent_dir in os.listdir(head_dir):
            new_path = head_dir + "/" + parent_dir
            if os.path.isdir(new_path):

                for file in os.listdir(new_path):

                    if file.endswith(".dat"):
                        creation_time = datetime.datetime.fromtimestamp(os.path.getmtime(new_path + "/" + file))
                        name = file[:-4]

                        data = self.extract_lt_file(new_path + "/" + file)
                        dic = {"name": name, "type": head_dir.split("/")[-1], "measured_time": creation_time, "data": data}

                        out.append(dic)
        return out

    def extract_lt_file(self, file):
        """A function extracting the measurement data from a SFG spectral file """
        col_names = ["time", "area", "apm", "surface_pressure"]
        temp = pd.read_csv(file, comment="#", sep='\t', usecols=[1, 2, 3, 4], names=col_names)
        return temp

    def import_liftoffs(self, file):
        """Imports the lift-off points determined manually from a file"""
        col_names = ["name", "lift_off"]
        temp = pd.read_csv(file, sep=';', names=col_names, skiprows=1)
        return temp

    # other spectra
    def extract_other_spectra(self, file, sep='\t', skip=0):
        """Extracts the data from UV, IR and Raman measurements"""
        temp = pd.read_csv(file, sep=sep, skiprows=skip, names=["x", "y"])
        return temp

    # other
    def import_substances(self, file):
        """Extracts the information about predefined subtances out of the corresponding file"""
        with open(file) as infile:
            return json.load(infile)

    def import_tensions(self, file):
        col_names = ["name", "tension"]
        temp = pd.read_csv(file, sep=';', names=col_names)
        return temp

    def import_other_spectra(self):

        for file in os.listdir("C:/Users/lange/Desktop/CharmingSFG/SFG/newport/Raman"):
            temp_df = self.extract_other_spectra("C:/Users/lange/Desktop/CharmingSFG/SFG/newport/Raman/"+file)
            dic = {"name":file, "data": temp_df}
            self.raman.append(dic)

        for file in os.listdir("C:/Users/lange/Desktop/CharmingSFG/SFG/newport/IR"):
            temp_df = self.extract_other_spectra("C:/Users/lange/Desktop/CharmingSFG/SFG/newport/IR/"+file)
            dic = {"name": file, "data": temp_df}
            self.ir.append(dic)

        for file in os.listdir("C:/Users/lange/Desktop/CharmingSFG/SFG/newport/UV"):
            temp_df = self.extract_other_spectra("C:/Users/lange/Desktop/CharmingSFG/SFG/newport/UV/"+file, skip=2, sep=",")
            dic = {"name": file, "data": temp_df}
            self.uv.append(dic)

    def import_salinity(self):
        sal = pd.read_excel("C:/Users/lange/Desktop/CharmingSFG/SFG/newport/stationsplan.xls")
        salinities = []
        for row in range(len(sal)):
            sur_sal = sal.loc[row, "Salinity surface"]
            dep_sal = sal.loc[row, "Salinity depth"]
            sur_temp = sal.loc[row, "Temperature surface"]
            dep_temp = sal.loc[row, "Temperature depth"]
            _hash = '0' + str(sal.loc[row, "hash"])
            label = sal.loc[row, "Leg"] + "-" + str(sal.loc[row, "Station Number"])
            lat = sal.loc[row, "Latitude"]
            long = sal.loc[row, "Longitude"]

            dic = {"surface_salinity": sur_sal, "deep_salinity": dep_sal,
                   "label": label, "longitude": long, "latitude": lat, "hash": _hash,
                   "surface_temperature": sur_temp, "deep_temperature": dep_temp}
            salinities.append(dic)
        return salinities

    # auxiliary functions


# todo: decouple the Importer so it is not to tightly bound to boknis/gasex
# todo: dacopuling by using a factory method for the importer?