import datetime
import json
import os
import timeit
from pathlib import Path

import pandas as pd


class Importer:

    """The Importer class performs the import steps from the raw measurement data."""

    def __init__(self):
        # paths
        self.base_path = Path.cwd()
        self.paths = {
            "regular_sfg": self.base_path / 'newport' / 'regular',
            "gasex_sfg": self.base_path / 'newport' / 'gasex_sfg',
            "boknis": self.base_path / 'newport' / 'boknis',
            "lt": self.base_path / 'newport' / 'lt',
            "gasex_lt": self.base_path / 'newport' / 'gasex_lt',
            "gasex_lift_off": self.base_path / 'newport' / 'liftoff_points.csv',
            "gasex_tension": self.base_path / 'newport' / 'gasex_surftens.txt',
            "substances": self.base_path / 'newport' / 'substances.json',
            "ir": self.base_path / 'newport' / 'IR',
            "uv": self.base_path / 'newport' / 'UV',
            "raman": self.base_path / 'newport' / 'Raman',
            "station_plan": self.base_path / 'newport' / 'stationsplan.xls'

        }
        # sfgs
        self.regular_sfg = self.import_sfg(str(self.paths["regular_sfg"]))
        self.gasex_sfg = self.import_sfg(str(self.paths["gasex_sfg"]))
        self.boknis = self.import_sfg(str(self.paths["boknis"]))

        # lts
        self.lt = self.import_lt(str(self.paths["lt"]))
        self.gasex_lt = self.import_lt(str(self.paths["gasex_lt"]))
        self.gasex_lift_off = self.import_liftoffs(str(self.paths["gasex_lift_off"]))

        # spectra
        self.ir = []
        self.uv = []
        self.raman = []
        self.import_other_spectra()

        # surface tension and salinity
        self.gasex_tension = self.import_tensions(str(self.paths["gasex_tension"]))
        self.salinity = self.import_salinity()

        #substances
        self.substances = self.import_substances(str(self.paths["substances"]))

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

                    dic = {"name": name, "type": parent_dir.split("/")[-1], "measured_time": creation_time,
                           "measurer": measurer,
                           "data": data}

                    if "dppc" in dic["name"] or "DPPC" in dic["name"]:

                        if "boknis" in parent_dir:
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
                        dic = {"name": name, "type": head_dir.split("/")[-1], "measured_time": creation_time,
                               "data": data}

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
        """Extracts the information about predefined substances out of the corresponding file"""
        with open(file) as infile:
            return json.load(infile)

    def import_tensions(self, file):
        col_names = ["name", "tension"]
        temp = pd.read_csv(file, sep=';', names=col_names)
        return temp

    def import_other_spectra(self):

        for file in os.listdir(str(self.paths["raman"])):
            temp_df = self.extract_other_spectra(str(self.paths["raman"] / file))
            dic = {"name": file, "data": temp_df}
            self.raman.append(dic)

        for file in os.listdir(str(self.paths["ir"])):
            temp_df = self.extract_other_spectra(str(self.paths["ir"] / file))
            dic = {"name": file, "data": temp_df}
            self.ir.append(dic)

        for file in os.listdir(str(self.paths["uv"])):
            temp_df = self.extract_other_spectra(str(self.paths["uv"] / file), skip=2, sep=",")
            dic = {"name": file, "data": temp_df}
            self.uv.append(dic)

    def import_salinity(self):
        sal = pd.read_excel(self.paths["station_plan"])
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
# todo: implement async io

if __name__ == "__main__":
    os.chdir("..")
    start = timeit.default_timer()
    i = Importer()
    end = timeit.default_timer()
    print(end - start)
