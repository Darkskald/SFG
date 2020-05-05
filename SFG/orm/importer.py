import datetime
import json
import os
import timeit
from functools import partial
from multiprocessing import Pool
from pathlib import Path

import pandas as pd

os.chdir("..")


class Importer:
    """The Importer class performs the import steps from the raw measurement data."""

    def __init__(self):
        # paths
        self.base_path = Path.cwd() / 'newport'
        # todo dieser Pfadspaß ist miserabel, das sollte durch ein yaml oder json geladen werden
        self.paths = {
            "regular_sfg": self.base_path / 'regular',
            "gasex_sfg": self.base_path / 'gasex_sfg',
            "boknis": self.base_path / 'boknis',
            "lt": self.base_path / 'lt',
            "gasex_lt": self.base_path / 'gasex_lt',
            "gasex_lift_off": self.base_path / 'liftoff_points.csv',
            "gasex_tension": self.base_path / 'gasex_surftens.txt',
            "substances": self.base_path / 'substances.json',
            "ir": self.base_path / 'newport' / 'IR',
            "uv": self.base_path / 'UV',
            "raman": self.base_path / 'Raman',
            "station_plan": self.base_path / 'stationsplan.xls',
            "be_database_parameters": self.base_path / 'be_data.csv',
            "water_samples": self.base_path / 'Wasserproben_komplett.xlsx'
        }
        # sfgs
        p = Pool(3)
        data = p.map(self.import_sfg, [self.paths["regular_sfg"], self.paths["gasex_sfg"], self.paths["boknis"]])
        p.close()
        self.regular_sfg, self.gasex_sfg, self.boknis = data

        # lts
        self.lt = self.import_lt(self.paths["lt"])
        self.gasex_lt = self.import_lt(self.paths["gasex_lt"])
        self.gasex_lift_off = self.import_liftoffs(self.paths["gasex_lift_off"])

        # spectra
        self.ir = self.import_xy_spectra("ir")
        self.uv = self.import_xy_spectra("uv")
        self.raman = self.import_xy_spectra("raman")

        # surface tension, GasEx cruise data and Boknis Eck parameters
        self.gasex_tension = pd.read_csv(self.paths["gasex_tension"], sep=';', names=["name", "tension"])
        self.station_plan = pd.read_excel(self.paths["station_plan"])
        self.be_database_parameters = pd.read_csv(self.paths["be_database_parameters"], sep=",", header=0, engine='c')
        self.water_samples = pd.read_excel(self.paths["water_samples"])

        # substances
        self.substances = self.import_substances(self.paths["substances"])

    # SFG
    def import_sfg(self, parent_dir):
        """SFG import function, calls extract_sfg_file() to collect the actual
        measurement file, extracts the file creation timestamp, measurer and name
        and packs it in a form suitable for persistence"""
        out = []

        for directory in parent_dir.iterdir():
            files = directory.rglob('*.sfg')
            data = list(map(partial(self.aux, func=self.extract_sfg_file, sfg=True), files))
            out.extend(data)
        return out

    def aux(self, file, func, sfg):

        creation_time = datetime.datetime.fromtimestamp(file.stat().st_ctime)
        data = func(file)
        dic = {"name": file.name, "type": file.parent.parent, "measured_time": creation_time, "data": data}

        if sfg:
            if "dppc" in dic["name"] or "DPPC" in dic["name"]:
                if "boknis" in file.parent.parent.name:
                    dic["type"] = "boknis_ref"
                else:
                    dic["type"] = "regular"
        return dic

    @staticmethod
    def extract_sfg_file(file):
        """A function extracting the measurement data from a SFG spectral file """
        col_names = ['wavenumbers', 'sfg', 'ir', 'vis']
        # todo: changed this to c engine for reasons of performance
        temp = pd.read_csv(file, sep="\t", usecols=[0, 1, 3, 4], names=col_names, encoding='utf8', engine='c')
        return temp

    # LT
    def import_lt(self, head_dir):
        """Imports LT data and collects the raw data by calling the extract_lt_function().
        It adds addtional metadata like creation timestamp, measured and name."""

        out = []

        for parent_dir in head_dir.iterdir():

            if parent_dir.is_dir():
                files = parent_dir.rglob('*.dat')
                data = list(map(partial(self.aux, func=self.extract_lt_file, sfg=False), files))
                out.extend(data)
        return out

    @staticmethod
    def extract_lt_file(file):
        """A function extracting the measurement data from a SFG spectral file """
        col_names = ["time", "area", "apm", "surface_pressure"]
        temp = pd.read_csv(file, comment="#", sep='\t', usecols=[1, 2, 3, 4], names=col_names, engine="c")
        return temp

    @staticmethod
    def import_liftoffs(file):
        """Imports the lift-off points determined manually from a file"""
        col_names = ["name", "lift_off"]
        temp = pd.read_csv(file, sep=';', names=col_names, skiprows=1)
        return temp

    # other spectra
    @staticmethod
    def extract_xy_spectra(file, sep='\t', skip=0):
        """Extracts the data from UV, IR and Raman measurements"""
        temp = pd.read_csv(file, sep=sep, skiprows=skip, names=["x", "y"])
        return {"name": file.name, "data": temp}

    def import_xy_spectra(self, spec_type):
        if spec_type != "uv":
            data = map(self.extract_xy_spectra, self.paths[spec_type].iterdir())
        else:
            data = map(partial(self.extract_xy_spectra, skip=2, sep=","), self.paths[spec_type].iterdir())
        return data

    # other
    @staticmethod
    def import_substances(file):
        """Extracts the information about predefined substances out of the corresponding file"""
        with open(file) as infile:
            return json.load(infile)


start = timeit.default_timer()
i = Importer()
end = timeit.default_timer()
print(end - start)
# print(i.regular_sfg)
