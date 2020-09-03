import datetime
import json
import os
import timeit
from functools import partial
from multiprocessing import Pool
from pathlib import Path
import re

import pandas as pd
import yaml


class Importer:
    """The Importer class performs the import steps from the raw measurement data."""

    def __init__(self):
        # paths
        self.base_path = Path.cwd() / 'newport'

        with open(Path.cwd() / "orm" / "import_paths.json") as infile:
            d = json.load(infile)
            self.paths = {key: (self.base_path / d[key]) for key in d}

        # sfgs
        p = Pool(3)
        data = p.map(self.import_sfg, [self.paths["regular_sfg"], self.paths["gasex_sfg"], self.paths["boknis"]])
        p.close()
        self.regular_sfg, self.gasex_sfg, self.boknis = data

        # lts
        p = Pool(2)
        data = p.map(self.import_lt, [self.paths["lt"], self.paths["gasex_lt"]])
        p.close()
        self.lt, self.gasex_lt = data
        self.gasex_lift_off = self.import_liftoffs(self.paths["gasex_lift_off"])

        # spectra
        p = Pool(3)
        data = map(self.import_xy_spectra, ["ir", "uv", "raman"])
        p.close()
        self.ir, self.uv, self.raman = data

        # surface tension, GasEx cruise data and Boknis Eck parameters
        self.gasex_tension = pd.read_csv(self.paths["gasex_tension"], sep=';', names=["name", "surface_tension"])

        self.station_plan = pd.read_excel(self.paths["station_plan"], dtype={'hash': str}).rename(
            columns={"Time [UTC]": "time", "Salinity surface": "salinity_surface", "Salinity depth": "salinity_depth",
                     "Temperature surface": "temperature_surface", "Temperature depth": "temperature_depth",
                     "Station Number": "station_number"})

        self.be_database_parameters = pd.read_csv(self.paths["be_database_parameters"],
                                                  sep=",", header=0, engine='c').rename(
            columns={"Depth [m]": "depth", "chlora": "chlorophyll_a"}).drop(
            columns=['chlora_flag', 'Name', 'Longitude', 'Latitude'])

        self.water_samples = self.import_water_samples()

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
        creation_time = datetime.datetime.fromtimestamp(os.path.getmtime(str(file)))
        data = func(file)
        dic = {"name": file.stem, "type": file.parent.parent.name, "measured_time": creation_time, "data": data}

        if sfg:
            date = file.parent.name.split(" ")[0]
            dic["name"] = f'{date}_{dic["name"]}'

            for key in ("wavenumbers", "sfg", "ir", "vis"):
                dic[key] = dic["data"][key]

            if "dppc" in dic["name"] or "DPPC" in dic["name"]:
                if "boknis" in file.parent.parent.name:
                    dic["type"] = "boknis_ref"
                else:
                    dic["type"] = "regular"
        else:
            for key in ("time", "area", "apm", "surface_pressure"):
                dic[key] = dic["data"][key]

        del dic['data']
        return dic

    @staticmethod
    def extract_sfg_file(file):
        """A function extracting the measurement data from a SFG spectral file """
        col_names = ['wavenumbers', 'sfg', 'ir', 'vis']
        temp = pd.read_csv(file, sep="\t", usecols=[0, 1, 3, 4], names=col_names, encoding='utf8', engine='c').apply(
            Importer.nparray_to_str)
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
        temp = pd.read_csv(file, comment="#", sep='\t', usecols=[1, 2, 3, 4], names=col_names, engine="c").apply(
            Importer.nparray_to_str)
        return temp

    @staticmethod
    def import_liftoffs(file):
        """Imports the lift-off points determined manually from a file"""
        with open(file) as infile:
            return pd.DataFrame(yaml.load(infile, Loader=yaml.FullLoader).items(), columns=["name", "lift_off"])

    # other spectra
    @staticmethod
    def extract_xy_spectra(file, sep='\t', skip=0):
        """Extracts the data from UV, IR and Raman measurements"""
        temp = pd.read_csv(file, sep=sep, skiprows=skip, names=["x", "y"]).apply(Importer.nparray_to_str)
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

    def import_water_samples(self):
        temp = pd.read_excel(self.paths["water_samples"],
                             header=2, sheet_name="Samples").drop(
            columns=['Quick description: Result', 'Remarks', 'Additional information', 'Quadrant']).rename(
            columns={'Sampler no.': 'sampler_no', 'Dips per sample': 'dips_per_sample',
                     'Sample container type': 'sample_container_type',
                     'Drainage time per dip (sec)\nSampling depth (m)': 'drainage_time_or_depth',
                     'Volume collected (ml)': 'volume_collected',
                     'Sea surface observational codes': 'sea_surface_observational_codes', 'Ship No.': 'ship_no',
                     'Location No.': 'location_no', 'Speed (m/s)': 'wind_speed', 'Direction (°)': 'wind_direction',
                     'Sample taken by': 'sampled_by', 'Period (m)': 'wave_period', 'Height (m)': 'wave_height',
                     'Air': 'air_temperature', 'Water': 'water_temperature', 'Experiment dates': 'experiment_date',
                     'Type of pollutant collected': 'pollutant_type', 'Date': 'sampling_date', 'Time': 'sampling_time'})
        temp["experiment_date"] = temp["experiment_date"].apply(Importer.to_measurement_date)
        return temp

    @staticmethod
    def nparray_to_str(array):
        return ",".join(array.values.astype(str))

    @staticmethod
    def to_measurement_date(string):
        rg = re.compile("^SFG \\(\d{4}-\d{2}-\d{2}\\)$")
        if re.match(rg, str(string)):
            temp = string.strip("SFG (").strip(")").split("-")
            d = datetime.date(int(temp[0]), int(temp[1]), int(temp[2]))
        else:
            d = None
        return d


if __name__ == "__main__":
    os.chdir("..")
    start = timeit.default_timer()
    i = Importer()
    end = timeit.default_timer()
    print(end - start)
    print(i.gasex_lift_off.keys())
