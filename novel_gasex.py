from spectrum import SfgSpectrum, LtIsotherm, DummyPlotter
import pandas as pd
import sqlite3
import numpy as np
from datetime import datetime

conn = sqlite3.connect("test.db")


class Station:

    def __init__(self, data):
        self.data = data
        self.samples = []

        self.stats = {
            "plate_coverage": None,
            "plate_lift_off": None,
            "plate_tension": None,
            "plate_max_pressure": None,
            "screen_coverage": None,
            "screen_lift_off": None,
            "screen_tension": None,
            "screen_max_pressure": None,
            "sml_coverage": None,
            "sml_lift_off": None,
            "sml_tension": None,
            "sml_max_pressure": None,
            "deep_coverage": None,
            "deep_lift_off": None,
            "deep_tension": None,
            "deep_max_pressure": None
        }
        # todo: calculate average coverage, tension, max_pressure, lift_off
        # todo: for sml, screen, plate and deep

    def get_date(self):
        date = datetime.strptime(self.data["date"], '%Y-%m-%d').date()
        return date


class Sample:

    def __init__(self, data):
        self.data = data
        self.sfg = []
        self.lt = []
        self.tension = None

        self.sfg_spectra = []
        self.lt_isotherms = []

    # todo: invoke the calculations necessary for the station

    def get_tension(self):
        if len(self.tension) > 1:
            raise ValueError("Invalid number of surface tension values!")

        elif len(self.tension) == 0:
            return None

        else:
            return self.tension.loc[0, "surface_tension"]

    def make_sfg(self):

        if len(self.sfg) > 1:
            raise ValueError(f'Invalid number of SFG spectra for {self.sfg["name"]}!')

        else:
            if len(self.sfg["name"]) > 0:

                name = (self.sfg["name"].values[0][0])
                time = datetime.strptime((self.sfg.loc[0, "measured_time"]), '%Y-%m-%d %H:%M:%S')
                meta = {"name": name, "time": time}

                # pass a tuple of column names to the map function, fetching the corresponding data from the datframe
                # and convert them to numpy arrays. Pass them to the SfgSpectrum constructor and unpack them
                data = ("wavenumbers", "sfg", "ir", "vis")
                tup = map(lambda x: np.fromstring(self.sfg[x].values[0], sep=";"), data)
                s = SfgSpectrum(*tup, meta)
                self.sfg_spectra.append(s)

    def make_lt(self):

        for i in range(len(self.lt["name"])):
            name = self.lt.loc[i, "name"].values[0]
            time = datetime.strptime((self.lt.loc[i, "measured_time"]), '%Y-%m-%d %H:%M:%S')

            lift_off = self.lt.loc[i, "lift_off"]
            if lift_off is not None:
                lift_off = float(lift_off)

            data = ("time", "area", "apm", "surface_pressure")
            tup = map(lambda x: np.fromstring(self.lt.loc[i, x], sep=";"), data)
            l = LtIsotherm(name, time, *tup, lift_off=lift_off)
            self.lt_isotherms.append(l)


# todo: encapsulate SQL interaction to class --> replace session controll manager
# todo: reimplement a top-level function which provides DPPC calibration data for each day of measurement


# functional part of the module to operate on the database

def get_stations():
    cmd = "SELECT * FROM stations"
    stations = pd.read_sql(cmd, conn)

    out = []
    counter = 0

    for item in stations["id"]:
        stat = Station(stations.loc[counter,])
        cmd2 = f'SELECT * from samples WHERE station_id = {item}'
        samples = pd.read_sql(cmd2, conn)

        counter2 = 0
        for element in samples["id"]:
            samp = Sample(samples.loc[counter2,])

            cmd3 = f"""SELECT * from gasex_sfg INNER JOIN sfg on sfg.name = gasex_sfg.name WHERE sample_id = {element}"""
            samp.sfg = pd.read_sql(cmd3, conn)

            cmd4 = f"""SELECT * from gasex_surftens WHERE sample_id = {element}"""
            samp.tension = pd.read_sql(cmd4, conn)

            cmd5 = cmd3.replace("sfg", "lt")
            samp.lt = pd.read_sql(cmd5, conn)
            stat.samples.append(samp)
            counter2 += 1

        out.append(stat)
        counter += 1

    return out

def get_dppc_data():
    pass

q = get_stations()

for stat in q:
    for sa in stat.samples:
        sa.make_lt()
        print(sa.data["type"])
