import pandas as pd
import sqlite3
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime

conn = sqlite3.connect("test.db")


class Station:

    def __init__(self, data):
        self.data = data
        self.samples = []

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

        self.sfg_spectra = None
        self.lt_isotherms = None

    # todo: generate the objects from the dataframes
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
            raise ValueError("Invalid number of SFG spectra!")

        else:
            # inject metainfo and values to sfg spec
            pass

    def make_lt(self):
        pass



# todo: encapsulate SQL interaction to class --> replace session controll manager
# todo: generate sfg and lt objects from sql query results, be careful about creation time!
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

        out.append(stat)
        counter += 1

    return out


q = get_stations()


# testcode section
for stat in q:
    for sample in stat.samples:
        print(sample.get_tension())
        print(sample.make_sfg())