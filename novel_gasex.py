from spectrum import SfgSpectrum, LtIsotherm, DummyPlotter
import pandas as pd
import sqlite3
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats
from datetime import datetime

conn = sqlite3.connect("test.db")


class GasexManager:

    def __init__(self, database):

        self.conn = sqlite3.connect(database)
        self.dates = self.get_dppc_reference()
        self.stations, self.station_table = self.get_stations()
        self.calculate_sample_values()
        self.set_new_pd_columns()

    def get_dppc_reference(self):
        cmd = "SELECT * FROM sfg WHERE name GLOB '*DPPC_*.*' AND measured_time BETWEEN '2018-01-01' AND '2018-12-31'"
        dppc_specs = pd.read_sql(cmd, self.conn)
        dates = {}
        for i in range(len(dppc_specs)):
            name = dppc_specs.loc[i, "name"]
            time = datetime.strptime((dppc_specs.loc[i, "measured_time"]), '%Y-%m-%d %H:%M:%S')
            meta = {"name": name, "time": time}
            data = ("wavenumbers", "sfg", "ir", "vis")
            tup = map(lambda x: np.fromstring(dppc_specs.loc[i, x], sep=";"), data)
            s = SfgSpectrum(*tup, meta)
            sr = s.yield_spectral_range()

            # if sr[0] < 2800 and sr[1] > 3010 and name.split("_")[-1] != "ppp":
            if name.split("_")[-1] != "ppp":
                q = s.calculate_ch_integral(average="gernot")
                if time.date() not in dates:
                    dates[time.date()] = [q]
                else:
                    dates[time.date()].append(q)

        for item in dates:
            dates[item] = np.average(np.array(dates[item]))

        return dates

    def get_stations(self):
        cmd = "SELECT * FROM stations"
        stations = pd.read_sql(cmd, self.conn)

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

        return out, stations

    def calculate_sample_values(self):
        for station in self.stations:
            for sample in station.samples:
                sample.calc_values(self.dates)
                if sample.tension is not None:
                    if sample.data["type"] == "deep":
                        sample.tension += GasexManager.correct_salinity(station.data["deep_salinity"])
                    else:
                        sample.tension += GasexManager.correct_salinity(station.data["surface_salinity"])
                    sample.tension = round(sample.tension, 2)

            station.get_all_stats()

    def set_new_pd_columns(self):

        cols = [i for i in self.stations[0].stats]
        self.station_table = self.station_table.reindex(columns=self.station_table.columns.tolist() + cols)
        for station in self.stations:
            for stat in station.stats:
                self.station_table.loc[self.station_table.id == station.data["id"], stat] = station.stats[stat]

    def to_sql(self, name):
        pass

    def to_excel(self):
        self.station_table.to_excel("stations.xlsx")

    @staticmethod
    def correct_salinity(salinity):
        return 0.5225-0.0391*salinity


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


    def get_date(self):
        date = datetime.strptime(self.data["date"], '%Y-%m-%d').date()
        return date

    def calc_stat(self, var, var2):
        value = [getattr(i, var) for i in self.samples if i.data["type"] in var2]
        value = [float(i) for i in value if i is not None]
        n = int(len(value))
        av = np.average(value)
        std = np.std(value)
        return av, std, n

    def get_all_stats(self):

        dic1 = {"deep": ("deep",), "sml": ("p", "s"), "plate": ("p",), "screen": ("s",)}
        dic2 = {"tension": "tension", "lift": "lift_off", "max": "max_pressure", "coverage": "coverage"}
        newdic = {}

        for item in self.stats:
            temp = item.split("_")
            value = self.calc_stat(dic2[temp[1]], dic1[temp[0]])
            newdic[item] = value[0]
            newdic[item + "_std"] = value[1]
            newdic[item + "_n"] = value[2]

        self.stats = newdic


class Sample:

    def __init__(self, data):
        self.data = data
        self.sfg = []
        self.lt = []

        self.tension = None
        self.max_pressure = None
        self.lift_off = None
        self.coverage = None

        self.sfg_spectra = []
        self.lt_isotherms = []

        self.get_tension()

    # todo: invoke the calculations necessary for the station

    def get_tension(self):
        if self.tension is None:
            pass

        elif len(self.tension) > 1:
            raise ValueError("Invalid number of surface tension values!")

        elif len(self.tension) == 0:
            self.tension = None

        else:
            self.tension = float(self.tension.loc[0, "surface_tension"])

    def make_sfg(self):

        if len(self.sfg) > 1:
            raise ValueError(f'Invalid number of SFG spectra for {self.sfg["name"]}!')

        else:
            if len(self.sfg["name"]) > 0:

                name = (self.sfg["name"].values[0][0])
                time = datetime.strptime((self.sfg.loc[0, "measured_time"]), '%Y-%m-%d %H:%M:%S')
                meta = {"name": name, "time": time}

                # pass a tuple of column names to the map function, fetching the corresponding data from the dataframe
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
            # pass a tuple of column names to the map function, fetching the corresponding data from the dataframe
            # and convert them to numpy arrays. Pass them to the LtIsotherm constructor and unpack them
            data = ("time", "area", "apm", "surface_pressure")
            tup = map(lambda x: np.fromstring(self.lt.loc[i, x], sep=";"), data)
            l = LtIsotherm(name, time, *tup, lift_off=lift_off)
            self.lt_isotherms.append(l)

    def calc_pressure(self):

        try:
            temp = sorted(self.lt_isotherms, key=lambda x: x.measured_time)
            pressure = round(temp[0].get_maximum_pressure(), 1)
            self.max_pressure = pressure
        except IndexError:
            pass

    def calc_coverage(self, dates):

        try:
            factor = dates[self.sfg_spectra[0].meta["time"].date()]
            integral = self.sfg_spectra[0].calculate_ch_integral(average="gernot")
            if integral < 0:
                integral = 0
            self.coverage = round(np.sqrt(integral / factor), 4)

        except IndexError:
            pass

    def get_lift_off(self):
        try:
            temp = sorted(self.lt_isotherms, key=lambda x: x.measured_time)
            # divide the raw lift-off point by the start area of compression
            self.lift_off = round(temp[0].lift_off/(np.max(temp[0].area)), 3)

        except (IndexError, TypeError):
            pass

    def calc_values(self, dates):
        self.get_tension()
        self.make_lt()
        self.make_sfg()

        self.calc_pressure()
        self.calc_coverage(dates)
        self.get_lift_off()

G = GasexManager("test.db")

p1 = (G.station_table["sml_lift_off"].corr(G.station_table["sml_max_pressure"]))
s1 = (G.station_table["sml_lift_off"].corr(G.station_table["sml_max_pressure"], method="spearman"))

p2 = (G.station_table["sml_coverage"].corr(G.station_table["sml_max_pressure"]))
s2 = (G.station_table["sml_coverage"].corr(G.station_table["sml_max_pressure"], method="spearman"))

p3 = (G.station_table["sml_coverage"].corr(G.station_table["sml_tension"]))
s3 = (G.station_table["sml_coverage"].corr(G.station_table["sml_tension"], method="spearman"))

p4 = (G.station_table["sml_coverage"].corr(G.station_table["sml_lift_off"]))
s4 = (G.station_table["sml_coverage"].corr(G.station_table["sml_lift_off"], method="spearman"))

p5 = (G.station_table["sml_tension"].corr(G.station_table["sml_lift_off"]))
s5 = (G.station_table["sml_tension"].corr(G.station_table["sml_lift_off"], method="spearman"))

p6 = (G.station_table["sml_tension"].corr(G.station_table["sml_max_pressure"]))
s6 = (G.station_table["sml_tension"].corr(G.station_table["sml_max_pressure"], method="spearman"))


out1 = f'The correlation between lift_off and max_pressure is {p1} (linear) and {s1} (nonlinear)'
out2 = f'The correlation between coverage and max_pressure is {p2} (linear) and {s2} (nonlinear)'
out3 = f'The correlation between coverage and tension is {p3} (linear) and {s3} (nonlinear)'
out4 = f'The correlation between coverage and lift_off is {p4} (linear) and {s4} (nonlinear)'
out5 = f'The correlation between tension and lift_off is {p5} (linear) and {s5} (nonlinear)'
out6 = f'The correlation between tension and max_pressure is {p6} (linear) and {s6} (nonlinear)'


pst = scipy.stats.ttest_ind(G.station_table["plate_tension"].values, G.station_table["screen_tension"].values,
                             equal_var=False, nan_policy='omit')

psc = scipy.stats.ttest_ind(G.station_table["plate_coverage"].values, G.station_table["screen_coverage"].values,
                             equal_var=False, nan_policy='omit')

psl = scipy.stats.ttest_ind(G.station_table["plate_lift_off"].values, G.station_table["screen_lift_off"].values,
                             equal_var=False, nan_policy='omit')

psp = scipy.stats.ttest_ind(G.station_table["plate_max_pressure"].values, G.station_table["screen_max_pressure"].values,
                             equal_var=False, nan_policy='omit')

sdt = scipy.stats.ttest_ind(G.station_table["deep_tension"].values, G.station_table["sml_tension"].values,
                             equal_var=False, nan_policy='omit')

sdc = scipy.stats.ttest_ind(G.station_table["deep_coverage"].values, G.station_table["sml_coverage"].values,
                             equal_var=False, nan_policy='omit')

sdl = scipy.stats.ttest_ind(G.station_table["deep_lift_off"].values, G.station_table["sml_lift_off"].values,
                             equal_var=False, nan_policy='omit')

sdp = scipy.stats.ttest_ind(G.station_table["deep_max_pressure"].values, G.station_table["sml_max_pressure"].values,
                             equal_var=False, nan_policy='omit')

print(f'plate screen tension: {pst}')
print(f'plate screen coverage: {psc}')
print(f'plate screen lift_off: {psl}')
print(f'plate screen pressure: {psp}')

print(f'sml deep tension: {sdt}')
print(f'sml deep coverage: {sdc}')
print(f'sml deep lift_off: {sdl}')
print(f'sml deep pressure: {sdp}')

print(scipy.stats.ttest_ind(G.station_table["deep_lift_off"].values, G.station_table["sml_lift_off"].values,
                             equal_var=False, nan_policy='omit'))



G.station_table["date"] = pd.to_datetime(G.station_table["date"])
mask1 = (G.station_table['date'] > '2018-6-1') & (G.station_table['date'] <= '2018-7-1')
cruise1 = G.station_table.loc[mask1]

mask2 = (G.station_table['date'] > '2018-9-1') & (G.station_table['date'] <= '2018-10-1')
cruise2 = G.station_table.loc[mask2]

cct = scipy.stats.ttest_ind(cruise1["sml_tension"].values, cruise2["sml_tension"].values,
                             equal_var=False, nan_policy='omit')
ccc = scipy.stats.ttest_ind(cruise1["sml_coverage"].values, cruise2["sml_coverage"].values,
                             equal_var=False, nan_policy='omit')
ccl = scipy.stats.ttest_ind(cruise1["sml_lift_off"].values, cruise2["sml_lift_off"].values,
                             equal_var=False, nan_policy='omit')
ccp = scipy.stats.ttest_ind(cruise1["sml_max_pressure"].values, cruise2["sml_max_pressure"].values,
                             equal_var=False, nan_policy='omit')

dcct = scipy.stats.ttest_ind(cruise1["deep_tension"].values, cruise2["deep_tension"].values,
                             equal_var=False, nan_policy='omit')
dccc = scipy.stats.ttest_ind(cruise1["deep_coverage"].values, cruise2["deep_coverage"].values,
                             equal_var=False, nan_policy='omit')
dccl = scipy.stats.ttest_ind(cruise1["deep_lift_off"].values, cruise2["deep_lift_off"].values,
                             equal_var=False, nan_policy='omit')
dccp = scipy.stats.ttest_ind(cruise1["deep_max_pressure"].values, cruise2["deep_max_pressure"].values,
                             equal_var=False, nan_policy='omit')

print(f'sml tension c1c2: {cct}')
print(f'sml coverage c1c2: {ccc}')
print(f'sml lift_off c1c2: {ccl}')
print(f'sml pressure c1c2: {ccp}')
print(f'deep tension c1c2: {dcct}')
print(f'deep coverage c1c2: {dccc}')
print(f'deep lift_off c1c2: {dccl}')
print(f'deep pressure c1c2: {dccp}')

small = G.station_table[G.station_table["type"] == "small"]
big = G.station_table[G.station_table["type"] == "big"]

bst = scipy.stats.ttest_ind(big["sml_tension"].values, small["sml_tension"].values,
                             equal_var=False, nan_policy='omit')

bsc = scipy.stats.ttest_ind(big["sml_coverage"].values, small["sml_coverage"].values,
                             equal_var=False, nan_policy='omit')

bsl = scipy.stats.ttest_ind(big["sml_lift_off"].values, small["sml_lift_off"].values,
                             equal_var=False, nan_policy='omit')

bsp = scipy.stats.ttest_ind(big["sml_max_pressure"].values, small["sml_max_pressure"].values,
                             equal_var=False, nan_policy='omit')

print(f'sml tension bs: {bst}')
print(f'sml coverage bs: {bsc}')
print(f'sml lift_off bs: {bsl}')
print(f'sml pressure bs: {bsp}')

#G.station_table.to_excel("gasex_data.xlsx", na_rep="-")
G.station_table.plot.scatter(x="deep_tension", y="sml_tension")
plt.show()