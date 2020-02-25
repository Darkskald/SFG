from SFG.spectrum.spectrum import SfgSpectrum, LtIsotherm, SfgAverager
import pandas as pd
import sqlite3
import numpy as np


from datetime import datetime


class GasexManager:

    def __init__(self, database, dppc_flag=False):

        self.dppc_flag = dppc_flag
        self.dppc = []

        self.conn = sqlite3.connect(database)
        self.dates = self.get_dppc_reference()
        self.stations, self.station_table = self.get_stations()
        self.calculate_sample_values()
        self.set_new_pd_columns()
        self.station_table["date"] = pd.to_datetime(self.station_table["date"])

        print("Initialization of GasEx Manager successful!")

    def get_dppc_reference(self):
        """Returns a dictionary with each day of measurement and the corresponding DPPC integrals"""
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
                # q = s.calculate_ch_integral()
                if time.date() not in dates:
                    # dates[time.date()] = [q]
                    dates[time.date()] = [s]
                else:
                    # dates[time.date()].append(q)
                    dates[time.date()].append(s)

        for item in dates:
            # dates[item] = np.average(np.array(dates[item]))
            dates[item] = SfgAverager(dates[item]).integral

        dates = {k: v for k, v in dates.items() if not np.isnan(v)}

        return dates

    def get_stations(self):
        """Get the station data from the SQL database and transform them to the corresponding objects and a suitable
        Pandas dataframe which is used for the further data processing"""
        cmd = "SELECT * FROM stations"
        stations = pd.read_sql(cmd, self.conn)

        out = []
        counter = 0

        for item in stations["id"]:
            stat = Station(stations.loc[counter,])
            cmd2 = f'SELECT * from samples WHERE station_id = {item}'
            samples = pd.read_sql(cmd2, self.conn)

            counter2 = 0
            for element in samples["id"]:
                samp = Sample(samples.loc[counter2,])

                cmd3 = f"""SELECT * from gasex_sfg INNER JOIN sfg on sfg.name = gasex_sfg.name WHERE sample_id = {element}"""
                samp.sfg = pd.read_sql(cmd3, self.conn)

                cmd4 = f"""SELECT * from gasex_surftens WHERE sample_id = {element}"""
                samp.tension = pd.read_sql(cmd4, self.conn)

                cmd5 = cmd3.replace("sfg", "lt")
                samp.lt = pd.read_sql(cmd5, self.conn)
                stat.samples.append(samp)
                counter2 += 1

            out.append(stat)
            counter += 1

        return out, stations

    def calculate_sample_values(self):
        """Convenience function to call the stations (and, therefore, samples) functions and collect the average data
        from the raw single measurements."""
        for station in self.stations:
            for sample in station.samples:
                #important: info about DPPC is passed here
                sample.calc_values(self.dates)
                
                if sample.tension is not None:
                    if sample.data["type"] == "deep":
                        sample.tension += GasexManager.correct_salinity(station.data["deep_salinity"])
                    else:
                        sample.tension += GasexManager.correct_salinity(station.data["surface_salinity"])
                    sample.tension = round(sample.tension, 2)

            station.get_all_stats()
            station.calc_coverage(self.dates)

    def set_new_pd_columns(self):
        """Add the calculated per-station average values to the Pandas station table dataframe."""

        cols = [i for i in self.stations[0].stats]
        cols.append("doy")
        self.station_table = self.station_table.reindex(columns=self.station_table.columns.tolist() + cols)
        for station in self.stations:
            self.station_table.loc[self.station_table.id == station.data["id"], "doy"] = station.doy
            for stat in station.stats:
                self.station_table.loc[self.station_table.id == station.data["id"], stat] = station.stats[stat]

    def to_sql(self, name):
        """Futute function to export the final station table to SQL"""
        raise NotImplementedError

    def to_excel(self):
        """Converts the final station table to an Excel file."""
        self.station_table.to_excel("stations.xlsx")

    def fetch_dppc_spec(self, spec):
        """A function which returns the suitable SFG spectrum object for the reference measurement of the input spectrum
        spec."""
        for item in self.dppc:
            if item.meta["time"].date() == spec.meta["time"].date():
                return spec, item
        raise ValueError("No matching spectrum found!")

    @staticmethod
    def correct_salinity(salinity):
        """A function yielding a salinity-dependent correction factor for the surface tension. The reference salinity where
        the factor equals zero is 17 PSU."""
        return 0.5225-0.0391*salinity


class Cruise:
    """The class taking care of statistics between the two cruises"""
    def __init__(self, station_frame):
        # calculate sml, deep, screen, plate averages and standard deviations
        pass


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
        self.doy = self.get_date().timetuple().tm_yday + (1-self.data["number"])*0.2

    def get_date(self):
        """Generate a datetime.date object from the raw date string."""
        date = datetime.strptime(self.data["date"], '%Y-%m-%d').date()
        return date

    def calc_stat(self, var, var2):
        """A generic function to calculate all necessary station attributes (eg. SML coverage, deep tension..)."""
        value = [getattr(i, var) for i in self.samples if i.data["type"] in var2]
        value = [float(i) for i in value if i is not None]
        n = int(len(value))
        av = np.nanmean(value)
        std = np.nanstd(value)

        return av, std, n

    def get_all_stats(self):
        """Function to collect the station's values including standard deviation and number of samples included
        in the average."""

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

    def calc_coverage(self, dppc_ref):

        average_coverages = {"sml": {"spec": None, "coverage": None},
                             "plate": {"spec": None, "coverage": None},
                             "screen": {"spec": None, "coverage": None},
                             "deep": {"spec": None, "coverage": None}
                             }

        sml_specs = []
        deep_specs = []
        plate_specs = []
        screen_specs = []
        # todo: append to the lists above based on the sample properties
        # todo: calculate coverages for all of this values
        # todo: new dictionary value in the Station's dic with the new attributes (av_sml_coverage)
        # todo: apply the averaging procedure to DPPC as well

        for samp in self.samples:
            if len(samp.sfg_spectra) > 0:
                if samp.data["type"] == "p":
                    plate_specs.append(samp.sfg_spectra[0])
                    sml_specs.append(samp.sfg_spectra[0])

                elif samp.data["type"] == "s":
                    screen_specs.append(samp.sfg_spectra[0])
                    sml_specs.append(samp.sfg_spectra[0])

                elif samp.data["type"] == "deep":
                    deep_specs.append(samp.sfg_spectra[0])

        #print(len(deep_specs))
        self.av_sml_coverage = SfgAverager(sml_specs, references=dppc_ref).coverage
        self.av_deep_coverage = SfgAverager(deep_specs, references=dppc_ref).integral
        #print(self.av_deep_coverage)
        self.av_plate_coverage = SfgAverager(plate_specs, references=dppc_ref).coverage
        self.av_screen_coverage = SfgAverager(screen_specs, references=dppc_ref).coverage

        av_sml_coverage = SfgAverager(sml_specs, references=dppc_ref)
        av_deep_coverage = SfgAverager(deep_specs, references=dppc_ref)
        av_plate_coverage = SfgAverager(plate_specs, references=dppc_ref)
        av_screen_coverage = SfgAverager(screen_specs, references=dppc_ref)


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

    def get_tension(self):
        """Calculates the sample's surface tension. Sets it to None if no tension measurement is available for
        this sample."""
        if self.tension is None:
            pass

        elif len(self.tension) > 1:
            raise ValueError("Invalid number of surface tension values!")

        elif len(self.tension) == 0:
            self.tension = None

        else:
            self.tension = float(self.tension.loc[0, "surface_tension"])

    def make_sfg(self):
        """Generate a list of SFG objects from the dataframe stored in self.sfg"""

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
        """Generate a list of LTbjects from the dataframe stored in self.lt"""

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
        """Calculate the maximum surface pressure for the LT isotherm measurement of the sample."""

        try:
            temp = sorted(self.lt_isotherms, key=lambda x: x.measured_time)
            pressure = round(temp[0].get_maximum_pressure(), 1)
            self.max_pressure = pressure
        except IndexError:
            pass

    def calc_coverage(self, dates):

        try:
            factor = dates[self.sfg_spectra[0].meta["time"].date()]
            integral = self.sfg_spectra[0].calculate_ch_integral()
            if integral < 0:
                integral = 0
            self.coverage = round(np.sqrt(integral / factor), 4)


        except IndexError:
            pass

    def get_lift_off(self):
        """Extract the lift-off value from the dataframe and normalize it to the maximum initial area in order to make
        it comparable."""
        try:
            temp = sorted(self.lt_isotherms, key=lambda x: x.measured_time)
            # divide the raw lift-off point by the start area of compression
            self.lift_off = round(temp[0].lift_off/(np.max(temp[0].area)), 3)

        except (IndexError, TypeError):
            pass

    def calc_values(self, dates):
        """Top level function to call all data-collecting routines."""
        self.get_tension()
        self.make_lt()
        self.make_sfg()

        self.calc_pressure()
        self.calc_coverage(dates)
        self.get_lift_off()

