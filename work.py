from orm import DatabaseWizard
from spectrum import SfgSpectrum, LtIsotherm, SfgAverager

import numpy as np
import matplotlib.pyplot as plt


class WorkDatabaseWizard(DatabaseWizard):

    def __init__(self):
        super().__init__()
        self.session.commit()

    def get_dppc_references(self):
        """A function querying the sfg table for DPPC reference spectra, generating the corresponding objects,
        calculating the CH integral making use of the SfgAverager class and returning a dictionary of date objects
        with the corresponding intensities."""
        dates = {}

        q_dppc = self.session.query(self.sfg). \
            filter(self.sfg.name.op('GLOB')('*DPPC_*.*')). \
            filter(self.sfg.measured_time.between('2018-01-01', '2018-12-31')) \
            .filter(~self.sfg.name.contains('ppp'))

        for item in q_dppc:
            s = WorkDatabaseWizard.construct_sfg(item)
            _date = s.meta["time"].date()
            if _date not in dates:
                dates[_date] = [s]
            else:
                dates[_date].append(s)

        for item in dates:
            dates[item] = SfgAverager(dates[item]).integral

        # get rid of days where no DPPC spectra were recorded
        dates = {k: v for k, v in dates.items() if not np.isnan(v)}

        return dates

    # auxiliary functions

    @staticmethod
    def to_array(string):
        """Converts the raw data stored as strings back to numpy float ndarrays."""
        try:
            return np.fromstring(string, sep=",")
        except TypeError:
            return np.nan

    @staticmethod
    def construct_sfg(or_object):
        """A function constructing the SFG object from the orm declarative class."""
        meta = {"name": or_object.name, "time": or_object.measured_time}
        args = ("wavenumbers", "sfg", "vis", "ir")
        args = [WorkDatabaseWizard.to_array(getattr(or_object, i)) for i in args]
        s = SfgSpectrum(*args, meta)
        return s

    @staticmethod
    def construct_lt(or_object):
        """A function constructing the LT object from the orm declarative class."""
        args = (or_object.name, or_object.measured_time)
        add_args = ["time", "area", "apm", "surface_pressure", "lift_off"]
        add_args = [WorkDatabaseWizard.to_array(getattr(or_object, i)) for i in add_args]
        l = LtIsotherm(args[0], args[1], *add_args)
        return l


class SomeManager:

    def __init__(self):

        self.wdw = WorkDatabaseWizard()
        self.stations = self.get_stations()
        self.calculate_sample_values()
        self.persist_stations()

    def get_stations(self):
        stations = []
        or_stations = self.wdw.session.query(self.wdw.stations)
        for item in or_stations:

            station = Station(item)

            or_samples = self.wdw.session.query(self.wdw.samples) \
                .filter(self.wdw.samples.station_id == item.id)

            for or_sample in or_samples:

                sample = Sample(or_sample)

                or_sfg = self.wdw.session.query(self.wdw.sfg) \
                    .join(self.wdw.gasex_sfg, self.wdw.gasex_sfg.name == self.wdw.sfg.name). \
                    filter(self.wdw.gasex_sfg.sample_id == or_sample.id)
                for or_spec in or_sfg:
                    s = self.wdw.construct_sfg(or_spec)
                    sample.sfg_spectra.append(s)

                or_lt = self.wdw.session.query(self.wdw.lt) \
                    .join(self.wdw.gasex_sfg, self.wdw.gasex_lt.name == self.wdw.lt.name). \
                    filter(self.wdw.gasex_lt.sample_id == or_sample.id)
                for or_lt in or_lt:
                    lt = self.wdw.construct_lt(or_lt)
                    sample.lt_isotherms.append(lt)

                or_tension = self.wdw.session.query(self.wdw.gasex_surftens) \
                    .filter(self.wdw.gasex_surftens.sample_id == or_sample.id). \
                    one_or_none()
                if or_tension is not None:
                    tension = float(or_tension.surface_tension)
                    sample.tension = tension

                sample.calc_values(self.wdw.get_dppc_references())
                station.samples.append(sample)

            stations.append(station)
        return stations

    def calculate_sample_values(self):
        """Convenience function to call the stations (and, therefore, samples) functions and collect the average data
        from the raw single measurements."""
        for station in self.stations:
            for sample in station.samples:

                if sample.tension is not None:
                    if sample.data.type == "deep":
                        sample.raw_tension = sample.tension
                        sample.tension += SomeManager.correct_salinity(float(station.data.deep_salinity))
                    else:
                        sample.raw_tension = sample.tension
                        sample.tension += SomeManager.correct_salinity(float(station.data.surface_salinity))

                    sample.tension = round(sample.tension, 2)

            station.get_all_stats()
            station.reference_tension_temperature()

    def persist_stations(self):

        for station in self.stations:
            station.persist_stats(self.wdw)

    @staticmethod
    def correct_salinity(salinity):
        """A function yielding a salinity-dependent correction factor for the surface tension. The reference salinity where
        the factor equals zero is 17 PSU."""
        return 0.52552 - 0.0391 * salinity


class Sample:

    def __init__(self, data):

        # instead of the pandas dataframe, the or_class is used here
        self.data = data

        self.tension = None
        self.raw_tension = None
        self.max_pressure = None
        self.lift_off = None
        self.coverage = None

        self.sfg_spectra = []
        self.lt_isotherms = []

    def get_lift_off(self):
        """Extract the lift-off value from the dataframe and normalize it to the maximum initial area in order to make
        it comparable."""

        if len(self.lt_isotherms) > 0:
            temp = sorted(self.lt_isotherms, key=lambda x: x.measured_time)
            # divide the raw lift-off point by the start area of compression
            self.lift_off = np.round(temp[0].lift_off/(np.max(temp[0].area)), 3)

    def calc_coverage(self, dates):
        """Calculates the surface coverage by normalizing the CH integral to the
        DPPC reference intensity which is provided as a dictionary of dates."""

        try:
            factor = dates[self.sfg_spectra[0].meta["time"].date()]
            integral = self.sfg_spectra[0].calculate_ch_integral()
            if integral < 0:
                integral = 0
            self.coverage = np.round(np.sqrt(integral / factor), 4)

        except IndexError:
            pass

    def calc_pressure(self):
        """Calculate the maximum surface pressure for the LT isotherm measurement of the sample."""

        try:
            temp = sorted(self.lt_isotherms, key=lambda x: x.measured_time)
            pressure = round(temp[0].get_maximum_pressure(), 1)
            self.max_pressure = pressure
        except IndexError:
            pass

    def calc_values(self, dates):
        """Top level function to call all data-collecting routines. Needs a dictionary
        of DPPC reference intensities for the corresponding dates."""
        self.calc_pressure()
        self.calc_coverage(dates)
        self.get_lift_off()


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
            "deep_max_pressure": None,
            "sml_rawtension": None,
            "deep_rawtension": None,
        }
        self.doy = self.data.date.timetuple().tm_yday + (1-self.data.number)*0.2

        self.temp_tension_deep = None
        self.temp_tension_sml = None

    def calc_stat(self, value, sample_type):
        """A generic function to calculate all necessary station attributes (eg. SML coverage, deep tension..)."""
        value = [getattr(i, value) for i in self.samples if i.data.type in sample_type]
        value = [float(i) for i in value if i is not None]
        len_ = [i for i in value if np.isnan(i) == False]

        n = len(len_)
        av = np.round(np.nanmean(value), decimals=2)
        std = np.round(np.nanstd(value), decimals=4)

        return av, std, n

    def get_all_stats(self):
        """Function to collect the station's values including standard deviation and number of samples included
        in the average."""

        types = {"deep": ("deep",), "sml": ("p", "s"), "plate": ("p",), "screen": ("s",)}
        values = {"tension": "tension", "lift": "lift_off", "max": "max_pressure", "coverage": "coverage",
                  "rawtension": "raw_tension"}
        new_stats = {}

        for item in self.stats:

            temp = item.split("_")
            value = self.calc_stat(values[temp[1]], types[temp[0]])
            new_stats[item] = value[0]
            new_stats[item + "_std"] = value[1]
            new_stats[item + "_n"] = value[2]

        self.stats = new_stats

    def persist_stats(self, wizard):

        stats = wizard.station_stats()
        stats.station_id = self.data.id
        for key in self.stats:
            if self.stats[key] is not None:
                setattr(stats, key, self.stats[key])

        wizard.session.add(stats)
        wizard.session.commit()

    def reference_tension_temperature(self):
        try:
            temp1_hi = Station.calc_sal_tension(float(self.data.surface_temperature), float(self.data.surface_salinity))
            temp1_lo = Station.calc_sal_tension(21, float(self.data.surface_salinity))
            temp1 = self.stats["sml_rawtension"] + (temp1_hi - temp1_lo)
            # print(f"""t: {self.data.surface_temperature} S: {self.data.surface_salinity}
            # raw: {self.stats["sml_rawtension"]} after: {temp1}""")
            self.stats["sml_rawtension"] = temp1
        except:
            pass

        try:
            temp2_hi = Station.calc_sal_tension(float(self.data.deep_temperature), float(self.data.deep_salinity))
            temp2_lo = Station.calc_sal_tension(21, float(self.data.deep_salinity))
            temp2 = self.stats["deep_rawtension"] + (temp2_hi - temp2_lo)
            #print(f"""t: {self.data.deep_temperature} S: {self.data.deep_salinity}
            #raw: {self.stats["deep_rawtension"]} after: {temp2}""")
            self.stats["deep_rawtension"] = temp2
        except:
            pass

    # auxiliary functions
    @staticmethod
    def calc_pure_tension(t):
        temp = (1 - ((t + 273.15) / 647.096))
        out = (235.8 * temp ** 1.256) * (1 - 0.625 * temp)
        return out

    @staticmethod
    def calc_sal_tension(t, s):
        temp = Station.calc_pure_tension(t)
        factor = 1 + (3.766 * 10 ** -4) * s + (2.347 * 10 ** -6) * s * t
        return temp * factor


S = SomeManager()

#for s in S.stations:
    #print(f'sml: {s.stats["sml_tension"]} {s.stats["sml_rawtension"]} deep: {s.stats["deep_tension"]}\n {s.stats["deep_rawtension"]}')

# Benchmark the new version
# todo: compare the dictionary with DPPC references for gasex values produced by the old and the new routines
# todo: compare the output values

# todo: find a suitable name for the data attribute of samples and stations
# todo: make a new table with measurement days and the corresponding maximum dppc intensities
# todo: the DPPC spectra have to be averaged on the station level and not on the sample level

