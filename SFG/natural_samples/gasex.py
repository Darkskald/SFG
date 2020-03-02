from SFG.orm.orm import WorkDatabaseWizard

import numpy as np
import pandas as pd
from scipy import stats
from typing import Tuple


class GasExManager:
    """The GasExManager classes operates on the data recorded during and in connection with the Baltic GasEx
    cruise in 2018. If instantiated with new=True, Samples and Stations classes are created from the raw data
    stored in the database, connected to each other, values calculated and persisted. For everyday-use, the
    new=False option is recommended. In this case, a Pandas dataframe containing the data of all Stations
    is generated. This is suitable for plotting and data analysis."""

    def __init__(self, new=True):
        self.wdw = WorkDatabaseWizard()

        if new:
            self.stations = self.get_stations()
            self.calculate_sample_values()
            self.persist_stations()

        self.station_table = self.generate_station_table()
        self.station_table["date"] = pd.to_datetime(self.station_table["date"])

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
                        sample.tension += GasExManager.correct_salinity(float(station.data.deep_salinity))
                    else:
                        sample.raw_tension = sample.tension
                        sample.tension += GasExManager.correct_salinity(float(station.data.surface_salinity))

                    sample.tension = round(sample.tension, 2)
                sample.set_values()
                self.wdw.session.commit()

            station.get_all_stats()
            station.reference_tension_temperature()

    def persist_stations(self):
        """Invokes the persist_stat function of each station, ensuring that the calculated values are written to
        the station_stats table in the database."""

        for station in self.stations:
            station.persist_stats(self.wdw)

    @staticmethod
    def correct_salinity(salinity):
        """A function yielding a salinity-dependent correction factor for the surface tension. The reference salinity where
        the factor equals zero is 17 PSU."""
        return 0.52552 - 0.0391 * salinity

    def generate_station_table(self):
        """Generates a pandas dataframe from the joint columns of stations and station_stats table."""
        command = """
        SELECT *
        FROM stations
        INNER JOIN station_stats
        on stations.id = station_stats.station_id;
        """
        temp = pd.read_sql(command, self.wdw.session.bind)
        return temp


class Sample:
    """The sample class unites all measurement data corresponding to a specific sample conveniently. It is used
    by the Station class to calculate the corresponding values"""
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

    def __str__(self):
        temp = f'\nhash: {self.data.sample_hash}\nnumber_sfg: {len(self.sfg_spectra)}\nnumber_lt: {len(self.lt_isotherms)}'
        return temp

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
            with open("errorlog.txt", "w") as outfile:
                outfile.write(str(self))

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

    def set_values(self):
        self.data.surface_tension = self.tension
        self.data.max_pressure = self.max_pressure
        self.data.lift_off = self.lift_off
        self.data.coverage = self.coverage


class Station:
    """This class represents a cruise station. It takes care about its samples, calculates average values and
    makes them available for plotting purposes via its station table (a Pandas dataframe)"""
    def __init__(self, data):
        self.data = data
        self.doy = self.data.date.timetuple().tm_yday + (1 - self.data.number) * 0.2
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
        stats.doy = self.doy
        wizard.session.add(stats)
        wizard.session.commit()

    def reference_tension_temperature(self):
        try:
            temp1_hi = Station.calc_sal_tension(float(self.data.surface_temperature), float(self.data.surface_salinity))
            temp1_lo = Station.calc_sal_tension(21, float(self.data.surface_salinity))
            temp1 = self.stats["sml_rawtension"] + (temp1_hi - temp1_lo)
            self.stats["sml_rawtension"] = round(temp1, 2)
        except:
            pass

        try:
            temp2_hi = Station.calc_sal_tension(float(self.data.deep_temperature), float(self.data.deep_salinity))
            temp2_lo = Station.calc_sal_tension(21, float(self.data.deep_salinity))
            temp2 = self.stats["deep_rawtension"] + (temp2_hi - temp2_lo)
            self.stats["deep_rawtension"] = round(temp2, 2)
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


class GasExWorkDatabaseWizard(WorkDatabaseWizard):

    def __init__(self):
        super().__init__()

    def load_station_data(self) -> pd.DataFrame:
        """Extracts all station data as Pandas dataframe"""

        cmd = """
            SELECT 
            
            date AS 'sampling date',
            doy AS 'day of the year',
            
            plate_coverage AS 'surface coverage plate',
            plate_coverage_std,
            plate_tension AS 'surface tension plate',
            plate_tension_std,
            plate_lift_off AS 'lift-off compression ratio plate',
            plate_lift_off_std,
            plate_max_pressure AS 'max. surface pressure plate',
            plate_max_pressure_std, 
            
            screen_coverage AS 'surface coverage screen',
            screen_coverage_std,
            screen_tension AS 'surface tension screen',
            screen_tension_std,
            screen_lift_off AS 'lift-off compression ratio screen',
            screen_lift_off_std,
            screen_max_pressure AS 'max. surface pressure screen',
            screen_max_pressure_std, 
            
            sml_coverage AS 'surface coverage SML',
            sml_coverage_std,
            sml_tension AS 'surface tension SML',
            sml_tension_std,
            sml_lift_off AS 'lift-off compression ratio SML',
            sml_lift_off_std,
            sml_max_pressure AS 'max. surface pressure SML',
            sml_max_pressure_std,
            
            
            deep_coverage AS 'surface coverage bulk',
            deep_coverage_std,
            deep_tension AS 'surface tension bulk',
            deep_tension_std,
            deep_lift_off AS 'lift-off compression ratio bulk',
            deep_lift_off_std,
            deep_max_pressure AS 'max. surface pressure bulk',
            deep_max_pressure_std 
            
            FROM stations
            INNER JOIN station_stats
            on stations.id = station_stats.station_id;
            """
        return GasExWorkDatabaseWizard.yield_df_from_sql(cmd, self.engine)

    def load_samples_by_type(self, _type) -> pd.DataFrame:

        base_cmd = """SELECT
                    sample_hash AS 'name', 
                    coverage AS 'surface coverage',
                    max_pressure AS 'max. surface pressure',
                    lift_off AS 'lift-off compression ratio',
                    surface_tension AS 'surface tension'
                    FROM samples
                    """
        end_cmd = {
            "all": ";",
            "sml": """WHERE type='p' OR type='s';""",
            "bulk": """WHERE type='deep';""",
            "p": """WHERE type='p';""",
            "s": """WHERE type='s';""",
            "alkor": """WHERE location='a';""",
            "rubber": """WHERE location='r';""",
        }
        cmd = f'{base_cmd}{end_cmd[_type]}'
        return GasExWorkDatabaseWizard.yield_df_from_sql(cmd, self.engine)

    @staticmethod
    def yield_df_from_sql(sql_statement: str, db) -> pd.DataFrame:
        """Auxiliary function to load data from an SQL database yielding a dataframe"""
        return pd.read_sql(sql_statement, db)

    @staticmethod
    def apply_test(df: pd.DataFrame, column: str, test_func) -> Tuple[float, float]:
        """Function to apply a statistical test function to a column of a dataframe"""
        temp = df[column].to_numpy()
        temp = temp[np.logical_not(np.isnan(temp))]
        return test_func(temp)

if __name__ == "__main__":
    GasExManager()


# todo: make a new table with measurement days and the corresponding maximum dppc intensities
# todo: the SFG spectra have to be averaged on the station level and not on the sample level
# todo: map parameters like tension, lt_measures and coverage to samples

