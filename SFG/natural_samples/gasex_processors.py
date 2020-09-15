import itertools as ito
from string import Template
from typing import List, Dict, Tuple

import numpy as np
import pandas as pd
from scipy import stats
from specsnake.sfg_spectrum import AverageSpectrum, SfgAverager

from SFG.orm.gasex_dtos import GasexSamples, GasexLt, GasexStations
from SFG.orm.interact import DbInteractor


class SampleProcessor:
    """This class has the purpose to map a list of sample objects to a list of property dictionaries that
    might be converted to a pandas dataframe for data analysis."""

    # todo: add the salinity correction to the surface tension!

    def __init__(self, samples: List[GasexSamples], interactor: DbInteractor):
        self.samples = samples
        self.ia = interactor

    def convert_to_dict(self, sample: GasexSamples):
        # create the basic dictionary
        temp = sample.to_basic_dict()

        # add the surface tension
        tension = SampleProcessor.get_corrected_salinity(sample)
        temp["surface_tension"] = tension if tension is not None else None

        # add the surface pressure of the first measured LT
        first_measured: List[GasexLt] = SampleProcessor.sort_by_measured_time([s.lt for s in sample.lt])
        if len(first_measured) > 0:
            lt_spec_object = DbInteractor.construct_lt(first_measured[0])
            temp_surface_pressure = lt_spec_object.get_maximum_pressure()
            # remove obvious outliers with a surface pressure bigger than pure water
            temp["max_surface_pressure"] = temp_surface_pressure if temp_surface_pressure < 72 else None
            temp["lift_off_compression_ratio"] = (first_measured[0].lift_off.lift_off) / np.max(lt_spec_object.area) if \
                first_measured[0].lift_off is not None else None
        else:
            temp["max_surface_pressure"] = None
            temp["lift_off_compression_ratio"] = None

        # SFG coverage
        temp["coverage"] = self.ia.get_coverage(sample.sfg.sfg) if sample.sfg is not None else None
        return temp

    def get_list_of_sample_dicts(self):
        return list(map(self.convert_to_dict, self.samples))

    @staticmethod
    def get_corrected_salinity(sample: GasexSamples):
        temp = sample.tension.surface_tension if sample.tension is not None else None
        if temp is not None:
            # first step: correct the surface tension value for the difference between calibration (20 °C)
            # and real lab temperature(21 °C)
            raw = float(sample.tension.surface_tension) * 0.99703

            # second step: correct for the salinity according to the type of the sample (surface or bulk)
            if sample.type == "deep":
                raw += SampleProcessor.correct_salinity(float(sample.station.station_plan.salinity_depth))
            else:
                raw += SampleProcessor.correct_salinity(float(sample.station.station_plan.salinity_surface))
            return round(raw, 2)

    @staticmethod
    def sort_by_measured_time(spec_list):
        return sorted(spec_list, key=lambda x: x.measured_time)

    @staticmethod
    def map_samples_to_category(samples: List[GasexSamples]) -> Dict[str, List[GasexSamples]]:
        """Match the samples to their corresponding types."""
        category_map = {category: list(group) for category, group in ito.groupby(samples, lambda x: x.type)}
        category_map["sml"] = [i for i in samples if i.type in ("p", "s")]
        return category_map

    @staticmethod
    def correct_salinity(salinity: float) -> float:
        """A function yielding a salinity-dependent correction factor for the surface tension. The reference salinity
        where the factor equals zero is 17 PSU."""
        return 0.52552 - 0.0391 * salinity


class StationProcessor:
    """This class has the purpose to map a list of station objects to a list of property dictionaries that
    might be converted to a pandas dataframe for data analysis. Note that this class has the ability to average multiple
    SFG spectra and other measurements to get an average value for plate, screen and SML, bulk samples."""

    def __init__(self, stations: List[GasexStations], interactor: DbInteractor):
        self.stations = stations
        self.interactor = interactor
        self.sample_processor = SampleProcessor([], self.interactor)

    def convert_to_dict(self, station: GasexStations):
        temp = station.to_basic_dict()
        # categories to check: plate, screen, deep, SML
        #  properties to get: coverage, max.pressure, lift_off, tension
        samples_to_categories = SampleProcessor.map_samples_to_category(station.samples)

        for cat in samples_to_categories:
            if len(samples_to_categories[cat]) > 0:
                value_dict_list = self.get_values_from_sample_list(samples_to_categories[cat])
                for value_dict in value_dict_list:
                    value = value_dict["value"]
                    for key in value_dict:
                        if key != "value":
                            temp[f'{cat}_{value}_{key}'] = value_dict[key]

            average_spec, coverage = self.average_sfg_spectra(samples_to_categories[cat])
            temp["coverage_averaged"] = coverage

        return temp

    def get_values_from_sample_list(self, samples: List[GasexSamples]):
        """Calculate the required values (coverage, tension, surface pressure and lift_off-point) from a
        given list of samples."""
        self.sample_processor.samples = samples
        value_dict_list = pd.DataFrame(self.sample_processor.get_list_of_sample_dicts())

        values = ["surface_tension", "max_surface_pressure", "lift_off_compression_ratio", "coverage"]
        return [StationProcessor.average_value(value_dict_list, i) for i in values]

    def average_sfg_spectra(self, samples: List[GasexSamples]) -> Tuple[AverageSpectrum, float]:
        """Averages the SFG spectra belonging to the provided list of samples. Returns the average spectrum
        and the calculated coverage."""
        spectra = list(
            map(self.interactor.construct_sfg, [i.sfg.sfg for i in samples if (i is not None and i.sfg is not None)]))

        # calculate the average spectrum, enforce default wavenumber scaling and baseline correction
        averager = SfgAverager(spectra, enforce_scale=True, baseline=True,
                               references=self.interactor.get_reference_integrals())

        # get the coverage
        coverage = averager.coverage
        return averager.average_spectrum, coverage

    def get_list_of_station_dicts(self):
        """Convenience function applying the dictionary calculation method on each station."""
        return list(map(self.convert_to_dict, self.stations))

    def get_station_data_frame(self) -> pd.DataFrame:
        """Generate the Dataframe containing all the station-averaged values."""
        return pd.DataFrame(self.get_list_of_station_dicts())

    @staticmethod
    def average_value(df: pd.DataFrame, value: str):
        temp = df[value].astype('float64').describe()[["count", "mean", "std"]].to_dict()
        temp["value"] = value
        return temp


class SamplePlotProcessor:

    def __init__(self, df: pd.DataFrame):
        self.df = df
        self.tension_unit = "mNm<sup>-1</sup>"
        self.properties = {
            "coverage": "surface coverage",
            "surface_tension": f'surface tension/ {self.tension_unit}',
            "max_surface_pressure": f'max. surface pressure/ {self.tension_unit}',
            "lift_off_compression_ratio": "lift-off compression ratio",
        }
        self.label_map = {
            "s": "screen",
            "p": "plate",
            "a": "Alkor",
            "r": "zodiac",
            "sml": "SML",
            "deep": "bulk"
        }

        self.stat_titles = {
            ""
        }

    def split_dataset(self, category: str, variants: Tuple[str, str]):
        """Splits the dataframe by the given categories into two."""
        return (*(self.filter_by_category(category, i) for i in variants), *variants)

    def filter_by_category(self, category: str, value: str) -> pd.DataFrame:
        """A convenience function to filter for a specific column value"""
        return self.df[self.df[category] == value]

    def gather_sml(self) -> pd.Series:
        """Convert the plate and screen tag to SML in order to put them all together."""
        return self.df["type"].apply(lambda x: "sml" if x in ('s', 'p') else x)

    def unwrap_properties_for_plotting(self, df1, df2, label_1, label_2):
        """Prepare a dictionary that is suitable for the low-level plotly plotting API."""
        dicts_for_plotly_trace = {}
        for prop in self.properties:
            temp1 = self.get_plotly_dict_from_dataframe(df1, prop, label_1)
            temp2 = self.get_plotly_dict_from_dataframe(df2, prop, label_2)
            dicts_for_plotly_trace[self.properties[prop]] = sorted([*temp1, *temp2], key=lambda x: x['name'])

        return dicts_for_plotly_trace

    def get_plotly_dict_from_dataframe(self, df: pd.DataFrame, _property: str, label: str):
        cruise_1 = df[df['cruise'] == 1]
        cruise_2 = df[df['cruise'] == 2]
        return ({'name': f'C1, {self.label_map[label]}', 'y': cruise_1[_property]},
                {'name': f'C2, {self.label_map[label]}', 'y': cruise_2[_property]})


class SampleLatexProcessor:
    table_template = Template(
        """\centering\n \caption[SFG spectral regionss]{Spectral regions probed during VSFG 
        measurements of the SML samples}\n\label{tab:sfgreg}\n$tabular\end{table}""")

    def render_tabular(self, values):
        table = """\\begin{tabular}{cccc}\n
        \\hline\n 
        category & t-statistics & p-Value & equal mean\Tstrut\Bstrut \\\\ \n
        \\hline\n"""
        for value in values:
            table += f'{value["category"]} & {value["tstat"]} & {value["pval"]}  & {value["equal"]} \\\\ \n'
        return table + "\end{tabular}"

    def render_table(self, values):
        return self.table_template.substitute(tabular=self.render_tabular(values))


# helper functions for statistics

def by_category(df, category, value) -> pd.DataFrame:
    """Filter a dataframe according to a column having a specific value"""
    return df[df[category] == value]


def split_dataset(df: pd.DataFrame, category: str, variants: Tuple[str, str]):
    """Splits the dataframe by the given categories into two."""
    return tuple(by_category(df, category, i) for i in variants)


def first_split_then_cruise(df: pd.DataFrame, category: str, variants: Tuple[str, str]):
    """Splits the dataframe first according to the category and variants and each of this variants by cruise"""
    var1, var2 = split_dataset(df, category, variants)
    return {variants[0]: split_dataset(var1, 'cruise', (1, 2)), variants[1]: split_dataset(var2, 'cruise', (1, 2))}


def first_cruise_then_split(df: pd.DataFrame, category: str, variants: Tuple[str, str]):
    cruise1, cruise2 = split_dataset(df, 'cruise', (1, 2))
    return {variants[0]: split_dataset(cruise1, category, variants),
            variants[1]: split_dataset(cruise2, category, variants)}


def compare_columns_by_ttest(col1, col2):
    return stats.ttest_ind(col1, col2, equal_var=False, nan_policy='omit')


def apply_ttest_along_properties(df1: pd.DataFrame, df2: pd.DataFrame):
    """Apply the two-sided t test to all the properties for convenience"""
    properties = ('surface_tension', 'coverage', 'max_surface_pressure', 'lift_off_compression_ratio')
    return {i: compare_columns_by_ttest(df1[i], df2[i]) for i in properties}


def apply_ttest_by_category_and_cruise(df: pd.DataFrame, category: str, variants: Tuple[str, str]):
    separated = first_split_then_cruise(df, category, variants)
    print(variants[0])
    print(separated[variants[0]][0])
    print(variants[1])
    print(separated[variants[1]][1])
    cruise1 = apply_ttest_along_properties(separated[variants[0]][0], separated[variants[1]][0])
    cruise2 = apply_ttest_along_properties(separated[variants[0]][1], separated[variants[1]][1])
    return {"cruise1": cruise1, "Cruise2": cruise2}


def apply_ttest_by_cruise_and_category(df: pd.DataFrame, category: str, variants: Tuple[str, str]):
    separated = first_split_then_cruise(df, category, variants)
    variant1 = apply_ttest_along_properties(separated[variants[0]][0], separated[variants[0]][1])
    variant2 = apply_ttest_along_properties(separated[variants[1]][0], separated[variants[1]][1])

    return {variants[0]: variant1, variants[1]: variant2}


def gather_sml(df: pd.DataFrame) -> pd.Series:
    """Convert the plate and screen tag to SML in order to put them all together."""
    return df["type"].apply(lambda x: "sml" if x in ('s', 'p') else x)
