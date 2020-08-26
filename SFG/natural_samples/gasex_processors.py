import itertools as ito
from typing import List, Dict, Tuple

import numpy as np
import pandas as pd

from SFG.orm.gasex_dtos import GasexSamples, GasexLt, GasexStations
from SFG.orm.interact import DbInteractor
from SFG.spectrum.sfg_spectrum import SfgAverager, AverageSpectrum


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
            temp["max_surface_pressure"] = lt_spec_object.get_maximum_pressure()
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

    def split_dataset(self, category: str, variants: Tuple[str, str]):
        """Splits the dataframe by the given categories into two."""
        return (self.filter_by_category(category, i) for i in variants)

    def filter_by_category(self, category: str, value: str) -> pd.DataFrame:
        """A convenience function to filter for a specific column value"""
        return self.df[self.df[category] == value]

    def gather_sml(self) -> pd.Series:
        """Convert the plate and screen tag to SML in order to put them all together."""
        return self.df["type"].apply(lambda x: "sml" if x in ('s', 'p') else x)

    def unwrap_properties_for_plotting(self, df1, df2):
        """Prepare a dictionary that is suitable for the low-level plotly plotting API."""
        dicts_for_plotly_trace = {}
        row_index = 1
        col_index = 1
        for prop in self.properties:
            temp = [SamplePlotProcessor.get_plotly_dict_from_dataframe(i, prop, row_index, col_index) for i in
                    (df1, df2)]
            dicts_for_plotly_trace[self.properties[prop]] = temp
            row_index += 1
            col_index += 1
        return dicts_for_plotly_trace

    @staticmethod
    def get_plotly_dict_from_dataframe(df: pd.DataFrame, _property: str, row: int, col: int):
        return {'x': df['cruise'], 'y': df[_property], "row": row, "col": col}
