from functools import partial

import numpy as np

from SFG.orm.gasex_dtos import GasexSamples, GasexStations
from SFG.orm.interact import DbInteractor


class SampleWrapper:

    # todo: mit Gernot abklären was mit den anderen Lts passieren soll
    def __init__(self, sample: GasexSamples):
        self.sample = sample

    def get_tension(self):
        try:
            # correct the surface tension value for the difference between calibration (20) and real lab temperature(21)
            raw = float(self.sample.tension.surface_tension) * 0.99703

            if self.sample.type == "deep":
                raw += SampleWrapper.correct_salinity(float(self.sample.station.station_plan.salinity_depth))
            else:
                raw += SampleWrapper.correct_salinity(float(self.sample.station.station_plan.salinity_surface))

            return round(raw, 2)

        except AttributeError:
            return None

    def get_max_pressure(self):
        try:
            lt_isotherms = list(map(DbInteractor.construct_lt, [i.lt for i in self.sample.lt]))
            temp = sorted(lt_isotherms, key=lambda x: x.measured_time)
            pressure = round(temp[0].get_maximum_pressure(), 1)
            return pressure
        except IndexError:
            return None

    def get_lift_off(self):
        temp = sorted([i.lt for i in self.sample.lt], key=lambda x: x.measured_time)
        if len(temp) > 0:
            try:
                lo = temp[0].lift_off.lift_off
                lt = DbInteractor.construct_lt(temp[0])
                return np.round(lo / (np.max(lt.area)), 3)
            except AttributeError:
                return None
        else:
            return None

    @DeprecationWarning
    def get_coverage(self, references):
        """ATTENTION: THIS IS APPARENTLY WRONG BECAUSE THE SQUARE ROOT IS MISSING"""
        reference = references[self.sample.sfg.sfg.measured_time.date()]
        integral = DbInteractor.construct_sfg(self.sample.sfg.sfg).calculate_ch_integral()
        return integral/reference


    @staticmethod
    def correct_salinity(salinity):
        """A function yielding a salinity-dependent correction factor for the surface tension. The reference salinity
        where the factor equals zero is 17 PSU."""
        return 0.52552 - 0.0391 * salinity


class StationWrapper:

    def __init__(self, station: GasexStations, references):
        self.station = station
        self.cruise = self.get_cruise()
        self.references = references

        self.sample_types = {"deep": ("deep",), "sml": ("p", "s"), "plate": ("p",), "screen": ("s",)}
        self.calculation_function = {'tension': 'get_tension', 'lift_off': 'get_lift_off',
                                     'max_pressure': 'get_max_pressure', 'coverage': 'get_coverage'}

    def select_samples(self, sample_type):
        return list(filter(lambda x: x.type in self.sample_types[sample_type], self.station.samples))

    def calculate_data(self, sample_type, target):
        func_name = self.calculation_function[target]
        selection = [SampleWrapper(i) for i in self.select_samples(sample_type)]
        selection_values = np.array(list(map(self.get_function(func_name), selection)), dtype=np.float64)
        len_ = selection_values[~np.isnan(selection_values)]

        n = len(len_)
        av = np.round(np.nanmean(selection_values), decimals=2)
        std = np.round(np.nanstd(selection_values), decimals=4)
        return av, std, n

    def get_cruise(self):
        if self.station.date.month == 6:
            return 1
        elif self.station.date.month == 9:
            return 2
        else:
            raise ValueError(f'Invalid cruise month in {self.station}')

    def get_function(self, func_name):
        if func_name != 'coverage':
            return getattr(SampleWrapper, func_name)
        else:
            return partial(getattr(SampleWrapper, func_name), references=self.references)