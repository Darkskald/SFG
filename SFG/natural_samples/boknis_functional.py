"""
tasks:

final_value <- relate_to_reference <- get_reference <- get_value

this can be implemented in a class conveniently!

attributes:

- input_sfg (this may be a list)
- property_name
- reference_relator
- reference_calculator
- reference_getter
- value_getter

methods:

- compose_calculator
- calculate value
- output values in a suitable fo
"""
from __future__ import annotations
from typing import List, Callable

import numpy as np
import pandas as pd

from SFG.orm.boknis_dtos import BoknisEckSamplingDay
from SFG.orm.interact import DbInteractor
from specsnake.sfg_spectrum import SfgSpectrum, SfgAverager, AverageSpectrum

from SFG.orm.base_dtos import SFG


class SamplingDayAnalyzer:

    def __init__(self, bes: BoknisEckSamplingDay, svc: SfgValueCalculator):
        self.sampling_day = bes
        self.value_calculator = svc
        self.sampling_date = self.sampling_day.sampling_date

    def apply_svc(self):
        data = self.value_calculator([i.sfg for i in self.sampling_day.spectra if i.sfg.measurement_day is not None])
        return {f'{self.value_calculator.property_name}': data, 'date': self.sampling_date}


class SfgValueCalculator:

    def __init__(self, property_name: str,
                 reference_relator: Callable[[float, float], float],
                 reference_calculator: Callable[[List[SfgSpectrum]], float],
                 value_getter: Callable[[SfgSpectrum], float],
                 reference_getter: Callable[[SFG], List[SfgSpectrum]] = None):
        """

        :param property_name:
        :param reference_relator:
        :param reference_calculator:
        :param value_getter:
        :param reference_getter:
        """
        self.property_name = property_name
        self.reference_relator = reference_relator
        self.reference_calculator = reference_calculator
        self.value_getter = value_getter
        self.reference_getter = SfgValueCalculator.default_get_references \
            if reference_getter is None else reference_getter

        self.calculator = self.compose_calculator()

    def __call__(self, spectra: List[SFG]):
        return list(map(self.calculator, spectra))

    def compose_calculator(self):
        def calculator(s: SFG) -> float:
            # reference value
            references: List[SfgSpectrum] = self.reference_getter(s)
            blubb = [DbInteractor.construct_sfg(i) for i in references]
            ref_to_use = self.reference_calculator(blubb)
            ref_value = self.value_getter(ref_to_use)

            # spectrum
            spectrum_value = self.value_getter(DbInteractor.construct_sfg(s))
            print(spectrum_value, ref_value)
            return self.reference_relator(spectrum_value, ref_value)

        return calculator

    @staticmethod
    def default_get_references(s: SFG):
        return [i.sfg for i in s.measurement_day.references]


# Reference relators

def root_ratio(value, reference):
    return np.sqrt(value / reference)


def ratio_of_roots(value, reference):
    return np.sqrt(value) / np.sqrt(reference)


def ratio(value, reference):
    return value / reference


# Reference getters

def get_average_spectrum(spectra: List[SfgSpectrum], enforce_scale=True) -> AverageSpectrum:
    return SfgAverager(spectra, enforce_scale=enforce_scale).average_spectrum


# Property calculators (functions that extract a certain value from SFG spectra, e.g. the CH integral)

def ch_integral(s: SfgSpectrum):
    return s.calculate_ch_integral()


def ch_integral_root(s: SfgSpectrum):
    return s.calculate_ch_integral(s.root_baseline_correction)


def max_ch_intensity(s: SfgSpectrum):
    indices = s.get_xrange_indices(2750, 3000)
    return np.max(s.y[indices[0]:indices[1] + 1])


# Testcode section
itc = DbInteractor()
bes = itc.session.query(itc.be_sampling_day).all()
# temp = itc.session.query(itc.boknis_eck).all()[45:46]
# to_test = [i.sfg for i in temp]

svc = SfgValueCalculator('test',
                         reference_relator=root_ratio,
                         reference_calculator=get_average_spectrum,
                         value_getter=ch_integral)

master = pd.DataFrame([SamplingDayAnalyzer(day, svc).apply_svc() for day in bes])
