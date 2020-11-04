import logging
from datetime import date
from typing import Dict, List, Set, Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import rcParams
from matplotlib.dates import MonthLocator, DateFormatter
from matplotlib.ticker import MaxNLocator
from pandas import DatetimeIndex
from specsnake.sfg_spectrum import SfgSpectrum, SfgAverager
from sqlalchemy import extract

from SFG.orm.base_dtos import MeasurementDay
from SFG.orm.boknis_dtos import BoknisEckSamplingDay, BoknisEck
from SFG.orm.interact import DbInteractor


class NoMeasurementDayError(Exception):
    pass


class StatisticSampleAnalyzer:
    years = (2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2017, 2018)
    qmap = {1: "Jan-Mar", 2: "Apr-Jun", 3: "Jul-Sep", 4: "Oct-Dec"}

    def __init__(self, itc: DbInteractor):
        self.itc = itc
        self.sampling_day_analyzers = [SamplingDayAnalyzer(i, self.itc) for i in
                                       itc.session.query(itc.be_sampling_day).all()]
        keys = ('sml', 'bulk', 'one', 'greater_one', 'date')
        self.data = pd.DataFrame(
            [{key: i.get_coverage_by_category_or_date(key) for key in keys} for i in
             self.sampling_day_analyzers]).replace(
            np.inf, np.nan).sort_values(by=['date'])

    def plot_scatter_matrix(self, path):
        rcParams['figure.figsize'] = (10.7, 8.2)
        rcParams['lines.markersize'] = 4
        rcParams['lines.markeredgewidth'] = 0.3
        rcParams['legend.fontsize'] = 'x-small'
        rcParams['xtick.minor.visible'] = True
        rcParams['ytick.minor.visible'] = True

        rcParams['xtick.major.size'] = 6
        rcParams['ytick.major.size'] = 6
        rcParams['xtick.minor.size'] = 3
        rcParams['ytick.minor.size'] = 3

        rcParams['xtick.major.width'] = 1.5
        rcParams['ytick.major.width'] = 1.5

        histprops = {'alpha': 0.6, 'edgecolor': 'black', 'facecolor': 'green'}
        scatterprops = {'edgecolor': 'black', 'aa': True, 'alpha': 0.7, 'color': 'blue'}

        fig = plt.figure(constrained_layout=True)
        gs = fig.add_gridspec(4, 4)

        # create histogram axis
        f1_1 = fig.add_subplot(gs[0, 0])
        f2_2 = fig.add_subplot(gs[1, 1])
        f3_3 = fig.add_subplot(gs[2, 2])
        f4_4 = fig.add_subplot(gs[3, 3])

        # create scatter axis
        f2_1 = fig.add_subplot(gs[1, 0])
        f3_1 = fig.add_subplot(gs[2, 0])
        f3_2 = fig.add_subplot(gs[2, 1])
        f4_1 = fig.add_subplot(gs[3, 0])
        f4_2 = fig.add_subplot(gs[3, 1])
        f4_3 = fig.add_subplot(gs[3, 2])

        # histograms
        f1_1.hist(self.data["sml"], **histprops, label="SML coverage")
        f2_2.hist(self.data["bulk"], **histprops, label="bulk coverage")
        f3_3.hist(self.data["one"], **histprops, label="1 m depth coverage")
        f4_4.hist(self.data["greater_one"], **histprops, label=" > 1 m depth coverage")
        for ax in (f1_1, f2_2, f3_3, f4_4):
            ax.legend()

        # scatter
        f2_1.scatter(self.data['sml'], self.data['bulk'], **scatterprops)
        f2_1.set_ylabel('bulk')
        pearson = self.data["sml"].corr(self.data["bulk"])
        t = f'R={pearson:.3f}'
        f2_1.text(0.5, 0.32, t, color='red')

        f3_1.scatter(self.data['sml'], self.data['one'], **scatterprops)
        f3_1.set_ylabel('1 m depth')
        pearson = self.data["sml"].corr(self.data["one"])
        t = f'R={pearson:.3f}'
        f3_1.text(0.43, 0.13, t, color='red')

        f3_2.scatter(self.data['bulk'], self.data['one'], **scatterprops)
        pearson = self.data["bulk"].corr(self.data["one"])
        t = f'R={pearson:.3f}'
        f3_2.text(0.46, 0.15, t, color='red')

        f4_1.scatter(self.data['sml'], self.data['greater_one'], **scatterprops)
        f4_1.set_xlabel('sml')
        f4_1.set_ylabel(' > 1 m depth')
        pearson = self.data["sml"].corr(self.data["greater_one"])
        t = f'R={pearson:.3f}'
        f4_1.text(0.47, 0.14, t, color='red')

        f4_2.scatter(self.data['bulk'], self.data['greater_one'], **scatterprops)
        f4_2.set_xlabel('bulk')
        pearson = self.data["bulk"].corr(self.data["greater_one"])
        t = f'R={pearson:.3f}'
        f4_2.text(0.58, 0.14, t, color='red')

        f4_3.scatter(self.data['one'], self.data['greater_one'], **scatterprops)
        f4_3.set_xlabel('1 m depth')
        pearson = self.data["one"].corr(self.data["greater_one"])
        t = f'R={pearson:.3f}'
        f4_3.text(0.16, 0.75, t, color='red')
        plt.savefig(f'{path}scatter_matrix.png')
        # plt.show()

    def plot_by_time(self):
        """BE time series plotting with rectangular stripes to help with orientation in the plot"""
        fig, ax = plt.subplots()
        ax2 = ax.twinx()

        months = MonthLocator(range(1, 13), bymonthday=1, interval=6)
        monthsFmt = DateFormatter("%b '%y")

        ax.scatter(self.data["date"], self.data["sml"], color="red")

        ts = pd.Series(self.data['sml'].rolling(window=2).mean().values, index=DatetimeIndex(self.data['date']))
        ax.plot(ts.interpolate(method='spline', order=3), color='red')
        # ax.plot(roll, test, color='red')

        ax.scatter(self.data["date"], self.data["bulk"], color="blue")
        # ax2.scatter(df["sampling_date"], df["chlorophyll"], color="green", marker="+")

        # ax.plot(df["sampling_date"], df["sml_"+param]*scale, color="red")
        # ax.plot(df["sampling_date"], df["bulk_"+param]*scale, color="blue")
        # ax2.plot(df["sampling_date"], df["chlorophyll"], color="green")

        # stripes to indicate certain times of the year
        for i in range(8, 20, 1):
            lower = str(date(2000 + i, 6, 1))
            upper = str(date(2000 + i, 12, 1))
            ax.axvspan(lower, upper, color="gray", alpha=0.4)

        # average lines
        # ax.axhline(df["sml_" + param].mean(), color="red", linestyle="dashed")
        # ax.axhline(df["bulk_" + param].mean(), color="blue", linestyle="dashed")
        """
        legend_elements = [Line2D([0], [0], marker='^', label='Bulk water',
                                  markerfacecolor='blue', mew=0.3, mec="blue", aa=True, linestyle=''),
                           Line2D([0], [0], marker='o', label='Surface microlayer',
                                  markerfacecolor='red', mew=0.3, mec="red", aa=True, linestyle=''),
                           Line2D([0], [0], marker='+', label='Chlorophyll a',
                                  markerfacecolor='green', mew=2, mec="green", aa=True, linestyle='', markersize=10)
                           ]
        """
        # ax.legend(handles=legend_elements)
        ax.xaxis.set_major_locator(months)
        ax.xaxis.set_major_formatter(monthsFmt)
        ax.xaxis.set_minor_locator(MonthLocator())
        ax.autoscale_view()
        fig.autofmt_xdate()

        ax.set_xlabel("time ")
        # ax.set_ylabel(leg)
        # ax2.set_ylabel("Chlorophyll a concentration/\n µg/L")
        rcParams['xtick.labelsize'] = 'small'

    def sampling_day_plot(self):
        for sda in self.sampling_day_analyzers:
            try:
                sml = sda.setup_averager(sda.type_map['sml'])
                bulk = sda.setup_averager(sda.type_map['bulk'])
                dppc = sda.get_reference_spectra()

                fig: plt.Figure = plt.figure()
                gs = fig.add_gridspec(3, 1)
                ax3 = fig.add_subplot(gs[2, 0])
                ax1 = fig.add_subplot(gs[0, 0], sharex=ax3)
                ax2 = fig.add_subplot(gs[1, 0], sharex=ax3)
                gs.update(hspace=0, wspace=0)

                for ax in (ax1, ax2):
                    ax.tick_params(axis='x', which='both', bottom=False, labelbottom=False)
                    ax.yaxis.set_major_locator(MaxNLocator(prune='lower', nbins=3))

                fig.text(0.5, 0.015, s="wavenumber/ cm$^{-1}$", ha="center", va="center")
                fig.text(0.014, 0.5, s="SFG intensity/ arb. u.", ha="center", va="center", rotation=90)
                # add plots

                ## SML
                if len(sml.spectra) > 0:
                    for s in sml.spectra:
                        ax1.plot(s.x, s.y, label=sda.name_transformer(s.name), alpha=0.6)
                        label = f'average sc:{sml.coverage:.3f}, int:{sml.average_spectrum.calculate_ch_integral():.4f}'
                    ax1.plot(sml.average_spectrum.x, sml.average_spectrum.y, label=label, color='red', marker='o')
                    boundaries = sml.average_spectrum.slice_by_borders(2750, 3000)
                    ax1.fill_between(sml.average_spectrum.x[boundaries[0]:boundaries[1] + 1],
                                     sml.average_spectrum.y[boundaries[0]:boundaries[1] + 1], color='red', alpha=1,
                                     zorder=9)
                    ax1.axvline(2750, linestyle='--', color='red')
                    ax1.axvline(3000, linestyle='--', color='red')

                ## Bulk
                if len(bulk.spectra) > 0:
                    for s in bulk.spectra:
                        ax2.plot(s.x, s.y, label=sda.name_transformer(s.name), alpha=0.6)
                        label = f'average sc:{bulk.coverage:.3f}, int:{bulk.average_spectrum.calculate_ch_integral():.4f}'
                    ax2.plot(bulk.average_spectrum.x, bulk.average_spectrum.y, label=label, color='blue', marker='^')
                    boundaries = bulk.average_spectrum.slice_by_borders(2750, 3000)
                    ax2.fill_between(bulk.average_spectrum.x[boundaries[0]:boundaries[1] + 1],
                                     bulk.average_spectrum.y[boundaries[0]:boundaries[1] + 1], color='blue', alpha=1,
                                     zorder=9)
                    ax2.axvline(2750, linestyle='--', color='blue')
                    ax2.axvline(3000, linestyle='--', color='blue')

                for s in dppc:
                    integral = s.calculate_ch_integral()
                    ax3.plot(s.x, s.y, label=f'{sda.name_transformer(s.name)}, int: {integral:.3f}')

                for ax in (ax1, ax2, ax3):
                    ax.grid(axis='x')
                    ax.legend(fontsize=8)

                fig.tight_layout()
                fig.suptitle(f'Boknis Eck sampling day {sda.sampling_day.sampling_date}')
                fig.savefig(f'/home/flo/Schreibtisch/Boknis/{sda.sampling_day.sampling_date}.png')
                plt.close(fig)
            except NoMeasurementDayError:
                pass

    # filter functions
    def maximum_of_year(self, year: int, key: str) -> float:
        """Return the maximum value of a property (key) of a certain year."""
        return self.data[self.data['date'].map(lambda x: x.year == year)][key].max()

    def normalize_to_maximum(self):
        temp = []
        for y in self.years:
            year_df = self.data[self.data['date'].map(lambda x: x.year == y)]
            for key in ('sml', 'bulk', 'one', 'greater_one'):
                new_col = f'avg_{key}'
                year_df[new_col] = year_df[key].map(lambda x: x / self.maximum_of_year(y, key))
            temp.append(year_df)
        return pd.concat(temp)

    def get_chlorophyll_data(self) -> pd.DataFrame:
        """Generate a dataframe of the chlorophyll values and their corresponding dates"""
        temp = []
        chlorophyll = pd.DataFrame([{"date": i.Time.date(), "chlorophyll": float(i.chlorophyll_a)} for i in
                                    self.itc.session.query(self.itc.be_chlorophyll).all()])
        for y in self.years:
            year_df = chlorophyll[chlorophyll['date'].map(lambda x: x.year == y)]
            new_col = f'avg_chlorophyll'
            max = year_df['chlorophyll'].max()
            year_df[new_col] = year_df['chlorophyll'].map(lambda x: x / max)
            temp.append(year_df)
        return pd.concat(temp)

    def plot_trimester(self):
        sml_samples = self.itc.session.query(self.itc.boknis_eck).join(self.itc.be_water_samples).filter(
            self.itc.be_water_samples.sampler_no.in_((1, 2))).all()
        sml_samples = [i for i in sml_samples if i.sfg.measurement_day is not None]

        trimester_map = {
            quartal: [self.map_and_normalize_to_reference(i) for i in sml_samples if self.map_to_quartal(i) == quartal]
            for quartal in (1, 2, 3, 4)}

        for key in trimester_map:
            average_spectrum = SfgAverager(trimester_map[key], enforce_scale=True, baseline=True).average_spectrum
            plt.plot(average_spectrum.x, average_spectrum.y, label=f'{self.qmap[key]} ({len(trimester_map[key])} spectra)')
        plt.legend()
        plt.xlabel(average_spectrum.x_unit)
        plt.ylabel(average_spectrum.y_unit)
        plt.savefig('/home/flo/Dropbox/PhdTex/latexthesistemplate/plots/boknis/quartal.png')

    @staticmethod
    def map_to_quartal(spectrum: BoknisEck):
        if spectrum.sampling_date.month in (1, 2, 3):
            return 1
        elif spectrum.sampling_date.month in (4, 5, 6):
            return 2
        elif spectrum.sampling_date.month in (7, 8, 9):
            return 3
        elif spectrum.sampling_date.month in (10, 11, 12):
            return 4

    def map_and_normalize_to_reference(self, spectrum: BoknisEck):
        reference = spectrum.sfg.measurement_day.dppc_integral
        sfg_spectrum = self.itc.construct_sfg(spectrum.sfg)
        sfg_spectrum.y = sfg_spectrum.normalize(external=reference)
        return sfg_spectrum


class YearAnalyzer:
    """This class is responsible to get the year's maxima for SML, bulk and 1 meter sampling."""
    years = (2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2017, 2018)

    def __init__(self, itc: DbInteractor):
        self.itc = itc

    def get_sampling_day_by_year(self, year):
        return self.itc.session.query(self.itc.be_sampling_day).filter(
            extract('year', self.itc.be_sampling_day.sampling_date) == year).all()


# todo:  the operation on a list of Boknis spectra may be handled on a higher order of hierarchy
class SamplingDayAnalyzer:

    def __init__(self, sampling_day: BoknisEckSamplingDay, itc: DbInteractor):
        self.sampling_day = sampling_day
        self.itc = itc
        self.type_map = self.map_to_type()

    def map_to_type(self) -> Dict[str, List[BoknisEck]]:
        """Map the sfg spectra contained in the spectrum list to the types 'sml' and 'bulk'."""
        # sampler 1,2: SML
        sml = {'sml': [spec for spec in self.sampling_day.spectra if spec.table_entry.sampler_no in (1, 2)]}
        # sampler 3,4: bulk
        bulk = {'bulk': [spec for spec in self.sampling_day.spectra if spec.table_entry.sampler_no in (3, 4)]}
        # samples from 1 meter depth
        one = {'one': [spec for spec in bulk['bulk'] if spec.table_entry.drainage_time_or_depth == 1]}
        # samples from depth > 1m
        greater_one = {'greater_one': [spec for spec in bulk['bulk'] if spec.table_entry.drainage_time_or_depth > 1]}
        sml.update(bulk)
        sml.update(one)
        sml.update(greater_one)
        return sml

    def convert_to_sfg(self, spectra: List[BoknisEck]) -> List[SfgSpectrum]:
        """Converts a list of BoknisEck dtos to the corresponding list of SfgSpectrum objects."""
        return [self.itc.construct_sfg(i.sfg) for i in spectra]

    def setup_averager(self, spectra: List[BoknisEck], baseline=True) -> SfgAverager:
        """Convert a list of BoknisEck dtos via their SfgSpectrum objects to an SfgAverager"""
        return SfgAverager(self.convert_to_sfg(spectra), self.itc.references, baseline=baseline, enforce_scale=True)

    def get_measurement_days(self) -> Set[MeasurementDay]:
        all_sfg = [i.sfg for i in self.sampling_day.spectra]
        return set(i.measurement_day for i in all_sfg)

    def get_reference_spectra(self) -> List[SfgSpectrum]:
        measurement_days = self.get_measurement_days()
        if None in measurement_days:
            raise NoMeasurementDayError
        measurement_day_spectra = [i.references for i in measurement_days]
        return [self.itc.construct_sfg(y.sfg) for x in measurement_day_spectra for y in x]

    def get_coverage_by_category_or_date(self, category: str, threshold: float = 1) -> Any:
        """Calculates the coverage for a given type of samples via an SfgAverager or returns the sampling date."""
        if category == 'date':
            return self.sampling_day.sampling_date
        else:
            result = self.setup_averager(self.type_map[category]).coverage
            if result is not None:
                return result if result < threshold else None

    @staticmethod
    def name_transformer(input: str):
        return input.replace('Wasserprobe', 'WP').replace('DPPC-Referenz', 'Ref')
    # plot all SML (plus average)
    # plot all bulk (plus average)
    # get samples from depth = 1 m
    # get bulk samples without 1m
    # calculate the coverages for each category
    # todo: is it possible to make the averager plot itself?


# testcode section

logging.basicConfig(level=logging.DEBUG, filename='/home/flo/Schreibtisch/flotest.log')

plt.style.use('../mpl_config/origin.mpltstyle')


# rcParams['legend.fontsize'] = 'xx-small'

def map_to_quartal(month: int):
    quartmap = {
        1: (1, 2, 3),
        2: (4, 5, 6),
        3: (7, 8, 9),
        4: (10, 11, 12)}
    for key in quartmap:
        if month in quartmap[key]:
            return key


itc = DbInteractor()
test = itc.session.query(itc.be_sampling_day).all()
q = StatisticSampleAnalyzer(itc)
# test = q.normalize_to_maximum()
# df = test #[test['date'] < date(2013, 1, 1)]
# df['month'] = df['date'].map(lambda x: x.month)
# df['quartal'] = df['date'].map(lambda x: map_to_quartal(x.month))

# fig = go.Figure()
# fig.add_trace(go.Box(x=df['quartal'], y=df['  sml'], boxmean=True))
# fig.add_trace(go.Box(x=df['date']+0.5, y=df['avg_bulk']))

# fig.show()
# q.data.to_csv('/home/flo/Schreibtisch/test.csv')
plotpath = "/home/flo/Dropbox/PhdTex/latexthesistemplate/plots/boknis/"
# q.normalize_to_maximum().to_csv("/home/flo/Dropbox/all_be_data.csv")
# q.get_chlorophyll_data().to_csv("/home/flo/Dropbox/chlorophyll_data.csv")
print(q.plot_trimester())
