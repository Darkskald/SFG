from functools import partial

import peakutils
import yaml
import matplotlib.pyplot as plt
import numpy as np
import logging

from sqlalchemy import extract

from SFG.orm.interact import DbInteractor
from specsnake.sfg_spectrum import SfgAverager, SfgSpectrum


class InvalidSelectorTypeException(Exception):
    pass


class SelectionRuleApplicator:
    allowed_types = ("boknis", "regular", "gasex")
    default_beginning = "2008-01-01"
    default_end = "2019-12-31"

    def __init__(self, interactor: DbInteractor, ruleset):
        logging.error(ruleset)
        self.check_type(ruleset['type'])
        self.interactor = interactor
        self.ruleset = ruleset

        self.query_map = {
            "boknis": self.interactor.session.query(self.interactor.boknis_eck),
            "regular": self.interactor.session.query(self.interactor.regular_sfg),
            "gasex": self.interactor.session.query(self.interactor.gasex_sfg)
        }

        self.type_map = {
            "regular": "regular_sfg",
            "boknis": "boknis_eck",
            "gasex": "gasex_sfg"
        }

        self.refine_map = {
            "regular": self.apply_regular_filter,
            "boknis": self.apply_boknis_filter,
            "gasex": self.apply_regular_filter
        }

    def check_type(self, type_to_check: str):
        if type_to_check not in self.allowed_types:
            raise InvalidSelectorTypeException

    def get_spectra(self):
        temp = self.get_query(self.ruleset['type'])
        temp = self.generate_temporal_filter(temp)
        temp = self.generate_special_filters(temp)

        return [i.sfg for i in temp.all()]

    def generate_temporal_filter(self, query):
        if 'beginning' in self.ruleset:
            beginning = self.ruleset['beginning']
            if 'end' in self.ruleset:
                # return query.filter(self.interactor.sfg.measured_time.between(self.ruleset['beginning'], self.ruleset['end']))
                end = self.ruleset['end']
            else:
                end = self.default_end
        else:
            beginning = self.default_beginning
            if 'end' in self.ruleset:
                end = self.ruleset['end']
            else:
                end = self.default_end
        logging.error(f'applying temporal filter with beginning {beginning}: end {end}')
        return query.filter(self.interactor.sfg.measured_time.between(beginning, end))

    def generate_special_filters(self, query):
        if 'special' in self.ruleset:
            return self.refine_map[self.ruleset['type']](query)
        else:
            return query

    def get_query(self, spec_type):
        return self.query_map[spec_type]

    def apply_regular_filter(self, query):
        for entry in self.ruleset['special']:
            entity = getattr(self.interactor, self.type_map[self.ruleset['type']])
            attribute = getattr(entity, entry)
            logging.error(f' entity:{entity}, attribute:{attribute}')
            query = query.filter(attribute == self.ruleset['special'][entry])
        return query

    def apply_boknis_filter(self, query):
        for entry in self.ruleset['special']:
            entity = getattr(self.interactor, self.type_map[self.ruleset['type']])
            if entry in ('year', 'month'):
                attribute = getattr(entity, 'sampling_date')
                query = query.filter(extract(entry, attribute) == self.ruleset['special'][entry])
                logging.error(f' entity:{entity}, attribute:{attribute}')
            else:
                attribute = getattr(entity, entry)
                query = query.filter(attribute == self.ruleset['special'][entry])
                logging.error(f' entity:{entity}, attribute:{attribute}')
        return query


class BatchProcessor:

    def __init__(self, interactor: DbInteractor, config_path):
        self.interactor = interactor
        self.config = BatchProcessor.load_config(config_path)
        logging.error(self.config)
        if 'spectra' in self.config:
            self.spectra = [self.unwrap_spectrum_config(i) for i in self.config["spectra"]]
        else:
            applicator = SelectionRuleApplicator(self.interactor, self.config['rule'])
            spectra = applicator.get_spectra()
            self.spectra = [{'object': self.interactor.construct_sfg(i), 'label': "", 'properties': {}} for i in
                            spectra]
            print(self.spectra)

        if "mpl_config" in self.config:
            plt.style.use(self.config["mpl_config"])

        if "average" in self.config:
            efc = self.config['average']['enforce_scale']
            self.spectra.append({
                'object': SfgAverager([i['object'] for i in self.spectra], enforce_scale=efc).average_spectrum,
                'label': self.config['average']['label'],
                'properties': self.config['average']['properties']
            })

    def get_spectrum_by_name(self, name: str):
        return self.interactor.session.query(i.sfg).filter(i.sfg.name == name).one()

    def plot_batch(self, plot_func):
        # call plotting stuff here
        plot_func(self.spectra)
        if self.config['interactive']:
            plt.show()
        else:
            path = f"{self.config['output_path']}/{self.config['name']}"
            plt.savefig(path)

    def unwrap_spectrum_config(self, spectrum_config):
        return {
            "object": self.interactor.construct_sfg(self.get_spectrum_by_name(spectrum_config['spectrum']['name'])),
            "label": spectrum_config['spectrum']['label'],
            "properties": spectrum_config['spectrum']['properties']
        }

    @staticmethod
    def load_config(path):
        with open(path) as infile:
            return yaml.load(infile, Loader=yaml.FullLoader)


# Plotting functions
def basic_plotting(spectra):
    for spec in spectra:
        sfg = spec['object']
        plt.plot(sfg.x, sfg.y, label=spec['label'], **spec['properties'])

    plt.xlabel(sfg.x_unit)
    plt.ylabel(sfg.y_unit)
    plt.legend()


def generate_n_shared_x(n):
    return plt.subplots(n, 1, sharex=True)


# with and without baseline
# ir/vis/raw demo
def demonstrate_normalization(spectra):
    fig, axis = generate_n_shared_x(2)
    ax1 = axis[0]
    ax2 = axis[1]
    for s in spectra:
        sfg = s['object']
        ax1.plot(sfg.x, sfg.raw_intensity)
        ax1.plot(sfg.x, sfg.ir_intensity)

        ax2.plot(sfg.x, sfg.y)


def create_broken_axis(limits):
    fig = plt.figure()
    gs = fig.add_gridspec(1, 2)

    # First pair
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1], sharey=ax1)
    ax1.set_xlim(limits[0], limits[1])
    ax2.set_xlim(limits[2], limits[3])

    ax1.spines['right'].set_visible(False)
    ax2.spines['left'].set_visible(False)
    ax1.yaxis.tick_left()
    # don't put tick labels at the top
    ax2.yaxis.tick_right()
    ax2.tick_params(labelright=False)

    d = .012  # how big to make the diagonal lines in axes coordinates
    # arguments to pass plot, just so we don't keep repeating them
    kwargs = dict(transform=ax1.transAxes, color='k', clip_on=False)
    ax1.plot((1 - d, 1 + d), (-d, +d), **kwargs)
    ax1.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)

    kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
    ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)
    ax2.plot((-d, +d), (-d, +d), **kwargs)

    # ax1.set_ylabel("norm. intensity/ arb. u.")

    # plot before
    # spec = by_name(itc, "20170901_BX9_5_x2_#1_5mM")
    # ax1.plot(spec.x, spec.y, label="initial", marker="o")
    # ax2.plot(spec.x, spec.y, label="initial", marker="o")

    # plot after
    # spec2 = by_name(itc, "20170901_BX9_5_x2_#2_30p_5mM")
    # ax1.plot(spec2.x, spec2.y, label="after 30 min.", marker="o")
    # ax2.plot(spec2.x, spec2.y, label="after 30 min.", marker="o")

    # ax2.legend()
    # gs.tight_layout(fig, rect=[0, 0.05, 1, 1])
    # fig.text(0.55, 0.04, s="wavenumber/ cm$^{-1}$", ha="center", va="center"
    return ax1, ax2


def plot_broken_axis(spectra, limits):
    ax1, ax2 = create_broken_axis(limits)
    for spec in spectra:
        sfg = spec['object']
        ax1.plot(sfg.x, sfg.y, **spec['properties'])
        ax2.plot(sfg.x, sfg.y, label=spec['label'], **spec['properties'])
    ax2.legend()
