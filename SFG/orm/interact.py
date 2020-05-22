import datetime
import os
from typing import Dict, List
import numpy as np

from SFG.orm.orm import DatabaseWizard
from SFG.spectrum.averagers import SfgAverager, DummyPlotter
from SFG.spectrum.spectrum import SfgSpectrum, LtIsotherm, BaseSpectrum


class DbInteractor(DatabaseWizard):

    def __init__(self):
        super().__init__()
        self.session.commit()

    def get_dppc_references(self) -> Dict[datetime.date, float]:
        """A function querying the sfg table for DPPC reference spectra, generating the corresponding objects,
        calculating the CH integral making use of the SfgAverager class and returning a dictionary of date objects
        with the corresponding intensities."""
        dates = {}

        q_dppc = self.session.query(self.sfg). \
            filter(self.sfg.name.op('GLOB')('*DPPC_*.*')). \
            filter(self.sfg.measured_time.between('2018-01-01', '2018-12-31')) \
            .filter(~self.sfg.name.contains('ppp'))

        for item in q_dppc:
            s = DbInteractor.construct_sfg(item)
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

    def fetch_by_specid(self, specid, sfg=True) -> BaseSpectrum:
        """Fetches the SFG spectrum with the id specid from the database and retuns it as an SFG spectrum object."""
        if sfg:
            spec = self.session.query(self.sfg).filter(self.sfg.id == specid).one()
            return self.construct_sfg(spec)
        else:
            lt = self.session.query(self.lt).filter(self.lt.id == specid).one()
            return self.construct_lt(lt)

    def get_spectrum_by_name(self, name) -> SfgSpectrum:
        """Returns the SFG spectrum object for a given file name"""
        temp = self.session.query(self.sfg).filter(self.sfg.name == name).one()
        return self.construct_sfg(temp)

    def get_spectrum_by_property(self, property_, target) -> SfgSpectrum:
        """A convenience function to collect spectra based on properties like surfactant, sensitizer etc."""
        temp = self.session.query(self.regular_sfg). \
            filter(getattr(self.regular_sfg, property_) == target).all()
        out = []
        for item in temp:
            out.append(self.get_spectrum_by_name(item.name))
        return out

    def convert_regular_to_lt(self, reg_lt) -> LtIsotherm:
        """Converts a RegularLt object directly into the Lt object of the spectrum module."""
        lt = self.session.query(self.lt).filter(self.lt.id == reg_lt.ltid).one()
        return DbInteractor.construct_lt(lt)

    def convert_regular_to_sfg(self, reg_sfg) -> SfgSpectrum:
        """Converts a RegularSfg object directly into the Sfg object of the spectrum module.
        It remains the former regular_sfg object as part of the new spectrum's meta attribute
        for access of the metadata stored in the regular_sfg object."""
        sfg = self.session.query(self.sfg).filter(self.sfg.id == reg_sfg.specid).one()
        temp = DbInteractor.construct_sfg(sfg)
        temp.meta["regular"] = reg_sfg
        return temp

    def map_data_to_dates(self, data) -> Dict[datetime.date, List[BaseSpectrum]]:
        """This function maps a list of spectra to their corresponding date of measurement
        """
        dates = {}
        for item in data:
            sampling_date = item.measured_time.date()
            if sampling_date not in dates:
                dates[sampling_date] = [item]
            else:
                dates[sampling_date].append(item)
        return dates

    def origin_preview_date(self, surfacant="NA", out_dir="out", max_size=6):
        temp = self.session.query(self.regular_sfg).filter(self.regular_sfg.surfactant == surfacant).all()
        temp = [self.session.query(self.sfg).filter(self.sfg.id == i.specid).one() for i in temp]
        dates = self.map_data_to_dates(temp)
        for key in dates:
            dir_name = out_dir + "/" + str(key)
            os.mkdir(dir_name)
            sfg_spectra = [self.construct_sfg(i) for i in dates[key]]

            for spec in sfg_spectra:
                df = spec.convert_to_export_dataframe()
                df.to_csv(f'{dir_name}/' + spec.name + ".csv", index=False, sep=";")

            sub_speclist = [sfg_spectra[i:i + max_size] for i in range(0, len(sfg_spectra), max_size)]
            for index, item in enumerate(sub_speclist):
                DummyPlotter(item, save=True, savedir=dir_name, savename=f'preview{index}').plot_all()

    # TODO: direct conversion from regular to Sfg orm object

    @staticmethod
    def to_array(string) -> np.ndarray:
        """Converts the raw data stored as strings back to numpy float ndarrays."""
        try:
            return np.fromstring(string, sep=",")
        except TypeError:
            return np.nan

    @staticmethod
    def construct_sfg(or_object) -> SfgSpectrum:
        """A function constructing the SFG object from the orm declarative class."""
        meta = {"name": or_object.name, "time": or_object.measured_time}
        args = ("wavenumbers", "sfg", "vis", "ir")
        args = [DbInteractor.to_array(getattr(or_object, i)) for i in args]
        s = SfgSpectrum(*args, meta)
        return s

    @staticmethod
    def construct_lt(or_object) -> LtIsotherm:
        """A function constructing the LT object from the orm declarative class."""
        args = (or_object.name, or_object.measured_time)
        add_args = ["time", "area", "apm", "surface_pressure", "lift_off"]
        add_args = [DbInteractor.to_array(getattr(or_object, i)) for i in add_args]
        l = LtIsotherm(args[0], args[1], *add_args)
        return l
