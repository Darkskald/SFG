# -*- coding: utf-8 -*-
"""
Created on Tue Jun 25 11:42:19 2019

@author: lange

ORM-part of SQlalchemy
"""

import os
import timeit
from datetime import timedelta
from functools import partial
import itertools as ito
from typing import List

from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from sqlalchemy.orm.exc import NoResultFound

from SFG.orm.boknis_dtos import BoknisEck, BoknisDatabaseParameters, BoknisWaterSamples, MeasurementDay, BoknisEckData
from SFG.orm.gasex_dtos import GasExSfg, GasexLt, GasexSurftens, GasexStations, StationStat, GasexSamples, \
    GasexStationPlan, LiftOff
from SFG.orm.importer import Importer
from SFG.orm.base_dtos import *
from SFG.orm.parsing import *

import numpy as np

from SFG.spectrum.spectrum import SfgSpectrum


class DatabaseWizard:
    """The basic class for interaction with databases via sqlalchemy."""

    def __init__(self):

        # SQL init
        # todo: this path is not reliable
        self.engine = create_engine(f'sqlite:///{os.getcwd()}/orm.db', echo=False)
        self.Base = Base

        # new stuff
        # todo: order the orm classes and put the corresponding headings

        # SFG orm
        self.sfg = SFG
        self.regular_sfg = RegularSfg
        self.boknis_eck = BoknisEck

        # LT orm
        self.lt = Lt
        self.regular_lt = RegularLt

        # BE specific
        # todo measurement days for all SFG
        self.measurement_days = MeasurementDay
        self.be_data = BoknisEckData
        self.be_params = BoknisDatabaseParameters
        self.be_water_samples = BoknisWaterSamples

        # GasEx specific
        self.stations = GasexStations
        self.station_plan = GasexStationPlan
        self.station_stats = StationStat
        self.samples = GasexSamples
        self.gasex_surftens = GasexSurftens
        self.gasex_lt = GasexLt
        self.gasex_sfg = GasExSfg
        self.gasex_station_plan = GasexStationPlan
        self.gasex_lift_off = LiftOff

        # Misc
        self.substances = Substances
        self.ir = IR
        self.raman = Raman
        self.uv = UV

        # SQL creation
        self.Base.metadata.create_all(self.engine)
        Session = sessionmaker(bind=self.engine)
        self.session = Session()

    @staticmethod
    def to_array(string) -> np.ndarray:
        """Converts the raw data stored as strings back to numpy float ndarrays."""
        try:
            return np.fromstring(string, sep=",")
        except TypeError:
            return np.nan

    @staticmethod
    def construct_sfg(or_object, time_correction=False) -> SfgSpectrum:
        """A function constructing the SFG object from the orm declarative class."""
        meta = {"name": or_object.name, "time": or_object.measured_time}
        args = ("wavenumbers", "sfg", "vis", "ir")
        args = [ImportDatabaseWizard.to_array(getattr(or_object, i)) for i in args]
        s = SfgSpectrum(*args, meta)

        # as some spectra are recorded after midnight, it is necessary to date it back to get the correct DPPC reference
        if time_correction:
            if 0 <= s.meta["time"].hour < 8:
                s.meta["time"] -= timedelta(days=1)

        return s


# todo: correction factor of old surface tension function


class ImportDatabaseWizard(DatabaseWizard):
    """This class is designed to take care about the initial raw data import from the measurement files.
    It provides the same declarative classes as the base class DatabaseWizard and will integrate with
    the class for the actual data interaction WorkDatabaseWizard. It makes use of the Importer class,
    delegating all raw data import and preprocessing to it. Usually, the function performing the object-
    relational mapping is called "persist_x", but in case of complex import operations there may be
    two cooperatig functions write_x and persist_x."""

    def __init__(self):
        super().__init__()
        # initial import from raw data
        self.importer = Importer()
        self.persist_substances()
        self.persist_xy_spectra()
        self.persist_sfg()
        self.persist_lt()

        # station plan of the GasEx cruise, surface tensions and lift-off points
        self.importer.station_plan.to_sql('gasex_station_plan', con=self.engine, if_exists='append', index=False)
        self.importer.gasex_tension.to_sql('gasex_surftens', con=self.engine, if_exists='append', index=False)
        self.importer.gasex_lift_off.to_sql('gasex_lift_off', con=self.engine, if_exists='append', index=False)

        # import BoknisEck metadata
        self.importer.be_database_parameters.to_sql('boknis_database_parameters', con=self.engine, if_exists='append',
                                                    index=False)

        # import the Boknis Eck table by Kristian
        self.importer.water_samples.to_sql('boknis_water_samples', con=self.engine, if_exists='append',
                                           index=False)

        self.persist_stations(self.get_stations_and_samples())

        # metainfo refinement of regular samples
        self.add_regular_sfg()
        self.add_regular_lt()

        # GasEx-specific
        self.map_lift_off()
        self.map_station_plan()
        self.map_surface_tensions()
        self.populate_gasex_lt_sfg()

        # commit
        self.get_measurement_days()
        self.generate_boknis()
        self.match_boknis_table()
        self.populate_be_days()

        self.session.commit()

    # import of the raw data

    def persist_sfg(self):
        """Iterates over all folders containing SFG files. Calls the persist_sfg() for each file."""
        folders = [*self.importer.regular_sfg, *self.importer.gasex_sfg, *self.importer.boknis]
        self.engine.execute(
            self.sfg.__table__.insert(),
            folders
        )

    def persist_lt(self):
        """Iterates over all folders containing LT files. Calls the write_lt() for each file."""
        folders = [*self.importer.gasex_lt, *self.importer.lt]

        self.engine.execute(
            self.lt.__table__.insert(),
            folders
        )

    def persist_tensions(self):
        """Processes the surface tension dataset by iterating over the table rows,
         creating instances of the declarative class and filling it with the data. """
        # todo: what does this factor mean?tension.surface_tension = row["tension"] * 0.99703
        # self.session.add_all(tensions)

    def persist_substances(self):
        """Processes the subtances dataset by generating instances of the declarative substances
        class."""
        substances_orm = []
        for substance in self.importer.substances:
            s = Substances()
            for key in substance:
                setattr(s, key, substance[key])
            substances_orm.append(s)
        self.session.add_all(substances_orm)

    def persist_xy_spectra(self):
        """Generates declarative class instances for  UV, Raman and IR data."""
        spec_properties = {"raman": [self.raman, "wavenumbers", "intensity"],
                           "ir": [self.ir, "wavenumbers", "transmission"],
                           "uv": [self.uv, "wavelength", "absorbance"]}
        for key in spec_properties:
            # get the desired spectra from the importer via the name
            spec_models = [ImportDatabaseWizard.convert_xy_spectra(i, spec_properties[key]) for i in
                           getattr(self.importer, key)]

            self.session.add_all(spec_models)
            self.session.commit()

    # population of the GasEx-related tables

    def get_stations_and_samples(self) -> Dict[str, List]:
        """This function creates the station and sample tables from GasEx from the names of the SFG spectra
        and LT isotherms"""

        # Step1: fetching the name columns and creating sets for stations and samples by parsing the names
        gasex_sfg_names = self.session.query(self.sfg.name).filter(self.sfg.type == 'gasex_sfg').all()
        gasex_sfg_names = [i for i, in gasex_sfg_names]
        gasex_lt_names = self.session.query(self.lt.name).filter(self.lt.type == 'gasex_lt').all()
        gasex_lt_names = [i for i, in gasex_lt_names]
        names = list(map(get_hashes, (gasex_sfg_names + gasex_lt_names)))

        stations = set(i["station_hash"] for i in names)
        stations = {i: [] for i in stations}

        # Step 2: map the samples to the stations
        samples = set(i["sample_hash"] for i in names)
        for s in samples:
            stations[sample_hash_to_station_hash(s)].append(s)

        return stations

    def persist_stations(self, station_dict):
        id_counter = 1
        for key in station_dict:
            station_object = self.stations()
            station_object.id = id_counter
            station_object.hash = key
            meta_data = process_station_hash(key)
            station_object.type = get_station_type(station_dict[key][0])
            for meta in meta_data:
                setattr(station_object, meta, meta_data[meta])

            sample_objects = list(
                map(partial(self.create_sample_object, station_id=station_object.id), station_dict[key]))
            self.session.add_all(sample_objects)
            self.session.add(station_object)
            self.session.commit()
            id_counter += 1

    def create_sample_object(self, sample, station_id) -> GasexSamples:
        sample_object = self.samples()
        sample_object.station_id = station_id
        sample_object.sample_hash = sample
        meta_data = process_sample_hash(sample)
        for meta in meta_data:
            setattr(sample_object, meta, meta_data[meta])
        return sample_object

    # population of regular LT and SFG tables

    def add_regular_sfg(self):
        """This function iterates over all regular type SFG objects, performs the name refinement
        provided in refine_regular() and persists this information in the regular_sfg table."""
        names = self.session.query(self.sfg.id, self.sfg.name).filter(self.sfg.type == 'regular').all()
        regular_sfg_objects = []
        for id, name in names:
            data = refine_regular(name, self.get_substances())
            regular_sfg_objects.append(self.create_regular_object(name, data, id))

        self.session.add_all(regular_sfg_objects)
        self.session.commit()

    def add_regular_lt(self):
        names = self.session.query(self.lt.id, self.lt.name).filter(self.lt.type == 'lt').all()
        regular_lt_objects = []
        for id, name in names:
            data = refine_regular_lt(name)
            regular_lt_objects.append(self.create_regular_object(name, data, id, sfg=False))

        self.session.add_all(regular_lt_objects)
        self.session.commit()

    def create_regular_object(self, name, data, parent_id, sfg=True):
        if sfg:
            regular_object = self.regular_sfg(**data)
            regular_object.specid = parent_id
        else:
            regular_object = self.regular_lt(**data)
            regular_object.ltid = parent_id

        regular_object.name = name

        return regular_object

    # additional GasEx tables

    def map_lift_off(self):
        lift_off = self.session.query(self.gasex_lift_off).all()
        for l in lift_off:
            id = self.session.query(self.lt.id).filter(self.lt.name == l.name).one()[0]
            l.lt_id = id
        self.session.commit()

    def map_station_plan(self):
        stations = self.session.query(self.stations.hash, self.stations.id).all()
        for tup in stations:
            plan_entry = self.session.query(self.station_plan).filter(self.station_plan.hash == tup[0]).one()
            plan_entry.station_id = tup[1]
        self.session.commit()

    def map_surface_tensions(self):
        tensions = self.session.query(self.gasex_surftens).all()
        for t in tensions:
            sample_hash = get_hashes("d_" + t.name)["sample_hash"]
            id = self.session.query(self.samples.id).filter(self.samples.sample_hash == sample_hash).one()[0]
            t.sample_id = id
        self.session.commit()

    def populate_gasex_lt_sfg(self):
        sfg = self.session.query(self.sfg).filter(self.sfg.type == "gasex_sfg").all()
        lt = self.session.query(self.lt).filter(self.lt.type == "gasex_lt").all()
        orm_objects = []

        for s in sfg:
            gasex_sfg = self.gasex_sfg()
            hashes = get_hashes(s.name)
            sample_id = \
                self.session.query(self.samples.id).filter(self.samples.sample_hash == hashes["sample_hash"]).one()[0]
            gasex_sfg.sfg_id = s.id
            gasex_sfg.sample_id = sample_id
            orm_objects.append(gasex_sfg)

        for l in lt:
            gasex_lt = self.gasex_lt()
            hashes = get_hashes(l.name)
            sample_id = \
                self.session.query(self.samples.id).filter(self.samples.sample_hash == hashes["sample_hash"]).one()[0]
            gasex_lt.lt_id = l.id
            gasex_lt.sample_id = sample_id
            orm_objects.append(gasex_lt)

        self.session.add_all(orm_objects)
        self.session.commit()

    # misc
    def get_substances(self):
        """This function queries the db_wizard's session for the substances to use this information
        for the regular_refine() method."""
        substances = {}
        q = self.session.query(self.substances)
        for item in q:
            substances[item.abbreviation] = item.sensitizing
        return substances

    def get_measurement_days(self):
        boknis = self.session.query(self.sfg).filter(self.sfg.type == 'boknis_ref').all()
        other_dppc = self.session.query(self.regular_sfg).filter(self.regular_sfg.surfactant == "DPPC").filter(
            self.regular_sfg.sensitizer is not None).all()
        other_dppc = list(map(lambda x: x.sfg, other_dppc))

        temp = ito.groupby(boknis + other_dppc, key=lambda x: x.measured_time.date())
        measurement_days = {key: list(specs) for key, specs in temp}
        # todo: create an auxiliary dto mapping the DPPC reference to this entries
        for key in measurement_days:
            specs = measurement_days[key]
            md = self.measurement_days()
            md.date = key
            integrals = []
            for s in specs:
                s_t = self.construct_sfg(s)
                integral = round(s_t.calculate_ch_integral(), 4)
                integrals.append(integral)
            md.dppc_integral = np.sum(integrals) / len(specs)
            md.dppc_no = len(specs)
            self.session.add(md)
        self.session.commit()

    # todo: create a relation between the SFG measurements and the corresponding reference integral

    # boknis

    def generate_boknis(self):
        """Populates the SQL table for BoknisEck spectra with metadata obtained
        from the systematic names of the SFG table."""

        samples = self.session.query(self.sfg).filter(self.sfg.type == "boknis")
        boknis_specs = []

        for item in samples:
            boknis_spec = self.boknis_eck()
            boknis_spec.specid = item.id
            boknis_spec.name = item.name

            for ex in numreg:
                temp = re.search(ex, item.name)

                if temp is not None:
                    refined = item.name[temp.start():temp.end()]
                    temp = re.search("\d{1,2}", refined)
                    number = refined[temp.start():temp.end()]
                    number_match = True
                    break

            if number_match is False:
                repstr = f'Number parsing was not possible for {item.name}\n'
                self.log += repstr

            for ex in datereg:

                temp = re.search(ex, item.name)

                if temp is not None:
                    refined = item.name[temp.start():temp.end()]
                    temp = re.search("\d{8}", refined)
                    sampling_date = refined[temp.start():temp.end()]
                    date_match = True
                    break

            if date_match is False:
                print(f'Date parsing was not possible for {item.name}\n')

            if number is not None:
                number = int(number)

            boknis_spec.sample_number = number
            boknis_spec.sampling_date = ImportDatabaseWizard.convert_date(sampling_date)
            boknis_spec.is_mapped = False

            boknis_specs.append(boknis_spec)

        self.session.add_all(boknis_specs)
        self.session.commit()

    def match_boknis_table(self):
        entries = self.session.query(self.be_water_samples).all()
        list(map(self.map_boknis, entries))
        self.session.commit()

    def map_boknis(self, entry: BoknisWaterSamples):
        try:
            temp: BoknisEck = self.session.query(self.boknis_eck).join(self.sfg).filter(
                self.boknis_eck.sampling_date == entry.sampling_date).filter(
                self.boknis_eck.sample_number == entry.Sample).order_by(self.sfg.measured_time).all()[0]

            temp.is_mapped = True
            temp.table_entry_id = entry.id
            entry.sfg_id = temp.id
            if entry.sampler_no in (4, 3):
                temp.sample_type = "deep"
            elif entry.sampler_no in (1, 2):
                temp.sample_type = "sml"

        # todo: establish logging
        except (NoResultFound, IndexError) as e:
            # print(f'{entry} {e} has no suitable Spectrum!')
            pass

    def populate_be_days(self):
        dates = set([i for i, in self.session.query(self.boknis_eck.sampling_date).filter(
            self.boknis_eck.table_entry_id != None).all()])

        for d in dates:
            new_day = self.be_data()
            new_day.sampling_date = d
            self.session.add(new_day)
        self.session.commit()

        spectra = self.session.query(self.boknis_eck).filter(self.boknis_eck.table_entry_id != None).all()
        for s in spectra:
            day_id = self.session.query(self.be_data).filter(self.be_data.sampling_date == s.sampling_date).one().id
            s.boknis_sampling_day_id = day_id

        entries = self.session.query(self.be_water_samples).filter(self.be_water_samples.sfg_id != None).all()
        for e in entries:
            day_id = self.session.query(self.be_data).filter(self.be_data.sampling_date == e.sampling_date).one().id
            e.sampling_day_id = day_id

        self.session.commit()

    # auxiliary functions (static)
    @staticmethod
    def convert_xy_spectra(spectrum_dict, properties):
        model = properties[0]()
        # set the x value
        setattr(model, 'name', spectrum_dict["name"])
        setattr(model, properties[1], spectrum_dict["data"]["x"])
        setattr(model, properties[2], spectrum_dict["data"]["y"])
        return model

    @staticmethod
    def convert_date(date):
        temp = str(date)
        year = int(temp[0:4])
        month = int(temp[4:6])
        day = int(temp[6:])

        return datetime.date(year, month, day)


if __name__ == "__main__":
    start = timeit.default_timer()
    D = ImportDatabaseWizard()
    end = timeit.default_timer()
    print(end - start)
