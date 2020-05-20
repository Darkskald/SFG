# -*- coding: utf-8 -*-
"""
Created on Tue Jun 25 11:42:19 2019

@author: lange

ORM-part of SQlalchemy
"""
import datetime
import os
import re
import timeit
from functools import partial
from typing import Dict, List

from sqlalchemy import create_engine
from sqlalchemy.exc import IntegrityError
from sqlalchemy.orm import sessionmaker

from SFG.orm.importer import Importer
from SFG.orm.orm_classes import *
from SFG.orm.parsing import refine_regular_lt, get_hashes, sample_hash_to_station_hash, process_station_hash, \
    process_sample_hash, get_station_type, refine_regular


class DatabaseWizard:
    """The basic class for interaction with databases via sqlalchemy."""

    def __init__(self):

        # SQL init
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
        self.stations = Stations
        self.station_stats = StationStat
        self.samples = Samples
        self.gasex_surftens = GasexSurftens
        self.gasex_lt = GasexLt
        self.gasex_sfg = GasExSfg
        self.gasex_station_plan = GasexStationPlan

        # Misc
        self.substances = Substances
        self.ir = IR
        self.raman = Raman
        self.uv = UV

        # SQL creation
        self.Base.metadata.create_all(self.engine)
        Session = sessionmaker(bind=self.engine)
        self.session = Session()

    def factory(self, dic):
        """A convenience class to create the classes for sqlalchemy ORM. The information to create
        the tables are stored in dictionaries (class attributes of the DatabaseWizard)"""
        name = dic["__tablename__"]

        if name == "sfg":
            dic["__table_args__"] = ((UniqueConstraint("name", "type")),)

        elif name == "lt":
            dic["__table_args__"] = ((UniqueConstraint("name", "measured_time")),)
        return type(name, (self.Base,), dic)


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

        # todo this is strange
        self.persist_stations(self.get_stations_and_samples())

        self.add_regular_sfg()
        self.add_regular_lt()

        # commit
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
        #self.session.add_all(tensions)

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
        spec_properties = {"raman": [self.raman, "wavenumbers", "intensities"],
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

    def create_sample_object(self, sample, station_id) -> Samples:
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

    # todo: regular LT table poplulation
    # todo GasEx LT table population
    # todo GasEx SFG table population
    # todo add rest of PostProcessor functions an delete itf

    # misc

    def get_substances(self):
        """This function queries the db_wizard's session for the substances to use this information
        for the regular_refine() method."""
        substances = {}
        q = self.session.query(self.substances)
        for item in q:
            substances[item.abbreviation] = item.sensitizing
        return substances

    # auxiliary functions (static)
    @staticmethod
    def convert_xy_spectra(spectrum_dict, properties):
        model = properties[0]()
        # set the x value
        setattr(model, 'name', spectrum_dict["name"])
        setattr(model, properties[1], spectrum_dict["data"]["x"])
        setattr(model, properties[2], spectrum_dict["data"]["y"])
        return model


"""
letzlich passieren hier mehrere Dinge: die Namen werden geparst, die Metadaten extrahiet und die Metadaten
anschließend persistiert. Es wäre eine Überlegung wert, diesen Schritt in den Importer zu verlegen, da die Daten
dort ohnehin im RAM vorgehalten werden. Das Tokenizen mus auf jeen Fall von dem Persistieren klar getrennt werden

--> Hier muss eindeutig eine klare uns vernünftige Architektur her

--> das Parsing braucht Zugriff auf die Substances

--> das Persisting kann weiterhin im PostProcessor passieren

--> das Parsing kann irgendwie standardisiert weren
"""


class PostProcessor:
    """The PostProcessor class is designed to work with the ImportDatabaseWizard. While the IDW just
    initially populates the database with the raw data, the PP creates relationships (eg. between
    samples, stations and measurements) and adds additional data like Lift-off points."""

    def __init__(self, database_wizard, new=True):

        self.db_wizard = database_wizard
        self.substances = self.get_substances()
        if new:
            self.add_regular_info()
            self.add_regular_lt_info()

            self.populate_gasex_tables()
            self.map_samples()
            self.map_tensions()
            self.add_salinity()
            self.add_lift_off()

    # todo: add
    def add_regular_lt_info(self):
        """This function iterates over all regular type SFG objects, performs the name refinement
        provided in refine_regular() and persists this information in the regular_sfg table."""
        q = self.db_wizard.session.query(self.db_wizard.lt).filter \
            (self.db_wizard.lt.type == "lt")
        reg_specs = []

        for item in q:
            name = item.name
            meta_info = refine_regular_lt(name)
            reg_spec = self.db_wizard.regular_lt()
            reg_spec.name = name
            reg_spec.ltid = item.id

            for key in meta_info:
                setattr(reg_spec, key, meta_info[key])

            # ensure that BXn isotherms are interpreted correctly
            if reg_spec.surfactant is None and reg_spec.sensitizer is not None:
                reg_spec.surfactant = reg_spec.sensitizer
                reg_spec.sensitizer = None
            reg_specs.append(reg_spec)

        self.db_wizard.session.add_all(reg_specs)
        self.db_wizard.session.commit()


    # natural sample name processing

    # gasex management
    def populate_gasex_tables(self):
        """This function iterates over the gasex_sfg spectra and gasex_lt isotherms, extracts information
        about corresponding samples and stations and matches them to their ids. In the course of this
        method, the stations, samples, gasex_lt and gasex_sfg tables are populated."""

        q_sfg = self.db_wizard.session.query(self.db_wizard.sfg) \
            .filter(self.db_wizard.sfg.type == "gasex_sfg")

        q_lt = self.db_wizard.session.query(self.db_wizard.lt) \
            .filter(self.db_wizard.lt.type == "gasex_lt")

        for dataset in q_sfg, q_lt:

            for item in dataset:

                hashdic = PostProcessor.generate_hashdic(item.name)

                # todo: das ist absoluter Schwachsinn und frisst Ressourcen ohne Ende! --> Vorsortierung
                try:
                    station = self.db_wizard.stations()
                    station.hash = hashdic["station_hash"]
                    station.type = hashdic["station_type"]
                    station.date = hashdic["date"]
                    station.number = hashdic["station_number"]
                    self.db_wizard.session.add(station)
                    self.db_wizard.session.commit()
                except IntegrityError:
                    self.db_wizard.session.rollback()

                try:
                    sample = self.db_wizard.samples()
                    sample.sample_hash = hashdic["sample_hash"]
                    sample.location = hashdic["location"]
                    sample.type = hashdic["sample_type"]
                    sample.number = hashdic["sample_number"]
                    self.db_wizard.session.add(sample)
                    self.db_wizard.session.commit()
                except IntegrityError:
                    self.db_wizard.session.rollback()

                finally:
                    sample_id = self.db_wizard.session.query(self.db_wizard.samples) \
                        .filter(self.db_wizard.samples.sample_hash == hashdic["sample_hash"])

                    if dataset == q_lt:
                        temp = self.db_wizard.gasex_lt()
                    elif dataset == q_sfg:
                        temp = self.db_wizard.gasex_sfg()

                    temp.sample_id = sample_id.all()[0].id
                    temp.sample_hash = hashdic["sample_hash"]
                    temp.name = item.name
                    self.db_wizard.session.add(temp)
                    self.db_wizard.session.commit()

    def map_tensions(self):
        """Maps the surface tension table (GasEx) to the corresponding samples"""
        q = self.db_wizard.session.query(self.db_wizard.gasex_surftens)

        for tension in q:
            try:
                sample_hash = self.get_hashes(tension.name)["sample_hash"]
                q_stat = self.db_wizard.session.query(self.db_wizard.samples) \
                    .filter(self.db_wizard.samples.sample_hash == sample_hash).all()[0]
                tension.sample_id = q_stat.id
            except IndexError:
                print(sample_hash)
        self.db_wizard.session.commit()

    def add_lift_off(self):
        """A function adding the manually determined lift-off points to the gasex_lt table."""
        for _, lift_off in self.db_wizard.importer.gasex_lift_off.iterrows():
            try:
                q = self.db_wizard.session.query(self.db_wizard.lt) \
                    .filter(self.db_wizard.lt.name == lift_off["name"]).all()[0]
                q.lift_off = lift_off["lift_off"]
            except IndexError:
                with open("corrupt.txt", "a") as logfile:
                    logfile.write(lift_off["name"] + "\n")

        self.db_wizard.session.commit()

    """
    DEPRECATED
    
        def add_salinity(self):
        This function adds information about the "official" station label, geographical position and
        salinity (surface and depth) to the station's table

        for dic in self.db_wizard.importer.salinity:
            try:
                q = self.db_wizard.session.query(self.db_wizard.stations) \
                    .filter(self.db_wizard.stations.hash == dic["hash"]).all()[0]

                for key in dic:
                    if key != "hash":
                        setattr(q, key, dic[key])

            except IndexError:
                print(f'ERROR: Processing {dic["hash"]}')
    
        def map_samples(self):
        Maps the sample table to the corresponding station id (foreign key

        q = self.db_wizard.session.query(self.db_wizard.samples)
        for sample in q:
            station_hash = self.get_station_from_sample(sample.sample_hash)

            q_stat = self.db_wizard.session.query(self.db_wizard.stations) \
                .filter(self.db_wizard.stations.hash == station_hash).all()[0]
            sample.station_id = q_stat.id
        self.db_wizard.session.commit()
    
    def disconnect(self):
        self.db_wizard.session.close()
    """



if __name__ == "__main__":
    # P = PostProcessor(D)
    # P.disconnect()
    start = timeit.default_timer()
    D = ImportDatabaseWizard()
    end = timeit.default_timer()
    print(end - start)
