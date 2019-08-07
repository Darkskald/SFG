# -*- coding: utf-8 -*-
"""
Created on Tue Jun 25 11:42:19 2019

@author: lange

ORM-part of SQlalchemy
"""

from importer import Importer
from spectrum import SfgSpectrum, SfgAverager, LtIsotherm

from sqlalchemy import create_engine, Column, Integer, Text, ForeignKey, UniqueConstraint, TIMESTAMP, Float
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.exc import IntegrityError
from sqlalchemy.orm import sessionmaker
import numpy as np

import re
import datetime


class DatabaseWizard:

    Sfg = {
        "__tablename__": 'sfg',
        "id": Column(Integer, primary_key=True),
        "name": Column(Text),
        "measured_time": Column(TIMESTAMP),
        "measurer": Column(Text),
        "type": Column(Text),
        "wavenumbers": Column(Text),
        "sfg": Column(Text),
        "ir": Column(Text),
        "vis": Column(Text),
    }

    BoknisEck = {
        "__tablename__": 'boknis_eck',
        "id": Column(Integer, primary_key=True),
        "name": Column(Text, unique=True),
        "specid": Column(Integer, ForeignKey('sfg.id'))
    }

    GasExSfg = {
        "__tablename__": 'gasex_sfg',
        "id": Column(Integer, primary_key=True),
        "name": Column(Text, ForeignKey('sfg.name'), unique=True),
        "sample_id": Column(Integer, ForeignKey('samples.id')),
        "sample_hash": Column(Text)
    }

    RegularSfg = {
        "__tablename__": 'regular_sfg',
        "id": Column(Integer, primary_key=True),
        "name": Column(Text),
        "specid": Column(Integer, ForeignKey('sfg.id')),
        "surfactant": Column(Text),
        "surfactant_vol": Column(Text),
        "surfactant_conc": Column(Text),
        "sensitizer": Column(Text),
        "sensitizer_vol": Column(Text),
        "sensitizer_conc": Column(Text),
        "photolysis": Column(Text),
        "sample_no": Column(Text),
        "measurement_no": Column(Text),
        "comment": Column(Text)
    }

    Lt = {
        "__tablename__": 'lt',
        "id": Column(Integer, primary_key=True),
        "name": Column(Text),
        "type": Column(Text),
        "measured_time": Column(TIMESTAMP),
        "measurer": Column(Text),
        "time": Column(Text),
        "area": Column(Text),
        "apm": Column(Text),
        "surface_pressure": Column(Text),
        "lift_off": Column(Text),
    }

    RegularLt = {
        "__tablename__": 'regular_lt',
        "id": Column(Integer, primary_key=True),
        "name": Column(Text),
        "ltid": Column(Integer, ForeignKey('lt.id')),
        "surfactant": Column(Text),
        "sensitizer": Column(Text),
        "ratio": Column(Text),
        "conc": Column(Text),
        "spreading_volume": Column(Text),
        "speed": Column(Text),
        "sample_no": Column(Text),
        "measurement_no": Column(Text),
    }

    GasexLt = {
        "__tablename__": 'gasex_lt',
        "id": Column(Integer, primary_key=True),
        "name": Column(Text),
        "sample_id": Column(Integer, ForeignKey('samples.id')),
        "sample_hash": Column(Text, ForeignKey('gasex_lt.name')),
    }

    GasexSurftens = {
        "__tablename__": 'gasex_surftens',
        "id": Column(Integer, primary_key=True),
        "name": Column(Text),
        "surface_tension": Column(Text),
        "sample_id": Column(Integer, ForeignKey('samples.id')),
    }

    Substances = {
        "__tablename__": 'substances',
        "id": Column(Integer, primary_key=True),
        "name": Column(Text, unique=True),
        "abbreviation": Column(Text),
        "molar_mass": Column(Text),
        "sensitizing": Column(Text),
    }

    Stations = {
        "__tablename__": 'stations',
        "id": Column(Integer, primary_key=True),
        "hash": Column(Text, unique=True),
        "type": Column(Text),
        "date": Column(TIMESTAMP),
        "number": Column(Integer),
        "label": Column(Text),
        "longitude": Column(Text),
        "latitude": Column(Text),
        "surface_salinity": Column(Text),
        "deep_salinity": Column(Text),
        "surface_temperature": Column(Text),
        "deep_temperature": Column(Text),
    }

    StationStat = {

        "__tablename__": 'station_stats',
        "id": Column(Integer, primary_key=True),
        "station_id": Column(Integer, ForeignKey("stations.id"), unique=True),
        "plate_coverage": Column(Float),
        "plate_coverage_std": Column(Float),
        "plate_coverage_n": Column(Integer),

        "plate_lift_off": Column(Float),
        "plate_lift_off_std": Column(Float),
        "plate_lift_off_n": Column(Integer),

        "plate_tension": Column(Float),
        "plate_tension_std": Column(Float),
        "plate_tension_n": Column(Integer),

        "plate_max_pressure": Column(Float),
        "plate_max_pressure_std": Column(Float),
        "plate_max_pressure_n": Column(Integer),

        "screen_coverage": Column(Float),
        "screen_coverage_std": Column(Float),
        "screen_coverage_n": Column(Integer),

        "screen_lift_off": Column(Float),
        "screen_lift_off_std": Column(Float),
        "screen_lift_off_n": Column(Integer),

        "screen_tension": Column(Float),
        "screen_tension_std": Column(Float),
        "screen_tension_n": Column(Integer),

        "screen_max_pressure": Column(Float),
        "screen_max_pressure_std": Column(Float),
        "screen_max_pressure_n": Column(Integer),

        "sml_coverage": Column(Float),
        "sml_coverage_std": Column(Float),
        "sml_coverage_n": Column(Integer),

        "sml_lift_off": Column(Float),
        "sml_lift_off_std": Column(Float),
        "sml_lift_off_n": Column(Integer),

        "sml_tension": Column(Float),
        "sml_tension_std": Column(Float),
        "sml_tension_n": Column(Integer),

        "sml_max_pressure": Column(Float),
        "sml_max_pressure_std": Column(Float),
        "sml_max_pressure_n": Column(Integer),

        "deep_coverage": Column(Float),
        "deep_coverage_std": Column(Float),
        "deep_coverage_n": Column(Integer),

        "deep_lift_off": Column(Float),
        "deep_lift_off_std": Column(Float),
        "deep_lift_off_n": Column(Integer),

        "deep_tension": Column(Float),
        "deep_tension_std": Column(Float),
        "deep_tension_n": Column(Integer),

        "deep_max_pressure": Column(Float),
        "deep_max_pressure_std": Column(Float),
        "deep_max_pressure_n": Column(Integer),

        "sml_rawtension":  Column(Float),
        "deep_rawtension": Column(Float),
    }

    Samples = {
        "__tablename__": 'samples',
        "id": Column(Integer, primary_key=True),
        "station_id": Column(Integer, ForeignKey('stations.id')),
        "sample_hash": Column(Text, unique=True),
        "location": Column(Text),
        "type": Column(Text),
        "number": Column(Integer),
        "coverage": Column(Float),
        "max_pressure": Column(Float),
        "lift_off": Column(Float),
        "surface_tension": Column(Float),
    }

    IR = {
        "__tablename__": 'ir',
        "id": Column(Integer, primary_key=True),
        "name": Column(Text),
        "wavenumbers": Column(Text),
        "transmission": Column(Text),
    }

    Raman = {
        "__tablename__": 'raman',
        "id": Column(Integer, primary_key=True),
        "name": Column(Text),
        "wavenumbers": Column(Text),
        "intensity": Column(Text),
    }

    uv = {
        "__tablename__": 'uv',
        "id": Column(Integer, primary_key=True),
        "name": Column(Text),
        "wavelength": Column(Text),
        "absorbance": Column(Text),
    }

    def __init__(self):

        # SQL init
        self.engine = create_engine('sqlite:///orm.db', echo=False)
        self.Base = declarative_base()

        # ORM object creation
        self.sfg = self.factory(DatabaseWizard.Sfg)
        self.regular_sfg = self.factory(DatabaseWizard.RegularSfg)

        self.lt = self.factory(DatabaseWizard.Lt)
        self.regular_lt = self.factory(DatabaseWizard.RegularLt)

        self.boknis_eck = self.factory(DatabaseWizard.BoknisEck)

        self.substances = self.factory(DatabaseWizard.Substances)

        self.stations = self.factory(DatabaseWizard.Stations)
        self.station_stats = self.factory(DatabaseWizard.StationStat)
        self.samples = self.factory(DatabaseWizard.Samples)
        self.gasex_surftens = self.factory(DatabaseWizard.GasexSurftens)
        self.gasex_lt = self.factory(DatabaseWizard.GasexLt)
        self.gasex_sfg = self.factory(DatabaseWizard.GasExSfg)

        self.ir = self.factory(DatabaseWizard.IR)
        self.raman = self.factory(DatabaseWizard.Raman)
        self.uv = self.factory(DatabaseWizard.uv)

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
        self.persist_tensions()
        self.persist_substances()
        self.persist_spectra()

        self.write_sfg()
        self.write_lt()

        # commit
        self.session.commit()

    def persist_sfg(self, dic):
        """Generates an instance of the declarative sfg class and fills it with the data
        obtained from the importer."""
        temp_sfg = self.sfg()

        for key in ("name", "measured_time", "measurer", "type"):
            setattr(temp_sfg, key, dic[key])

        for key in("wavenumbers", "sfg", "ir", "vis"):
            _data = dic["data"][key]
            temp = ImportDatabaseWizard.nparray_to_str(_data)
            setattr(temp_sfg, key, temp)

        return temp_sfg

    def write_sfg(self):
        """Iterates over all folders containing SFG files. Calls the persist_sfg() for each file."""
        folders = self.importer.regular_sfg, self.importer.gasex_sfg, self.importer.boknis

        for folder in folders:
            sfgs = []
            for item in folder:
                sfgs.append(self.persist_sfg(item))

            self.session.add_all(sfgs)
            self.session.commit()

    def persist_lt(self, dic):
        """Generates an instance of the declarative lt class and fills it with the data
        obtained from the importer."""
        temp_lt = self.lt()
        for key in ("name", "measured_time", "type"):
            setattr(temp_lt, key, dic[key])

        for key in ("time", "area", "apm", "surface_pressure"):
            _data = dic["data"][key]
            temp = ImportDatabaseWizard.nparray_to_str(_data)
            setattr(temp_lt, key, temp)

        return temp_lt

    def write_lt(self):
        """Iterates over all folders containing LT files. Calls the write_lt() for each file."""
        folders = self.importer.gasex_lt, self.importer.lt

        for folder in folders:
            lts = []
            for item in folder:
                lts.append(self.persist_lt(item))

            self.session.add_all(lts)
            self.session.commit()

    def persist_tensions(self):
        """Processes the surface tension dataset by iterating over the table rows,
         creating instances of the declarative class and filling it with the data. """
        tensions = []

        for _, row in self.importer.gasex_tension.iterrows():
            tension = self.gasex_surftens()
            tension.name = row["name"]
            tension.surface_tension = row["tension"]*0.99703
            tensions.append(tension)

        self.session.add_all(tensions)

    def persist_substances(self):
        """Processes the subtances dataset by generating instances of the declarative substances
        class."""
        substances = []
        for _, row in self.importer.substances.iterrows():
            substance = self.substances()
            substance.name = row["name"]
            substance.abbreviation = row["short"]
            substance.molar_mass = row["molar_mass"]
            substance.sensitizing = row["sensitizing"]
            substances.append(substance)

        self.session.add_all(substances)

    def persist_spectra(self):
        """Generates declarative class instances for  UV, Raman and IR data."""
        raman = []
        ir = []
        uv = []

        for item in self.importer.raman:
            spec = self.raman()
            spec.name = item["name"]
            spec.wavenumbers = ",".join(item["data"]["x"].values.astype(str))
            spec.intensity = ",".join(item["data"]["y"].values.astype(str))
            raman.append(spec)
        self.session.add_all(raman)

        for item in self.importer.ir:
            spec = self.ir()
            spec.name = item["name"]
            spec.wavenumbers = ",".join(item["data"]["x"].values.astype(str))
            spec.transmission = ",".join(item["data"]["y"].values.astype(str))
            ir.append(spec)
        self.session.add_all(ir)

        for item in self.importer.uv:
            spec = self.uv()
            spec.name = item["name"]
            spec.wavelength = ",".join(item["data"]["x"].values.astype(str))
            spec.absorbance = ",".join(item["data"]["y"].values.astype(str))
            uv.append(spec)
        self.session.add_all(uv)

    # auxiliary functions
    @staticmethod
    def nparray_to_str(array):
        return ",".join(array.values.astype(str))


class PostProcessor:
    """The PostProcessor class is designed to work with the ImportDatabaseWizard. While the IDW just
    initially populates the database with the raw data, the PP creates relationships (eg. between
    samples, stations and measurements) and adds additional data like Lift-off points."""

    def __init__(self, database_wizard):

        self.db_wizard = database_wizard
        self.substances = self.get_substances()

        self.add_regular_info()
        self.add_regular_lt_info()

        self.populate_gasex_tables()
        self.map_samples()
        self.map_tensions()
        self.add_salinity()
        self.add_lift_off()

    # sfg refinement
    def refine_regular(self, namestring):
        """This function parses the names of the spectra in the sfg table with the type "regular" in
        order to extract additional metainformation. It makes use of regular expressions."""

        process_list = namestring.split("_")
        sample = re.compile('x\d')
        measurement = re.compile('#\d')
        photolysis = re.compile('\d{1,3}p')
        spread_vol = re.compile('\d{1,2}.\d{1,2}')
        spread_vol2 = re.compile('\d{1,2}')
        conc = re.compile('\d{1,2}mM')

        date = datetime.date(int(process_list[0][0:4]), int(process_list[0][4:6]),
                             int(process_list[0][6:]))
        surf = None
        sens = None
        surf_v = None
        sens_v = None
        surf_c = None
        sens_c = None
        sample_nr = None
        measurement_nr = None
        photo = None
        comment = None

        for item in process_list[1:]:

            if item in self.substances:
                if surf is None:
                    surf = item
                else:
                    if self.substances[item] == "y":
                        sens = item

            elif re.match(sample, item):
                sample_nr = item

            elif re.match(measurement, item):
                measurement_nr = item

            elif re.match(photolysis, item):
                photo = item

            elif re.match(conc, item):

                if surf_c is None:
                    surf_c = item
                else:
                    sens_c = item

            elif re.match(spread_vol, item) or re.match(spread_vol2, item):

                if surf_v is None:
                    surf_v = item
                else:
                    sens_v = item

            else:
                comment = item

        return {"sensitizer": sens, "date": date, "sample_no": sample_nr, "measurement_no": measurement_nr,
                "surfactant_conc": surf_c, "sensitizer_conc": sens_c, "surfactant": surf, "surfactant_vol": surf_v,
                "sensitizer_vol": sens_v, "comment": comment, "full_name": namestring, "photolysis": photo}

    def add_regular_info(self):
        """This function iterates over all regular type SFG objects, performs the name refinement
        provided in refine_regular() and persists this information in the regular_sfg table."""
        q = self.db_wizard.session.query(self.db_wizard.sfg).filter\
            (self.db_wizard.sfg.type == "regular")
        reg_specs = []

        for item in q:
            name = item.name
            meta_info = self.refine_regular(name)
            reg_spec = self.db_wizard.regular_sfg()
            reg_spec.name = name
            reg_spec.specid = item.id

            for key in meta_info:
                setattr(reg_spec, key, meta_info[key])

            reg_specs.append(reg_spec)

        self.db_wizard.session.add_all(reg_specs)
        self.db_wizard.session.commit()

    def add_regular_lt_info(self):
        """This function iterates over all regular type SFG objects, performs the name refinement
        provided in refine_regular() and persists this information in the regular_sfg table."""
        q = self.db_wizard.session.query(self.db_wizard.lt).filter\
            (self.db_wizard.lt.type == "lt")
        reg_specs = []

        for item in q:
            name = item.name
            meta_info = LtTokenizer.process(name)
            reg_spec = self.db_wizard.regular_lt()
            reg_spec.name = name
            reg_spec.ltid = item.id

            for key in meta_info:
                setattr(reg_spec, key, meta_info[key])

            reg_specs.append(reg_spec)

        self.db_wizard.session.add_all(reg_specs)
        self.db_wizard.session.commit()

    # auxiliary functions
    def get_substances(self):
        """This function queries the db_wizard's session for the substances to use this information
        for the regular_refine() method."""
        substances = {}
        q = self.db_wizard.session.query(self.db_wizard.substances)
        for item in q:
            substances[item.abbreviation] = item.sensitizing
        return substances

    # natural sample name processing
    @staticmethod
    def get_hashes(name):
        """ Extract the hashes from file names. Remember to strip the file ending as well as  the leading LT or date
        before passing namestring to this function."""

        temp = name.split("_")

        station_hash = temp[0] + temp[1][1]

        try:
            if temp[1][0] != "c":
                sample_hash = temp[0] + temp[1] + temp[2] + temp[3]

            else:
                sample_hash = temp[0] + temp[1] + "deep"

            dic = PostProcessor.process_sample_hash(sample_hash)

            dic["name"] = name
            dic["station_hash"] = station_hash
            dic["sample_hash"] = sample_hash
        except IndexError:
            print(temp)

        return dic

    @staticmethod
    def process_station_hash(hash):
        date = hash[:4]
        date = datetime.date(2018, int(date[0:2]), int(date[2:]))
        temp = hash[4]

        if temp in ("c", "r"):
            stype = "big"
        else:
            stype = "small"

        number = hash[5]

        return {"date": date, "station_type": stype, "station_number": number}

    @staticmethod
    def process_sample_hash(hash):
        dic = PostProcessor.process_station_hash(hash)
        location = hash[4]
        dic["location"] = location

        if location == "c":
            if "deep" in hash or "low" in hash:
                dic["sample_type"] = "deep"
                dic["sample_number"] = "-"

        elif location == "r":
            dic["sample_type"] = hash[6]
            if dic["sample_type"] in ("p","P"):
                dic["sample_number"] = hash[7]
                dic["sample_type"] = "p"

            elif dic["sample_type"] == "s":
                dic["sample_number"] = hash[8]

        elif location == "a":
            dic["sample_type"] = hash[6]
            dic["sample_number"] = hash[8]

        else:
            raise ValueError(f"Invalid sample hash {hash}")

        return dic

    @staticmethod
    def generate_hashdic(namestring):

        temp = namestring.split("_")[1:]
        temp = "_".join(temp)

        dic = PostProcessor.get_hashes(temp)
        return dic

    @staticmethod
    def get_station_from_sample(sample_hash):
        temp = sample_hash[:4]+sample_hash[5]
        return temp

    # gasex management
    def populate_gasex_tables(self):
        """This function iterates over the gasex_sfg spectra and gasex_lt isotherms, extracts information
        about corresponding samples and stations and matches them to their ids. In the course of this
        method, the stations, samples, gasex_lt and gasex_sfg tables are populated."""

        q_sfg = self.db_wizard.session.query(self.db_wizard.sfg)\
            .filter(self.db_wizard.sfg.type == "gasex_sfg")

        q_lt = self.db_wizard.session.query(self.db_wizard.lt)\
            .filter(self.db_wizard.lt.type == "gasex_lt")

        for dataset in q_sfg, q_lt:

            for item in dataset:

                hashdic = PostProcessor.generate_hashdic(item.name)

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

                sample_id = self.db_wizard.session.query(self.db_wizard.samples)\
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

    def map_samples(self):
        """Maps the sample table to the corresponding station id (foreign key"""

        q = self.db_wizard.session.query(self.db_wizard.samples)
        for sample in q:
            station_hash = self.get_station_from_sample(sample.sample_hash)

            q_stat = self.db_wizard.session.query(self.db_wizard.stations)\
                .filter(self.db_wizard.stations.hash == station_hash).all()[0]
            sample.station_id = q_stat.id
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

    def add_salinity(self):
        """This function adds information about the "official" station label, geographical position and
        salinity (surface and depth) to the station's table"""

        for dic in self.db_wizard.importer.salinity:
            try:
                q = self.db_wizard.session.query(self.db_wizard.stations)\
                    .filter(self.db_wizard.stations.hash == dic["hash"]).all()[0]

                for key in dic:
                    if key != "hash":
                        setattr(q, key, dic[key])

            except IndexError:
                print(f'ERROR: Processing {dic["hash"]}')

        self.db_wizard.session.commit()

    def add_lift_off(self):
        """A function adding the manually determined lift-off points to the gasex_lt table."""
        for _, lift_off in self.db_wizard.importer.gasex_lift_off.iterrows():
            try:
                q = self.db_wizard.session.query(self.db_wizard.lt)\
                    .filter(self.db_wizard.lt.name == lift_off["name"]).all()[0]
                q.lift_off = lift_off["lift_off"]
            except IndexError:
                with open("corrupt.txt", "a") as logfile:
                    logfile.write(lift_off["name"]+"\n")

        self.db_wizard.session.commit()


class WorkDatabaseWizard(DatabaseWizard):

    def __init__(self):
        super().__init__()
        self.session.commit()

    def get_dppc_references(self):
        """A function querying the sfg table for DPPC reference spectra, generating the corresponding objects,
        calculating the CH integral making use of the SfgAverager class and returning a dictionary of date objects
        with the corresponding intensities."""
        dates = {}

        q_dppc = self.session.query(self.sfg). \
            filter(self.sfg.name.op('GLOB')('*DPPC_*.*')). \
            filter(self.sfg.measured_time.between('2018-01-01', '2018-12-31')) \
            .filter(~self.sfg.name.contains('ppp'))

        for item in q_dppc:
            s = WorkDatabaseWizard.construct_sfg(item)
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

    @staticmethod
    def to_array(string):
        """Converts the raw data stored as strings back to numpy float ndarrays."""
        try:
            return np.fromstring(string, sep=",")
        except TypeError:
            return np.nan

    @staticmethod
    def construct_sfg(or_object):
        """A function constructing the SFG object from the orm declarative class."""
        meta = {"name": or_object.name, "time": or_object.measured_time}
        args = ("wavenumbers", "sfg", "vis", "ir")
        args = [WorkDatabaseWizard.to_array(getattr(or_object, i)) for i in args]
        s = SfgSpectrum(*args, meta)
        return s

    @staticmethod
    def construct_lt(or_object):
        """A function constructing the LT object from the orm declarative class."""
        args = (or_object.name, or_object.measured_time)
        add_args = ["time", "area", "apm", "surface_pressure", "lift_off"]
        add_args = [WorkDatabaseWizard.to_array(getattr(or_object, i)) for i in add_args]
        l = LtIsotherm(args[0], args[1], *add_args)
        return l


class LtTokenizer:

    regex = {
        "sample_no": re.compile(r'x\d'),
        "measurement_no": re.compile(r'#\d'),
        "ratio": re.compile(r'\dto\d'),
        "surfactant": re.compile(r'\wA'),
        "sensitizer": re.compile(r'BX\d'),
        "speed": re.compile(r'\d.\d*'),
        "spreading_volume": re.compile(r'\d\d'),
        "conc": re.compile(r'\d+')
    }

    @staticmethod
    def process(string):
        out = {}
        for token in string.split("_"):
            for item in LtTokenizer.regex:
                temp = re.match(LtTokenizer.regex[item], token)
                if temp is not None:
                    out[item] = token
        return out


if __name__ == "__main__":
    D = ImportDatabaseWizard()
    P = PostProcessor(D)


# todo: take care about the "measurer" field in LT
