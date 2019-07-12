# -*- coding: utf-8 -*-
"""
Created on Tue Jun 25 11:42:19 2019

@author: lange

ORM-part of SQlalchemy
"""

from importer import Importer

from sqlalchemy import create_engine, Column, Integer, Text, DateTime, ForeignKey
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker



class DatabaseWizard:

    Sfg = {
        "__tablename__": 'sfg',
        "id": Column(Integer, primary_key=True),
        "name": Column(Text, unique=True),
        "measured_time": Column(DateTime),
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
        "name": Column(Text, unique=True),
        "specid": Column(Integer),
        "surfactant": Column(Text, ForeignKey('sfg.id')),
        "surfactant_vol": Column(Text),
        "sensitizer": Column(Text),
        "sensitizer_vol": Column(Text),
        "photolysis": Column(Text),
        "comment": Column(Text)
    }

    Lt = {
        "__tablename__": 'lt',
        "id": Column(Integer, primary_key=True),
        "name": Column(Text, unique=True),
        "type": Column(Text),
        "measured_time": Column(DateTime),
        "time": Column(Text),
        "area": Column(Text),
        "apm": Column(Text),
        "surface_pressure": Column(Text),
        "lift_off": Column(Text),
    }

    GasexLt = {
        "__tablename__": 'gasex_lt',
        "id": Column(Integer, primary_key=True),
        "name": Column(Text, unique=True),
        "sample_id": Column(Integer, ForeignKey('samples.id')),
        "sample_hash": Column(Text, ForeignKey('gasex_lt.name')),
    }

    GasexSurftens = {
        "__tablename__": 'gasex_surftens',
        "id": Column(Integer, primary_key=True),
        "name": Column(Text, unique=True),
        "surface_tension": Column(Text),
        "sample_id": Column(Integer, ForeignKey('samples.id')),
    }

    Substances = {
        "__tablename__": 'substances',
        "id": Column(Integer, primary_key=True),
        "name": Column(Text, unique=True),
        "long_name:": Column(Text),
        "molar_mass": Column(Text),
        "sensitizing": Column(Text),
    }

    Stations = {
        "__tablename__": 'stations',
        "id": Column(Integer, primary_key=True),
        "hash": Column(Text, unique=True),
        "type": Column(Text),
        "date": Column(DateTime),
        "number": Column(Integer),
        "label": Column(Text),
        "longitude": Column(Text),
        "latidude": Column(Text),
        "surface_salinity": Column(Text),
        "deep_salinity": Column(Text),
    }

    Samples = {
        "__tablename__": 'samples',
        "id": Column(Integer, primary_key=True),
        "station_id": Column(Integer, ForeignKey('stations.id')),
        "sample_hash": Column(Text, unique=True),
        "location": Column(Text),
        "type": Column(Text),
        "number": Column(Integer),
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
        engine = create_engine('sqlite:///orm.db', echo=True)
        self.Base = declarative_base()

        # ORM object creation
        self.sfg = self.factory(DatabaseWizard.Sfg)
        self.regular_sfg = self.factory(DatabaseWizard.RegularSfg)
        self.lt = self.factory(DatabaseWizard.Lt)
        self.boknis_eck = self.factory(DatabaseWizard.BoknisEck)

        self.substances = self.factory(DatabaseWizard.Substances)

        self.stations = self.factory(DatabaseWizard.Stations)
        self.samples = self.factory(DatabaseWizard.Samples)
        self.gasex_surftens = self.factory(DatabaseWizard.GasexSurftens)
        self.gasex_lt = self.factory(DatabaseWizard.GasexLt)
        self.gasex_sfg = self.factory(DatabaseWizard.GasExSfg)
        self.lift_off = self.factory(DatabaseWizard.lift_off)

        self.ir = self.factory(DatabaseWizard.IR)
        self.raman = self.factory(DatabaseWizard.Raman)
        self.uv = self.factory(DatabaseWizard.uv)

        # SQL creation
        self.Base.metadata.create_all(engine)
        Session = sessionmaker(bind=engine)
        self.session = Session()

        # initial import from raw data
        self.importer = Importer()
        self.persist_tensions()

        # commit
        self.session.commit()

    def factory(self, dic):
        """A convenience class to create the classes for sqlalchemy ORM. The information to create
        the tables are stored in dictionaries (class attributes of the DatabaseWizard)"""
        name = dic["__tablename__"]
        return type(name, (self.Base,), dic)

    # mapping functions

    def persist_sfg(self):
        pass

    def persist_lt(self):
        pass

    def persist_tensions(self):
        tensions = []

        for _, row in self.importer.gasex_tension.iterrows():
            tension = self.gasex_surftens()
            tension.name = row["name"]
            tension.surface_tension = row["tension"]
            tensions.append(tension)

        self.session.add_all(tensions)

    def persist_lift_offs(self):
        # todo: this can be put together with map tensions, very redundant code
        # todo: no lift_off table so far
        pass



DatabaseWizard()




