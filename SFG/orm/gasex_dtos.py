from sqlalchemy import Column, Integer, ForeignKey, Text, TIMESTAMP, Float
from sqlalchemy.orm import relationship

from SFG.orm.base_dtos import Base


class GasExSfg(Base):
    __tablename__ = 'gasex_sfg'
    id = Column(Integer, primary_key=True)
    sample_id = Column(Integer, ForeignKey('gasex_samples.id'))
    sfg_id = Column(Integer, ForeignKey('sfg.id'))

    sfg = relationship('SFG', back_populates='gasex')
    sample = relationship('GasexSamples', back_populates='sfg')

    def __repr__(self):
        return f'GasEx_SFG_{self.sample.sample_hash}'


class GasexLt(Base):
    __tablename__ = 'gasex_lt'
    id = Column(Integer, primary_key=True)
    lt_id = Column(Integer, ForeignKey('lt.id'))
    sample_id = Column(Integer, ForeignKey('gasex_samples.id'))

    lt = relationship("Lt", back_populates="gasex")
    sample = relationship('GasexSamples', back_populates='lt')

    def __repr__(self):
        return f'GasEx_LT_{self.sample.sample_hash}'


class GasexSurftens(Base):
    __tablename__ = 'gasex_surftens'
    id = Column(Integer, primary_key=True)
    name = Column(Text)
    surface_tension = Column(Text)
    sample_id = Column(Integer, ForeignKey('gasex_samples.id'))

    sample = relationship('GasexSamples', back_populates='tension')

    def __repr__(self):
        return f'GasEx_tension_{self.sample.sample_hash}'


class GasexStations(Base):
    __tablename__ = 'gasex_stations'
    id = Column(Integer, primary_key=True)
    hash = Column(Text, unique=True)
    type = Column(Text)
    date = Column(TIMESTAMP)
    number = Column(Integer)

    samples = relationship('GasexSamples', back_populates='station')
    station_plan = relationship('GasexStationPlan', uselist=False, back_populates='station')

    def __repr__(self):
        return f'GasEx_Station_{self.hash}_{self.type}'

    def to_basic_dict(self):
        """Return a basic dictionary of the station's properties"""
        return {"station_hash": self.hash, "date": self.date, "number": self.number, "doy": self.get_corrected_doy(),
                "station_type": self.type}

    def get_doy(self) -> int:
        """Convert the datetime of the station to the day of the year"""
        return self.date.timetuple().tm_yday

    def get_corrected_doy(self) -> float:
        """As there exist max. four stations a day and it is in general not desirable to plot all datapoints of one
         day of a year at exaclty the same location, this hack adds a fraction to the day of the year based on the sample
         number: corrected_doy = doy + (station_number-1)/4"""
        return self.get_doy() + (1 - self.number) / 4


class GasexSamples(Base):
    __tablename__ = 'gasex_samples'
    id = Column(Integer, primary_key=True)
    station_id = Column(Integer, ForeignKey('gasex_stations.id'))
    sample_hash = Column(Text, unique=True)
    location = Column(Text)
    type = Column(Text)
    number = Column(Integer)

    station = relationship('GasexStations', back_populates='samples')
    tension = relationship('GasexSurftens', uselist=False, back_populates='sample')
    sfg = relationship('GasExSfg', uselist=False, back_populates='sample')
    lt = relationship('GasexLt', back_populates='sample')

    def __repr__(self):
        return f'GasEx_Sample_{self.sample_hash}_object'

    def to_basic_dict(self):
        return {"sample_hash": self.sample_hash, "cruise": 1 if self.sample_hash.startswith('06') else 2,
                "doy": self.station.get_corrected_doy(), "type": self.type, "number": self.number,
                "location": self.location}


class GasexStationPlan(Base):
    __tablename__ = 'gasex_station_plan'
    id = Column(Integer, primary_key=True)
    station_id = Column(Integer, ForeignKey('gasex_stations.id'))
    Leg = Column(Text)
    time = Column(TIMESTAMP)
    station_number = Column(Text)
    Latitude = Column(Text)
    Longitude = Column(Text)
    hash = Column(Text)
    salinity_surface = Column(Text)
    salinity_depth = Column(Text)
    temperature_surface = Column(Text)
    temperature_depth = Column(Text)

    station = relationship('GasexStations', back_populates='station_plan')


class LiftOff(Base):
    __tablename__ = 'gasex_lift_off'
    id = Column(Integer, primary_key=True)
    lt_id = Column(Integer, ForeignKey('lt.id'))
    name = Column(Text)
    lift_off = Column(Float)
    sample = relationship("Lt")

    lt = relationship('Lt', uselist=False, back_populates='lift_off')
