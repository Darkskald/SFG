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


class GasexLt(Base):
    __tablename__ = 'gasex_lt'
    id = Column(Integer, primary_key=True)
    lt_id = Column(Integer, ForeignKey('lt.id'))
    sample_id = Column(Integer, ForeignKey('gasex_samples.id'))

    lt = relationship("Lt", back_populates="gasex")
    sample = relationship('GasexSamples', back_populates='lt')


class GasexSurftens(Base):
    __tablename__ = 'gasex_surftens'
    id = Column(Integer, primary_key=True)
    name = Column(Text)
    surface_tension = Column(Text)
    sample_id = Column(Integer, ForeignKey('gasex_samples.id'))

    sample = relationship('GasexSamples', back_populates='tension')


class GasexStations(Base):
    __tablename__ = 'gasex_stations'
    id = Column(Integer, primary_key=True)
    hash = Column(Text, unique=True)
    type = Column(Text)
    date = Column(TIMESTAMP)
    number = Column(Integer)

    samples = relationship('GasexSamples', back_populates='station')
    station_plan = relationship('GasexStationPlan', uselist=False, back_populates='station')


@DeprecationWarning
class StationStat(Base):
    __tablename__ = 'station_stats'
    id = Column(Integer, primary_key=True)
    doy = Column(Float)

    station_id = Column(Integer, ForeignKey("gasex_stations.id"), unique=True)
    plate_coverage = Column(Float)
    plate_coverage_std = Column(Float)
    plate_coverage_n = Column(Integer)

    plate_lift_off = Column(Float)
    plate_lift_off_std = Column(Float)
    plate_lift_off_n = Column(Integer)

    plate_tension = Column(Float)
    plate_tension_std = Column(Float)
    plate_tension_n = Column(Integer)

    plate_max_pressure = Column(Float)
    plate_max_pressure_std = Column(Float)
    plate_max_pressure_n = Column(Integer)

    screen_coverage = Column(Float)
    screen_coverage_std = Column(Float)
    screen_coverage_n = Column(Integer)

    screen_lift_off = Column(Float)
    screen_lift_off_std = Column(Float)
    screen_lift_off_n = Column(Integer)

    screen_tension = Column(Float)
    screen_tension_std = Column(Float)
    screen_tension_n = Column(Integer)

    screen_max_pressure = Column(Float)
    screen_max_pressure_std = Column(Float)
    screen_max_pressure_n = Column(Integer)

    sml_coverage = Column(Float)
    sml_coverage_std = Column(Float)
    sml_coverage_n = Column(Integer)

    sml_lift_off = Column(Float)
    sml_lift_off_std = Column(Float)
    sml_lift_off_n = Column(Integer)

    sml_tension = Column(Float)
    sml_tension_std = Column(Float)
    sml_tension_n = Column(Integer)

    sml_max_pressure = Column(Float)
    sml_max_pressure_std = Column(Float)
    sml_max_pressure_n = Column(Integer)

    deep_coverage = Column(Float)
    deep_coverage_std = Column(Float)
    deep_coverage_n = Column(Integer)

    deep_lift_off = Column(Float)
    deep_lift_off_std = Column(Float)
    deep_lift_off_n = Column(Integer)

    deep_tension = Column(Float)
    deep_tension_std = Column(Float)
    deep_tension_n = Column(Integer)

    deep_max_pressure = Column(Float)
    deep_max_pressure_std = Column(Float)
    deep_max_pressure_n = Column(Integer)

    sml_rawtension = Column(Float)
    deep_rawtension = Column(Float)


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

    def to_basic_dict(self):
        return {"sample_hash": self.sample_hash, "cruise": 1 if self.sample_hash.startswith('06') else 2,
                "type": self.type, "number": self.number, "location": self.location}


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
