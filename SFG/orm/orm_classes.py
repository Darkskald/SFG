from sqlalchemy import Column, Integer, Text, ForeignKey, UniqueConstraint, TIMESTAMP, \
    Float, Date, Boolean

from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship

Base = declarative_base()


# independent


class SFG(Base):
    __tablename__ = "sfg"
    id = Column(Integer, primary_key=True)
    name = Column(Text)
    measured_time = Column(TIMESTAMP)
    type = Column(Text)
    wavenumbers = Column(Text)
    sfg = Column(Text)
    ir = Column(Text)
    vis = Column(Text)

    __table_args__ = (UniqueConstraint("name", "type"),)


class Lt(Base):
    __tablename__ = 'lt'
    id = Column(Integer, primary_key=True)
    name = Column(Text)
    type = Column(Text)
    measured_time = Column(TIMESTAMP)
    time = Column(Text)
    area = Column(Text)
    apm = Column(Text)
    surface_pressure = Column(Text)
    lift_off = Column(Text)

    __table_args__ = (UniqueConstraint("name", "measured_time"),)
    children = relationship("RegularLt", uselist=False,
                            back_populates="parent")

    children2 = relationship("GasexLt", uselist=False,
                             back_populates="parent")


class IR(Base):
    __tablename__ = 'ir'
    id = Column(Integer, primary_key=True)
    name = Column(Text)
    wavenumbers = Column(Text)
    transmission = Column(Text)


class Raman(Base):
    __tablename__ = 'raman'
    id = Column(Integer, primary_key=True)
    name = Column(Text)
    wavenumbers = Column(Text)
    intensity = Column(Text)


class UV(Base):
    __tablename__ = 'uv'
    id = Column(Integer, primary_key=True)
    name = Column(Text)
    wavelength = Column(Text)
    absorbance = Column(Text)


class BoknisDatabaseParameters(Base):
    __tablename__ = 'boknis_database_parameters'
    id = Column(Integer, primary_key=True)
    Name = Column(Text)
    Latidude = Column(Text)
    Longitude = Column(Text)
    Time = Column(TIMESTAMP)
    Cast = Column(Text)
    Label = Column(Text)
    depth = Column(Text)
    chlorophyll_a = Column(Text)


class GasexStationPlan(Base):
    __tablename__ = 'gasex_station_plan'
    id = Column(Integer, primary_key=True)
    Leg = Column(Text)
    time = Column(TIMESTAMP)
    station_number = Column(Text)
    Latidude = Column(Text)
    Longitude = Column(Text)
    hash = Column(Text)
    salinity_surface = Column(Text)
    salinity_depth = Column(Text)
    temperature_surface = Column(Text)
    temerature_depth = Column(Text)


class Substances(Base):
    __tablename__ = 'substances'
    id = Column(Integer, primary_key=True)
    name = Column(Text, unique=True)
    abbreviation = Column(Text)
    molar_mass = Column(Text)
    sensitizing = Column(Text)

    def __repr__(self):
        out = f"""compound {self.name}, abbreviation {self.abbreviation}, molar mass: {self.molar_mass}, sensitizing: {self.sensitizing}"""
        return out


# calculated


class BoknisEck(Base):
    __tablename__ = 'boknis_eck'
    id = Column(Integer, primary_key=True)
    name = Column(Text, unique=True)
    specid = Column(Integer, ForeignKey('sfg.id'))
    sample_type = Column(Text)
    sampling_date = Column(Date)
    sample_number = Column(Integer)
    location_number = Column(Integer)
    is_mapped = Column(Boolean)
    depth = Column(Integer)


class GasExSfg(Base):
    __tablename__ = 'gasex_sfg'
    id = Column(Integer, primary_key=True)
    name = Column(Text, ForeignKey('sfg.name'), unique=True)
    sample_id = Column(Integer, ForeignKey('samples.id'))
    sample_hash = Column(Text)


class RegularSfg(Base):
    __tablename__ = 'regular_sfg'
    id = Column(Integer, primary_key=True)
    name = Column(Text)
    specid = Column(Integer, ForeignKey('sfg.id'))
    surfactant = Column(Text)
    surfactant_vol = Column(Text)
    surfactant_conc = Column(Text)
    sensitizer = Column(Text)
    sensitizer_vol = Column(Text)
    sensitizer_conc = Column(Text)
    photolysis = Column(Text)
    sample_no = Column(Text)
    measurement_no = Column(Text)
    ratio = Column(Text)
    comment = Column(Text)


class MeasurementDay(Base):
    __tablename__ = 'measurement_days'
    id = Column(Integer, primary_key=True)
    date = Column(Date, unique=True)
    dppc_integral = Column(Float)


class BoknisEckData(Base):

    __tablename__ = 'be_data'
    id = Column(Integer, primary_key=True)
    sampling_date = Column(Date, unique=True)
    bulk_no = Column(Integer)
    sml_no = Column(Integer)

    sml_coverage = Column(Float)
    one_coverage = Column(Float)

    sml_ch = Column(Float)
    sml_oh1 = Column(Float)
    sml_oh2 = Column(Float)
    sml_dangling = Column(Float)

    bulk_coverage = Column(Float)
    bulk_ch = Column(Float)
    bulk_oh1 = Column(Float)
    bulk_oh2 = Column(Float)
    bulk_dangling = Column(Float)

    chlorophyll = Column(Float)


class RegularLt(Base):
    __tablename__ = 'regular_lt'
    id = Column(Integer, primary_key=True)
    name = Column(Text)
    ltid = Column(Integer, ForeignKey('lt.id'))
    surfactant = Column(Text)
    sensitizer = Column(Text)
    ratio = Column(Text)
    conc = Column(Text)
    spreading_volume = Column(Text)
    speed = Column(Text)
    sample_no = Column(Text)
    measurement_no = Column(Text)
    parent = relationship("Lt", uselist=False, back_populates="children")

    def __repr__(self):
        representation = f"""RegularLT object (id {self.id}) {self.name}"""
        return representation


class GasexLt(Base):
    __tablename__ = 'gasex_lt'
    id = Column(Integer, primary_key=True)
    ltid = Column(Integer, ForeignKey('lt.id'))
    name = Column(Text)
    sample_id = Column(Integer, ForeignKey('samples.id'))
    sample_hash = Column(Text)
    parent = relationship("Lt",  uselist=False, back_populates="children2")


class GasexSurftens(Base):
    __tablename__ = 'gasex_surftens'
    id = Column(Integer, primary_key=True)
    name = Column(Text)
    surface_tension = Column(Text)
    sample_id = Column(Integer, ForeignKey('samples.id'))


class Stations(Base):
    __tablename__ = 'stations'
    id = Column(Integer, primary_key=True)
    hash = Column(Text, unique=True)
    type = Column(Text)
    date = Column(TIMESTAMP)
    number = Column(Integer)
    label = Column(Text)
    longitude = Column(Text)
    latitude = Column(Text)
    surface_salinity = Column(Text)
    deep_salinity = Column(Text)
    surface_temperature = Column(Text)
    deep_temperature = Column(Text)


class StationStat(Base):
    __tablename__ = 'station_stats'
    id = Column(Integer, primary_key=True)
    doy = Column(Float)

    station_id = Column(Integer, ForeignKey("stations.id"), unique=True)
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


class Samples(Base):
    __tablename__ = 'samples'
    id = Column(Integer, primary_key=True)
    station_id = Column(Integer, ForeignKey('stations.id'))
    sample_hash = Column(Text, unique=True)
    location = Column(Text)
    type = Column(Text)
    number = Column(Integer)
    coverage = Column(Float)
    max_pressure = Column(Float)
    lift_off = Column(Float)
    surface_tension = Column(Float)
