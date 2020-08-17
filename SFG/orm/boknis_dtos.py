from sqlalchemy import Column, Integer, Text, ForeignKey, Date, Boolean, TIMESTAMP, Time, Float
from sqlalchemy.orm import relationship

from SFG.orm.base_dtos import Base


class BoknisEck(Base):
    __tablename__ = 'boknis_eck'
    id = Column(Integer, primary_key=True)
    name = Column(Text, unique=True)
    specid = Column(Integer, ForeignKey('sfg.id'))
    table_entry_id = Column(Integer, ForeignKey('boknis_water_samples.id'))
    boknis_sampling_day_id = Column(Integer, ForeignKey('be_data.id'))

    sample_type = Column(Text)
    sampling_date = Column(Date)
    sample_number = Column(Integer)
    location_number = Column(Integer)
    is_mapped = Column(Boolean)
    depth = Column(Integer)

    sfg = relationship('SFG', back_populates='boknis')
    table_entry = relationship('BoknisWaterSamples', uselist=False, back_populates='spectrum')
    boknis_sampling_day = relationship('BoknisEckData', back_populates='spectra')

    def __repr__(self):
        temp = f'{self.sfg.name} | BoknisEck, sampled {self.sampling_date}, measured {self.sfg.measured_time.date()}'
        return temp


class BoknisDatabaseParameters(Base):
    __tablename__ = 'boknis_database_parameters'
    id = Column(Integer, primary_key=True)
    sampling_date_id = Column(Integer, ForeignKey('be_data.id'))
    Time = Column(TIMESTAMP)
    Cast = Column(Text)
    Label = Column(Text)
    depth = Column(Text)
    chlorophyll_a = Column(Text)
    boknis_sampling_day = relationship('BoknisEckData', uselist=False, back_populates='parameters')


@DeprecationWarning
class BoknisWaterSamples(Base):
    __tablename__ = 'boknis_water_samples'
    id = Column(Integer, primary_key=True)
    sfg_id = Column(Integer, ForeignKey('sfg.id'))
    sampling_day_id = Column(Integer, ForeignKey('be_data.id'))

    Sample = Column(Integer)
    sampler_no = Column(Integer)
    dips_per_sample = Column(Integer)
    sample_container_type = Column(Integer)
    drainage_time_or_depth = Column(Integer)
    volume_collected = Column(Integer)
    sea_surface_observational_codes = Column(Text)
    sampling_date = Column(Date)
    sampling_time = Column(Time)
    ship_no = Column(Integer)
    location_no = Column(Integer)
    Latitude = Column(Float)
    Longitude = Column(Float)
    wind_speed = Column(Float)
    wind_direction = Column(Text)
    wave_period = Column(Text)
    wave_height = Column(Text)
    air_temperature = Column(Float)
    water_temperature = Column(Float)
    pollutant_type = Column(Float)
    sampled_by = Column(Text)
    Storage = Column(Text)
    experiment_date = Column(Date)

    spectrum = relationship('BoknisEck', back_populates='table_entry')
    boknis_sampling_day = relationship('BoknisEckData', back_populates='samples')

    def __repr__(self):
        temp = f'BoknisEck table entry sampled: {self.sampling_date} measured: {self.experiment_date} number: {self.Sample}'
        return temp


class BoknisEckData(Base):
    __tablename__ = 'be_data'
    id = Column(Integer, primary_key=True)
    sampling_date = Column(Date, unique=True)

    spectra = relationship('BoknisEck', back_populates='boknis_sampling_day')
    samples = relationship('BoknisWaterSamples', back_populates='boknis_sampling_day')
    parameters = relationship('BoknisDatabaseParameters', back_populates='boknis_sampling_day')

    # bulk_no = Column(Integer)
    # sml_no = Column(Integer)

    # sml_coverage = Column(Float)
    # one_coverage = Column(Float)

    # sml_ch = Column(Float)
    # sml_oh1 = Column(Float)
    # sml_oh2 = Column(Float)
    # sml_dangling = Column(Float)

    # bulk_coverage = Column(Float)
    # bulk_ch = Column(Float)
    # bulk_oh1 = Column(Float)
    # bulk_oh2 = Column(Float)
    # bulk_dangling = Column(Float)

    # chlorophyll = Column(Float)
