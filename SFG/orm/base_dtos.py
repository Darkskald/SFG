from sqlalchemy import Column, Integer, Text, ForeignKey, UniqueConstraint, TIMESTAMP, Date, Float

from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship

Base = declarative_base()


class SFG(Base):
    __tablename__ = "sfg"
    id = Column(Integer, primary_key=True)
    name = Column(Text)
    measured_time = Column(TIMESTAMP)
    measurement_day_id = Column(Integer, ForeignKey('measurement_days.id'))
    type = Column(Text)
    wavenumbers = Column(Text)
    sfg = Column(Text)
    ir = Column(Text)
    vis = Column(Text)

    __table_args__ = (UniqueConstraint("name", "type"),)

    regular = relationship('RegularSfg', uselist=False, back_populates='sfg')
    gasex = relationship('GasExSfg', uselist=False, back_populates='sfg')
    boknis = relationship('BoknisEck', uselist=False, back_populates='sfg')
    measurement_day = relationship('MeasurementDay', back_populates='spectra')

    def __repr__(self):
        return (f'SFG spectrum DTO {self.name}')


class Lt(Base):
    __tablename__ = 'lt'
    __table_args__ = (UniqueConstraint("name", "measured_time"),)
    id = Column(Integer, primary_key=True)
    name = Column(Text)
    type = Column(Text)
    measured_time = Column(TIMESTAMP)
    time = Column(Text)
    area = Column(Text)
    apm = Column(Text)
    surface_pressure = Column(Text)

    regular = relationship('RegularLt', uselist=False, back_populates='lt')
    gasex = relationship('GasexLt', uselist=False, back_populates='lt')
    lift_off = relationship('LiftOff', uselist=False, back_populates='lt')

    def __repr__(self):
        return f'Langmuir trough isotherm DTO {self.name}'


class IR(Base):
    __tablename__ = 'ir'
    id = Column(Integer, primary_key=True)
    name = Column(Text)
    wavenumbers = Column(Text)
    transmission = Column(Text)

    def __repr__(self):
        return f'Infrared spectrum DTO {self.name}'


class Raman(Base):
    __tablename__ = 'raman'
    id = Column(Integer, primary_key=True)
    name = Column(Text)
    wavenumbers = Column(Text)
    intensity = Column(Text)

    def __repr__(self):
        return f'Raman spectrum DTO {self.name}'


class UV(Base):
    __tablename__ = 'uv'
    id = Column(Integer, primary_key=True)
    name = Column(Text)
    wavelength = Column(Text)
    absorbance = Column(Text)

    def __repr__(self):
        return f'UV spectrum DTO {self.name}'


class Substances(Base):
    __tablename__ = 'substances'
    id = Column(Integer, primary_key=True)
    name = Column(Text, unique=True)
    abbreviation = Column(Text)
    molar_mass = Column(Text)
    sensitizing = Column(Text)

    def __repr__(self):
        return f'Infrared spectrum DTO {self.name}'

    def __str__(self):
        out = f"""compound {self.name}, abbreviation {self.abbreviation}, molar mass: {self.molar_mass}, sensitizing: {self.sensitizing}"""
        return out


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

    sfg = relationship('SFG', back_populates='regular')

    def __repr__(self):
        return f'Regular SFG DTO {self.name}'


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

    lt = relationship("Lt", uselist=False, back_populates="regular")

    def __repr__(self):
        return f'Regular LT DTO {self.name}'


class MeasurementDay(Base):
    __tablename__ = 'measurement_days'
    id = Column(Integer, primary_key=True)
    date = Column(Date, unique=True)
    dppc_integral = Column(Float)
    dppc_no = Column(Integer)

    references = relationship('ReferenceSpectrum', back_populates="measurement_day")
    spectra = relationship('SFG', back_populates='measurement_day')


class ReferenceSpectrum(Base):
    __tablename__ = 'reference_spectra'
    id = Column(Integer, primary_key=True)
    sfg_id = Column(Integer, ForeignKey('sfg.id'))
    measurement_day_id = Column(Integer, ForeignKey('measurement_days.id'))

    measurement_day = relationship('MeasurementDay', uselist=False, back_populates='references')
    sfg = relationship('SFG')
