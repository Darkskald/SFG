from sqlalchemy import Column, Integer, Text, ForeignKey, UniqueConstraint, TIMESTAMP

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

    regular = relationship('RegularSfg', uselist=False, back_populates='sfg')
    gasex = relationship('GasExSfg', uselist=False, back_populates='sfg')
    boknis = relationship('BoknisEck', uselist=False, back_populates='sfg')


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
        representation = f"""RegularLT object (id {self.id}) {self.name}"""
        return representation
