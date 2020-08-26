from datetime import timedelta

from sqlalchemy import func

from SFG.natural_samples.gasex_processors import StationProcessor, SampleProcessor
from SFG.orm.interact import DbInteractor

import pandas as pd

interactor = DbInteractor()
stations = interactor.session.query(interactor.stations).all()

temp = sorted([s.get_corrected_doy() for s in stations])
print(len(temp))


"""
stations = interactor.session.query(interactor.stations).all()
samples = interactor.session.query(interactor.samples).all()

s = SampleProcessor(samples, interactor)

for sa in s.samples:
    print(s.get_corrected_salinity(sa))

for s in temp:
    if 0 <= s.measured_time.hour < 8:
        s.measured_time -= timedelta(days=1)
    temp = interactor.session.query(interactor.measurement_days).filter(func.DATE(s.measured_time) == interactor.measurement_days.date).all()
    print(temp)
temp = interactor.session.query(interactor.sfg).all()

stations = interactor.session.query(interactor.stations).all()
for s in stations:
    print(s.get_doy())
    
from SFG.orm import interact
from SFG.natural_samples.gasex_processors import StationProcessor
import itertools as ito

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


interactor = interact.DbInteractor()
stations = interactor.session.query(interactor.stations).all()
sp = StationProcessor(stations, interactor)
df = sp.get_station_data_frame()


"""

