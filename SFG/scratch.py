from datetime import timedelta

from sqlalchemy import func

from SFG.natural_samples.gasex_processors import StationProcessor, SampleProcessor
from SFG.orm.interact import DbInteractor

import pandas as pd

interactor = DbInteractor()

"""
stations = interactor.session.query(interactor.stations).all()
samples = interactor.session.query(interactor.samples).all()

s = SampleProcessor(samples, interactor)

for sa in s.samples:
    print(s.get_corrected_salinity(sa))
"""

temp = interactor.session.query(interactor.sfg).all()

for s in temp:
    if 0 <= s.measured_time.hour < 8:
        s.measured_time -= timedelta(days=1)
    temp = interactor.session.query(interactor.measurement_days).filter(func.DATE(s.measured_time) == interactor.measurement_days.date).all()
    print(temp)
