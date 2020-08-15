from SFG.natural_samples.gasex_processors import StationProcessor, SampleProcessor
from SFG.orm.interact import DbInteractor

import pandas as pd

interactor = DbInteractor()
stations = interactor.session.query(interactor.stations).all()
samples = interactor.session.query(interactor.samples).all()

s = SampleProcessor(samples, interactor)

for sa in s.samples:
    print(s.get_corrected_salinity(sa))
