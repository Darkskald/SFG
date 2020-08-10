from SFG.natural_samples.gasex_processors import StationProcessor
from SFG.orm.interact import DbInteractor

import pandas as pd

interactor = DbInteractor()
stations = interactor.session.query(interactor.stations).all()
samples = interactor.session.query(interactor.samples).all()

statpro = StationProcessor(stations, interactor)
test = samples[0:3]


pd.DataFrame(statpro.get_list_of_station_dicts()).to_csv("~/Schreibtisch/dideldum.csv")

