from pandas import DataFrame, notnull

from SFG.orm.interact import DbInteractor, StationProcessor, SampleProcessor

interactor = DbInteractor()
stations = interactor.session.query(interactor.stations).all()
samples = interactor.session.query(interactor.samples).all()

statpro = StationProcessor(stations, interactor)
print(SampleProcessor.map_samples_to_category(samples))