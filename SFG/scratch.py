from pandas import DataFrame, notnull

from SFG.orm.interact import DbInteractor, SampleProcessor

interactor = DbInteractor()
samples = interactor.session.query(interactor.samples).all()

test = SampleProcessor(samples, interactor)
temp = test.get_list_of_sample_dicts()
temp = DataFrame(temp)

print(temp[notnull(temp["lift_off_compression_ration"])])