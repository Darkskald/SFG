from pandas import DataFrame, notnull

from SFG.orm.interact import DbInteractor, SampleProcessor

interactor = DbInteractor()
samples = interactor.session.query(interactor.samples).all()

test = SampleProcessor(samples, interactor).get_list_of_sample_dicts()
df = DataFrame(test)
print(df[["sample_hash", "cruise"]])