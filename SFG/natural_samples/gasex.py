import os
from functools import partial

from SFG.orm.interact import DbInteractor
from SFG.orm.boknis_dtos import BoknisEckData
from SFG.spectrum.averagers import SfgAverager


class BoknisWrapper:

    def __init__(self, day: BoknisEckData):
        self.day = day

    def get_by_type(self, type="sml"):
        temp = filter(lambda x: x.sample_type == type, self.day.spectra)
        return list(temp)

    @staticmethod
    def get_coverage(references, subset):
        sfgs = [i.sfg for i in subset]
        sfgs = list(map(partial(DbInteractor.construct_sfg, time_correction=True), sfgs))
        averager = SfgAverager(sfgs, references, enforce_scale=True, baseline=True)
        return averager.coverage


if __name__ == "__main__":
    import matplotlib.pyplot as plt

    os.chdir("..")
    w = DbInteractor()
    t = w.session.query(w.be_data).all()
    for t_x in t:
        b = BoknisWrapper(t_x)
        print(len(t))
        plt.scatter(b.day.sampling_date, b.get_coverage(w.references, b.get_by_type("sml")), color="red")
        plt.scatter(b.day.sampling_date, b.get_coverage(w.references, b.get_by_type("deep")), color="blue")
    plt.show()




    """
    t = w.session.query(w.sfg).filter(w.sfg.type == 'boknis_ref').all()
    temp = ito.groupby(t, key=lambda x: x.measured_time.date())
    # todo: include SfgAverager class in benchmark
    measurement_days = {key: list(specs) for key, specs in temp}
    for key in measurement_days:
        specs = measurement_days[key]
        for s in specs:
            s_t = DbInteractor.construct_sfg(s)
            integral = round(s_t.calculate_ch_integral(), 4)
            plt.plot(s_t.x, s_t.y, label=f'{s_t.name} {str(integral)}')

        plt.title(str(key))
        plt.legend()
        plt.savefig("/home/flo/Schreibtisch/dppc_check/" + str(key) + ".png")
        plt.cla()

    print(measurement_days)
    """
