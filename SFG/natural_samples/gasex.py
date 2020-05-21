import os
import numpy as np

from SFG.orm.interact import WorkDatabaseWizard
from SFG.orm.orm_classes import Samples


class SampleWrapper:

    def __init__(self, sample: Samples):
        self.sample = sample

    # todo: proper temperature + salinity correction
    def get_tension(self):
        try:
            return self.sample.tension.surface_tension
        except AttributeError:
            return None

    def get_max_pressure(self):
        try:
            lt_isotherms = list(map(WorkDatabaseWizard.construct_lt, [i.lt for i in self.sample.lt]))
            temp = sorted(lt_isotherms, key=lambda x: x.measured_time)
            pressure = round(temp[0].get_maximum_pressure(), 1)
            return pressure
        except IndexError:
            return None

    def get_lift_off(self):
        temp = sorted([i.lt for i in self.sample.lt], key=lambda x: x.measured_time)
        if len(temp) > 0:
            try:
                lo = temp[0].lift_off.lift_off
                lt = WorkDatabaseWizard.construct_lt(temp[0])
                return np.round(lo / (np.max(lt.area)), 3)
            except AttributeError:
                return None
        else:
            return None

    # todo: das kann man auch elegant über eine relation lösen
    # todo: fix the problems with the integration and baseline correction
    def get_coverage(self, dates):
        pass


if __name__ == "__main__":
    os.chdir("..")
    w = WorkDatabaseWizard()
    s = w.session.query(w.samples).all()
    for j in s:
        sw = SampleWrapper(j)
        sw.get_lift_off()
