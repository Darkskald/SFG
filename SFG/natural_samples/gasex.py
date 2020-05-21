import os
import numpy as np

from SFG.orm.interact import WorkDatabaseWizard
from SFG.orm.orm_classes import Samples


class SampleWrapper:

    # todo: mit Gernot abklären was mit den anderen Lts passieren soll
    def __init__(self, sample: Samples):
        self.sample = sample

    def get_tension(self):
        try:
            # correct the surface tension value for the difference between calibration (20) and real lab temperature(21)
            raw = float(self.sample.tension.surface_tension) * 0.99703

            if self.sample.type == "deep":
                raw += SampleWrapper.correct_salinity(float(self.sample.station.station_plan.salinity_depth))
            else:
                raw += SampleWrapper.correct_salinity(float(self.sample.station.station_plan.salinity_surface))

            return round(raw, 2)

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

    # todo: fix the problems with the integration and baseline correction
    def get_coverage(self, dates):
        pass

    @staticmethod
    def correct_salinity(salinity):
        """A function yielding a salinity-dependent correction factor for the surface tension. The reference salinity
        where the factor equals zero is 17 PSU."""
        return 0.52552 - 0.0391 * salinity


if __name__ == "__main__":
    os.chdir("..")
    w = WorkDatabaseWizard()
    s = w.session.query(w.samples).all()
    for j in s:
        sw = SampleWrapper(j)
        print(sw.get_tension())
