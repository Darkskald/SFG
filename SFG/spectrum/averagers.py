import datetime

import numpy as np

from SFG.spectrum import SfgSpectrum, DummyPlotter, LtIsotherm


class AverageSpectrum(SfgSpectrum):

    def __init__(self, wavenumbers, intensities, meta):
        self.wavenumbers = wavenumbers
        self.normalized_intensity = intensities
        self.meta = meta
        self.baseline_corrected = None
        self.regions = None

        # ensure nan-values in intensity and their corresponding wavenumbers are removed
        mask = np.isfinite(self.normalized_intensity)
        self.normalized_intensity = self.normalized_intensity[mask]
        self.wavenumbers = self.wavenumbers[mask]
        super().setup_spec()
        super().set_regions()


class AverageLt(LtIsotherm):

    def __init__(self, area, pressure, meta):
        # todo: invert the wavenunmber/ intensity scale to make it ascending
        self.area = area
        self.pressure = pressure
        self.meta = meta
        self._name = self.meta["name"]
        super().setup_spec()


class SfgAverager:
    # todo: throw an error and plot the spectra if the integral is NAN or zero!
    # todo: the benchmark function MUST display the integral value and the baseline
    """This class takes a list of SFG spectra and generates an average spectrum of them by interpolation and
    averaging. It is possible to pass a dictionary of date:dppc_integral key-value-pairs in order to calculate
    the coverage."""
    def __init__(self, spectra, references=None, enforce_scale=False, name="default", debug=False):
        self.failure_count = 0
        self.log = ""
        self.log += "Log file for averaging spectra\n"
        self.spectra = spectra
        self.references = references
        self.enforce_scale = enforce_scale
        self.name = name

        if len(self.spectra) == 0:
            print("Warning: zero spectra to average in SfgAverager!")
            self.average_spectrum = None
            self.integral = None
            self.coverage = None

        else:
            self.day_counter = {}
            self.average_spectrum = self.average_spectra()

            self.integral = self.average_spectrum.calculate_ch_integral()
            self.coverage = self.calc_coverage()

            if debug:
                if self.integral < 0:
                    self.benchmark()
                    print("Warning: negative integral value in SfgAverager!")
                    self.integral = 0
                    self.coverage = 0

    def average_spectra(self):
        """Function performing the averaging: it ensures that all spectra are interpolated to have the same shape,
        then they are averaged. A AverageSpectrum  object is constructed and returned."""
        to_average = []

        # sort spectra by length of the wavenumber array (lambda)
        if self.enforce_scale is False:
            self.spectra.sort(key=lambda x: x.yield_wn_length(), reverse=True)
            root_x_scale = self.spectra[0].x
        else:
            root_x_scale = SfgAverager.enforce_base()

        # get y values by interpolation and collect the y values in a list
        # collect the dates of measurement for DPPC referencing
        for item in self.spectra:

            if 0 <= item.meta["time"].hour < 8:
                item.meta["time"] -= datetime.timedelta(days=1)

            date = item.meta["time"].date()

            if date not in self.day_counter:
                self.day_counter[date] = 1
            else:
                self.day_counter[date] += 1

            new_intensity = np.interp(root_x_scale, item.x, item.y)
            mask = (root_x_scale > np.max(item.x)) | (root_x_scale < np.min(item.x))
            new_intensity[mask] = np.nan
            to_average.append(new_intensity)

        to_average = np.array(to_average)
        average = np.nanmean(to_average, axis=0)
        std = np.nanstd(to_average, axis=0)

        # prepare meta data for average spectrum
        if self.name == "default":
            newname = self.spectra[0].name + "baseAV"
        else:
            newname = self.name
        in_new = [n.name for n in self.spectra]
        s_meta = {"name": newname, "made_from": in_new, "std": std}

        #with open("blabla.txt", "a") as outfile:
           #outfile.write(f'name: {newname} x: {root_x_scale}, y: {average}\n')

        s = AverageSpectrum(root_x_scale, average, s_meta)

        return s

    def calc_reference_part(self):
        """Calculate the participation of each DPPC references. This is important if the spectra to average are
        measured on different sampling days. If, for example,  5 samples are to average and 3 of them are measured
        on one day, 2 on another, the final coverage is calculated by dividing the AveragedSpectrum integral by the
        weighted sum of the DPPC integrals of the corresponding days, in our example (2/5 * DPPC_1) + (3/5 * DPPC_2)"""

        spec_number = len(self.spectra)
        total = 0
        self.log += f'Start of the reference calculating section: \n'

        for date in self.day_counter:
            # divide by total number of spectra in the average
            self.log += f'date {date} divided by the number of spectra {spec_number}\n'
            self.day_counter[date] /= spec_number

            # multiply the weighting factor by the integral of the day
            try:
                self.log += f"""Now multiplying the factor {self.day_counter[date]} 
                by the reference integral {self.references[date]}\n"""

                self.day_counter[date] *= self.references[date]
                total += self.day_counter[date]
            except KeyError:
                self.failure_count += 1
                self.log += f'Error: no suitable DPPC reference found four date {date}\n'

        self.log += f'Finalizing calculation. The total factor now is {total}.\n'

        return total

    def calc_coverage(self):
        """A convenience function  to calculate the surface coverage"""

        if self.references is not None:
            dppc_factor = self.calc_reference_part()
            coverage = np.sqrt(self.integral / dppc_factor)
            return coverage

        else:
            #print(f'Coverage not available for reference samples, integral is {self.integral}!')
            pass

    def benchmark(self):
        self.create_log()
        l = [i for i in self.spectra]
        l.append(self.average_spectrum)
        p = DummyPlotter(l, save=True, savedir="benchmark", savename=self.spectra[0].name, special="AV")
        p.plot_all()

    def create_log(self):

        name = "benchmark/" + self.spectra[0].name + ".log"

        s = f'This average contains {len(self.spectra)} SFG spectra:\n'

        self.log += 80 * "-" + "\n"
        self.log += s
        for i in self.spectra:
            self.log += (i.name + "\n")

        self.log += 80 * "-" + "\n"
        s = f'integral: {self.integral}\ncoverage: {self.coverage}\n'

        #TODO replace this by proper logging
        #with open(name, "w") as outfile:
            #outfile.write(self.log)

    @staticmethod
    def enforce_base():
        reg1 = np.arange(2750, 3055, 5)
        reg2 = np.arange(3050, 3670, 20)
        reg3 = np.arange(3650, 3845, 5)
        new = np.concatenate((reg1, reg2, reg3), axis=None)
        return new


class LtAverager:

    def __init__(self, spectra):

        if len(spectra) == 0:
            raise ValueError("Warning: zero spectra to average in LtAverager!")

        self.spectra = spectra

    def average_lt(self, apm=True):
        """Function performing the averaging: it ensures that all spectra are interpolated to have the same shape,
        then they are averaged. A AverageSpectrum  object is constructed and returned."""
        to_average = []

        # sort spectra by length of the wavenumber array (lambda)
        self.spectra.sort(key=lambda x: np.min(x.area), reverse=True)

        if apm:
            area_var = 'apm'

        else:
            area_var = 'area'

        root_x_scale, c = self.spectra[0].cut_away_decay(getattr(self.spectra[0], area_var))
        # get y values by interpolation and collect the y values in a list
        for item in self.spectra:
            x_array = getattr(item, area_var)
            area, pressure = item.cut_away_decay(x_array)
            new_pressure = np.interp(root_x_scale[::-1], area[::-1],
                                      pressure[::-1])
            to_average.append(new_pressure)

        to_average = np.array(to_average)
        average = np.nanmean(to_average, axis=0)
        std = np.nanstd(to_average, axis=0)

        # prepare meta data for average spectrum
        newname = self.spectra[0].name + "baseAV"
        in_new = [n.name for n in self.spectra]
        s_meta = {"name": newname, "made_from": in_new, "std": std}
        s = AverageLt(root_x_scale[::-1], average, s_meta)

        return s