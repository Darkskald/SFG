import numpy as np
import pandas as pd

from SFG.spectrum.base_spectrum import BaseSpectrum


class LtIsotherm(BaseSpectrum):
    """A class to represent experimental Langmuir trough isotherms, handling time, area, area per molecule and
    surface pressure"""

    def __init__(self, name, measured_time, time, area, apm, pressure, lift_off=None, correct=True):
        self.name = name
        self.measured_time = measured_time
        self.time = time
        self.area = area
        self.apm = apm
        self.pressure = pressure
        self.compression_factor = self.calc_compression_factor()

        # todo: remove this hack which ensures the spectrum complies with the metaclass
        self.x = None
        self.y = None
        self.x_unit = "area/ cm$^{2}$"
        self.y_unit = "surface pressure/ mNm$^{-1}$"

        self.speed = None
        self.measurement_number = None

        self.lift_off = lift_off

        if correct is True:
            p = np.min(self.pressure)
            if p < 0:
                self.pressure += np.abs(p)

        self.setup_spec()

    def __lt__(self, other):
        if self.get_maximum_pressure() < other.get_maximum_pressure():
            return True

    def setup_spec(self) -> None:
        self._x = self.area
        self._y = self.pressure

    def drop_ascii(self) -> None:
        """Drops an ascii file with semikolon-separated data in the form time;area;surface pressure. Intention
        is easy interfacing with external software like Excel or Origin"""

        with open(self.name + ".out", "w") as outfile:
            for a, b, c in zip(self.time, self.area, self.pressure):
                outfile.write(str(a) + ";" + str(b) + ";" + str(c) + "\n")

    def convert_to_export_dataframe(self):
        """This function returns a Pandas dataframe, suitable for data export to
        origin and similar other programs"""
        data = {
            "time": self.time,
            "area": self.area,
            "area per molecule": self.apm,
            "surface pressure": self.pressure
        }
        return pd.DataFrame(data=data)

    def get_maximum_pressure(self, shrinked=None):
        """Returns the maximum measured surface pressure. Note: This property is uesd for the less-then
        operator implementation of this class!"""
        if shrinked == None:
            return np.max(self.pressure)
        else:
            try:
                return np.max(shrinked)
            except:
                # todo specify the type of error numpy will throw
                raise TypeError("Can not calc maximum for this operand")

    def calc_compression_factor(self):
        max = np.max(self.area)
        return self.area / max

    def derive_pressure(self):
        """Calculates the difference quotient of the surface pressure with respect to the area.
        Useful for calculation of surface elasticity"""
        return np.diff(self.pressure) / np.diff(self.area)

    def create_pointlist(self, x_array):
        """Returns a list containing the index, the x_array (usually area or time) and the surface pressure.
        This is used for example by the GUI functions to find the closest datapoint to a mouse-defined position
        in the plot"""

        output = []
        for i, (a, b) in enumerate(zip(x_array, self.pressure)):
            output.append((a, b, i))
        return output

    def get_slice(self, x_array, lower, upper, smooth=False):
        """Returns a slice of the x_array (usually time or area) defined by the lower and upper integer
        index."""

        x_out = x_array[lower:upper + 1]
        if smooth == False:
            y_out = self.pressure[lower:upper + 1]
        else:
            y_out = self.smooth()[lower:upper + 1]

        return x_out, y_out

    def calculate_elasticity(self) -> np.ndarray:
        """Returns the surface elasticity of the isotherm"""

        xdata = self.area[::-1]
        ydata = self.pressure[::-1]
        out = []

        for i in range(len(self.pressure) - 1):
            p = abs(ydata[i + 1] - ydata[i])
            a = abs((xdata[i + 1] - xdata[i]))
            A = (xdata[i + 1] - xdata[i]) / 2

            out.append(p / a * A)

        return np.array(out[::-1])

    def cut_away_decay(self, x_array):
        max_area = np.argmax(self.area < 31.2)
        return self.get_slice(x_array, 0, max_area)

    # todo: replace with numpy argmax function
    @staticmethod
    def get_closest_index(array_datapoints, check) -> int:
        d = 1000000000000
        index = None
        for point in array_datapoints:

            d_temp = np.sqrt((check[0] - point[0]) ** 2 + (check[1] - point[1]) ** 2)
            if d_temp < d:
                d = d_temp
                index = point[2]

        return index


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


class AverageLt(LtIsotherm):

    def __init__(self, area, pressure, meta):
        self.area = area
        self.pressure = pressure
        self.meta = meta
        self._name = self.meta["name"]
        super().setup_spec()