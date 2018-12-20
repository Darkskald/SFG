import numpy as np
from scipy import stats
from scipy.signal import savgol_filter
from scipy.interpolate import interp1d

class LtIsotherm:
    """A class to represent experimental Langmuir trough isotherms, handling time, area, area per molecule and
    surface pressure"""

    def __init__(self, *args, correct=True):

        self.name = args[0]
        self.measured_time = args[1]
        self.time = np.array(args[2].split(";")).astype(np.float)
        self.area = np.array(args[3].split(";")).astype(np.float)
        self.apm = np.array(args[4].split(";")).astype(np.float)
        self.pressure = np.array(args[5].split(";")).astype(np.float)
        self.compression_factor = self.calc_compression_factor()

        self.sample_hash = None
        self.speed = None
        self.measurement_number = None

        if correct is True:
            p = np.min(self.pressure)
            if p < 0:
                self.pressure += np.abs(p)

        self.monotonic = True
        try:
            self.correct_increase()
        except:
            self.monotonic = False


    def __str__(self):
        return self.name + " LtIsotherm Object"

    def __repr__(self):
        return self.name + " LtIsotherm Object"

    def __lt__(self, other):
        if self.get_maximum_pressure() < other.get_maximum_pressure():
            return True

    def find_increase(self):
        start = None
        end = None
        for i,j in enumerate(self.area[:-3]):
            try:
                if self.area[i] > self.area[i+1]:
                    if self.area[i] > self.area[i + 2]:
                        if self.area[i] > self.area[i + 3]:
                            start = i
                            break

            except IndexError:
                raise IndexError("Monotony first step failed!")

        for i,j in enumerate(self.area[start:-1]):
            try:
                if self.area[i] < self.area[i+1]:
                    end = i
            except IndexError:
                raise IndexError("Monotony second step failed!")

        return start, end


    def correct_increase(self):

        start,end = self.find_increase()
        self.area = self.area[start:end+1]
        self.pressure = self.pressure[start:end+1]
        self.time = self.time[start:end+1]
        self.compression_factor = self.compression_factor[start:end+1]
        self.apm = self.apm[start:end+1]








    def drop_ascii(self):
        """Drops an ascii file with semikolon-separated data in the form time;area;surface pressure. Intention
        is easy interfacing with external software like Excel or Origin"""

        with open(self.name + ".out", "w") as outfile:
            for a, b, c in zip(self.time, self.area, self.pressure):
                outfile.write(str(a) + ";" + str(b) + ";" + str(c) + "\n")

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
        return (self.area/max)

    def derive_pressure(self):
        """Calculates the difference quotient of the surface pressure with respect to the area.
        Useful for calculation of surface elasticity"""
        return np.diff(self.pressure) / np.diff(self.area)

    def same_sample(self, other):
        """Checks wether two isotherms belong to the same sample. This is the case if the same sample
        is measured several times in a row. Returns a bool."""
        if self.sample_hash == other.sample_hash:
            return True
        else:
            return False

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

    def calculate_elasticity(self):
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

    def smooth(self):
        """Performs a smooth operation of the measured pressure involving a Savitzky-Golay-filter"""
        return savgol_filter(self.pressure, 9, 3)

    def cut_away_decay(self, x_array):

        max = self.get_maximum_pressure()
        index = int(np.where(self.pressure==max)[0][0])
        return self.get_slice(x_array, 0, index)

    def find_lift_off(self, x_array, inc=0.1, anchor=False):
        x, y = self.cut_away_decay(x_array)
        if anchor is True:
            x,y = x[self.x_anchor:], y[self.x_anchor:]

        if len(x) < 100:
            new_x = np.linspace(np.min(x), np.max(x), 500)
            interpol = interp1d(x, y)
            new_y = interpol(new_x)

            x = new_x
            y = new_y

        tup = LtIsotherm.process_lift_off(x, y, inc)

        if tup is False:
            return False

        if tup[1] is None or tup[2]is None:
            print(self.name, len(x), len(y), tup[1], tup[2])

        lift_off = LtIsotherm.calculate_intersection(tup[1], tup[2])
        index = LtIsotherm.get_closest_index(self.create_pointlist(self.area), lift_off)
        lift_off = self.area[index], self.pressure[index]
        return lift_off

    def optimize_lift_off(self):

        initial = 0.02
        maxiter = 25
        anc = False

        init = self.find_lift_off(self.area, initial)

        if init[0] > self.area[self.x_anchor]:
            anc = True
            init = self.find_lift_off(self.area, initial, anchor=anc)

        if init[1] > 0.8:

            pressure = init[1]

            while pressure > 0.8 and maxiter > 0:

                    point2 = self.find_lift_off(self.area, initial+0.05, anchor=anc)

                    if point2 is False:
                        break

                    if point2[0] > self.area[self.x_anchor]:
                        break
                    else:
                        initial += 0.01
                        maxiter -= 1

        return self.find_lift_off(self.area, initial)


    @staticmethod
    def process_lift_off(x, y, inc=0.1):

        start = 0
        diff = 0
        increment = int(inc * len(x))
        line1 = None
        line2 = None
        x_out = None

        if (len(x) - increment * 2 < 0):
            return False



        while start < (len(x) - increment * 2):

            left_x = x[start:start + increment]
            left_y = y[start:start + increment]

            right_x = x[start + increment:start + (2 * increment)]
            right_y = y[start + increment:start + (2 * increment)]

            slope_l, intercept_l, r_value, p_value, std_err = \
                stats.linregress(left_x, left_y)

            slope_r, intercept_r, r_value, p_value, std_err = \
                stats.linregress(right_x, right_y)

            temp = np.abs(slope_l - slope_r)

            if temp >= diff:
                diff = temp
                line1 = [slope_l, intercept_l]
                line2 = [slope_r, intercept_r]
                x_out = x[start:start + (2 * increment)]

            start += 1

        return diff, line1, line2, x_out

    @staticmethod
    def calculate_intersection(line1, line2):

        x = (line2[1] - line1[1]) / (line1[0] - line2[0])

        y = line1[0] * x + line1[1]

        return x, y

    @staticmethod
    def get_closest_index(array_datapoints, check):

        d = 1000000000000
        index = None
        for point in array_datapoints:

            d_temp = np.sqrt((check[0] - point[0]) ** 2 + (check[1] - point[1]) ** 2)
            if d_temp < d:
                d = d_temp
                index = point[2]

        return index




    # todo: cut away decay before lift-off detection, smooth, check for very round isotherms