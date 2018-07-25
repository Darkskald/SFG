import numpy as np
import matplotlib.pyplot as plt


def generate_baseline(point_1, point_2):
    slope = (point_2[1] - point_1[1]) / (point_2[0] - point_1[0])
    intercept = point_2[1] - slope * point_2[0]

    def baseline(x):
        return (slope * x + intercept)

    return baseline


class vibration:


    def __init__(self, *args):
        if len(args) != 4:
            raise ValueError("Wrong number of peak parameters!")

        for argument in args:
            self.frequency = args[0]
            self.tau = args[1]
            self.amplitude = args[2]
            self.phasefaktor = args[3]


    def yield_lorentzian(self, wavenumbers):
        intensity = (self.frequency - wavenumbers) - 1j * self.tau
        intensity = 1 / intensity
        return intensity


    def calculate_phase(self, wavenumbers):
        intensity = self.yield_lorentzian(wavenumbers)
        phase = np.arccos(np.real(intensity) / np.abs(intensity)) + self.phasefaktor
        return phase


    def calc_total_function(self, wavenumbers):
        intensity = self.yield_lorentzian(wavenumbers)
        phase = self.calculate_phase(wavenumbers)
        # phase = 1
        total = np.abs(self.amplitude * intensity) * np.e ** (phase * 1j)
        return total


testarray = np.linspace(2450, 3840, 5000)

v = vibration(2882., 1, 1, 0)
v2 = vibration(2875., 1, 1, 0)

p = v.calc_total_function(testarray)
p2 = v2.calc_total_function(testarray)

p = np.absolute(p + p2) ** 2

plt.plot(testarray, p ** 2)

plt.show()