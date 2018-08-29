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

#bar plot safe
to_plot = []
fix, ax1 = plt.subplots()
ax1.set_xlabel("day of cruise")
ax1.set_ylabel("average surface pressure/ mN/m")
ax2 = ax1.twinx()
ax2.set_ylabel("positive samples/ percent")
ax1.set_title("Surfactant occurence in GasEx 1 (June '18), \nmeasured by Langmuir Trough")
ax1.grid(True)
day_percent = []
for daylist in q.days.values():

    average = 0
    positive = 0
    negative = 0
    s = set([i.create_sample_hash() for i in daylist])
    count=(len(s))

    for isotherm in daylist:
            max_pres = isotherm.get_maximum_pressure()
            average += max_pres

            if max_pres < 2:
                negative += 1
            elif max_pres > 2 < 70:
                positive += 1

    average/=len(daylist)
    ratio = (positive/len(daylist))*100
    day_percent.append((isotherm.day, ratio, count))
    p=ax1.scatter(isotherm.day, average, color="green", s=72)


#ax1.legend(handles=legend_elements, scatterpoints=1).draggable()
rects = ax2.bar([a[0] for a in day_percent],[a[1] for a in day_percent], alpha=0.45)

for rect, day in zip(rects, day_percent):
    height = rect.get_height()
    width = rect.get_width()
    ax2.text(rect.get_x()+0.4*width, rect.get_y()+height*0.4, str((day[2])), color="blue")


ax2.axhline(50, color="red", alpha=0.4, antialiased=True,linewidth=2)
ax1.legend((p,), ["average surface pressure"], scatterpoints=1)
ax2.set_yticks([0, 25, 50, 75, 100])
plt.show()