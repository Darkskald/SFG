from ipy_interpreter import IpyInterpreter
from Classes import Plotter
import numpy as np
import matplotlib.pyplot as plt
"""A supporting module with useful functions for all the data management and analysis in the SFG module"""


def simple_analysis():
    Substances = {}
    i = Ipy_Interpreter()

    with open("name_info/Surfactants.txt", "r") as infile:
        for line in infile:
            collect = line.split(":")
            Substances[collect[0]] = collect[1].strip()

    with open("name_info/Sensitizers.txt", "r") as infile:
        for line in infile:
            collect = line.split(":")
            Substances[collect[0]] = collect[1].strip()

    for key in Substances:
        if key != "DPPC":
            dates = []
            # erster Schritt: alle Daten des Surfactants holen
            i.get("su " + key)

            # zweiter Schritt: Datum extrahieren
            for j in i.subset:
                dates.append(j.name.date)

            # an jedem Tag alle Probennummern extrahieren
            for d in dates:
                i.get("su " + key)
                daysamples = []
                i.refine("d " + d)
                for spectrum in i.subset:
                    daysamples.append(spectrum.name.sample_number)

                    # die einzelnen Samples zusammen plotten
                for sample in daysamples:
                    plotllist = [q for q in i.subset if q.name.sample_number == sample]
                    if len(plotllist) > 1:
                        plot_title = key + " " + str(d) + " " + str(sample)
                        P = Plotter(plotllist, title=plot_title)
                        P.custom_plot()

                i.clear()

        # subset wieder säubern
        i.clear()
"""
def yield_peaklist():
    yield a defined list of peaks separated from each other by minimum the threshold value in wavenumber
    pair_get = []
    out = []
    num = num
    if intensity == "default":
        intensity = self.normalized_intensity

    for i in range(len(self.wavenumbers)):
        pair_get.append([intensity[i], self.wavenumbers[i]])

    while len(out) < (num):

        if len(out) != 0:

            for i in range(len(out)):
                k = max(pair_get)

                if np.abs(out[i][1] - k[1]) < threshold:
                    pair_get.remove(k)
                    break
                if i == (len(out) - 1):
                    out.append(k)
                    pair_get.remove(k)


        else:
            k = max(pair_get)
            out.append(k)
"""
def calc_dish_area(diameter):
    """A auxialiary function to calculate the area of a teflon dish in square angstroms. Diameter given in cm."""
    radius = diameter * 0.5
    area = np.pi * radius ** 2
    area = area * 10 ** 16  # conversion to square angstroms
    return area


def calc_area_per_molecule(area, concentration, volume):
    """The function calculates the area per molecule. The area should be given in square angstroms, the
    concentration in milimole per liter and the volume in microleter"""

    concentration = concentration * 10 ** -3  # conversion in mol per liter
    volume = volume * 10 ** -6  # conversion in liter
    amount = volume * concentration
    molecules = (6.022 * 10 ** 23) * amount  # number of molecules
    area_per_molecule = area / molecules

    return area_per_molecule

#advanced plotting functionality







# Test code section
I = IpyInterpreter()
I.update()
I.get("se PA")
set = I.subset
