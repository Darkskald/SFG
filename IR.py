
import csv
import numpy as np
import matplotlib.pyplot as plt
from Classes import *

class IR_Spectrum:
    
    def __init__(self,name,wavenumbers,intensities):
        self.name = name
        self.wavenumbers = wavenumbers
        self.intensities = intensities

class Raman_Spectrum:
    
    def __init__(self,name,wavenumbers,intensities):
        self.name = name
        self.wavenumbers = wavenumbers
        self.intensities = intensities

class IR_Plotter:
    """This class takes a list of IR,Raman or SFG Spectra as input and provides plotting functionality"""
    def __init__(self,speclist,title="Default"):
        self.speclist = speclist
        self.title = title


    def simple_plot(self):
        """This function is made to overlap spectra of a common type and does not provide double y axis!"""
        
        fig = plt.figure()
        ax = plt.subplot(111)
        for spectrum in self.speclist:
            ax.plot(spectrum.wavenumbers,spectrum.intensities)
        
        
        
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])      
        ax.grid()
        ax.set_title(self.title)
        ax.set_xlabel("Wavenumber/ $cm^{-1}$")
        ax.set_ylabel("Intensity/ a.u.")
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        plt.show()






def get_iraman_data(filename):
    with open (filename,"r") as infile:
        #function can handle IR files as well as Raman files
        c = csv.reader(infile, delimiter="\t")
        wavenumbers = []
        intensities = []
        for line in c:
            wavenumbers.append(line[0])
            intensities.append(line[1])
        return(wavenumbers,intensities)








#test code section
a = get_iraman_data("BX9")
b = get_iraman_data("BX9.dpt")


d = FileFetcher("20170901_BX9_5_x1_#2_5mM.sfg").sfg
c = (d.wavenumbers,d.normalized_intensity)


fig, ax1 = plt.subplots()

ax1.plot(a[0],a[1])
ax1.set_xlabel("Wavenumbers")
ax1.set_ylabel("Raman intensity")
ax1.set_xlim(1000,4000)

ax2 = ax1.twinx()   #creating another axis object sharing x axis, twiny() also possible
ax2.plot(b[0],b[1],color="r")
ax2.set_xlim(1000,4000)

"""
ax3 = ax1.twinx()
ax3.plot(c[0],c[1],color="black")
ax2.set_xlim(1000,4000)
"""
ax2.set_ylabel("Transmission")

#fig.tight_layout()
plt.show()

    