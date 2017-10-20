import csv
import matplotlib.pyplot as plt

class UV_Spectrum:
    
    def __init__ (self, name, wavelengths, intensities):
        
        self.name = name
        self.wavelengths = wavelengths
        self.intensities = intensities

def get_data(filename):
    with open(filename,"r") as infile:
        next(infile)
        next(infile)
        c = csv.reader(infile,delimiter="\t")
        wavelength = []
        intensity = []
        for line in c:
             l= line[0].split(",") 
             if len(l) == 2:
                 print(l)
                 wavelength.append(float(l[0]))
                 intensity.append(float(l[1]))
    return (wavelength,intensity)
    
a = get_data("BX6_2.csv")
b = get_data("BX9_1.csv")
c = get_data("BX12.csv")

plt.plot(a[0],a[1],label="BX6")
plt.plot(b[0],b[1],label="BX9")
plt.plot(c[0],c[1],label="BX12")
plt.legend()
plt.show()