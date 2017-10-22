from Classes import *
from ipy_interpreter import Ipy_Interpreter
import numpy as np
import matplotlib.pyplot as plt

colors = {
    0 : "b",
    1 : "g",
    2: "c",
    3: "m",
    4 : "y",
    5 : "k",
    }


def detailed_analysis(x_array,y_array):
    x_array = x_array[::-1]
    y_array = y_array[::-1]
    
    slopes = [(y_array[i+1]-y_array[i])/(x_array[i+1]-x_array[i]) for i in range((len(x_array)-1))]
    possible_peaks = []
    peak_tuples = []
    
    
    for i in range(1,len(y_array)-1):
        if slopes[i-1] > 0 and slopes[i] < 0:
            possible_peaks.append(i)
            
           
            
    average_intensity = np.average(y_array)
    
    confirmed_peaks = [i for i in possible_peaks if (y_array[i] > average_intensity*5)]
  
    for i in confirmed_peaks:
        
  
        left = 0
        right = 0
        center = i
        k = i-2
        

        # check for left border
        while slopes[k] > 0:
            k -= 1
        left = k
        
        # check for right border
        k = i+1
  
        while ((slopes[k] < 0) and (k+1 <= len(slopes)-1)):
             k += 1
             
        right = k
        
        
        peak_tuples.append((center,left,right))
     
    data_out = []
    for i in peak_tuples:
        indices = (i[0],i[1],i[2])
        print("center: "+str(x_array[i[0]]))
        center = x_array[i[0]]
        print("left: "+str(x_array[i[1]]))
        left = x_array[i[1]]
        print("right: "+str(x_array[i[2]]))
        right = x_array[i[2]]
        peak_slice_x = x_array[i[1]:i[2]+1]
        peak_slice_y = y_array[i[1]:i[2]+1]
        area = integrate_peak(peak_slice_x,peak_slice_y)
        datapoints = len(peak_slice_x)
        
        data_out.append((center,left,right,peak_slice_x,peak_slice_y,datapoints,area,indices))
    return data_out

def integrate_peak(x_array,y_array):
 
    area = 0
    for i in range(len(x_array)-1):
        dx = abs(x_array[i+1]-x_array[i])
        print(dx)
        square = dx*y_array[i]
        triangle = dx*abs(y_array[i+1]-y_array[i])*0.5
        total = square + triangle
        area += total
    return area    
    
#test code section
"""
i = Ipy_Interpreter()
i.get("su DPPC")
i.keep("80")
i.subset[0].smooth()
a = i.subset[0].wavenumbers
b = i.subset[0].normalized_intensity
q=detailed_analysis(a,b)
print(q)
i.plot()


i = Ipy_Interpreter()
i.get("su DPPC")

for j in range(len(i.subset)):
    i.subset[j].smooth()
    a = i.subset[j].wavenumbers
    b = i.subset[j].normalized_intensity
    try:
        q=detailed_analysis(a,b)
        print(q)
    except LookupError:
        print("Something happened with indexing")  
x = np.linspace(-0.5,3.3,500)
y = 0.8*np.sin(2*x)**2

print("Hello")
"""
i = Ipy_Interpreter()
i.get("su DPPC")
i.refine("d 20170914")
sfg = i.subset[0]

x= sfg.wavenumbers
y=sfg.normalized_intensity

analyse = detailed_analysis(x,y)

x= sfg.wavenumbers[::-1]
y=sfg.normalized_intensity[::-1]


fig = plt.figure()
ax = plt.subplot(111)
ax.plot(x,y)
ax.set_ylim(0)
j = 0
for i in analyse:
  
    print("Test:")
    print(i)
    center = i[0]
    left = i[1]
    right = i[2]
    area = i[6]
    indices = i[7]
    print(indices)
    centerx = x[indices[0]]
    centery = y[indices[0]]
    print(centerx,center)

    ax.fill_between(i[3],i[4])
    ax.axvline(left,color=colors[j],linewidth=2)
    ax.axvline(right,color=colors[j],linewidth=2)

    ax.annotate((float("{0:.2f}".format(centerx))),(centerx, centery)
   
                )
    ax.annotate("area: "+str((float("{0:.2f}".format(area)))),(centerx-0.2,centery*0.3))
    j +=1
   

plt.show()
