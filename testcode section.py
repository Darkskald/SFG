import numpy as np


def integrate_peak(x_array, y_array):
    area = 0
    for i in range(len(x_array) - 1):
        dx = abs(x_array[i + 1] - x_array[i])
        square = dx * y_array[i]
        triangle = dx * abs(y_array[i + 1] - y_array[i]) * 0.5
        total = square + triangle
        area += total
    return area

x = np.linspace(0,5,10000)
y = (2*x**2)+(0.4*x)-3
y2 = (6*x**3)+(2*x**2)+(0.4*x)-3
y3 = 0.5*np.sin(x)+2*x

als = integrate_peak(x,y)
aln = np.trapz(y,x)

print("My Integrate")
print(als)
print("\n")

print("Numpy")
print(aln)
print("\n")

a2i = integrate_peak(x,y2)
a3i = integrate_peak(x,y3)
a2n = integrate_peak(y2,x)
a3n = integrate_peak(y3,x)
print("My Integrate")
print(a2i)
print("\n")

print("Numpy")
print(a2n)
print("\n")

print("My Integrate")
print(a3i)
print("\n")

print("Numpy")
print(a3n)
print("\n")