from ipy_interpreter import IpyInterpreter
from Classes import Analyzer, Plotter

I = IpyInterpreter()
I.get("su PA")
A = Analyzer(I.subset)
P = Plotter(I.subset)
P.bar_peaks()
