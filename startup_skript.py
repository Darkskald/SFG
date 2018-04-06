from ipy_interpreter import IpyInterpreter
from Classes import Analyzer, Plotter, TexInterface

I = IpyInterpreter()
I.get("su DPPC")

s = I.subset[4]
P = Plotter([s])
P.marked_peaks()

