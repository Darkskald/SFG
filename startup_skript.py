from ipy_interpreter import IpyInterpreter
from Classes import Analyzer, Plotter, TexInterface, SqlWizard, SqlExtractor
import os
I = IpyInterpreter()
I.get("su DPPC")


d = I.database.database

Sw = SqlWizard(d)
Sq = SqlExtractor("sfg.db")
print(os.getcwd())
s = Sq.fetch_single("single_sub_sfg",2)

P = Plotter ([s])
P.simple_plot()