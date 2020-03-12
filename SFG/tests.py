import SFG as sf

import matplotlib.pyplot as plt
from matplotlib import gridspec
import os

from SFG.natural_samples import BEDatabaseWizard
from SFG.orm.orm import WorkDatabaseWizard
from SFG.spectrum.averagers import SfgAverager

dirname = os.path.dirname(__file__)
p = os.path.join(dirname, 'mpl_config/origin.mpltstyle')

plt.style.use(p)

# w = sf.orm.WorkDatabaseWizard()
"""
temp = w.session.query(w.sfg).filter(w.sfg.type =="gasex_sfg").all()
specs = [w.construct_sfg(i) for i in temp]

for spec in specs:
     plt.plot(spec.x, spec.y)
     plt.xlabel(spec.x_unit)
     plt.ylabel(spec.y_unit)
     plt.title(spec.name)
     plt.savefig(spec.name+ ".png")
     plt.cla()
"""
#specs = [w.fetch_by_specid(i) for i in [1256, 1107]]
"""
for spec in specs:
     p = spec.convert_to_export_dataframe()
     p.to_csv(spec.name + ".csv", index=False, sep=";")

import peakutils
w = WorkDatabaseWizard()
testspec = w.construct_sfg(w.session.query(w.sfg).all()[1792])
base = peakutils.baseline(testspec.y, deg=3)
plt.plot(testspec.x, testspec.y)
plt.plot(testspec.x, base)
plt.show()
"""
#for spec in specs:
     #sf.plotting.baseline_demo_dppc(spec)
b = BEDatabaseWizard()
trimesters = b.fetch_by_quartal(selection="t", refine="sml")

for t in trimesters:
    spectra = b.normalize_to_reference_integral(trimesters[t])
    avg = SfgAverager(spectra, baseline=True, enforce_scale=True).average_spectrum
    plt.plot(avg.x, avg.y, label=t)
plt.legend()
plt.savefig("trimesters_base_extra.png")







