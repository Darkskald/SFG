import SFG as sf

import matplotlib.pyplot as plt
from matplotlib import gridspec
import os

from SFG.natural_samples import BEDatabaseWizard
from SFG.orm.orm import WorkDatabaseWizard
from SFG.spectrum.averagers import SfgAverager

dirname = os.path.dirname(__file__)
p = os.path.join(dirname, 'mpl_config/qt.mpltstyle')

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

months = [i for i in range(1, 13)]
years = [2009, 2010, 2011, 2012, 2013, 2014, 2015, 2017, 2018]


b = BEDatabaseWizard()
for year in years:
    t = b.query_by_year(year, "deep")

    fig = plt.figure(figsize=(11.69,8.27))
    spec = gridspec.GridSpec(ncols=4, nrows=3, figure=fig)

    month_counter = 1
    for row in range(3):
        for col in range(4):

            samples = b.refine_by_month(t, month_counter).all()
            spectra = [b.convert_be_to_sfg(i) for i in samples]

            if len(spectra) > 0:
                avg = SfgAverager(spectra, baseline=True)
                temp = fig.add_subplot(spec[row, col])
                temp.set_xlim(2700, 3900)
                for s in avg.spectra:
                    temp.plot(s.x, s.y, color="black", alpha=0.6, linewidth=0.5)
                temp.plot(avg.average_spectrum.x, avg.average_spectrum.y, color="red")
                temp.axhline(0, linestyle="--", color="black")
                temp.set_title(str(month_counter))

            month_counter += 1

    fig.suptitle(str(year))
    plt.tight_layout()
    plt.savefig(f'{year}_deep.png')





