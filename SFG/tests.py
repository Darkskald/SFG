import SFG as sf

import matplotlib.pyplot as plt
import os

dirname = os.path.dirname(__file__)
p = os.path.join(dirname, 'mpl_config/origin.mpltstyle')

plt.style.use(p)

w = sf.orm.WorkDatabaseWizard()
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
specs = [w.fetch_by_specid(i) for i in [1256, 1107]]
"""
for spec in specs:
     p = spec.convert_to_export_dataframe()
     p.to_csv(spec.name + ".csv", index=False, sep=";")
"""
for spec in specs:
     sf.plotting.baseline_demo_dppc(spec)


