import os
import sqlite3

import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd

import SFG as sf
from SFG.orm.import_db_controller import WorkDatabaseWizard
from SFG.spectrum.sfg_spectrum import DummyPlotter

# p = pathlib.Path().cwd() / "SFG" / "mpl_config" / "origin.mpltstyle"
# plt.style.use(str(p))

mpl.rcParams["figure.subplot.hspace"] = 0.3

import numpy as np

db = sqlite3.connect("orm.db")
command = """SELECT 
date,
doy,
type, 
label,
sml_tension,
deep_tension,
sml_coverage,
deep_coverage,
sml_max_pressure,
deep_max_pressure,
sml_lift_off,
deep_lift_off,
sml_rawtension,
deep_rawtension
FROM stations
INNER JOIN station_stats
on stations.id = station_stats.station_id;
"""

command2 = """
SELECT
sample_hash AS 'name', 
coverage AS 'surface coverage',
max_pressure AS 'max. surface pressure',
lift_off AS 'lift-off compression ratio',
surface_tension AS 'surface tension'
FROM samples;
"""

command4 = """
SELECT 
sample_hash,
location,
type,
coverage AS 'surface coverage',
surface_tension AS 'surface tension',
max_pressure AS 'max. surface pressure',
lift_off AS ' lift-off compression ratio'
FROM samples
WHERE type='p'
OR type='s';
"""

command3 = """
SELECT 

date AS 'sampling date',
doy AS 'day of the year',

plate_coverage AS 'surface coverage plate',
plate_coverage_std,
plate_tension AS 'surface tension plate',
plate_tension_std,
plate_lift_off AS 'lift-off compression ratio plate',
plate_lift_off_std,
plate_max_pressure AS 'max. surface pressure plate',
plate_max_pressure_std, 

screen_coverage AS 'surface coverage screen',
screen_coverage_std,
screen_tension AS 'surface tension screen',
screen_tension_std,
screen_lift_off AS 'lift-off compression ratio screen',
screen_lift_off_std,
screen_max_pressure AS 'max. surface pressure screen',
screen_max_pressure_std, 

sml_coverage AS 'surface coverage SML',
sml_coverage_std,
sml_tension AS 'surface tension SML',
sml_tension_std,
sml_lift_off AS 'lift-off compression ratio SML',
sml_lift_off_std,
sml_max_pressure AS 'max. surface pressure SML',
sml_max_pressure_std,


deep_coverage AS 'surface coverage bulk',
deep_coverage_std,
deep_tension AS 'surface tension bulk',
deep_tension_std,
deep_lift_off AS 'lift-off compression ratio bulk',
deep_lift_off_std,
deep_max_pressure AS 'max. surface pressure bulk',
deep_max_pressure_std 

FROM stations
INNER JOIN station_stats
on stations.id = station_stats.station_id;
"""

command5 = """
SELECT
sampling_date AS 'sampling_date',
sml_coverage AS 'SML coverage',
bulk_coverage AS 'bulk coverage',
sml_no,
bulk_no
FROM be_data;

"""
"""
g = sf.gasex.GasExWorkDatabaseWizard()
samples = g.session.query(g.samples).all()
exsample = samples[191]

e_hash = exsample.id
sfg_name = g.session.query(g.gasex_sfg).filter(g.gasex_sfg.sample_id == e_hash).one().name
inter = g.session.query(g.sfg).filter(g.sfg.name == sfg_name).one()
sfg = g.construct_sfg(inter)

print(inter.measured_time.date())

lt_names = g.session.query(g.gasex_lt).filter(g.gasex_lt.sample_id == e_hash).all()
lts = [g.load_lt_by_name(i) for i in lt_names]

w = WorkDatabaseWizard()
w.origin_preview_date(surfacant="HA")
b = BEDatabaseWizard()
quartals = b.fetch_by_quartal(refine="sml", selection="t")
to_demo = [b.fetch_by_specid(i.specid) for i in quartals["t1"]]

test = to_demo[110]

lower = test.get_xrange_indices(np.min(test.x), 3030)
upper = test.get_xrange_indices(3030, np.max(test.x))

left = test.y[lower[0]:lower[1]+1]
right = (test.y[upper[0]+1:upper[1]+1])

base = peakutils.baseline(test.y[lower[0]:lower[1]+1], deg=1)
base2 = peakutils.baseline(test.y[upper[0]+1:upper[1]+1], deg=1)

peaks = peakutils.indexes(left-base, min_dist=3)

#plt.plot(test.x, test.y)
bs_corrected = np.concatenate([left-base, right-base2])
plt.axvline(test.x[lower[1]], color="red")
plt.axhline(0, color="blue")
for i in peaks:
    plt.axvline(test.x[i], color="green")
#plt.plot(test.x[lower[0]:lower[1]], left-base, color="black")
#plt.plot(test.x[upper[0]:upper[1]], right-base2, color="black")
plt.plot(test.x, bs_corrected)
plt.savefig("peak.png")
"""
def origin_preview_date():
    w = WorkDatabaseWizard()
    temp = w.session.query(w.regular_sfg).filter(w.regular_sfg.surfactant == "NA").all()
    temp = [w.session.query(w.sfg).filter(w.sfg.id == i.specid).one() for i in temp]
    dates = w.map_data_to_dates(temp)
    for key in dates:
        dir_name = str(key)
        os.mkdir(dir_name)
        sfg_spectra = [w.construct_sfg(i) for i in dates[key]]
        for spec in sfg_spectra: # type:import SFG.spectrum.sfg_spectrum
SFG.spectrum.sfg_spectrum.SfgSpectrum
            df = spec.convert_to_export_dataframe()
            df.to_csv(f'{dir_name}/'+spec.name + ".csv", index=False, sep=";")
        DummyPlotter(sfg_spectra, save=True, savedir=dir_name).plot_all()


def plot_sample(sfg, lts, name):
    fig, axs = plt.subplots(2)
    fig.suptitle(name)
    axs[0].plot(sfg.x, sfg.y, color="black")
    axs[0].set_xlabel(sfg.x_unit)
    axs[0].set_ylabel(sfg.y_unit.replace("/", "/\n"))

    for lt in lts:
        axs[1].plot(lt.x, lt.y, label=lt.name)
    axs[1].set_xlabel(lts[0].x_unit)
    axs[1].set_ylabel(lts[0].y_unit.replace("/", "/\n"))
    axs[1].legend()
    plt.savefig("test.png")


temp = pd.read_sql(command5, db)
temp = temp.replace([np.inf, -np.inf], np.nan)
temp = temp[(pd.notna(temp["SML coverage"]))]
temp["total"] = ((temp["SML coverage"] * temp["sml_no"]) + (temp["bulk coverage"] * temp["bulk_no"])) / (
            temp["sml_no"] + temp["bulk_no"])

mask = pd.isna(temp["total"])
test = temp["SML coverage"].loc[mask].copy(deep=True)
temp["total"].loc[mask] = test
temp.drop(labels=["bulk_no", "sml_no"], axis=1, inplace=True)
print(temp)
temp.to_csv("~/Schreibtisch/new_be_data.csv", sep=";", index=False)
