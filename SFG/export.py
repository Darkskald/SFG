import SFG as sf

import pandas as pd
import matplotlib.pyplot as plt
import pathlib
import sqlite3

from SFG.orm.orm import WorkDatabaseWizard

p = pathlib.Path().cwd() / "SFG" / "mpl_config" / "origin.mpltstyle"
plt.style.use(str(p))

import matplotlib as mpl

mpl.rcParams["figure.subplot.hspace"] = 0.3

import numpy as np
from scipy import stats

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
"""

w = WorkDatabaseWizard()
temp = w.session.query(w.sfg).all()

print(w.map_data_to_dates(temp))

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



#