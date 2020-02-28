import pandas as pd
import matplotlib.pyplot as plt
#plt.style.use("mpl_config/origin.mpltstyle")
import sqlite3

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
#df = pd.read_sql(command, db).to_excel("surfactant_data.xlsx")

df = pd.read_sql(command2, db)

q = df.corr()
print(q)


