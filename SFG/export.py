import pandas as pd
import sqlite3

db = sqlite3.connect("orm.db")
command = """SELECT 
date, 
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
df = pd.read_sql(command, db).to_excel("surfactant_data.xlsx")

