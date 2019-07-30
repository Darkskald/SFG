import pandas as pd
import sqlite3

db = sqlite3.connect("orm.db")
command = """SELECT 
hash, 
type, 
date, 
label,
longitude,
latitude,
surface_salinity,
deep_salinity,
plate_tension,
screen_tension,
sml_tension,
deep_tension

FROM stations
INNER JOIN station_stats
on stations.id = station_stats.station_id;
"""
df = pd.read_sql(command, db).to_excel("tension_data.xlsx")

