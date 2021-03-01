from SFG.orm import interact
from SFG.natural_samples.gasex_processors import *

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

from plotly.subplots import make_subplots
import plotly.graph_objects as go


# Section for SQL queries

def setup_gasex(itc):
    # Sample fun

    stations = itc.session.query(itc.stations).all()
    samples = itc.session.query(itc.samples).all()

    sp = StationProcessor(stations, itc)
    samp = SampleProcessor(samples, itc)

    # fun with stations
    stations = sp.get_station_data_frame()

    # fun with samples
    samples = pd.DataFrame(samp.get_list_of_sample_dicts())
    stations_cleaned = stations[stations.columns.drop(list(stations.filter(regex='_(std)|(count$)')))]

    return stations_cleaned


def query_to_frame(query: str, itc: DbInteractor) -> pd.DataFrame:
    return pd.read_sql(query, itc.engine)


def query_to_shapiro(query: str, prop: str, itc: DbInteractor):
    """
    Get a record from the databsae using a specific SQL query string and test for normal
    distribution of the data using Shapiro-Wilkens test.

    :param prop: column name to retrieve from the output dataframe
    :param query: SQL query string to retrieve the relevant data
    :param itc: Database interactor
    :return: test statistics for the Shapiro-Wilk test for normal distribution
    """
    df = query_to_frame(query, itc)
    return stats.shapiro(df[prop])


itc = interact.DbInteractor()
# print(query_to_frame("SELECT * FROM sample_stats WHERE location='r';", itc))
# setup_gasex(itc).to_sql("station_values", itc.engine, index=False)
"""
query = "SELECT surface_tension FROM sample_stats WHERE type IN ('p', 's') AND surface_tension IS NOT NULL;"
query2 = "SELECT surface_tension FROM sample_stats WHERE type='deep' AND surface_tension IS NOT NULL"
print(query_to_shapiro(query, "surface_tension", itc))

print(query_mwu(query, query2, "surface_tension", "surface_tension", itc))
"""

ef_compare = ["tension", "coverage", "pressure"]


def calc_enrichment_by_prop_and_type(prop: str, sample_type: str, itc):
    q = f'SELECT {sample_type}_{prop}, d_{prop} FROM station_enrichment_test WHERE d_{prop} ' \
        f'IS NOT NULL AND {sample_type}_{prop} IS NOT NULL'
    df = pd.read_sql(q, itc.engine)
    print(q)
    return df[f'{sample_type}_{prop}'] / df[f'd_{prop}']


def calculate_enrichment_factors(itc: DbInteractor):
    indicators = ["tension", "coverage", "pressure"]
    result = {f'{i}_{j}': calc_enrichment_by_prop_and_type(i, j, itc) for j in ('p', 's') for i in indicators}
    return pd.DataFrame(result)


def plot_enrichment(itc: DbInteractor):
    df = pd.read_sql("SELECT * from enrichment_factors", itc.engine)
    fig = make_subplots(rows=3, cols=2)
    for row, indicator in enumerate(("tension", "pressure", "coverage")):
        if row == 0:
            fig.add_trace(go.Box(y=df[f'{indicator}_s'], marker_color="green",
                                 boxmean=True, name="screen " + indicator, legendgroup="screen", showlegend=False),
                          row=row + 1, col=1)
            fig.add_trace(go.Box(y=df[f'{indicator}_p'], marker_color="red",
                                 boxmean=True, name="plate " + indicator, legendgroup="plate", showlegend=False),
                          row=row + 1, col=2)

            fig.add_trace(go.Box(y=df[f'{indicator}_s'], marker_color="green",
                                 boxmean=True, name="screen", visible="legendonly"), row=row + 1, col=1)
            fig.add_trace(go.Box(y=df[f'{indicator}_p'], marker_color="red",
                                 boxmean=True, name="plate", visible="legendonly"), row=row + 1, col=2)
            # hack for marker and mean
            fig.add_trace(
                go.Scatter(x=[0], y=[-1], name="mean", mode="lines", line=dict(dash='dot', color='red'),
                           visible="legendonly"))
            fig.add_trace(
                go.Scatter(x=[0], y=[-1], name="median", mode="lines", line=dict(color='red'), visible="legendonly"))

        fig.add_trace(go.Box(y=df[f'{indicator}_s'], marker_color="green",
                             boxmean=True, showlegend=False, name="screen " + indicator, legendgroup="screen"),
                      row=row + 1, col=1)
        fig.add_trace(go.Box(y=df[f'{indicator}_p'], marker_color="red",
                             boxmean=True, showlegend=False, name="plate " + indicator, legendgroup="plate"),
                      row=row + 1, col=2)
        fig.update_yaxes(title_text="EF", row=row + 1, col=1)

    fig.update_layout(height=1900, width=3000, font_size=44, legend={'itemsizing': 'constant'})
    fig.show()


# calc_enrichment_by_prop_and_type("coverage", "p", itc)
# calculate_enrichment_factors(itc).to_sql("enrichment_factors", itc.engine)
# plot_enrichment(itc)
# df = pd.read_sql("SELECT * FROM enrichment_factors;", itc.engine)
# cols = [i for i in df if i != "index"]

# pairs = [f'{i}_{j}' for i in ('tension', 'pressure', 'coverage') for j in ('p', 's')]
# pairs = [pairs[:2], pairs[2:4], pairs[4:]]
# for p in pairs:
#    print(
#        f'{p[0]} to {p[1]}: {stats.ttest_ind(df[p[0]].replace("", np.nan).dropna(), df[p[1]].replace("", np.nan).dropna())}')

# categories


# Plate cruise 1
# Plate cruise 2
# Screen cruise 1
# Screen cruise 2
# Alkor cruise 1
# Alkor cruise 2
# rubber cruise 1
# rubber cruise 2

query_map = {
    "sml_cruise_1": "SELECT surface_tension, max_surface_pressure, lift_off_compression_ratio,coverage FROM sample_stats WHERE cruise=1 AND type IN ('s', 'p');",
    "sml_cruise_2": "SELECT surface_tension, max_surface_pressure, lift_off_compression_ratio,coverage FROM sample_stats WHERE cruise=2 AND type IN ('s', 'p');",
    "bulk_cruise_1": "SELECT surface_tension, max_surface_pressure, lift_off_compression_ratio,coverage FROM sample_stats WHERE cruise=1 AND type='deep';",
    "bulk_cruise_2": "SELECT surface_tension, max_surface_pressure, lift_off_compression_ratio,coverage FROM sample_stats WHERE cruise=2 AND type='deep';",
    "plate_cruise_1": "SELECT surface_tension, max_surface_pressure, lift_off_compression_ratio,coverage FROM sample_stats WHERE cruise=1 AND type='p';",
    "plate_cruise_2": "SELECT surface_tension, max_surface_pressure, lift_off_compression_ratio,coverage FROM sample_stats WHERE cruise=2 AND type='p';",
    "screen_cruise_1": "SELECT surface_tension, max_surface_pressure, lift_off_compression_ratio,coverage FROM sample_stats WHERE cruise=1 AND type='s';",
    "screen_cruise_2": "SELECT surface_tension, max_surface_pressure, lift_off_compression_ratio,coverage FROM sample_stats WHERE cruise=2 AND type='s';",
    "alkor_cruise_1": "SELECT surface_tension, max_surface_pressure, lift_off_compression_ratio,coverage FROM sample_stats WHERE cruise=1 AND location='a';",
    "alkor_cruise_2": "SELECT surface_tension, max_surface_pressure, lift_off_compression_ratio,coverage FROM sample_stats WHERE cruise=2 AND location='a';",
    "rubber_cruise_1": "SELECT surface_tension, max_surface_pressure, lift_off_compression_ratio,coverage FROM sample_stats WHERE cruise=1 AND location='r';",
    "rubber_cruise_2": "SELECT surface_tension, max_surface_pressure, lift_off_compression_ratio,coverage FROM sample_stats WHERE cruise=2 AND location='r';"
}

properties = ("surface_tension", "max_surface_pressure", "lift_off_compression_ratio", "coverage")

pairs = (("sml_cruise_1", "sml_cruise_2"),
         ("bulk_cruise_1", "bulk_cruise_2"),
         ("plate_cruise_1", "screen_cruise_1"),
         ("plate_cruise_2", "screen_cruise_2"),
         ("rubber_cruise_1", "alkor_cruise_1"),
         ("rubber_cruise_2", "alkor_cruise_2"),
         ("sml_cruise_1", "bulk_cruise_1"),
         ("sml_cruise_2", "bulk_cruise_2"),
         )
mwu_result = {"set_1": [], "set_2": [], "sai": [], "statistics": [], "p_value": [], "h0_accepted": []}


def query_mwu(query1: str, query2: str, prop1: str, prop2: str, itc: DbInteractor):
    v1 = query_to_frame(query1, itc)[prop1]
    v2 = query_to_frame(query2, itc)[prop2]
    return stats.mannwhitneyu(v1, v2)


for pair in pairs:
    print(pair)
    q1 = query_map[pair[0]]
    q2 = query_map[pair[1]]
    for prop in properties:
        try:
            temp = query_mwu(q1, q2, prop, prop, itc)
            mwu_result['set_1'].append(pair[0])
            mwu_result['set_2'].append(pair[1])
            mwu_result['sai'].append(prop)
            mwu_result['statistics'].append(temp[0])
            mwu_result['p_value'].append(temp[1])
            mwu_result['h0_accepted'].append('y' if temp[1] > 0.05 else 'n')
        except TypeError:
            mwu_result['set_1'].append(pair[0])
            mwu_result['set_2'].append(pair[1])
            mwu_result['sai'].append(prop)
            mwu_result['statistics'].append(None)
            mwu_result['p_value'].append(None)
            mwu_result['h0_accepted'].append(None)

pd.DataFrame(mwu_result).to_sql("mwu_results", itc.engine)


def get_shapiro_results(itc: DbInteractor, query_map):
    shapiro_results = {"name": [], "statistics": [], "p_value": []}
    for key in query_map:
        df = pd.read_sql(query_map[key], itc.engine)
        for col in df:
            try:
                temp = stats.shapiro(df[col].dropna())
                shapiro_results['name'].append(f'{key}_{col}')
                shapiro_results['statistics'].append(temp[0])
                shapiro_results['p_value'].append(temp[1])
            except:
                shapiro_results['name'].append(f'{key}_{col}')
            shapiro_results['statistics'].append(None)
            shapiro_results['p_value'].append(None)

    return shapiro_results
