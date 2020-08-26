from datetime import timedelta

from sqlalchemy import func

from SFG.natural_samples.gasex_processors import StationProcessor, SampleProcessor, SamplePlotProcessor
from SFG.orm.interact import DbInteractor

import pandas as pd

from plotly.subplots import make_subplots
import plotly.graph_objects as go

interactor = DbInteractor()
stations = interactor.session.query(interactor.stations).all()
samples = interactor.session.query(interactor.samples).all()
samp = SampleProcessor(samples, interactor)
smps = pd.DataFrame(samp.get_list_of_sample_dicts())

sps = SamplePlotProcessor(smps)
sps.df["type"] = sps.gather_sml()
out = sps.unwrap_properties_for_plotting(*sps.split_dataset("type", ('sml', 'deep')))



def plot_sample_by_plotly(dic):
    fig = make_subplots(rows=2, cols=2)
    coord = [(1, 1), (1, 2), (2, 1), (2, 2)]
    count = 0

    for key in dic:
        row, col = coord[count]
        if count == 0:
            fig.add_trace(go.Box(**dic[key][0], marker_color="red", boxmean=True), row=row, col=col)
            fig.add_trace(go.Box(**dic[key][1], marker_color="green", boxmean=True), row=row, col=col)
        else:
            fig.add_trace(go.Box(**dic[key][0], marker_color="red", boxmean=True, showlegend=False), row=row, col=col)
            fig.add_trace(go.Box(**dic[key][1], marker_color="green", boxmean=True, showlegend=False), row=row, col=col)
        fig.add_trace(go.Box(**dic[key][2], marker_color="red", boxmean=True, showlegend=False), row=row, col=col)
        fig.add_trace(go.Box(**dic[key][3], marker_color="green", boxmean=True, showlegend=False), row=row, col=col)
        fig.update_yaxes(title_text=key, row=row, col=col)
        count += 1

    fig.update_layout(height=1000, width=1400, font_size=18)
    fig.show()
plot_sample_by_plotly(out)
"""
stations = interactor.session.query(interactor.stations).all()
samples = interactor.session.query(interactor.samples).all()

s = SampleProcessor(samples, interactor)

for sa in s.samples:
    print(s.get_corrected_salinity(sa))

for s in temp:
    if 0 <= s.measured_time.hour < 8:
        s.measured_time -= timedelta(days=1)
    temp = interactor.session.query(interactor.measurement_days).filter(func.DATE(s.measured_time) == interactor.measurement_days.date).all()
    print(temp)
temp = interactor.session.query(interactor.sfg).all()

stations = interactor.session.query(interactor.stations).all()
for s in stations:
    print(s.get_doy())
    
from SFG.orm import interact
from SFG.natural_samples.gasex_processors import StationProcessor
import itertools as ito

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


interactor = interact.DbInteractor()
stations = interactor.session.query(interactor.stations).all()
sp = StationProcessor(stations, interactor)
df = sp.get_station_data_frame()


"""
