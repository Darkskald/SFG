from datetime import timedelta, date

from sqlalchemy import func

from SFG.natural_samples.gasex_processors import StationProcessor, SampleProcessor, SamplePlotProcessor, \
    SampleLatexProcessor, first_split_then_cruise, apply_ttest_along_properties, apply_ttest_by_category_and_cruise, \
    gather_sml
from SFG.orm.interact import DbInteractor

import pandas as pd

from plotly.subplots import make_subplots
import plotly.graph_objects as go
import plotly.express as px

interactor = DbInteractor()
stations = interactor.session.query(interactor.stations).all()
samples = interactor.session.query(interactor.samples).all()
stp = StationProcessor(stations, interactor)
samp = SampleProcessor(samples, interactor)
smps = pd.DataFrame(samp.get_list_of_sample_dicts())

sps = SamplePlotProcessor(smps)
#

#sml, deep, _, _ = sps.split_dataset("type", ('sml', 'deep'))


#test = apply_ttest_by_category_and_cruise(smps, 'type', ('sml', 'deep'))
#print(test)
#print(sml[sml["cruise"] == 1])
# print(sps.compare_columns_by_ttest(sml["max_surface_pressure"], deep["max_surface_pressure"]))
#test = [{"category": "surace_tension", "tstat": -1.34, "pval": 1.2, "equal": "accepted"}]


def plot_sample_by_plotly(dic):
    fig = make_subplots(rows=2, cols=2)
    coord = [(1, 1), (1, 2), (2, 1), (2, 2)]
    count = 0

    for key in dic:
        row, col = coord[count]
        if count == 0:
            fig.add_trace(go.Box(**dic[key][0], marker_color="red", boxmean=True), row=row, col=col)
            fig.add_trace(go.Box(**dic[key][1], marker_color="green", boxmean=True), row=row, col=col)
            # hack for marker and mean
            fig.add_trace(
                go.Scatter(x=[0], y=[-1], name="mean", mode="lines", line=dict(dash='dot', color='red'), visible="legendonly"))
            fig.add_trace(
                go.Scatter(x=[0], y=[-1], name="median", mode="lines", line=dict(color='red'), visible="legendonly"))
            # fig.add_trace(go.Scatter(name="mean"))
        else:
            fig.add_trace(go.Box(**dic[key][0], marker_color="red", boxmean=True, showlegend=False), row=row, col=col)
            fig.add_trace(go.Box(**dic[key][1], marker_color="green", boxmean=True, showlegend=False), row=row, col=col)
        fig.add_trace(go.Box(**dic[key][2], marker_color="red", boxmean=True, showlegend=False), row=row, col=col)
        fig.add_trace(go.Box(**dic[key][3], marker_color="green", boxmean=True, showlegend=False), row=row, col=col)
        fig.update_yaxes(title_text=key, row=row, col=col)
        count += 1

    fig.update_layout(height=1000, width=1400, font_size=18)
    fig.show()


# this code is used to produce box plots
#smps["type"] = gather_sml(smps)
out = sps.unwrap_properties_for_plotting(*sps.split_dataset("location", ('a', 'r')))
plot_sample_by_plotly(out)
