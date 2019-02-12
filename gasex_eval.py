import csv
import time
import datetime
import copy
import traceback
import logging

# scientific libraries
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.lines import Line2D
from mpl_toolkits.mplot3d import Axes3D

import sqlite3
from scipy.signal import savgol_filter
from scipy.integrate import simps as sp
from scipy.interpolate import interp1d, splrep
from scipy.spatial import distance


def set_plot_properties(big=True):
    rcParams['xtick.labelsize'] = 12
    rcParams['ytick.labelsize'] = 12
    rcParams['legend.frameon'] = False

    rcParams['axes.labelsize'] = 14
    rcParams['axes.labelweight'] = 'normal'
    rcParams['axes.labelpad'] = 6.0

    rcParams['figure.autolayout'] = True
    rcParams['figure.figsize'] = 6.8, 4.25

    rcParams['lines.linewidth'] = 1
    rcParams['lines.markersize'] = 4

    rcParams['legend.fontsize'] = 10

    if big is False:
        rcParams['xtick.labelsize'] = 10
        rcParams['ytick.labelsize'] = 10

        rcParams['figure.figsize'] = 3.15, 1.96

        rcParams['lines.markersize'] = 4
        rcParams['lines.linewidth'] = 0.5

        rcParams['axes.labelsize'] = 10
        rcParams['axes.labelweight'] = 'normal'
        rcParams['figure.autolayout'] = True
        rcParams['legend.fontsize'] = 8


def scatter_maxpressure_day(isothermlist):
    """Create a scatter plot (day of cruise vs. maximum surface pressure) of the LtIsotherms
    provided in the isothermlist."""
    for isotherm in isothermlist:

        if isotherm.type == "p":
            color = "red"
        else:
            color = "blue"
        if "1" in isotherm.station:
            marker = "o"
        elif "2" in isotherm.station:
            marker = "x"
        elif "3" in isotherm.station:
            marker = "*"
        else:
            marker = "8"
        plt.scatter(isotherm.day, isotherm.get_maximum_pressure(), color=color, s=40, marker=marker)
        plt.text(isotherm.day + 0.1, isotherm.get_maximum_pressure() + 0.1,
                 isotherm.station + isotherm.type + isotherm.number,
                 fontsize=7)
        # plt.text(isotherm.day+0.2, isotherm.get_maximum_pressure()+0.1, isotherm.name, fontsize=7)

    plt.xlabel("days of cruise")
    plt.ylabel("maximum pressure/ mNm$^{-1}$")
    plt.grid()
    # legend_elements = [Line2D([0], [0], marker='o', color='r', label='glass plate', markerfacecolor='g', markersize=40)]
    # plt.legend(handles=legend_elements)
    plt.show()


def plot_vs_time(isothermlist):
    """Plots surface pressure vs time for each of the LtIsotherms provided in isothermlist."""
    for i in isothermlist:
        plt.plot(i.time, i.pressure)

    plt.xlabel("time")
    plt.ylabel("maximum pressure/ mNm$^{-1}$")
    plt.grid()
    plt.show()


def plot_per_sample(isothermlist):
    """Matches the items provided in the list of LtIsotherms according to their samples and exports all
     plots of a sample as pdf."""

    hashes = []
    for isotherm in isothermlist:

        if isotherm.create_sample_hash() not in hashes:
            hashes.append(isotherm.create_sample_hash())
            leg = isotherm.name
            plt.plot(isotherm.time, isotherm.pressure, label=leg)

            for partner in isotherm.partners:
                leg = partner.name
                plt.plot(partner.time, partner.pressure, label=leg)

            plt.xlabel("time/ s")
            plt.ylabel("maximum pressure/ mNm$^{-1}$")
            plt.title(str(isotherm.day) + " " + isotherm.station + " " + isotherm.type + " " + str(isotherm.number))
            plt.grid()
            plt.legend()
            plt.savefig(isotherm.create_sample_hash() + ".png")
            plt.cla()

        else:
            print("Already processed!")


def lt_plot_stats_new(S):
    """Takes Session controll manager as argument. Performs the bar/max surface pressure plot plot for the LT
    isotherms"""
    S.setup_for_gasex()

    for s in S.stations.values():
        s.join_samples()
        s.count_per_type()

    fig, (ax1, u_ax) = plt.subplots(nrows=2, ncols=1)
    ax1.set_xlabel("station number")
    ax1.set_ylabel("average surface pressure/ mN/m")
    ax2 = ax1.twinx()

    ax2.set_xlabel("day of cruise")
    ax2.set_ylabel("positive samples/ percent")

    ax1.set_title("Surfactant occurence in GasEx 1 (June '18), \nmeasured by Langmuir Trough\n\n\n")
    ax1.grid(True)

    days = [i for i in range(1, 15)]
    ax3 = u_ax.twiny()

    ax3.xaxis.set_ticks_position("bottom")
    ax3.xaxis.set_label_position("bottom")
    ax3.spines["bottom"].set_position(("axes", -0.35))
    ax3.set_xlabel("day of cruise")

    u_ax.set_xlabel("station number")
    u_ax.set_ylabel("Norm. SFG intensity/ arb.u.")
    u_ax.grid(True)
    u_ax.set_ylim(-0.0001, 0.0006)

    p_percentages = []
    s_percentages = []
    t_percentages = []

    stations = []
    plates = []
    totals = []
    screens = []

    p_std = []
    s_std = []
    t_std = []

    for s in S.stations.values():

        if len(s.lt_isotherms) > 0:
            station = s.station_number
            stations.append(station)
            screens.append(s.stats["screen_av"])
            plates.append(s.stats["plate_av"])
            totals.append(s.stats["total_av"])

            t_std.append(s.stats["std_total"])
            s_std.append(s.stats["std_screen"])
            p_std.append(s.stats["std_plate"])

            ax3.scatter(s.cruise_day, s.stats["screen_av"], s=0)
            s_percentages.append([station, s.stats["percent_screen"]])
            p_percentages.append([station, s.stats["percent_plate"]])
            t_percentages.append([station, s.stats["total_percent"]])

        max_ins = []

        for spectrum in s.sfg_spectra:  # type: SfgSpectrum

            station = s.station_number
            b = spectrum.slice_by_borders(3000, 2800)
            max = np.max(spectrum.normalized_intensity[b[0]:b[1] + 1])
            max_ins.append(max)

        max_ins = np.array(max_ins)
        av = np.average(max_ins)
        std = np.std(max_ins)
        u_ax.errorbar(station, av, std, alpha=0.9, fmt="o", color="b", label="average", barsabove="true", capsize=5,
                      capthick=2)

    ax1.errorbar(stations, totals, yerr=t_std, fmt="o", color="r", barsabove="true", capsize=5, capthick=2)
    rects2 = ax2.bar([a[0] + 0.2 for a in t_percentages], [a[1] for a in t_percentages], alpha=0.45, label="total")

    ax2.legend()
    # u_ax.legend()
    plt.tight_layout()
    plt.show()


def sfg_with_lt(S):
    """Takes Session controll manager as argument. Performs the bar/max surface pressure plot plot for the LT
        isotherms"""
    S.setup_for_gasex()

    for s in S.stations.values():
        s.join_samples()
        s.count_per_type()

    fig, ax1 = plt.subplots()
    ax1.set_xlabel("station number")
    ax1.set_ylabel("norm. SFG intensity/ arb. u.")
    ax2 = ax1.twinx()

    days = [i for i in range(1, 15)]
    ax3 = ax1.twiny()
    ax3.set_xlabel("day of cruise")
    ax3.xaxis.set_ticks_position("bottom")
    ax3.xaxis.set_label_position("bottom")
    ax3.spines["bottom"].set_position(("axes", -0.25))

    ax2.set_ylabel("positive samples/ percent")
    ax1.set_title("Surfactant occurence in GasEx 1 (June '18), \nmeasured by Langmuir Trough and SFG\n\n")
    ax1.grid(True)

    p_percentages = []
    s_percentages = []
    t_percentages = []

    for s in S.stations.values():

        station = s.station_number

        for spectrum in s.sfg_spectra:  # type: SfgSpectrum

            b = spectrum.slice_by_borders(3000, 2800)
            max = np.max(spectrum.normalized_intensity[b[0]:b[1] + 1])
            mapping = {"p": ["g", "plate"],
                       "s1": ["b", "screen"],
                       "s2": ["b", "screen"],
                       "s3": ["b", "screen"],
                       "deep": "black",
                       "low": "black"}
            try:
                ax1.scatter(station, max, alpha=0.9, color=mapping[spectrum.name.type][0],
                            label=mapping[spectrum.name.type][1],
                            s=60)
            except KeyError:
                ax1.scatter(station, max, alpha=0.9, color="black", label="CTD", s=60)

        ax3.scatter(s.cruise_day, max, s=0)
        ax1.set_ylim(0, 0.0008)
        # ax1.legend()

        s_percentages.append([station, s.stats["percent_screen"]])
        p_percentages.append([station, s.stats["percent_plate"]])
        t_percentages.append([station, s.stats["total_percent"]])

    rects1 = ax2.bar([a[0] - 0.2 for a in p_percentages], [a[1] for a in p_percentages], alpha=0.35, width=0.2,
                     color="g", label="plate, trough (positive)")
    rects2 = ax2.bar([a[0] for a in s_percentages], [a[1] for a in s_percentages], alpha=0.35, width=0.2, color="b",
                     label="screen, trough (positive)")
    rects2 = ax2.bar([a[0] + 0.2 for a in t_percentages], [a[1] for a in t_percentages], alpha=0.35, width=0.2,
                     color="r", label="total, trough (positive)")

    # ax2.legend()

    h1 = [Line2D([0], [0], marker='o', color='w', markerfacecolor='g', markersize=10),
          Line2D([0], [0], marker='o', color='w', markerfacecolor='b', markersize=10),
          Line2D([0], [0], marker='o', color='w', markerfacecolor='black', markersize=10)]

    l1 = ["plate, SFG", "screen, SFG", "CTD, SFG"]
    h2, l2 = ax2.get_legend_handles_labels()
    ax1.legend(h1 + h2, l1 + l2, loc="upper center", ncol=2).draggable()
    plt.sca(ax1)
    plt.tight_layout()
    plt.show()


def lt_sfg_integral(S):
    """Takes Session controll manager as argument. Performs the bar/max surface pressure plot plot for the LT
    isotherms"""

    S.setup_for_gasex()

    for s in S.stations.values():
        s.join_samples()
        s.count_per_type()

    fig, (ax1, u_ax) = plt.subplots(nrows=2, ncols=1)
    ax1.set_xlabel("station number")
    ax1.set_ylabel("average surface pressure/ mN/m")
    ax2 = ax1.twinx()

    ax2.set_xlabel("day of cruise")
    ax2.set_ylabel("positive samples/ percent")

    ax1.set_title("Surfactant occurrence in GasEx 1 (June '18)\n\n")
    ax1.grid(True)

    days = [i for i in range(1, 15)]
    ax3 = u_ax.twiny()

    ax3.xaxis.set_ticks_position("bottom")
    ax3.xaxis.set_label_position("bottom")
    ax3.spines["bottom"].set_position(("axes", -0.35))
    ax3.set_xlabel("day of cruise")

    u_ax.set_xlabel("station number")
    u_ax.set_ylabel("Integrated SFG CH intensity/ arb.u.")
    u_ax.grid(True)
    u_ax.set_ylim(-0.00025, 0.014)

    p_percentages = []
    s_percentages = []
    t_percentages = []

    stations = []
    plates = []
    totals = []
    screens = []

    p_std = []
    s_std = []
    t_std = []

    for s in S.stations.values():

        if len(s.lt_isotherms) > 0:
            station = s.station_number
            stations.append(station)
            screens.append(s.stats["screen_av"])
            plates.append(s.stats["plate_av"])
            totals.append(s.stats["total_av"])

            t_std.append(s.stats["std_total"])
            s_std.append(s.stats["std_screen"])
            p_std.append(s.stats["std_plate"])

            ax3.scatter(s.cruise_day, s.stats["screen_av"], s=0)
            s_percentages.append([station, s.stats["percent_screen"]])
            p_percentages.append([station, s.stats["percent_plate"]])
            t_percentages.append([station, s.stats["total_percent"]])

        max_ins = []

        if len(s.sfg_spectra) > 0:

            temp = []

            for spec in s.sfg_spectra:
                temp.append(spec.calculate_ch_integral())

            average = np.average(temp)
            std = np.std(temp)
            u_ax.errorbar(s.station_number, average, yerr=std, fmt="o", color="b", barsabove="true", capsize=5,
                          capthick=2)

    ax1.errorbar(stations, totals, yerr=t_std, fmt="o", color="r", barsabove="true", capsize=5, capthick=2)
    rects2 = ax2.bar([a[0] + 0.2 for a in t_percentages], [a[1] for a in t_percentages], alpha=0.45, label="positive")

    ax2.legend()
    # u_ax.legend()
    plt.tight_layout()
    plt.show()


def lt_sfg_integral_dppc(S):
    """Takes Session controll manager as argument. Performs the bar/max surface pressure plot plot for the LT
    isotherms"""

    S.setup_for_gasex()
    figure, plots = plt.subplots(3, 1)
    t_percentages = []

    ltplot = plots[0]
    # ltplot.set_title("Surfactant occurrence in GasEx 1 (June '18)\n\n", fontweight='bold')
    ltplot.set_ylabel("Average surface pressure/\n mN/m")
    ltplot.xaxis.set_ticklabels([])

    bars = ltplot.twinx()
    bars.set_ylabel("Positive samples/\n percent")

    sfg = plots[1]
    sfg.set_ylim(-0.00025, 0.014)
    sfg.set_ylabel("Integrated SFG \nCH intensity/ arb.u.")
    sfg.xaxis.set_ticklabels([])

    dppc = plots[2]
    dppc.set_ylabel("Surface coverage/\n %")
    dppc.set_xlabel("\nStation number")

    base = ltplot.twiny()
    base.xaxis.set_ticks_position("top")
    base.xaxis.set_label_position("top")
    # base.spines["top"].set_position(("axes", +0.7))
    base.set_xlabel("Day of cruise\n")

    for station in S.stations.values():  # type: Station

        if station.cruise_day < 13:
            base.scatter(station.cruise_day, station.stats["total_av"], s=0)

        if len(station.lt_isotherms) > 0:
            ltplot.scatter(station.station_number, station.stats["total_av"], color="green")
            ltplot.text(station.station_number + 0.15, station.stats["total_av"] + 0.15, str(station.isotherm_count))
            t_percentages.append([station.station_number, station.stats["total_percent"]])

        if len(station.sfg_spectra) > 0:
            if station.make_average_sfg() is not None:
                sfg.scatter(station.station_number, station.make_average_sfg(), color="blue")

            if station.make_average_sfg(dppc=S.dppc_ints) is not None:
                dppc.scatter(station.station_number, station.make_average_sfg(dppc=S.dppc_ints)[0] * 100, color="red")

    bars.bar([a[0] + 0.2 for a in t_percentages], [a[1] for a in t_percentages], alpha=0.45, label="positive")
    plt.show()


def lt_integral_average(S):
    """Takes Session controll manager as argument. Performs the bar/max surface pressure plot plot for the LT
    isotherms"""
    S.setup_for_gasex()

    for s in S.stations.values():
        s.join_samples()
        s.count_per_type()

    fig, (ax1, u_ax) = plt.subplots(nrows=2, ncols=1)
    ax1.set_xlabel("station number")
    ax1.set_ylabel("average surface pressure/ mN/m")
    ax2 = ax1.twinx()

    ax2.set_xlabel("day of cruise")
    ax2.set_ylabel("positive samples/ percent")

    ax1.set_title("Surfactant occurrence in GasEx 1 (June '18)\n\n")
    ax1.grid(True)

    days = [i for i in range(1, 15)]
    ax3 = u_ax.twiny()

    ax3.xaxis.set_ticks_position("bottom")
    ax3.xaxis.set_label_position("bottom")
    ax3.spines["bottom"].set_position(("axes", -0.35))
    ax3.set_xlabel("day of cruise")

    u_ax.set_xlabel("station number")
    u_ax.set_ylabel("Integrated SFG CH intensity/ arb.u.")
    u_ax.grid(True)
    u_ax.set_ylim(-0.00025, 0.014)

    p_percentages = []
    s_percentages = []
    t_percentages = []

    stations = []
    plates = []
    totals = []
    screens = []

    p_std = []
    s_std = []
    t_std = []

    for s in S.stations.values():  # type: Station

        if len(s.lt_isotherms) > 0:
            station = s.station_number
            stations.append(station)
            screens.append(s.stats["screen_av"])
            plates.append(s.stats["plate_av"])
            totals.append(s.stats["total_av"])

            t_std.append(s.stats["std_total"])
            s_std.append(s.stats["std_screen"])
            p_std.append(s.stats["std_plate"])

            ax3.scatter(s.cruise_day, s.stats["screen_av"], s=0)
            s_percentages.append([station, s.stats["percent_screen"]])
            p_percentages.append([station, s.stats["percent_plate"]])
            t_percentages.append([station, s.stats["total_percent"]])

        max_ins = []

        if len(s.sfg_spectra) > 0:

            av_spec = s.make_average_sfg()
            if isinstance(av_spec, AddedSpectrum):
                integral = av_spec.calculate_ch_integral()
                u_ax.scatter(s.station_number, integral)

    ax1.errorbar(stations, totals, yerr=t_std, fmt="o", color="r", barsabove="true", capsize=5, capthick=2)
    rects2 = ax2.bar([a[0] + 0.2 for a in t_percentages], [a[1] for a in t_percentages], alpha=0.45, label="positive")

    ax2.legend()
    # u_ax.legend()
    plt.tight_layout()
    plt.show()


def baseline_demo(spectrum, name="default"):
    spectrum.correct_baseline(average="gernot")

    test = np.linspace(2750, 3050, 10000)
    func = spectrum.make_ch_baseline(average="gernot")
    borders = spectrum.slice_by_borders(3000, np.min(spectrum.wavenumbers))

    f, axarr = plt.subplots(2, sharex=True)
    axarr[0].plot(spectrum.wavenumbers, spectrum.normalized_intensity, label=spectrum.name.full_name, linewidth=1.5,
                  marker="o", markersize=3)
    axarr[0].plot(test, func(test), color="r", label="Baseline")
    axarr[0].set_xlabel(str(spec.calculate_ch_integral(average="gernot")))
    axarr[0].set_ylabel("Norm. SFG intensity/ arb. u.")
    # axarr[0].set_title("Demonstration of the automatic baseline subtraction and integration")
    axarr[0].legend()

    axarr[1].plot(spectrum.wavenumbers, spectrum.baseline_corrected, label=spectrum.name.full_name, linewidth=1.5,
                  marker="o", markersize=3)
    axarr[1].fill_between(spectrum.wavenumbers[borders[0]:borders[1] + 1],
                          spectrum.baseline_corrected[borders[0]:borders[1] + 1])
    axarr[1].set_xlabel("Wavenumber/ cm$^{-1}$")
    axarr[1].set_ylabel("Norm. SFG intensity/ arb. u.")
    axarr[1].legend()
    name = spectrum.name.full_name
    plt.savefig(name + ".png")
    # plt.show()


def baseline_demo_dppc(spectrum, ref, name="default"):
    rcParams['xtick.labelsize'] = 18
    rcParams['ytick.labelsize'] = 18
    spectrum.correct_baseline(average="gernot")

    test = np.linspace(2750, 3050, 10000)
    func = spectrum.make_ch_baseline(average="gernot")
    borders = spectrum.slice_by_borders(3000, np.min(spectrum.wavenumbers))

    f, axarr = plt.subplots(3, sharex=True)

    axarr[0].plot(ref.wavenumbers, ref.normalized_intensity, label="reference", linewidth=1.5,
                  marker="o", markersize=3)

    axarr[0].legend(frameon=False)

    axarr[1].plot(spectrum.wavenumbers, spectrum.normalized_intensity, label="spectrum", linewidth=1.5,
                  marker="o", markersize=3)
    axarr[1].plot(test, func(test), color="r", label="baseline")

    axarr[1].legend(frameon=False)

    axarr[2].plot(spectrum.wavenumbers, spectrum.baseline_corrected, label="spectrum", linewidth=1.5,
                  marker="o", markersize=3)
    axarr[2].fill_between(spectrum.wavenumbers[borders[0]:borders[1] + 1],
                          spectrum.baseline_corrected[borders[0]:borders[1] + 1])
    axarr[2].set_xlabel("wavenumber/ cm$^{-1}$", fontsize=26)

    axarr[2].legend(frameon=False)
    f.text(0.025, 0.5, 'norm. intensity/ arb. u.', ha='center', va='center', rotation='vertical', fontsize=26)
    plt.savefig(spectrum.name.full_name + ".png")
    # plt.show()


def benchmark_baseline(speclist):
    for spectrum in speclist:  # type: SfgSpectrum

        standard = spectrum.calculate_ch_integral()
        regress = spectrum.calculate_ch_integral(average="min_reg")
        gernot = spectrum.calculate_ch_integral(average="gernot")
        return standard, regress, gernot


def plot_lt_isotherm(isotherm):  # type: LtIsotherm
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.set_xlabel("Total area (cm$^{-2}$)\n", fontsize=20)
    ax.set_ylabel("Surface pressure (mN/m)\n", fontsize=20)
    ax.set_ylim(-0.1, 2)
    ax.plot(isotherm.area, isotherm.pressure, linewidth=3)
    plt.show()


def broken_axis(x, y, lim):
    fig, (ax, ax2) = plt.subplots(1, 2, sharey=True)

    ax.set_ylabel("Surface tension/ $mN \cdot m^{-1}$")
    fig.text(0.5, 0.015, s="Day of the year", ha="center", va="center", size=14)

    ax.set_xlim(lim[0], lim[1])
    ax2.set_xlim(lim[2], lim[3])

    ax.spines['right'].set_visible(False)
    ax2.spines['left'].set_visible(False)
    ax.yaxis.tick_left()
    ax.tick_params(labeltop='off')  # don't put tick labels at the top
    ax2.yaxis.tick_right()

    d = .012  # how big to make the diagonal lines in axes coordinates
    # arguments to pass plot, just so we don't keep repeating them
    kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
    ax.plot((1 - d, 1 + d), (-d, +d), **kwargs)
    ax.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)

    kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
    ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)
    ax2.plot((-d, +d), (-d, +d), **kwargs)

    ax.scatter(x, y, marker="o")
    ax2.scatter(x, y, marker="o")


def broken_axis_errorbar(lim):
    fig, (ax, ax2) = plt.subplots(1, 2, sharey=True)

    ax.set_ylabel("Surface tension/ $mN \cdot m^{-1}$")
    fig.text(0.5, 0.015, s="day of the year", ha="center", va="center", size=16)

    ax.set_xlim(lim[0], lim[1])
    ax2.set_xlim(lim[2], lim[3])

    ax.spines['right'].set_visible(False)
    ax2.spines['left'].set_visible(False)
    ax.yaxis.tick_left()
    ax.tick_params(labeltop='off')  # don't put tick labels at the top
    ax2.yaxis.tick_right()

    d = .012  # how big to make the diagonal lines in axes coordinates
    # arguments to pass plot, just so we don't keep repeating them
    kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
    ax.plot((1 - d, 1 + d), (-d, +d), **kwargs)
    ax.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)

    kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
    ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)
    ax2.plot((-d, +d), (-d, +d), **kwargs)

    return ax, ax2


def plot_stats(stations, stats, ylabel="Surface tension/ $mN \cdot m^{-1}$"):
    axes = broken_axis_errorbar([153, 166, 254, 266])
    counter = 0
    colormap = {0: "ro", 1: "bo", 2: "go", 3: "purple"}
    axes[0].set_ylabel(ylabel)

    for s in stats:

        doy = []
        out = []
        err = []

        for station in stations:

            if station.stats[s] is not None:
                doy.append(station.get_doy())
                out.append(station.stats[s][0])
                err.append(station.stats[s][1])

            axes[0].errorbar(doy, out, yerr=err, fmt=colormap[counter], color="r", barsabove="true", capsize=5,
                             capthick=1, ecolor="black", elinewidth=1.0,
                             markeredgecolor="black", markeredgewidth=0.4, antialiased=True)
            axes[1].errorbar(doy, out, yerr=err, fmt=colormap[counter], color="r", barsabove="true", capsize=5,
                             capthick=1, ecolor="black", elinewidth=1.0,
                             markeredgecolor="black", markeredgewidth=0.4, antialiased=True, label=s)
            counter += 1

    axes[1].legend()

    plt.show()


def plot_stats_scatter(stations, stats, ylabel="Surface tension/ $mN \cdot m^{-1}$", scatter=False):
    rcParams['xtick.labelsize'] = 18
    rcParams['ytick.labelsize'] = 18
    rcParams['axes.labelsize'] = 'x-large'

    axes = broken_axis_errorbar([153, 166, 254, 266])
    counter = 0
    colormap = {0: "r", 1: "b", 2: "g", 3: "p"}
    legendmap = {"tension_deep": "bulkwater",
                 "tension_sml": "SML",
                 "coverage_sml": "SML",
                 "coverage_deep": "bulkwater",
                 "pressure_sml": "SML",
                 "pressure_deep": "bulkwater",
                 "liftoff_sml": "SML",
                 "liftoff_"
                 "deep": "bulkwater"}

    axes[0].set_ylabel(ylabel, fontsize=22)

    legend_elements = []
    labels = []

    for s in stats:

        av1 = []
        av2 = []

        for station in stations:

            if station.stats[s] is not None:
                d = station.get_doy()
                y = station.stats[s][0]
                axes[0].scatter(d, y, color=colormap[counter])
                axes[1].scatter(d, y, color=colormap[counter])

                if d < 165:
                    av1.append(y)
                else:
                    av2.append(y)
        av1 = np.average(av1)
        av2 = np.average(av2)
        axes[0].axhline(y=av1, color=colormap[counter], alpha=0.5, linestyle="--")
        axes[1].axhline(y=av2, color=colormap[counter], alpha=0.5, linestyle="--")

        labels.append(legendmap[s])
        labels.append(legendmap[s] + " average")
        counter += 1

    for i in range(counter):
        legend_elements.append(Line2D([0], [0], marker='o', color="w", markerfacecolor=colormap[i]))

        legend_elements.append(Line2D([0], [0], color=colormap[i],
                                      alpha=0.5, linestyle="--"))

    axes[1].legend(legend_elements, labels, frameon=False, )
    plt.show()


def plot_stats_scatter_small(stations, stats, ylabel="Surface tension/ $mN \cdot m^{-1}$", scatter=False):
    axes = broken_axis_errorbar([153, 166, 254, 266])
    counter = 0
    colormap = {0: "r", 1: "b", 2: "g", 3: "p"}
    legendmap = {"tension_deep": "bulkwater",
                 "tension_sml": "SML",
                 "coverage_sml": "SML",
                 "coverage_deep": "bulkwater",
                 "coverage_plate": "glass plate",
                 "coverage_screen": "screen",
                 "pressure_sml": "SML",
                 "pressure_deep": "bulkwater",
                 "liftoff_sml": "SML",
                 "liftoff_deep": "bulkwater"
                 }

    axes[0].set_ylabel(ylabel, fontsize=14)

    legend_elements = []
    labels = []

    for s in stats:

        av1 = []
        av2 = []

        for station in stations:

            if station.stats[s] is not None:
                d = station.get_doy()
                y = station.stats[s][0]
                error = station.stats[s][1]
                axes[0].errorbar(d, y, yerr=error, marker="o", mfc=colormap[counter],
                                 color="r", barsabove="true", capsize=5,
                                 capthick=1, ecolor="black", elinewidth=1.0,
                                 markeredgecolor="black", markeredgewidth=0.4,
                                 antialiased=True)

                axes[1].errorbar(d, y, yerr=error, marker="o", mfc=colormap[counter],
                                 color="r", barsabove="true", capsize=5,
                                 capthick=1, ecolor="black", elinewidth=1.0,
                                 markeredgecolor="black", markeredgewidth=0.4,
                                 antialiased=True)

                if d < 165:
                    av1.append(y)
                else:
                    av2.append(y)
        av1 = np.average(av1)
        av2 = np.average(av2)
        axes[0].axhline(y=av1, color=colormap[counter], alpha=0.5, linestyle="--")
        axes[1].axhline(y=av2, color=colormap[counter], alpha=0.5, linestyle="--")

        labels.append(legendmap[s])
        labels.append(legendmap[s] + " average")
        counter += 1

    for i in range(counter):
        legend_elements.append(Line2D([0], [0], marker='o', color="w", markerfacecolor=colormap[i]))

        legend_elements.append(Line2D([0], [0], color=colormap[i],
                                      alpha=0.5, linestyle="--"))

    axes[1].legend(legend_elements, labels, frameon=False, )
    plt.tight_layout()
    plt.show()


def tension_average(station):
    average = []
    x = station.get_doy()

    for tension in station.tensions:
        average.append(tension[1])

    average = np.array(average)
    av_out = np.average(average)
    std = np.std(average)

    return x, av_out, std


def correlation_plot(stations, value1, value2, spec1, spec2, average=True):
    dic = {
        "t": "surface tension/ $mN \cdot m^{-1}$",
        "s": "SFG CH integral/ arb u.",
        "c": "surface coverage",
        "p": "surface pressure/ $mN \cdot m^{-1}$"
    }
    if average is True:
        for station in stations:

            x = station.stats[value1]
            y = station.stats[value2]

            if x is not None and y is not None:
                plt.errorbar(x[0], y[0], xerr=x[1], yerr=y[1], fmt="ro",
                             color="r", barsabove="true", capsize=5,
                             capthick=1, ecolor="black", elinewidth=1.0,
                             markeredgecolor="black", markeredgewidth=0.4,
                             antialiased=True)

    else:
        for station in stations:
            plt.plot()

    plt.xlabel(dic[spec1])
    plt.ylabel(dic[spec2])
    plt.show()

from scm import SessionControlManager

S = SessionControlManager("sfg.db", "test")
S.setup_for_gasex()
set_plot_properties()

stats = [s for s in S.stations.values()]

small = [s for s in stats if s.type == "small"]
big = [s for s in stats if s.type == "big"]

import csv


def export_station_data(stations, big=True):
    if big:

        with open("big_stations.csv", "w") as outfile:
            writer = csv.writer(outfile, delimiter=";")
            header = "Station", "Date", "Type", "coverage_screen", "coverage_screen_error", "coverage_plate", \
                     "coverage_plate_error", "coverage_deep", "coverage_deep_error", "coverage_sml", \
                     "coverage_sml_error", "max_pressure_screen", "max_pressure_error", "lift_off_screen", \
                     "lift_off_screen_error", "max_pressure_plate", "max_pressure_plate", "lift_off_plate", \
                     "lift_off_screen_plate_error", "max_pressure_deep", "max_pressure_deep_error", "lift_off_deep", \
                     "lift_off_deep_error", "max_pressure_sml", "max_pressure_sml_error", "lift_off_sml", \
                     "lift_off_sml_error"

            writer.writerow(header)
            for st in stations:
                name = st.station_hash
                date = st.date
                try:
                    coverage_screen, coverage_screen_error, n = st.stats["coverage_screen"]
                except:
                    coverage_screen, coverage_screen_error = "None", "None"

                try:
                    coverage_plate, coverage_plate_error, n = st.stats["coverage_plate"]
                except:
                    coverage_plate, coverage_plate_error = "None", "None"

                try:
                    coverage_deep, coverage_deep_error, n = st.stats["coverage_deep"]
                except:
                    coverage_deep, coverage_deep_error = "None", "None"

                try:
                    coverage_sml, coverage_sml_error, n = st.stats["coverage_sml"]
                except:
                    coverage_sml, coverage_sml_error = "None", "None"

                #
                try:
                    pressure_screen, pressure_screen_error, n = st.stats["pressure_screen"]
                except:
                    pressure_screen, pressure_screen_error = "None", "None"

                try:
                    pressure_plate, pressure_plate_error, n = st.stats["pressure_plate"]
                except:
                    pressure_plate, pressure_plate_error = "None", "None"

                try:
                    pressure_deep, pressure_deep_error, n = st.stats["pressure_deep"]
                except:
                    pressure_deep, pressure_deep_error = "None", "None"

                try:
                    pressure_sml, pressure_sml_error, n = st.stats["pressure_sml"]
                except:
                    pressure_sml, pressure_sml_error = "None", "None"
                #
                try:
                    liftoff_screen, liftoff_screen_error, n = st.stats["liftoff_screen"]
                except:
                    liftoff_screen, liftoff_screen_error = "None", "None"

                try:
                    liftoff_plate, liftoff_plate_error, n = st.stats["liftoff_plate"]
                except:
                    liftoff_plate, liftoff_plate_error = "None", "None"

                try:
                    liftoff_deep, liftoff_deep_error, n = st.stats["liftoff"]
                except:
                    liftoff_deep, liftoff_deep_error = "None", "None"

                try:
                    liftoff_sml, liftoff_sml_error, n = st.stats["liftoff_sml"]
                except:
                    liftoff_sml, liftoff_sml_error = "None", "None"

                writer.writerow([name, date, str(st.type), coverage_screen, coverage_screen_error, coverage_plate,
                                coverage_plate_error, coverage_deep, coverage_deep_error, coverage_sml,
                                coverage_sml_error, pressure_screen, pressure_screen_error, liftoff_screen, liftoff_screen_error,
                                pressure_plate, pressure_plate_error, liftoff_plate, liftoff_plate_error, pressure_deep,
                                pressure_deep_error, liftoff_deep, liftoff_deep_error, pressure_sml, pressure_sml_error,
                                liftoff_sml, liftoff_sml_error])


with open("ltdata.csv", "w") as outfile:

    writer = csv.writer(outfile, delimiter=";")
    writer.writerow(["name", "station_type", "type", "max_surface_pressure", "lift_off",
                     "initial_area"])
    stats = sorted(stats)
    for s in stats:
        for isotherm in s.lt_isotherms:
            if isotherm.lift_off is not None:
                writer.writerow([isotherm.name, isotherm.sample_hash.station_type,
                                isotherm.sample_hash.sample_type, isotherm.get_maximum_pressure(),
                                isotherm.lift_off, np.max(isotherm.area)])

