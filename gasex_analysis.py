from novel_gasex import GasexManager

# scientific libraries
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator
import scipy.stats
import seaborn as sns
plt.style.use(['seaborn-ticks', 'seaborn-poster'])


def baseline_demo_dppc(spectrum, ref, name="default"):

    spectrum.correct_baseline(average="gernot")

    test = np.linspace(2750, 3050, 10000)
    func = spectrum.make_ch_baseline(average="gernot")
    borders = spectrum.slice_by_borders(3000, np.min(spectrum.wavenumbers))


    f, axarr = plt.subplots(3, sharex=True)

    axarr[0].plot(ref.wavenumbers, ref.normalized_intensity, label="reference", linewidth=1.5,
                  marker="o", markersize=3)

    ref_borders = ref.slice_by_borders(3000, np.min(spectrum.wavenumbers))
    axarr[0].fill_between(ref.wavenumbers[ref_borders[0]:ref_borders[1] + 1],
                          ref.normalized_intensity[ref_borders[0]:ref_borders[1] + 1])

    axarr[0].legend(frameon=False)

    axarr[1].plot(spectrum.wavenumbers, spectrum.normalized_intensity, label="spectrum", linewidth=1.5,
                  marker="o", markersize=3)
    axarr[1].plot(test, func(test), color="r", label="baseline")

    axarr[1].legend(frameon=False)

    axarr[2].plot(spectrum.wavenumbers, spectrum.baseline_corrected, label="spectrum", linewidth=1.5,
                  marker="o", markersize=3)
    axarr[2].fill_between(spectrum.wavenumbers[borders[0]:borders[1] + 1],
                          spectrum.baseline_corrected[borders[0]:borders[1] + 1])
    axarr[2].set_xlabel("wavenumber/ cm$^{-1}$")

    axarr[2].legend(frameon=False)
    f.text(0.025, 0.5, 'norm. intensity/ arb. u.', ha='center', va='center', rotation='vertical', size=16)
    for ax in axarr.flat:
        ax.set_xlim(2750, 3300)
    #plt.savefig(spectrum.name.full_name + ".png")
    plt.show()


def broken_axis_errorbar(lim, label="default"):
    fig, (ax, ax2) = plt.subplots(1, 2, sharey=True)

    ax.set_ylabel(label)
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


def plot_doy_by_attribute(manager, names):

    labels = {
        "tension": "Surface tension/ $mN \cdot m^{-1}$",
        "pressure": "Max. surface pressure/ $mN \cdot m^{-1}$",
        "coverage": "% surface coverage",
        "lift_off": "lift-off compression ration"
    }
    for item in labels:
        if item in names[0]:
            axes = broken_axis_errorbar([153, 167, 254, 266], label=labels[item])

    if not axes:
        axes = broken_axis_errorbar([153, 167, 254, 266], label=labels[item])

    axes[1].xaxis.set_major_locator(MaxNLocator(integer=True))
    i = 1
    for name in names:
        axes[0].scatter(manager.station_table["date"].dt.dayofyear, manager.station_table[name])

        axes[1].scatter(manager.station_table["date"].dt.dayofyear, manager.station_table[name], label=name)
        i+=1


    axes[1].legend().draggable()
    plt.show()


if __name__ == "__main__":
    G = GasexManager("test.db", dppc_flag=True)
    plot_doy_by_attribute(G, ("sml_max_pressure", "deep_max_pressure"))
    # axes = broken_axis_errorbar([153, 167, 254, 266])
    # axes[1].xaxis.set_major_locator(MaxNLocator(integer=True))
    # axes[0].scatter(G.station_table["date"].dt.dayofyear, G.station_table["sml_tension"])
    # axes[1].scatter(G.station_table["date"].dt.dayofyear, G.station_table["sml_tension"],)
    # plt.show()
    #
    # spec = G.stations[1].samples[1].sfg_spectra[0]
    # q = G.fetch_dppc_spec(spec)
    # baseline_demo_dppc(q[0], q[1])
    # todo: -0.2155, -73.11
    # todo: stationwise difference between sml_coverage and deep_coverage
    # todo: what about lt isotherms for deep water in cruise 1?
    # todo: change aspect ration
    # todo: fourth row in sfg plot with doy vs coverage
    # todo: bechriftung größer bei theresa


    # mask1 = (G.station_table['date'] > '2018-6-1') & (G.station_table['date'] <= '2018-7-1')
    # cruise1 = G.station_table.loc[mask1]
    #
    # mask2 = (G.station_table['date'] > '2018-9-1') & (G.station_table['date'] <= '2018-10-1')
    # cruise2 = G.station_table.loc[mask2]
    #
    # cct = scipy.stats.ttest_ind(cruise1["sml_tension"].values, cruise2["sml_tension"].values,
    #                             equal_var=False, nan_policy='omit')
    # ccc = scipy.stats.ttest_ind(cruise1["sml_coverage"].values, cruise2["sml_coverage"].values,
    #                             equal_var=False, nan_policy='omit')
    # ccl = scipy.stats.ttest_ind(cruise1["sml_lift_off"].values, cruise2["sml_lift_off"].values,
    #                             equal_var=False, nan_policy='omit')
    # ccp = scipy.stats.ttest_ind(cruise1["sml_max_pressure"].values, cruise2["sml_max_pressure"].values,
    #                             equal_var=False, nan_policy='omit')
    #
    # dcct = scipy.stats.ttest_ind(cruise1["deep_tension"].values, cruise2["deep_tension"].values,
    #                              equal_var=False, nan_policy='omit')
    # dccc = scipy.stats.ttest_ind(cruise1["deep_coverage"].values, cruise2["deep_coverage"].values,
    #                              equal_var=False, nan_policy='omit')
    # dccl = scipy.stats.ttest_ind(cruise1["deep_lift_off"].values, cruise2["deep_lift_off"].values,
    #                              equal_var=False, nan_policy='omit')
    # dccp = scipy.stats.ttest_ind(cruise1["deep_max_pressure"].values, cruise2["deep_max_pressure"].values,
    #                              equal_var=False, nan_policy='omit')
    #
    # print(f'sml tension c1c2: {cct}')
    # print(f'sml coverage c1c2: {ccc}')
    # print(f'sml lift_off c1c2: {ccl}')
    # print(f'sml pressure c1c2: {ccp}')
    # print(f'deep tension c1c2: {dcct}')
    # print(f'deep coverage c1c2: {dccc}')
    # print(f'deep lift_off c1c2: {dccl}')
    # print(f'deep pressure c1c2: {dccp}')
    #
    # small = G.station_table[G.station_table["type"] == "small"]
    # big = G.station_table[G.station_table["type"] == "big"]
    #
    # bst = scipy.stats.ttest_ind(big["sml_tension"].values, small["sml_tension"].values,
    #                             equal_var=False, nan_policy='omit')
    #
    # bsc = scipy.stats.ttest_ind(big["sml_coverage"].values, small["sml_coverage"].values,
    #                             equal_var=False, nan_policy='omit')
    #
    # bsl = scipy.stats.ttest_ind(big["sml_lift_off"].values, small["sml_lift_off"].values,
    #                             equal_var=False, nan_policy='omit')
    #
    # bsp = scipy.stats.ttest_ind(big["sml_max_pressure"].values, small["sml_max_pressure"].values,
    #                             equal_var=False, nan_policy='omit')
    #
    # print(f'sml tension bs: {bst}')
    # print(f'sml coverage bs: {bsc}')
    # print(f'sml lift_off bs: {bsl}')
    # print(f'sml pressure bs: {bsp}')
