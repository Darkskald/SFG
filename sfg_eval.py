from Classes import SessionControlManager as scm

# scientific libraries
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.lines import Line2D

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



def sfg_plot_broken_axis(speclist, lower, upper, title="default", normalized="false"):
    """Produces a pre-formatted SFG plot from a list of SFG spectrum objects"""

    fig, (ax, ax2) = plt.subplots(1, 2, sharey=True)
    #ax.set_xlabel("Wavenumber/ cm$^{-1}$", fontsize=10)
    #ax.set_ylabel("Norm. SFG intensity/ arb. u.", fontsize=10)

    ax.set_xlim(lower[0], lower[1])
    ax2.set_xlim(upper[0], upper[1])
    fig.text(0.025, 0.5, 'Norm. intensity/ arb. u.', ha='center', va='center', rotation='vertical')
    fig.text(0.5, 0.025, 'Wavenumber/ cm$^{-1}$', ha='center', va='center')

    ax.spines['right'].set_visible(False)
    ax2.spines['left'].set_visible(False)
    ax.yaxis.tick_left()
    ax.tick_params(labeltop='off')  # don't put tick labels at the top
    ax2.yaxis.tick_right()

    d = .01  # how big to make the diagonal lines in axes coordinates
    # arguments to pass plot, just so we don't keep repeating them
    kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
    ax.plot((1 - d, 1 + d), (-d, +d), **kwargs)
    ax.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)

    kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
    ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)
    ax2.plot((-d, +d), (-d, +d), **kwargs)

    inc = 0.25 / len(speclist)
    counter = 0
    for spectrum in speclist:
        eff_alpha = 0.75 + inc * counter
        if normalized == "false":
            ax.plot(spectrum.wavenumbers, spectrum.normalized_intensity, linewidth=1.5, marker="s",
                    alpha=eff_alpha, label=spectrum.name.full_name)
            ax2.plot(spectrum.wavenumbers, spectrum.normalized_intensity, linewidth=1.5, marker="s",
                     alpha=eff_alpha, label=spectrum.name.full_name)
        elif normalized == "true":
            ax.plot(spectrum.wavenumbers, spectrum.normalize_to_highest(), linewidth=1.5, marker="s", markersize=4,
                    alpha=eff_alpha, label=spectrum.name.full_name)
        counter += 1

    ax2.legend()
    return fig


def sfg_pub_plot(speclist, title="default", normalized="false"):
    """Produces a pre-formatted SFG plot from a list of SFG spectrum objects"""

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.set_xlabel("Wavenumber/ cm$^{-1}$")
    ax.set_ylabel("Norm. SFG intensity/ arb. u.")

    inc = 0.25 / len(speclist)
    counter = 0
    for spectrum in speclist:
        eff_alpha = 0.75 + inc * counter

        if normalized == "false":
            ax.plot(spectrum.wavenumbers, spectrum.normalized_intensity, linewidth=1.5, marker="s",
                    alpha=eff_alpha, label=spectrum.name.full_name)
        elif normalized == "true":

            ax.plot(spectrum.wavenumbers, spectrum.normalize_to_highest(), linewidth=1.5, marker="s",
                    alpha=eff_alpha, label=spectrum.name.full_name)

        counter += 1
    ax.legend(frameon=False)
    #fig.tight_layout()
    return fig


def sfg_stack_plot(speclist):
    """Stacks the provided SfgSpectrum objects after normalizing them to 1 internally."""
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.set_xlabel("Wavenumber/ cm$^{-1}$", fontsize=10)
    ax.set_ylabel("Norm. SFG intensity/ arb. u.", fontsize=10)
    ax.set_yticks([])
    ax.set_yticklabels([])

    inc = 0.25 / len(speclist)
    counter = 0
    offset = 0
    for spectrum in speclist:
        eff_alpha = 0.75 + inc * counter
        ax.plot(spectrum.wavenumbers, spectrum.normalize_to_highest() + offset, linewidth=1.5, marker="o", markersize=3,
                alpha=eff_alpha, label=spectrum.name.full_name)
        counter += 1
        offset += 0.2

    #size = fig.get_size_inches()
    #ratio = size[0] / size[1]
    #fig.set_size_inches(3.2 * ratio, 3.2)
    #fig.tight_layout()
    return fig


def sfg_doublestack_plot(speclist1, speclist2):
    """Two stacked-by-offset plots, one on the bottom, on on the top of the figure."""

    fig = plt.figure()
    ax = fig.add_subplot(111) # for label only
    ax1 = fig.add_subplot(211)
    ax2 = fig.add_subplot(212)
    ax2.set_xlabel("Wavenumber/ cm$^{-1}$", fontsize=22)
    ax.set_ylabel("Norm. SFG intensity/ arb. u.", fontsize=22)


    ax.spines['top'].set_color('none')
    ax.spines['bottom'].set_color('none')
    ax.spines['left'].set_color('none')
    ax.spines['right'].set_color('none')
    ax.tick_params(labelcolor='w', top=False, bottom=False, left=False, right=False)

    ax2.set_yticks([])
    ax2.set_yticklabels([])
    ax1.set_yticks([])
    ax1.set_yticklabels([])

    inc = 0.25 / len(speclist1)
    counter = 0
    offset = 0
    for spectrum in speclist1:
        eff_alpha = 0.75 + inc * counter
        ax1.plot(spectrum.wavenumbers, spectrum.normalize_to_highest() + offset, linewidth=1.5, marker="o",
                 markersize=3,
                 alpha=eff_alpha, label=spectrum.name.full_name)
        counter += 1
        offset += 0.4

    inc = 0.25 / len(speclist1)
    counter = 0
    offset = 0
    for spectrum in speclist2:
        eff_alpha = 0.75 + inc * counter
        ax2.plot(spectrum.wavenumbers, spectrum.normalize_to_highest() + offset, linewidth=1.5, marker="o",
                 markersize=3,
                 alpha=eff_alpha, label=spectrum.name.full_name)
        counter += 1
        offset += 0.2
    ax1.legend()
    ax2.legend()
    size = fig.get_size_inches()
    ratio = size[0] / size[1]
    fig.set_size_inches(3.2 * ratio, 3.2)
    fig.tight_layout()
    return fig


def sfg_doublestack_plot_broken(speclist1, speclist2, lower, upper, labels1=None, labels2=None):
    """Two stacked-by-offset plots, one on the bottom, on on the top of the figure."""

    fig, (row_1, row_2) = plt.subplots(nrows=2, ncols=2, sharey=True)

    left_top = row_1[0]
    right_top = row_1[1]

    left_bottom = row_2[0]
    right_bottom = row_2[1]




    left_top.spines['right'].set_visible(False)
    right_top.spines['left'].set_visible(False)
    left_bottom.spines['right'].set_visible(False)
    right_bottom.spines['left'].set_visible(False)
    left_top.yaxis.tick_left()
    left_bottom.yaxis.tick_left()
    left_top.tick_params(labeltop='off')  # don't put tick labels at the top
    left_bottom.tick_params(labeltop='off')
    right_top.yaxis.tick_right()
    right_bottom.yaxis.tick_right()

    left_top.set_yticks([])
    left_top.set_yticklabels([])
    left_bottom.set_yticks([])
    left_bottom.set_yticklabels([])

    right_top.set_xticklabels([])
    left_top.set_xticklabels([])


    d = .015  # how big to make the diagonal lines in axes coordinates
    d2 = 0.02

    # arguments to pass plot, just so we don't keep repeating them
    kwargs = dict(transform=left_top.transAxes, color='k', clip_on=False)
    left_top.plot((1 - d, 1 + d), (-d, +d), **kwargs)
    left_top.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)

    kwargs.update(transform=right_top.transAxes)  # switch to the bottom axes
    right_top.plot((-d, +d), (1 - d, 1 + d), **kwargs)
    right_top.plot((-d, +d), (-d, +d), **kwargs)

    # arguments to pass plot, just so we don't keep repeating them
    kwargs = dict(transform=left_bottom.transAxes, color='k', clip_on=False)
    left_bottom.plot((1 - d, 1 + d), (-d, +d), **kwargs)
    left_bottom.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)

    kwargs.update(transform=right_bottom.transAxes)  # switch to the bottom axes
    right_bottom.plot((-d, +d), (1 - d, 1 + d), **kwargs)
    right_bottom.plot((-d, +d), (-d, +d), **kwargs)


    offset = 0
    legend_counter = 0

    for spectrum in speclist1:
        if labels1 == None:
            label = spectrum.name.full_name
        else:
            label = labels1[legend_counter]
            legend_counter += 1

        left_top.plot(spectrum.wavenumbers, spectrum.normalize_to_highest() + offset, marker="s",
                 markersize=3, label=label)

        right_top.plot(spectrum.wavenumbers, spectrum.normalize_to_highest() + offset, marker="s",
                      markersize=3, label=label)

        offset += 0.4

    offset = 0
    legend_counter = 0
#
    for spectrum in speclist2:

        if labels1 == None:
            label = spectrum.name.full_name
        else:
            label = labels1[legend_counter]
            legend_counter += 1

        left_bottom.plot(spectrum.wavenumbers, spectrum.normalize_to_highest() + offset, marker="s",
                 markersize=3, label=label)

        right_bottom.plot(spectrum.wavenumbers, spectrum.normalize_to_highest() + offset, marker="s",
                      markersize=3, label=label)

        offset += 0.4

    fig.text(0.025, 0.5, 'Norm. intensity/ arb. u.', ha='center', va='center', rotation='vertical', fontsize=12)
    fig.text(0.5, 0.015, 'Wavenumber/ cm$^{-1}$', ha='center', va='center', fontsize=12)
    bl = left_bottom.legend(fontsize=8)
    tl = left_top.legend(fontsize=8)

    bl.draggable()
    tl.draggable()

    left_top.set_xlim(lower[0], lower[1])
    right_top.set_xlim(upper[0], upper[1])

    left_bottom.set_xlim(lower[0], lower[1])
    right_bottom.set_xlim(upper[0], upper[1])
    print(upper, lower)

    plt.show()


def finalize_figure(fig, title="test2"):
    """Makes figures publication-ready and exports them as pdf. Takes the title for the output file as argument"""
    #size = fig.get_size_inches()
    #ratio = size[0] / size[1]
    #fig.set_size_inches(3.2 * ratio, 3.2)
    #fig.tight_layout()
    fig.savefig(title + ".pdf", dpi=600)

def get_by_name(names, scm):

    out = []
    for name in names:
        scm.get("name "+name)
        out.append(scm.subset[0])
        scm.clear()

    return out


def decay_double(speclist):

    fig, axes = plt.subplots(nrows=1, ncols=2, sharey=True)

    sfg_to_axes(axes[0], speclist)

    for i, spectrum in enumerate(speclist):
        axes[1].scatter(i, np.max(spectrum.normalized_intensity), color="r", edgecolor="black", alpha=0.7, s=22)

    axes[0].legend()
    axes[1].legend()
    axes[1].set_xlabel("time of photolysis/ min.")

    return fig


def sfg_to_axes(axes, speclist):

    for spectrum in speclist: # type: Classes.SfgSpectrum

        axes.plot(spectrum.wavenumbers, spectrum.normalized_intensity, marker="s")

    axes.set_xlabel("Wavenumber/ cm$^{-1}$")
    axes.set_ylabel("Norm. SFG intensity/ arb. u.")


S = scm("sfg.db", "test")
li = ["20180323_PA_BX12_10_x2_#1_0.5zu1mM", "20180323_PA_BX12_10_x2_#2_80p_0.5zu1mM", "20180323_PA_10_x2_#1_45p_5mM", "20180323_PA_10_x1_#2_5mM"]
#li = ["20170901_BX9_5_x2_#1_5mM", "20170901_BX9_5_x2_#2_30p_5mM"]
# li = ["20180130_BX12_3.5_x1_#6_4p_Peak1", "20180130_BX12_3.5_x1_#9_7p_Peak1",
#       "20180130_BX12_3.5_x1_#5_3p_Peak1", "20180130_BX12_3.5_x1_#8_6p_Peak1", "20180130_BX12_3.5_x1_#4_2p_Peak1",
#       "20180130_BX12_3.5_x1_#7_5p_Peak1", "20180130_BX12_3.5_x1_#3_1p_Peak1", "20180130_BX12_3.5_x1_#2_Peak1init"]








# S.get("name 20180323_PA_10_x1_#2_5mM")
#
# pa = S.subset[0]
# S.clear()
#
# S.get("name 20180315_SA_BX12_10_x1_#2_mixedlayer1zu2mM")
# sa_bx = S.subset[0]
# S.clear()
#
# S.get("name 20180320_PA_BX12_10_x1_#1_1.5zu1.5mM")
# pa_bx = S.subset[0]
# S.clear()
#
# S.get("name 20180524_BX12_8_x2_#2_1mM")
# bx = S.subset[0]
# S.clear()
#
# S.get("name 20160316_SA_5_x1_#1_5mM")
# sa = S.subset[0]
# S.clear()
#
# stack1 = [sa, sa_bx, bx]
# stack2 = [pa, pa_bx, bx]
# labels1 = ["Stearic Acid, pure", "Stearic Acid/ BX12 (1:2)", "BX12, pure"]
# labels2 = ["Palmitic Acid, pure", "Palmitic Acid/ BX12 (1:1)", "BX12, pure"]
# set_plot_properties()
# f = sfg_doublestack_plot_broken(stack1, stack2, [1600,1850], [2700,3300], labels1, labels2)
set_plot_properties()
li = get_by_name(li, S)
sfg_pub_plot(li)
plt.show()
# test