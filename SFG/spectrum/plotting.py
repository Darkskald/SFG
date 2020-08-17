import matplotlib.pyplot as plt
import numpy as np

# todo: remove all plotting from the core spectrum module

class BatchPlotter:

    def __init__(self, specdic, save=False, savedir="", savename="default", title="default"):
        self.specdic = specdic
        self.save = save
        self.savedir = savedir
        self.savename = savename
        self.title = title

    def plot_all(self, xlim=None, ylim=None, legend=True, line=False):
        for key in self.specdic:
                spectrum = self.specdic[key]
                if not line:
                    plt.plot(spectrum.x, spectrum.y, label=key, marker="^", linestyle="-")
                else:
                    plt.plot(spectrum.x, spectrum.y, label=key)

        plt.xlabel(spectrum.x_unit)
        plt.ylabel(spectrum.y_unit)
        plt.minorticks_on()
        plt.title(self.title)
        if xlim is not None:
            plt.xlim(*xlim)
        if ylim is not None:
            plt.ylim(*ylim)
        if legend:
            plt.legend()

        if self.save is False:
            plt.show()

        else:
            path = self.savedir + "/" + self.savename + ".png"
            plt.tight_layout()
            plt.savefig(path)
            plt.close()

    def plot_broken_axis(self, lim, line=False):
        fig, (ax, ax2) = plt.subplots(1, 2, sharey=True)

        ax.set_ylabel("intensity/ arb.u.")
        fig.text(0.5, 0.03, s="wavenumber/ cm$^{⁻1}$", ha="center", va="center")

        ax.set_xlim(lim[0], lim[1])
        ax2.set_xlim(lim[2], lim[3])

        ax.spines['right'].set_visible(False)
        ax2.spines['left'].set_visible(False)
        ax.yaxis.tick_left()

        ax2.yaxis.tick_right()

        d = .012  # how big to make the diagonal lines in axes coordinates
        # arguments to pass plot, just so we don't keep repeating them
        kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
        ax.plot((1 - d, 1 + d), (-d, +d), **kwargs)
        ax.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)

        kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
        ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)
        ax2.plot((-d, +d), (-d, +d), **kwargs)

        for key in self.specdic:
            spectrum = self.specdic[key]
            if not line:
                ax.plot(spectrum.x, spectrum.y, label=key, marker="^", linestyle="-")
                ax2.plot(spectrum.x, spectrum.y, label=key, marker="^", linestyle="-")
            else:
                ax.plot(spectrum.x, spectrum.y, label=key)
                ax2.plot(spectrum.x, spectrum.y, label=key)
        #l = fig.legend(*ax2.get_legend_handles_labels(), 'best')
        l = ax2.legend()
        l.set_draggable(True)
        if self.save is False:
            plt.show()
        else:
            path = self.savedir + "/" + self.savename + ".png"
            plt.savefig(path)
            plt.close()


def baseline_demo_dppc(spectrum, integral= "", coverage= ""):
    spectrum.correct_baseline()

    test = np.linspace(2750, 3050, 10000)
    func = spectrum.make_ch_baseline()


    borders = spectrum.slice_by_borders(np.min(spectrum.wavenumbers), 3000)

    f, axarr = plt.subplots(2, sharex=True)
    axarr[0].ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
    axarr[1].ticklabel_format(style='sci', axis='y', scilimits=(0, 0))

    axarr[0].plot(spectrum.wavenumbers, spectrum.normalized_intensity, label=spectrum.meta["name"],
                  marker="v", markersize=5, color="blue")
    axarr[0].plot(test, func(test), color="r", label="red")


    axarr[1].plot(spectrum.wavenumbers, spectrum.baseline_corrected, label=spectrum.meta["name"],
                  marker="v", markersize=5, color="blue")
    axarr[1].fill_between(spectrum.wavenumbers[borders[0]:borders[1] + 1],
                          spectrum.baseline_corrected[borders[0]:borders[1] + 1], facecolor="blue")
    axarr[1].set_xlabel("wavenumber/ cm$^{-1}$")

    if type(integral) != str:
        integral = "{0:.4e}".format(integral)
        axarr[1].plot([], [], label=f'integral: {integral}, coverage: {coverage}')

    #f.text(0.07, 0.5, 'norm. intensity/ arb. u.', ha='center', va='center', rotation='vertical')
    f.text(0.03, 0.5, 'norm. intensity/ arb. u.', ha='center', va='center', rotation='vertical')

    for ax in axarr.flat:
        ax.set_xlim(2750, 3800)

    #plt.tight_layout()
    plt.savefig(spectrum.meta["name"] + ".png")
    plt.close()


def advanced_baseline_demo_dppc(spectrum, integral="", coverage= ""):
    spectrum.correct_baseline()

    test = np.linspace(2750, 3050, 10000)
    func = spectrum.make_ch_baseline()

    borders = spectrum.slice_by_borders(2800, 3000)
    dangling_borders = spectrum.slice_by_borders(3670, 3760)
    oh2_borders = spectrum.slice_by_borders(3350, 3670)
    oh_borders = spectrum.slice_by_borders(3000, 3350)

    f, axarr = plt.subplots(2, sharex=True)
    axarr[0].ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
    axarr[1].ticklabel_format(style='sci', axis='y', scilimits=(0, 0))

    axarr[0].plot(spectrum.wavenumbers, spectrum.normalized_intensity, label=spectrum.meta["name"], linewidth=1.5,
                  marker="o", markersize=3)
    axarr[0].plot(test, func(test), color="r", label="baseline")

    axarr[1].plot(spectrum.wavenumbers, spectrum.baseline_corrected, label=spectrum.meta["name"], linewidth=1.5,
                  marker="o", markersize=3, color="black")

    axarr[1].fill_between(spectrum.wavenumbers[borders[0]:borders[1] + 1],
                          spectrum.baseline_corrected[borders[0]:borders[1] + 1], color="r")

    axarr[1].fill_between(spectrum.wavenumbers[dangling_borders[0]:dangling_borders[1] + 1],
                          spectrum.baseline_corrected[dangling_borders[0]:dangling_borders[1] + 1], color="blue")

    axarr[1].fill_between(spectrum.wavenumbers[oh_borders[0]:oh_borders[1] + 1],
                          spectrum.baseline_corrected[oh_borders[0]:oh_borders[1] + 1], color="g")

    axarr[1].fill_between(spectrum.wavenumbers[oh2_borders[0]:oh2_borders[1] + 1],
                          spectrum.baseline_corrected[oh2_borders[0]:oh2_borders[1] + 1], color="purple")


    axarr[1].set_xlabel("wavenumber/ cm$^{-1}$")

    axarr[1].axhline(0, linestyle="--", color="r")

    # OH  dangling
    axarr[1].axvline(3670, color="b")
    axarr[1].axvline(3760, color="b")

    # OH
    axarr[1].axvline(3005, color="g")
    axarr[1].axvline(3350, color="g")


    # OH2
    axarr[1].axvline(3670, color="purple")
    axarr[1].axvline(3350, color="purple")

    if type(integral) != str:
        integral = "{0:.4e}".format(integral)
        axarr[1].plot([], [], label=f'integral: {integral}, coverage: {coverage}')

    f.text(0.025, 0.5, 'norm. intensity/ arb. u.', ha='center', va='center', rotation='vertical')

    for ax in axarr.flat:
        ax.set_xlim(2750, 3800)

    plt.savefig(spectrum.meta["name"] + ".png")
    plt.close()