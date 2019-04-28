from spectrum import SfgSpectrum

import re
import sqlite3
from datetime import datetime, date

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.dates import MonthLocator, DateFormatter
from matplotlib.lines import Line2D


class SampleNameParser:
    # sample numbers:
    n1 = re.compile("\D\d{1,2}\D")
    n2 = re.compile("[ a-zA-Z_-]\d{1,2}$")
    n3 = re.compile("\d{1,2}-?#$")
    n4 = re.compile("-#\d{1,2}$")

    numreg = [n1, n2, n3, n4]
    # sampling dates:
    d1 = re.compile("^\d{8}_\d{1,2}\D")
    d2 = re.compile("_[a-zA-z -]*\d{8}")
    d3 = re.compile("^\d{8}_[a-zA-Z]*[ -]")
    d4 = re.compile("^\d{8}_\d{1,2}$")
    d5 = re.compile("^\d{8}_[a-zA-Z]{2}\d{1,2}-")

    datereg = [d1, d2, d3, d4, d5]
    sep = "*"*90+"\n"

    def __init__(self):
        self.log = ""
        self.sample_counter = 0
        self.df = pd.read_excel("Wasserproben_komplett.xlsx", header=2, sheet_name="Samples")
        self.df = self.df[self.df["Location No."] == 3]
        self.gasex_deep = {date(2018, 6, 1): (0.019299999999999998, 0.029944320879612154),
                           date(2018, 9, 1): (0.02368125, 0.03131569532099678)}
        self.gasex_sml = {date(2018, 6, 1): (0.02919075757575757, 0.02888070463447273),
                          date(2018, 9, 1): (0.04120840909090909, 0.024149130954924634)}

        self.deep = self.df[(self.df["Sampler no."] == 4) | (self.df["Sampler no."] == 3)]
        self.df = self.df[(self.df["Sampler no."] == 1) | (self.df["Sampler no."] == 2)]

        self.df = self.df.reset_index(drop=True)
        self.deep = self.deep.reset_index(drop=True)

        self.conn = sqlite3.connect("test.db")
        self.samples = pd.read_sql("SELECT * from sfg where type = \"boknis\"", self.conn)
        self.sfg = []
        self.reference_per_date = (self.get_dppc_intensities())
        self.new_df_cols()
        self.process()
        self.samples["sampling_date"] = pd.to_datetime(self.samples["sampling_date"])

        self.coverages = self.match_test(self.df)
        self.deep_coverages = self.match_test(self.deep)

        self.coverages.update(self.gasex_sml)
        self.deep_coverages.update(self.gasex_deep)

        self.plot_coverage()
        #self.write_log()

    def new_df_cols(self):

        cols = ["sampling_date", "number"]
        self.samples = self.samples.reindex(columns=self.samples.columns.tolist() + cols)

    def process(self):

        for index in range(len(self.samples)):

            number = None
            sampling_date = None
            name = self.samples.loc[index, "name"]
            number_match = False
            date_match = False

            for ex in SampleNameParser.numreg:

                temp = re.search(ex, name)

                if temp is not None:
                    refined = name[temp.start():temp.end()]
                    temp = re.search("\d{1,2}", refined)
                    number = refined[temp.start():temp.end()]
                    number_match = True
                    break

            if number_match is False:
                repstr = f'Number parsing was not possible for {name}\n'
                self.log += repstr

            for ex in SampleNameParser.datereg:

                temp = re.search(ex, name)

                if temp is not None:
                    refined = name[temp.start():temp.end()]
                    temp = re.search("\d{8}", refined)
                    sampling_date = refined[temp.start():temp.end()]
                    date_match = True
                    break

            if date_match is False:
                repstr = f'Date parsing was not possible for {name}\n'
                self.log += repstr

            if number is not None:
                number = int(number)

            self.samples.loc[index, "number"] = number
            self.samples.loc[index, "sampling_date"] = sampling_date

    def match_test(self, target):
        out = {}

        for i in range(len(target)):

            try:
                mask1 = (self.samples["sampling_date"] == target.loc[i, "Date"])
                mask2 = (self.samples["number"] == target.loc[i, "Sample"])
                sampling_date = target.loc[i, "Date"]
                match = self.samples[mask1 & mask2].reset_index(drop=True)
                self.sample_counter += 1

                for j in range(len(match)):
                    name = match.loc[j, "name"]
                    time = datetime.strptime((match.loc[j, "measured_time"]), '%Y-%m-%d %H:%M:%S')
                    meta = {"name": name, "time": time}
                    data = ("wavenumbers", "sfg", "ir", "vis")
                    tup = map(lambda x: np.fromstring(match.loc[j, x], sep=";"), data)
                    s = SfgSpectrum(*tup, meta)


                    info = f"""sample number {target.loc[i, "Sample"]} from sampling date {target.loc[i, "Date"]}
                    matched to {match.head()} \n"""
                    self.log += SampleNameParser.sep
                    self.log += f'Sample counter: {self.sample_counter}\n'
                    self.log += info

                    # now get the coverage
                    try:
                        measure_date = s.meta["time"].date()
                        factor = self.reference_per_date[measure_date]
                        integral = s.calculate_ch_integral(average="gernot")

                        if integral < 0:
                            integral = 0
                        coverage = round(np.sqrt(integral / factor), 4)

                        if sampling_date not in out:
                            out[sampling_date] = [coverage]
                        else:
                            out[sampling_date].append(coverage)
                        #baseline_demo_dppc(s, integral, coverage)

                    except IndexError:
                        errstr = f'Error occurred in inner block of match_test with {s.meta["name"]}\n'
                        self.log += errstr

                    except KeyError:
                        errstr = f'No DPPC reference found for spectrum {s.meta["name"]}\n'
                        self.log += errstr

            except TypeError:
                errstr = f'Type Error occurred in outer block of match_test for {target.iloc[i]}\n'
                self.log += errstr

        for i in out:
            av = np.average(np.array(out[i]))
            std = np.std(np.array(out[i]))
            out[i] = av, std

        return out

    def get_dppc_intensities(self):
        """Returns a dictionary with each day of measurement and the corresponding DPPC integrals"""
        cmd = "SELECT * FROM sfg WHERE type = \"boknis_ref\""
        dppc_specs = pd.read_sql(cmd, self.conn)
        dates = {}
        for i in range(len(dppc_specs)):
            name = dppc_specs.loc[i, "name"]

            time = datetime.strptime((dppc_specs.loc[i, "measured_time"]), '%Y-%m-%d %H:%M:%S')
            meta = {"name": name, "time": time}
            data = ("wavenumbers", "sfg", "ir", "vis")
            tup = map(lambda x: np.fromstring(dppc_specs.loc[i, x], sep=";"), data)
            s = SfgSpectrum(*tup, meta)
            #baseline_demo_dppc(s)
            sr = s.yield_spectral_range()

            # if sr[0] < 2800 and sr[1] > 3010 and name.split("_")[-1] != "ppp":
            if name.split("_")[-1] != "ppp":
                q = s.calculate_ch_integral(average="gernot")
                if time.date() not in dates:
                    dates[time.date()] = [q]
                else:
                    dates[time.date()].append(q)


        for item in dates:
            dates[item] = np.average(np.array(dates[item]))

        return dates

    def plot_coverage(self):
        fig, ax = plt.subplots()
        months = MonthLocator(range(1, 13), bymonthday=1, interval=3)
        monthsFmt = DateFormatter("%b '%y")

        for item in self.coverages:

            ax.errorbar(item, self.coverages[item][0], yerr=self.coverages[item][1],
                        markerfacecolor="red", marker="o", ecolor="red", capsize=5,
                        capthick=10, mec="black", mew=0.3, aa=True, elinewidth=1)

        for item in self.deep_coverages:

            ax.errorbar(item, self.deep_coverages[item][0], yerr=self.deep_coverages[item][1],
                        markerfacecolor="blue",  marker="o", ecolor="blue", capsize=5,
                        capthick=10, mec="black", mew=0.3, aa=True, elinewidth=1)

        for i in range(8, 20, 1):
            lower = date(2000+i, 3, 1)
            upper = date(2000 + i, 9, 1)
            ax.axvspan(lower, upper, color="g", alpha=0.4)

        ax.xaxis.set_major_locator(months)
        ax.xaxis.set_major_formatter(monthsFmt)
        ax.autoscale_view()
        fig.autofmt_xdate()
        legend_elements = [Line2D([0], [0], marker='o', label='Bulk water',
                          markerfacecolor='blue', mew=0.3,  mec="black", aa=True, linestyle=''),
                           Line2D([0], [0], marker='o', label='Surface microlayer',
                                   markerfacecolor='red', mew=0.3, mec="black", aa=True, linestyle='')
                            ]
        ax.legend(handles=legend_elements)
        ax.grid(True)
        ax.set_xlabel("time ")
        ax.set_ylabel("Surface coverage/ %")

        plt.title("Boknis Eck time series surface coverage measured by SFG",
                  fontdict={'fontweight': "bold", 'fontsize': 18})

        plt.show()

    def write_log(self):

        with open("log.txt", "w") as outfile:
            outfile.write(self.log)


def baseline_demo_dppc(spectrum, integral= "", coverage= ""):
    spectrum.correct_baseline(average="gernot")

    test = np.linspace(2750, 3050, 10000)
    func = spectrum.make_ch_baseline(average="gernot")
    borders = spectrum.slice_by_borders(3000, np.min(spectrum.wavenumbers))

    f, axarr = plt.subplots(2, sharex=True)
    axarr[0].ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
    axarr[1].ticklabel_format(style='sci', axis='y', scilimits=(0, 0))

    axarr[0].plot(spectrum.wavenumbers, spectrum.normalized_intensity, label=spectrum.meta["name"], linewidth=1.5,
                  marker="o", markersize=3)
    axarr[0].plot(test, func(test), color="r", label="baseline")

    axarr[1].plot(spectrum.wavenumbers, spectrum.baseline_corrected, label=spectrum.meta["name"], linewidth=1.5,
                  marker="o", markersize=3)
    axarr[1].fill_between(spectrum.wavenumbers[borders[0]:borders[1] + 1],
                          spectrum.baseline_corrected[borders[0]:borders[1] + 1])
    axarr[1].set_xlabel("wavenumber/ cm$^{-1}$")
    if type(integral) != str:
        integral = "{0:.4e}".format(integral)
        axarr[1].plot([], [], label=f'integral: {integral}, coverage: {coverage}')

    f.text(0.025, 0.5, 'norm. intensity/ arb. u.', ha='center', va='center', rotation='vertical')

    for ax in axarr.flat:
        ax.set_xlim(2750, 3300)

    plt.savefig("plots/"+spectrum.meta["name"] + ".png")
    plt.close()

plt.style.use(['seaborn-ticks', 'seaborn-notebook'])

rcParams['figure.figsize'] = 16.4, 10.8
rcParams['axes.labelsize'] = 18
rcParams['font.size'] = 18
rcParams['figure.subplot.bottom'] = 0.12


s = SampleNameParser()

# todo: remove odd DPPC references
# todo: search and remove odd-looking dppc samples
# todo: include the samples measured by Ann-Kathrin by renaming
# todo: remove outliers (with more than 100 % coverage)
# todo: document the module
# todo: write a logfile with all samples included in this analysis, check it
# todo: make a list with the most abundant things I do with the mouse to enhance productivity
# todo: search and remove odd-looking water samples
