from spectrum import SfgSpectrum

import pandas as pd
import sqlite3
import numpy as np

import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.dates import MonthLocator, WeekdayLocator, DateFormatter

from datetime import datetime
import re
plt.style.use(['seaborn-ticks', 'seaborn-notebook'])
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

    def __init__(self):
        self.df = pd.read_excel("Wasserproben_komplett.xlsx", header=2, sheet_name="Samples")
        self.df = self.df[self.df["Location No."] == 3]

        self.deep = self.df[(self.df["Sampler no."] == 4)]
        self.df = self.df[(self.df["Sampler no."] == 1) | (self.df["Sampler no."] == 2)]

        self.df = self.df.reset_index(drop=True)
        self.deep = self.deep.reset_index(drop=True)

        self.conn = sqlite3.connect("test.db")
        self.samples = pd.read_sql("SELECT * from sfg where type = \"boknis\"", self.conn)
        self.sfg = []
        self.reference_per_date = (self.get_dppc_intensities())
        self.coverages = {}

        self.new_df_cols()
        self.process()
        self.samples["sampling_date"] = pd.to_datetime(self.samples["sampling_date"])

        self.match_test()
        self.plot_coverage()


    def new_df_cols(self):

        cols = ["sampling_date", "number"]
        self.samples = self.samples.reindex(columns=self.samples.columns.tolist() + cols)

    def process(self):

        for index in range(len(self.samples)):

            number = None
            sampling_date = None
            name = self.samples.loc[index, "name"]

            for ex in SampleNameParser.numreg:

                temp = re.search(ex, name)

                if temp is not None:
                    refined = name[temp.start():temp.end()]
                    temp = re.search("\d{1,2}", refined)
                    number = refined[temp.start():temp.end()]
                    break

            for ex in SampleNameParser.datereg:

                temp = re.search(ex, name)

                if temp is not None:
                    refined = name[temp.start():temp.end()]
                    temp = re.search("\d{8}", refined)
                    sampling_date = refined[temp.start():temp.end()]
                    break

            if number is not None:
                number = int(number)

            self.samples.loc[index, "number"] = number
            self.samples.loc[index, "sampling_date"] = sampling_date

    def match_test(self):

        for i in range(len(self.df)):

            try:
                mask1 = (self.samples["sampling_date"] == self.df.loc[i, "Date"])
                mask2 = (self.samples["number"] == self.df.loc[i, "Sample"])
                sampling_date = self.df.loc[i, "Date"]
                match = self.samples[mask1 & mask2].reset_index(drop=True)

                for j in range(len(match)):
                    name = match.loc[j, "name"]
                    time = datetime.strptime((match.loc[j, "measured_time"]), '%Y-%m-%d %H:%M:%S')
                    meta = {"name": name, "time": time}
                    data = ("wavenumbers", "sfg", "ir", "vis")
                    tup = map(lambda x: np.fromstring(match.loc[j, x], sep=";"), data)
                    s = SfgSpectrum(*tup, meta)

                    # now get the coverage
                    try:
                        measure_date = s.meta["time"].date()
                        factor = self.reference_per_date[measure_date]
                        integral = s.calculate_ch_integral(average="gernot")

                        if integral < 0:
                            integral = 0
                        coverage = round(np.sqrt(integral / factor), 4)

                        if sampling_date not in self.coverages:
                            self.coverages[sampling_date] = [coverage]
                        else:
                            self.coverages[sampling_date].append(coverage)
                        #baseline_demo_dppc(s, integral, coverage)

                    except IndexError:
                        pass

                    except KeyError:
                        print(measure_date, "No DPPC referece")

            except TypeError:
                pass

        for i in self.coverages:
            self.coverages[i] = np.average(np.array(self.coverages[i]))

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
            ax.plot_date(item, self.coverages[item], xdate=True, color="red")

        ax.xaxis.set_major_locator(months)
        ax.xaxis.set_major_formatter(monthsFmt)
        ax.autoscale_view()
        fig.autofmt_xdate()
        ax.grid(True)
        plt.show()


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



rcParams['axes.labelsize'] = 14
rcParams['font.size'] = 14
rcParams['figure.subplot.bottom'] = 0.12

s = SampleNameParser()

