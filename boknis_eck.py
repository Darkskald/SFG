from orm import WorkDatabaseWizard
from spectrum import SfgAverager, DummyPlotter

from sqlalchemy import extract
from sqlalchemy.exc import IntegrityError
from sqlalchemy.orm.exc import NoResultFound
# todo: fix the double import
from datetime import datetime, date, timedelta
import datetime
import functools
import pandas as pd
import re
import numpy as np

import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.dates import MonthLocator, DateFormatter
from matplotlib.lines import Line2D
import matplotlib.ticker as ticker


class BoknisEckExtension:
    # sample number regular expressions to extract the number of a sample from the filename:
    n1 = re.compile("\D\d{1,2}\D")
    n2 = re.compile("[ a-zA-Z_-]\d{1,2}$")
    n3 = re.compile("\d{1,2}-?#$")
    n4 = re.compile("-#\d{1,2}$")

    numreg = [n1, n2, n3, n4]
    # sampling dates to extract the sampling date of a sample from the filename:
    d1 = re.compile("^\d{8}_\d{1,2}\D")
    d2 = re.compile("_[a-zA-z -]*\d{8}")
    d3 = re.compile("^\d{8}_[a-zA-Z]*[ -]")
    d4 = re.compile("^\d{8}_\d{1,2}$")
    d5 = re.compile("^\d{8}_[a-zA-Z]{2}\d{1,2}-")

    datereg = [d1, d2, d3, d4, d5]
    sep = "*" * 90 + "\n"

    def __init__(self, new=False):
        self.wz = WorkDatabaseWizard()
        # read the master excel sheet
        self.df = pd.read_excel("Wasserproben_komplett.xlsx", header=2, sheet_name="Samples")

        """
        don't forget the GasEx data!
        """
        if new:
            self.write_references()
            self.persist_gasex_references()
            self.references = self.fetch_dppc_integrals()
            self.process()
            self.match_to_table()
            self.chlorophyll = BoknisEckExtension.prepare_chorophyll_data()
            data = self.map_to_sampling_date(self.references)
            self.persist_average_data(data[0], data[1], data[2])
            self.include_gasex_average()

        self.references = self.fetch_dppc_integrals()

    # setting up the data
    def get_references(self):
        """Retrieve the DPPC integral per day of measurement"""

        results = self.wz.session.query(self.wz.sfg).filter(self.wz.sfg.type == "boknis_ref").all()
        dates = {}

        for spec in results:
            spec = self.wz.construct_sfg(spec)

            if 0 <= spec.meta["time"].hour < 8:
                spec.meta["time"] -= timedelta(days=1)

            date = spec.meta["time"].date()

            if date not in dates:
                dates[date] = [spec]
            else:
                dates[date].append(spec)

        for date in dates:
            av = SfgAverager(dates[date])
            integral = av.integral
            dates[date] = integral

        return dates

    def write_references(self):
        """Write the reference integrals to the database"""

        refs = self.get_references()
        temp = []

        for date in refs:
            measurement_day = self.wz.measurement_days()
            measurement_day.dppc_integral = refs[date]
            measurement_day.date = date
            temp.append(measurement_day)

        self.wz.session.add_all(temp)
        self.wz.session.commit()

    def fetch_dppc_integrals(self):
        """Get the DPPC integrals stored in the database"""
        measurement_days = self.wz.session.query(self.wz.measurement_days).all()
        out = {i.date: i.dppc_integral for i in measurement_days}
        return out

    def fetch_gasex_references(self):
        """Get the DPPC references for the GasEx samples  from the database"""
        dates = {}
        q = self.wz.session.query(self.wz.sfg).filter(self.wz.sfg.measured_time.between("2018-06-01", "2018-12-31"))
        q = q.filter(self.wz.sfg.name.like("%DPPC%"))

        for item in q.all():
            date = item.measured_time.date()
            if date not in dates:
                dates[date] = [self.wz.construct_sfg(item)]
            else:
                dates[date].append(self.wz.construct_sfg(item))

        for key in dates:
            dates[key] = SfgAverager(dates[key]).integral

        return dates

    def persist_gasex_references(self):
        """Add the DPPC reference date of the GasEx period to the databse"""

        d = self.fetch_gasex_references()
        for key in d:
            try:
                measurement_day = self.wz.measurement_days()
                measurement_day.date = key
                measurement_day.dppc_integral = d[key]
                self.wz.session.add(measurement_day)
                self.wz.session.commit()
            except IntegrityError:
                self.wz.session.rollback()

    def process(self):
        """Populates the SQL table for BoknisEck spectra with metadata obtained
        from the systematic names of the SFG table."""

        samples = self.wz.session.query(self.wz.sfg).filter(self.wz.sfg.type == "boknis")
        boknis_specs = []

        for item in samples:
            boknis_spec = self.wz.boknis_eck()
            boknis_spec.specid = item.id
            boknis_spec.name = item.name

            for ex in BoknisEckExtension.numreg:
                temp = re.search(ex, item.name)

                if temp is not None:
                    refined = item.name[temp.start():temp.end()]
                    temp = re.search("\d{1,2}", refined)
                    number = refined[temp.start():temp.end()]
                    number_match = True
                    break

            if number_match is False:
                repstr = f'Number parsing was not possible for {item.name}\n'
                self.log += repstr

            for ex in BoknisEckExtension.datereg:

                temp = re.search(ex, item.name)

                if temp is not None:
                    refined = item.name[temp.start():temp.end()]
                    temp = re.search("\d{8}", refined)
                    sampling_date = refined[temp.start():temp.end()]
                    date_match = True
                    break

            if date_match is False:
                repstr = f'Date parsing was not possible for {item.name}\n'
                self.log += repstr

            if number is not None:
                number = int(number)

            boknis_spec.sample_number = number
            boknis_spec.sampling_date = BoknisEckExtension.convert_date(sampling_date)
            boknis_spec.is_mapped = False

            boknis_specs.append(boknis_spec)

        self.wz.session.add_all(boknis_specs)
        self.wz.session.commit()

    def match_to_table(self):
        """Traverses trough the lines of the master Excel sheet and maps the spectra of
        the Boknis Eck SQL table to them."""

        for index, row in self.df.iterrows():
            date = row["Date"].date()
            number = row["Sample"]

            try:
                q = self.wz.session.query(self.wz.boknis_eck).\
                    filter(self.wz.boknis_eck.sampling_date == date)\
                    .filter(self.wz.boknis_eck.sample_number == number).all()

                q = q[0]

                if row["Sampler no."] in (4, 3):
                    q.sample_type = "deep"
                    q.depth = int(row["Drainage time per dip (sec)\nSampling depth (m)"])

                elif row["Sampler no."] in (1, 2):
                    q.sample_type = "sml"
                    q.depth = 0

                q.is_mapped = True
                q.location_number = int(row["Location No."])
                self.wz.session.commit()

            except:
                pass

    # operating on the data
    def get_boknis_specs(self):
        """This function retrieves all Boknis Eck samples from the boknis_eck SQL table."""
        q = self.wz.session.query(self.wz.boknis_eck).\
            filter(self.wz.boknis_eck.is_mapped == 1).\
            filter(self.wz.boknis_eck.location_number == 3)
        return q

    def add_spectrum(self, query):
        """This function retrieves via Foreign-Key-relationship the spectral data of the SFG
        spectrum from the database and constructs the actual SFG object, setting it as attribute
        of the corresponding boknis_eck object."""

        out = []
        for item in query.all():
            temp = self.wz.session.query(self.wz.sfg).filter(self.wz.sfg.id == item.specid).one()
            item.sfg_spectrum = self.wz.construct_sfg(temp)
            out.append(item)

        return out

    def retrieve_reference(self):
        """This function creates a dictionary of DPPC integrals mapped to days of measurement
        from the measurement_days SQL table."""

        reference_integrals = {}
        query = self.wz.session.query(self.wz.measurement_days).all()
        for item in query:
            reference_integrals[item.date] = item.dppc_integral
        return reference_integrals

    def map_to_sampling_date(self, references):

        sml = self.add_spectrum(self.get_boknis_specs().filter(self.wz.boknis_eck.sample_type == 'sml'))
        deep = self.add_spectrum(self.get_boknis_specs().filter(self.wz.boknis_eck.sample_type == 'deep').
                                 filter(self.wz.boknis_eck.depth > 1))

        one = self.add_spectrum(self.get_boknis_specs().filter(self.wz.boknis_eck.sample_type == 'deep').
                                filter(self.wz.boknis_eck.depth == 1))


        out = []
        for item in (sml, deep, one):
            dates = {}
            for spectrum in item:
                date = spectrum.sampling_date
                if date not in dates:
                    dates[date] = [spectrum]
                else:
                    dates[date].append(spectrum)
            dates = self.average_sampling_dates(dates, references)
            out.append(dates)

        return out

    def average_sampling_dates(self, date_dict, references):
        new_dict = {}
        for key in date_dict:
            to_average = [i.sfg_spectrum for i in date_dict[key]]
            av = SfgAverager(to_average, references, enforce_scale=True)
            if av.coverage == np.inf:
                av.benchmark()
            new_dict[key] = av
        return new_dict

    def persist_average_data(self, sml_dates, bulk_dates, one_dates):
        """Write the data obtained by the average spectra to the database"""
        # todo: this is the right place for intense debugging and plotting

        for date in sml_dates:
            be_data_orm = self.wz.be_data()
            be_data_orm.sampling_date = date

            # sml
            be_data_orm.sml_no = len(sml_dates[date].spectra)
            be_data_orm.sml_coverage = sml_dates[date].coverage
            be_data_orm.sml_ch = sml_dates[date].average_spectrum.calc_region_integral("CH")
            be_data_orm.sml_oh1 = sml_dates[date].average_spectrum.calc_region_integral("OH")
            be_data_orm.sml_oh2 = sml_dates[date].average_spectrum.calc_region_integral("OH2")
            be_data_orm.sml_dangling = sml_dates[date].average_spectrum.calc_region_integral("dangling")

            be_data_orm.chlorophyll = BoknisEckExtension.get_mean_by_date(self.chlorophyll, date)

            # bulk
            try:
                be_data_orm.bulk_no = len(bulk_dates[date].spectra)
                be_data_orm.bulk_coverage = bulk_dates[date].coverage
                be_data_orm.bulk_ch = bulk_dates[date].average_spectrum.calc_region_integral("CH")
                be_data_orm.bulk_oh1 = bulk_dates[date].average_spectrum.calc_region_integral("OH")
                be_data_orm.bulk_oh2 = bulk_dates[date].average_spectrum.calc_region_integral("OH2")
                be_data_orm.bulk_dangling = bulk_dates[date].average_spectrum.calc_region_integral("dangling")
            except KeyError:
                pass

            # one meter
            try:
                be_data_orm.one_coverage = one_dates[date].coverage
            except KeyError:
                pass

            finally:
                self.wz.session.add(be_data_orm)
                self.wz.session.commit()

        # todo: bacpopulate spectra table

    def provide_dataframe(self):
        """A convenience function providing the compiled Boknis Eck data SQL table as Pandas dataframe"""
        query = self.wz.session.query(self.wz.be_data)
        df = pd.read_sql(query.statement, query.session.bind)
        return df

    # include gasex measurements
    def calculate_gasex_average(self):
        """Calculates the average SML and bulk spectra for the GasEx cruise, both for June and September"""

        sml = self.wz.session.query(self.wz.samples).filter(self.wz.samples.type.in_(("s", "p")))
        deep = self.wz.session.query(self.wz.samples).filter(self.wz.samples.type == "deep")

        sml_june = sml.filter(self.wz.samples.sample_hash.like('06%'))
        sml_september = sml.filter(self.wz.samples.sample_hash.like('09%'))

        deep_june = deep.filter(self.wz.samples.sample_hash.like('06%'))
        deep_september = deep.filter(self.wz.samples.sample_hash.like('09%'))

        input = [self.construct_gasex_sfgs(i) for i in (sml_june, sml_september, deep_june, deep_september)]
        func = functools.partial(BoknisEckExtension.get_average_spectrum, references=self.references)
        output = map(func, input)
        return list(output)

    def construct_gasex_sfgs(self, query):
        """Collects the spectra specified in query and converts them to actual SfgSpectrum objects for averaging"""

        samples = query.all()
        to_average = []
        for sample in samples:
            try:
                spec = self.wz.session.query(self.wz.gasex_sfg).filter(self.wz.gasex_sfg.sample_hash == sample.sample_hash).one()
                spec_raw = self.wz.session.query(self.wz.sfg).filter(self.wz.sfg.name == spec.name).one()
                spec = self.wz.construct_sfg(spec_raw)
                to_average.append(spec)
            except NoResultFound:
                pass
        return to_average

    def include_gasex_average(self):
        """Creates an artificial set of Boknis Eck sampling days to store the information calculated by the GasEx
        average spectra"""
        average_spectra = self.calculate_gasex_average()

        # setup new sampling dates for database
        june = self.wz.be_data()
        september = self.wz.be_data()
        june.sampling_date = datetime.date(2018, 6, 21)
        september.sampling_date = datetime.date(2018, 9, 20)

        # set the necessary values for june
        june.sml_no = 1
        june.sml_coverage = average_spectra[0].coverage
        june.sml_ch = average_spectra[0].integral
        june.sml_oh1 = average_spectra[0].calc_region_integral("OH")
        june.sml_oh2 = average_spectra[0].calc_region_integral("OH2")
        june.sml_dangling = average_spectra[0].calc_region_integral("dangling")
        june.bulk_no = 1
        june.bulk_coverage = average_spectra[2].coverage
        june.bulk_ch = average_spectra[2].integral
        june.bulk_oh1 = average_spectra[2].calc_region_integral("OH")
        june.bulk_oh2 = average_spectra[2].calc_region_integral("OH2")
        june.bulk_dangling = average_spectra[2].calc_region_integral("dangling")
        june.chlorophyll = BoknisEckExtension.get_mean_by_date(self.chlorophyll, june.sampling_date)
        
        # set the values for september
        september.sml_no = 1
        september.sml_coverage = average_spectra[1].coverage
        september.sml_ch = average_spectra[1].integral
        september.sml_oh1 = average_spectra[1].calc_region_integral("OH")
        september.sml_oh2 = average_spectra[1].calc_region_integral("OH2")
        september.sml_dangling = average_spectra[1].calc_region_integral("dangling")
        september.bulk_no = 1
        september.bulk_coverage = average_spectra[3].coverage
        september.bulk_ch = average_spectra[3].integral
        september.bulk_oh1 = average_spectra[3].calc_region_integral("OH")
        september.bulk_oh2 = average_spectra[3].calc_region_integral("OH2")
        september.bulk_dangling = average_spectra[3].calc_region_integral("dangling")
        september.chlorophyll = BoknisEckExtension.get_mean_by_date(self.chlorophyll, september.sampling_date)

        self.wz.session.add_all([june, september])
        self.wz.session.commit()

    @staticmethod
    def prepare_chorophyll_data():
        be = pd.read_csv("newport/be_data.csv", sep=",", header=0)
        be = be[be["Depth [m]"] == 1]
        be["Time"] = pd.to_datetime(be["Time"])
        be["Time"] = be["Time"].apply(lambda x: x.date())
        return be

    @staticmethod
    def get_mean_by_date(df, date):
        temp = df[df["Time"] == date]
        return temp["chlora"].mean()

    # auxiliary functions

    @staticmethod
    def convert_date(date):
        temp = str(date)
        year = int(temp[0:4])
        month = int(temp[4:6])
        day = int(temp[6:])

        return datetime.date(year, month, day)

    @staticmethod
    def get_average_spectrum(speclist, references):
        s = SfgAverager(speclist, references=references, enforce_scale=True)
        av = s.average_spectrum
        av.integral = s.integral
        av.coverage = s.coverage
        return av


class Plotter:

    def __init__(self):
        self.fig = plt.figure()

    def plot_sfg_list(self, speclist, title="default", save=False):
        ax = self.fig.add_subplot(1, 1, 1)
        ax.set_xlabel(speclist[0].x_unit)
        ax.set_ylabel(speclist[0].y_unit)

        for item in speclist:
            ax.plot(item.x, item.y, label=item.name)

        ax.set_title(title)

        self.fig.legend()
        if save:
            plt.style.use("output.mplstyle")
            plt.tight_layout()
            plt.savefig("boknis_dates/"+title+".png")
            plt.close()

    def plot_raw_ir_vis(self, spec, save=False):
        # task1: SFG plot with IR, Vis, Raw
        ax = self.fig.add_subplot(1, 1, 1)

        # set limit for x axis to 4 digits (for wavenumber display)
        formatter = ticker.ScalarFormatter()
        formatter.set_powerlimits((-3, 4))
        ax.xaxis.set_major_formatter(formatter)
        ax.set_xlabel(spec.x_unit)

        # SFG
        ax.plot(spec.x, spec.raw_intensity, color="blue", label="SFG", marker="s")
        ax.set_ylabel("raw SFG intensity/ counts")

        # IR
        ax2 = ax.twinx()
        ax2.plot(spec.x, spec.ir_intensity, color="red", label="IR")
        ax2.plot(spec.x, spec.vis_intensity, color="green", label="VIS")
        ax2.set_ylabel("intensity/ counts")

        ax.grid(True)
        ax.set_title(spec.name)
        if save:
            plt.style.use("output.mplstyle")
            plt.tight_layout()
            self.fig.legend()
            plt.savefig(f'boknis_spectra/{spec.name}_raw.png')
            plt.close()

    def plot_raw_ir_vis_norm(self, spec, baseline=True, integral=False, save=False):
        # task1: SFG plot with IR, Vis, Raw
        ax = self.fig.add_subplot(2, 1, 1)

        # set limit for x axis to 4 digits (for wavenumber display)
        formatter = ticker.ScalarFormatter()
        formatter.set_powerlimits((-3, 4))
        ax.xaxis.set_major_formatter(formatter)
        ax.set_xlabel(spec.x_unit)

        # SFG
        ax.plot(spec.x, spec.raw_intensity, color="blue", label="SFG", marker="s")
        ax.set_ylabel("raw SFG intensity/\n counts")

        # IR
        ax2 = ax.twinx()
        ax2.plot(spec.x, spec.ir_intensity, color="red", label="IR")
        ax2.plot(spec.x, spec.vis_intensity, color="green", label="VIS")
        ax2.set_ylabel("intensity/\n counts")

        ax.grid(True)
        ax.set_title(spec.name)

        # normalized
        ax3 = self.fig.add_subplot(2, 1, 2, sharex=ax)
        ax3.set_ylabel("norm. SFG intensity/\narb. u.")
        ax3.plot(spec.x, spec.y, label="norm.", color="black", marker="o")
        ax3.grid(True)

        # baseline
        if baseline:
            base_data = np.linspace(2750, 3050, 10000)
            func = spec.make_ch_baseline()
            ax3.plot(base_data, func(base_data), color="purple", label="baseline", linewidth=2.5, alpha=0.75)

        # integral
        if integral:
            ax3.text(3050, np.max(spec.normalized_intensity)/2, f'integral: {spec.calculate_ch_integral():.4f}')


        plt.style.use("output.mplstyle")

        plt.tight_layout()
        # remove vertical gap between subplots
        plt.subplots_adjust(hspace=.0)


        if save:

            self.fig.legend(loc=4)
            plt.savefig(f'boknis_spectra/{spec.name}_raw.png')
            plt.close()

    def plot_sfg_averager(self, averager):

        plt.style.use("output.mplstyle")
        ax = self.fig.add_subplot(1, 1, 1)

        for spectrum in averager.spectra:
            ax.plot(spectrum.x, spectrum.y, alpha=0.35)

        ax.plot(averager.average_spectrum.x, averager.average_spectrum.y, marker="s", color="r", linewidth=2.5,
                label="average")

        ax.set_xlabel(averager.average_spectrum.x_unit)
        ax.set_ylabel(averager.average_spectrum.y_unit)



        ax.legend()


class BEDatabaseWizard(WorkDatabaseWizard):

    def __init__(self):
        super().__init__()
        self.df = BoknisEckExtension().provide_dataframe()

        # convert to datetime
        self.df["sampling_date"] = pd.to_datetime(self.df["sampling_date"])

        # remove outliers
        self.df["sml_coverage"] = self.df["sml_coverage"].mask(self.df["sml_coverage"] > 1)
        self.df["bulk_coverage"] = self.df["bulk_coverage"].mask(self.df["bulk_coverage"] > 1)
        self.df["one_coverage"] = self.df["one_coverage"].mask(self.df["one_coverage"] > 1)
        self.df["sml_coverage"] = self.df["sml_coverage"].replace([np.inf, -np.inf], np.nan)
        self.df["bulk_coverage"] = self.df["bulk_coverage"].replace([np.inf, -np.inf], np.nan)

        # remove nans in coverage
        # self.df = self.df[self.df["sml_coverage"].notna()]
        # self.df = self.df[self.df["bulk_coverage"].notna()]

    def filter_date(self, beginning, end):
        """Get all BE samples taken between beginning and end and returns them as a dataframe"""
        beginning = pd.Timestamp(beginning)
        end = pd.Timestamp(end)
        mask = (self.df["sampling_date"] > beginning) & (self.df["sampling_date"] < end)
        return self.df[mask]

    def filter_year(self, year):
        """Returns a dataframe containing all samples taken in a specified year"""
        beginning = datetime.date(year, 1, 1)
        end = datetime.date(year, 12, 31)
        return self.filter_date(beginning, end)

    def filter_month(self, month):
        """Returns a dataframe containing all samples taken in a specified month"""
        mask = (self.df["sampling_date"].dt.month == month)
        return self.df[mask]

    def get_be_spectra_monthwise(self, sample_type="all"):
        """Returns a dict of spectra mapped to their month of sampling"""
        out = {}
        if sample_type == "all":
            temp = self.session.query(self.boknis_eck).filter(self.boknis_eck.is_mapped == 1).all()

        elif sample_type == "sml":
            temp = self.session.query(self.boknis_eck).filter(self.boknis_eck.sample_type == "sml").all()
        elif sample_type == "1":
            temp = self.session.query(self.boknis_eck).filter(self.boknis_eck.sample_type == "deep")\
                .filter(self.boknis_eck.depth == 1).all()
        elif sample_type == "bulk":
            temp = self.session.query(self.boknis_eck).filter(self.boknis_eck.sample_type == "deep") \
                .filter(self.boknis_eck.depth > 1).all()

        for item in temp:
            if item.sampling_date.month not in out:
                out[item.sampling_date.month] = [self.fetch_by_specid(item.specid)]
            else:
                out[item.sampling_date.month].append(self.fetch_by_specid(item.specid))

        return out

    def process_years(self, lower=2009, upper=2018):
        """Returns a dictionary of dataframes with the BE data and additional normalized to the year's maximum value
        columns"""
        out = {}
        years = (i for i in range(lower, upper+1) if i != 2016)
        for year in years:
            out[year] = BEDatabaseWizard.normalize_year_records(self.filter_year(year))
        return out

    # convenience query methods
    def get_data_per_sampling_date(self):
        """Returns a dictionary of BoknisEck spectra objects mapped to their date of sampling. Note that in that stadium
        they are not plottable SFG objects."""
        dates = self.session.query(self.be_data).all()
        dates = [i.sampling_date for i in dates]

        date_dic = {}

        for item in dates:
            temp = self.session.query(self.boknis_eck).filter(self.boknis_eck.sampling_date == item)\
                .filter(self.boknis_eck.is_mapped == 1).\
                filter(self.boknis_eck.location_number == 3).all()
            date_dic[item] = temp

        return date_dic

    def convert_be_to_sfg(self, be):
        """Accepts a BoknisEck spectrum object. Fetches the corresponding SFG data from the database and creates the
        actual object"""
        temp = self.session.query(self.sfg).filter(self.sfg.id == be.specid).one()
        return self.construct_sfg(temp)

    def fetch_by_month(self, month, refine="all"):
        """Fetches all spectra sampled on a special month (passed as integer). The refine kwarg makes it possible to
        filter sml and bulk samples"""
        t = self.session.query(self.boknis_eck).filter(self.boknis_eck.is_mapped == 1).\
            filter(self.boknis_eck.location_number == 3)\
            .filter(extract('month', self.boknis_eck.sampling_date) == month)

        if refine == "sml":
            t = t.filter(self.boknis_eck.sample_type == "sml")

        elif refine == "deep":
            t = t.filter(self.boknis_eck.sample_type == "deep").filter(self.boknis_eck.depth > 1)

        elif refine == 1:
            t = t.filter(self.boknis_eck.sample_type == "deep").filter(self.boknis_eck.depth == 1)

        print(f'debug: refine is {refine}')
        return t.all()

    def fetch_by_quartal(self, refine="all", selection="q"):
        """Map the BoknisEck spectra to their corresponding quartals of the year. The refine kwarg makes it possible to
        filter sml and bulk samples"""
        # todo: implement 1m
        # todo: write a generic funtion to refine any query by those criteria
        # todo: find appropriate word for "quartal" and refactor

        quartals = {"q1": (1, 2, 3), "q2": (4, 5, 6), "q3": (7, 8, 9), "q4": (10, 11, 12)}
        trimester = {"t1": (3, 4, 5, 6), "t2": (7, 8, 9, 10), "t3": (11, 12, 1, 2)}

        to_get = None
        out = {}

        if selection == "q":
            to_get = quartals
        elif selection == "t":
            to_get = trimester

        for q in to_get:

            temp = []
            for month in to_get[q]:
                specs = self.fetch_by_month(month, refine=refine)
                temp.extend(specs)
            out[q] = temp

        return out

    def normalize_to_reference_integral(self, speclist):
        """This function takes a list of BoknisEck spectra, converts them to SfgSpectrum objects
         and normalizes them to the DPPC integral of ther measurement day. Note that spectra
         without suitable reference are dropped and are not included in the output"""
        out = []
        for spectrum in speclist:
            spec = self.convert_be_to_sfg(spectrum)
            if 0 <= spec.meta["time"].hour < 8:
                spec.meta["time"] -= timedelta(days=1)
            try:
                ref = self.session.query(self.measurement_days).filter(
                self.measurement_days.date == spec.meta["time"].date()).one().dppc_integral
                spec.normalized_intensity = spec.normalize_to_highest(external_norm=ref)
                spec.setup_spec()
                out.append(spec)

            except NoResultFound:
                print(f'{spec} has no DPPC reference')
        return out

    # plotting

    def plot_by_sampling_date(self):
        """Fetch all spectra of all sampling days, match them to the sampling dates, plot them and save them as png"""
        dates = self.get_data_per_sampling_date()

        for d in dates:
            #dates[d] = [self.convert_be_to_sfg(item) for item in dates[d]]
            dates[d] = self.normalize_to_reference_integral(dates[d])
        for item in dates:
            try:
                P = Plotter()
                P.plot_sfg_list(dates[item], title=str(item), save=True)
            except IndexError:
                pass

    def plot_all_raw_be(self):
        """Creates plots of all BE data showing raw intensity, IR and vis"""
        temp = self.session.query(self.boknis_eck).\
            filter(self.boknis_eck.location_number == 3).filter(self.boknis_eck.is_mapped == 1).all()
        for item in temp:
            spec = self.convert_be_to_sfg(item)
            #Plotter().plot_raw_ir_vis(spec, save=True)
            Plotter().plot_raw_ir_vis_norm(spec, save=True)

    def plot_all_be_references(self):
        """Creates plots of all BE data showing raw intensity, IR and vis"""
        temp = self.session.query(self.sfg).filter(self.sfg.type == "boknis_ref").all()
        for item in temp:
            spec = self.construct_sfg(item)
            Plotter().plot_raw_ir_vis_norm(spec, save=True, integral=True)

    @staticmethod
    def normalize_year_records(df):
        """Normalize the dataset of a particular year to the corresponding maximum values"""
        # todo: implement this for all the other params
        for i in ("sml_coverage", "bulk_coverage", "sml_ch", "bulk_ch"):
            df[f'norm_{i}'] = df[i] / df[i].mean()
        return df

    @staticmethod
    def convert_to_origin_date(df):
        """Sets the dataframe'S sampling_date column to an origin-friendly setting"""
        df['sampling_date'] = df["sampling_date"].dt.strftime('%d.%m.%Y')

    @staticmethod
    def compare_to_gernot(df):

        df['origin'] = (df['sml_no'] * df['sml_coverage']/(df['bulk_no'] + df['sml_no'])) + \
        (df['bulk_no'] * df['bulk_coverage'] / (df['bulk_no'] + df['sml_no']))

        df['origin'] = (df['origin']**2)*50
        BEDatabaseWizard.convert_to_origin_date(df)
        new = df[['sampling_date', 'origin']]
        new.to_csv("origin_boknis_out.csv")


# auxiliary functions not contained in classes
def scale_to_polar(array):
    """Normalize the input array(day of the year) to the interval [0, 2 PI]. Useful for polar plots."""
    min_ = 1
    max_ = 365
    scale = np.pi * 2 * ((array - min_) / (max_ - min_))
    return scale


def plot_by_time(dataframes, param, leg, scale=1):
    """BE time series plotting with rectangular stripes to help with orientation in the plot"""
    fig, ax = plt.subplots()
    ax2 = ax.twinx()

    months = MonthLocator(range(1, 13), bymonthday=1, interval=6)
    monthsFmt = DateFormatter("%b '%y")

    for df in dataframes:
        ax.scatter(df["sampling_date"], df["sml_" + param] * scale, color="red")
        ax.scatter(df["sampling_date"], df["bulk_" + param] * scale, color="blue", marker="^")
        ax2.scatter(df["sampling_date"], df["chlorophyll"], color="green", marker="+")

        # ax.plot(df["sampling_date"], df["sml_"+param]*scale, color="red")
        # ax.plot(df["sampling_date"], df["bulk_"+param]*scale, color="blue")
        # ax2.plot(df["sampling_date"], df["chlorophyll"], color="green")

    # stripes to indicate certain times of the year
    for i in range(8, 20, 1):
        lower = str(date(2000 + i, 6, 1))
        upper = str(date(2000 + i, 12, 1))
        ax.axvspan(lower, upper, color="gray", alpha=0.4)

    # average lines
    ax.axhline(df["sml_" + param].mean(), color="red", linestyle="dashed")
    ax.axhline(df["bulk_" + param].mean(), color="blue", linestyle="dashed")

    legend_elements = [Line2D([0], [0], marker='^', label='Bulk water',
                              markerfacecolor='blue', mew=0.3, mec="blue", aa=True, linestyle=''),
                       Line2D([0], [0], marker='o', label='Surface microlayer',
                              markerfacecolor='red', mew=0.3, mec="red", aa=True, linestyle=''),
                       Line2D([0], [0], marker='+', label='Chlorophyll a',
                              markerfacecolor='green', mew=2, mec="green", aa=True, linestyle='', markersize=10)
                       ]

    ax.legend(handles=legend_elements)
    ax.xaxis.set_major_locator(months)
    ax.xaxis.set_major_formatter(monthsFmt)
    ax.xaxis.set_minor_locator(MonthLocator())
    ax.autoscale_view()
    fig.autofmt_xdate()

    ax.set_xlabel("time ")
    ax.set_ylabel(leg)
    ax2.set_ylabel("Chlorophyll a concentration/\n µg/L")
    rcParams['xtick.labelsize'] = 'small'


def plot_ratio_by_time(dataframes, params):
    """Plotting function to deal with the ratio of different region integrals"""
    fig, ax = plt.subplots()
    ax2 = ax.twinx()

    months = MonthLocator(range(1, 13), bymonthday=1, interval=6)
    monthsFmt = DateFormatter("%b '%y")

    for df in dataframes:

        if params[1] == "oh":
            sml_total_oh = df["sml_oh1"] + df["sml_oh2"]
            bulk_total_oh = df["bulk_oh1"] + df["bulk_oh2"]
            sml_ratio = df["sml_" + params[0]] / sml_total_oh
            bulk_ratio = df["bulk_" + params[0]] / bulk_total_oh

        else:
            sml_ratio = df["sml_" + params[0]] / df["sml_" + params[1]]
            bulk_ratio = df["bulk_" + params[0]] / df["bulk_" + params[1]]

        ax.scatter(df["sampling_date"], sml_ratio, color="red")
        ax.scatter(df["sampling_date"], bulk_ratio, color="blue", marker="^")
        ax2.scatter(df["sampling_date"], df["chlorophyll"], color="green", marker="+")

        # ax.plot(df["sampling_date"], sml_ratio, color="red")
        # ax.plot(df["sampling_date"], bulk_ratio, color="blue")
        # ax2.plot(df["sampling_date"], df["chlorophyll"], color="green")

    for i in range(8, 20, 1):
        lower = date(2000 + i, 3, 1)
        upper = date(2000 + i, 9, 1)
        ax.axvspan(lower, upper, color="gray", alpha=0.4)

    legend_elements = [Line2D([0], [0], marker='^', label='Bulk water',
                              markerfacecolor='blue', mew=0.3, mec="blue", aa=True, linestyle=''),
                       Line2D([0], [0], marker='o', label='Surface microlayer',
                              markerfacecolor='red', mew=0.3, mec="red", aa=True, linestyle=''),
                       Line2D([0], [0], marker='+', label='Chlorophyll a',
                              markerfacecolor='green', mew=2, mec="green", aa=True, linestyle='', markersize=10)
                       ]

    ax.legend(handles=legend_elements)
    ax.xaxis.set_major_locator(months)
    ax.xaxis.set_major_formatter(monthsFmt)
    ax.autoscale_view()
    fig.autofmt_xdate()

    ax.set_xlabel("time ")
    ax.set_ylabel(f'{params[0]}/{params[1]} ratio')
    ax2.set_ylabel("Chlorophyll a concentration/\n µg/L")
    rcParams['xtick.labelsize'] = 'small'

    plt.show()


def polar_plot_coverage(records):
    """Plots the coverage as function of the day of the year projected on a polar axis."""
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='polar')

    # ensure 0° is pointing upwards
    ax.set_theta_zero_location("N")
    # make theta grow clockwise
    ax.set_theta_direction(-1)

    # labels of the theta gridlines set as month abbreviations
    angles = (0, 31, 59.25, 90.25, 120.25, 151.25, 181.25, 212.25, 243.25, 273.25, 304.25, 334.25)
    months = ("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")

    # customize theta gridlines and labels
    lines, label = plt.thetagrids(angles, months)

    for key in records:
        records[key].sort_values("sampling_date", inplace=True)
        dates = records[key]["sampling_date"].apply(lambda x: x.timetuple().tm_yday)
        dates = scale_to_polar(dates)

        # todo: make this function more generic accepting an argument for the plot params
        c = ax.plot(dates, records[key]["norm_sml_coverage"], marker="o")
        #c = ax.plot(dates, records[key]["norm_bulk_coverage"], marker="x")


if __name__ == "__main__":
    BoknisEckExtension(new=True)
