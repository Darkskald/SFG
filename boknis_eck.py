from orm import WorkDatabaseWizard
from spectrum import SfgAverager, DummyPlotter

from sqlalchemy.exc import IntegrityError
from sqlalchemy.orm.exc import NoResultFound
from datetime import timedelta
import datetime
import functools
import pandas as pd
import re
import numpy as np

import matplotlib.pyplot as plt


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
            self.persist_average_data(data[0], data[1])
            self.include_gasex_average()

        self.references = self.fetch_dppc_integrals()

    # setting up the data
    def get_references(self):

        results = self.wz.session.query(self.wz.sfg).filter(self.wz.sfg.type == "boknis_ref").all()
        dates = {}

        for spec in results:
            spec = self.wz.construct_sfg(spec)

            if 0 < spec.meta["time"].hour < 8:
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
        measurement_days = self.wz.session.query(self.wz.measurement_days).all()
        out = {i.date: i.dppc_integral for i in measurement_days}
        return out

    def fetch_gasex_references(self):
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

                elif row["Sampler no."] in (1, 2):
                    q.sample_type = "sml"

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
        deep = self.add_spectrum(self.get_boknis_specs().filter(self.wz.boknis_eck.sample_type == 'deep'))

        out = []
        for item in (sml, deep):
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

    def persist_average_data(self, sml_dates, bulk_dates):

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

            finally:
                self.wz.session.add(be_data_orm)
                self.wz.session.commit()

        # todo: bacpopulate spectra table

    def provide_dataframe(self):
        query = self.wz.session.query(self.wz.be_data)
        df = pd.read_sql(query.statement, query.session.bind)
        return df

    # include gasex measurements
    def calculate_gasex_average(self):

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
        be["chlora"] = be["chlora"].mask(be["chlora"] < 0)
        be = be[be["chlora"].notnull()]
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


class BEDatabaseWizard(WorkDatabaseWizard):
    def __init__(self):
        super().__init__()
        self.df = BoknisEckExtension().provide_dataframe()

        # convert to datetime
        self.df["sampling_date"] = pd.to_datetime(self.df["sampling_date"])

        # remove outliers
        self.df["sml_coverage"] = self.df["sml_coverage"].mask(self.df["sml_coverage"] > 1)
        self.df["bulk_coverage"] = self.df["bulk_coverage"].mask(self.df["bulk_coverage"] > 1)
        self.df["sml_coverage"] = self.df["sml_coverage"].replace([np.inf, -np.inf], np.nan)
        self.df["bulk_coverage"] = self.df["bulk_coverage"].replace([np.inf, -np.inf], np.nan)

        # remove nans in coverage
        self.df = self.df[self.df["sml_coverage"].notna()]
        self.df = self.df[self.df["bulk_coverage"].notna()]

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

    def get_be_spectra_monthwise(self):
        out = {}
        temp = self.session.query(self.boknis_eck).all()
        for item in temp:
            if item.sampling_date.month not in out:
                out[item.sampling_date.month] = [self.fetch_by_specid(item.specid)]
            else:
                out[item.sampling_date.month].append(self.fetch_by_specid(item.specid))

        return out


if __name__ == "__main__":
    BoknisEckExtension(new=True)