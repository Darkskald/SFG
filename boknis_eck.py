from orm import WorkDatabaseWizard
from spectrum import SfgAverager, DummyPlotter

from datetime import timedelta
import datetime
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
        self.references = self.get_references()

        """
        don't forget the GasEx data!
        """

        if new:
            self.write_references()
            self.process()
            self.match_to_table()
            self.chlorophyll = BoknisEckExtension.prepare_chorophyll_data()
            data = self.map_to_sampling_date(self.references)
            self.persist_average_data(data[0], data[1])


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
            av = SfgAverager(to_average, references)
            if av.coverage == np.inf:
                av.benchmark()
            new_dict[key] = av
        return new_dict

    def persist_average_data(self, sml_dates, bulk_dates):

        for date in sml_dates:
            be_data_orm = self.wz.be_data()
            be_data_orm. sampling_date = date

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

        # todo: bacpopulate spectra table
        self.wz.session.commit()

    def provide_dataframe(self):
        query = self.wz.session.query(self.wz.be_data)
        df = pd.read_sql(query.statement, query.session.bind)
        return df

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

if __name__ == "__main__":
    BoknisEckExtension(new=True)