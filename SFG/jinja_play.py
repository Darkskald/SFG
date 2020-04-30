import datetime
import os
import subprocess

import matplotlib.pyplot as plt
from jinja2 import Template
from sqlalchemy import func
from sqlalchemy.orm.exc import NoResultFound

from SFG.natural_samples import BEDatabaseWizard

dirname = os.path.dirname(__file__)
p = os.path.join(dirname, 'mpl_config/origin.mpltstyle')

plt.style.use(p)

tex_plate = """
\\documentclass{article}
\\usepackage[utf8]{inputenc}
\\renewcommand{\\familydefault}{\sfdefault}
\\begin{document}
\\section*{{sec1}}

\\begin{center}
\\begin{tabular}{ ccc } 
 \hline
 SML & bulk & total \\\\
 \hline
 \hline
 cell4 & cell5 & cell6 \\\\ 
 cell7 & cell8 & cell9 \\\\ 
 \hline
\end{tabular}
\end{center}




\\end{document}
"""


class SamplingDay:

    def __init__(self, sampling_day: datetime.datetime, be_spectra, be_manager: BEDatabaseWizard):
        self.sampling_day = sampling_day
        self.be_manager = be_manager
        self.be_spectra = be_spectra

        self.measurement_days = {}

        self.repr = """"""

        self.get_measurement_days()
        self.process_measurement_days()
        self.get_measurement_tex_string()

    def get_measurement_days(self):
        for spec in self.be_spectra:
            spec = b.convert_be_to_sfg(spec)
            if 0 <= spec.meta["time"].hour < 8:
                spec.meta["time"] -= datetime.timedelta(days=1)
            measurement_date = spec.meta["time"].date()

            if measurement_date not in self.measurement_days:
                self.measurement_days[measurement_date] = {"spectra": None, "references": None, "integral": None}

    def process_measurement_days(self):
        for key in self.measurement_days:
            self.measurement_days[key]["references"] = b.session.query(b.sfg).filter(b.sfg.type == "boknis_ref").filter(
                func.date(b.sfg.measured_time) == key).all()
            try:
                self.measurement_days[key]["integral"] = b.session.query(b.measurement_days).filter(
                    b.measurement_days.date == key).one().dppc_integral
            except NoResultFound:
                with open("no_dppc_reference.txt", "a") as outfile:
                    outfile.write(f'{key} has no DPPC!\n')
            self.measurement_days[key]["spectra"] = b.session.query(b.sfg).filter(b.sfg.type == "boknis").filter(
                func.date(b.sfg.measured_time) == key).all()

    def plot_references(self):
        pass

    def get_measurement_tex_string(self):
        for m in self.measurement_days:
            integral_sql = self.measurement_days[m]["integral"]
            sample_integrals = {i.name: BEDatabaseWizard.construct_sfg(i).calculate_ch_integral() for i in
                                self.measurement_days[m]["spectra"]}
            try:
                single_integrals = {i.name: BEDatabaseWizard.construct_sfg(i).calculate_ch_integral() for i in
                                    self.measurement_days[m]["references"]}

            except KeyError:
                print("NO SUCH ENTRY")

            self.repr += f'Measurement day {m}\n'
            self.repr += f'integral: {integral_sql}\n'
            if single_integrals:
                for item in single_integrals:
                    self.repr += f'{item}: {single_integrals[item]}\n'

            print(self.repr)


def produce_pdf(prefix):
    file = 'p{prefix}.tex'
    with open(file, "w") as outfile:
        outfile.write(t)
    subprocess.run(["pdflatex", "test.tex"])
    subprocess.run(["rm", "*.log", "*.aux"])


# references = BoknisEckExtension().fetch_dppc_integrals()

def process_date_dic(dic, references):
    out = []
    for key in dic:

        measurement_days = {}
        spectra = [b.convert_be_to_sfg(i) for i in dic[key]]

        for spec in dic[key]:
            measurement_date = b.convert_be_to_sfg(spec).meta["time"].date()
            if measurement_date not in measurement_days:
                measurement_days[measurement_date] = []

        for m in measurement_days:
            temp = b.session.query(b.sfg).filter(b.sfg.type == "boknis_ref").filter(
                func.date(b.sfg.measured_time) == m).all()
            measurement_days[m] = temp
            try:
                a = references[m]
                out.append(SamplingDay(key, spectra, measurement_days))
            except:
                pass

    return out


t = Template(tex_plate).render(sec1="{blubb}")

b = BEDatabaseWizard()
dic = b.get_data_per_sampling_date()
for key in dic:
    s = SamplingDay(key, dic[key], b)
    # print(s.measurement_days)
