import sqlite3
import os
import shutil
import datetime
import numpy as np
import csv


class SqlWizard:

    def __init__(self, name):
        self.db = sqlite3.connect(name)
        self.cur = self.db.cursor()
        self.databases = {
            "sfg": """
            id INTEGER PRIMARY KEY,
            name TEXT,
            measured_time TIMESTAMP,
            measurer TEXT,
            wavenumbers TEXT,
            sfg TEXT,
            ir TEXT,
            vis TEXT,
            surfactant TEXT,
            sensitizer TEXT,
            photolysis TEXT,
            CONSTRAINT unique_name UNIQUE(name)""",

            "lt":
            """
            id INTEGER PRIMARY KEY,
            name TEXT,
            measured_time TIMESTAMP,
            measurer TEXT,
            time TEXT,
            area TEXT,
            apm TEXT,
            surface_pressure TEXT,
            lift_off TEXT,
            CONSTRAINT unique_name UNIQUE(name)
            """,

            "ir":
                """
                id INTEGER PRIMARY KEY,
                name TEXT,
                wavenumbers TEXT,
                transmission TEXT,
                CONSTRAINT unique_name UNIQUE(name)
                 """,

            "raman":
                """
                id INTEGER PRIMARY KEY,
                name TEXT,
                wavenumbers TEXT,
                intensity TEXT,
                CONSTRAINT unique_name UNIQUE(name)
                """,

            "uv":
                """
                id INTEGER PRIMARY KEY,
                name TEXT,
                wavelength TEXT,
                absorbance TEXT,
                CONSTRAINT unique_name UNIQUE(name)
                """,

            "gasex_surftens":
                """
                id INTEGER PRIMARY KEY,
                name TEXT,
                surface_tension TEXT,
                CONSTRAINT unique_name UNIQUE(name)
                """

        }
        self.create_db()

    def create_db(self):
        commands = []
        for key in self.databases:
            command = f"""
            CREATE TABLE IF NOT EXISTS {key} (
            {self.databases[key]});"""
            commands.append(command)

        print(commands)
        self.exec_sql(commands)

    def exec_sql(self, commands):

        """A function to perform a number of SQL commands and commit them."""
        for c in commands:
            self.cur.execute(c)

        self.db.commit()

    # sfg
    def write_sfg_data(self, diclist):
        commands = []
        pass

    # lt
    def write_lt_data(self, diclist):
        for dic in diclist:
            command = f"""
             INSERT INTO lt_gasex
                (
                name,
                measured_time,
                time,
                area,
                apm,
                surface_pressure,
                measurer
                )
                VALUES(?,?,?,?,?,?);
            """
            try:

                self.cur.execute(command, (dic["name"], dic["measured_time"], dic["time"], dic["area"], dic["apm"],
                                           dic["surface_pressure"], dic["measurer"]))

            except sqlite3.IntegrityError as e:
                print("Spectrum already in database!")
        self.db.commit()

    # surface tension
    def write_surface_tension(self):
        out = []

        with open("gasex_surftens.txt", "r") as infile:

            for line in infile:

                try:
                    temp = line.strip().split(";")
                    out.append([temp[0], temp[1]])

                except IndexError:
                    pass

        db = self.db
        cur = self.cur
        command = f'INSERT INTO gasex_surftens(name, surface_tension)VALUES(?,?);'
        for tup in out:
            try:

                cur.execute(command, (tup[0], tup[1]))

            except sqlite3.IntegrityError as e:
                print(f'Surface tension {tup} already in database!')
        db.commit()

    # ir raman uv
    def write_spectra(self, target_folder):
        db = self.db
        cur = self.cur

        for file in os.listdir(target_folder):

            path = target_folder + "/" + file

            if target_folder in ("IR", "Raman"):
                data = Importer.fetch_iraman_data(path)
                x_data = "wavenumbers"

                if target_folder == "Raman":
                    y_data = "intensity"
                    database = "raman"

                elif target_folder == "IR":
                    y_data = "transmission"
                    database = "ir"

            elif target_folder == "UV":
                data = Importer.fetch_uv_dat(path)
                x_data = "wavelength"
                y_data = "absorbance"
                database = "ir"

            command = \
                f"""
                INSERT INTO {database}
                (
                name,
                {x_data},
                {y_data}
                )
                VALUES(?,?,?);
                """
            try:

                cur.execute(command, (str(file), str(data[0]), str(data[1])))

            except sqlite3.IntegrityError as e:
                print("Spectrum already in database!")

        db.commit()


    def write_liftoffs(self):

        names = []
        liftoffs = []
        with open("liftoff_points.csv") as csvfile:
            csvreader = csv.reader(csvfile, delimiter=";")
            next(csvfile)
            for row in csvreader:
                names.append(row[0])
                liftoffs.append(float(row[1]))

        command = \
            f"""
            UPDATE lt
            SET lift_off = ?
            WHERE
            name = ?;
            """
        for i, k in enumerate(names):
            try:

                self.cur.execute(command, (names[i], liftoffs[i]))

            except sqlite3.IntegrityError as e:
                print("Lift-off already in database!")
        self.db.commit()




class Importer:

    def __init__(self):
        self.wizard = SqlWizard("test.db")

        self.sfg = []
        self.lt = []
        self.tensions = []
        self.liftoffs = []

        self.ir = []
        self.raman = []
        self.uv = []

        self.station_meta = []

        self.import_sfg()
        self.import_lt()

    def import_sfg(self):

        for file in os.listdir("test"):
            measurer = file.split(" ")[1]

            for data in os.listdir("test/" + file):
                if data.endswith(".sfg"):
                    creation_time = datetime.datetime.fromtimestamp(os.path.getmtime("test/" + file + "/" + data))
                    name = data[:-4]
                    wavenumbers, sfg, ir, vis = self.extract_sfg_data("test/" + file + "/" + data)

                    dic = {"name": name, "measured_time": creation_time, "measurer": measurer,
                           "wavenumbers": wavenumbers, "sfg": sfg, "ir": ir, "vis": vis}
                    self.sfg.append(dic)

    def extract_sfg_data(self, file):
        data_collect = []

        with open(file, "r") as infile:

            readCSV = csv.reader(infile, delimiter="\t")
            for row in readCSV:
                data_collect.append(row)

            data_package = [0] * len(data_collect[0])
            for i in range(len(data_collect[0])):
                # mind this complex list comprehension
                data_package[i] = [j[i] for j in data_collect]
                # remove useless column
            del data_package[2]

            # convert strings to float and generate numpy array
            convert = []
            for i in data_package:
                q = [float(j) for j in i]
                convert.append(q)
            convert = [np.array(i) for i in convert]

            return convert

    def import_lt(self):
        for file in os.listdir("test"):
            measurer = file.split(" ")[1]

            for data in os.listdir("test/" + file):
                if data.endswith(".dat"):
                    creation_time = datetime.datetime.fromtimestamp(os.path.getmtime("test/" + file + "/" + data))
                    name = data[:-4]
                    time, area, apm, pressure = self.extract_lt_data("test/" + file + "/" + data)

                    dic = {"name": name, "measured_time": creation_time, "measurer": measurer, "time": time,
                           "area": area,
                           "apm": apm, "pressure": pressure}
                    print(dic)
                    self.lt.append(dic)

    def extract_lt_data(self, file):

        with open(file, "r") as infile:

            collector = []
            for line in infile:

                if line[0] != "#":
                    temp = line.strip().split("\t")
                    collector.append(temp)
            time = []
            area = []
            apm = []
            surface_pressure = []

            for i in collector:
                time.append(i[1])
                area.append(i[2])
                apm.append(i[3])
                surface_pressure.append(i[4])

        return np.array(time), np.array(area), np.array(apm), np.array(surface_pressure)

    # ir raman uv
    @staticmethod
    def fetch_iraman_data(filename):

        with open(filename, "r") as infile:
            # function can handle IR files as well as Raman files
            c = csv.reader(infile, delimiter="\t")
            wavenumbers = []
            intensities = []
            for line in c:
                wavenumbers.append(line[0])
                intensities.append(line[1])
            return wavenumbers, intensities

    @staticmethod
    def fetch_uv_data(filename):

        with open(filename, "r") as infile:
            try:
                next(infile)
                next(infile)
            except:
                pass
            c = csv.reader(infile, delimiter="\t")
            wavelength = []
            intensity = []
            for line in c:
                l = line[0].split(",")
                if len(l) == 2:
                    wavelength.append(float(l[0]))
                    intensity.append(float(l[1]))
        return wavelength, intensity

    def write_spectra(self, target_folder):
        db = self.wizard.db
        cur = self.wizard.curc

        for file in os.listdir(target_folder):

            path = target_folder + "/" + file

            if target_folder in ("IR", "Raman"):
                data = Importer.fetch_iraman_data(path)
                x_data = "wavenumbers"

                if target_folder == "Raman":
                    y_data = "intensity"

                elif target_folder == "IR":
                    y_data = "transmission"

            elif target_folder == "UV":
                data = Importer.fetch_uv_dat(path)
                x_data = "wavelength"
                y_data = "absorbance"

            command = \
                f"""
                INSERT INTO {database}
                (
                name,
                {x_data},
                {y_data}
                )
                VALUES(?,?,?);
                """
            try:

                cur.execute(command, (str(file), str(data[0]), str(data[1])))

            except sqlite3.IntegrityError as e:
                print("Spectrum already in database!")

        db.commit()
        db.close()


# testcode section
# todo: insert data of lt and sfg in database


I = Importer()
