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
            
            CONSTRAINT unique_name UNIQUE(name)""",

            "boknis_eck": """
                    id INTEGER PRIMARY KEY,
                    name TEXT,
                    specid INTEGER,
                    FOREIGN KEY (specid) REFERENCES sfg(id)
                    """,

            "gasex_sfg": """
                    id INTEGER PRIMARY KEY,
                    name TEXT,
                    spec_id INTEGER,
                    sample_id INTEGER,
                    FOREIGN KEY (spec_id) REFERENCES sfg(id),
                    FOREIGN KEY(sample_id) REFERENCES samples(id)
                   """,

            "regular_sfg": """
                    id INTEGER PRIMARY KEY,
                    name TEXT,
                    specid INTEGER,
                    surfactant TEXT,
                    sensitizer TEXT,
                    photolysis TEXT,
                    comment TEXT,
                    FOREIGN KEY (specid) REFERENCES sfg(id)
                """,

            "lt": """
                id INTEGER PRIMARY KEY,
                name TEXT,
                type TEXT,
                measured_time TIMESTAMP,
                measurer TEXT,
                time TEXT,
                area TEXT,
                apm TEXT,
                surface_pressure TEXT,
                lift_off TEXT,
                CONSTRAINT unique_name UNIQUE(name)
                """,

            "gasex_lt": """
                id INTEGER PRIMARY KEY,
                station_id INTEGER,
                sample_id INTEGER,
                name TEXT,
                type TEXT,
                measured_time TIMESTAMP,
                measurer TEXT,
                time TEXT,
                area TEXT,
                apm TEXT,
                surface_pressure TEXT,
                lift_off TEXT,
                sample_hash TEXT,
                sample_id INTEGER,
                station_hash TEXT,
                CONSTRAINT unique_name UNIQUE(name),
                FOREIGN KEY(sample_id) REFERENCES samples(id)
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

            "uv": """
                id INTEGER PRIMARY KEY,
                name TEXT,
                wavelength TEXT,
                absorbance TEXT,
                CONSTRAINT unique_name UNIQUE(name)
                """,

            "gasex_surftens": """
                id INTEGER PRIMARY KEY,
                sample_id INTEGER,
                name TEXT,
                surface_tension TEXT,
                CONSTRAINT unique_name UNIQUE(name),
                FOREIGN KEY(sample_id) REFERENCES samples(id)
                """,

            "substances": """
                id INTEGER PRIMARY KEY,
                name TEXT,
                long_name TEXT,
                molar_mass TEXT,
                sensitizing TEXT,
                CONSTRAINT unique_name UNIQUE(name)
                """,

            "stations": """
                id INTEGER PRIMARY KEY,
                hash TEXT,
                type TEXT,
                number INTEGER
                
            """,

            "samples": """
                    id INTEGER PRIMARY KEY,
                    station_id INTEGER,
                    sample_hash TEXT,
                    location TEXT,
                    type TEXT,
                    number INTEGER,
                    FOREIGN KEY(station_id) REFERENCES stations(id)
                    """

        }
        self.create_db()

        self.write_substances()

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
            try:
                self.cur.execute(c)
            except sqlite3.OperationalError as e:
                print(c)

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

    def write_substances(self):
        """A function to extract the sensitizer/surfactant metadata and pass them to the database"""
        substances = Importer.extract_substances()

        command = f'INSERT INTO substances(name, long_name, molar_mass, sensitizing)VALUES(?,?,?,?)'

        for dic in substances:
            try:
                self.cur.execute(command, (dic["abbr"], dic["name"], dic["mass"], dic["photoactive"]))
            except sqlite3.IntegrityError as e:
                print(f'Substance {dic} already in database!')

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

        # todo: create the import folder final structure
        # todo: import sfg data from the folders setting MEASURER and TYPE
        # todo: strip DATE_....._.sfg und LT_...._.dat and gain station and sample hashes
        # todo: create station/sample hashes for the tensions as well
        # todo: create a dictionary of the form "stationhash" : {stationhashdata} to collect all stations /samples from LT SFG Tension
        # todo: create the station and sample tables, check if they are consistent
        # todo: populate the tables with sfg,lt and tension using the station and sample IDs


    def import_sfg(self, folder):

        for file in os.listdir(folder):
            date, measurer = file.split(" ")

            for data in os.listdir(folder + "/" + file):
                if data.endswith(".sfg"):
                    creation_time = datetime.datetime.fromtimestamp(os.path.getmtime(folder + "/" + file + "/" + data))
                    name = date + "_" + data[:-4]
                    # tuple contains wavenumber, sfg, ir, vis in this order
                    extract = list(self.extract_sfg_data(folder + "/" + file + "/" + data))

                    for var, j in enumerate(extract):
                        extract[j] = ";".join(extract[j].astype(str))

                    dic = {"name": name, "type": folder, "measured_time": creation_time, "measurer": measurer,
                           "wavenumbers": extract[0], "sfg": extract[1], "ir": extract[2], "vis": extract[3]}

                    if "dppc" in dic["name"] or "DPPC" in dic["name"]:
                        dic["type"] = "regular"

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

    def import_lt(self, folder):
        for file in os.listdir(folder):
            measurer = file.split(" ")[1]

            for data in os.listdir(folder + "/" + file):
                if data.endswith(".dat"):
                    creation_time = datetime.datetime.fromtimestamp(os.path.getmtime(folder + "/" + file + "/" + data))
                    name = data[:-4]

                    # list contains time, area, area per molecule and pressure in this order
                    extract = list(self.extract_lt_data(folder + "/" + file + "/" + data))

                    for var, j in enumerate(extract):
                        extract[j] = ";".join(extract[j].astype(str))

                    dic = {"name": name, "type": folder, "measured_time": creation_time, "measurer": measurer,
                           "time": extract[0],
                           "area": extract[1],
                           "apm": extract[2], "pressure": extract[3]}
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

    @staticmethod
    def extract_substances():

        substances = []

        with open("substances.txt", "r") as infile:
            try:
                for line in infile:
                    name, abbreviation, mass, photoactive = line.split(";")

                    substances.append({"name": name.strip(),
                                       "abbr": abbreviation.strip(),
                                       "mass": mass.strip(),
                                       "photoactive": photoactive.strip()})
            except:
                print(line.split(";") + " is not a valid substance file line")

        return substances


class NaturalSampleExtension:
    """This class will handle station and sample hashes, generate the station SQL tables"""

    def __init__(self):
        pass

    @staticmethod
    def get_hashes(name):
        """ Extract the hashes from file names. Remember to strip the file ending as well as  the leading LT or date
        before passing namestring to this function."""

        temp = name.split("_")

        station_hash = temp[0] + temp[1]

        if temp[1][0] != "c":
            sample_hash = station_hash + temp[2] + temp[3]
        else:
            sample_hash = station_hash + temp[2]

        return {"name": name, "station_hash": station_hash, "sample_hash": sample_hash}

    @staticmethod
    def process_station_hash(hash):
        date = hash[:4]
        temp = hash[4]

        if temp in ("c", "r"):
            stype = "big"
        else:
            stype = "small"

        number = hash[5]

        return {"date": date, "station_type": stype, "station_number": number}

    @staticmethod
    def process_sample_hash(hash):
        dic = NaturalSampleExtension.process_station_hash(hash)
        location = hash[4]
        dic["location"] = location

        if location == "c":
            if "deep" in hash or "low" in hash:
                dic["sample_type"] = "deep"

        elif location == "r":
            dic["sample_type"] = hash[6]
            if dic["sample_type"] == "p":
                dic["sample_number"] = hash[7]
            elif dic["sample_type"] == "s":
                dic["sample_number"] = hash[8]

        elif location == "a":
            dic["sample_type"] = hash[6]
            dic["sample_number"] = hash[8]

        else:
            raise ValueError(f"Invalid sample hash {hash}")

        return dic







# testcode section
# todo: insert data of lt and sfg in database


# I = Importer()

ex1 = "0922_r1_p_4_1_#1"
ex2 = "0922_r1_p_4"
ex3 = "0922_c1_low_1_#1"
ex4 = "0922_a1_s1_2"

a = NaturalSampleExtension.get_hashes(ex1)
b = a["sample_hash"]

o = NaturalSampleExtension.get_hashes(ex4)["sample_hash"]
print(o)
print(NaturalSampleExtension.process_sample_hash(o))
