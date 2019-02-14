import sqlite3
import os
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
            type TEXT,
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
                    sample_id INTEGER,
                    sample_hash TEXT,
                    FOREIGN KEY (name) REFERENCES sfg(name),
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
                time TEXT,
                area TEXT,
                apm TEXT,
                surface_pressure TEXT,
                lift_off TEXT,
                CONSTRAINT unique_name UNIQUE(name)
                """,

            "gasex_lt": """
                id INTEGER PRIMARY KEY,
                sample_id INTEGER,
                sample_hash TEXT,
                name TEXT,
                CONSTRAINT unique_name UNIQUE(name),
                FOREIGN KEY(sample_id) REFERENCES samples(id),
                FOREIGN KEY (name) REFERENCES gasex_lt(name)
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
                date TIMESTAMP,
                number INTEGER,
                CONSTRAINT unique_name UNIQUE(hash)
                
            """,

            "samples": """
                    id INTEGER PRIMARY KEY,
                    station_id INTEGER,
                    sample_hash TEXT,
                    location TEXT,
                    type TEXT,
                    number INTEGER,
                    FOREIGN KEY(station_id) REFERENCES stations(id),
                    CONSTRAINT unique_name UNIQUE(sample_hash)
                    """

        }
        self.create_db()

        self.write_substances()
        spec_folders = ("UV", "IR", "Raman")
        for folder in spec_folders:
            self.write_spectra(folder)


    def create_db(self):

        commands = []
        for key in self.databases:
            command = f"""
            CREATE TABLE IF NOT EXISTS {key} (
            {self.databases[key]});"""
            commands.append(command)

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
    def write_sfg_data(self, dic):
        command  = """
        INSERT INTO sfg(name, measured_time, measurer,type, wavenumbers, sfg, ir, vis)
        values(?,?,?,?,?,?,?,?)
        """
        try:
            self.cur.execute(command, (dic["name"], dic["measured_time"], dic["measurer"],
                             dic["type"], dic["wavenumbers"], dic["sfg"], dic["ir"], dic["vis"]))
        except sqlite3.IntegrityError as e:
            #print(f'Spectrum {dic["name"]} already in database!')
            pass
        self.db.commit()

    # lt
    def write_lt_data(self, dic):
            command = f"""
             INSERT INTO lt
                (
                name,
                type,
                measured_time,
                time,
                area,
                apm,
                surface_pressure
                )
                VALUES(?,?,?,?,?,?,?)
            """
            try:

                self.cur.execute(command, (dic["name"], dic["type"], dic["measured_time"], dic["time"], dic["area"], dic["apm"],
                                           dic["pressure"]))

            except sqlite3.IntegrityError as e:
                #print("Spectrum already in database!")
                pass
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
                #print(f'Surface tension {tup} already in database!')
                pass
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
                data = Importer.fetch_uv_data(path)
                x_data = "wavelength"
                y_data = "absorbance"
                database = "uv"

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
                #print("Spectrum already in database!")
                pass

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
        print(names, liftoffs)
        command = \
            f"""
            UPDATE lt
            SET lift_off=?
            WHERE
            name=?    
            """
        for i, k in enumerate(names):
            try:
                print(names[i])
                self.cur.execute(command, (liftoffs[i], names[i]))

            except sqlite3.IntegrityError as e:
                #print("Lift-off already in database!")
                pass
        self.db.commit()

    def write_substances(self):
        """A function to extract the sensitizer/surfactant metadata and pass them to the database"""
        substances = Importer.extract_substances()

        command = f'INSERT INTO substances(name, long_name, molar_mass, sensitizing)VALUES(?,?,?,?)'

        for dic in substances:
            try:
                self.cur.execute(command, (dic["abbr"], dic["name"], dic["mass"], dic["photoactive"]))
            except sqlite3.IntegrityError as e:
                #print(f'Substance {dic} already in database!')
                pass

        self.db.commit()


class Importer:

    def __init__(self):
        os.chdir("newport")
        self.wizard = SqlWizard("test.db")

        self.sfg = []
        self.lt = []
        self.tensions = []
        self.liftoffs = []

        # import sfg
        gasex = self.import_sfg("gasex_sfg")
        boknis = self.import_sfg("boknis")
        regular = self.import_sfg("regular")
        for j in (gasex, boknis, regular):
            for dic in j:
                self.wizard.write_sfg_data(dic)

        # import lt
        gasex_lt = self.import_lt("gasex_lt")
        lt = self.import_lt("lt")
        for j in (gasex_lt, lt):
            for dic in j:
                self.wizard.write_lt_data(dic)

        gasex = gasex_lt + gasex
        for dic in gasex:
            if dic["type"] in ("gasex_sfg", "gasex_lt"):
                hashdic = NaturalSampleExtension.generate_hashdic(dic["name"])

                command = """
                          INSERT INTO stations(hash,type,date,number)VALUES(?,?,?,?)
                """
                try:
                    self.wizard.cur.execute(command, (hashdic["station_hash"], hashdic["station_type"],
                                                      hashdic["date"], hashdic["station_number"]))
                except sqlite3.IntegrityError:
                    pass
                self.wizard.db.commit()

                command2 = """
                           INSERT INTO samples(sample_hash, location, type, number)VALUES(?,?,?,?)
                """
                try:
                    self.wizard.cur.execute(command2, (hashdic["sample_hash"], hashdic["location"],
                                                      hashdic["sample_type"], hashdic["sample_number"]))
                except sqlite3.IntegrityError:
                    pass

                command3 = f'SELECT id FROM samples WHERE sample_hash="{hashdic["sample_hash"]}"'
                self.wizard.cur.execute(command3)
                record = self.wizard.cur.fetchall()[0][0]
                command4 = f'INSERT INTO {dic["type"]}(sample_id, sample_hash, name)VALUES(?,?,?)'#
                try:
                    self.wizard.cur.execute(command4, (record, hashdic["sample_hash"], dic["name"]))
                except sqlite3.IntegrityError as e:
                    pass

        self.wizard.db.commit()

        self.wizard.write_liftoffs()
        self.wizard.write_surface_tension()
        self.map_samples()
        self.map_tensions()

    def import_sfg(self, folder):
        out = []

        for file in os.listdir(folder):
            date, measurer = file.split(" ")

            for data in os.listdir(folder + "/" + file):
                if data.endswith(".sfg"):
                    creation_time = datetime.datetime.fromtimestamp(os.path.getmtime(folder + "/" + file + "/" + data))
                    name = date + "_" + data[:-4]
                    # tuple contains wavenumber, sfg, ir, vis in this order
                    extract = list(self.extract_sfg_data(folder + "/" + file + "/" + data))

                    for j, var in enumerate(extract):
                        extract[j] = ";".join(extract[j].astype(str))

                    dic = {"name": name, "type": folder, "measured_time": creation_time, "measurer": measurer,
                           "wavenumbers": extract[0], "sfg": extract[1], "ir": extract[2], "vis": extract[3]}

                    if "dppc" in dic["name"] or "DPPC" in dic["name"]:
                        dic["type"] = "regular"

                    out.append(dic)

        return out

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
        out = []
        for file in os.listdir(folder):
            if file.endswith(".dat"):

                    creation_time = datetime.datetime.fromtimestamp(os.path.getmtime(folder + "/" + file))
                    name = file[:-4]

                    # list contains time, area, area per molecule and pressure in this order
                    extract = list(self.extract_lt_data(folder + "/" + file))

                    for j, var in enumerate(extract):
                        extract[j] = ";".join(extract[j].astype(str))

                        dic = {"name": name, "type": folder, "measured_time": creation_time,
                               "time": extract[0],
                               "area": extract[1],
                               "apm": extract[2], "pressure": extract[3]}
                    out.append(dic)

        return out

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
                if len(i) == 7:

                    time.append(i[1])
                    area.append(i[2])
                    apm.append(i[3])
                    surface_pressure.append(i[4])

        return np.array(time), np.array(area), np.array(apm), np.array(surface_pressure)

    def map_samples(self):
        """A function to connect the samples with their corresponding station ID"""
        command = "SELECT id,sample_hash FROM samples"
        self.wizard.cur.execute(command)
        tup = self.wizard.cur.fetchall()

        for record in tup:
            station_hash = NaturalSampleExtension.get_station_from_sample(record[1])
            command2 = f'SELECT id FROM stations WHERE hash="{station_hash}"'
            self.wizard.cur.execute(command2)
            station_id = self.wizard.cur.fetchall()[0][0]
            command3 = f'UPDATE samples SET station_id={station_id} WHERE id={record[0]}'
            self.wizard.cur.execute(command3)

        self.wizard.db.commit()

    def map_tensions(self):
        command = f'SELECT id, name FROM gasex_surftens'
        self.wizard.cur.execute(command)
        records = self.wizard.cur.fetchall()
        for record in records:
            sample_hash = NaturalSampleExtension.get_hashes(record[1])["sample_hash"]
            command2 = f'SELECT id FROM samples WHERE sample_hash="{sample_hash}"'
            self.wizard.cur.execute(command2)
            try:

                sample_id = self.wizard.cur.fetchall()[0][0]
                command3 = f'UPDATE gasex_surftens SET sample_id={sample_id} WHERE name="{record[1]}"'
                self.wizard.cur.execute(command3)
                self.wizard.db.commit()
            except IndexError:
                pass




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
                pass

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

    @staticmethod
    def refine_regular(namestring):
        """A function getting the metainformation out of the name of a regular spectrum"""
        pass


class NaturalSampleExtension:
    """This class will handle station and sample hashes, generate the station SQL tables"""

    def __init__(self):
        pass

    @staticmethod
    def get_hashes(name):
        """ Extract the hashes from file names. Remember to strip the file ending as well as  the leading LT or date
        before passing namestring to this function."""

        temp = name.split("_")

        station_hash = temp[0] + temp[1][1]

        try:
            if temp[1][0] != "c":
                sample_hash = temp[0] + temp[1] + temp[2] + temp[3]

            else:
                sample_hash = temp[0] + temp[1] + "deep"

            dic = NaturalSampleExtension.process_sample_hash(sample_hash)

            dic["name"] = name
            dic["station_hash"] = station_hash
            dic["sample_hash"] = sample_hash
        except IndexError:
            print(temp)

        return dic

    @staticmethod
    def process_station_hash(hash):
        date = hash[:4]
        date = datetime.date(2019, int(date[0:2]), int(date[2:]))
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
                dic["sample_number"] = "-"

        elif location == "r":
            dic["sample_type"] = hash[6]
            if dic["sample_type"] in ("p","P"):
                dic["sample_number"] = hash[7]
                dic["sample_type"] = "p"

            elif dic["sample_type"] == "s":
                dic["sample_number"] = hash[8]

        elif location == "a":
            dic["sample_type"] = hash[6]
            dic["sample_number"] = hash[8]

        else:
            raise ValueError(f"Invalid sample hash {hash}")

        return dic

    @staticmethod
    def generate_hashdic(namestring):

        temp = namestring.split("_")[1:]
        temp = "_".join(temp)

        dic = NaturalSampleExtension.get_hashes(temp)
        return dic

    def get_station_from_sample(sample_hash):
        temp = sample_hash[:4]+sample_hash[5]
        return temp

# Friday:
# todo: add the molar masses of the bxn species
# todo: sort the boknis eck spetra in the folders
# todo: sort the regular spectra in the folders
# todo: sort the lt isotherms in the folders
# todo: build a parser for the metainformation of the regular sfg spectra
# todo: comment all functions written so far, set todos for ugly parts of code
# todo: reorganize the scm: clear functions for getting data from database
# todo: rebase the refine operations on views
# todo: separation between the data in sql queries and the final cast into objects














# testcode section
# todo: insert data of lt and sfg in database
t = "20180628_0608_r4_p_2"

I = Importer()
I.map_samples()