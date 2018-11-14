class SpectraManager:
    """A class to handle UV, IR, and Raman specta"""

    def __init__(self, database):
        self.raman = []
        self.ir = []
        self.uv = []

    def construct_spectrum(self, item):
        """Constructs a spectrum object from the SQL query results"""
        pass

    def query(self, table):
        """Extract the data from the table"""
        pass

class SessionControlManager:
    """Successor of the IpyInterpreter. A class to access all experimental data, search data by certain match criteria
    and produce plots efficiently. Designed to work with a sqllite database containing the spectral raw data. Especially
    useful from interactive python environments (eg IPy)."""

    def __init__(self, database, id):
        # todo: handling of UV/Ir/Raman-Data has to be included
        self.db = sqlite3.connect(database, detect_types=sqlite3.PARSE_DECLTYPES | sqlite3.PARSE_COLNAMES)
        self.cur = self.db.cursor()

        self.session_id = id
        self.Surfactants = {}
        self.Sensitizers = {}
        self.Makros = {}
        self.stations = None
        self.get_senssurf_names()

        self.lt_manager = None
        self.spec_manager = None

        # former IpyInterpreter functionality, tracking the primary key in parallel
        self.subset_ids = []
        self.subset = []

        self.recover_ids = []
        self.recover = []

        self.dppc_ints = None

        self.set_spec_manager()


    # database methods

    def fetch_single(self, number, database, condition3=None, condition4=None):
        """Fetches the data of a single SFG spectrum from the database. Returns an SfgSpectrum object. The
        conditional parameters allow a selection of the spectrum according to specified properties. The
        default_data kwarg controlls which database the spectrum should be extracted from."""

        # todo: A separate table for BoknisEck data is desirable

        command = "SELECT * FROM " + database + " WHERE ROWID=" + str(number)

        if condition3 is not None:
            command += " AND " + str(condition3) + "=" + str(condition4)

        self.cur.execute(command)

        result = self.cur.fetchall()[0]
        spec = self.construct_sfg(result, database)

        return spec

    def construct_sfg(self, query_result, database):
        """A function to create an SFG spectrum from the result of a query in the SQL database"""

        result = query_result
        creationtime = result[2]

        if database == "sfg_database":

            sysname = SystematicName(result[1], creationtime)
            wavenumber = np.array(result[-7].split(";")).astype(np.float)
            sfg = np.array(result[-6].split(";")).astype(np.float)
            ir = np.array(result[-5].split(";")).astype(np.float)
            vis = np.array(result[-4].split(";")).astype(np.float)

        else:

            if "DPPC" not in result[1]:
                sysname = SystematicGasExName(result[1], creationtime)
                wavenumber = np.array(result[-4].split(";")).astype(np.float)
                sfg = np.array(result[-3].split(";")).astype(np.float)
                ir = np.array(result[-2].split(";")).astype(np.float)
                vis = np.array(result[-1].split(";")).astype(np.float)

            else:
                sysname = SystematicName(result[1], creationtime)
                wavenumber = np.array(result[-4].split(";")).astype(np.float)
                sfg = np.array(result[-3].split(";")).astype(np.float)
                ir = np.array(result[-2].split(";")).astype(np.float)
                vis = np.array(result[-1].split(";")).astype(np.float)

        spec = SfgSpectrum(wavenumber, sfg, ir, vis, sysname)

        return spec

    def general_fetch(self, condition_1=None, condition_2=None, database="sfg_database"):
        """A function to fetch spectra from the database according to specified match criteria.
        The default_data kwarg controlls the datatable that is used."""

        command = "SELECT * from " + database

        if condition_1 is not None and condition_2 is not None:
            command += (" WHERE " + condition_1 + "=" + condition_2)

        self.cur.execute(command)

        for item in self.cur.fetchall():
            id = item[0]
            self.subset_ids.append(id)

            try:
                self.subset.append(self.construct_sfg(item, database))
            except:
                logging.error(traceback.format_exc())

    def general_refine(self, condition1, condition2, database):
        """Refinement of the actual subset by applying further match criteria. This is the
        actual implementation of the get method abstracting away the database access routine from the
        user."""
        temp = []
        temp_id = []

        for id in self.subset_ids:

            try:
                s = self.fetch_single(id, database, condition1, condition2)
                temp.append(s)
                temp_id.append(id)
            except IndexError:
                pass

        self.recover = self.subset
        self.recover_ids = self.subset_ids

        self.subset = temp
        self.subset_ids = temp_id

    # user interaction functions

    def plot(self):
        """Runs an external PyQt application plotting all current spectra in the subset. The external
        app gives access to further manual analysis functionality"""
        run_app(self.subset, self.session_id)

    def show(self):
        """Prints a list of the spectra of the current subset including their subset index by which they
        can be accessed easily"""
        for i, spectrum in enumerate(self.subset):
            print(str(i) + " " + str(spectrum))

    def clear(self):
        """Deletes all spectra from the current subset"""
        self.subset = []
        self.subset_ids = []

    def get(self, flagstring, database="sfg_database", ref=False):
        """Fetches spectra according to the desired properties and adds them to the subset"""
        t = self.flagstring_split(flagstring)

        if t[0] == "su" or t[0] == "se":
            condition1 = self.retranslate_name(t[0])
            condition2 = "\"" + self.retranslate_name(t[1]) + "\""

            if ref is False:
                try:
                    self.general_fetch(condition1, condition2, database=database)

                except:
                    pass

            if ref is True:
                try:
                    self.general_refine(condition1, condition2, database)

                except:
                    pass

        elif t[0] == "name":
            self.general_fetch(t[0], "\"" + t[1] + ".sfg" + "\"")

    def ref(self, flagstring):
        """Refines the current subset by applying further match criteria"""
        self.get(flagstring, database="sfg_database", ref=True)

    def by_time(self, time1, time2, database, refine=False):
        """Fetch or  refine the spectral data by time of measurement. The Number has to be given as a string, embraced
        in '' quotation marks to pass it to the SQL query. Note this function is the only user function that directly
        accesses the SQL database"""

        if refine == False:

            command = "SELECT * from " + database + " WHERE measured_time between " + time1 + " and " + time2

            self.cur.execute(command)
            keys = []
            for item in self.cur.fetchall():
                id = item[0]
                self.subset_ids.append(id)
                self.subset.append(self.fetch_single(id))
        else:
            temp = []
            temp_id = []

            for id in self.subset_ids:

                try:
                    command = "SELECT * FROM " + database + " WHERE ROWID=" + str(
                        id) + " AND measured_time between " + time1 + " and " + time2
                    self.cur.execute(command)
                    result = self.cur.fetchall()[0]
                    s = self.construct_sfg(result, database)
                    temp.append(s)
                    temp_id.append(id)
                except IndexError:
                    pass

            self.recover = self.subset
            self.recover_ids = self.subset_ids

            self.subset = temp
            self.subset_ids = temp_id

    def rec(self):
        """A recover function that resets the subset to the state before the last refinement or keep call"""
        self.subset = self.recover
        self.subset_ids = self.recover_ids

    def remove(self, numbers):
        """Removes items by (a list of) index(indices from the subset"""

        options = numbers.split(",")
        to_remove = [self.subset[int(i)] for i in options]
        newlist = [i for i in self.subset if i not in to_remove]
        new_indexlist = [self.subset_ids[int(i)] for i in options]
        self.recover = self.subset
        self.recover_ids = self.subset_ids

        self.subset = newlist
        self.subset_ids = new_indexlist

    def keep(self, flagstring):
        """Removes everything but the specified items. Specification is done by a list of indices"""
        f = flagstring
        options = f.split(",")
        new_list = [self.subset[int(i)] for i in options]
        new_indexlist = [self.subset_ids[int(i)] for i in options]
        recover = [i for i in self.subset if i not in new_list]

        self.recover = recover
        self.recover_ids = self.subset_ids

        self.subset = new_list
        self.subset_ids = new_indexlist

    # specific functionality for GasEx

    def set_lt_manager(self):
        """Creates a LtManager object to give access to analysis tools for Langmuir trough isotherms besides
        the built-in SFG functionality. It can be seen as an add-in to the SessionControlManager to extend
        his capabilities."""

        self.lt_manager = LtManager(self.db)


    def collect_stations(self):
        """Creates Station objects from the station information of the currently available isotherms and
        spectra. This helps to keep track of the set of samples taken within one GasEx sampling station"""
        # todo: include SFG stations here as well as in the match to station function
        # todo: include error handling (LtManager defined?)
        stations = []
        for isotherm in self.lt_manager.isotherms:

            temp = isotherm.long_station
            if temp not in stations:
                stations.append(Station(temp))
        self.stations = {s.station_hash: s for s in stations}

    def collect_samples(self):


    def match_to_stations(self):
        """Matches the LtIsotherms in the LtManager to the generated station list. The
        stations attribute of the SessionControllManager is a dictionary with the station names
        as keys and a list with LtIsotherms matching the station as value"""

        for isotherm in self.lt_manager.isotherms:
            self.stations[isotherm.station_hash].lt_isotherms.append(isotherm)

        for spectrum in self.subset:
            if isinstance(spectrum.name, SystematicGasExName):

                try:
                    self.stations[spectrum.name.station_hash].sfg_spectra.append(spectrum)
                except KeyError:
                    self.stations[spectrum.name.station_hash] = Station(spectrum.name.station)
                    self.stations[spectrum.name.station_hash].sfg_spectra.append(spectrum)

    def get_station_numbers(self):
        """Traverse the stations and assign them a number in chronological order (1-n)"""

        stations = [i for i in self.stations.values()]
        stations.sort()
        for i, s in enumerate(stations):
            s.station_number = i + 1

    def fetch_gasex_sfg(self):
        """Fetches all SfgSpectra from the GasEx cruise and puts them in the subset attribute"""
        self.general_fetch(database="sfg_gasex")

    def setup_for_gasex(self):
        """Convenience function, calls all necessary functions to setup the SCM for gasex data
        processing"""
        self.set_lt_manager()
        self.fetch_gasex_sfg()
        self.collect_stations()
        self.match_to_stations()
        self.get_station_numbers()
        self.dppc_ints = self.get_dppc_average()
        for s in self.stations.values():
            s.join_samples()
            s.count_per_type()

    # handling DPPC normalization
    def get_dppc_average(self):
        """Takes a list of SfgSpectra as arguments. Returns a dictionary of dates as keys and the average SFG CH
        integral as values."""

        temp = {}
        out = {}
        #todo: erase the error which arises if the spectrum is not measured till 3000 cm
        for spec in self.subset:  # type: SfgSpectrum

            if not (isinstance(spec.name, SystematicGasExName)):

                if spec.name.surfactant == "DPPC":

                    if np.max(spec.wavenumbers >= 3000):

                        try:
                            temp[spec.name.creation_time.date()].append(spec)
                        except KeyError:
                            temp[spec.name.creation_time.date()] = [spec]



        for key in temp:

            intens = 0
            counter = 0

            for spec in temp[key]:
                intens += spec.calculate_ch_integral()
                counter += 1

            intens /= counter
            out[key] = intens

        return out

    # auxiliary functions

    def flagstring_split(self, flagstring):
        """This function processes the given flagstring and returns a list of SFG objects
        which are passed through the Finder methods utilizing the flags and options"""
        f = flagstring.split(" ")
        return f

    def get_senssurf_names(self):
        """Loads allowed names of surfactants and sensitizers from the control file"""

        with open("name_info/Surfactants.txt", "r") as infile:
            for line in infile:
                collect = line.split(":")
                self.Surfactants[collect[0]] = collect[1].strip()

        with open("name_info/Sensitizers.txt", "r") as infile:
            for line in infile:
                collect = line.split(":")
                self.Sensitizers[collect[0]] = collect[1].strip()

        with open("name_info/makros.txt", "r") as infile:
            for line in infile:
                collect = line.split(":")
                self.Makros[collect[0]] = collect[1].strip()

    def retranslate_name(self, stri):
        """Extracts name information from the string, returning the full name of the substance abbreviation"""

        if stri in self.Surfactants:
            return self.Surfactants[stri]
        elif stri in self.Sensitizers:
            return self.Sensitizers[stri]
        elif stri in self.Makros:
            return self.Makros[stri]
        else:
            print("Retranslation failed. Unknown expression.")

    def set_spec_manager(self):
        self.spec_manager = SpectraManager(self.db)

class LtManager:
    """A class to perform operations on a set of LtIsotherm objects. It is an extension to the SessionControllManager
    to extend his features with isotherm handling. It relies on sqlite databases as well."""

    def __init__(self, database, table="lt_gasex"):

        self.database = database
        self.cursor = database.cursor()
        self.table = table
        self.isotherms = []
        self.days = None
        self.ordered_days = {}

        self.get_all_isotherms()
        self.join_days_isotherms()
        self.order_by_sample()
        self.join_same_measurement()

    def get_all_isotherms(self):
        """Fetches all isotherm data from the database and transforms them to LtIsotherm objects."""
        command = "SELECT * from " + self.table
        self.cursor.execute(command)
        result = self.cursor.fetchall()
        for i in result:
            lt = LtIsotherm(i[1], i[2], i[4], i[5], i[6], i[7])
            self.isotherms.append(lt)

    def get_days(self):
        """Extract the days on which the samples were taking from the isotherms."""
        days = []
        for isotherm in self.isotherms:
            if isotherm.day not in days:
                days.append(isotherm.day)
        return days

    def join_days_isotherms(self):
        """Matches the isotherms to the days of measurement, storing this mapping as a dictionary
        in the days attribute."""

        days = self.get_days()
        self.days = {i: [j for j in self.isotherms if j.day == i] for i in days}

    def order_by_sample(self):
        """Matches isotherms to a certain sample, ensuring that all consecutive measurements of
        the same sample are stored together."""

        for item in self.days:
            stations = []

            for isotherm in self.days[item]:
                if isotherm.station not in stations:
                    stations.append(isotherm.station)

            self.ordered_days[item] = {i: [j for j in self.days[item] if j.station == i] for i in stations}

        for day in self.ordered_days:
            for station in self.ordered_days[day]:

                types = []
                for isotherm in self.ordered_days[day][station]:
                    if isotherm.type not in types:
                        types.append(isotherm.type)
                self.ordered_days[day][station] = {i: [j for j in self.ordered_days[day][station] if j.type == i] for i
                                                   in types}

    def join_same_measurement(self):
        """Stores information about the other measurements of a single sample in the isotherms partner
        attribute."""
        for i, isotherm in enumerate(self.isotherms):

            for isotherm2 in self.isotherms[i + 1:]:
                if isotherm.same_sample(isotherm2) and isotherm2 not in isotherm.partners:
                    isotherm.partners.append(isotherm2)