# module-internal imports
from SFG.legacy.new_gui import run_app
from SFG.legacy.sfg import SystematicName
from spectools import SpectraManager
from SFG.natural_samples.gasex import SampleHash, SystematicGasExName, LtManager, SfgSpectrum, Station

# standard utilities
import traceback
import logging

# scientific libraries
import numpy as np

import sqlite3


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
        self.tensions = None

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

    # todo: put this to importing process
    def collect_stations(self):
        """Creates Station objects from the station information of the currently available isotherms and
        spectra. This helps to keep track of the set of samples taken within one GasEx sampling station"""
        station_hashes = []
        types = []

        for spec in self.subset: #type: SfgSpectrum

            if isinstance(spec.name, SystematicGasExName) is True:
                h = spec.get_sample_hash()
                _hash = h.station_hash

                if _hash not in station_hashes:
                    station_hashes.append(_hash)

                    if h.station_type == "a":
                        types.append("small")

                    elif h.station_type in ("r", "c"):
                        types.append("big")

                    else:
                        print("Unknwonwn station type!")

        for isotherm in self.lt_manager.isotherms:#type: LtIsotherm
            _hash = isotherm.sample_hash.station_hash
            if _hash not in station_hashes:
                station_hashes.append(_hash)
                if isotherm.sample_hash.station_type == "a":
                    types.append("small")

                elif isotherm.sample_hash.station_type in ("r", "c"):
                    types.append("big")

                else:
                    print("Unknwonwn station type!")

        for tension in self.tensions:
            try:
                h = SampleHash(tension)
                _hash = h.station_hash

                if _hash not in station_hashes:
                    station_hashes.append(_hash)
                    if h.station_type == "a":
                        types.append("small")

                    elif h.station_type in ("r", "c"):
                        types.append("big")

                    else:
                        print("Unknwonwn station type!")

            except IndexError:
                print(tension)

        self.stations = {}
        for s in zip(station_hashes, types):
            self.stations[s[0]] = Station(s[0], s[1], parent=self)

    #todo: put this to importing process and map the samples in the database to the stations
    def map_to_stations(self):

        for spectrum in self.subset:
            if isinstance(spectrum.name, SystematicGasExName) is True:
                self.stations[spectrum.get_sample_hash().station_hash].sfg_spectra.append(spectrum)

        for tension in self.tensions:
            try:
                self.stations[SampleHash(tension).station_hash].tensions.append([tension, self.tensions[tension]])
            except IndexError:
                print(tension)


        for lt_isotherm in self.lt_manager.isotherms:
            self.stations[lt_isotherm.sample_hash.station_hash].lt_isotherms.append(lt_isotherm)

    def fetch_gasex_sfg(self):
        """Fetches all SfgSpectra from the GasEx cruise and puts them in the subset attribute"""
        self.general_fetch(database="sfg_gasex")

    #todo: reimplement this to realize it via database access only
    def setup_for_gasex(self):
        """Convenience function, calls all necessary functions to setup the SCM for gasex data
        processing"""
        self.set_lt_manager()
        self.fetch_gasex_sfg()
        self.tensions = self.fetch_tension_data()
        self.collect_stations()
        self.map_to_stations()
        self.dppc_ints = self.get_dppc_average()
        liftoffs = self.fetch_liftoff_points()

        for station in self.stations.values():
                station.arange_to_sample()
                station.only_slowest()
                for isotherm in station.lt_isotherms:
                    if isotherm.name in liftoffs:
                        isotherm.lift_off = float(liftoffs[isotherm.name])
                station.analyze_station_data()

    def fetch_tension_data(self):
        out = {}
        command = f'SELECT * from gasex_surftens'
        self.cur.execute(command)
        for item in self.cur.fetchall():
            out[item[1]] = float(item[2])
        return out

    def fetch_liftoff_points(self):
        out = {}
        command = f'SELECT * from gasex_liftoff'
        self.cur.execute(command)
        for item in self.cur.fetchall():
            out[item[1]] = float(item[2])
        return out


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

    #todo: remove all name processing stuff and auxialiary file dependencies
    def get_senssurf_names(self):
        """Loads allowed names of surfactants and sensitizers from the control file"""

        with open("../../name_info/Surfactants.txt", "r") as infile:
            for line in infile:
                collect = line.split(":")
                self.Surfactants[collect[0]] = collect[1].strip()

        with open("../../name_info/Sensitizers.txt", "r") as infile:
            for line in infile:
                collect = line.split(":")
                self.Sensitizers[collect[0]] = collect[1].strip()

        with open("../../name_info/makros.txt", "r") as infile:
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

    #todo should be obsolete
    def set_spec_manager(self):
        self.spec_manager = SpectraManager(self.db)

