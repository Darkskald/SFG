# standard utilities
import os
import shutil
import csv
import time
import datetime
import copy
import traceback
import logging

# scientific libraries
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.lines import Line2D

import sqlite3
from scipy.signal import savgol_filter
from scipy.integrate import simps as sp


class SystematicGasExName(SystematicName):
    """Special modification of the systematic name in order to fit the requirements for the GasEx
    cruise samples from June and September 2018"""

    def __init__(self, namestring, creation_time="unknown"):

        self.full_name = namestring
        self.name = self.full_name[:-4]
        self.creation_time = creation_time

        temp = self.name.split("_")

        if len(temp) > 3:
            self.measurement_date = temp[0]
            self.station = temp[1] + "_" + temp[2]
            self.type = temp[3]

            try:
                self.number = temp[4]

            except IndexError:
                self.number = 1


        else:
            self.name.date = temp[0]
            self.station = temp[1] + "_" + temp[2]
            self.type = temp[3][0]

        self.station_hash = temp[1] + "_" + temp[2][1]
        self.sample_hash_string = "_".join(temp[1:])

    def __str__(self):

        return self.name

    def __repr__(self):
        return self.__str__()


class SampleHash:

    def __init__(self, namestring):
        self.namestring = namestring
        self.process_list = namestring.split("_")

        if len(self.process_list) not in (3, 4):
            print(namestring)
            #raise IndexError("Invalid process list length!")

        self.station_type = None
        self.station_number = None

        self.sample_type = None
        self.sample_number = None

        self.date = None
        self.station_hash = None
        self.station_construct = None

        self.date_from_name()
        self.set_station_hash()
        self.get_type()


    def __eq__(self, other):

        if self.namestring == other.namestring:
            return True
        else:
            return False

    def date_from_name(self, year=2018):
        month = int(self.process_list[0][0:2])
        day = int(self.process_list[0][2:])
        self.date = datetime.date(year, month, day)

    def set_station_hash(self):
        try:
            self.station_hash = self.process_list[0] + self.process_list[1][1]
        except IndexError:
            print(f'Invalid sample name {self.namestring}')

    def get_type(self):
        self.station_type = self.process_list[1][0].lower()
        self.station_number = self.process_list[1][1]
        self.sample_type = self.process_list[2]

        if len(self.process_list) == 3:
            self.sample_number = 1
        elif len(self.process_list) == 4:
            self.sample_number = self.process_list[-1]

        if "s" in self.sample_type:
            self.sample_type = "s"

        elif self.sample_type in ("deep", "low"):
            self.sample_type = "c"

    def get_doy(self):

        return self.date.timetuple().tm_yday


class Sample:

    def __init__(self, sample_hash):
        self.sample_hash = sample_hash
        self.station_hash = self.sample_hash.station_hash  # todo:
        self.lt_isotherms = []  # list of corresponding Isotherms
        self.sfg_spectra = []  # list of sfg spectra

        self.ch_integral = None
        self.max_pressure = None
        self.surface_tension = None

    def __repr__(self):
        pass

    def __str__(self):
        pass


class LtManager:
    """A class to perform operations on a set of LtIsotherm objects. It is an extension to the SessionControllManager
    to extend his features with isotherm handling. It relies on sqlite databases as well."""

    # todo: obsolete. The small remaining amount of code does not justify a class
    def __init__(self, database, table="lt_gasex"):

        self.database = database
        self.cursor = database.cursor()
        self.table = table
        self.isotherms = []
        self.days = None
        self.ordered_days = {}

        self.get_all_isotherms()
        #todo : the functions below have to be replaced by the sample/station hash system
        #self.join_days_isotherms()
        #self.order_by_sample()
        #self.join_same_measurement()

    def get_all_isotherms(self):
        """Fetches all isotherm data from the database and transforms them to LtIsotherm objects."""
        command = "SELECT * from " + self.table
        self.cursor.execute(command)
        result = self.cursor.fetchall()
        for i in result:
            lt = LtIsotherm(i[1], i[2], i[4], i[5], i[6], i[7])
            self.isotherms.append(lt)