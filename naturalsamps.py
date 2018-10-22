from Classes import SfgSpectrum, LtIsotherm


class Sample:
    def __init__(self, sample_hash):
        pass

    def __repr__(self):
        pass

    def __str__(self):
        pass

class Station:
    """A class carrying all information and data for a given cruise station, especially SFG and isotherms"""

    def __init__(self, name):

        self.name = name
        self.station_hash = None
        self.date, self.id = self.name.split("_")
        self.cruise_day = int(self.date[2:]) - 2
        self.sfg_spectra = []
        self.lt_isotherms = []
        self.station_number = None
        self.lt_joined = {}
        self.stats = {
            "positive_plate": 0,
            "positive_screen": 0,
            "total_screen": 0,
            "total_plate": 0,
            "screen_av": 0,
            "plate_av": 0,
            "total": 0,
            "std_plate": 0,
            "std_screen": 0,
            "std_total": 0,
            "percent_plate": 0,
            "percent_screen": 0,
            "total_percent": 0,
            "total_av": 0,
        }
        self.isotherm_count = None

        self.set_station_hash()

    def __lt__(self, other):
        """Determines which of two stations took place earlier."""

        daytag = int(self.name.split("_")[0][2:])
        monthtag = int(self.name.split("_")[0][0:2])
        stationtag = int(self.id[1])

        daytag2 = int(other.name.split("_")[0][2:])
        monthtag2 = int(other.name.split("_")[0][0:2])
        stationtag2 = int(other.id[1])

        outbool = None

        if monthtag <= monthtag2:

            if daytag == daytag2:
                if stationtag < stationtag2:
                    outbool = True

            elif daytag < daytag2:
                outbool = True

            else:
                outbool = False
        else:
            outbool = False

        return outbool

    def __rpr__(self):
        return f'Station {self.id[1]} on date {self.date}'

    def __str__(self):
        return f'Station {self.id[1]} on date {self.date}'

    def set_station_hash(self):
        temp = self.name.split("_")
        try:
            self.station_hash = temp[0] + "_" + temp[1][1]
        except IndexError:
            print("*" * 80 + "\n" + str(temp))
            print(self.name)

    def count_per_type(self):
        """Counts the occurence of plate and screen samples as well as how many of those are positive
        with respect to surfactants (high surface pressure in Langmuir trough measurement). Stores all the
        statistical information in the stats dictionary."""

        # todo: handle CTD samples!
        # todo: implement a similar sample for SFG

        isos = [i[0] for i in self.lt_joined.values()]
        plate_array = []
        screen_array = []

        for isotherm in isos:  # type: LtIsotherm

            p_max = isotherm.get_maximum_pressure()

            if p_max < 73:
                if isotherm.type == "p":

                    if p_max > 1.5:
                        self.stats["positive_plate"] += 1

                    self.stats["total_plate"] += 1
                    self.stats["plate_av"] += p_max
                    plate_array.append(p_max)

                elif isotherm.type[0] == "s":

                    if 72 > p_max > 1.5:
                        self.stats["positive_screen"] += 1

                    self.stats["total_screen"] += 1
                    self.stats["screen_av"] += p_max
                    screen_array.append(p_max)

        self.stats["total"] = (self.stats["total_screen"] + self.stats["total_plate"])

        try:
            self.stats["total_av"] = (self.stats["screen_av"] + self.stats["plate_av"]) / self.stats["total"]
        except ZeroDivisionError:
            pass
        try:
            self.stats["screen_av"] /= self.stats["total_screen"]
        except ZeroDivisionError:
            pass
        try:
            self.stats["plate_av"] /= self.stats["total_plate"]
        except ZeroDivisionError:
            pass
        try:
            self.stats["percent_screen"] = self.stats["positive_screen"] / self.stats["total_screen"] * 100
        except ZeroDivisionError:
            pass
        try:
            self.stats["percent_plate"] = self.stats["positive_plate"] / self.stats["total_plate"] * 100
        except ZeroDivisionError:
            pass
        try:
            self.stats["total_percent"] = (self.stats["positive_screen"] + self.stats["positive_plate"]) / self.stats[
                "total"] * 100
        except ZeroDivisionError:
            pass
        total_array = plate_array + screen_array

        self.stats["std_plate"] = np.std(plate_array)
        self.stats["std_screen"] = np.std(screen_array)
        self.stats["std_total"] = np.std(total_array)

    def join_samples(self):
        """Joins Langmuir trough measurements of the same sample. Much more comprehensive than
        the approach in the LtManager."""

        for isotherm in self.lt_isotherms:

            if isotherm.create_sample_hash() not in self.lt_joined:
                self.lt_joined[isotherm.create_sample_hash()] = [isotherm]
            else:
                self.lt_joined[isotherm.create_sample_hash()].append(isotherm)

        for isolist in self.lt_joined.values():
            isolist.sort(reverse=True)

        self.isotherm_count = len(self.lt_joined)

    def print_stats(self):
        """Formatted output of the stations stats, calculated from the LtIsotherms belonging to the
        station."""

        for item in self.stats:
            s = f'{item} : {self.stats[item]}\n'
            print(s)

    def make_average_sfg(self, dppc=False):

        to_av = []
        for spec in self.sfg_spectra:
            if isinstance(spec.name, SystematicGasExName):
                if "deep" not in spec.name.type and "low" not in spec.name.type:
                    to_av.append(spec)

        if dppc == False:

            if len(to_av) != 0:
                out = to_av[0]
                for spec in to_av[1:]:
                    out += spec
                out.name.full_name = self.name
                out = out.calculate_ch_integral(average="gernot")

            else:
                out = None



        else:
            if len(to_av) != 0:
                dates = dppc
                temp= []

                for spec in to_av: # type: SfgSpectrum
                    dppc_integral = dates[spec.name.creation_time.date()]
                    ch_integral = spec.calculate_ch_integral(average="gernot")
                    temp.append(np.sqrt((ch_integral/dppc_integral)))

                out = np.average(temp), np.std(temp)
            else:
                out = None

        return out

    def get_overview(self, sfg=True, lt=True):
        """A function returning a formatted string of information about the station. Can be used to
        check wether everything was imported correctly"""
        separator = "-" * 80
        outstring = f'Station number {self.station_number} with name {self.name} '
        outstring += f'on {self.date[0:2]}.{self.date[2:]}.2018\n'
        outstring += f'The station contains {len(self.sfg_spectra)} SFG measurements '
        outstring += f' and  {len(self.lt_isotherms)} compression isotherm numbers.\n'

        if sfg == True:
            outstring += f'{separator}\n List of the SFG measurements:\n'
            for i, spectrum in enumerate(self.sfg_spectra):
                temp = f'{i+1}: {spectrum.name.full_name}\n'
                outstring += temp

        if lt == True:
            outstring += f'{separator}\n List of the LT measurements:\n'
            for i, spectrum in enumerate(self.lt_isotherms):
                temp = f'{i+1}: {spectrum.name}\n'
                outstring += temp
        outstring += separator + "\n"

        return outstring

class Cruise:
    def __init__(self, hash):
        self.hash = hash

    def __repr__(self):
        pass

    def __str__(self):
        pass


class SmlSfgSpectrum(SfgSpectrum):
    pass

class SmlLtIsotherm(LtIsotherm):
    pass

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

    def __str__(self):

        return self.name

    def __repr__(self):
        return self.__str__()