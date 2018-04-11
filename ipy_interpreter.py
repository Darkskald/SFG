from Classes import *
import Spec_utilities as spc
import time

class IpyInterpreter:
    """A class for interaction between the spectra database and the user by
    fetching single spectra upon matching criteria adding them to an internal
    list of SFG spectra."""

    def __init__(self):
        self.database = Finder()
        self.subset = []
        self.recover = []
        self.uvdata = spc.SpecDatabase("UV")
        self.irdata = spc.SpecDatabase("IR")
        self.ramandata = spc.SpecDatabase("Raman")
        self.planer = Planer()

        self.log_file = self.init_session_control()
    # SFG Management

    def show(self):
        """Prints all spectra contained in the current subset"""
        self.write_log_entry("show")

        tabstring = "Nr."+"\t\t"+"Date"+"\t" + "Surf." + "\t" + "spectral range" + "\t\t\t\t" + "full spectrum name"
        print(tabstring)
        print("\n")
        #todo this code snippet hast to be replaced by enumarate
        for i in range(len(self.subset)):
            print(str(i) + " : " + str(self.subset[i]))
            print("\n")

    def clear(self):
        """Removes everything from the subset"""
        self.write_log_entry("show")
        self.subset = []

    def get(self, flagstring):
        self.write_log_entry("get "+flagstring)

        """Gets spectra according to specified criteria from the database and passes them into the subset"""

        if flagstring == "bo":
            self.subset = self.database.comment_based()

        elif flagstring == "p":
            self.subset = self.database.photo_based()

        else:

            f = self.flagstring_split(flagstring)
            flag = f[0]
            options = f[1].split(",")

            if flag == "su":

                for i in options:
                    option = self.retranslate_name(i)
                    print(option)
                    collector = self.database.surfactant_based(option)
                    for j in collector:
                        self.subset.append(j)

            elif flag == "se":
                for i in options:
                    option = self.retranslate_name(i)
                    print(option)
                    collector = self.database.sensitizer_based(option)
                    for j in collector:
                        self.subset.append(j)

            elif flag == "s":

                for i in options:
                    collector = self.database.sample_based(i)
                    for j in collector:
                        self.subset.append(j)

            elif flag == "d":

                self.subset = self.database.date_based(options[0])

            elif flag == "suv":

                for i in options:
                    collector = self.database.sur_volume_based(i)
                    for j in collector:
                        self.subset.append(j)

            elif flag == "m":
                for i in options:
                    collector = self.database.measurement_based(i)
                    for j in collector:
                        self.subset.append(j)

            elif flag == "ml":
                months = options
                begin = int(months[0])
                end = int(months[1])
                self.subset = self.database.by_monthrange(begin, end)

            elif flag == "y":
                years = options
                self.subset = self.database.by_year(years)

    def flagstring_split(self, flagstring):
        """This function processes the given flagstring and returns a list of SFG objects
        which are passed through the Finder methods utilizing the flags and options"""
        f = flagstring.split(" ")
        flag = f[0]
        options = f[1]
        print("options: " + options)

        allowed = ["su", "se", "p", "c", "d", "s", "suv", "m", "ml", "y"]

        if flag not in allowed:
            print("Unknown flag")
        else:
            return (flag, options)

    def remove(self, numbers):
        """Removes items by (a list of) index(indices from the subset"""
        self.write_log_entry("remove"+str(numbers))

        options = numbers.split(",")
        to_remove = [self.subset[int(i)] for i in options]
        newlist = [i for i in self.subset if i not in to_remove]
        self.subset = newlist
        self.recover = to_remove

    def keep(self, flagstring):
        """Removes everything but the specified items. Specification is done by a list of indices"""
        self.write_log_entry("keep " + str(flagstring))
        f = flagstring
        options = f.split(",")
        new_list = [self.subset[int(i)] for i in options]
        recover = [i for i in self.subset if i not in new_list]
        self.recover = recover
        self.subset = new_list

    def plot(self, flag=False):
        """Plots an overlay of all spectra currently present in the subset"""

        s1 = "plot "+str(flag)+"\n"
        s2 = "Specs included in plot: \n"
        for spec in self.subset:
            s2 += "\t"+spec.name.full_name+"\n"
        s = s1+s2
        self.write_log_entry(s)

        p = Plotter(self.subset)
        # noinspection PySimplifyBooleanCheck
        if flag == False:
            p.simple_plot()
        elif flag == "raw":
            p.raw_plot()
        elif flag == "rawi":
            p.raw_plot_plus_ir()
        elif flag == "stack":
            p.stack_plot()

    def ref(self, flagstring):
        """Keeps a a subset of the current subset according to specified selection criteria. This function is totally
        equivalent to the get function, but is applied to the subset and not the overall database"""
        self.write_log_entry("refine "+flagstring)
        self.recover = self.subset
        f = self.flagstring_split(flagstring)
        flag = f[0]
        options = f[1]

        if flag == "su":
            options = self.retranslate_name(options)
            print(options)
            self.subset = self.database.surfactant_based(options, subset=self.subset)

        elif flag == "se":
            options = self.retranslate_name(options)
            print(options)
            self.subset = self.database.sensitizer_based(options, subset=self.subset)

        elif flag == "p":
            if options == "t":
                self.subset = self.database.photo_based(subset=self.subset)
            elif options == "f":
                self.subset = self.database.photo_based(subset=self.subset, photolyzed=False)

        elif flag == "s":
            self.subset = self.database.sample_based(options, subset=self.subset)

        elif flag == "d":
            self.subset = self.database.date_based(options, self.subset)

        elif flag == "suv":

            self.subset = self.database.sur_volume_based(options, self.subset)

        elif flag == "m":
            self.subset = self.database.measurement_based(options, self.subset)

        elif flag == "ml":
            months = options.split(",")
            print("options: ", options, type(options))
            begin = int(months[0])
            end = int(months[1])
            self.subset = self.database.by_monthrange(begin, end, self.subset)

        elif flag == "y":
            years = options.split(",")
            self.subset = self.database.by_year(years, self.subset)

    def update(self):
        """This is called when new spectra are added from the lab. It makes them available, creates the final systematic
        name and copies them to the Library"""
        self.write_log_entry("update")
        Importer()
        LibraryManager().update()

    def retranslate_name(self, stri):
        """This auxiliary function creates the long form of the name of a surfactant or sensitizer, using a text file as
        database"""
        # load the allowed surfactans and sensitizers from files
        self.Surfactants = {}
        self.Sensitizers = {}

        with open("name_info/Surfactants.txt", "r") as infile:
            for line in infile:
                collect = line.split(":")
                self.Surfactants[collect[0]] = collect[1].strip()

        with open("name_info/Sensitizers.txt", "r") as infile:
            for line in infile:
                collect = line.split(":")
                self.Sensitizers[collect[0]] = collect[1].strip()
        if stri in self.Surfactants:
            return self.Surfactants[stri]
        elif stri in self.Sensitizers:
            return self.Sensitizers[stri]
        else:
            print("Retranslation failed. Unknown expression.")

    def rec(self):
        self.write_log_entry("recover")
        """Restores the subset to one step before the last refinement(keep, refine)"""
        for i in self.recover:
            self.subset.append(i)

    def analyze_peaks(self, number=4):

        a = Analyzer(self.subset)
        a.list_peaks(number)

    # further analytics (IR/Raman/UV)

    def show_uv(self):
        """Prints the names of stored UV/Vis spectra to the console"""
        counter = 1
        for spectrum in self.uvdata.database:
            print(str(counter) + " : " + spectrum.name)
            counter += 1

    def show_ir(self):
        """Prints the names of stored IR spectra to the console"""
        counter = 1
        for spectrum in self.irdata.database:
            print(str(counter) + " : " + spectrum.name)
            counter += 1

    def show_raman(self):
        """Prints the names of stored Raman spectra to the console"""
        counter = 1
        for spectrum in self.ramandata.database:
            print(str(counter) + " : " + spectrum.name)
            counter += 1

    # section for planer

    def dn(self, number):
        """dn is abbreviation vor done"""
        self.planer.done(number)

    def ata(self, string):
        """ata is abbreviation for add_task"""
        self.planer.add_task(string)

    def st(self):
        """st is abbreviation for show tasks"""
        self.planer.show_tasks()

    def note(self, note):
        timestamp = time.strftime("%b %d %Y %H:%M:%S")
        with open("notes", "a") as outfile:
            outfile.write(timestamp+"\t"+note+"\n")

    def init_session_control(self):
        timestamp = time.strftime("%b_%d_%Y_%H_%M")
        with open("logs/"+timestamp, "w") as outfile:
            outfile.write("Logfile of IpyInterpreter session\n")
        return timestamp

    def write_log_entry(self, keyword):
        timestamp = time.strftime("%b %d %Y %H:%M:%S")
        with open("logs/"+self.log_file, "a") as outfile:
            outfile.write(timestamp + "\t" + keyword + "\n")
