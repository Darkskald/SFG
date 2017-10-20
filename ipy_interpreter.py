from Classes import *


class Ipy_Interpreter:
    """A class for interaction between the spectra datebase and the user by
    fetching single spectra upon matching criteria adding them to an internal
    list of SFG spectra."""
    def __init__(self):
        self.database = Finder()
        self.subset = []
        self.recover = []
        
    def show(self):
        for i in range(len(self.subset)):
            print (str(i)+" : "+self.subset[i].name.full_name)         
    def clear(self):
        self.subset = []
    def get(self,flagstring):
        
        if flagstring =="bo":
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
                    for i in collector:
                        self.subset.append(i)
            
            elif flag == "se":
                for i in options:
                    option = self.retranslate_name(i)
                    print(option)            
                    collector = self.database.sensitizer_based(option)
                    for i in collector:
                        self.subset.append(i)
            
            
            elif flag =="s":
                for i in options:
                    collector = self.database.sample_based(i)
                    for i in collector:
                        self.subset.append(i)
            
            elif flag =="d":
                self.subset = self.database.date_based(options[0])
            
            elif flag =="suv":
                for i in options:
                 collector = self.database.sur_volume_based(i)
                 for i in collector:
                        self.subset.append(i)
            
            elif flag == "m":
                for i in options:
                    collector = self.database.measurement_based(i)
                    for i in collector:
                        self.subset.append(i)
        
            
    def flagstring_split(self,flagstring):
        """This function processes the given flagstring and returns a list of SFG objects
        which are passed through the Finder methods utilizing the flags and options"""
        f = flagstring.split(" ")
        flag = f[0]
        options = f[1]
        
        allowed = ["su","se","p","c","d","s","suv","m"]
        
        if flag not in allowed:
            print("Unknown flag")
        else:
            return (flag,options)
        
    def remove(self,numbers):   
        options = numbers.split(",")
        to_remove = [self.subset[int(i)] for i in options]
        newlist = [i for i in self.subset if i not in to_remove]
        self.subset = newlist
        self.recover = to_remove
        
    def keep(self,flagstring):
        f = flagstring
        options = f.split(",")
        new_list = [self.subset[int(i)] for i in options]
        recover = [i for i in self.subset if i not in new_list]
        self.recover = recover
        self.subset = new_list
    def plot(self,flag=False):
        p = Plotter(self.subset)
        if flag == False:
            p.simple_plot()
        elif flag == "raw":
            p.raw_plot()
        elif flag == "rawi":
            p.raw_plot_plus_ir()
        
    def refine (self,flagstring):
        f = self.flagstring_split(flagstring)
        flag = f[0]
        options = f[1]
        
        if flag == "su":
            options = self.retranslate_name(options)
            print(options)
            self.subset = self.database.surfactant_based(options,subset=self.subset)
        
        elif flag == "se":
            options = self.retranslate_name(options)
            print(options)
            self.subset = self.database.sensitizer_based(options,subset=self.subset)
        
        elif flag == "p":
            if options == "t":
                self.subset = self.database.photo_based(subset=self.subset)
            elif options == "f":
                self.subset = self.database.photo_based(subset=self.subset,photolyzed=False)
        
        elif flag =="s":
            self.subset = self.database.sample_based(options,subset=self.subset)
        
        elif flag =="d":
            self.subset = self.database.date_based(options,self.subset)
            
        elif flag =="suv":
            
            self.subset = self.database.sur_volume_based(options,self.subset)
            
        elif flag == "m":
            self.subset = self.database.measurement_based(options,self.subset)        
            
    def update(self):
        Importer()
        Library_Manager().update()
        
        
    def retranslate_name(self,stri):
            #load the allowed surfactans and sensitizers from files
            self.Surfactants = {}
            self.Sensitizers = {}
        
            with open ("name_info/Surfactants.txt","r") as infile:
                for line in infile:
                    collect = line.split(":")
                    self.Surfactants[collect[0]] = collect[1].strip()
        
            with open ("name_info/Sensitizers.txt","r") as infile:
                for line in infile:
                    collect = line.split(":")
                    self.Sensitizers[collect[0]] = collect[1].strip()
            if stri in self.Surfactants:
                return self.Surfactants[stri]
            elif stri in self.Sensitizers:
                return self.Sensitizers[stri]
            else:
                print("Retranslation failed. Unknown expression.")    
        
    def recovery(self):
        for i in self.recover:
            self.subset.append(i)
    def analyze_peaks(self,number=4):
        a = Analyzer(self.subset)
        a.list_peaks(number)
            
            
            
      