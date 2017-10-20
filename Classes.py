import os
import shutil
import csv
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
from scipy.signal import savgol_filter
rcParams['mathtext.default'] = 'regular'


#tools to import new spectra v
class Importer:
    #first: make a list of day folders in the archive directory
    def __init__(self):
        self.new_folders = []
        for folder in os.listdir("archive"):
            self.new_folders.append(folder)

        #now: get the files of each folder, strip FL off, add daytag to filename
        self.new_files = []
        for folder in self.new_folders:
            daytag = (folder.split("FL"))[0].strip()
            #now creating the day_information_file
            day_info_file = daytag+ ".dif"
            if self.gate_keeper("library/day_information",day_info_file) == False:
                with open ("library/day_information/"+day_info_file,"w") as outfile:
                    outfile.write(daytag+"\n")
                    files = []
                    for file in os.listdir("archive/"+folder):
                        files.append(file)
                    outfile.write(str(len(files))+"\n")
                    for file in files:
                        outfile.write(file+"\n")
                    outfile.write("#")

            #next: collect all spectra per day and copy them to the library 
            for file in os.listdir("archive/"+folder):
                new_filename = daytag+"_"+ file
                #the gate_keeper function prohibits double import
                if self.gate_keeper("library",new_filename) == False:
                    shutil.copy2("archive/"+folder+"/"+file,"library/"+new_filename)


    def gate_keeper(self,checkdir,filename):
        #checks if filename is already present in the directory chekdir
        checklist = []
        for file in os.listdir(checkdir):
           checklist.append(file)
        if filename in checklist:
            print("File "+filename+" already exists")
            return True
        else:
            return False
            
#tools to maintain the library
class Day_meta:
    pass
class Library_Manager:
    #controll and maintain the library management file
    def __init__(self):
        #self.entries is a list of lists extracted from the lines of the library file and splitted at the ;
        self.entries = []
        with open("library_management.txt","r") as infile:
            for line in infile:
                buffer = line.split(";")
                self.entries.append(buffer)
         

    def update(self):
        files = []
        for file in os.listdir("library"):
            if file.endswith(".sfg"):
                files.append(file)
        

        newfiles = []
        for file in files:
            for checkfile in self.entries:
                if file == checkfile[1].strip():
                    print("File already in library!")
            else:
                    newfiles.append(file)
        
        
        with open ("library_management.txt","a") as outfile:
            counter = len(self.entries)
            for file in newfiles:
                counter += 1
                sfg = FileFetcher(file).sfg
                specrange = sfg.yield_spectral_range()
                outfile.write(str(counter)+";"+file+";"+str(specrange[0])+";"+str(specrange[1])+";"+str(specrange[2])+"\n")
 
#tools for direct spectra management
class SFG_Spectrum:

    def __init__(self,wavenumbers,intensity,ir_intensity,vis_intensity,systematic_name):
        self.wavenumbers = wavenumbers
        self.raw_intensity = intensity
        self.vis_intensity = vis_intensity
        self.ir_intensity = ir_intensity
        self.name = systematic_name
        self.normalized_intensity = self.raw_intensity/(self.vis_intensity*self.ir_intensity)
     

    #normalize an given array to its maximum, typically the normalized or raw intensity
    def normalize_to_highest(self,intensity="default"):
        if intensity == "default":
            intensity = self.normalized_intensity
        norm_factor = np.max(intensity)
        return (intensity/norm_factor)

    #yield a defined list of peaks separated from each other by minimum the threshold value in wavenumber
    def yield_peaklist(self,intensity="default",num=6,threshold=25):
         pair_get=[]
         out=[]
         num=num
         if intensity == "default":
             intensity = self.normalized_intensity

         for i in range(len(self.wavenumbers)):
            pair_get.append([intensity[i],self.wavenumbers[i]])
        
         while len(out) < (num):
        
        
              if len(out) != 0:
            
                 for i in range(len(out)):
                    k=max(pair_get)
                    
                    if np.abs(out[i][1]-k[1]) < threshold:
                       pair_get.remove(k) 
                       break
                    if i==(len(out)-1):
                        out.append(k)
                        pair_get.remove(k)
           
              
              else:
                   k=max(pair_get)    
                   out.append(k)
            
           
         return out

    #returns a list containing maximum and minimum wavenumer and the number of data points
    def yield_spectral_range(self):
        return [min(self.wavenumbers),max(self.wavenumbers),len(self.wavenumbers)]

    def yield_increment(self):
        """Calculates stepsize and wavenumbers where the stepsize is changed"""
        borders = []
        stepsize = []
        current = self.wavenumbers[0]
        currentstep = abs(current-self.wavenumbers[1])
        borders.append(current)

        for wavenumber in self.wavenumbers[1:]:
            s = abs(wavenumber - current )
            if s != currentstep:
                stepsize.append(currentstep)
                borders.append(current)
                currentstep = s
                current = wavenumber
            else:
                current = wavenumber
        borders.append(current)
        stepsize.append(s)
        increment = [borders,stepsize]
        return increment
    
    
    def smooth(self,points=9,order=5):
        y = savgol_filter(self.normalized_intensity,points,order)
        self.normalized_intensity = y
class Systematic_Name:

    def __init__(self,namestring):
        #load the allowed surfactans and sensitizers from files
        self.Surfactants = {}
        self.Sensitizers = {}

        with open ("../name_info/Surfactants.txt","r") as infile:
            for line in infile:
                collect = line.split(":")
                self.Surfactants[collect[0]] = collect[1].strip()

        with open ("../name_info/Sensitizers.txt","r") as infile:
            for line in infile:
                collect = line.split(":")
                self.Sensitizers[collect[0]] = collect[1].strip()

        #now the actual processing begins
        self.full_name = namestring
        processing_list =self.full_name[:-4]#removes .sfg
        self.processing_list = processing_list.split("_") 
        self.date = self.processing_list [0]
      
        #setting defaults
        self.surfactant = "unknown"
        self.surfactant_spread_volume = "unknown"
        self.sensitizer = "-"
        self.sensitizer_spread_volume = "-"
        self.photolysis ="none"
        self.sample_number = 1
        self.measurement = "#1"         #das ist nicht optimal, sample_number als int und measurement als string zu haben
        self.comment = "none"
        if len(self.processing_list) == 2:
            self.check_boknis()
        #traversing the processing list 
        for i in self.processing_list[1:]:
           

            if i in self.Surfactants:
                self.surfactant = self.Surfactants[i]

            elif i in self.Sensitizers:
                self.sensitizer = self.Sensitizers[i]

            #the following line could result in a problem with the date (processing_list[0])
            elif self.is_number(i) == True:
                if self.surfactant_spread_volume == "unknown":
                    self.surfactant_spread_volume = i
                else:
                    sensitizer_spread_volume = i
            #pH value info will be comment
            elif "p" in i and "H" in i:
                self.comment = i

            #photolysis handling
            elif "p" in i:
                 if self.photolysis =="none":
                     time=i.strip("p")
                     if self.is_number(time)==True:
                        self.photolysis= time+" min."    
                     elif "h" in time:
                         if self.is_number(time.strip("h"))==True:                   
                           self.photolysis= str(int(time.strip("h"))*60)+" min."
                         else:
                             self.comment=i
                 else:
                    self.comment=i

            #check for sample number
            elif "x" in i:
               if self.is_number(i.strip("x"))==True:
                 self.sample_number=i.strip("x")
               else:
                   self.comment=i

            #check for measurement number
            elif "#" in i:
              self.measurement= i  #hier auch ein String! ggf. aendern
            
            else:
                if self.comment == "none":
                    self.comment = i
                else:
                    self.comment += i
        
        #postprocessing
       
        if self.surfactant == "unknown" and self.sensitizer != "-":
            self.surfactant = self.sensitizer
            self.sensitizer = "-"

        if self.sensitizer=="-":
            self.sensitizer_spread_volume ="-"
       
    def is_number(self,s):
        try:
            float(s)
            return True
        except ValueError:
            return False 
    
    def yield_info(self):
        return(self.full_name,self.date,self.surfactant,self.surfactant_spread_volume,
               self.sensitizer,self.sensitizer_spread_volume,self.photolysis,
               self.sample_number,self.measurement,self.comment)

    def check_boknis(self):
 
        if self.is_number(self.processing_list[0]) == True and self.is_number(self.processing_list[1]) == True:
            self.comment ="BoknisEckSample"
            self.sample_number = self.processing_list[1]
            self.surfactant="Nat. surface sample"
    def date_split(self):
        year = self.date[:3]
        month = self.date[4:5]
        day = self.date[6:7]
        return(year,month,day)
        
class FileFetcher:
    """The file fetcher class interconnects different subdirectories of the folder and simplifies
    the handling of filenames and filepaths. It changes the working directory to the file storage directory,
    usually the library, creates a Data_Collector obect and returns a SFG object"""

    def __init__(self,filename,destination="library"):
        self.filename = filename
        initial_wd = os.getcwd()
        os.chdir(destination)

        self.collector = Data_Collector(filename)
        self.sfg = self.collector.yield_SFG_Spectrum()
        os.chdir(initial_wd)

#extract spectral data from .sfg file
class Data_Collector:
    def __init__(self,filename):
        #the filename is ONLY the filename, not containing any directory information. 
        self.file = filename


        #now start extraction
        data_collect=[]
            
        with open(self.file, "r") as infile:

            readCSV=csv.reader(infile, delimiter="\t")    
            for row in readCSV:
                     data_collect.append(row)
                    
            data_package=[0]*len(data_collect[0])
            for i in range(len(data_collect[0])):
                #mind this complex list comprehension
                data_package[i]=[j[i] for j in data_collect] 
            #remove useles column
            del data_package[2] 

            #convert strings to float and generate numpy array
            convert=[]
            for i in data_package:
                q=[float(j) for j in i]
                convert.append(q)
            convert=[np.array(i) for i in convert]
            self.data_package = convert

    def yield_SFG_Spectrum(self):
        #this function will need exception handling. It retuns an SFG spectrum object
        sysname=Systematic_Name(self.file)
        data=self.data_package
        sfg=SFG_Spectrum(data[0],data[1],data[3],data[2],sysname)
        return sfg

class Finder:
    """Powerfull class generating sfg objects of ALL files in library and traversing them for match
    criteria. Contains a list of SFG objects whitch mach the criteria. Later features to extract 
    information from libary management file and day information file will be added"""

    def __init__(self):      
        self.database = []
        for file in os.listdir("library"):
            if file.endswith(".sfg"):
                sfg = FileFetcher(file).sfg
                self.database.append(sfg)

    #Basic search funcions
    def date_based(self,dateflag,subset="default"):
        if subset == "default":
            subset = self.database
        matches = []
        dateflag = dateflag.strip("d")
        for spectrum in subset:
            if spectrum.name.date == dateflag:
                matches.append(spectrum)
        return matches

    def sample_based(self,sampleflag,subset="default"):
        if subset == "default":
            subset = self.database
        matches = []
        for spectrum in subset:
            if spectrum.name.sample_number == sampleflag or str(spectrum.name.sample_number) == str(sampleflag):
                matches.append(spectrum)
        return matches

    def surfactant_based(self,surfactantflag,subset="default"):
        if subset == "default":
            subset = self.database
        matches = []
        for spectrum in subset:
            if spectrum.name.surfactant == surfactantflag:
                matches.append(spectrum)
        return matches

    def sensitizer_based(self,sensitizerflag,subset="default"):
        if subset == "default":
            subset = self.database
        matches=[]
        for spectrum in subset:
            if spectrum.name.surfactant == sensitizerflag:
                matches.append(spectrum)
        return matches
    
    def photo_based(self,subset="default",photolyzed=True):
        if photolyzed == True:
            if subset == "default":
                subset = self.database
            matches = []
            for spectrum in subset:
                if spectrum.name.photolysis != "none":
                    matches.append(spectrum)
        else:
            if subset == "default":
                subset = self.database
            matches = []
            for spectrum in subset:
                if spectrum.name.photolysis == "none":
                    matches.append(spectrum)            
        return matches
    def measurement_based(self,measurementflag,subset="default"):
        if subset == "default":
            subset = self.database
        measures = measurementflag.split(",")
        matches = [ i for i in subset if i.name.measurement in measures]
        return matches
    def comment_based(self,option="BoknisEckSample",subset="default"):
        if subset == "default":
            subset = self.database
        matches = []
        for spectrum in subset:
            if spectrum.name.comment == option:
                matches.append(spectrum)
        return matches
    def sur_volume_based(self,amount, subset= "default"):
        if subset == "default":
            subset = self.database
        matches = [i for i in subset if i.name.surfactant_spread_volume == amount]
        return matches
        
    def unclassified(self):
            subset = self.database
            matches = [i for i in subset if i.name.surfactant == "unknown"]
            return matches

   
class Interpreter:

    """Class providing the functionality for the command line of the plotting routine"""

    def __init__(self,command):
       
       commandlist = command.split(" ")
       if len(commandlist) == 1:
           if commandlist[0] == "ud":
               self.update()
           elif commandlist[0] == "x":
               sys.exit()

       elif len(command) < 3:
           print("Wrong number of parameters! Try again")
       else:
           self.type = commandlist[0]
           self.flags = commandlist [1]
           self.options = commandlist[2]

           if self.type == "plot" and self.flags == "f":
                files = self.options.split(",")
                self.plothandler(files)

           elif self.type == "plot" and self.flags =="d,su,sa":
               self.sample_surfactant_date(self.options)

           elif self.type == "plot" and self.flags ==  "bo":
               self.plot_boknis()

           elif self.type =="list":
               self.list_()
           

    def plothandler(self,filelist):

        sfg_objects = []
        for file in filelist:
            sfg = FileFetcher(file).sfg
            sfg_objects.append(sfg)
        P = Plotter(sfg_objects)
        P.simple_plot()

    def sample_surfactant_date(self,options):
        options = options.split(",")
        if len(options) != 3:
            print("sa_su_da calling. Invalid number of flags!")
        else:
            dateflag = options[0]
            surflag = self.retranslate_name(options[1])
            sampleflag = int(options[2])

            f = Finder()
            f1 = f.date_based(dateflag)
            f2 = f.surfactant_based(surflag,f1)
            f3 = f.sample_based(sampleflag,f2)

            Plotter(f3).simple_plot()

    def plot_boknis(self):
        if self.options == "none":
            f = Finder()
            sfg = f.comment_based()
            Plotter(sfg).simple_plot()

        elif "s" in self.options:
            number = int(self.options[1])
            f = Finder()
            sfg = f.comment_based()
            sfg = f.sample_based(number,sfg)
            Plotter(sfg).simple_plot()

        else:
            print("Not yet implemented")
            print(self.options, len(self.options))

    def update(self):
        Importer()

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

    def list_(self):
        f = Finder()

        if self.flags == "su":
            surfoptions = self.options.split(",") #if comma-separated multiplit sensitizers are given
            subsets = []
            for  i in surfoptions:
                i= self.retranslate_name(i)
                f1 = f.surfactant_based(i)
                for j in f1:
                    subsets.append(j)
            self.answerset = subsets
        
        elif self.flags =="d":

            subsets = []
            dateoptions = self.options.split(",")
            for i in dateoptions:
                subset = f.date_based(self.options)
                for j in subset:
                    subsets.append(j)
            self.answerset = subsets

 
#tools for plotting
class Plotter:

    def __init__(self,speclist,raw=False,title="Default"):
        self.speclist = speclist
        self.raw = raw
        self.title = title
        
    def simple_plot(self):
        
        fig = plt.figure()
        ax = plt.subplot(111)
        for spectrum in self.speclist:
            wl = spectrum.wavenumbers
            if self.raw == False:
                intensity = spectrum.normalized_intensity
            else:
                intensity = spectrum.raw_intensity
            ax.plot(wl,intensity,label=spectrum.name.full_name,linestyle='--',markersize=4,marker="o")
        
        
        
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])      
        ax.grid()
        ax.set_title(self.title)
        ax.set_xlabel("Wavenumber/ $cm^{-1}$")
        ax.set_ylabel("Intensity/ a.u.")
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        plt.show()
        
    def raw_plot(self):
        for spectrum in self.speclist:
            wl = spectrum.wavenumbers
            intensity = spectrum.raw_intensity
            plt.plot(wl,intensity,label=spectrum.name.full_name,linestyle='--',markersize=4,marker="o")
    
        plt.grid(True)
        plt.title(self.title)
        plt.xlabel("Wavenumber/ $cm^{-1}$")
        plt.ylabel("Raw intensity/ a.u.")
        if len(self.speclist < 6):
            plt.legend(loc="upper right")
        else:
            plt.legend(bbox_to_anchor(1,0.5),loc='center left', numpoints=1)
        plt.show() 
    def raw_plot_plus_ir(self):
        for spectrum in self.speclist:
            wl = spectrum.wavenumbers
            intensity = spectrum.raw_intensity
            ir = spectrum.ir_intensity
            plt.plot(wl,intensity,label=spectrum.name.full_name,linestyle='--',markersize=4,marker="o")
            plt.plot(wl,ir,label=spectrum.name.full_name+" IR_intensity",linestyle='--')
        plt.grid(True)
        plt.title(self.title)
        plt.xlabel("Wavenumber/ $cm^{-1}$")
        plt.ylabel("Raw intensity/ a.u.")
        plt.legend(loc="upper right")
        plt.show() 
    def custom_plot(self):
        
        fig = plt.figure()
        ax = plt.subplot(111)
        for spectrum in self.speclist:
            wl = spectrum.wavenumbers
            if self.raw == False:
                intensity = spectrum.normalized_intensity
            else:
                intensity = spectrum.raw_intensity
            stitle = str(spectrum.name.sample_number)+" "+spectrum.name.comment+" "+str(spectrum.name.measurement)
            ax.plot(wl,intensity,label=stitle,linestyle='--',markersize=4,marker="o")
        
        
        
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])      
        ax.grid()
        ax.set_title(self.title)
        ax.set_xlabel("Wavenumber/ $cm^{-1}$")
        ax.set_ylabel("Intensity/ a.u.")
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        plt.savefig("analysis_out/"+self.title+".pdf")
        plt.close()
        with open("analysis_out/"+self.title+".log","w") as outfile:

            outfile.write("Logfile of Spectal anaylsis routine")
            outfile.write("List of included specta/n")
            for spectrum in self.speclist:
                outfile.write(spectrum.name.full_name+"/n")

class Analyzer:
    """This class takes, what a surprise, a list of SFG spectra as constructor argument. Its purpose
    is to perform analytical tasks with the spectral data, e.g. compare peaks, integral, datapoints
    and will be extendet to handle statistics in the future"""
    def __init__(self,speclist):
        
        self.speclist = speclist
    def list_peaks(self,number):
        intensities = []
        wavenumbers = []
        
        for i in self.speclist:
            j = i.yield_peaklist(num=number)
            for k in j:
                intensities.append(k[0])
                wavenumbers.append(k[1])
        #plt.scatter(wavenumbers,intensities)
        plt.hist(wavenumbers,rwidth=0.02,normed=True)
        plt.show()
#test code section
"""dc = Data_Collector("20170811_DPPC_9.81.sfg")
dc2 = Data_Collector("20170811_SA_10_x1_#1.sfg")

sg1 = dc.yield_SFG_Spectrum()
sg2= dc2.yield_SFG_Spectrum()


p = Plotter([sg1,sg2])
p.simple_plot()
def testfunc():
    d = os.getcwd()
    print(d)
    a = FileFetcher("20170811_DPPC_9.81.sfg")
    s=a.sfg
    d = os.getcwd()
    p=Plotter([s])
    p.simple_plot()
    print(s.yield_increment())
    print(d)
I = Importer()
L = Library_Manager()
L.update()

f = Finder()
f1=f.date_based("20170811")
f2=f.surfactant_based("Stearic Acid",f1)

f3=f.sample_based(1,f2)

p = Plotter(f3)
p.simple_plot()"""

#testfunc()
