from ipy_interpreter import Ipy_Interpreter
from Classes import Plotter

"""A supporting module with useful functions for all the data management and analysis in the SFG module"""

def simple_analysis():

    Substances = {}
    i = Ipy_Interpreter()

    with open ("name_info/Surfactants.txt","r") as infile:
                for line in infile:
                    collect = line.split(":")
                    Substances[collect[0]] = collect[1].strip()
        
    with open ("name_info/Sensitizers.txt","r") as infile:
                for line in infile:
                    collect = line.split(":")
                    Substances[collect[0]] = collect[1].strip()

    for key in Substances:
        if key != "DPPC":
            dates = []
            #erster Schritt: alle Daten des Surfactants holen
            i.get("su "+key)

            #zweiter Schritt: Datum extrahieren
            for j in i.subset:
                dates.append(j.name.date)

            #an jedem Tag alle Probennummern extrahieren
            for d in dates:
                i.get("su "+key)
                daysamples = []
                i.refine("d "+d)
                for spectrum in i.subset:
                    daysamples.append(spectrum.name.sample_number)

             #die einzelnen Samples zusammen plotten   
                for sample in daysamples:
                    plotllist = [q for q in i.subset if q.name.sample_number == sample]
                    if len(plotllist) > 1:
                        plot_title = key+" "+str(d)+" "+str(sample)
                        P = Plotter(plotllist,title=plot_title)
                        P.custom_plot()
                    


                i.clear()

             


       


        #subset wieder säubern
        i.clear()



#Test code section

simple_analysis()

