# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.

"""

from Classes import SessionControlManager, baseline_demo
import matplotlib.pyplot as plt


class StationTexGenerator:
    
    def __init__(self, station):
        
        self.station = station
        self.sfg_labels = []
        self.lt_labels = []
        self.make_sfg_plots()
        self.make_lt_plots()
    
    def make_namestring(self):
        s = f'Station {self.station.station_number} on {self.station.date.strftime("%m/%d")} of type {self.station.type}'
        return s

    def make_infostring(self):
        
        
        s = f"""\nIt contains:\n\\begin{{itemize}}\n"""
        s += f"""\\item {len(self.station.sfg_spectra)} SFG spectra\n"""
        s += f"""\\item {len(self.station.lt_isotherms)} compression isotherms\n""" 
        s += f"""\\item {len(self.station.tensions)} surface tension measurements\n"""
        s += f"""\\end{{itemize}}\n"""
                       
        return s

    def make_sfg_infostring(self):
        
        pass
    
    def make_tension_table(self):
        
        
        s = """\\begin{table}\n\\centering\n\\begin{tabular}{cc}Sample  & surface tension /\\\\&\\(\si{\\milli\\newton\\per\\meter}\\)\\\\\n\\hline\n\\hline"""
            
        
#        for tens in self.station.tensions:
#            if tens is not None:
#                
#                texname = tens[0].split("_")
#                texname = " ".join(texname)
#                
#                s += f'{texname} & {tens[1]} \\\\\n'
        try:
            #s += f'\\hline\n'
            s += f'\n'
            s += f'total average & {self.station.stats["tension_average"][0]:.2f} \\\\\n'
            s += f'SML samples & {self.station.stats["tension_sml"][0]:.2f} \\\\\n'
            s += f'deep water & {self.station.stats["tension_deep"][0]:.2f} \\\\\n'
            s += f'plate average & {self.station.stats["tension_plate"][0]:.2f} \\\\\n'
            s += f'screen average & {self.station.stats["tension_screen"][0]:.2f} \\\\\n'
            
        except TypeError:
            pass
        
    
        s += "\\hline\n\\end{tabular}\n\\end{table}"
        
        return s
    
    def make_sfg_table(self):
        s = """\\begin{table}\n\\centering\n\\begin{tabular}{cc}Samples & CH integral /\\\\&arb. units\\\\\n\\hline\n\\hline"""
        
        try:
            s += f'\n'
            s += f'average & {self.station.stats["ch_average"][0]:.7f} \\\\\n'
            s += f'deep water & {self.station.stats["ch_deep"][0]:.7f} \\\\\n'
            s += f'SML sample & {self.station.stats["ch_sml"][0]:.7f} \\\\\n'
            s += f'plate average & {self.station.stats["ch_plate"][0]:.7f} \\\\\n'
            s += f'screen average & {self.station.stats["ch_screen"][0]:.7f} \\\\\n'
        
        except TypeError:
            pass
            
            
        s += "\\hline\n\\end{tabular}\n\\end{table}"
        
        return s
    
    def make_lt_table(self):
        
         s = """\\begin{table}\n\\centering\n\\begin{tabular}{cc}Sample  & max. surface pressure /\\\\&\\(\si{\\milli\\newton\\per\\meter}\\)\\\\\n\\hline\n\\hline"""
         
         try:
            s += f'\n'
            s += f'average & {self.station.stats["pressure_average"][0]:.2} \\\\\n'
            s += f'deep water & {self.station.stats["pressure_deep"][0]:.2f} \\\\\n'
            s += f'SML sample & {self.station.stats["pressure_sml"][0]:.2f} \\\\\n'
            s += f'plate average & {self.station.stats["pressure_plate"][0]:.2f} \\\\\n'
            s += f'screen average & {self.station.stats["pressure_screen"][0]:.2f} \\\\\n'
        
         except TypeError:
            pass
         
         
         s += "\\hline\n\\end{tabular}\n\\end{table}"
        
         return s

    def set_plot(self, path):

        s = f"""\\begin{{figure}}\n\\centering\n
            \\includegraphics[width=0.7\\linewidth]{{{path}}}\n
            \\end{{figure}}\n"""
        return s
        

    def generate_tex(self):
        
        s = f'\\subsubsection* {{{self.make_namestring()}}}\n'
        
        s += self.make_infostring()
        
        if len(self.station.sfg_spectra) > 0:
            
            s += f'Overview about the SFG measurements of the sample: \n'
            s += self.make_sfg_table()
            
            s += f'The following SFG spectra belong to the station:\n'
            for spec in self.sfg_labels:
                s += self.set_plot(spec)

        if len(self.station.tensions) > 0:
            s += f'Overview about the surface tension measurements of the sample: \n'
            s += self.make_tension_table()
        
        if len(self.station.lt_isotherms) > 0:
            s += f'Overview about the compression isotherm measurements of the sample: \n'
            s += self.make_lt_table()
            s += f'\nThe following compression isotherms belong to the station:\n'
            for spec in self.lt_labels:
                s += self.set_plot(spec)
             
        s += "\n"
        return s
    
    def make_sfg_plots(self):
        label = self.station.station_hash+"_sfg"+".png"
        plt.xlabel("Wavenumber/ cm$^{-1}$")
        plt.ylabel("Norm. SFG intensity/ arb. u.")
        for spec in self.station.sfg_spectra:
            plt.plot(spec.wavenumbers, spec.normalized_intensity, label=spec.name.full_name[:-4])
        
        plt.legend()
        plt.savefig(label)

        label = ("fig/" + label)
        self.sfg_labels.append(label)
        plt.cla()
    
    def make_lt_plots(self):
        label = self.station.station_hash + "_lt" + ".png"
        plt.xlabel("area/ cm$^{-1}$")
        plt.ylabel("Surface pressure/ mN $\cdot m^{-1}$")
        for iso in self.station.lt_isotherms:
            plt.plot(iso.area, iso.pressure, label=iso.name)
            
        plt.legend()
        plt.savefig(label)
        label = ("fig/" + label)
        self.lt_labels.append(label)
        plt.cla()
        


S = SessionControlManager("sfg.db", "test")
S.setup_for_gasex()

stats = []

for s in S.stations.values():
    stats.append(s)

stats.sort()

T = StationTexGenerator(stats[0])

s = T.generate_tex()
with open("out.tex", "w") as outfile:
    outfile.write(s)



