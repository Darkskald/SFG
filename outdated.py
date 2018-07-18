class Planer:
    """A class to perform task management and plan the research tasks for a period of measurements"""

    def __init__(self):

        self.tasks = []
        self.count = 0
        self.refresh()
        print("Planer initialized. The current task count is " + str(self.count))
        self.show_tasks()

    def show_tasks(self):

        index = 1
        for i in self.tasks:
            print(str(index) + " : " + i)
            index += 1

    def add_task(self, string):

        with open("tasks", "a") as infile:
            infile.write(string + "\n")
        self.refresh()

    def refresh(self):
        self.tasks = []
        with open("tasks", "r") as infile:
            for line in infile:
                self.tasks.append(line)
                self.count = len(self.tasks)

    def done(self, number):

        print(str(self.tasks[number - 1]) + "was removed!")
        del self.tasks[number - 1]
        self.update_file()

    def update_file(self):

        with open("tasks", "w") as outfile:
            for task in self.tasks:
                outfile.write(task)


class TexInterface:
    def __init__(self, name, sfg_list):
        self.name = name
        self.blocks = []
        self.out_dir = "tex_out/"
        self.fig_dir = "fig/"
        self.sfg_list = sfg_list

        # section control
        self.create_header()

    def add_figure(self, filename, label):
        figstring = "\\begin{figure}[h!]\n\centering\n\includegraphics[scale=0.4]"
        figstring += "{" + self.fig_dir + filename + "}\n"
        figstring += "\caption{" + label + ".}\n\end{figure}\n"
        self.blocks.append(figstring)

    def add_tabular(self, spectrum):
        corr_name = " ".join(spectrum.name.full_name.split("_"))
        corr_name = corr_name.replace("#", "no")
        tabular_string = "\\begin{table}[h!]\n"
        tabular_string += "\caption{" + "Peaklist of " + corr_name + "}" + "\n"
        tabular_string += "\\begin{center}\n"

        tabular_string += spectrum.drop_tex_peaktable()

        tabular_string += "\\end{center}\n"
        tabular_string += "\\end{table}\n"
        print("DEBUUUUUUUUUG")
        self.blocks.append(tabular_string)

    def load_template(self, name):
        pass

    def create_header(self):
        with open("texhead", "r") as infile:
            string = infile.read()
            self.blocks.append(string)

    def write_file(self):
        self.blocks.append("\n" + "\end{document}")
        outstring = "\n".join(self.blocks)
        with open(self.out_dir + self.name + ".tex", "w") as outfile:
            outfile.write(outstring)

    def compile(self):
        pass

    def add_section(self, name):
        self.blocks.append("\section{" + name + "}" + "\n")

    def traverse_sfg_list(self):
        for spectrum in self.sfg_list:
            corr_name = " ".join(spectrum.name.full_name.split("_"))
            corr_name = corr_name.replace("#", "NR. ")
            self.add_section(corr_name)
            P = Plotter([spectrum])
            s = spectrum.name.full_name[:-4].replace(".", "x")
            s = s.replace("#", "no")
            P.title = s
            P.simple_plot(mode="save")
            self.add_figure(s + ".pdf", "Simple plot with normalized intensity")
            self.add_tabular(spectrum)

    def join_speclist(self):
        P = Plotter(self.sfg_list)
        P.title = self.name
        P.simple_plot(mode="save")
        self.add_section(self.name)
        self.add_figure(self.name + ".pdf", "Joint plot of a list of spectra")


class Plotter:
    """Simple class providing an interface between the matplotlib API and the SFG objects. Depending on the specific method
    and kwargs, a variety of different methods is accessible"""
    def __init__(self, speclist, raw=False, title="Default"):
        self.speclist = speclist
        self.raw = raw
        self.title = title
        self.save_dir = "tex_out/fig"

    def simple_plot(self, mode="show"):

        fig = plt.figure()
        ax = plt.subplot(111)
        for spectrum in self.speclist:
            wl = spectrum.wavenumbers
            if self.raw == False:
                intensity = spectrum.normalized_intensity
            else:
                intensity = spectrum.raw_intensity
            ax.plot(wl, intensity, label=spectrum.name.full_name, linestyle='--', markersize=4, marker="o")

        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
        ax.grid()
        ax.set_title(self.title)
        ax.set_xlabel("Wavenumber/ $cm^{-1}$")
        ax.set_ylabel("Intensity/ a.u.")
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        if mode == "show":
            plt.show()
        elif mode =="save":
            plt.savefig(self.save_dir+"/"+self.title+".pdf")

    def raw_plot(self, mode="show"):
        for spectrum in self.speclist:
            wl = spectrum.wavenumbers
            intensity = spectrum.raw_intensity
            plt.plot(wl, intensity, label=spectrum.name.full_name, linestyle='--', markersize=4, marker="o")

        plt.grid(True)
        plt.title(self.title)
        plt.xlabel("Wavenumber/ $cm^{-1}$")
        plt.ylabel("Raw intensity/ a.u.")
        if len(self.speclist) < 6:
            plt.legend(loc="upper right")
        else:
            plt.legend(bbox_to_anchor=(1, 0.5), loc='center left', numpoints=1)
        if mode == "show":
            plt.show()
        elif mode =="save":
            plt.savefig(self.save_dir+"/"+self.title+".pdf")

    def raw_plot_plus_ir(self, mode="show"):
        for spectrum in self.speclist:
            wl = spectrum.wavenumbers
            intensity = spectrum.raw_intensity
            ir = spectrum.ir_intensity
            plt.plot(wl, intensity, label=spectrum.name.full_name, linestyle='dotted', markersize=4, marker="o")
            plt.plot(wl, ir, label=spectrum.name.full_name + " IR_intensity", linestyle='--')
        plt.grid(True)
        plt.title(self.title)
        plt.xlabel("Wavenumber/ $cm^{-1}$")
        plt.ylabel("Raw intensity/ a.u.")
        plt.legend(loc="upper right")
        if mode == "show":
            plt.show()
        elif mode =="save":
            plt.savefig(self.save_dir+"/"+self.title+".pdf")

    def custom_plot(self):

        fig = plt.figure()
        ax = plt.subplot(111)
        for spectrum in self.speclist:
            wl = spectrum.wavenumbers
            if self.raw == False:
                intensity = spectrum.normalized_intensity
            else:
                intensity = spectrum.raw_intensity
            photo = ""

            if spectrum.name.photolysis != "none":
                photo = spectrum.name.photolysis

            stitle = str(spectrum.name.sample_number) + " " + spectrum.name.comment + " " + str(
                spectrum.name.measurement) + " " + photo
            ax.plot(wl, intensity, label=stitle, linestyle='dotted', markersize=2, marker="o")

        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
        ax.grid()
        ax.set_title(self.title)
        ax.set_xlabel("Wavenumber/ $cm^{-1}$")
        ax.set_ylabel("Intensity/ a.u.")
        lgd = ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        plt.savefig("analysis_out/" + self.title + ".pdf", bbox_extra_artists=(lgd,), bbox_inches='tight')
        plt.close()
        with open("analysis_out/" + self.title + ".log", "w") as outfile:

            outfile.write("Logfile of Spectral analysis routine")
            outfile.write("List of included spectra\n")
            for spectrum in self.speclist:
                outfile.write(spectrum.name.full_name + "\n")

    def stack_plot(self, mode="show"):
        spacer = 0
        fig = plt.figure()
        ax = plt.subplot(111)
        for spectrum in self.speclist:
            fig = plt.figure()
            ax = plt.subplot(111)
            for spectrum in self.speclist:
                wl = spectrum.wavenumbers
                intensity = spectrum.normalize_to_highest()+spacer
                ax.plot(wl, intensity, label=spectrum.name.full_name, linestyle='--', markersize=4, marker="o")
                spacer += 0.5

            box = ax.get_position()
            ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
            ax.grid()
            ax.set_title(self.title)
            ax.set_xlabel("Wavenumber/ $cm^{-1}$")
            ax.set_ylabel("Intensity/ a.u.")
            ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
            ax.axes.get_yaxis().set_ticks([])
            self.log_plot("stack")
            if mode == "show":
                plt.show()
            elif mode == "save":
                plt.savefig(self.save_dir + "/" + self.title + ".pdf")

    def log_plot(self, mode):
        # Part 1: Generate logfile
        with open(self.title+"_"+mode, "w") as outfile:
            outfile.write("Log file for plot " + self.title + "\n")
            outfile.write("The mode of the plot is: " + mode + "\n")
            outfile.write("Plot contains the following spectra: \n")
            for spectrum in self.speclist:
                outfile.write(spectrum.name.full_name + "\n")

    def bar_peaks(self, number=4):

        A = Analyzer(self.speclist)
        dic = A.count_peak_abundance(number)


        for key in dic:
            plt.bar(key, dic[key], width=3.)

        plt.xlabel("wavenumber")
        plt.ylabel("rel. peak abundance")
        plt.grid(True)
        plt.minorticks_on()
        plt.show()

    def marked_peaks(self, mode="show"):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax2 = ax.twiny()

        for spectrum in self.speclist:
            wl = spectrum.wavenumbers
            if self.raw == False:
                intensity = spectrum.normalized_intensity
            else:
                intensity = spectrum.raw_intensity
            ax.plot(wl, intensity, label=spectrum.name.full_name, linestyle='--', markersize=4, marker="o")

            peakticks = []
            for peak in spectrum.detailed_analysis(threshold=1.2):
                ax2.axvline(peak[0], color="red")
                ax2.axvline(peak[1], ls="dashed", color="blue")
                ax2.axvline(peak[2], ls="dashed", color="blue")
                peakticks.append(peak[0])

        ax.grid()
        ax.set_xlabel("Wavenumber/ $cm^{-1}$")
        ax.set_ylabel("Intensity/ a.u.")
        ax.legend()


        ax2.set_xticks(peakticks)
        ax2.set_xticklabels(peakticks, rotation=45, color='blue')
        ax2.set_xlim(ax.get_xlim())

        if mode == "show":
            plt.show()
        elif mode == "save":
            plt.savefig(self.save_dir + "/" + self.title + ".pdf")

    def relation_coverage_peaks(self):

        A = Analyzer(self.speclist)
        fig = plt.figure()
        ax = plt.subplot(111)
        for entry in A.intensity_vs_surfcoverage():
            ax.scatter(entry[0], entry[1])
        ax.set_xlabel("area per molecule/ A$^{2}$")
        ax.set_ylabel("maximum intensity")
        plt.show()

    def relation_3d(self):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        wavelengths = []
        intensities = []
        surface_concentrations = []

        for spectrum in self.speclist:
            for peak in spectrum.peaks:
                wavelengths.append(peak[0])
                intensities.append(peak[3])
                surface_concentrations.append(spectrum.calc_area_per_molecule())

        ax.scatter(wavelengths, intensities, surface_concentrations)
        ax.set_xlabel("wavelength")
        ax.set_ylabel("intensity")
        ax.set_zlabel("area per molecule")
        plt.show()