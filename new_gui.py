from plotgui import Ui_MainWindow

import time
import sys
import matplotlib
# Make sure that we are using QT5
matplotlib.use('Qt5Agg')
from PyQt5 import QtCore, QtWidgets, QtGui

from numpy import arange, sin, pi
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
from matplotlib import rcParams
rcParams['mathtext.default'] = 'regular'




class MyMplCanvas(FigureCanvas):
    """Ultimately, this is a QWidget (as well as a FigureCanvasAgg, etc.)."""

    def __init__(self, parent=None, width=5, height=4, dpi=100):

        fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = fig.add_subplot(111)

        FigureCanvas.__init__(self, fig)
        self.setParent(parent)

        FigureCanvas.setSizePolicy(self,
                                   QtWidgets.QSizePolicy.Expanding,
                                   QtWidgets.QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)
        self.init_axes()
        self.compute_initial_figure()

    def init_axes(self,title="default"):
        self.axes.grid()
        self.axes.set_title(title)
        self.axes.set_xlabel("Wavenumber/ $cm^{-1}$")
        self.axes.set_ylabel("Intensity/ arb.u.")

    def compute_initial_figure(self):
        pass


class MyStaticMplCanvas(MyMplCanvas):
    """Simple canvas with a sine plot."""

    def compute_initial_figure(self):
        pass

    def plot_sfg(self, sfg_list,title="default", flag="none"):
        self.axes.cla()
        self.init_axes(title)
        for spectrum in sfg_list:
            if flag == "IR":
                self.axes.plot(spectrum.wavenumbers, spectrum.ir_intensity, label=(spectrum.name.full_name + "_IR"), linestyle='dotted', markersize=2, marker="o")
                self.axes.plot(spectrum.wavenumbers, spectrum.raw_intensity, label=spectrum.name.full_name, linestyle='dotted', markersize=4, marker="o")

            elif flag == "Vis":
                self.axes.plot(spectrum.wavenumbers, spectrum.vis_intensity, label=spectrum.name.full_name, linestyle='dotted', markersize=2, marker="o")
                self.axes.plot(spectrum.wavenumbers, spectrum.raw_intensity, label=(spectrum.name.full_name + "_vis"), linestyle='dotted', markersize=4, marker="o")
            elif flag == "b":
                self.axes.plot(spectrum.wavenumbers, spectrum.ir_intensity, label=(spectrum.name.full_name + "_IR"), linestyle='dotted', markersize=2, marker="o")
                self.axes.plot(spectrum.wavenumbers, spectrum.raw_intensity, label=spectrum.name.full_name, linestyle='dotted', markersize=4, marker="o")
                self.axes.plot(spectrum.wavenumbers, spectrum.vis_intensity, label=(spectrum.name.full_name + "_vis"), linestyle='dotted', markersize=2, marker="o")
            elif flag == "r":
                self.axes.plot(spectrum.wavenumbers, spectrum.raw_intensity, label=spectrum.name.full_name, linestyle='dotted', markersize=4, marker="o")
            else:
                self.plot_dataset(spectrum.wavenumbers, spectrum.normalized_intensity, spectrum.name.full_name, "")

        self.axes.legend()
        self.draw()

    def plot_dataset(self, x_dataset, y_dataset, fullname, label):
        self.axes.plot(x_dataset, y_dataset, label=(fullname + label),
                       linestyle='dotted', markersize=2, marker="o")
        self.axes.legend()

    def refresh(self):
        self.draw()

    def clear(self):
        self.axes.cla()

    def onclick(self, event):
        print('%s click: button=%d, x=%d, y=%d, xdata=%f, ydata=%f' %
              ('double' if event.dblclick else 'single', event.button,
               event.x, event.y, event.xdata, event.ydata))

    def draw_peakline(self, coordinates, color="blue"):

        #self.axes.cla()
        self.axes.axvline(coordinates[0], color=color)
        self.draw()
        print(coordinates[0])


class UiWindow(QtWidgets.QMainWindow, Ui_MainWindow):

    def __init__(self):
        QtWidgets.QMainWindow.__init__(self)
        self.setupUi(self)
        l = QtWidgets.QVBoxLayout(self.plotframe)
        self.sc = MyStaticMplCanvas(self.plotframe, width=5, height=4, dpi=100)

        self.toolbar = NavigationToolbar(self.sc, self)
        l.addWidget(self.sc)
        #l.addWidget(dc)
        l.addWidget(self.toolbar)

        self.setWindowTitle("Manifestum 1.0 Plot Gui")
        self.setWindowIcon(QtGui.QIcon("logo.png"))
        self.plotLabelLineEdit.setText("Default")


class UiManager:

    def __init__(self, window, speclist, session_id):

        self.window = window
        self.speclist = speclist
        self.session_id = session_id


        cid = self.window.sc.mpl_connect('button_press_event', self.mouse_handler)
        self.window.setLabel.clicked.connect(self.set_title)
        self.window.bApply.clicked.connect(self.apply)
        self.window.clearNotes.clicked.connect(self.clear_notes)
        self.window.bClear.clicked.connect(self.clear_peaks)
        self.window.exportNotes.clicked.connect(self.export_notes)
        self.window.exportPeaks.clicked.connect(self.export_peaks)

        #self.window.sc.plot_sfg(self.subset, flag=self.flag)
        self.window.show()

        self.window.bIntegrate.clicked.connect(self.integrate)
        self.modus = "normal"
        self.counter = 0
        self.leftborder = None
        self.rightborder = None

        #Plotting system state
        self.normalized = False
        self.ir = False
        self.vis = False
        self.smooth = False
        self.root = False
        self.normalized_to_highest = False
        self.plot_handler()


    def test(self):
        normalized = self.window.checkBox.checkState()
        show_IR = self.window.checkBox_2.checkState()
        show_Vis = self.window.checkBox_3.checkState()
        ntitle = self.window.plotLabelLineEdit.text()
        mode = self.window.cStacked.checkState()

        if mode != 0:
            if show_IR != 0:
                if show_Vis != 0:
                    self.flag = "b"

                else:
                    self.flag ="IR"
            elif show_Vis != 0:
                self.flag = "Vis"

            elif normalized == 0:
                self.flag = "r"
            else:
                #todo add the possbility of pure raw plotting
                self.flag= "none"
            print(self.flag)
            self.window.sc.plot_sfg(self.subset, title=ntitle, flag=self.flag)

        else:
            if show_IR != 0:
                if show_Vis != 0:
                    self.flag = "b"

                else:
                    self.flag = "IR"
            elif show_Vis != 0:
                self.flag = "Vis"

            elif normalized == 0:
                self.flag = "r"
            else:
                self.flag = "none"

            self.window.sc.plot_sfg([self.active_spectrum], title=ntitle, flag=self.flag)
        self.window.sc.draw()
        print("boing blöööök bumm miau")
        print(normalized, show_IR, show_Vis)

    def pick(self, event):
            x = event.xdata
            y = event.ydata

            if self.modus == "normal":
             self.window.Peaks.insertPlainText("Peak: "+str(x)+","+str(y)+"\n")
             self.window.sc.draw_peakline((x, y))

            else:
                self.window.sc.draw_peakline((x, y), color="red")

    def clear_peaks(self):
        self.window.Peaks.clear()

    def clear_notes(self):
        self.window.Peaks_2.clear()

    def export_peaks(self):

        text = self.window.Peaks.toPlainText()
        print(self.session_id)

        with open(self.session_id+".txt", "a") as outfile:
            outfile.write("Peak export:\n")
            outfile.write(text)
            outfile.write("\n"+"*"*80+"\n")

    def export_notes(self):

        text = self.window.Peaks_2.toPlainText()

        with open(self.session_id+".txt", "a") as outfile:
            outfile.write("Notes export:\n")
            outfile.write(text)
            outfile.write("\n" + "*" * 80+"\n")

    def integrate(self):

        if len(self.speclist == 1):

            if self.counter == 0:
                print("Select first peak: ")
                self.modus = "integrate"

            elif self.counter == 1:
                print("Select second peak: ")
                self.modus = "next"
        else:
            print("Integration not possible for mor than one active spectrum!\n")

    def set_title(self):
        text = self.window.plotLabelLineEdit.text()
        self.window.sc.axes.set_title(text)
        self.window.sc.draw()
        self.window.plotLabelLineEdit.setText("")

    def mouse_handler(self, event):

        if self.modus == "normal":
            if event.button == 3 and self.window.pickingModeBox.isChecked():
                self.pick(event)

        elif self.modus == "integrate":
            if event.button == 3:
                self.counter = 1
                self.integrate(event)
                self.pick(event)
                self.leftborder = (event.xdata, event.ydata)

        elif self.modus == "next" and event.button == 3:
            self.counter += 1
            self.integrate(event)
            self.pick(event)
            self.rightborder = (event.xdata, event.ydata)

            print("Borders selected",self.leftborder, self.rightborder)

            spectrum = self.speclist[0]  # type: SfgSpectrum
            point_list = spectrum.create_pointlist(spectrum.normalized_intensity[::-1])
            left_index = self.get_closest_index(point_list, self.leftborder)
            right_index = self.get_closest_index(point_list, self.rightborder)
            x_array = (spectrum.wavenumbers[::-1])[left_index:right_index+1]

            if self.normalized is True:
                y_array = (spectrum.normalized_intensity[::-1])[left_index:right_index+1]
            else:
                y_array = (spectrum.raw_intensity[::-1])[left_index:right_index + 1]

            area = spectrum.integrate_peak(x_array, y_array)
            self.window.Peaks.insertPlainText("Integral between "+ str(self.leftborder)
                                              + " and " + str(self.rightborder) + ": "
                                              + str(area) + "\n")

            self.modus = "normal"
            self.counter = 0

    def plot_handler(self):

        self.window.sc.clear()
        self.window.sc.init_axes()
        norm_factor = 0

        if self.normalized_to_highest is True:
            for spectrum in self.speclist:
                temp = spectrum.yield_maximum()
                if temp > norm_factor:
                    norm_factor = temp

        for spectrum in self.speclist:
            x_data = spectrum.wavenumbers
            name = spectrum.name.full_name[:-4]

            if self.normalized_to_highest is True:
                y_data = spectrum.normalize_to_highest(external_norm=norm_factor)
                self.window.sc.plot_dataset(x_data, y_data, name, "")

            else:
                if self.normalized is True:
                    self.window.sc.plot_dataset(x_data, spectrum.normalized_intensity, name, "")
                elif self.normalized is False:

                    if self.smooth is False:
                        self.window.sc.plot_dataset(x_data, spectrum.raw_intensity, name, "_raw")
                    else:
                        self.window.sc.plot_dataset(x_data, spectrum.smooth(), name, "_smoothed")
                if self.ir is True:
                    self.window.sc.plot_dataset(x_data, spectrum.ir_intensity, name, "_IR")

                if self.vis is True:
                    self.window.sc.plot_dataset(x_data, spectrum.vis_intensity, name, "_Vis")

        self.window.sc.refresh()

    def apply(self):

        self.normalized = self.window.normalizedBox.isChecked()
        self.ir = self.window.showIrBox.isChecked()
        self.vis = self.window.showVisBox.isChecked()
        self.smooth = self.window.smoothBox.isChecked()
        self.root = self.window.rootBox.isChecked()
        self.normalized_to_highest = self.window.checkBox_6.isChecked()

        self.plot_handler()

    def get_closest_index(self, array_datapoints, check):

        d = 1000000000000
        index = None
        for point in array_datapoints:

            d_temp = np.sqrt((check[0] - point[0]) ** 2 + (check[1] - point[1]) ** 2)
            if d_temp < d:
                d = d_temp
                index = point[2]

        return index


def run_app(speclist, session_id):
    qApp = QtWidgets.QApplication(sys.argv)
    ui = UiWindow()
    M = UiManager(ui,speclist, session_id)
    sys.exit(qApp.exec_())



#run_app([])

