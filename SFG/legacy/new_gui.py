from plotgui import Ui_MainWindow
from SFG.legacy.ltgui import Ui_MainWindow as LtW

import time
import sys
import matplotlib
import traceback
import logging
# Make sure that we are using QT5
matplotlib.use('Qt5Agg')
from PyQt5 import QtWidgets, QtGui


import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
from matplotlib import rcParams
import numpy as np
from scipy.optimize import curve_fit

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

    def draw_peakline(self, coordinates, color="blue"):

        #self.axes.cla()
        self.axes.axvline(coordinates[0], color=color)
        self.draw()
        print(coordinates[0])

    def refresh(self):
        self.draw()

    def clear(self):
        self.axes.cla()


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



    def onclick(self, event):
        print('%s click: button=%d, x=%d, y=%d, xdata=%f, ydata=%f' %
              ('double' if event.dblclick else 'single', event.button,
               event.x, event.y, event.xdata, event.ydata))


class LtCanvas(MyMplCanvas):
    """The graphical represantion of the LtIsotherm plot to be embedded in the Qt window."""

    def init_axes(self,title="default"):
        self.axes.grid()
        self.axes.set_title(title)
        self.axes.set_xlabel("time/ s")
        self.axes.set_ylabel("surface pressure/ mNm$^{-1}$")

    def plot_dataset(self, x_dataset, y_dataset, label):
        self.axes.plot(x_dataset, y_dataset, label=label)


class UiWindow(QtWidgets.QMainWindow, Ui_MainWindow):
    """The window to handle SFG plots."""

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


class LtWindow(QtWidgets.QMainWindow, LtW):
    """The GUI to handle Langmuir trough isotherms."""

    def __init__(self):
        QtWidgets.QMainWindow.__init__(self)
        self.setupUi(self)
        l = QtWidgets.QVBoxLayout(self.plotframe)
        self.sc = LtCanvas(self.plotframe)
        self.toolbar = NavigationToolbar(self.sc, self)
        l.addWidget(self.sc)
        l.addWidget(self.toolbar)

        self.sc.init_axes()


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

        if len(self.speclist) == 1:

            if self.counter == 0:

                print("Select first peak1: ")
                self.modus = "integrate"

            elif self.counter == 1:
                print("Select second peak: ")
                self.modus = "next"

            else:
                print("Undefined counter state!")
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
                self.pick(event)
                self.leftborder = (event.xdata, event.ydata)
                self.counter += 1
                self.integrate()

        elif self.modus == "next" and event.button == 3:

            self.pick(event)
            self.rightborder = (event.xdata, event.ydata)

            print("Borders selected",self.leftborder, self.rightborder)

            spectrum = self.speclist[0]  # type: SfgSpectrum
            point_list = spectrum.create_pointlist(spectrum.normalized_intensity[::-1])

            print(point_list, self.leftborder)
            left_index = self.get_closest_index(point_list, self.leftborder)
            right_index = self.get_closest_index(point_list, self.rightborder)

            x_array = (spectrum.wavenumbers[::-1])[left_index:right_index+1]

            if self.normalized is True:
                if self.root is False:
                 y_array = (spectrum.normalized_intensity[::-1])[left_index:right_index+1]
                else:
                    y_array = (spectrum.root()[::-1])[left_index:right_index+1]
            else:
                y_array = (spectrum.raw_intensity[::-1])[left_index:right_index + 1]

            area = spectrum.integrate_peak(x_array, y_array)


            l = (f'{self.leftborder[0]:.2f}')+","+(f'{self.leftborder[1]:.2f}')
            r = (f'{self.rightborder[0]:.2f}')+","+(f'{self.rightborder[1]:.2f}')


            self.window.Peaks.insertPlainText("Integral between " + l
                                              + " and " +  r + ": "
                                              + str(area) + "\n")

            self.window.sc.axes.fill_between(x_array,y_array)
            self.window.sc.refresh()

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

                    if self.root is True:

                        self.window.sc.plot_dataset(x_data, spectrum.root(), name, "")
                    else:
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


class LtUiManager:

    def __init__(self, window, isotherm_list, session_id):

        self.mode_dic = {
            "apm":"area per molecule/ $\AA^{2}$ ",
            "time":"time/ s",
            "area":"area/ cm$^{-2}$"
        }

        self.window = window # type: LtW
        self.isotherm_list = isotherm_list
        self.session_id = session_id
        self.mode = time
        self.modus = "normal"
        self.borders = []

        cid = self.window.sc.mpl_connect('button_press_event', self.mouse_handler)

        self.window.pbsave.clicked.connect(self.export)
        self.window.pbapply.clicked.connect(self.apply)
        self.window.pbselect.clicked.connect(self.select_borders)

        self.apply()
        self.window.show()

    def get_closest_index(self, array_datapoints, check):

        d = 1000000000000
        index = None
        for point in array_datapoints:

            d_temp = np.sqrt((check[0] - point[0]) ** 2 + (check[1] - point[1]) ** 2)
            if d_temp < d:
                d = d_temp
                index = point[2]

        return index

    def export(self):

        text = self.window.teout.toPlainText()

        with open(self.session_id + ".txt", "a") as outfile:
            outfile.write("Notes export:\n")
            outfile.write(text)
            outfile.write("\n" + "*" * 80 + "\n")

    def clear(self):
        self.window.teout.clear()

    def apply(self):
        self.mode = self.window.ycb.currentText()

        if self.mode == "area per molecule":
            self.mode = "apm"

        if self.mode == "total area":
            self.mode = "area"

        self.plot_handler()

    def plot_handler(self):

        self.window.sc.clear()
        self.window.sc.init_axes()

        for isotherm in self.isotherm_list:
            self.window.sc.plot_dataset(getattr(isotherm, self.mode), isotherm.pressure, isotherm.name)

        try:
            if self.window.cbshow.isChecked() and self.mode == "area":
                self.window.sc.plot_dataset(isotherm.area[:-1], isotherm.calculate_elasticity(), label="elastic")
        except Exception as e:
            logging.error(traceback.format_exc())

        self.window.sc.axes.set_xlabel(self.mode_dic[self.mode])
        self.window.sc.axes.legend()
        self.window.sc.refresh()

    def mouse_handler(self, event):

        if self.modus == "normal":
            if event.button == 3 and self.window.cbpick.isChecked():
                self.pick(event)

        elif self.modus == "left":

            if event.button == 3:

                self.borders.append((event.xdata, event.ydata))
                self.pick(event)
                self.modus = "right"

        elif self.modus == "right":

            if event.button == 3:
                self.borders.append((event.xdata, event.ydata))
                self.pick(event)
                self.fit()
                self.borders = []
                self.modus = "normal"

    def pick(self, event):
            x = event.xdata
            y = event.ydata

            if self.modus == "normal":
             x_ = f'x_data: {x:.2f}'
             y_ = f'y_data: {y:.2f}'

             self.window.teout.insertPlainText(x_+" "+y_+"\n")
             self.window.sc.draw_peakline((x, y))

            else:
                self.window.sc.draw_peakline((x, y), color="red")

    def select_borders(self):
        if self.modus == "normal":
            self.modus = "left"

    def extract_indices(self):

        if len(self.isotherm_list) != 1:
            raise ValueError("Operation not defined for more than one isotherm")

        else:
            x_data = getattr(self.isotherm_list[0], self.mode)
            values = self.isotherm_list[0].create_pointlist(x_data)


            index_left = self.get_closest_index(values, self.borders[0])
            index_right = self.get_closest_index(values, self.borders[1])

            return (index_left, index_right)

    def extract_slice(self):

        indices = self.extract_indices()
        iso = self.isotherm_list[0]
        xdata = getattr(iso, self.mode)
        sliced = iso.get_slice(xdata, indices[0], indices[1])

        return sliced

    def fit(self):

        data = self.extract_slice()
        amplitude = data[1][0]



        fit_function = fit_function_mono
        fit_function2 = fit_function_di



        try:
         offset = data[0][0]
         xdata = data[0]



         popt1, pcov = curve_fit(fit_function, xdata, data[1],p0=(60, 1/307, 7))
         popt2, pcov = curve_fit(fit_function2, xdata, data[1], p0=(1, 1e-2, 1, 1e-2, 1, 1))


         plt.plot(xdata, data[1])
         plt.scatter(testx, testy)


         mono_label = f'${popt1[0]:.2f}*e^{{-{popt1[1]:.6f}*x}}+{popt1[2]:.4f}$'
         di_label = f'$({popt1[0]:.2f}*e^{{-{popt2[1]:.6f}*x}}+{popt2[2]:.4f})+({popt2[3]:.2f}*e^{{-{popt2[4]:.6f}*x}}+{popt2[5]:.4f})$'

         plt.plot(xdata, fit_function(data[0], *popt1), label=mono_label)
         plt.plot(xdata, fit_function2(data[0], *popt2), label=di_label)
         plt.plot(testx, fit_function(testx, *pre_opt), label="pre")


         plt.legend(prop={'size':7})
         plt.show()
         self.window.teout.insertPlainText(f'rate constant: {popt1[0]:.6f}\n')
        except Exception as e:
            logging.error(traceback.format_exc())


def fit_function_mono(x, a, b, c):
    return a*np.exp(-b*x)+ c


def fit_function_di(x, a, b, c, d, e, f):
    return a*np.exp(-b*x)+ c + f*np.exp(-d*x)+e


def run_app(speclist, session_id):
    qApp = QtWidgets.QApplication(sys.argv)
    ui = UiWindow()
    M = UiManager(ui, speclist, session_id)
    sys.exit(qApp.exec_())


def run_lt_app(isotherm_list, session_id):
    qApp = QtWidgets.QApplication(sys.argv)
    L = LtWindow()
    M = LtUiManager(L, isotherm_list, session_id)
    sys.exit(qApp.exec_())


#run_app([])

