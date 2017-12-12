from __future__ import unicode_literals
import numpy as np
from scipy.integrate import simps as sp

import time
import sys
import os
import random
import matplotlib
# Make sure that we are using QT5
matplotlib.use('Qt5Agg')
from PyQt5 import QtCore, QtWidgets

from numpy import arange, sin, pi
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
from matplotlib import rcParams
rcParams['mathtext.default'] = 'regular'


from SFGPlot import Ui_MainWindow
from ipy_interpreter import IpyInterpreter

def integrate_peak(x_array, y_array):
    area = 0
    for i in range(len(x_array) - 1):
        dx = abs(x_array[i + 1] - x_array[i])
        square = dx * y_array[i]
        triangle = dx * abs(y_array[i + 1] - y_array[i]) * 0.5
        total = square + triangle
        area += total
    return area
"""
x = np.linspace(0,5,2)
y = (2*x**2)+(0.4*x)-3
y2 = ((6*x**3)+(2*x**2)+(0.4*x)-3)**2
y3 = (0.5*np.sin(x)+2*x)**2

als = integrate_peak(x,y)
aln = np.trapz(y,x)

a2i = integrate_peak(x,y2)
a3i = integrate_peak(x,y3)
a2n = integrate_peak(y2,x)
a3n = integrate_peak(y3,x)

a1s = sp(y, x)
a2s = sp(y2, x)
a3s = sp(y3, x)

print("My Integrate")
print(als)
print("\n")

print("Numpy")
print(aln)
print("\n")

print("Simpson")
print(a1s)
print("\n")

print("My Integrate")
print(a2i)
print("\n")

print("Numpy")
print(a2n)
print("\n")

print("Simpson")
print(a2s)
print("\n")

print("My Integrate")
print(a3i)
print("\n")

print("Numpy")
print(a3n)
print("\n")

print("Simpson")
print(a3s)
print("\n")
"""
# embedding_in_qt5.py --- Simple Qt5 application embedding matplotlib canvases
#
# Copyright (C) 2005 Florent Rougon
#               2006 Darren Dale
#               2015 Jens H Nielsen
#
# This file is an example program for matplotlib. It may be used and
# modified with no restriction; raw copies as well as modified versions
# may be distributed without limitation.



progname = os.path.basename(sys.argv[0])
progversion = "0.1"


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
        self.axes.set_ylabel("Intensity/ a.u.")

    def compute_initial_figure(self):
        pass


class MyStaticMplCanvas(MyMplCanvas):
    """Simple canvas with a sine plot."""

    def compute_initial_figure(self):
        t = arange(0.0, 3.0, 0.01)
        s = sin(2*pi*t)
        self.axes.plot(t, s)

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
                self.axes.plot(spectrum.wavenumbers, spectrum.normalized_intensity, label=spectrum.name.full_name, linestyle='dotted', markersize=4, marker="o")


        #box = self.axes.get_position()
        #self.axes.set_position([box.x0, box.y0, box.width * 0.8, box.height])

        #self.axes.legend(loc='center left', bbox_to_anchor=(1, 0.5))


    def clear(self):
        print("Debuggedidooo")
        self.axes.cla()

    def onclick(self, event):
        print('%s click: button=%d, x=%d, y=%d, xdata=%f, ydata=%f' %
              ('double' if event.dblclick else 'single', event.button,
               event.x, event.y, event.xdata, event.ydata))




class MyDynamicMplCanvas(MyMplCanvas):
    """A canvas that updates itself every second with a new plot."""

    def __init__(self, *args, **kwargs):
        MyMplCanvas.__init__(self, *args, **kwargs)
        timer = QtCore.QTimer(self)
        timer.timeout.connect(self.update_figure)
        timer.start(1000)

    def compute_initial_figure(self):
        self.axes.plot([0, 1, 2, 3], [1, 2, 0, 4], 'r')

    def update_figure(self):
        # Build a list of 4 random integers between 0 and 10 (both inclusive)
        l = [random.randint(0, 10) for i in range(4)]
        self.axes.cla()
        self.axes.plot([0, 1, 2, 3], l, 'r')
        self.draw()


class ApplicationWindow(QtWidgets.QMainWindow):
    def __init__(self):
        QtWidgets.QMainWindow.__init__(self)
        self.setAttribute(QtCore.Qt.WA_DeleteOnClose)
        self.setWindowTitle("application main window")

        self.file_menu = QtWidgets.QMenu('&File', self)
        self.file_menu.addAction('&Quit', self.fileQuit,
                                 QtCore.Qt.CTRL + QtCore.Qt.Key_Q)
        self.menuBar().addMenu(self.file_menu)

        self.help_menu = QtWidgets.QMenu('&Help', self)
        self.menuBar().addSeparator()
        self.menuBar().addMenu(self.help_menu)

        self.help_menu.addAction('&About', self.about)

        self.main_widget = QtWidgets.QWidget(self)

        l = QtWidgets.QVBoxLayout(self.main_widget)
        sc = MyStaticMplCanvas(self.main_widget, width=5, height=4, dpi=100)
        dc = MyDynamicMplCanvas(self.main_widget, width=5, height=4, dpi=100)
        self.toolbar = NavigationToolbar(sc, self)
        l.addWidget(sc)
        l.addWidget(dc)
        l.addWidget(self.toolbar)

        self.main_widget.setFocus()
        self.setCentralWidget(self.main_widget)

        self.statusBar().showMessage("All hail matplotlib!", 2000)

    def fileQuit(self):
        self.close()

    def closeEvent(self, ce):
        self.fileQuit()

    def about(self):
        QtWidgets.QMessageBox.about(self, "About",
                                    """embedding_in_qt5.py example
Copyright 2005 Florent Rougon, 2006 Darren Dale, 2015 Jens H Nielsen

This program is a simple example of a Qt5 application embedding matplotlib
canvases.

It may be used and modified with no restriction; raw copies as well as
modified versions may be distributed without limitation.

This is modified from the embedding in qt4 example to show the difference
between qt4 and qt5"""
                                )

""""
qApp = QtWidgets.QApplication(sys.argv)

aw = ApplicationWindow()
aw.setWindowTitle("%s" % progname)
aw.show()
sys.exit(qApp.exec_())
#qApp.exec_()
"""


class UiWindow(QtWidgets.QMainWindow, Ui_MainWindow):

    def __init__(self):
        QtWidgets.QMainWindow.__init__(self)
        self.setupUi(self)
        l = QtWidgets.QVBoxLayout(self.plotframe)
        self.sc = MyStaticMplCanvas(self.plotframe, width=5, height=4, dpi=100)

        #dc = MyDynamicMplCanvas(self.plotframe, width=5, height=4, dpi=100)
        self.toolbar = NavigationToolbar(self.sc, self)
        l.addWidget(self.sc)
        #l.addWidget(dc)
        l.addWidget(self.toolbar)

        self.setWindowTitle("SFG Plot Gui")
        self.plotLabelLineEdit.setText("Default")





class UiManager:

    def __init__(self, subset, window):
        self.subset = subset
        #todo calling the constructor of the window here, removing it from arglist
        self.window = window
        self.flag = "none"
        if self.subset[0]:
            self.active_spectrum = self.subset[0]
        else:
            self.active_spectrum = "none"
        self.index = 0
        self.peakcounter = 1
        self.exports = 1

        #todo implementation of toolbar and PlotCanvas missing here


        self.fill_form()
        self.window.bRefresh.clicked.connect(self.test)
        self.window.bClear.clicked.connect(self.clear_peaks)
        self.window.bExport.clicked.connect(self.export)
        #self.window.pushButton.clicked.connect(self.next)
        cid = self.window.sc.mpl_connect('button_press_event', self.pick)

        self.window.sc.plot_sfg(self.subset, flag=self.flag)
        self.window.show()


    def fill_form(self):
        #individual properties
        self.window.fullNameLineEdit.setText(self.active_spectrum.name.full_name)
        self.window.dateLineEdit.setText(self.active_spectrum.name.date)
        self.window.surfactantLineEdit.setText(self.active_spectrum.name.surfactant)
        self.window.spreadingVolumeLineEdit.setText(self.active_spectrum.name.surfactant_spread_volume)
        self.window.spreadingVolumeLineEdit_2.setText(self.active_spectrum.name.sensitizer_spread_volume)
        self.window.sensitizerLineEdit.setText(self.active_spectrum.name.sensitizer)
        self.window.photolysisLineEdit.setText(self.active_spectrum.name.photolysis)
        self.window.specralRangeLineEdit.setText(str(self.active_spectrum.yield_spectral_range()))
        self.window.Comment.insertPlainText(self.active_spectrum.name.comment)
        #collective properties
        self.window.totalSpectraCountLineEdit.setText(str(len(self.subset)))
        self.window.plotLabelLineEdit.setText("default")

    def test(self):
        normalized = self.window.checkBox.checkState()
        show_IR = self.window.checkBox_2.checkState()
        show_Vis = self.window.checkBox_3.checkState()
        ntitle = self.window.plotLabelLineEdit.text()


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
        self.window.sc.draw()

        print("boing blöööök bumm miau")
        print(normalized, show_IR, show_Vis)

    def printy(self):
        print("blööök"*1000)

    def next(self):

        self.index += 1
        try:
            self.active_spectrum = self.subset[self.index]
        except IndexError:
            self.active_spectrum= self.subset[0]
            self.index = 0
            self.window.sc.plot_sfg(self.subset, title="default", flag=self.flag)
            self.window.sc.draw()

        self.fill_form()

    def pick(self, event):

        if event.button == 3:
            x = event.xdata
            y = event.ydata
            self.window.Peaks.insertPlainText(str(x)+"\n"+str(y)+"\n\n")
        print('%s click: button=%d, x=%d, y=%d, xdata=%f, ydata=%f' %
              ('double' if event.dblclick else 'single', event.button,
               event.x, event.y, event.xdata, event.ydata))

    def clear_peaks(self):
        self.window.Peaks.clear()

    def export(self):
        timestamp = time.strftime("%b %d %Y %H:%M:%S")
        title = self.window.plotLabelLineEdit.text()
        with open(str(self.exports)+".txt", "w") as outfile:
            text = self.window.Peaks.toPlainText()
            outfile.write(timestamp+" "+title+"\n\n")
            outfile.write(text)
            self.exports += 1


#test code section
qApp = QtWidgets.QApplication(sys.argv)
ui = UiWindow()

I = IpyInterpreter()
I.get("bo")
#I.subset=[I.subset[0]]

um = UiManager(I.subset, ui)

sys.exit(qApp.exec_())
