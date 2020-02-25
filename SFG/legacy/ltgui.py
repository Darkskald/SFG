# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file '.\ltgui.ui'
#
# Created by: PyQt5 UI code generator 5.6
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(634, 677)
        self.centralwidget = QtWidgets.QWidget(MainWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.verticalLayout = QtWidgets.QVBoxLayout(self.centralwidget)
        self.verticalLayout.setObjectName("verticalLayout")
        self.plotframe = QtWidgets.QWidget(self.centralwidget)
        self.plotframe.setObjectName("plotframe")
        self.verticalLayout.addWidget(self.plotframe)
        self.toolframe = QtWidgets.QWidget(self.centralwidget)
        self.toolframe.setObjectName("toolframe")
        self.gridLayout = QtWidgets.QGridLayout(self.toolframe)
        self.gridLayout.setContentsMargins(0, 0, 0, 0)
        self.gridLayout.setObjectName("gridLayout")
        self.pbintegrate = QtWidgets.QPushButton(self.toolframe)
        self.pbintegrate.setObjectName("pbintegrate")
        self.gridLayout.addWidget(self.pbintegrate, 1, 5, 1, 1)
        self.cbpick = QtWidgets.QCheckBox(self.toolframe)
        self.cbpick.setObjectName("cbpick")
        self.gridLayout.addWidget(self.cbpick, 1, 4, 1, 1)
        self.teout = QtWidgets.QPlainTextEdit(self.toolframe)
        self.teout.setObjectName("teout")
        self.gridLayout.addWidget(self.teout, 1, 1, 3, 2)
        spacerItem = QtWidgets.QSpacerItem(10, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.gridLayout.addItem(spacerItem, 3, 0, 1, 1)
        self.cbshow = QtWidgets.QCheckBox(self.toolframe)
        self.cbshow.setObjectName("cbshow")
        self.gridLayout.addWidget(self.cbshow, 2, 4, 1, 1)
        self.pbselect = QtWidgets.QPushButton(self.toolframe)
        self.pbselect.setObjectName("pbselect")
        self.gridLayout.addWidget(self.pbselect, 2, 5, 1, 1)
        self.pbclear = QtWidgets.QPushButton(self.toolframe)
        self.pbclear.setObjectName("pbclear")
        self.gridLayout.addWidget(self.pbclear, 4, 2, 1, 1)
        spacerItem1 = QtWidgets.QSpacerItem(136, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.gridLayout.addItem(spacerItem1, 1, 3, 1, 1)
        spacerItem2 = QtWidgets.QSpacerItem(136, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.gridLayout.addItem(spacerItem2, 2, 3, 1, 1)
        self.pbsave = QtWidgets.QPushButton(self.toolframe)
        self.pbsave.setObjectName("pbsave")
        self.gridLayout.addWidget(self.pbsave, 4, 1, 1, 1)
        self.ycb = QtWidgets.QComboBox(self.toolframe)
        self.ycb.setObjectName("ycb")
        self.ycb.addItem("")
        self.ycb.addItem("")
        self.ycb.addItem("")
        self.gridLayout.addWidget(self.ycb, 0, 5, 1, 1)
        self.pbapply = QtWidgets.QPushButton(self.toolframe)
        self.pbapply.setObjectName("pbapply")
        self.gridLayout.addWidget(self.pbapply, 0, 4, 1, 1)
        self.verticalLayout.addWidget(self.toolframe)
        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 634, 21))
        self.menubar.setObjectName("menubar")
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(MainWindow)
        self.statusbar.setObjectName("statusbar")
        MainWindow.setStatusBar(self.statusbar)

        self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "MainWindow"))
        self.pbintegrate.setText(_translate("MainWindow", "Integrate"))
        self.cbpick.setText(_translate("MainWindow", "Pick"))
        self.cbshow.setText(_translate("MainWindow", "Show Elasticity"))
        self.pbselect.setText(_translate("MainWindow", "Select Borders"))
        self.pbclear.setText(_translate("MainWindow", "Clear"))
        self.pbsave.setText(_translate("MainWindow", "Save"))
        self.ycb.setItemText(0, _translate("MainWindow", "time"))
        self.ycb.setItemText(1, _translate("MainWindow", "area per molecule"))
        self.ycb.setItemText(2, _translate("MainWindow", "total area"))
        self.pbapply.setText(_translate("MainWindow", "Apply"))

