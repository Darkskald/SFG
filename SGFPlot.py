# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file '.\SFGPlot.ui'
#
# Created by: PyQt5 UI code generator 5.6
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(892, 856)
        self.centralwidget = QtWidgets.QWidget(MainWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.formLayoutWidget = QtWidgets.QWidget(self.centralwidget)
        self.formLayoutWidget.setGeometry(QtCore.QRect(700, 20, 191, 421))
        self.formLayoutWidget.setObjectName("formLayoutWidget")
        self.formLayout = QtWidgets.QFormLayout(self.formLayoutWidget)
        self.formLayout.setContentsMargins(0, 0, 0, 0)
        self.formLayout.setObjectName("formLayout")
        self.fullNameLabel = QtWidgets.QLabel(self.formLayoutWidget)
        self.fullNameLabel.setObjectName("fullNameLabel")
        self.formLayout.setWidget(0, QtWidgets.QFormLayout.LabelRole, self.fullNameLabel)
        self.fullNameLineEdit = QtWidgets.QLineEdit(self.formLayoutWidget)
        self.fullNameLineEdit.setObjectName("fullNameLineEdit")
        self.formLayout.setWidget(0, QtWidgets.QFormLayout.FieldRole, self.fullNameLineEdit)
        self.dateLabel = QtWidgets.QLabel(self.formLayoutWidget)
        self.dateLabel.setObjectName("dateLabel")
        self.formLayout.setWidget(1, QtWidgets.QFormLayout.LabelRole, self.dateLabel)
        self.dateLineEdit = QtWidgets.QLineEdit(self.formLayoutWidget)
        self.dateLineEdit.setObjectName("dateLineEdit")
        self.formLayout.setWidget(1, QtWidgets.QFormLayout.FieldRole, self.dateLineEdit)
        self.surfactantLineEdit = QtWidgets.QLineEdit(self.formLayoutWidget)
        self.surfactantLineEdit.setObjectName("surfactantLineEdit")
        self.formLayout.setWidget(2, QtWidgets.QFormLayout.FieldRole, self.surfactantLineEdit)
        self.spreadingVolumeLabel = QtWidgets.QLabel(self.formLayoutWidget)
        self.spreadingVolumeLabel.setObjectName("spreadingVolumeLabel")
        self.formLayout.setWidget(3, QtWidgets.QFormLayout.LabelRole, self.spreadingVolumeLabel)
        self.spreadingVolumeLineEdit = QtWidgets.QLineEdit(self.formLayoutWidget)
        self.spreadingVolumeLineEdit.setObjectName("spreadingVolumeLineEdit")
        self.formLayout.setWidget(3, QtWidgets.QFormLayout.FieldRole, self.spreadingVolumeLineEdit)
        self.sensitizerLabel = QtWidgets.QLabel(self.formLayoutWidget)
        self.sensitizerLabel.setObjectName("sensitizerLabel")
        self.formLayout.setWidget(4, QtWidgets.QFormLayout.LabelRole, self.sensitizerLabel)
        self.sensitizerLineEdit = QtWidgets.QLineEdit(self.formLayoutWidget)
        self.sensitizerLineEdit.setObjectName("sensitizerLineEdit")
        self.formLayout.setWidget(4, QtWidgets.QFormLayout.FieldRole, self.sensitizerLineEdit)
        self.spreadingVolumeLabel_2 = QtWidgets.QLabel(self.formLayoutWidget)
        self.spreadingVolumeLabel_2.setObjectName("spreadingVolumeLabel_2")
        self.formLayout.setWidget(5, QtWidgets.QFormLayout.LabelRole, self.spreadingVolumeLabel_2)
        self.spreadingVolumeLineEdit_2 = QtWidgets.QLineEdit(self.formLayoutWidget)
        self.spreadingVolumeLineEdit_2.setObjectName("spreadingVolumeLineEdit_2")
        self.formLayout.setWidget(5, QtWidgets.QFormLayout.FieldRole, self.spreadingVolumeLineEdit_2)
        self.photolysisLabel = QtWidgets.QLabel(self.formLayoutWidget)
        self.photolysisLabel.setObjectName("photolysisLabel")
        self.formLayout.setWidget(6, QtWidgets.QFormLayout.LabelRole, self.photolysisLabel)
        self.photolysisLineEdit = QtWidgets.QLineEdit(self.formLayoutWidget)
        self.photolysisLineEdit.setObjectName("photolysisLineEdit")
        self.formLayout.setWidget(6, QtWidgets.QFormLayout.FieldRole, self.photolysisLineEdit)
        self.specralRangeLabel = QtWidgets.QLabel(self.formLayoutWidget)
        self.specralRangeLabel.setObjectName("specralRangeLabel")
        self.formLayout.setWidget(7, QtWidgets.QFormLayout.LabelRole, self.specralRangeLabel)
        self.specralRangeLineEdit = QtWidgets.QLineEdit(self.formLayoutWidget)
        self.specralRangeLineEdit.setObjectName("specralRangeLineEdit")
        self.formLayout.setWidget(7, QtWidgets.QFormLayout.FieldRole, self.specralRangeLineEdit)
        spacerItem = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.formLayout.setItem(8, QtWidgets.QFormLayout.LabelRole, spacerItem)
        self.commentLabel = QtWidgets.QLabel(self.formLayoutWidget)
        self.commentLabel.setObjectName("commentLabel")
        self.formLayout.setWidget(9, QtWidgets.QFormLayout.LabelRole, self.commentLabel)
        self.Comment = QtWidgets.QTextEdit(self.formLayoutWidget)
        self.Comment.setObjectName("Comment")
        self.formLayout.setWidget(10, QtWidgets.QFormLayout.SpanningRole, self.Comment)
        self.surfactantLabel = QtWidgets.QLabel(self.formLayoutWidget)
        self.surfactantLabel.setObjectName("surfactantLabel")
        self.formLayout.setWidget(2, QtWidgets.QFormLayout.LabelRole, self.surfactantLabel)
        self.plotLabelLabel = QtWidgets.QLabel(self.formLayoutWidget)
        self.plotLabelLabel.setObjectName("plotLabelLabel")
        self.formLayout.setWidget(11, QtWidgets.QFormLayout.LabelRole, self.plotLabelLabel)
        self.plotLabelLineEdit = QtWidgets.QLineEdit(self.formLayoutWidget)
        self.plotLabelLineEdit.setObjectName("plotLabelLineEdit")
        self.formLayout.setWidget(11, QtWidgets.QFormLayout.FieldRole, self.plotLabelLineEdit)
        self.label = QtWidgets.QLabel(self.centralwidget)
        self.label.setGeometry(QtCore.QRect(690, 0, 47, 13))
        self.label.setObjectName("label")
        self.formLayoutWidget_2 = QtWidgets.QWidget(self.centralwidget)
        self.formLayoutWidget_2.setGeometry(QtCore.QRect(710, 450, 174, 51))
        self.formLayoutWidget_2.setObjectName("formLayoutWidget_2")
        self.formLayout_2 = QtWidgets.QFormLayout(self.formLayoutWidget_2)
        self.formLayout_2.setContentsMargins(0, 0, 0, 0)
        self.formLayout_2.setObjectName("formLayout_2")
        self.totalSpectraCountLabel = QtWidgets.QLabel(self.formLayoutWidget_2)
        self.totalSpectraCountLabel.setObjectName("totalSpectraCountLabel")
        self.formLayout_2.setWidget(0, QtWidgets.QFormLayout.LabelRole, self.totalSpectraCountLabel)
        self.totalSpectraCountLineEdit = QtWidgets.QLineEdit(self.formLayoutWidget_2)
        self.totalSpectraCountLineEdit.setObjectName("totalSpectraCountLineEdit")
        self.formLayout_2.setWidget(0, QtWidgets.QFormLayout.FieldRole, self.totalSpectraCountLineEdit)
        self.pushButton = QtWidgets.QPushButton(self.formLayoutWidget_2)
        self.pushButton.setObjectName("pushButton")
        self.formLayout_2.setWidget(1, QtWidgets.QFormLayout.FieldRole, self.pushButton)
        self.horizontalLayoutWidget = QtWidgets.QWidget(self.centralwidget)
        self.horizontalLayoutWidget.setGeometry(QtCore.QRect(0, 0, 691, 501))
        self.horizontalLayoutWidget.setObjectName("horizontalLayoutWidget")
        self.horizontalLayout = QtWidgets.QHBoxLayout(self.horizontalLayoutWidget)
        self.horizontalLayout.setContentsMargins(0, 0, 0, 0)
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.plotframe = QtWidgets.QWidget(self.horizontalLayoutWidget)
        self.plotframe.setObjectName("plotframe")
        self.horizontalLayout.addWidget(self.plotframe)
        self.gridLayoutWidget = QtWidgets.QWidget(self.centralwidget)
        self.gridLayoutWidget.setGeometry(QtCore.QRect(590, 540, 101, 71))
        self.gridLayoutWidget.setObjectName("gridLayoutWidget")
        self.gridLayout = QtWidgets.QGridLayout(self.gridLayoutWidget)
        self.gridLayout.setContentsMargins(0, 0, 0, 0)
        self.gridLayout.setObjectName("gridLayout")
        self.checkBox_2 = QtWidgets.QCheckBox(self.gridLayoutWidget)
        self.checkBox_2.setObjectName("checkBox_2")
        self.gridLayout.addWidget(self.checkBox_2, 2, 0, 1, 1)
        self.checkBox_3 = QtWidgets.QCheckBox(self.gridLayoutWidget)
        self.checkBox_3.setObjectName("checkBox_3")
        self.gridLayout.addWidget(self.checkBox_3, 3, 0, 1, 1)
        self.checkBox = QtWidgets.QCheckBox(self.gridLayoutWidget)
        self.checkBox.setObjectName("checkBox")
        self.gridLayout.addWidget(self.checkBox, 0, 0, 1, 1)
        self.horizontalLayoutWidget_2 = QtWidgets.QWidget(self.centralwidget)
        self.horizontalLayoutWidget_2.setGeometry(QtCore.QRect(420, 530, 160, 41))
        self.horizontalLayoutWidget_2.setObjectName("horizontalLayoutWidget_2")
        self.horizontalLayout_2 = QtWidgets.QHBoxLayout(self.horizontalLayoutWidget_2)
        self.horizontalLayout_2.setContentsMargins(0, 0, 0, 0)
        self.horizontalLayout_2.setObjectName("horizontalLayout_2")
        self.pushButton_4 = QtWidgets.QPushButton(self.horizontalLayoutWidget_2)
        self.pushButton_4.setObjectName("pushButton_4")
        self.horizontalLayout_2.addWidget(self.pushButton_4)
        self.pushButton_3 = QtWidgets.QPushButton(self.horizontalLayoutWidget_2)
        self.pushButton_3.setObjectName("pushButton_3")
        self.horizontalLayout_2.addWidget(self.pushButton_3)
        self.formLayoutWidget_3 = QtWidgets.QWidget(self.centralwidget)
        self.formLayoutWidget_3.setGeometry(QtCore.QRect(710, 520, 171, 111))
        self.formLayoutWidget_3.setObjectName("formLayoutWidget_3")
        self.formLayout_3 = QtWidgets.QFormLayout(self.formLayoutWidget_3)
        self.formLayout_3.setContentsMargins(0, 0, 0, 0)
        self.formLayout_3.setObjectName("formLayout_3")
        self.xLabelLabel = QtWidgets.QLabel(self.formLayoutWidget_3)
        self.xLabelLabel.setObjectName("xLabelLabel")
        self.formLayout_3.setWidget(0, QtWidgets.QFormLayout.LabelRole, self.xLabelLabel)
        self.xLabelLineEdit = QtWidgets.QLineEdit(self.formLayoutWidget_3)
        self.xLabelLineEdit.setObjectName("xLabelLineEdit")
        self.formLayout_3.setWidget(0, QtWidgets.QFormLayout.FieldRole, self.xLabelLineEdit)
        self.yLabelLabel = QtWidgets.QLabel(self.formLayoutWidget_3)
        self.yLabelLabel.setObjectName("yLabelLabel")
        self.formLayout_3.setWidget(1, QtWidgets.QFormLayout.LabelRole, self.yLabelLabel)
        self.yLabelLineEdit = QtWidgets.QLineEdit(self.formLayoutWidget_3)
        self.yLabelLineEdit.setObjectName("yLabelLineEdit")
        self.formLayout_3.setWidget(1, QtWidgets.QFormLayout.FieldRole, self.yLabelLineEdit)
        self.xMinLabel = QtWidgets.QLabel(self.formLayoutWidget_3)
        self.xMinLabel.setObjectName("xMinLabel")
        self.formLayout_3.setWidget(2, QtWidgets.QFormLayout.LabelRole, self.xMinLabel)
        self.xMinLineEdit = QtWidgets.QLineEdit(self.formLayoutWidget_3)
        self.xMinLineEdit.setObjectName("xMinLineEdit")
        self.formLayout_3.setWidget(2, QtWidgets.QFormLayout.FieldRole, self.xMinLineEdit)
        self.xMaxLabel = QtWidgets.QLabel(self.formLayoutWidget_3)
        self.xMaxLabel.setObjectName("xMaxLabel")
        self.formLayout_3.setWidget(3, QtWidgets.QFormLayout.LabelRole, self.xMaxLabel)
        self.xMaxLineEdit = QtWidgets.QLineEdit(self.formLayoutWidget_3)
        self.xMaxLineEdit.setObjectName("xMaxLineEdit")
        self.formLayout_3.setWidget(3, QtWidgets.QFormLayout.FieldRole, self.xMaxLineEdit)
        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 892, 21))
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
        self.fullNameLabel.setText(_translate("MainWindow", "Full name"))
        self.dateLabel.setText(_translate("MainWindow", "Date"))
        self.spreadingVolumeLabel.setText(_translate("MainWindow", "Surf. volume"))
        self.sensitizerLabel.setText(_translate("MainWindow", "Sensitizer"))
        self.spreadingVolumeLabel_2.setText(_translate("MainWindow", "Sens volume"))
        self.photolysisLabel.setText(_translate("MainWindow", "Photolysis"))
        self.specralRangeLabel.setText(_translate("MainWindow", "Specral range"))
        self.commentLabel.setText(_translate("MainWindow", "Comment"))
        self.surfactantLabel.setText(_translate("MainWindow", "Surfactant"))
        self.plotLabelLabel.setText(_translate("MainWindow", "Plot label"))
        self.label.setText(_translate("MainWindow", "current"))
        self.totalSpectraCountLabel.setText(_translate("MainWindow", "total spectra count"))
        self.pushButton.setText(_translate("MainWindow", "Select"))
        self.checkBox_2.setText(_translate("MainWindow", "show IR"))
        self.checkBox_3.setText(_translate("MainWindow", "show Vis"))
        self.checkBox.setText(_translate("MainWindow", "normalized"))
        self.pushButton_4.setText(_translate("MainWindow", "Refresh"))
        self.pushButton_3.setText(_translate("MainWindow", "Integrate"))
        self.xLabelLabel.setText(_translate("MainWindow", "X label"))
        self.yLabelLabel.setText(_translate("MainWindow", "Y label"))
        self.xMinLabel.setText(_translate("MainWindow", "X min,max"))
        self.xMaxLabel.setText(_translate("MainWindow", "Y min,max"))
