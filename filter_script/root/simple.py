# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'simple.ui'
#
# Created: Thu Feb  2 16:34:23 2012
#      by: PyQt4 UI code generator 4.8.4
#
# WARNING! All changes made in this file will be lost!

from PyQt4 import QtCore, QtGui

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    _fromUtf8 = lambda s: s

class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName(_fromUtf8("MainWindow"))
        MainWindow.resize(1024, 580)
        MainWindow.setProperty(_fromUtf8("QmodelIndex"), _fromUtf8(""))
        self.centralwidget = QtGui.QWidget(MainWindow)
        self.centralwidget.setFocusPolicy(QtCore.Qt.TabFocus)
        self.centralwidget.setObjectName(_fromUtf8("centralwidget"))
        self.tabWidget = QtGui.QTabWidget(self.centralwidget)
        self.tabWidget.setGeometry(QtCore.QRect(420, 60, 591, 481))
        self.tabWidget.setObjectName(_fromUtf8("tabWidget"))
        self.fEPSPanalyser = QtGui.QWidget()
        self.fEPSPanalyser.setObjectName(_fromUtf8("fEPSPanalyser"))
        self.fEPSP_button = QtGui.QPushButton(self.fEPSPanalyser)
        self.fEPSP_button.setGeometry(QtCore.QRect(440, 380, 131, 71))
        self.fEPSP_button.setObjectName(_fromUtf8("fEPSP_button"))
        self.debugBox = QtGui.QCheckBox(self.fEPSPanalyser)
        self.debugBox.setGeometry(QtCore.QRect(20, 140, 84, 20))
        self.debugBox.setObjectName(_fromUtf8("debugBox"))
        self.layoutWidget = QtGui.QWidget(self.fEPSPanalyser)
        self.layoutWidget.setGeometry(QtCore.QRect(20, 30, 191, 101))
        self.layoutWidget.setObjectName(_fromUtf8("layoutWidget"))
        self.gridLayout_2 = QtGui.QGridLayout(self.layoutWidget)
        self.gridLayout_2.setMargin(0)
        self.gridLayout_2.setObjectName(_fromUtf8("gridLayout_2"))
        self.frequency_line = QtGui.QLineEdit(self.layoutWidget)
        self.frequency_line.setObjectName(_fromUtf8("frequency_line"))
        self.gridLayout_2.addWidget(self.frequency_line, 0, 2, 1, 1)
        self.substance_line = QtGui.QLineEdit(self.layoutWidget)
        self.substance_line.setObjectName(_fromUtf8("substance_line"))
        self.gridLayout_2.addWidget(self.substance_line, 1, 2, 1, 1)
        self.label_3 = QtGui.QLabel(self.layoutWidget)
        self.label_3.setObjectName(_fromUtf8("label_3"))
        self.gridLayout_2.addWidget(self.label_3, 1, 1, 1, 1)
        self.label_2 = QtGui.QLabel(self.layoutWidget)
        self.label_2.setObjectName(_fromUtf8("label_2"))
        self.gridLayout_2.addWidget(self.label_2, 0, 1, 1, 1)
        self.database_checkBox = QtGui.QCheckBox(self.fEPSPanalyser)
        self.database_checkBox.setGeometry(QtCore.QRect(20, 170, 131, 20))
        self.database_checkBox.setChecked(True)
        self.database_checkBox.setObjectName(_fromUtf8("database_checkBox"))
        self.clusterizationBox = QtGui.QCheckBox(self.fEPSPanalyser)
        self.clusterizationBox.setGeometry(QtCore.QRect(20, 200, 241, 20))
        self.clusterizationBox.setObjectName(_fromUtf8("clusterizationBox"))
        self.tabWidget.addTab(self.fEPSPanalyser, _fromUtf8(""))
        self.tab = QtGui.QWidget()
        self.tab.setObjectName(_fromUtf8("tab"))
        self.scrollArea = QtGui.QScrollArea(self.tab)
        self.scrollArea.setGeometry(QtCore.QRect(270, 0, 311, 291))
        self.scrollArea.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAsNeeded)
        self.scrollArea.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAsNeeded)
        self.scrollArea.setObjectName(_fromUtf8("scrollArea"))
        self.scrollAreaWidgetContents = QtGui.QWidget()
        self.scrollAreaWidgetContents.setGeometry(QtCore.QRect(0, 0, 307, 287))
        self.scrollAreaWidgetContents.setObjectName(_fromUtf8("scrollAreaWidgetContents"))
        self.imageLabel = QtGui.QLabel(self.scrollAreaWidgetContents)
        self.imageLabel.setGeometry(QtCore.QRect(0, 0, 8, 16))
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.imageLabel.sizePolicy().hasHeightForWidth())
        self.imageLabel.setSizePolicy(sizePolicy)
        self.imageLabel.setFrameShape(QtGui.QFrame.Box)
        self.imageLabel.setScaledContents(True)
        self.imageLabel.setObjectName(_fromUtf8("imageLabel"))
        self.scrollArea.setWidget(self.scrollAreaWidgetContents)
        self.layoutWidget1 = QtGui.QWidget(self.tab)
        self.layoutWidget1.setGeometry(QtCore.QRect(270, 320, 301, 161))
        self.layoutWidget1.setObjectName(_fromUtf8("layoutWidget1"))
        self.gridLayout = QtGui.QGridLayout(self.layoutWidget1)
        self.gridLayout.setMargin(0)
        self.gridLayout.setObjectName(_fromUtf8("gridLayout"))
        self.createListButton = QtGui.QPushButton(self.layoutWidget1)
        self.createListButton.setObjectName(_fromUtf8("createListButton"))
        self.gridLayout.addWidget(self.createListButton, 0, 0, 1, 1)
        self.graphButton = QtGui.QPushButton(self.layoutWidget1)
        self.graphButton.setObjectName(_fromUtf8("graphButton"))
        self.gridLayout.addWidget(self.graphButton, 0, 1, 1, 1)
        self.clearList = QtGui.QPushButton(self.layoutWidget1)
        self.clearList.setObjectName(_fromUtf8("clearList"))
        self.gridLayout.addWidget(self.clearList, 1, 0, 1, 1)
        self.rmItemButton = QtGui.QPushButton(self.layoutWidget1)
        self.rmItemButton.setObjectName(_fromUtf8("rmItemButton"))
        self.gridLayout.addWidget(self.rmItemButton, 1, 1, 1, 1)
        self.processingLabel = QtGui.QLabel(self.tab)
        self.processingLabel.setGeometry(QtCore.QRect(0, 190, 91, 16))
        self.processingLabel.setObjectName(_fromUtf8("processingLabel"))
        self.processedList = QtGui.QListWidget(self.tab)
        self.processedList.setGeometry(QtCore.QRect(0, 210, 256, 291))
        self.processedList.setAcceptDrops(True)
        self.processedList.setDragEnabled(True)
        self.processedList.setDragDropMode(QtGui.QAbstractItemView.InternalMove)
        self.processedList.setSelectionMode(QtGui.QAbstractItemView.ExtendedSelection)
        self.processedList.setObjectName(_fromUtf8("processedList"))
        self.layoutWidget2 = QtGui.QWidget(self.tab)
        self.layoutWidget2.setGeometry(QtCore.QRect(20, 10, 221, 171))
        self.layoutWidget2.setObjectName(_fromUtf8("layoutWidget2"))
        self.gridLayout_3 = QtGui.QGridLayout(self.layoutWidget2)
        self.gridLayout_3.setMargin(0)
        self.gridLayout_3.setObjectName(_fromUtf8("gridLayout_3"))
        self.label = QtGui.QLabel(self.layoutWidget2)
        self.label.setObjectName(_fromUtf8("label"))
        self.gridLayout_3.addWidget(self.label, 0, 0, 1, 1)
        self.startLine = QtGui.QLineEdit(self.layoutWidget2)
        self.startLine.setObjectName(_fromUtf8("startLine"))
        self.gridLayout_3.addWidget(self.startLine, 0, 1, 1, 1)
        self.label_4 = QtGui.QLabel(self.layoutWidget2)
        self.label_4.setObjectName(_fromUtf8("label_4"))
        self.gridLayout_3.addWidget(self.label_4, 1, 0, 1, 1)
        self.stopLine = QtGui.QLineEdit(self.layoutWidget2)
        self.stopLine.setObjectName(_fromUtf8("stopLine"))
        self.gridLayout_3.addWidget(self.stopLine, 1, 1, 1, 1)
        self.label_5 = QtGui.QLabel(self.layoutWidget2)
        self.label_5.setObjectName(_fromUtf8("label_5"))
        self.gridLayout_3.addWidget(self.label_5, 2, 0, 1, 1)
        self.rstrideLine = QtGui.QLineEdit(self.layoutWidget2)
        self.rstrideLine.setObjectName(_fromUtf8("rstrideLine"))
        self.gridLayout_3.addWidget(self.rstrideLine, 2, 1, 1, 1)
        self.label_6 = QtGui.QLabel(self.layoutWidget2)
        self.label_6.setObjectName(_fromUtf8("label_6"))
        self.gridLayout_3.addWidget(self.label_6, 3, 0, 1, 1)
        self.cstrideLine = QtGui.QLineEdit(self.layoutWidget2)
        self.cstrideLine.setObjectName(_fromUtf8("cstrideLine"))
        self.gridLayout_3.addWidget(self.cstrideLine, 3, 1, 1, 1)
        self.debugBox2 = QtGui.QCheckBox(self.layoutWidget2)
        self.debugBox2.setObjectName(_fromUtf8("debugBox2"))
        self.gridLayout_3.addWidget(self.debugBox2, 4, 0, 1, 2)
        self.tabWidget.addTab(self.tab, _fromUtf8(""))
        self.tab_2 = QtGui.QWidget()
        self.tab_2.setObjectName(_fromUtf8("tab_2"))
        self.frame = QtGui.QFrame(self.tab_2)
        self.frame.setGeometry(QtCore.QRect(20, 10, 291, 181))
        self.frame.setFrameShape(QtGui.QFrame.Panel)
        self.frame.setFrameShadow(QtGui.QFrame.Raised)
        self.frame.setObjectName(_fromUtf8("frame"))
        self.clearDbButton = QtGui.QPushButton(self.frame)
        self.clearDbButton.setGeometry(QtCore.QRect(10, 150, 121, 25))
        self.clearDbButton.setObjectName(_fromUtf8("clearDbButton"))
        self.saveDbButton = QtGui.QPushButton(self.frame)
        self.saveDbButton.setGeometry(QtCore.QRect(160, 150, 121, 25))
        self.saveDbButton.setObjectName(_fromUtf8("saveDbButton"))
        self.layoutWidget3 = QtGui.QWidget(self.tab_2)
        self.layoutWidget3.setGeometry(QtCore.QRect(30, 20, 271, 137))
        self.layoutWidget3.setObjectName(_fromUtf8("layoutWidget3"))
        self.gridLayout_4 = QtGui.QGridLayout(self.layoutWidget3)
        self.gridLayout_4.setMargin(0)
        self.gridLayout_4.setObjectName(_fromUtf8("gridLayout_4"))
        self.label_7 = QtGui.QLabel(self.layoutWidget3)
        self.label_7.setObjectName(_fromUtf8("label_7"))
        self.gridLayout_4.addWidget(self.label_7, 2, 0, 1, 1)
        self.dbNameLine = QtGui.QLineEdit(self.layoutWidget3)
        self.dbNameLine.setObjectName(_fromUtf8("dbNameLine"))
        self.gridLayout_4.addWidget(self.dbNameLine, 2, 1, 1, 1)
        self.label_8 = QtGui.QLabel(self.layoutWidget3)
        self.label_8.setObjectName(_fromUtf8("label_8"))
        self.gridLayout_4.addWidget(self.label_8, 3, 0, 1, 1)
        self.dbUserLine = QtGui.QLineEdit(self.layoutWidget3)
        self.dbUserLine.setObjectName(_fromUtf8("dbUserLine"))
        self.gridLayout_4.addWidget(self.dbUserLine, 3, 1, 1, 1)
        self.label_9 = QtGui.QLabel(self.layoutWidget3)
        self.label_9.setObjectName(_fromUtf8("label_9"))
        self.gridLayout_4.addWidget(self.label_9, 4, 0, 1, 1)
        self.dbPassLine = QtGui.QLineEdit(self.layoutWidget3)
        self.dbPassLine.setObjectName(_fromUtf8("dbPassLine"))
        self.gridLayout_4.addWidget(self.dbPassLine, 4, 1, 1, 1)
        self.label_10 = QtGui.QLabel(self.layoutWidget3)
        self.label_10.setObjectName(_fromUtf8("label_10"))
        self.gridLayout_4.addWidget(self.label_10, 0, 0, 1, 2)
        self.dbServerIpLine = QtGui.QLineEdit(self.layoutWidget3)
        self.dbServerIpLine.setObjectName(_fromUtf8("dbServerIpLine"))
        self.gridLayout_4.addWidget(self.dbServerIpLine, 1, 1, 1, 1)
        self.label_11 = QtGui.QLabel(self.layoutWidget3)
        self.label_11.setObjectName(_fromUtf8("label_11"))
        self.gridLayout_4.addWidget(self.label_11, 1, 0, 1, 1)
        self.tabWidget.addTab(self.tab_2, _fromUtf8(""))
        self.sourceList = QtGui.QListView(self.centralwidget)
        self.sourceList.setGeometry(QtCore.QRect(20, 60, 351, 481))
        self.sourceList.setAcceptDrops(True)
        self.sourceList.setDragEnabled(True)
        self.sourceList.setDragDropOverwriteMode(True)
        self.sourceList.setDragDropMode(QtGui.QAbstractItemView.DragDrop)
        self.sourceList.setDefaultDropAction(QtCore.Qt.CopyAction)
        self.sourceList.setSelectionMode(QtGui.QAbstractItemView.ExtendedSelection)
        self.sourceList.setObjectName(_fromUtf8("sourceList"))
        self.pathLine = QtGui.QLineEdit(self.centralwidget)
        self.pathLine.setGeometry(QtCore.QRect(20, 30, 941, 25))
        self.pathLine.setMaxLength(32765)
        self.pathLine.setDragEnabled(True)
        self.pathLine.setObjectName(_fromUtf8("pathLine"))
        self.pathLabel = QtGui.QLabel(self.centralwidget)
        self.pathLabel.setGeometry(QtCore.QRect(20, 10, 138, 16))
        self.pathLabel.setObjectName(_fromUtf8("pathLabel"))
        self.exitButton = QtGui.QPushButton(self.centralwidget)
        self.exitButton.setGeometry(QtCore.QRect(860, 550, 146, 25))
        self.exitButton.setObjectName(_fromUtf8("exitButton"))
        self.upButton = QtGui.QPushButton(self.centralwidget)
        self.upButton.setGeometry(QtCore.QRect(970, 30, 31, 25))
        self.upButton.setObjectName(_fromUtf8("upButton"))
        MainWindow.setCentralWidget(self.centralwidget)
        self.label_3.setBuddy(self.substance_line)
        self.label_2.setBuddy(self.frequency_line)
        self.label.setBuddy(self.startLine)
        self.label_4.setBuddy(self.stopLine)
        self.label_5.setBuddy(self.rstrideLine)
        self.label_6.setBuddy(self.cstrideLine)
        self.label_7.setBuddy(self.dbNameLine)
        self.label_8.setBuddy(self.dbUserLine)
        self.label_9.setBuddy(self.dbPassLine)
        self.label_11.setBuddy(self.dbServerIpLine)
        self.pathLabel.setBuddy(self.pathLine)

        self.retranslateUi(MainWindow)
        self.tabWidget.setCurrentIndex(0)
        QtCore.QObject.connect(self.exitButton, QtCore.SIGNAL(_fromUtf8("clicked()")), MainWindow.close)
        QtCore.QObject.connect(self.clearList, QtCore.SIGNAL(_fromUtf8("clicked()")), self.processedList.clear)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        MainWindow.setWindowTitle(QtGui.QApplication.translate("MainWindow", "MainWindow", None, QtGui.QApplication.UnicodeUTF8))
        self.fEPSP_button.setText(QtGui.QApplication.translate("MainWindow", "Start fEPSP-analiser", None, QtGui.QApplication.UnicodeUTF8))
        self.debugBox.setText(QtGui.QApplication.translate("MainWindow", "debug", None, QtGui.QApplication.UnicodeUTF8))
        self.frequency_line.setText(QtGui.QApplication.translate("MainWindow", "200000", None, QtGui.QApplication.UnicodeUTF8))
        self.label_3.setText(QtGui.QApplication.translate("MainWindow", "substance", None, QtGui.QApplication.UnicodeUTF8))
        self.label_2.setText(QtGui.QApplication.translate("MainWindow", "frequency", None, QtGui.QApplication.UnicodeUTF8))
        self.database_checkBox.setText(QtGui.QApplication.translate("MainWindow", "write to database", None, QtGui.QApplication.UnicodeUTF8))
        self.clusterizationBox.setText(QtGui.QApplication.translate("MainWindow", "replace clusterization by stimDetect", None, QtGui.QApplication.UnicodeUTF8))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.fEPSPanalyser), QtGui.QApplication.translate("MainWindow", "fEPSP-analyser", None, QtGui.QApplication.UnicodeUTF8))
        self.createListButton.setText(QtGui.QApplication.translate("MainWindow", "создать список", None, QtGui.QApplication.UnicodeUTF8))
        self.graphButton.setText(QtGui.QApplication.translate("MainWindow", "graph", None, QtGui.QApplication.UnicodeUTF8))
        self.clearList.setText(QtGui.QApplication.translate("MainWindow", "Очистить список", None, QtGui.QApplication.UnicodeUTF8))
        self.rmItemButton.setText(QtGui.QApplication.translate("MainWindow", "remove item", None, QtGui.QApplication.UnicodeUTF8))
        self.processingLabel.setText(QtGui.QApplication.translate("MainWindow", "В обработку:", None, QtGui.QApplication.UnicodeUTF8))
        self.label.setText(QtGui.QApplication.translate("MainWindow", "start from", None, QtGui.QApplication.UnicodeUTF8))
        self.startLine.setText(QtGui.QApplication.translate("MainWindow", "0", None, QtGui.QApplication.UnicodeUTF8))
        self.label_4.setText(QtGui.QApplication.translate("MainWindow", "stop at", None, QtGui.QApplication.UnicodeUTF8))
        self.stopLine.setText(QtGui.QApplication.translate("MainWindow", "8000", None, QtGui.QApplication.UnicodeUTF8))
        self.label_5.setText(QtGui.QApplication.translate("MainWindow", "rstride", None, QtGui.QApplication.UnicodeUTF8))
        self.rstrideLine.setText(QtGui.QApplication.translate("MainWindow", "2", None, QtGui.QApplication.UnicodeUTF8))
        self.label_6.setText(QtGui.QApplication.translate("MainWindow", "cstride", None, QtGui.QApplication.UnicodeUTF8))
        self.cstrideLine.setText(QtGui.QApplication.translate("MainWindow", "20", None, QtGui.QApplication.UnicodeUTF8))
        self.debugBox2.setText(QtGui.QApplication.translate("MainWindow", "debug", None, QtGui.QApplication.UnicodeUTF8))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab), QtGui.QApplication.translate("MainWindow", "graph", None, QtGui.QApplication.UnicodeUTF8))
        self.clearDbButton.setText(QtGui.QApplication.translate("MainWindow", "Clear db config", None, QtGui.QApplication.UnicodeUTF8))
        self.saveDbButton.setText(QtGui.QApplication.translate("MainWindow", "Save db config", None, QtGui.QApplication.UnicodeUTF8))
        self.label_7.setText(QtGui.QApplication.translate("MainWindow", "Database Name", None, QtGui.QApplication.UnicodeUTF8))
        self.label_8.setText(QtGui.QApplication.translate("MainWindow", "User Name", None, QtGui.QApplication.UnicodeUTF8))
        self.label_9.setText(QtGui.QApplication.translate("MainWindow", "Password", None, QtGui.QApplication.UnicodeUTF8))
        self.label_10.setText(QtGui.QApplication.translate("MainWindow", "DataBase configuration", None, QtGui.QApplication.UnicodeUTF8))
        self.label_11.setText(QtGui.QApplication.translate("MainWindow", "Server address", None, QtGui.QApplication.UnicodeUTF8))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab_2), QtGui.QApplication.translate("MainWindow", "config", None, QtGui.QApplication.UnicodeUTF8))
        self.pathLine.setText(QtGui.QApplication.translate("MainWindow", "data/", None, QtGui.QApplication.UnicodeUTF8))
        self.pathLabel.setText(QtGui.QApplication.translate("MainWindow", "Путь до директории:", None, QtGui.QApplication.UnicodeUTF8))
        self.exitButton.setText(QtGui.QApplication.translate("MainWindow", "выход", None, QtGui.QApplication.UnicodeUTF8))
        self.upButton.setText(QtGui.QApplication.translate("MainWindow", "../", None, QtGui.QApplication.UnicodeUTF8))

