# -*- coding: utf-8 -*-
"""
Created on Sat May 25 08:19:13 2019

@author: KTMS458
"""

from PyQt5 import QtCore, QtWidgets
from PyQt5.QtWidgets import QDialog, QFileDialog,QColorDialog, QMessageBox

import os.path as path
import xml.dom.minidom as xml
from PyQt5.QtWidgets import QWidget
import seaborn as sns



class ProjectWizard(QDialog):
    def __init__(self,input_table):
        super(QWidget,self).__init__()
        self.DataFrame=input_table
        self.setupUi()
        self.configFile=""
        self.projectDirectory=""
#        self.exec_()

    def closeEvent(self,event):
        valid=self.validateInput()        
        if not valid:
            event.ignore()
            return
        output=self.createWarnBox('Are you sure to quit?\nHave you saved the configuration?')
        if not output: 
            event.ignore()
            return
        self.updateDataFrame()
        
    
    def setupUi(self):
        self._want_to_close = True
        self.DataFrame['basename']=[path.basename(x) for x in self.DataFrame['Location']]
        self.DataFrame['label']=[item.split(".")[0] for item in self.DataFrame['basename']]
        
        self.setObjectName("Dialog")
        self.resize(845, 553)
        self.setMinimumSize(QtCore.QSize(845, 553))
        self.setMaximumSize(QtCore.QSize(845, 553))
        self.centralwidget = QtWidgets.QWidget(self)
        self.centralwidget.setObjectName("centralwidget")
        
        
        self.gridLayout = QtWidgets.QGridLayout(self.centralwidget)
        self.gridLayout.setContentsMargins(-1, -1, 15, 2)
        self.gridLayout.setObjectName("gridLayout")

        self.saveButton=QtWidgets.QPushButton()
        self.saveButton.setMinimumSize(QtCore.QSize(100, 30))
        self.saveButton.setMaximumSize(QtCore.QSize(100, 30))
        self.saveAsButton=QtWidgets.QPushButton()
        self.saveAsButton.setMinimumSize(QtCore.QSize(100, 30))
        self.saveAsButton.setMaximumSize(QtCore.QSize(100, 30))
        self.loadButton=QtWidgets.QPushButton()
        self.loadButton.setMinimumSize(QtCore.QSize(100, 30))
        self.loadButton.setMaximumSize(QtCore.QSize(100, 30))
        self.okButton=QtWidgets.QPushButton()
        self.okButton.setMinimumSize(QtCore.QSize(100, 30))
        self.okButton.setMaximumSize(QtCore.QSize(100, 30))

        
        self.saveButton.clicked.connect(self.saveConfigFile)
        self.saveAsButton.clicked.connect(self.saveAsConfigFile)
        self.loadButton.clicked.connect(self.loadConfigFile)
        #self.okButton.clicked.connect(self.closeEvent)

        self.horizontalLayoutButton = QtWidgets.QHBoxLayout()
        self.horizontalLayoutButton.setObjectName("horizontalLayoutButton")
        
        spacerItem = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayoutButton.addItem(spacerItem)  
        self.horizontalLayoutButton.addWidget(self.saveButton)
        self.horizontalLayoutButton.addWidget(self.saveAsButton)
        self.horizontalLayoutButton.addWidget(self.loadButton)
        self.horizontalLayoutButton.addWidget(self.okButton)
        self.gridLayout.addLayout(self.horizontalLayoutButton,2,0,1,1)
        self.tabWidget = QtWidgets.QTabWidget(self)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.tabWidget.sizePolicy().hasHeightForWidth())
        self.tabWidget.setSizePolicy(sizePolicy)
        self.tabWidget.setMinimumSize(QtCore.QSize(821, 501))
        self.tabWidget.setMaximumSize(QtCore.QSize(821, 501))
        self.tabWidget.setStyleSheet("")
        self.tabWidget.setTabPosition(QtWidgets.QTabWidget.North)
        self.tabWidget.setTabShape(QtWidgets.QTabWidget.Triangular)
        self.tabWidget.setObjectName("tabWidget")
        self.ProjectTab = QtWidgets.QWidget()
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.ProjectTab.sizePolicy().hasHeightForWidth())
        self.ProjectTab.setSizePolicy(sizePolicy)
        self.ProjectTab.setStyleSheet("")
        self.ProjectTab.setObjectName("ProjectTab")
        


        self.gbNewProject = QtWidgets.QGroupBox(self.ProjectTab)
        self.gbNewProject.setGeometry(QtCore.QRect(10, 20, 791, 151))
        self.gbNewProject.setMinimumSize(QtCore.QSize(750, 151))
        self.gbNewProject.setMaximumSize(QtCore.QSize(812, 151))
        self.gbNewProject.setStyleSheet("border-color: rgb(0, 0, 0);")
        self.gbNewProject.setAlignment(QtCore.Qt.AlignLeading|QtCore.Qt.AlignLeft|QtCore.Qt.AlignVCenter)
        self.gbNewProject.setObjectName("gbNewProject")
        self.verticalLayout_4 = QtWidgets.QVBoxLayout(self.gbNewProject)
        self.verticalLayout_4.setObjectName("verticalLayout_4")
        self.CreateNewProjectGridlayout = QtWidgets.QGridLayout()
        self.CreateNewProjectGridlayout.setHorizontalSpacing(20)
        self.CreateNewProjectGridlayout.setVerticalSpacing(10)
        self.CreateNewProjectGridlayout.setObjectName("CreateNewProjectGridlayout")
        self.lbDLIMSProject = QtWidgets.QLabel(self.gbNewProject)
        self.lbDLIMSProject.setObjectName("lbDLIMSProject")
        self.CreateNewProjectGridlayout.addWidget(self.lbDLIMSProject, 1, 3, 1, 1)
        self.leDLIMSProject = QtWidgets.QLineEdit(self.gbNewProject)
        self.leDLIMSProject.setObjectName("leDLIMSProject")
        self.CreateNewProjectGridlayout.addWidget(self.leDLIMSProject, 1, 4, 1, 1)
        self.leProjectCode = QtWidgets.QLineEdit(self.gbNewProject)
        self.leProjectCode.setObjectName("leProjectCode")
        self.CreateNewProjectGridlayout.addWidget(self.leProjectCode, 0, 4, 1, 1)
        self.lbProjectCode = QtWidgets.QLabel(self.gbNewProject)
        self.lbProjectCode.setObjectName("lbProjectCode")
        self.CreateNewProjectGridlayout.addWidget(self.lbProjectCode, 0, 3, 1, 1)
        self.lbProjectFolder = QtWidgets.QLabel(self.gbNewProject)
        self.lbProjectFolder.setObjectName("lbProjectFolder")
        self.CreateNewProjectGridlayout.addWidget(self.lbProjectFolder, 2, 0, 1, 1)
        self.lbAnalyst = QtWidgets.QLabel(self.gbNewProject)
        self.lbAnalyst.setObjectName("lbAnalyst")
        self.CreateNewProjectGridlayout.addWidget(self.lbAnalyst, 1, 0, 1, 1)
        self.lbProjectName = QtWidgets.QLabel(self.gbNewProject)
        self.lbProjectName.setObjectName("lbProjectName")
        self.CreateNewProjectGridlayout.addWidget(self.lbProjectName, 0, 0, 1, 1)
        self.leProjectFolder = QtWidgets.QLineEdit(self.gbNewProject)
        self.leProjectFolder.setObjectName("leProjectFolder")
        self.CreateNewProjectGridlayout.addWidget(self.leProjectFolder, 2, 1, 1, 1)
        self.toolButton = QtWidgets.QToolButton(self.gbNewProject)
        self.toolButton.clicked.connect(self.selectDir)
        self.toolButton.setObjectName("toolButton")
        self.CreateNewProjectGridlayout.addWidget(self.toolButton, 2, 2, 1, 1)
        self.leAnalyst = QtWidgets.QLineEdit(self.gbNewProject)
        self.leAnalyst.setObjectName("leAnalyst")
        self.CreateNewProjectGridlayout.addWidget(self.leAnalyst, 1, 1, 1, 2)
        self.leProjectName = QtWidgets.QLineEdit(self.gbNewProject)
        self.leProjectName.setObjectName("leProjectName")
        self.CreateNewProjectGridlayout.addWidget(self.leProjectName, 0, 1, 1, 2)
        self.verticalLayout_4.addLayout(self.CreateNewProjectGridlayout)
        self.tabWidget.addTab(self.ProjectTab, "")
        self.TabFile = QtWidgets.QWidget()
        self.TabFile.setStyleSheet("")
        self.TabFile.setObjectName("TabFile")
        
        self.ProjectTab.setTabOrder(self.leProjectName,self.leAnalyst)
        self.ProjectTab.setTabOrder(self.leAnalyst,self.leProjectFolder)
        self.ProjectTab.setTabOrder(self.leProjectFolder,self.leProjectCode)
        self.ProjectTab.setTabOrder(self.leProjectCode,self.leDLIMSProject)
        
        self.widget = QtWidgets.QWidget(self.TabFile)
        self.widget.setGeometry(QtCore.QRect(2, 12, 811, 461))
        self.widget.setObjectName("widget")
        self.gridLayout_2 = QtWidgets.QGridLayout(self.widget)
        self.gridLayout_2.setContentsMargins(0, 0, 0, 0)
        self.gridLayout_2.setHorizontalSpacing(10)
        self.gridLayout_2.setObjectName("gridLayout_2")
        self.lbReferenceFile = QtWidgets.QLabel(self.widget)
        self.lbReferenceFile.setObjectName("lbReferenceFile")
        self.gridLayout_2.addWidget(self.lbReferenceFile, 0, 0, 1, 1)
        self.cbReferenceFile = QtWidgets.QComboBox(self.widget)
        self.cbReferenceFile.addItems(self.DataFrame['label'])
        
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Maximum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.cbReferenceFile.sizePolicy().hasHeightForWidth())
        self.cbReferenceFile.setSizePolicy(sizePolicy)
        self.cbReferenceFile.setSizeAdjustPolicy(QtWidgets.QComboBox.AdjustToMinimumContentsLength)
        self.cbReferenceFile.setObjectName("cbReferenceFile")
        self.gridLayout_2.addWidget(self.cbReferenceFile, 0, 1, 1, 1)
        self.frame_2 = QtWidgets.QFrame(self.widget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.frame_2.sizePolicy().hasHeightForWidth())
        self.frame_2.setSizePolicy(sizePolicy)
        self.frame_2.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame_2.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame_2.setObjectName("frame_2")
        self.gridLayout_2.addWidget(self.frame_2, 0, 2, 1, 1)
        

        
        self.tableWidget = QtWidgets.QTableWidget(self.widget)
        self.tableWidget.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAsNeeded)
        self.tableWidget.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAsNeeded)
        self.tableWidget.setSizeAdjustPolicy(QtWidgets.QAbstractScrollArea.AdjustToContentsOnFirstShow)
        self.tableWidget.setObjectName("tableWidget")
        self.tableWidget.setColumnCount(5)
        self.tableWidget.setRowCount(len(self.DataFrame))
        
        ls_order=[str(i+1) for i in range(len(self.DataFrame)) ]
        color_pallete=sns.color_palette("muted",len(ls_order))
        color_pallete=color_pallete.as_hex()
        for n,item in enumerate(self.DataFrame['label']):
            color=color_pallete[n]
            self.tableWidget.setItem(n, 0, QtWidgets.QTableWidgetItem(str(item)))
            comboBox = QtWidgets.QComboBox()
            comboBox.addItems(['Sample','Blank','System Suit'])
            self.tableWidget.setCellWidget(n,1,comboBox)
            comboBox = QtWidgets.QComboBox()
            comboBox.addItems(ls_order)
            self.tableWidget.setCellWidget(n,2,comboBox)
            colorButton=QtWidgets.QPushButton()
            colorButton.clicked.connect(self.getColor)
            colorButton.setStyleSheet("background-color: %s"%(color))
            colorButton.setProperty("row",n)
            self.tableWidget.setCellWidget(n,3,colorButton)
            self.tableWidget.setItem(n, 4, QtWidgets.QTableWidgetItem("%s"%(color)))
        
        self.tableWidget.itemChanged.connect(self.adjustColor)

                                                    
        
        self.tableWidget.setVerticalHeaderLabels(self.DataFrame['basename'])
        item = QtWidgets.QTableWidgetItem()
        self.tableWidget.setHorizontalHeaderItem(0, item)
        item = QtWidgets.QTableWidgetItem()
        self.tableWidget.setHorizontalHeaderItem(1, item)
        item = QtWidgets.QTableWidgetItem()
        self.tableWidget.setHorizontalHeaderItem(2, item)
        item = QtWidgets.QTableWidgetItem()
        self.tableWidget.setHorizontalHeaderItem(3, item)
        item = QtWidgets.QTableWidgetItem()
        self.tableWidget.setHorizontalHeaderItem(4, item)
        self.tableWidget.horizontalHeader().setCascadingSectionResizes(True)
        self.tableWidget.horizontalHeader().setDefaultSectionSize(170)
        self.tableWidget.horizontalHeader().setStretchLastSection(True)
        self.tableWidget.verticalHeader().setCascadingSectionResizes(True)
        self.tableWidget.verticalHeader().setSortIndicatorShown(False)
        self.tableWidget.verticalHeader().setStretchLastSection(False)        
        
        
        
        self.gridLayout_2.addWidget(self.tableWidget, 1, 0, 1, 3)
        self.tabWidget.addTab(self.TabFile, "")
        self.XICTab = QtWidgets.QWidget()
        self.XICTab.setObjectName("XICTab")
        self.gbSmoothXIC = QtWidgets.QGroupBox(self.XICTab)
        self.gbSmoothXIC.setGeometry(QtCore.QRect(10, 150, 791, 141))
        self.gbSmoothXIC.setObjectName("gbSmoothXIC")
        self.cbApplySmoothing = QtWidgets.QCheckBox(self.gbSmoothXIC)
        self.cbApplySmoothing.setGeometry(QtCore.QRect(20, 20, 151, 17))
        self.cbApplySmoothing.setObjectName("cbApplySmoothing")
        self.cbApplySmoothing.setChecked(True)
        self.cbApplySmoothing.stateChanged.connect(self.changeSmooth)
        
        
        self.layoutWidget = QtWidgets.QWidget(self.gbSmoothXIC)
        self.layoutWidget.setGeometry(QtCore.QRect(30, 40, 231, 81))
        self.layoutWidget.setObjectName("layoutWidget")
        self.smoothGridlayout = QtWidgets.QGridLayout(self.layoutWidget)
        self.smoothGridlayout.setSizeConstraint(QtWidgets.QLayout.SetDefaultConstraint)
        self.smoothGridlayout.setContentsMargins(0, 0, 0, 0)
        self.smoothGridlayout.setObjectName("smoothGridlayout")
        self.sbCycles = QtWidgets.QSpinBox(self.layoutWidget)
        self.sbCycles.setObjectName("sbCycles")
        self.smoothGridlayout.addWidget(self.sbCycles, 1, 1, 1, 1)
        self.cbSmoothType = QtWidgets.QComboBox(self.layoutWidget)
        self.cbSmoothType.setEditable(False)
        self.cbSmoothType.setDuplicatesEnabled(False)
        self.cbSmoothType.setObjectName("cbSmoothType")
        self.cbSmoothType.addItems(['Move average','Gaussian','Savitzky Golay'])
        
        
        self.smoothGridlayout.addWidget(self.cbSmoothType, 0, 1, 1, 1)
        self.dsbSmoothWindow = QtWidgets.QDoubleSpinBox(self.layoutWidget)
        self.dsbSmoothWindow.setDecimals(1)
        self.dsbSmoothWindow.setObjectName("dsbSmoothWindow")
        self.smoothGridlayout.addWidget(self.dsbSmoothWindow, 2, 1, 1, 1)
        self.lbSelectSmoothType = QtWidgets.QLabel(self.layoutWidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.lbSelectSmoothType.sizePolicy().hasHeightForWidth())
        self.lbSelectSmoothType.setSizePolicy(sizePolicy)
        self.lbSelectSmoothType.setAlignment(QtCore.Qt.AlignLeading|QtCore.Qt.AlignLeft|QtCore.Qt.AlignVCenter)
        self.lbSelectSmoothType.setObjectName("lbSelectSmoothType")
        self.smoothGridlayout.addWidget(self.lbSelectSmoothType, 0, 0, 1, 1)
        self.lbCycles = QtWidgets.QLabel(self.layoutWidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.lbCycles.sizePolicy().hasHeightForWidth())
        self.lbCycles.setSizePolicy(sizePolicy)
        self.lbCycles.setObjectName("lbCycles")
        self.smoothGridlayout.addWidget(self.lbCycles, 1, 0, 1, 1)
        self.lbWindow = QtWidgets.QLabel(self.layoutWidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.lbWindow.sizePolicy().hasHeightForWidth())
        self.lbWindow.setSizePolicy(sizePolicy)
        self.lbWindow.setObjectName("lbWindow")
        self.smoothGridlayout.addWidget(self.lbWindow, 2, 0, 1, 1)
        self.gbGeneralXIC = QtWidgets.QGroupBox(self.XICTab)
        self.gbGeneralXIC.setGeometry(QtCore.QRect(10, 20, 791, 111))
        self.gbGeneralXIC.setMinimumSize(QtCore.QSize(750, 111))
        self.gbGeneralXIC.setObjectName("gbGeneralXIC")
        self.layoutWidget_2 = QtWidgets.QWidget(self.gbGeneralXIC)
        self.layoutWidget_2.setGeometry(QtCore.QRect(30, 20, 471, 80))
        self.layoutWidget_2.setObjectName("layoutWidget_2")
        self.XICGeneralGridlayout = QtWidgets.QGridLayout(self.layoutWidget_2)
        self.XICGeneralGridlayout.setContentsMargins(0, 0, 0, 0)
        self.XICGeneralGridlayout.setObjectName("XICGeneralGridlayout")
        self.dsbTimeWindow = QtWidgets.QDoubleSpinBox(self.layoutWidget_2)
        self.dsbTimeWindow.setDecimals(1)
        self.dsbTimeWindow.setObjectName("dsbTimeWindow")
        self.XICGeneralGridlayout.addWidget(self.dsbTimeWindow, 1, 1, 1, 1)
        self.lbCreateBasedOn = QtWidgets.QLabel(self.layoutWidget_2)
        self.lbCreateBasedOn.setObjectName("lbCreateBasedOn")
        self.XICGeneralGridlayout.addWidget(self.lbCreateBasedOn, 0, 0, 1, 1)

        self.cbCreateBasedOn = QtWidgets.QComboBox(self.layoutWidget_2)
        self.cbCreateBasedOn.setObjectName("cbCreateBasedOn")
        self.cbCreateBasedOn.addItems(['monoisotopic peak','highest intense peak'])
        self.XICGeneralGridlayout.addWidget(self.cbCreateBasedOn, 0, 1, 1, 1)
        self.lbXICTimeWindow = QtWidgets.QLabel(self.layoutWidget_2)
        self.lbXICTimeWindow.setObjectName("lbXICTimeWindow")
        self.XICGeneralGridlayout.addWidget(self.lbXICTimeWindow, 1, 0, 1, 1)
        self.frame = QtWidgets.QFrame(self.layoutWidget_2)
        self.frame.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame.setObjectName("frame")
        self.XICGeneralGridlayout.addWidget(self.frame, 0, 2, 1, 1)
        self.tabWidget.addTab(self.XICTab, "")
        
        self.FullMSTab = QtWidgets.QWidget()
        self.FullMSTab.setObjectName("FullMSTab")
        self.gbFullMSGeneral_3 = QtWidgets.QGroupBox(self.FullMSTab)
        self.gbFullMSGeneral_3.setGeometry(QtCore.QRect(10, 100, 791, 80))
        self.gbFullMSGeneral_3.setObjectName("gbFullMSGeneral_3")
        self.lbErrorTolerance_2 = QtWidgets.QLabel(self.gbFullMSGeneral_3)
        self.lbErrorTolerance_2.setGeometry(QtCore.QRect(30, 20, 111, 20))
        self.lbErrorTolerance_2.setObjectName("lbErrorTolerance_2")
        self.leErrorTolerance_2 = QtWidgets.QSpinBox(self.gbFullMSGeneral_3)
        self.leErrorTolerance_2.setGeometry(QtCore.QRect(160, 20, 201, 20))
        self.leErrorTolerance_2.setObjectName("leErrorTolerance_2")
        self.tabWidget.addTab(self.FullMSTab, "")
        
        self.MSMSTab = QtWidgets.QWidget()
        self.MSMSTab.setObjectName("MSMSTab")
        self.gbGeneralMSMS = QtWidgets.QGroupBox(self.MSMSTab)
        self.gbGeneralMSMS.setGeometry(QtCore.QRect(10, 20, 791, 361))
        self.gbGeneralMSMS.setObjectName("gbGeneralMSMS")
        self.gridLayoutWidget = QtWidgets.QWidget(self.gbGeneralMSMS)
        self.gridLayoutWidget.setGeometry(QtCore.QRect(10, 260, 282, 81))
        self.gridLayoutWidget.setObjectName("gridLayoutWidget")
        self.gridLayout_5 = QtWidgets.QGridLayout(self.gridLayoutWidget)
        self.gridLayout_5.setSizeConstraint(QtWidgets.QLayout.SetMinimumSize)
        self.gridLayout_5.setContentsMargins(0, 0, 0, 0)
        self.gridLayout_5.setObjectName("gridLayout_5")
        self.cbFragmentLosses = QtWidgets.QCheckBox(self.gridLayoutWidget)
        self.cbFragmentLosses.setObjectName("cbFragmentLosses")
        self.gridLayout_5.addWidget(self.cbFragmentLosses, 0, 2, 1, 1)
        self.cbIonX = QtWidgets.QCheckBox(self.gridLayoutWidget)
        self.cbIonX.setObjectName("cbIonX")
        self.gridLayout_5.addWidget(self.cbIonX, 0, 1, 1, 1)
        self.cbIonY = QtWidgets.QCheckBox(self.gridLayoutWidget)
        self.cbIonY.setObjectName("cbIonY")
        self.gridLayout_5.addWidget(self.cbIonY, 1, 1, 1, 1)
        self.cbIonZ = QtWidgets.QCheckBox(self.gridLayoutWidget)
        self.cbIonZ.setObjectName("cbIonZ")
        self.gridLayout_5.addWidget(self.cbIonZ, 2, 1, 1, 1)
        self.cbIonC = QtWidgets.QCheckBox(self.gridLayoutWidget)
        self.cbIonC.setObjectName("cbIonC")
        self.gridLayout_5.addWidget(self.cbIonC, 2, 0, 1, 1)
        self.cbIonB = QtWidgets.QCheckBox(self.gridLayoutWidget)
        self.cbIonB.setObjectName("cbIonB")
        self.gridLayout_5.addWidget(self.cbIonB, 1, 0, 1, 1)
        self.cbIonA = QtWidgets.QCheckBox(self.gridLayoutWidget)
        self.cbIonA.setObjectName("cbIonA")
        self.gridLayout_5.addWidget(self.cbIonA, 0, 0, 1, 1)
        self.lbIonSeries = QtWidgets.QLabel(self.gbGeneralMSMS)
        self.lbIonSeries.setGeometry(QtCore.QRect(10, 240, 138, 13))
        self.lbIonSeries.setObjectName("lbIonSeries")
        self.layoutWidget_6 = QtWidgets.QWidget(self.gbGeneralMSMS)
        self.layoutWidget_6.setGeometry(QtCore.QRect(10, 130, 281, 74))
        self.layoutWidget_6.setObjectName("layoutWidget_6")
        self.gridLayout_7 = QtWidgets.QGridLayout(self.layoutWidget_6)
        self.gridLayout_7.setContentsMargins(0, 0, 0, 0)
        self.gridLayout_7.setObjectName("gridLayout_7")
        self.lbErrorTolerance = QtWidgets.QLabel(self.layoutWidget_6)
        self.lbErrorTolerance.setObjectName("lbErrorTolerance")
        self.gridLayout_7.addWidget(self.lbErrorTolerance, 0, 0, 1, 1)
        self.leErrorTolerance = QtWidgets.QLineEdit(self.layoutWidget_6)
        self.leErrorTolerance.setObjectName("leErrorTolerance")
        self.gridLayout_7.addWidget(self.leErrorTolerance, 0, 1, 1, 1)
        self.lbUnit = QtWidgets.QLabel(self.layoutWidget_6)
        self.lbUnit.setObjectName("lbUnit")
        self.gridLayout_7.addWidget(self.lbUnit, 1, 0, 1, 1)
        self.lbAssignTo = QtWidgets.QLabel(self.layoutWidget_6)
        self.lbAssignTo.setObjectName("lbAssignTo")
        self.gridLayout_7.addWidget(self.lbAssignTo, 2, 0, 1, 1)
        self.cbAssignTo = QtWidgets.QComboBox(self.layoutWidget_6)
        self.cbAssignTo.setObjectName("cbAssignTo")
        self.cbAssignTo.addItems(['closest','highest intensity'])
        
        self.gridLayout_7.addWidget(self.cbAssignTo, 2, 1, 1, 1)
        self.cbMSMSUnit = QtWidgets.QComboBox(self.layoutWidget_6)
        self.cbMSMSUnit.setObjectName("cbMSMSUnit")
        self.cbMSMSUnit.addItems(['ppm','Da'])
        self.gridLayout_7.addWidget(self.cbMSMSUnit, 1, 1, 1, 1)
        self.tabWidget.addTab(self.MSMSTab, "")
        self.gridLayout.addWidget(self.tabWidget, 0, 0, 1, 1)

        self.retranslateUi()
        self.tabWidget.setCurrentIndex(0)


    def retranslateUi(self):
        _translate = QtCore.QCoreApplication.translate
        #self.setWindowTitle(_translate("Dialog", "Dialog"))
        self.gbNewProject.setTitle(_translate("Dialog", "Create a new project"))
        self.lbDLIMSProject.setText(_translate("Dialog", "DLIMS project"))
        self.lbProjectCode.setText(_translate("Dialog", "Project code"))
        self.lbProjectFolder.setText(_translate("Dialog", "Project folder*"))
        self.lbAnalyst.setText(_translate("Dialog", "Analyst*"))
        self.lbProjectName.setText(_translate("Dialog", "Project name*"))
        self.toolButton.setText(_translate("Dialog", "..."))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.ProjectTab), _translate("Dialog", "Project"))
        item = self.tableWidget.horizontalHeaderItem(0)
        item.setText(_translate("Dialog", "Label"))
        item = self.tableWidget.horizontalHeaderItem(1)
        item.setText(_translate("Dialog", "Sample Type"))
        item = self.tableWidget.horizontalHeaderItem(2)
        item.setText(_translate("Dialog", "Order"))
        item = self.tableWidget.horizontalHeaderItem(3)
        item.setText(_translate("Dialog", "Color"))
        item = self.tableWidget.horizontalHeaderItem(4)
        item.setText(_translate("Dialog", "Color code"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.TabFile), _translate("Dialog", "Files"))
        self.gbSmoothXIC.setTitle(_translate("Dialog", "Smooth"))
        self.cbApplySmoothing.setText(_translate("Dialog", "Apply smooting"))
        self.lbSelectSmoothType.setText(_translate("Dialog", "Select type"))
        self.lbCycles.setText(_translate("Dialog", "cycles "))
        self.lbWindow.setText(_translate("Dialog", "window"))
        self.gbGeneralXIC.setTitle(_translate("Dialog", "General"))
        self.lbCreateBasedOn.setText(_translate("Dialog", "Create based on "))
        self.lbXICTimeWindow.setText(_translate("Dialog", "Time window +/-"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.XICTab), _translate("Dialog", "XIC"))
        self.gbFullMSGeneral_3.setTitle(_translate("Dialog", "PYQMS"))
        self.lbErrorTolerance_2.setText(_translate("Dialog", "error tolerance in ppm"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.FullMSTab), _translate("Dialog", "Full MS"))
        self.gbGeneralMSMS.setTitle(_translate("Dialog", "General"))
        self.cbFragmentLosses.setText(_translate("Dialog", "Fragment losses"))
        self.cbIonX.setText(_translate("Dialog", "x"))
        self.cbIonY.setText(_translate("Dialog", "y"))
        self.cbIonZ.setText(_translate("Dialog", "z"))
        self.cbIonC.setText(_translate("Dialog", "c"))
        self.cbIonB.setText(_translate("Dialog", "b"))
        self.cbIonA.setText(_translate("Dialog", "a"))
        self.lbIonSeries.setText(_translate("Dialog", "Define ion series to consider "))

        self.lbErrorTolerance.setText(_translate("Dialog", "error tolerance"))
        self.lbUnit.setText(_translate("Dialog", "unit"))
        self.lbAssignTo.setText(_translate("Dialog", "assign to"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.MSMSTab), _translate("Dialog", "MSMS"))
        
        self.saveButton.setText(_translate("Dialog", "Save"))
        self.saveAsButton.setText(_translate("Dialog", "Save as"))
        self.loadButton.setText(_translate("Dialog", "Load"))
        self.okButton.setText(_translate("Dialog", "OK"))
        self.okButton.clicked.connect(self.close)
        
        
        
    def selectDir(self):
        self.openDialog=QFileDialog()
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        self.projectDirectory= QFileDialog.getExistingDirectory(self.openDialog,"select project directory" )
        self.leProjectFolder.setText(self.projectDirectory)

    def getColor(self):
        color = QColorDialog.getColor()
        focused=self.tableWidget.focusWidget()
        row=focused.property("row")
        colorButton=self.tableWidget.cellWidget(row,3)
        colorButton.setStyleSheet("background-color:%s"%(color.name()))
        self.tableWidget.setItem(row, 4, QtWidgets.QTableWidgetItem(color.name()))        


    def adjustColor(self,item):
        col=item.column()
        if col!=4:
            return
        row=item.row()
        colorButton=self.tableWidget.cellWidget(row,3)
        color=self.tableWidget.item(row,4).text()
        import re
        match = re.search(r'^#(?:[0-9a-fA-F]{3}){1,2}$', color)       
        if match: 
            colorButton.setStyleSheet("background-color:%s"%(color))
        else: 
            color=colorButton.palette().button().color()
            self.tableWidget.setItem(row, 4, QtWidgets.QTableWidgetItem(str(color.name())))
            

    def createEmptyXML(self):
        XML="""<?xml version='1.0' encoding='UTF-8'?>
        <configuration></configuration>        """
        return(XML)
        
    def createConfigFile(self):
        emptyXML=self.createEmptyXML()
        dom= xml.parseString(emptyXML)
        root=dom.childNodes[0]
        
        ProjectSection=dom.createElement("section")
        ProjectSection.setAttribute("name","project")
        
        projectName=dom.createElement("param")
        projectName.setAttribute("name","project name")
        projectName.setAttribute("value",self.leProjectName.text())
        ProjectSection.appendChild(projectName)
        
        analyst=dom.createElement("param")
        analyst.setAttribute("name","analyst")
        analyst.setAttribute("value",self.leAnalyst.text())
        ProjectSection.appendChild(analyst)

        projectDir=dom.createElement("param")
        projectDir.setAttribute("name","project directory")
        projectDir.setAttribute("value",self.leProjectFolder.text())
        ProjectSection.appendChild(projectDir)        

        projectCode=dom.createElement("param")
        projectCode.setAttribute("name","project code")
        projectCode.setAttribute("value",self.leProjectCode.text())
        ProjectSection.appendChild(projectCode)         

        DLIMS=dom.createElement("param")
        DLIMS.setAttribute("name","DLIMS project")
        DLIMS.setAttribute("value",self.leDLIMSProject.text())
        ProjectSection.appendChild(DLIMS)
        
        root.appendChild(ProjectSection)
        
        XICSection=dom.createElement("section")
        XICSection.setAttribute("name","XIC")
        
        basedOn=dom.createElement("param")
        basedOn.setAttribute("name","based on")
        basedOn.setAttribute("value",self.cbCreateBasedOn.currentText())
        XICSection.appendChild(basedOn)
        
        timeWindow=dom.createElement("param")
        timeWindow.setAttribute("name","time window")
        timeWindow.setAttribute("value",str(self.dsbTimeWindow.value()))
        XICSection.appendChild(timeWindow)
        
        applySmooth=dom.createElement("param")
        applySmooth.setAttribute("name","apply smooth")
        applySmooth.setAttribute("value",str(self.cbApplySmoothing.isChecked()))
        XICSection.appendChild(applySmooth)
        
        smoothType=dom.createElement("param")
        smoothType.setAttribute("name","apply type")
        smoothType.setAttribute("value",self.cbSmoothType.currentText())
        XICSection.appendChild(smoothType)        
        
        smoothCycles=dom.createElement("param")
        smoothCycles.setAttribute("name","smooth cycles")
        smoothCycles.setAttribute("value",str(self.sbCycles.value()))
        XICSection.appendChild(smoothCycles)         

        window=dom.createElement("param")
        window.setAttribute("name","smooth window")
        window.setAttribute("value",str(self.dsbSmoothWindow.value()))
        XICSection.appendChild(window) 
        
        root.appendChild(XICSection)
        
        FullMSSection=dom.createElement("section")
        FullMSSection.setAttribute("name","Full MS")
               
        
        errorTolerance=dom.createElement("param")
        errorTolerance.setAttribute("name","error tolerance")
        errorTolerance.setAttribute("value",self.leErrorTolerance_2.text())
        FullMSSection.appendChild(errorTolerance)        
 
        PYQMSppm=dom.createElement("param")
        PYQMSppm.setAttribute("name","PYQMSppm")
        PYQMSppm.setAttribute("value",str(self.leErrorTolerance_2.value()))
        FullMSSection.appendChild(PYQMSppm) 
       
        root.appendChild(FullMSSection)
        
        MSMSSection=dom.createElement("section")
        MSMSSection.setAttribute("name","MSMS")
                
        errorTolerance=dom.createElement("param")
        errorTolerance.setAttribute("name","error tolerance")
        errorTolerance.setAttribute("value",self.leErrorTolerance.text())
        MSMSSection.appendChild(errorTolerance) 

        unit=dom.createElement("param")
        unit.setAttribute("name","error unit")
        unit.setAttribute("value",self.cbMSMSUnit.currentText())
        MSMSSection.appendChild(unit) 
        
        assign=dom.createElement("param")
        assign.setAttribute("name","assign to")
        assign.setAttribute("value",self.cbAssignTo.currentText())
        MSMSSection.appendChild(assign)

        aIon=dom.createElement("param")
        aIon.setAttribute("name","a-ion")
        aIon.setAttribute("value",str(self.cbIonA.isChecked()))
        MSMSSection.appendChild(aIon) 

        bIon=dom.createElement("param")
        bIon.setAttribute("name","b-ion")
        bIon.setAttribute("value",str(self.cbIonB.isChecked()))
        MSMSSection.appendChild(bIon)         

        cIon=dom.createElement("param")
        cIon.setAttribute("name","c-ion")
        cIon.setAttribute("value",str(self.cbIonC.isChecked()))
        MSMSSection.appendChild(cIon) 

        xIon=dom.createElement("param")
        xIon.setAttribute("name","x-ion")
        xIon.setAttribute("value",str(self.cbIonX.isChecked()))
        MSMSSection.appendChild(xIon)

        yIon=dom.createElement("param")
        yIon.setAttribute("name","y-ion")
        yIon.setAttribute("value",str(self.cbIonY.isChecked()))
        MSMSSection.appendChild(yIon)

        zIon=dom.createElement("param")
        zIon.setAttribute("name","z-ion")
        zIon.setAttribute("value",str(self.cbIonZ.isChecked()))
        MSMSSection.appendChild(zIon)

        fragIon=dom.createElement("param")
        fragIon.setAttribute("name","fragment losses")
        fragIon.setAttribute("value",str(self.cbFragmentLosses.isChecked()))
        MSMSSection.appendChild(fragIon)

        root.appendChild(MSMSSection)
        
        return dom

    def saveConfigFile(self):
        valid=self.validateInput()
        if not valid:
            return
        dom=self.createConfigFile()
        
        if not self.configFile:
            filePath=path.join(self.projectDirectory,"config.configXML")
            filePath=path.normpath(filePath)
        else:
            filePath=self.configFile
        file_handle = open(filePath,"w")
        dom.writexml(file_handle,addindent='\t',newl='\n')
        file_handle.close()
        self.createMessageBox("Config file is saved in:\n %s"%(filePath))
        self.configFile=filePath
        
        

    def saveAsConfigFile(self):
        valid=self.validateInput()
        if not valid:
            return        
        filePath=QFileDialog.getSaveFileName(self, 'Save File',"","config file (*.configXML)")
        filePath=filePath[0]
        if not filePath:
            return
        file_handle = open(filePath,"w")
        dom=self.createConfigFile()
        dom.writexml(file_handle,addindent='\t',newl='\n')
        file_handle.close()
        self.configFile=filePath
        
    def validateInput(self):        
        if not self.leProjectFolder.text() or not self.leAnalyst.text() or not self.leProjectFolder.text():
            self.createMessageBox("Set a value for all mandatory input fields.")
            return False
        elif not path.isdir(self.leProjectFolder.text()):
            self.createMessageBox("Project folder is not valid!")
            return False         
        else:
            return True

    def closeWindow(self):
        return
        output=self.createWarnBox('Are you sure to quit?')
        if not output: 
            return
        valid=self.validateInput()
        if not valid:
            return
        self.close()


    def createMessageBox(self,output):
        msg = QMessageBox()
        msg.setIcon(QMessageBox.Information)
        msg.setText(output)
        msg.setStandardButtons(QMessageBox.Ok )
        msg.setWindowTitle("Message Box")
        msg.exec_()
        del msg

    def createWarnBox(self,output):
        msg = QtWidgets.QMessageBox()
        msg.setIcon(QtWidgets.QMessageBox.Information)
        msg.setText(output)
        msg.setStandardButtons(QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No)
        msg.setWindowTitle("Warning")
        reply=msg.exec_()
        if reply == 16384:
            output=True
        else:
            output=False      
        del msg
        return output

        
    def loadConfigFile(self):
        filePath=QFileDialog.getOpenFileName(self, 'Open File',"","config file (*.configXML)")
        filePath=filePath[0]
        print filePath
        if not filePath:
            return
        doc = xml.parse(filePath)
        sections=doc.getElementsByTagName("section")
        for section in sections:
            if section.getAttribute("name")=="XIC":
                elements=section.getElementsByTagName("param")
                for element in elements:
                    if element.getAttribute("name")=="based on":
                        self.cbCreateBasedOn.setCurrentText(element.getAttribute("value"))
                    
                    elif element.getAttribute("name")=="time window":
                        self.dsbTimeWindow.setValue(float(element.getAttribute("value")))
                        
                    elif element.getAttribute("name")=="apply smooth":
                        if element.getAttribute("value")=="True":
                            self.cbApplySmoothing.setChecked(True)
                        else:
                            self.cbApplySmoothing.setChecked(False)

                    elif element.getAttribute("name")=="apply type":
                        self.cbSmoothType.setCurrentText(element.getAttribute("value"))
                        
                    elif element.getAttribute("name")=="smooth cycles":
                        self.sbCycles.setValue(int(element.getAttribute("value")))                        
                    
                    elif element.getAttribute("name")=="smooth window":
                        self.dsbSmoothWindow.setValue(float(element.getAttribute("value")))                   
                
            if section.getAttribute("name")=="Full MS":  
                elements=section.getElementsByTagName("param")
                for element in elements:                        
                    if element.getAttribute("name")=="PYQMSppm":
                        self.leErrorTolerance_2.setValue(int(element.getAttribute("value")))

            if section.getAttribute("name")=="MSMS":  
                elements=section.getElementsByTagName("param")
                for element in elements:                        
                    if element.getAttribute("name")=="error tolerance":
                        self.leErrorTolerance.setText(element.getAttribute("value"))
                    
                    elif element.getAttribute("name")=="error unit":
                        self.cbMSMSUnit.setCurrentText(element.getAttribute("value"))

                    elif element.getAttribute("name")=="assign to":
                        self.cbAssignTo.setCurrentText(element.getAttribute("value"))

                    
                    elif element.getAttribute("name")=="a-ion":
                        if element.getAttribute("value")=="True":
                            self.cbIonA.setChecked(True)
                        else:
                            self.cbIonA.setChecked(False)

                    elif element.getAttribute("name")=="b-ion":
                        if element.getAttribute("value")=="True":
                            self.cbIonB.setChecked(True)
                        else:
                            self.cbIonB.setChecked(False)

                    elif element.getAttribute("name")=="c-ion":
                        if element.getAttribute("value")=="True":
                            self.cbIonC.setChecked(True)
                        else:
                            self.cbIonC.setChecked(False)        

                    elif element.getAttribute("name")=="x-ion":
                        if element.getAttribute("value")=="True":
                            self.cbIonX.setChecked(True)
                        else:
                            self.cbIonX.setChecked(False)

                    elif element.getAttribute("name")=="y-ion":
                        if element.getAttribute("value")=="True":
                            self.cbIonY.setChecked(True)
                        else:
                            self.cbIonY.setChecked(False)

                    elif element.getAttribute("name")=="z-ion":
                        if element.getAttribute("value")=="True":
                            self.cbIonZ.setChecked(True)
                        else:
                            self.cbIonZ.setChecked(False)

                    elif element.getAttribute("name")=="fragment losses":
                        if element.getAttribute("value")=="True":
                            self.cbFragmentLosses.setChecked(True)
                        else:
                            self.cbFragmentLosses.setChecked(False)
        
    def changeSmooth(self):
        if self.cbApplySmoothing.isChecked():
            self.cbSmoothType.setEnabled(True)
            self.sbCycles.setEnabled(True)
            self.dsbSmoothWindow.setEnabled(True)
        else:
            self.cbSmoothType.setEnabled(False)
            self.sbCycles.setEnabled(False)
            self.dsbSmoothWindow.setEnabled(False)        
        
    def updateDataFrame(self):
        import pandas as pd
        rows=list()

        self.DataFrame.drop(['label','sample', 'type','order','color'],errors='ignore',axis=1,inplace=True)
        for n,item in enumerate(self.DataFrame.index):
            label=self.tableWidget.item(n,0).text()
            sampleType=self.tableWidget.cellWidget(n,1).currentText()
            order=int(self.tableWidget.cellWidget(n,2).currentText())
            color=self.tableWidget.item(n,4).text()
            baseName=self.tableWidget.verticalHeaderItem(n).text()
            row=[label,baseName,sampleType,order,color]
            rows.append(row)        
        output_table=pd.DataFrame(data=rows,columns=['label','basename','sample type','order','color'])
        output_table=output_table.merge(self.DataFrame,on="basename")
        self.DataFrame=output_table
# =============================================================================
#         if not self.configFile:
#             filePath=path.join(self.projectDirectory,"config.configXML")
#             filePath=path.normpath(filePath)
#             self.configFile=filePath
#             dom=self.createConfigFile()
#             file_handle = open(filePath,"w")
#             dom.writexml(file_handle,addindent='\t',newl='\n')
#             file_handle.close()
# =============================================================================


        
            

        


 
    





