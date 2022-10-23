# -*- coding: utf-8 -*-
"""
Created on Sat May 25 08:19:13 2019

@author: KTMS458
"""
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from PyQt5 import QtCore,  QtWidgets, QtGui
from PyQt5.QtWidgets import QMainWindow

import pandas as pd
import numpy as np
import re
from MedImmune.util import modTable
from matplotlib.widgets import SpanSelector
from MedImmune.Graphs4 import ScalarFormatterForceFormat
from matplotlib.figure import Figure
import seaborn as sns


class Ui_MainWindow(QMainWindow):
    def __init__(self,input_table):
        super(QMainWindow,self).__init__()
        self.input_table=input_table
        self.dict_condition={}
        self.setupUi()
        self.savedItems=dict()
        self.tempItems=dict()
        self.output_barplot=pd.DataFrame(columns=['label','Modification','value'])
        
        import sys
        sys._excepthook = sys.excepthook 
        def exception_hook(exctype, value, traceback):
            print(exctype, value, traceback)
            sys._excepthook(exctype, value, traceback) 
            sys.exit(1) 
        sys.excepthook = exception_hook

    def closeEvent(self,event):
        output=self.createMessageBox('Are you sure to quit?')
        if not output: 
            event.ignore()
            return
        self.output_table=pd.DataFrame(self.savedItems)
    
    def setupUi(self):        
        self.component=None
        self.componentNumber=0
        self.quantitationType=None
        self.currentSettings={}
        self.grouped=False
        self.groupedDF=None
        self.savedItems=pd.DataFrame()
        self.output_table=pd.DataFrame()
        _translate = QtCore.QCoreApplication.translate
        self.setObjectName("MainWindow")
        self.resize(1317, 819)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Maximum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.sizePolicy().hasHeightForWidth())
        self.setSizePolicy(sizePolicy)
        self.setAutoFillBackground(True)
        self.centralwidget = QtWidgets.QWidget(self)
        self.centralwidget.setObjectName("centralwidget")
        self.gridLayout_4 = QtWidgets.QGridLayout(self.centralwidget)
        self.gridLayout_4.setObjectName("gridLayout_4")


        self.splitter = QtWidgets.QSplitter(self.centralwidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.splitter.sizePolicy().hasHeightForWidth())
        self.splitter.setSizePolicy(sizePolicy)
        self.splitter.setOrientation(QtCore.Qt.Horizontal)
        self.splitter.setObjectName("splitter")
        self.layoutWidget = QtWidgets.QWidget(self.splitter)
        self.layoutWidget.setObjectName("layoutWidget")
        self.gridLayoutLeft = QtWidgets.QGridLayout(self.layoutWidget)
        self.gridLayoutLeft.setContentsMargins(0, 0, 0, 0)
        self.gridLayoutLeft.setObjectName("gridLayoutLeft")
        self.verticalLayoutTableView = QtWidgets.QVBoxLayout()
        self.verticalLayoutTableView.setObjectName("verticalLayoutTableView")
        self.ModTable = QtWidgets.QTableWidget(self.layoutWidget)
        self.ModTable.setMinimumSize(QtCore.QSize(300, 300))
        self.ModTable.setMaximumSize(QtCore.QSize(1280, 1281))
        self.ModTable.setObjectName("ModTable")
        
        self.ModTable.setColumnCount(3)
        item = QtWidgets.QTableWidgetItem()
        self.ModTable.setHorizontalHeaderItem(0, item)
        item.setText(_translate("Ui_MainWindow", "Protein")) 
        item = QtWidgets.QTableWidgetItem()
        self.ModTable.setHorizontalHeaderItem(1, item)
        item.setText(_translate("Ui_MainWindow", "Residue"))        
        item = QtWidgets.QTableWidgetItem()
        self.ModTable.setHorizontalHeaderItem(2, item)
        item.setText(_translate("Ui_MainWindow", "Modification"))

        
        header = self.ModTable.horizontalHeader()       
        header.setSectionResizeMode(0, QtWidgets.QHeaderView.ResizeToContents)
        header.setSectionResizeMode(1, QtWidgets.QHeaderView.ResizeToContents)
        header.setSectionResizeMode(2, QtWidgets.QHeaderView.Stretch)
        
        df=modTable(self.input_table)
        n =int(len(df))
        m =int(len(df.columns))
        self.ModTable.setRowCount(n)
        
        for a in range(n):
            for b in range(m):
                item=df.iloc[a,b]
                self.ModTable.setItem(a, b, QtWidgets.QTableWidgetItem(str(item)))
        
        self.ModTable.setSelectionBehavior(QtWidgets.QTableWidget.SelectRows)
        self.ModTable.setSelectionMode(QtWidgets.QTableWidget.SingleSelection)
        self.ModTable.clicked.connect(self.updatePeptideTable)

        

        
        self.verticalLayoutTableView.addWidget(self.ModTable)
        self.PeptideTable = QtWidgets.QTableWidget(self.layoutWidget)
        self.PeptideTable.setMinimumSize(QtCore.QSize(300, 300))
        self.PeptideTable.setMaximumSize(QtCore.QSize(1280, 1280))
        self.PeptideTable.setSelectionBehavior(QtWidgets.QTableWidget.SelectRows)
        self.PeptideTable.setSelectionMode(QtWidgets.QTableWidget.SingleSelection)
        
        self.PeptideTable.setObjectName("PeptideTable")
        
        
        self.PeptideTable.setColumnCount(14)
        header = self.PeptideTable.horizontalHeader()
        item = QtWidgets.QTableWidgetItem()
        self.PeptideTable.setHorizontalHeaderItem(0, item)
        item.setText(_translate("Ui_MainWindow", "Checked"))
        header.setSectionResizeMode(0, QtWidgets.QHeaderView.ResizeToContents)
        
        self.df_full=self.input_table.filter(['No.','Protein','Peptide Sequence','Modification','Site','Delta (ppm)',
                'Confidence Score','ID Type','RT','M/Z', 'Charge State', 'Mono Mass Exp.']).copy()
        
        for i,col in enumerate(self.df_full.columns):        
                item = QtWidgets.QTableWidgetItem()
                self.PeptideTable.setHorizontalHeaderItem(i+1, item)
                item.setText(_translate("Ui_MainWindow", col)) 
                header.setSectionResizeMode(i+1, QtWidgets.QHeaderView.ResizeToContents)
                header.setSectionResizeMode(i+1, QtWidgets.QHeaderView.Interactive)
        
        
        item = QtWidgets.QTableWidgetItem()
        self.PeptideTable.setHorizontalHeaderItem(1, item)
        item.setText(_translate("Ui_MainWindow", "#Component"))         
        
        item = QtWidgets.QTableWidgetItem()
        self.PeptideTable.setHorizontalHeaderItem(13, item)
        item.setText(_translate("Ui_MainWindow", "Group")) 
        
        
        self.verticalLayoutTableView.addWidget(self.PeptideTable)
        self.gridLayoutLeft.addLayout(self.verticalLayoutTableView, 0, 0, 1, 4)
        self.frame_2 = QtWidgets.QFrame(self.layoutWidget)
        self.frame_2.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame_2.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame_2.setObjectName("frame_2")
        self.gridLayoutLeft.addWidget(self.frame_2, 1, 0, 1, 1)
        self.formLayoutChargeState = QtWidgets.QFormLayout()
        self.formLayoutChargeState.setObjectName("formLayoutChargeState")
        self.groupBoxFilter = QtWidgets.QGroupBox(self.layoutWidget)
        self.groupBoxFilter.setObjectName("groupBoxFilter")
        self.gridLayout = QtWidgets.QGridLayout(self.groupBoxFilter)
        self.gridLayout.setObjectName("gridLayout")
        self.lbChargeState = QtWidgets.QLabel(self.groupBoxFilter)
        self.lbChargeState.setObjectName("lbChargeState")
        self.gridLayout.addWidget(self.lbChargeState, 1, 3, 1, 2)
        self.sbTop = QtWidgets.QSpinBox(self.groupBoxFilter)
        self.sbTop.setObjectName("sbTop")
        self.gridLayout.addWidget(self.sbTop, 1, 1, 1, 2)
        self.cbDefineChargeState = QtWidgets.QCheckBox(self.groupBoxFilter)
        self.cbDefineChargeState.setObjectName("cbDefineChargeState")
        self.gridLayout.addWidget(self.cbDefineChargeState, 3, 0, 1, 4)
        self.leChargeState = QtWidgets.QLineEdit(self.groupBoxFilter)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.leChargeState.sizePolicy().hasHeightForWidth())
        self.leChargeState.setSizePolicy(sizePolicy)
        self.leChargeState.setObjectName("leChargeState")
        self.gridLayout.addWidget(self.leChargeState, 3, 4, 1, 1)
        self.cbTop = QtWidgets.QCheckBox(self.groupBoxFilter)
        self.cbTop.setObjectName("cbTop")
        self.gridLayout.addWidget(self.cbTop, 1, 0, 1, 1)
        self.cbAllChargeStates = QtWidgets.QCheckBox(self.groupBoxFilter)
        self.cbAllChargeStates.setObjectName("cbAllChargeStates")
        self.gridLayout.addWidget(self.cbAllChargeStates, 2, 0, 1, 5)
        self.cbHighestIntenseChargeState = QtWidgets.QCheckBox(self.groupBoxFilter)
        self.cbHighestIntenseChargeState.setObjectName("cbHighestIntenseChargeState")
        self.cbHighestIntenseChargeState.stateChanged.connect(self.filterStateChanged)
        
        
        self.gridLayout.addWidget(self.cbHighestIntenseChargeState, 0, 0, 1, 5)
        self.formLayoutChargeState.setWidget(0, QtWidgets.QFormLayout.SpanningRole, self.groupBoxFilter)
        self.lbFilterBasedOnID = QtWidgets.QLabel(self.layoutWidget)
        self.lbFilterBasedOnID.setObjectName("lbFilterBasedOnID")
        self.formLayoutChargeState.setWidget(1, QtWidgets.QFormLayout.LabelRole, self.lbFilterBasedOnID)
        self.cbID = QtWidgets.QComboBox(self.layoutWidget)
        self.cbID.setMinimumSize(QtCore.QSize(201, 20))
        self.cbID.setMaximumSize(QtCore.QSize(201, 20))
        self.cbID.setObjectName("cbID")
        self.cbID.addItems(['any','both','MS2'])
        self.formLayoutChargeState.setWidget(1, QtWidgets.QFormLayout.FieldRole, self.cbID)
        self.gridLayoutLeft.addLayout(self.formLayoutChargeState, 1, 1, 1, 1)
        self.verticalLayoutLeftBottom = QtWidgets.QVBoxLayout()
        self.verticalLayoutLeftBottom.setObjectName("verticalLayoutLeftBottom")
        self.btGroup = QtWidgets.QPushButton(self.layoutWidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.btGroup.sizePolicy().hasHeightForWidth())
        self.btGroup.setSizePolicy(sizePolicy)
        self.btGroup.setMinimumSize(QtCore.QSize(130, 23))
        self.btGroup.setMaximumSize(QtCore.QSize(130, 23))
        self.btGroup.setObjectName("btGroup")
        self.btGroup.clicked.connect(self.groupComponent)
        
        self.verticalLayoutLeftBottom.addWidget(self.btGroup)
        self.leGroup=QtWidgets.QLineEdit(self.layoutWidget)
        self.leGroup.setSizePolicy(sizePolicy)
        self.leGroup.setMinimumSize(QtCore.QSize(130, 23))
        self.leGroup.setMaximumSize(QtCore.QSize(130, 23))
        self.leGroup.setObjectName("leGroup")
        self.verticalLayoutLeftBottom.addWidget(self.leGroup)
        self.btReset = QtWidgets.QPushButton(self.layoutWidget)
        self.btReset.setMinimumSize(QtCore.QSize(130, 23))
        self.btReset.setMaximumSize(QtCore.QSize(130, 23))
        self.btReset.setObjectName("btReset")
        self.btReset.clicked.connect(self.reset)
        self.verticalLayoutLeftBottom.addWidget(self.btReset)
        self.frame = QtWidgets.QFrame(self.layoutWidget)
        self.frame.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame.setObjectName("frame")
        self.verticalLayoutLeftBottom.addWidget(self.frame)
        self.gridLayoutLeft.addLayout(self.verticalLayoutLeftBottom, 1, 2, 1, 1)
        self.frame_3 = QtWidgets.QFrame(self.layoutWidget)
        self.frame_3.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame_3.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame_3.setObjectName("frame_3")
        self.gridLayoutLeft.addWidget(self.frame_3, 1, 3, 1, 1)
        self.layoutWidget1 = QtWidgets.QWidget(self.splitter)
        self.layoutWidget1.setObjectName("layoutWidget1")
        self.gridLayoutRight = QtWidgets.QGridLayout(self.layoutWidget1)
        self.gridLayoutRight.setContentsMargins(0, 0, 0, 0)
        self.gridLayoutRight.setObjectName("gridLayoutRight")
        self.horizontalLayoutRightMiddle = QtWidgets.QHBoxLayout()
        self.horizontalLayoutRightMiddle.setObjectName("horizontalLayoutRightMiddle")
                
        self.btProcess = QtWidgets.QPushButton(self.layoutWidget1)
        self.btProcess.setObjectName("btProcess")
        self.btProcess.clicked.connect(self.updateBarPlot)
        self.horizontalLayoutRightMiddle.addWidget(self.btProcess)
#        self.progressBar = QtWidgets.QProgressBar(self.layoutWidget1)
#        self.progressBar.setProperty("value", 0)
#        self.progressBar.setObjectName("progressBar")
#        self.horizontalLayoutRightMiddle.addWidget(self.progressBar)
        spacerItem = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayoutRightMiddle.addItem(spacerItem)
        
        self.FileSelection = CheckableComboBox(self.layoutWidget1)
        self.FileSelection.setObjectName("FileSelection")
        component=self.input_table['component'].iloc[0]
        
        fileNames=list()
        for fileName in component.rawfile_dict.keys():
            if component.rawfile_dict[fileName]['condition']:
                fileNames.append(component.rawfile_dict[fileName]['condition'])
                self.dict_condition[component.rawfile_dict[fileName]['condition']]=fileName
            else:
                fileNames.append(fileName)   
                self.dict_condition[fileName]=fileName
        self.FileSelection.addItems(fileNames)
        self.FileSelection.setMinimumSize(QtCore.QSize(400, 20))
        
        self.horizontalLayoutRightMiddle.addWidget(self.FileSelection)
        
        self.btSelectFiles = QtWidgets.QPushButton(self.layoutWidget1)
        self.btSelectFiles.setObjectName("btSelectFiles")
        self.btSelectFiles.clicked.connect(self.filterLines)
        self.horizontalLayoutRightMiddle.addWidget(self.btSelectFiles)
        
        

        self.gridLayoutRight.addLayout(self.horizontalLayoutRightMiddle, 2, 0, 1, 1)
        self.horizontalLayoutRightBottom = QtWidgets.QHBoxLayout()
        self.horizontalLayoutRightBottom.setObjectName("horizontalLayoutRightBottom")
        self.lbQuantType = QtWidgets.QLabel(self.layoutWidget1)
        self.lbQuantType.setMinimumSize(QtCore.QSize(115, 20))
        self.lbQuantType.setMaximumSize(QtCore.QSize(115, 20))
        self.lbQuantType.setObjectName("lbQuantType")
        self.horizontalLayoutRightBottom.addWidget(self.lbQuantType)
        self.cbQuantType = QtWidgets.QComboBox(self.layoutWidget1)
        self.cbQuantType.setMinimumSize(QtCore.QSize(151, 21))
        self.cbQuantType.setMaximumSize(QtCore.QSize(151, 21))
        self.cbQuantType.setObjectName("cbQuantType")
        self.cbQuantType.addItems(['MassArea (BiopharmaFinder)','mXIC','hXIC','PYQMS'])
        self.cbQuantType.currentTextChanged.connect(self.changeQuantType)
        
        
        self.horizontalLayoutRightBottom.addWidget(self.cbQuantType)
        spacerItem = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayoutRightBottom.addItem(spacerItem)        
        self.cbModName=QtWidgets.QComboBox(self.layoutWidget1)
        self.cbModName.setObjectName("cbModName")
        self.cbModName.setSizeAdjustPolicy(QtWidgets.QComboBox.AdjustToContents)
        self.horizontalLayoutRightBottom.addWidget(self.cbModName) 
        self.btDelete = QtWidgets.QPushButton(self.layoutWidget1)
        self.btDelete.setObjectName("Delete")
        self.btDelete.clicked.connect(self.deleteItem)
        self.horizontalLayoutRightBottom.addWidget(self.btDelete) 
        self.btSave = QtWidgets.QPushButton(self.layoutWidget1)
        self.btSave.setObjectName("btSave")
        self.btSave.clicked.connect(self.saveItem)
        
        self.horizontalLayoutRightBottom.addWidget(self.btSave)
        self.btLoad = QtWidgets.QPushButton(self.layoutWidget1)
        self.btLoad.setObjectName("btLoad")
        self.btLoad.clicked.connect(self.loadItem)
        self.horizontalLayoutRightBottom.addWidget(self.btLoad)                 
        self.gridLayoutRight.addLayout(self.horizontalLayoutRightBottom, 4, 0, 1, 1)


        self.gvBarplot=BarPlotCanvas()
        self.gvBarplot.setObjectName("gvBarplot")
        self.gridLayoutRight.addWidget(self.gvBarplot, 0, 0, 1, 1)
 
        self.gvQuantTrace=LinePlotCanvas()
        self.gvQuantTrace.setObjectName("gvQuantTrace")
        self.gridLayoutRight.addWidget(self.gvQuantTrace, 3, 0, 1, 1)       

        
        self.gridLayout_4.addWidget(self.splitter, 1, 0, 1, 1)
        self.setCentralWidget(self.centralwidget)


        self.retranslateUi(self)
        QtCore.QMetaObject.connectSlotsByName(self)

    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "MainWindow"))
        self.groupBoxFilter.setTitle(_translate("MainWindow", " Fiilter based on charge state"))
        self.lbChargeState.setText(_translate("MainWindow", "charge state based on intensity"))
        self.cbDefineChargeState.setText(_translate("MainWindow", "Define charge states"))
        self.cbTop.setText(_translate("MainWindow", "Top"))
        self.cbAllChargeStates.setText(_translate("MainWindow", "All charge states"))
        self.cbHighestIntenseChargeState.setText(_translate("MainWindow", "Highest intense charge state (unmodified)"))
        self.lbFilterBasedOnID.setText(_translate("MainWindow", "Filter based on id type"))
        self.btGroup.setText(_translate("MainWindow", "Group modified species"))
        self.btReset.setText(_translate("MainWindow", "Reset"))
        self.btProcess.setText(_translate("MainWindow", "Process"))
        self.btSave.setText(_translate("MainWindow", "Process and Save Item"))
        self.btDelete.setText(_translate("MainWindow", "Delete Item"))
        self.btLoad.setText(_translate("MainWindow", "Load"))
        self.btSelectFiles.setText(_translate("MainWindow", "Select Files"))
        self.lbQuantType.setText(_translate("MainWindow", "Select quantitation type"))

        #self.actionLoad_config_file.setText(_translate("MainWindow", "Load config file"))
    
    def updatePeptideTable(self,index,buff=None):
        self.grouped=False
        self.PeptideTable.setRowCount(0)
        self.gvBarplot.ax.clear()
        self.gvBarplot.draw()
        selectedRow=self.ModTable.currentRow()
        protein=self.ModTable.item(selectedRow,0).data(0)
        residue=self.ModTable.item(selectedRow,1).data(0)
        residue=int(re.findall(r'\d+',residue)[0])
        
        if buff:
            self.buff=buff            
        else:
            self.buff=self.df_full[self.input_table['Protein'].str.contains(protein) & (self.input_table['begin']<=residue)&(self.input_table['end']>=residue)].copy()
        self.buff.sort_values("Peptide Sequence", inplace=True)
        rows =int(len(self.buff))
        columns =int(len(self.buff.columns))
        self.PeptideTable.setRowCount(rows)
        for a in range(rows):
            checkBoxWidget=QtWidgets.QWidget()            
            layout=QtWidgets.QVBoxLayout(checkBoxWidget)            
            checkBox = QtWidgets.QCheckBox()
            checkBox.setProperty("row",a)
            layout.addWidget(checkBox)
            layout.setAlignment(QtCore.Qt.AlignCenter)
            layout.setContentsMargins(0,0,0,0)
            self.PeptideTable.setCellWidget(a,0,checkBoxWidget)
            checkBox.stateChanged.connect(self.state_changed)
            for b in range(columns):
                item=self.buff.iloc[a,b]
                self.PeptideTable.setItem(a, b+1, QtWidgets.QTableWidgetItem(str(item)))
        
        self.PeptideTable.clicked.connect(self.updateLinePlot)
        self.tempQuant={}        
        self.PeptideTable.setCurrentCell(0,0)
        index=self.PeptideTable.currentIndex()
        self.updateLinePlot(index)


    def state_changed(self,event):
        focused=self.PeptideTable.focusWidget()
        row=focused.property("row")
        checkBoxWidget=self.PeptideTable.cellWidget(row,0)
        checkBox=checkBoxWidget.children()[1]

        if checkBox.isChecked():
            font = QtGui.QFont()
            font.setBold(True)
            columns=self.PeptideTable.columnCount()
            for col in range(1,columns):
                try: self.PeptideTable.item(row, col).setFont(font)
                except: pass

        else:
            font = QtGui.QFont()
            font.setBold(False)
            columns=self.PeptideTable.columnCount()
            for col in range(1,columns):
                try: self.PeptideTable.item(row, col).setFont(font)
                except: pass
        if not self.quantitationType=='MassArea (BiopharmaFinder)':
            self.gvQuantTrace.span.set_active(True)
            
        else:
            self.gvQuantTrace.span.set_active(False)
            
    
    def updateLinePlot(self,clickedIndex):        
        row=clickedIndex.row()     
        componentNumber=self.buff['No.'].iloc[row]
        component=self.input_table['component'][self.input_table['No.']==componentNumber][0]
        kind=self.cbQuantType.currentText()
        if self.componentNumber and self.componentNumber==componentNumber and kind==self.quantitationType: return
        self.updateTempQuant()
        self.component=component
        self.componentNumber=componentNumber
        self.quantitationType=kind
        
        if kind=='MassArea (BiopharmaFinder)':
            kind="hXIC"
        
        self.gvQuantTrace.text=None
        
        if self.componentNumber in self.tempQuant.keys():
            tempQuant=self.tempQuant[self.componentNumber]
            self.gvQuantTrace.changeComponent(self.component,kind,tempQuant=tempQuant)
            self.filterLines()
        else:
            self.gvQuantTrace.changeComponent(self.component,kind)
            self.filterLines()
           
        
        buff=self.PeptideTable.cellWidget(row,0)
        if buff:
            buff=buff.children()[1]
        else:
            return
        if not self.quantitationType=='MassArea (BiopharmaFinder)':
            self.gvQuantTrace.span.set_active(True)
        else:
            self.gvQuantTrace.span.set_active(False)                        


    def updateTempQuant(self):
        kind=self.cbQuantType.currentText()
        for key in self.gvQuantTrace.artist_dict.keys():
            if  self.gvQuantTrace.artist_dict[key]['limit'] and self.componentNumber:
                if self.componentNumber in self.tempQuant.keys():
                    if key in self.tempQuant[self.componentNumber].keys():
                        self.tempQuant[self.componentNumber][key][kind]=self.gvQuantTrace.artist_dict[key]['limit']
                    else:    
                        self.tempQuant[self.componentNumber][key]={kind:self.gvQuantTrace.artist_dict[key]['limit']}
                else:
                    self.tempQuant[self.componentNumber]={key:{kind:self.gvQuantTrace.artist_dict[key]['limit']}}
    
    
    
    def getCheckedRows(self):
        rowCount=self.PeptideTable.rowCount()
        ls_checked=list()
        for row in range(rowCount):
            checkBoxWidget=self.PeptideTable.cellWidget(row,0)
            checkBox=checkBoxWidget.children()[1]
            if checkBox.isChecked():
                ls_checked.append(row)
        return(ls_checked)

  
    def updateBarPlot(self):
        self.tempItems=dict()
        ls_componentNumber=list()
        currentIndex=self.PeptideTable.currentIndex()
        self.updateLinePlot(currentIndex)
        self.updateTempQuant()
        ls_checked=self.getCheckedRows()
        if self.grouped:
            ls_group=list()   
            for row in ls_checked:
                componentNumber=int(self.PeptideTable.item(row,1).text())
                ls_componentNumber.append(componentNumber)
                group=self.PeptideTable.item(row,13).text()
                ls_group.append(group)
            
            self.df_barplot=self.input_table[self.input_table['No.'].isin(ls_componentNumber)].copy()
            buff=pd.DataFrame()
            buff['No.']=ls_componentNumber
            buff['Group']=ls_group
            self.df_barplot=pd.merge(self.df_barplot,buff,how='left',on='No.')
            self.groupedDF=self.df_barplot.groupby('Group').apply(self.getModQuant)
            self.groupedDF.reset_index(inplace=True)
            self.groupedDF.drop(['level_1'],inplace=True,axis=1)
            self.groupedDF.set_index('Group',inplace=True)
            maxMod=0
            maxCol=""
            for n,col in enumerate(self.groupedDF.columns):
                buffMax=self.groupedDF[col].max()
                if buffMax>maxMod:
                    maxMod=buffMax
                    maxCol=col
            
            self.groupedDF.sort_values(maxCol,ascending=False,inplace=True)
          
            self.groupedDF=self.groupedDF.transpose()
            for col in self.groupedDF.columns:
                if not col in ["None","unmodified"]:
                    colName="%"+col
                    self.groupedDF[colName]=self.groupedDF[col]/self.groupedDF[self.groupedDF.columns].sum(axis=1)*100
            
            ls_label=list()
            ls_color=list()
            for item in self.groupedDF.index:
                ls_label.append(self.gvQuantTrace.artist_dict[item]['label'])
                ls_color.append(self.gvQuantTrace.artist_dict[item]['color'])
            self.groupedDF['label']=ls_label
            self.groupedDF['color']=ls_color
            self.gvBarplot.updateBarplot(self.groupedDF)
            
         
        else: 
            for row in ls_checked:
                componentNumber=int(self.PeptideTable.item(row,1).text())
                ls_componentNumber.append(componentNumber)            
            self.df_barplot=self.input_table[self.input_table['No.'].isin(ls_componentNumber)].copy()                  
            self.groupedDF=self.df_barplot.groupby('Modification').apply(self.getModQuant)
            self.groupedDF.reset_index(inplace=True)
            self.groupedDF.drop(['level_1'],inplace=True,axis=1)
            self.groupedDF.set_index('Modification',inplace=True)
            
            maxMod=0
            maxCol=""
            for n,col in enumerate(self.groupedDF.columns):
                buffMax=self.groupedDF[col].max()
                if buffMax>maxMod:
                    maxMod=buffMax
                    maxCol=col
            
            self.groupedDF.sort_values(maxCol,ascending=False,inplace=True)
                        
            self.groupedDF=self.groupedDF.transpose()
            for col in self.groupedDF.columns:
                if not col in ["None","unmodified"]:
                    colName="%"+col
                    self.groupedDF[colName]=self.groupedDF[col]/self.groupedDF[self.groupedDF.columns].sum(axis=1)*100
            ls_label=list()
            ls_color=list()
            for item in self.groupedDF.index:
                ls_label.append(self.gvQuantTrace.artist_dict[item]['label'])
                ls_color.append(self.gvQuantTrace.artist_dict[item]['color'])
            self.groupedDF['label']=ls_label
            self.groupedDF['color']=ls_color
            self.gvBarplot.updateBarplot(self.groupedDF)
                
    
    def validatePeptideTable(self):
        pass

        
    def changeQuantType(self):
        index=self.PeptideTable.currentIndex()
        self.updateLinePlot(index)
        
        
    def groupComponent(self):
        group=self.leGroup.text()
        if not group:
            return

        ls_group=list()
        rowCount=self.PeptideTable.rowCount()
        for row in range(rowCount):
            checkBoxWidget=self.PeptideTable.cellWidget(row,0)
            checkBox=checkBoxWidget.children()[1]
            if checkBox.isChecked():
                self.PeptideTable.setItem(row, 13, QtWidgets.QTableWidgetItem(group))
                font = QtGui.QFont()
                font.setBold(True)        
                self.PeptideTable.item(row, 13).setFont(font)
                ls_group.append(group)
                
            elif not self.PeptideTable.item(row,13): 
                self.PeptideTable.setItem(row, 13, QtWidgets.QTableWidgetItem('None'))
                ls_group.append('None')
            elif self.PeptideTable.item(row,13).text()=='None':
                ls_group.append('None')
        
        self.buff['group']=ls_group
        self.grouped=True
            
    def reset(self):
        self.gvQuantTrace.artist_dict={}
        self.updatePeptideTable(0)
        self.gvQuantTrace.draw()
           
    
    def filterStateChanged(self):
        if self.cbHighestIntenseChargeState.isChecked():
            self.cbTop.setChecked(False)
            self.cbAllChargeStates.setChecked(False)
            self.cbDefineChargeState.setChecked(False)
        if self.cbTop.isChecked():
            self.cbHighestIntenseChargeState.setChecked(False)
            self.cbAllChargeStates.setChecked(False)
            self.cbDefineChargeState.setChecked(False)
        if self.cbAllChargeStates.isChecked():
            self.cbHighestIntenseChargeState.setChecked(False)
            self.cbTop.setChecked(False)
            self.cbDefineChargeState.setChecked(False)        
        if self.cbDefineChargeState.isChecked():
            self.cbHighestIntenseChargeState.setChecked(False)
            self.cbTop.setChecked(False)
            self.cbAllChargeStates.setChecked(False)

            
    def getModQuant(self,input_table):
        kind=self.quantitationType
        output_dict={}
        component=input_table['component'].iloc[0]
        for key in component.rawfile_dict.keys():
            output_dict[key]=[0]

        if kind=='MassArea (BiopharmaFinder)':
            for n,row in input_table.iterrows():
                componentNumber=row['No.']
                self.tempItems[componentNumber]={}               
                for key in output_dict.keys():
                    quantValue=row[key]                
                    output_dict[key]=[output_dict[key][0]+quantValue]
                    self.tempItems[componentNumber][key]=[round(quantValue,1),0.0,0.0]
            
            output_table=pd.DataFrame.from_dict(output_dict)
            return(output_table)            

       
        for n,row in input_table.iterrows():
            component=row['component']
            componentNumber=row['No.']
            self.tempItems[componentNumber]={}
            for key in component.rawfile_dict.keys():
                quantValue=0
                if componentNumber in self.tempQuant.keys():                    
                    if key in self.tempQuant[componentNumber].keys(): 
                        if kind in self.tempQuant[componentNumber][key].keys():
                            RT=component.rawfile_dict[key]['quant'][kind]['RT']  
                            intensity=component.rawfile_dict[key]['quant'][kind]['intensity']
                            buff=pd.DataFrame()
                            buff['RT']=RT
                            buff['intensity']=intensity
                            low, high=self.tempQuant[componentNumber][key][kind]
                            buff=buff[(buff['RT']>=low)& (buff['RT']<=high)]
                            quantValue=np.trapz(buff['intensity'],buff['RT'])
                            self.tempItems[componentNumber][key]=[round(quantValue,1),round(low,3),round(high,3)]
                
                if not quantValue:
                    quantValue=component.rawfile_dict[key]['quant'][kind]['value']
                    #low,high=component.rawfile_dict[key]['quant'][kind]['limit']
                    self.tempItems[componentNumber][key]=[round(quantValue,1),0.0,0.0]
                output_dict[key]=[output_dict[key][0]+quantValue]
        output_table=pd.DataFrame.from_dict(output_dict)
        return(output_table)
        
    def saveItem(self):
        self.updateBarPlot()
        selectedRow=self.ModTable.currentRow()
        if self.groupedDF is None:
            return
        protein=self.ModTable.item(selectedRow,0).data(0)
        residue=self.ModTable.item(selectedRow,1).data(0)
        itemName="%s_%s"%(protein,residue)
        for col in self.groupedDF.columns:
            if col[0]=="%":
                itemName+="_%s"%(col[1:])
        
        if self.quantitationType=='MassArea (BiopharmaFinder)':
            quantitationType="BiopharmaFinder"
        else:
            quantitationType=self.quantitationType
          
        
        itemName+="_%s"%(quantitationType)
        
        if itemName in self.savedItems.keys():
            overwrite=self.createMessageBox("Do you want to overwrite existing item?")
            if not overwrite:
                return
        
        Barplot_table=self.groupedDF.drop(labels=['color','None'],axis=1).melt(id_vars=['label']).copy()
        Barplot_table['quantID']=[itemName for i in Barplot_table.index]
        

                
        from io import BytesIO
        from PIL import Image
        ls_quantTrace=list()
        ls_componentNumber=list()
        for key in self.tempItems.keys():
            componentNumber=key
            ls_componentNumber.append(componentNumber)
            component=self.df_barplot['component'][self.df_barplot['No.']==componentNumber].iloc[0]
            output_fig=LinePlotCanvas(width=3,height=1.5,alpha=0.1)
            kind=self.quantitationType
            if kind=='MassArea (BiopharmaFinder)':
                kind="hXIC"
            output_fig.changeComponent(component,kind,self.tempQuant)
            output_fig.ax.axis("off")
            output_fig.ax.set_xticklabels=([])
            output_fig.ax.set_yticklabels=([])
            for line in output_fig.ax.get_lines():
                line.set(linewidth=0.5)
            output_fig.ax.set_xlabel("")
            for spine in output_fig.ax.spines.values():
                spine.set_visible(False)
            extent = output_fig.ax.get_window_extent().transformed(output_fig.figure.dpi_scale_trans.inverted())            
            fig_buffer=BytesIO()
            output_fig.figure.savefig(fig_buffer,bbox_inches = extent,pad_inches = 0)
            img=Image.open(fig_buffer)
            ls_quantTrace.append(img)
            quantTraceDict=dict(componentNumber=ls_componentNumber, images=ls_quantTrace)
        
        self.savedItems[itemName]=dict(out_barplot=Barplot_table, components=self.tempItems, 
                                       quantTrace=quantTraceDict, groupedDF=self.groupedDF.copy(),
                                       selectedRow=selectedRow, quantType=self.quantitationType)

        self.cbModName.clear()
        self.cbModName.addItems(self.savedItems.keys())    
        
        if self.grouped:
            ls_group=list(self.buff['group'])
            self.savedItems[itemName]['group']=ls_group


    
    def deleteItem(self):
        item=self.cbModName.currentText()
        if item:            
            delete=self.createMessageBox("Delet item: %s"%(item))
            if delete:
                del self.savedItems[item]
            else:
                return
        self.cbModName.clear()
        self.cbModName.addItems(self.savedItems.keys())
        return
        
        
    def loadItem(self):
        item=self.cbModName.currentText()
        print 1
        if not item:
            return
        item=self.savedItems[item]
        selectedRow=item['selectedRow']
        self.ModTable.setCurrentCell(selectedRow,0)
        self.updatePeptideTable(0)
        quantType=item['quantType']
        self.groupedDF=item['groupedDF']
        index=self.cbQuantType.findText(quantType)
        self.cbQuantType.setCurrentIndex(index)


        rowCount=self.PeptideTable.rowCount()
        #add group
        if "group" in item.keys():
            ls_group=item['group']            
            for row in range(rowCount):                  
                self.PeptideTable.setItem(row, 13, QtWidgets.QTableWidgetItem(ls_group[row]))

        for row in range(rowCount):
            componentNumber=int(self.PeptideTable.item(row,1).text())
            if componentNumber in item['components'].keys():         
                checkBoxWidget=self.PeptideTable.cellWidget(row,0)
                checkBox=checkBoxWidget.children()[1]
                checkBox.setFocus()
                checkBox.setChecked(True)
                 
        temp_dict=dict()
        for key in item['components'].keys():
            componentNumber=key
            buff=dict()
            for rawfile in item['components'][componentNumber].keys():                
                quantValue,low,high=item['components'][componentNumber][rawfile]
                if (low+high)>1:
                    buff[rawfile]={quantType:(low,high)}
            if buff:
                temp_dict[componentNumber]=buff
        if temp_dict:
            self.tempQuant=temp_dict
        
        self.updateBarPlot()
                    
                    
        


    def createMessageBox(self,output):
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


    def filterLines(self):
        selectedFiles,nSelectedFiles=self.FileSelection.checkedItems()
        for item in nSelectedFiles:
            fileName=self.dict_condition[item]
            self.gvQuantTrace.artist_dict[fileName]['line'].set_visible(False)
            self.gvQuantTrace.artist_dict[fileName]['line'].set_picker(0)
            self.gvQuantTrace.artist_dict[fileName]['fill'].set_visible(False)
            if self.gvQuantTrace.artist_dict[fileName]['markerline']:
                markerline=self.gvQuantTrace.artist_dict[fileName]['markerline']
                markerline.remove()
                del markerline
                self.gvQuantTrace.artist_dict[fileName]['markerline']=None

        
        for item in selectedFiles:
            fileName=self.dict_condition[item]
            self.gvQuantTrace.artist_dict[fileName]['line'].set_visible(True)
            self.gvQuantTrace.artist_dict[fileName]['line'].set_picker(5)
            self.gvQuantTrace.artist_dict[fileName]['fill'].set_visible(True)
            if self.gvQuantTrace.artist_dict[fileName]['markerline']:
                self.gvQuantTrace.artist_dict[fileName]['markerline'].set_visible(True)
        
        self.gvQuantTrace.draw() 
                    
class BarPlotCanvas(FigureCanvas):
    def __init__(self, parent=None, width=9, height=3.5):
        fig = Figure(figsize=(width, height))
        self.ax = fig.add_subplot(111)
        FigureCanvas.__init__(self, fig)
        import sys
        sys._excepthook = sys.excepthook 
        def exception_hook(exctype, value, traceback):
            print(exctype, value, traceback)
            sys._excepthook(exctype, value, traceback) 
            sys.exit(1) 
        sys.excepthook = exception_hook
    
    def updateBarplot(self,DataFrame):
        self.ax.clear()
        ls_mod=list()
        #maximum=0
        for col in DataFrame.columns:
            if col[0]=="%":
                ls_mod.append(col)
                #if DataFrame[col].max()>maximum:
                    #maximum=DataFrame[col].max()

        numberMod=len(set(ls_mod))
        barwidth=0.8/numberMod
        
        if numberMod<1:
            return
        elif numberMod==1:
            mod=ls_mod[0]
            height=DataFrame[mod]
            x=[x for x in range(len(height))]
            self.ax.bar(x,height,width=barwidth,edgecolor="white",align="edge",color=DataFrame['color'])
            for x_pos,y_pos in zip(x,height):
                x_pos=x_pos+barwidth/2                
                self.ax.text(x_pos,y_pos,str(round(y_pos,1)),verticalalignment="bottom",horizontalalignment="center",size="medium")
            self.ax.set_ylabel(mod)
        
        else: 
            colorPalette=sns.color_palette("muted",len(ls_mod))
            for n,mod in enumerate(ls_mod):
                height=DataFrame[mod]
                x=[(x+barwidth*n) for x in range(len(height))]
                self.ax.bar(x,height,width=barwidth,edgecolor="white",label=mod,align="edge",color=colorPalette[n])
                for x_pos,y_pos in zip(x,height):
                    x_pos=x_pos+barwidth/2                
                    self.ax.text(x_pos,y_pos,str(round(y_pos,1)),verticalalignment="bottom",horizontalalignment="center",size="small")
                self.ax.set_ylabel("%modified")
                self.ax.legend()
        
        self.ax.set_xticks([n+0.4 for n in  range(len(height))])            
        self.ax.set_xticklabels(DataFrame['label'])
        ylim=self.ax.get_ylim()[1]
        ylim=ylim+ylim*0.2
        self.ax.set_ylim(0,ylim)
        


        
        self.figure.tight_layout(pad=1, w_pad=0, h_pad=0.5)
        self.draw()

        

            



class LinePlotCanvas(FigureCanvas):
    def __init__(self, parent=None, width=9, height=3.5,alpha=0.3):
        fig = Figure(figsize=(width, height))
        self.ax = fig.add_subplot(111)
        self.alpha=alpha
        
        FigureCanvas.__init__(self, fig)
        import sys
        sys._excepthook = sys.excepthook 
        def exception_hook(exctype, value, traceback):
            print(exctype, value, traceback)
            sys._excepthook(exctype, value, traceback) 
            sys.exit(1) 
        sys.excepthook = exception_hook   
        self.setParent(parent)
        self.setMinimumSize(QtCore.QSize(600, 200))
        self.setMaximumSize(QtCore.QSize(1280, 1280))
        self.setFocusPolicy(QtCore.Qt.ClickFocus)
        self.setFocus()
        self.updateGeometry()
        yfmt = ScalarFormatterForceFormat()
        yfmt.set_powerlimits((0,0))
        self.ax.yaxis.set_major_formatter(yfmt) 
        self.ax.set_xlabel('RT in min')
        self.ax.set_ylabel('intensity')
        self.figure.tight_layout(pad=1, w_pad=0, h_pad=0.5)
        self.text=None
        self.artist_dict={}
        self.span = SpanSelector(self.ax, self.onselect, 'horizontal', useblit=True,
                    rectprops=dict(alpha=0.5, facecolor='red'),button=3)
        self.span.set_active(False)
        self.figure.canvas.mpl_connect('pick_event',self.onpick)
        self.figure.canvas.mpl_connect('key_press_event', self.press)
                     
    def press (self,event):
        print event.key


    def onpick(self,event):
        if not event.mouseevent.button==1:
            return
        thisline=event.artist
        label=thisline.get_label()
        markerline=self.artist_dict[label]['markerline']
    
                
        if markerline:
            markerline.remove()
            del markerline
            self.artist_dict[label]['markerline']=None
    
            if self.text:
                self.text.remove()        
                del self.text
                self.text=None
            self.draw()
            return
        
        if self.text:
            self.text.remove()        
            del self.text
            self.text=None
            
       
        ylow, yhigh=self.ax.get_ylim()
        xlow, xhigh=self.ax.get_xlim()
        xpadding=(xhigh-xlow)*0.02
        ypadding=(yhigh-ylow)*0.02
        textLabel=self.artist_dict[label]['label']
        self.text=self.ax.text((xhigh-xpadding),yhigh-ypadding,textLabel,
                size=9,verticalalignment='top',
                horizontalalignment='right')
        
                    
        xdata=thisline.get_xdata()
        ydata=thisline.get_ydata()
        markerline=self.ax.plot(xdata,ydata,color='yellow',linewidth=8,alpha=.5,label=label,zorder=1)        
        self.artist_dict[label]['markerline']=markerline[0]
        self.draw()    
        return (thisline)

    def onselect(self,xmin, xmax):
        output_list=list()
        for key in self.artist_dict.keys():
            if self.artist_dict[key]['markerline']:
                output_list.append(key)
        
        if not output_list:
            return
        
        for key in output_list:
            fill=self.artist_dict[key]['fill']
            fill.remove()
            del fill
            line=self.artist_dict[key]['line']
            xdata=line.get_xdata()
            ydata=line.get_ydata()
            buff=pd.DataFrame()
            buff['x']=xdata
            buff['y']=ydata
            color=self.artist_dict[key]['color']
            buff=buff[(buff['x']>=xmin)&(buff['x']<=xmax)]
            fill=self.ax.fill_between(buff['x'],buff['y'],alpha=0.3,label=key,color=color)
            self.artist_dict[key]['fill']=fill
            self.artist_dict[key]['limit']=(xmin,xmax)
        self.draw()

    def fill(self,key,limit):
        xmin,xmax=limit
        fill=self.artist_dict[key]['fill']
        fill.remove()
        del fill
        line=self.artist_dict[key]['line']
        xdata=line.get_xdata()
        ydata=line.get_ydata()
        buff=pd.DataFrame()
        buff['x']=xdata
        buff['y']=ydata
        color=self.artist_dict[key]['color']
        buff=buff[(buff['x']>=xmin)&(buff['x']<=xmax)]
        fill=self.ax.fill_between(buff['x'],buff['y'],alpha=self.alpha,label=key,color=color)
        self.artist_dict[key]['fill']=fill
        self.artist_dict[key]['limit']=(xmin,xmax)
        self.draw()


    def changeComponent(self,component,kind,tempQuant=None):
        self.ax.clear()
        for key in component.rawfile_dict.keys():        
            x=component.rawfile_dict[key]['quant'][kind]['RT']
            y=component.rawfile_dict[key]['quant'][kind]['intensity']
            
            if component.rawfile_dict[key]['condition']:
                label=component.rawfile_dict[key]['condition']
            else:
                label=key
            
            if component.rawfile_dict[key]['color']:
                color=component.rawfile_dict[key]['color']
                line=self.ax.plot(x,y,picker=5,label=key,color=color)[0]
                fill=self.ax.fill_between(x,y,alpha=self.alpha,label=key,color=color)
            else:
                line=self.ax.plot(x,y,picker=5,label=key)[0]
                fill=self.ax.fill_between(x,y,alpha=self.alpha,label=key)
                color=line.get_color()
                
                  
            self.artist_dict[key]={'line':line,
                                    'fill':fill,
                                    'markerline':None,
                                    'color':color,
                                    'limit':None,
                                    'label':label}  
            if tempQuant:
                if key in tempQuant.keys():
                    for item in tempQuant[key].keys():
                        if item == kind:
                            limit=tempQuant[key][item]
                            self.fill(key,limit)
    
        yfmt = ScalarFormatterForceFormat()
        yfmt.set_powerlimits((0,0))
        self.ax.yaxis.set_major_formatter(yfmt) 
        self.ax.set_xlabel('RT in min')
        self.ax.set_ylabel('intensity')
        self.figure.tight_layout(pad=1, w_pad=0, h_pad=0.5)
        self.draw()
        
class CheckableComboBox(QtWidgets.QComboBox):
    def __init__(self,input_table):
        super(QtWidgets.QComboBox,self).__init__()
        
    def addItem(self,item):
        super(CheckableComboBox,self).addItem(item)
        item = self.model().item(self.count()-1,0)
        item.setFlags(QtCore.Qt.ItemIsUserCheckable | QtCore.Qt.ItemIsEnabled)
        item.setCheckState(QtCore.Qt.Checked)
    
    def addItems(self,items):
        for item in items:
            self.addItem(item)
        
    def itemChecked(self,index):
        item = self.model().item(index,0)
        return item.checkState() == QtCore.Qt.Checked
    
    def checkedItems(self):
        checkedItems=list()
        uncheckedItems=list()
        for index in range(self.count()):
            item=self.model().item(index,0)
            if item.checkState() == QtCore.Qt.Checked:
                checkedItems.append(item.text())
            else:
                uncheckedItems.append(item.text())
                
        return([checkedItems,uncheckedItems])
        
    
    