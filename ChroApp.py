# -*- coding: utf-8 -*-
"""
Created on Sat Nov 16 10:56:57 2019

@author: KTMS458
"""

from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtWidgets import QMainWindow, QWidget

import numpy as np
from MedImmune.ThermoRawReader import ThermoRawReader
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(style="whitegrid",rc={'grid.color':'0.8','xtick.color': '.0','legend.fontsize': '8','ytick.color':'0','axes.labelcolor': '0.','grid.linewidth':1,'axes.linewidth':1,'axes.edgecolor': '0.0'})

from matplotlib.ticker import ScalarFormatter
from io import BytesIO
from PIL import Image
import pylab
import pandas as pd

class ScalarFormatterForceFormat(ScalarFormatter):
    def _set_format(self,vmin,vmax):  # Override function that finds format to use.
        self.format = "%1.1f"  # Give format here



class Ui_MainWindow(object):
    def setupUi(self, MainWindow,input_table,kind):
        self.DataFrame=input_table        
        MainWindow.setObjectName("Waters UV data extraction")
        MainWindow.resize(850, 612)
        self.centralwidget = QtWidgets.QWidget(MainWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.gridLayout_3 = QtWidgets.QGridLayout(self.centralwidget)
        self.gridLayout_3.setObjectName("gridLayout_3")
        self.formLayout_3 = QtWidgets.QFormLayout()
        self.formLayout_3.setObjectName("formLayout_3")
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(12)
        self.groupBoxLimits = QtWidgets.QGroupBox(self.centralwidget)
        self.groupBoxLimits.setMinimumSize(QtCore.QSize(212, 188))
        self.groupBoxLimits.setMaximumSize(QtCore.QSize(212, 220))
        self.groupBoxLimits.setObjectName("groupBoxLimits")
        self.gridLayout = QtWidgets.QGridLayout(self.groupBoxLimits)
        self.gridLayout.setObjectName("gridLayout")
        spacerItem2 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.gridLayout.addItem(spacerItem2, 1, 0, 1, 1)
        self.RefreshButton = QtWidgets.QPushButton(self.groupBoxLimits)
        self.RefreshButton.setObjectName("RefreshButton")
        self.gridLayout.addWidget(self.RefreshButton, 1, 1, 1, 1)
        self.RefreshButton.clicked.connect(self.refresh)
        
        spacerItem3 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.gridLayout.addItem(spacerItem3, 1, 2, 1, 1)
        self.formLayout = QtWidgets.QFormLayout()
        self.formLayout.setObjectName("formLayout")
        self.labelBegin = QtWidgets.QLabel(self.groupBoxLimits)
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(12)
        self.labelBegin.setFont(font)
        self.labelBegin.setObjectName("labelBegin")
        self.formLayout.setWidget(0, QtWidgets.QFormLayout.LabelRole, self.labelBegin)
        self.leBegin = QtWidgets.QLineEdit(self.groupBoxLimits)
        self.leBegin.setObjectName("leBegin")
        self.formLayout.setWidget(0, QtWidgets.QFormLayout.FieldRole, self.leBegin)
        self.labelBreak1 = QtWidgets.QLabel(self.groupBoxLimits)
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(12)
        self.labelBreak1.setFont(font)
        self.labelBreak1.setObjectName("labelBreak1")
        self.formLayout.setWidget(1, QtWidgets.QFormLayout.LabelRole, self.labelBreak1)
        self.leBreak1 = QtWidgets.QLineEdit(self.groupBoxLimits)
        self.leBreak1.setObjectName("leBreak1")
        self.formLayout.setWidget(1, QtWidgets.QFormLayout.FieldRole, self.leBreak1)
        self.labelBreak2 = QtWidgets.QLabel(self.groupBoxLimits)
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(12)
        self.labelBreak2.setFont(font)
        self.labelBreak2.setObjectName("labelBreak2")
        self.formLayout.setWidget(2, QtWidgets.QFormLayout.LabelRole, self.labelBreak2)
        self.leBreak2 = QtWidgets.QLineEdit(self.groupBoxLimits)
        self.leBreak2.setObjectName("leBreak2")
        self.formLayout.setWidget(2, QtWidgets.QFormLayout.FieldRole, self.leBreak2)
        self.labelBreak3 = QtWidgets.QLabel(self.groupBoxLimits)
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(12)
        self.labelBreak3.setFont(font)
        self.labelBreak3.setObjectName("labelBreak3")
        self.formLayout.setWidget(3, QtWidgets.QFormLayout.LabelRole, self.labelBreak3)
        self.leBreak3 = QtWidgets.QLineEdit(self.groupBoxLimits)
        self.leBreak3.setObjectName("leBreak3")
        self.formLayout.setWidget(3, QtWidgets.QFormLayout.FieldRole, self.leBreak3)
        self.labelEnd = QtWidgets.QLabel(self.groupBoxLimits)
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(12)
        self.labelEnd.setFont(font)
        self.labelEnd.setObjectName("labelEnd")
        self.formLayout.setWidget(4, QtWidgets.QFormLayout.LabelRole, self.labelEnd)
        self.leEnd = QtWidgets.QLineEdit(self.groupBoxLimits)
        self.leEnd.setObjectName("leEnd")
        self.formLayout.setWidget(4, QtWidgets.QFormLayout.FieldRole, self.leEnd)
        self.labelSpacing = QtWidgets.QLabel(self.groupBoxLimits)
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(12)
        self.labelSpacing.setFont(font)
        self.labelSpacing.setObjectName("labelSpacing")
        self.labelSpacing.setText("0")
        self.formLayout.setWidget(5, QtWidgets.QFormLayout.LabelRole, self.labelSpacing)
        self.leSpacing = QtWidgets.QLineEdit(self.groupBoxLimits)
        self.leSpacing.setObjectName("lineEdit")
        self.formLayout.setWidget(5, QtWidgets.QFormLayout.FieldRole, self.leSpacing)
        self.gridLayout.addLayout(self.formLayout, 0, 0, 1, 3)
        self.gridLayout_3.addWidget(self.groupBoxLimits, 2, 1, 1, 1)
        self.widget = QtWidgets.QWidget(self.centralwidget)
        self.widget.setMinimumSize(QtCore.QSize(0, 40))
        self.widget.setMaximumSize(QtCore.QSize(16777215, 40))
        self.widget.setObjectName("widget")
        self.gridLayout_3.addWidget(self.widget, 0, 0, 1, 1)
        self.label_3 = QtWidgets.QLabel(self.centralwidget)
        self.label_3.setMaximumSize(QtCore.QSize(212, 34))
        self.label_3.setText("")
        self.label_3.setObjectName("label_3")
        self.gridLayout_3.addWidget(self.label_3, 4, 1, 1, 1)
        self.frame = QtWidgets.QFrame(self.centralwidget)


        self.initFig(input_table,kind)
        self.canvas=FigureCanvas(self.fig) 
        self.canvas.setFocusPolicy(QtCore.Qt.ClickFocus )
        self.toolbar = NavigationToolbar(self.canvas, self.widget, coordinates=True)

        self.gridLayout_3.addWidget(self.toolbar, 0, 0, 1, 1)
        self.gridLayout_3.addWidget(self.canvas, 1, 0, 9, 1)
        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 850, 21))
        self.menubar.setObjectName("menubar")
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(MainWindow)
        self.statusbar.setObjectName("statusbar")
        MainWindow.setStatusBar(self.statusbar)

        self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "Chromatogram Traces"))
        self.groupBoxLimits.setTitle(_translate("MainWindow", "Limits"))
        self.RefreshButton.setText(_translate("MainWindow", "Refresh"))
        self.labelBegin.setText(_translate("MainWindow", "Begin"))
        self.labelBreak1.setText(_translate("MainWindow", "Break1"))
        self.labelBreak2.setText(_translate("MainWindow", "Break2"))
        self.labelBreak3.setText(_translate("MainWindow", "Break3"))
        self.labelEnd.setText(_translate("MainWindow", "End"))
        self.labelSpacing.setText(_translate("MainWindow", "Spacing"))


        
    def initFig(self,input_table,kind):
        ls_intensities=list()
        ls_RT=list()
        if kind=="UV":
            for n,row in input_table.iterrows():
                location=row['Location']
                rawfile=ThermoRawReader(location)
                df_UV=rawfile.getUV()
                ls_intensities.append(list(df_UV['intensity']))
                ls_RT.append(list(df_UV['RT']))
                rawfile.close()
        elif kind=="TIC":
            for n,row in input_table.iterrows():
                location=row['Location']
                rawfile=ThermoRawReader(location)
                df_TIC=rawfile.getTIC()
                ls_intensities.append(list(df_TIC['intensity']))
                ls_RT.append(list(df_TIC['RT']))
                rawfile.close()        
        elif kind=="BPC":
            for n,row in input_table.iterrows():
                location=row['Location']
                rawfile=ThermoRawReader(location)
                df_BPC=rawfile.getTIC()
                ls_intensities.append(list(df_BPC['intensity']))
                ls_RT.append(list(df_BPC['RT']))
                rawfile.close() 
        else:
            raise ValueError("Paramater kind must be either UV, TIC or BPC")
        
        self.DataFrame['RT']=ls_RT
        self.DataFrame['intensity']=ls_intensities
        
        self.fig=plt.figure(figsize=(9,3))
        self.ax=self.fig.add_subplot(3,1,1)
        plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        plt.tight_layout(pad=0.5, w_pad=0, h_pad=0)
        self.ax2=self.fig.add_subplot(3,2,3)

        self.ax3=self.fig.add_subplot(3,2,4)
        self.ax4=self.fig.add_subplot(3,2,5)
        self.ax5=self.fig.add_subplot(3,2,6)
        self.axDict={'ax2':{'ax':self.ax2,'lowerLimit':0,'higherLimit':0},
                     'ax3':{'ax':self.ax3,'lowerLimit':0,'higherLimit':0},
                     'ax4':{'ax':self.ax4,'lowerLimit':0,'higherLimit':0},
                     'ax5':{'ax':self.ax5,'lowerLimit':0,'higherLimit':0}}
        
        plt.tight_layout(pad=0.5, w_pad=2, h_pad=1)
                

        
        for n,row in self.DataFrame.iterrows():
            RT=row['RT']
            intensity=row['intensity']
            color=row['color']
            self.ax.plot(RT,intensity,linewidth=1.0,color=color)
                
          
        for item in ([self.ax.title, self.ax.xaxis.label, self.ax.yaxis.label] +
                     self.ax.get_xticklabels() + self.ax.get_yticklabels()):
            item.set_fontsize(12)
        
        self.ax.set_ylabel("AU")
        self.ax.set_xlabel('')

    
    def refresh(self):  
        self.Graph1lowerLimit=float(self.leBegin.text())
        self.Graph1higherLimit=float(self.leBreak1.text())
        self.axDict['ax2']['lowerLimit']=self.Graph1lowerLimit
        self.axDict['ax2']['higherLimit']=self.Graph1higherLimit        
        self.Graph2lowerLimit=float(self.leBreak1.text())
        self.Graph2higherLimit=float(self.leBreak2.text())
        self.axDict['ax3']['lowerLimit']=self.Graph2lowerLimit
        self.axDict['ax3']['higherLimit']=self.Graph2higherLimit
        self.Graph3lowerLimit=float(self.leBreak2.text())
        self.Graph3higherLimit=float(self.leBreak3.text())  
        self.axDict['ax4']['lowerLimit']=self.Graph3lowerLimit
        self.axDict['ax4']['higherLimit']=self.Graph3higherLimit
        self.Graph4lowerLimit=float(self.leBreak3.text())
        self.Graph4higherLimit=float(self.leEnd.text())
        self.axDict['ax5']['lowerLimit']=self.Graph4lowerLimit
        self.axDict['ax5']['higherLimit']=self.Graph4higherLimit
        
        try:self.spacing=float(self.leSpacing.text())
        except:self.spacing=0
        highest_value=0
        lowest_value=0
        
        yfmt = ScalarFormatterForceFormat()
        yfmt.set_powerlimits((0,0))
        
        m=0        
        for key in self.axDict: 
            ax=self.axDict[key]['ax']
            ax.clear()
        
        for n,row in self.DataFrame.iterrows():
            color=row['color']
            intensity=np.array(row['intensity'])
            RT=np.array(row['RT'])
            intensity=intensity+self.spacing*m            
            signal=np.dstack((RT,intensity))[0]
            tempMax=max(signal[:,1][(signal[:,0]>20)&(signal[:,0]<70)])
            tempMin=min(signal[:,1][(signal[:,0]>20)&(signal[:,0]<70)])
            if tempMax>highest_value:
                highest_value=tempMax
            if tempMin<lowest_value:
                lowest_value=tempMin
            m=m+1                    
            for key in self.axDict:
                ax=self.axDict[key]['ax']
                lowerLimit=self.axDict[key]['lowerLimit']
                higherLimit=self.axDict[key]['higherLimit']                
                tempSignal=signal[(signal[:,0]>lowerLimit)&(signal[:,0]<higherLimit)]
                tempRT=tempSignal[:,0]
                tempIntensity=tempSignal[:,1]                                                                
                ax.plot(tempRT,tempIntensity,linewidth=1.0,color=color) 
                ax.yaxis.set_major_formatter(yfmt)  
                plt.tight_layout(pad=0.5, w_pad=2, h_pad=0.5)

            
        #self.DataFrame['tempIntensities']=ls_tempIntensities
        self.lowest_value=lowest_value-0.05*highest_value
        self.highest_value=highest_value+0.05*highest_value
        self.ax2.set_ylim((self.lowest_value,self.highest_value)) 
        self.ax3.set_ylim((self.lowest_value,self.highest_value)) 
        self.ax4.set_ylim((self.lowest_value,self.highest_value)) 
        self.ax5.set_ylim((self.lowest_value,self.highest_value))

        self.canvas.draw()
        
     




class ChroApp(QMainWindow):
    def __init__(self,input_table,kind="UV"):
        super(QWidget,self).__init__()
        self.ui = Ui_MainWindow()
        self.ui.setupUi(self,input_table,kind) 
        

    def getOutputTables(self):
        sns.set_style("ticks")
        yfmt = ScalarFormatterForceFormat()
        yfmt.set_powerlimits((0,0))
        ls_images=list()
        ls_conditions=list()
        ls_lines=list()
        lowest_limit=self.ui.Graph1higherLimit
        highest_limit=self.ui.Graph4higherLimit
        m=0
        fig,ax=plt.subplots(1,figsize=(9,3))
        for n,row in self.ui.DataFrame.iterrows():
            color=row['color']
            condition=row['label']
            ls_conditions.append(condition)
            intensity=np.array(row['intensity'])
            RT=np.array(row['RT'])
            intensity=intensity+self.ui.spacing*m
            signal=np.dstack((RT,intensity))[0]
            m=m+1                            
            tempSignal=signal[(signal[:,0]>lowest_limit)&(signal[:,0]<highest_limit)]
            tempRT=tempSignal[:,0]
            tempIntensity=tempSignal[:,1]                                                                
            lines=ax.plot(tempRT,tempIntensity,linewidth=1.0,color=color) 
            ls_lines.extend(lines)
        
        ax.yaxis.set_major_formatter(yfmt)          
        ax.set_ylabel("intensity")
        ax.set_xlabel('time')
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.set_ylim((self.ui.lowest_value,self.ui.highest_value))
        plt.tight_layout(pad=0.5, w_pad=2, h_pad=0.5)
        fig_buffer=BytesIO()
        fig.savefig(fig_buffer,format='png',dpi=150)
        img=Image.open(fig_buffer)
        ls_images.append(img)



        for key in self.ui.axDict.keys():
            fig,ax=plt.subplots(1,figsize=(9,3))
            lowerLimit=self.ui.axDict[key]['lowerLimit']
            higherLimit=self.ui.axDict[key]['higherLimit']
            m=0
            for n,row in self.ui.DataFrame.iterrows():
                color=row['color']                
                intensity=np.array(row['intensity'])
                RT=np.array(row['RT'])
                intensity=intensity+self.ui.spacing*m            
                signal=np.dstack((RT,intensity))[0]                                                   
                tempSignal=signal[(signal[:,0]>lowerLimit)&(signal[:,0]<higherLimit)]
                tempRT=tempSignal[:,0]
                tempIntensity=tempSignal[:,1]                                                                
                ax.plot(tempRT,tempIntensity,linewidth=1.0,color=color)
                m=m+1
            ax.yaxis.set_major_formatter(yfmt)  
            ax.set_ylabel("intensity")
            ax.set_xlabel('time')
            ax.spines['right'].set_visible(False)
            ax.spines['top'].set_visible(False)
            ax.set_ylim((self.ui.lowest_value,self.ui.highest_value))
            plt.tight_layout(pad=0.5, w_pad=2, h_pad=0.5)
            fig_buffer=BytesIO()
            fig.savefig(fig_buffer,format='png',dpi=150)
            img=Image.open(fig_buffer)
            ls_images.append(img)
            

        figlegend = pylab.figure(figsize=(1.5,1))
        figlegend.legend(handles=ls_lines,labels=ls_conditions,loc='center')
        fig_buffer=BytesIO()
        pylab.tight_layout(pad=0.5, w_pad=0, h_pad=0.5)
        figlegend.savefig(fig_buffer,format='png',dpi=300, bbox_inches="tight",transparent = True, pad_inches=0)
        img=Image.open(fig_buffer)
        ls_images.append(img)
        
        output_table=pd.DataFrame(data=ls_images,columns=['images'],index=['all','firstQuartil','secondQuartil','thirdQuartil','fourthQuartil','legend'])
        return(output_table)
    

