# -*- coding: utf-8 -*-
"""
Created on Tue Nov 19 14:22:24 2019

@author: KTMS458
"""

from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtWidgets import QApplication, QMainWindow, QWidget
import sys
import re

import pylab


from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
import matplotlib.pyplot as plt
import seaborn as sns
from MedImmune.ThermoRawReader import ThermoRawReader
sns.set(style="whitegrid",rc={'grid.color':'0.8','xtick.color': '.0','legend.fontsize': '8','ytick.color':'0','axes.labelcolor': '0.','grid.linewidth':1,'axes.linewidth':1,'axes.edgecolor': '0.0'})

from io import BytesIO
from PIL import Image
from matplotlib.ticker import ScalarFormatter

import pandas as pd

input_table=pd.read_hdf(r"D:\Vahid\PepModQuantTestfiles\input_table_1.h5")


class ScalarFormatterForceFormat(ScalarFormatter):
    def _set_format(self,vmin,vmax):  # Override function that finds format to use.
        self.format = "%1.1f"  # Give format here


class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
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

        self.DataFrame=input_table
        self.initFig()
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
        MainWindow.setWindowTitle(_translate("MainWindow", "MainWindow"))
        self.groupBoxLimits.setTitle(_translate("MainWindow", "Limits"))
        self.RefreshButton.setText(_translate("MainWindow", "Refresh"))
        self.labelBegin.setText(_translate("MainWindow", "Begin"))
        self.labelBreak1.setText(_translate("MainWindow", "Break1"))
        self.labelBreak2.setText(_translate("MainWindow", "Break2"))
        self.labelBreak3.setText(_translate("MainWindow", "Break3"))
        self.labelEnd.setText(_translate("MainWindow", "End"))
        self.labelSpacing.setText(_translate("MainWindow", "Spacing"))


        
    def initFig(self):
        self.fig=plt.figure(figsize=(9,3))
        self.ax=self.fig.add_subplot(3,1,1)
        plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        plt.tight_layout(pad=0.5, w_pad=0, h_pad=0)
        self.ax2=self.fig.add_subplot(3,2,3)

        self.ax3=self.fig.add_subplot(3,2,4)
        self.ax4=self.fig.add_subplot(3,2,5)
        self.ax5=self.fig.add_subplot(3,2,6)
        plt.tight_layout(pad=0.5, w_pad=2, h_pad=1)
                
        #self.rawfile_dict=dict()
        for n,row in self.DataFrame.iterrows():
            location=row['Location']
            color=row['color']
            rawfile=ThermoRawReader(location)
            dfUV=rawfile.getUV()
            self.ax.plot(dfUV['RT'],dfUV['intensity'],linewidth=1.0,color=color)
            rawfile.close()
            #self.rawfile_dict[location]=dfUV
                
          
        for item in ([self.ax.title, self.ax.xaxis.label, self.ax.yaxis.label] +
                     self.ax.get_xticklabels() + self.ax.get_yticklabels()):
            item.set_fontsize(12)
        
        self.ax.set_ylabel("AU")
        self.ax.set_xlabel('')

    
    def refresh(self):  
        self.Graph1lowerLimit=float(self.leBegin.text())
        self.Graph1higherLimit=float(self.leBreak1.text())
        self.Graph2lowerLimit=float(self.leBreak1.text())
        self.Graph2higherLimit=float(self.leBreak2.text())
        self.Graph3lowerLimit=float(self.leBreak2.text())
        self.Graph3higherLimit=float(self.leBreak3.text())        
        self.Graph4lowerLimit=float(self.leBreak3.text())
        self.Graph4higherLimit=float(self.leEnd.text())
        try:self.spacing=float(self.leSpacing.text())
        except:self.spacing=0
        highest_value=0
        lowest_value=0
        
        yfmt = ScalarFormatterForceFormat()
        yfmt.set_powerlimits((0,0))
        

        self.ax2.clear()
        self.ax3.clear()
        self.ax4.clear()
        self.ax5.clear()
        
        
        for n,row in input_table.iterrows():
            location=row['Location']
            color=row['color']
            rawfile=ThermoRawReader(location)
            df_UV=rawfile.getUV()
            df_buff=df_UV[(df_UV['RT']>self.Graph1lowerLimit)&(df_UV['RT']<self.Graph1higherLimit)]
            df_buff['intensity']=df_buff['intensity']+n*self.spacing
            self.ax2.plot(df_buff['RT'],df_buff['intensity'],linewidth=1.0,color=color)
            self.ax2.yaxis.set_major_formatter(yfmt)  
            plt.tight_layout(pad=0.5, w_pad=2, h_pad=0.5)
            df_buff=df_UV[(df_UV['RT']>self.Graph2lowerLimit)&(df_UV['RT']<self.Graph2higherLimit)]
            df_buff['intensity']=df_buff['intensity']+n*self.spacing
            self.ax3.plot(df_buff['RT'],df_buff['intensity'],linewidth=1.0,color=color)
            self.ax3.yaxis.set_major_formatter(yfmt)  
            plt.tight_layout(pad=0.5, w_pad=2, h_pad=0.5)
            df_buff=df_UV[(df_UV['RT']>self.Graph3lowerLimit)&(df_UV['RT']<self.Graph3higherLimit)]
            df_buff['intensity']=df_buff['intensity']+n*self.spacing
            self.ax4.plot(df_buff['RT'],df_buff['intensity'],linewidth=1.0,color=color)
            self.ax4.yaxis.set_major_formatter(yfmt)  
            plt.tight_layout(pad=0.5, w_pad=2, h_pad=0.5)
            df_buff=df_UV[(df_UV['RT']>self.Graph4lowerLimit)&(df_UV['RT']<self.Graph4higherLimit)]
            df_buff['intensity']=df_buff['intensity']+n*self.spacing
            self.ax5.plot(df_buff['RT'],df_buff['intensity'],linewidth=1.0,color=color)
            self.ax5.yaxis.set_major_formatter(yfmt)  
            plt.tight_layout(pad=0.5, w_pad=2, h_pad=0.5)
        
        
        


        self.ax2.set_ylim((lowest_value,highest_value)) 
        self.ax3.set_ylim((lowest_value,highest_value)) 
        self.ax4.set_ylim((lowest_value,highest_value)) 
        self.ax5.set_ylim((lowest_value,highest_value))
        self.lowest_value=lowest_value-0.05*highest_value
        self.highest_value=highest_value+0.05*highest_value
        self.canvas.draw()
        
     




class AppWindow(QMainWindow):
    def __init__(self):
        super(QWidget,self).__init__()
        self.ui = Ui_MainWindow()
        self.ui.setupUi(self)
        self.show()  
    
    def getFigParameter(self):
        return ((self.ui.Graph1lowerLimit,self.ui.Graph1higherLimit),
                (self.ui.Graph2lowerLimit,self.ui.Graph2higherLimit),
                (self.ui.Graph3lowerLimit,self.ui.Graph3higherLimit),
                (self.ui.Graph4lowerLimit,self.ui.Graph4higherLimit),
                self.ui.spacing,
                self.ui.lowest_value,
                self.ui.highest_value)




app = QApplication(sys.argv)
w = AppWindow()
w.show()
app.exec_()
"""
sns.set_style("ticks")
ls_images=list()
figParameter=w.getFigParameter()
limits=figParameter[:4]
spacing=figParameter[-3]
lowest_value=figParameter[-2]-0.05*figParameter[-1]
highest_value=figParameter[-1]+0.05*figParameter[-1]
ls_columns=[x for x in input_table_1.columns if re.match(".+raw",x)]
for limit in limits:
    df_buff=input_table_1[(input_table_1['RT']>limit[0])&(input_table_1['RT']<limit[1])]
    time=df_buff['RT']    
    fig,ax=plt.subplots(1,figsize=(9,3))
    for n,item in enumerate(ls_columns):
       buff=df_buff[item]+n*spacing
       ax.plot(time,buff,linewidth=1.0,color=input_table_2['color'][n])
    ax.set_ylabel("intensity")
    ax.set_xlabel('time')
    ax.set_ylim((lowest_value,highest_value))
    #ax.legend(loc='upper right',ncol=len(input_table_2['condition'])/2)
    yfmt = ScalarFormatterForceFormat()
    yfmt.set_powerlimits((0,0))
    ax.yaxis.set_major_formatter(yfmt)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)     
    plt.tight_layout(pad=0.5, w_pad=0, h_pad=0.5)
    fig_buffer=BytesIO()
    fig.savefig(fig_buffer,format='png',dpi=150)
    img=Image.open(fig_buffer)
    ls_images.append(img)

lowest_limit=limits[0][0]
highest_limit=limits[-1][-1]
df_buff=input_table_1[(input_table_1['RT']>lowest_limit)&(input_table_1['RT']<highest_limit)]
time=df_buff['RT']
fig,ax=plt.subplots(1,figsize=(9,3))
ls_lines=list()
ls_conditions=list()
for n,item in enumerate(ls_columns):
   buff=df_buff[item]+n*spacing
   condition=input_table_2['label'][n]
   lines=ax.plot(time,buff,linewidth=1.0,color=input_table_2['color'][n])
   ls_lines.extend(lines)
   ls_conditions.append(condition)
ax.set_ylabel("intensity")
ax.set_xlabel('time')
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False) 
ax.set_ylim((lowest_value,highest_value))
yfmt = ScalarFormatterForceFormat()
yfmt.set_powerlimits((0,0))
ax.yaxis.set_major_formatter(yfmt)   
plt.tight_layout(pad=0.5, w_pad=0, h_pad=0.5)

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


output_table_1=pd.DataFrame()
output_table_2=pd.DataFrame()
output_table_1['output']=ls_images
output_table_2['output'] = [img]"""
