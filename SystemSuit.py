# -*- coding: utf-8 -*-
"""
Created on Mon Sep 23 07:30:20 2019

@author: KTMS458
"""
import pandas as pd
from MedImmune.ThermoRawReader import ThermoRawReader as RawReader
import seaborn as sns
from matplotlib.ticker import ScalarFormatter

sns.set(style="whitegrid",rc={'grid.color':'0.8','xtick.color': '.0','ytick.color':'0','axes.labelcolor': '0.','grid.linewidth':0.5,'axes.linewidth':1,'axes.edgecolor': '0.0'})
import matplotlib.pyplot as plt
import numpy as np
from mMass.mspy import mod_peakpicking as peak_picking
import os


dir_pepmix =os.path.abspath(os.path.dirname(__file__))
path_pepmix=os.path.join(dir_pepmix,'CalibrationStandards.xlsx')
class ScalarFormatterForceFormat(ScalarFormatter):
    def _set_format(self,vmin,vmax):  # Override function that finds format to use.
        self.format = "%1.1f"  # Give format here



class SystemSuit():
    def __init__(self,path):
        self.input_table=pd.read_excel(path_pepmix)
        self.path=path
        self.rawfile=RawReader(path)
        self.fig=plt.figure(figsize=(9,3.5))
        self.ax=plt.subplot(111)
        #self.ax2=plt.subplot(212)

    def __del__(self):
        self.rawfile.close()        

    def getAttributes(self):
        outputRows=list()
        for n,row in self.input_table.iterrows():
            peptideNumber=row['Peptide Number']
            MZ_low=row['MZ_low']
            MZ_high=row['MZ_high']
            massrange="%f-%f"%(MZ_low,MZ_high)
            dfXIC=self.rawfile.getXIC(massrange,8.5,25,smooth=False)
            try:
                signal=np.dstack((dfXIC['RT'],dfXIC['intensity']))[0]
                peaklist=peak_picking.labelscan(signal,relThreshold=0.01)
            except:
                continue
            if len(peaklist)<1:
                continue
            peakHighestIntensity=peaklist[0]        
            for peak in peaklist:
                if peak.ai>peakHighestIntensity.ai:
                    peakHighestIntensity=peak
            RT=peakHighestIntensity.mz
            intensity=peakHighestIntensity.ai
            fwhm=peakHighestIntensity.fwhm
            MSspectra=self.rawfile.getMS1spectraByRT(RT)
            MSspectra=MSspectra[(MSspectra['mz']<MZ_high) &  (MSspectra['mz']>MZ_low)]
            mz=MSspectra['mz'][MSspectra['intensity'].idxmax()]
            
            outputRow=[peptideNumber,RT,intensity,fwhm,mz]
            outputRows.append(outputRow)
        output_table=pd.DataFrame(data=outputRows,columns=["Peptide Number","RT","intensity","fwhm","MZ_observed"])
        output_table['relIntensity']=output_table['intensity']/(output_table['intensity'].max())
        self.input_table=pd.merge(self.input_table,output_table,on=['Peptide Number'])
            
  
    def createLinePlot(self):
        for n,row in self.input_table.iterrows():
            peptideNumber=row['Peptide Number']
            MZ_low=row['MZ_low']
            MZ_high=row['MZ_high']
            massrange="%f-%f"%(MZ_low,MZ_high)
            dfXIC=self.rawfile.getXIC(massrange,8.5,25,smooth=False)
            self.ax.plot(dfXIC['RT'],dfXIC['intensity'],linewidth=1.0,label=peptideNumber)
            continue
           
            if peptideNumber in [1,2,3,4,5,6,7,8]:
                dfXIC=self.rawfile.getXIC(massrange,8.5,25,smooth=False)
                self.ax.plot(dfXIC['RT'],dfXIC['intensity'],linewidth=1.0,label=peptideNumber)
            if peptideNumber in [8,9,10,11,12,13,14,15]:
                dfXIC=self.rawfile.getXIC(massrange,8.5,25,smooth=False)
                self.ax2.plot(dfXIC['RT'],dfXIC['intensity'],linewidth=1.0,label=peptideNumber)  
            self.ax.legend()
            #self.ax2.legend()
    
    def creatStackedLinePlot(self):
        fig, axes = plt.subplots(7, 2,figsize=(9,14),sharex=True)
        n=0
        for i in range(2):
            for ax in axes:
                MZ_low=self.input_table['MZ_low'][n]
                MZ_high=self.input_table['MZ_high'][n]
                massrange="%f-%f"%(MZ_low,MZ_high)
                dfXIC=self.rawfile.getXIC(massrange,8.5,25,smooth=False)
                ax[i].plot(dfXIC['RT'],dfXIC['intensity'])
                n=n+1
        
                

    