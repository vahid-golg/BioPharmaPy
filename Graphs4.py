# -*- coding: utf-8 -*-
"""
Created on Sat May 11 10:09:31 2019

@author: KTMS458
"""

import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.ticker import ScalarFormatter
import numpy as np


color_dict={'y':'#CF3234','b':'#3E7EB1','by':'#FF00FF'}
color=("#CF3234","#3E7EB1")
sns.set(style="whitegrid",rc={'grid.color':'0.8','xtick.color': '.0','ytick.color':'0','axes.labelcolor': '0.','grid.linewidth':0.5,'axes.linewidth':1,'axes.edgecolor': '0.0'},palette=color)

class ScalarFormatterForceFormat(ScalarFormatter):
    def _set_format(self,vmin,vmax):  # Override function that finds format to use.
        self.format = "%1.1f"  # Give format here


class Graphs():
    def __init__(self,component=None):
        self.fig=plt.figure(figsize=(9,3.5))
        self.ax=plt.subplot(111)
        self.component=component

    def getEmptyFig(self,ax,text="",x_min=0,x_max=1):
        x=x_min+(x_max-x_min)*0.5
        ax.set_ylim(0,1)
        if text:
            ax.annotate(text,(x,0.5),horizontalalignment='center',size=12,verticalalignment='center')
        ax.set_xlim(x_min,x_max)
        return(ax)

    def createMS1Fig(self,MS1Spectrum,fileName,title=""):
        if len(MS1Spectrum)>0:
            x_data_min,x_data_max=self.getMSspectraLimits(self.component.rawfile_dict[fileName]['matchedIsotopologues'])
            MS1Spectrum=MS1Spectrum[(MS1Spectrum['mz']>x_data_min) & (MS1Spectrum['mz']<x_data_max)]
            x_data,y_data=np.array(MS1Spectrum['mz']),np.array(MS1Spectrum['intensity'])
            compMSLabel=self.component.rawfile_dict[fileName]['matchedIsotopologues']
            ydata_max=self.component.rawfile_dict[fileName]['highestIntensityIsotopologue']
            ydata_max=ydata_max+0.15*ydata_max
            self.ax.plot(x_data,y_data,linewidth=1.0)
            self.ax.set_ylim(0,ydata_max)
        
            for item in compMSLabel:

                mz=item[0]
                if not mz:continue
                intensity=item[1]
                text='{:.4f}'.format(round(mz,4))
                self.ax.annotate(text,(mz,intensity+5),horizontalalignment='center',size=9,verticalalignment='bottom')            
            self.ax.set_ylabel('intensity')
        
            ylow, yhigh=self.ax.get_ylim()
            xlow, xhigh=self.ax.get_xlim()

 
        
            yfmt = ScalarFormatterForceFormat()
            yfmt.set_powerlimits((0,0))
            self.ax.yaxis.set_major_formatter(yfmt) 
            self.ax.set_xlabel('mz')
        else:
            self.ax=self.getEmptyFig(self.ax,text="No Signal")                
        
        if title:
            self.ax.set_title(title)
        plt.tight_layout(pad=1, w_pad=0, h_pad=0.5)
        plt.close()

    def getMS1subplot(self,component,ax,xmin,xmax,chopped=False):
        MSspectra=component.rawfile.getMS1spectraByRT(component.realRT)
        MSspectra=MSspectra[(MSspectra['mz']>xmin) & (MSspectra['mz']<xmax)]
        x_data,y_data=np.array(MSspectra['mz']),np.array(MSspectra['intensity'])
        compMSLabel=component.matchedIsotopologues
        ydata_max=component.highestIntensityIsotopologues
        
        ydata_max=ydata_max+0.15*ydata_max
        ax.plot(x_data,y_data,linewidth=1.0)
        ax.set_ylim(0,ydata_max)

        
        for item in compMSLabel:
            mz=item[0]
            intensity=item[1]
            text='{:.4f}'.format(round(mz,4))
            ax.annotate(text,(mz,intensity+5),horizontalalignment='center',size=9,verticalalignment='bottom')            
        ax.set_ylabel('intensity')
        
        ylow, yhigh=ax.get_ylim()
        xlow, xhigh=ax.get_xlim()
        xpadding=(xhigh-xlow)*0.02
        ypadding=(yhigh-ylow)*0.02
        ax.text((xhigh-xpadding),yhigh-ypadding,"file name = %s\nscan number = %i"%(component.rawfile.baseName,
                component.scanNumberHighestMS1Intensity),
                size=9,verticalalignment='top',
                horizontalalignment='right')
        if chopped:        
            plt.setp(ax.get_xticklabels(),visible=False)
            ax.set_yticks(ax.get_yticks()[1:])        
        
        yfmt = ScalarFormatterForceFormat()
        yfmt.set_powerlimits((0,0))
        ax.yaxis.set_major_formatter(yfmt) 
        return(ax)
            
        
    def getMSspectraLimits(self,cluster,refCluster=None):
        theoreticalIsotopes=self.component.theoIsotopologues
        #apply threshold
        theoreticalIsotopes=theoreticalIsotopes[theoreticalIsotopes[:,1]>0.001]
        cluster_xmin,cluster_xmax=min(theoreticalIsotopes[:,0]),max(theoreticalIsotopes[:,0])
        cluster_xdiffrence=cluster_xmax-cluster_xmin
        return(cluster_xmin-0.1*cluster_xdiffrence,cluster_xmax+0.1*cluster_xdiffrence)
        
        
    def getXICFig(self,xic_df,title="",legend=False):
        try:
            self.ax.cla()
        except:
            self.reset(referenceComponent=None)
        for n,row in xic_df.iterrows():
            x_data=row['XIC_RT']
            y_data=row['XIC_intens']
            color=row['color']
            label=row['label']
            self.ax.plot(x_data,y_data,linewidth=1.0,color=color,label=label)
        self.ax.legend(loc='best')
        self.ax.set_ylabel("intensity")
        #title=input_table['Peptide Sequence'][0]+', '+input_table['Site'][0]
        if title:
             self.ax.set_title(title,size=10)
        self.ax.set_xlabel('time')
    
        yfmt = ScalarFormatterForceFormat()
        yfmt.set_powerlimits((0,0))
        self.ax.yaxis.set_major_formatter(yfmt)   
        plt.tight_layout(pad=0.5, w_pad=2, h_pad=2)
        plt.close()

    def getMS2subplot(self,MS2spectrum,fragments,ax,relThreshould=0.05,label=""):
        maxIntens=max(MS2spectrum['mz'])
        x_data=MS2spectrum['mz']
        y_data=MS2spectrum['intensity']
        
        markerline, stemlines, baseline=ax.stem(x_data,y_data,linefmt='.5', markerfmt=' ', basefmt=' ')
        plt.setp(stemlines, 'linewidth', 1.5)
    
        for series in set(fragments['fragmentSerie']):
            color=color_dict[series]
            buff=fragments[fragments['fragmentSerie']==series]
            buff=buff[buff['ionMatchIntens']>(maxIntens*relThreshould)]
            if len(buff)==0:continue
            x_data=np.array(buff['ionMatches'])
            y_data=np.array(buff['ionMatchIntens'])
            markerline, stemlines, baseline=ax.stem(x_data,y_data,linefmt=color, markerfmt=' ', basefmt=' ')
            plt.setp(stemlines, 'linewidth',1.5)
            for n,ion in enumerate(buff['ionName']):
                ax.annotate(ion,(x_data[n],y_data[n]),horizontalalignment='center',size=10,verticalalignment='bottom',color=color)
            
            yfmt = ScalarFormatterForceFormat()
            yfmt.set_powerlimits((0,0))
            ax.yaxis.set_major_formatter(yfmt)  



            if ax.get_yticks()[0]==0:
                ax.set_yticks(ax.get_yticks()[1:])
            ax.text(0.98,0.95,label,verticalalignment='top',
                horizontalalignment='right',
                transform=ax.transAxes)
        minMS2Spectrum=min(MS2spectrum['mz'])-(min(MS2spectrum['mz'])*0.1)
        maxMS2Spectrum=max(MS2spectrum['mz'])+(max(MS2spectrum['mz'])*0.1)
        ax.set_xlim(minMS2Spectrum,maxMS2Spectrum)
        y_max=0.1*max(MS2spectrum['intensity'])+max(MS2spectrum['intensity'])
        ax.set_ylim(0,y_max)
        
        return(ax)
    
    def getMS2fig(self,MS2spectra,fragments,title=""):
        if len(MS2spectra)==0:
            self.ax=self.getEmptyFig(self.ax)
            plt.tight_layout(pad=1, w_pad=0, h_pad=0.5)
            plt.close()
            print "empty"
            return
            
        mz_values= MS2spectra['mz']
        if len(mz_values)==0:
            self.ax=self.getEmptyFig(self.ax,flipped=False,chopped=True)  
        
        self.ax=self.getMS2subplot(MS2spectra,fragments,self.ax)        
        self.ax.set_xlabel('mz')
        self.ax.set_ylabel('intensity')
        if title:
            self.ax.set_title(title)
        plt.tight_layout(pad=1, w_pad=0, h_pad=0.5)
        plt.close()




        
    def reset(self,referenceComponent):
        self.__init__(referenceComponent=referenceComponent)
        
