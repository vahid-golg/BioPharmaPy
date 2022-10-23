# -*- coding: utf-8 -*-
import bisect
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from  mMass.mspy.mod_proteo import fragment as fragment
from  mMass.mspy.mod_proteo import fragmentlosses as fragmentlosses
import mMass.mspy.blocks as blocks
import re
import xml.dom.minidom as xml
import os
import comtypes
from comtypes.client import GetModule, CreateObject
#from comtypes.gen import MSFileReaderLib
from ctypes import c_long,c_double,byref
#import peakutils
from collections import OrderedDict
#from mpl_toolkits.axisartist.anchored_artists import AnchoredText
from matplotlib.offsetbox import AnchoredText
from matplotlib.ticker import ScalarFormatter

class ScalarFormatterForceFormat(ScalarFormatter):
    def _set_format(self,vmin,vmax):  # Override function that finds format to use.
        self.format = "%1.1f"  # Give format here



color=("#CF3234","#3E7EB1")
sns.set(style="whitegrid",rc={'grid.color':'0.8','xtick.color': '.0','ytick.color':'0','axes.labelcolor': '0.','grid.linewidth':0.5,'axes.linewidth':1,'axes.edgecolor': '0.0'},palette=color)
blocks.loadModifications(r"C:\ProgramData\Anaconda2\Lib\site-packages\mMass\configs\modifications.xml")
color_dict={'y':'#CF3234','b':'#3E7EB1','by':'#FF00FF'}


def getMS2spectra(scan_df):
    rows=list()
    indexes=scan_df.index
    for ind in indexes:
        path=str(ind)
        scanNumber=int(scan_df.loc[ind])
        print(scanNumber)
        print(path)
        if scanNumber==-1:
            rows.append([path,[],[]])
            print("hier")
        else:
            mz_values,int_values=getMS2Label(path,scanNumber=scanNumber)               
            rows.append([path,mz_values,int_values])

    return (pd.DataFrame(rows, columns=['Location', 'raw_mz','raw_int']))





def movaver(signal, window, cycles=1, style='flat'):
    """Smooth signal by moving average filter. New array is returned.
        signal (numpy array) - signal data points
        window (float) - m/z window size for smoothing
        cycles (int) - number of repeating cycles
    """
    
    # approximate number of points within window
    window = int(window*len(signal)/(signal[-1][0]-signal[0][0]))
    window = min(window, len(signal))
    if window < 3:
        return signal.copy()
    if not window % 2:
        window -= 1
    
    # unpack mz and intensity
    xAxis, yAxis = np.hsplit(signal,2)
    xAxis = xAxis.flatten()
    yAxis = yAxis.flatten()
    
    # smooth the points
    while cycles:
        
        #CHECK_FORCE_QUIT()
        
        if style == 'flat':
            w = np.ones(window,'f')
        elif style == 'gaussian':
            r = np.array([(i-(window-1)/2.) for i in range(window)])
            w = np.exp(-(r**2/(window/4.)**2))
        else:
            w = eval('numpy.'+style+'(window)')        
        s = np.r_[yAxis[window-1:0:-1], yAxis, yAxis[-2:-window-1:-1]]
        y = np.convolve(w/w.sum(), s, mode='same')
        yAxis = y[window-1:-window+1]
        cycles -=1
    
    # return smoothed data
    xAxis.shape = (-1,1)
    yAxis.shape = (-1,1)
    data = np.concatenate((xAxis,yAxis), axis=1)
    
    return data.copy()

def getFragments(sequence, charge,series=['y','b'],fragment_losses=True,agentFormula='H'):
    if not agentFormula=='H': charge=1
    ls_mz=list()
    ls_fragmentSerie=list()
    ls_charge=list()
    ls_index=list()
    ls_fragmentLosses=list()
    charge=range(1,charge+1)
    fragments=fragment(sequence,series)
    if fragment_losses:
        fragments.extend(fragmentlosses(fragments, defined=True))
    for n in charge:
        for frag in fragments:
            ls_charge.append(n)
            ls_mz.append(frag.mz(n,agentFormula)[0])
            ls_fragmentSerie.append(frag.fragmentSerie)
            ls_index.append(frag.fragmentIndex)
            if frag.fragmentLosses:
                ls_fragmentLosses.append('-'+str(frag.fragmentLosses[0]))
            else: ls_fragmentLosses.append('')
                
                

    ls_ionNames=[frag+'%i%s%s'%(ls_index[n],ls_fragmentLosses[n],ls_charge[n]*'+') for n,frag in enumerate(ls_fragmentSerie)]
    buff=pd.DataFrame({'fragmentSerie':ls_fragmentSerie,'index':ls_index,'charge':ls_charge,'mz':ls_mz,'ionName':ls_ionNames,'fragmentLosses':ls_fragmentLosses})
    return(buff)




def getXICFig(df,title="",legend=False,smooth=False,window=0.05,cycles=2):
    fig,ax=plt.subplots(1,figsize=(9,3))
    for n,cond in enumerate(df['condition']):
        x_data=list(df['XIC_RT'][n])
        if x_data==[]:
            continue
        if not smooth:y_data=list(df['XIC_intens'][n])
        else:
            signal=np.dstack((x_data,list(df['XIC_intens'][n])))[0]
            signal=movaver(signal,window,cycles)
            y_data=np.hsplit(signal,2)[1]
        if not cond:
            label=os.path.basename(df.index[n])
        else:
            label=cond
        ax.plot(x_data,y_data,linewidth=1.0,color=df['color'][n],label=label)
    ax.legend(loc='best')
    ax.set_ylabel("intensity")
    #title=input_table['Peptide Sequence'][0]+', '+input_table['Site'][0]
    if title:
         ax.set_title(title,size=10)
    ax.set_xlabel('time')

    yfmt = ScalarFormatterForceFormat()
    yfmt.set_powerlimits((0,0))
    ax.yaxis.set_major_formatter(yfmt)   
    plt.tight_layout(pad=0.5, w_pad=2, h_pad=2)
    plt.close()
    return(fig)

def assignFragments(MS2spectra,fragments,error_tolerance=0.1):
    ls_ionMatch=list()
    ls_ionMatchIntens=list()
    for ion in fragments['mz']:
        minMZ=ion-error_tolerance
        maxMZ=ion+error_tolerance
        minpos=bisect.bisect(MS2spectra['raw_mz'][0],minMZ)
        maxpos=bisect.bisect(MS2spectra['raw_mz'][0],maxMZ)
        candidates=MS2spectra['raw_mz'][0][minpos:maxpos]
        if len(candidates)==0:
            ls_ionMatch.append(-1)
            ls_ionMatchIntens.append(-1)
            continue
        maxIonIntens=np.max(MS2spectra['raw_int'][0][minpos:maxpos])
        ionMatch=float(candidates[MS2spectra['raw_int'][0][minpos:maxpos]==np.max(MS2spectra['raw_int'][0][minpos:maxpos])])
        ls_ionMatch.append(ionMatch)
        ls_ionMatchIntens.append(maxIonIntens)
    fragments['ionMatches']=ls_ionMatch
    fragments['ionMatchIntens']=ls_ionMatchIntens
    return(fragments)


def modifyPeptide(sequence,modification,site,begin):
    sites=re.findall("\d+",site)
    sites=[int(item)-begin for item in sites]
    buff=list()
    modification=[x.strip() for x in modification.split(',')]
    for mod in modification:
        if mod=="nonspecific":continue            
        if mod=="Isomerization":continue
        if mod=="Unglycosylated":continue
        if mod=="K+":continue
        if mod=="Na+":continue
        if mod=="Gln->Pyro-Glu": mod="Gln2Pyro-Glu"
        buff.append(mod)
    for n,mod in enumerate(buff):
        try:
            sequence.setUnknownMod(sites[n],float(mod))
            return(sequence)
        except:pass
        try:
            output=sequence.modify(mod,sites[n],state='f')
        except: 
            return(False)
        if not output: 
            return(False)
    return(sequence)

def getMS2subplot(ax,MS2spectrum,fragments,flipped=False,relThreshould=0.1,label="",error_tolerance=0.1):
    fragments=assignFragments(MS2spectrum,fragments,error_tolerance)
    maxMZ=max(MS2spectrum['raw_mz'][0])
    minMZ=min(MS2spectrum['raw_mz'][0])
    maxIntens=max(MS2spectrum['raw_int'][0])
    width=(maxMZ-minMZ)/300
    x_data=MS2spectrum['raw_mz'][0]
    if flipped:y_data=MS2spectrum['raw_int'][0]*-1
    else: y_data=MS2spectrum['raw_int'][0]    
    ax.bar(x_data,y_data,width,color='grey')

    for series in set(fragments['fragmentSerie']):
        color=color_dict[series]
        buff=fragments[fragments['fragmentSerie']==series]
        buff=buff[buff['ionMatchIntens']>maxIntens*relThreshould]
        x_data=np.array(buff['ionMatches'])
        if flipped:y_data=np.array(buff['ionMatchIntens'])*-1
        else: y_data=np.array(buff['ionMatchIntens'])
        ax.bar(x_data,y_data,width,color=color)
        
        for n,ion in enumerate(buff['ionName']):
            if flipped: ax.annotate(ion,(x_data[n],y_data[n]),horizontalalignment='center',size=10,verticalalignment='top',color=color)
            else: ax.annotate(ion,(x_data[n],y_data[n]),horizontalalignment='center',size=10,verticalalignment='bottom',color=color)
        
        yfmt = ScalarFormatterForceFormat()
        yfmt.set_powerlimits((0,0))
        ax.yaxis.set_major_formatter(yfmt)  
        if flipped:
            y_max=-1*max(MS2spectrum['raw_int'][0])-0.1*max(MS2spectrum['raw_int'][0])
            ax.set_ylim(y_max,0)
            if ax.get_yticks()[-1]==0:
                ax.set_yticks(ax.get_yticks()[:-1])
            at=AnchoredText(label,loc=4,frameon=False)
            ax.add_artist(at)
        else:
            y_max=0.1*max(MS2spectrum['raw_int'][0])+max(MS2spectrum['raw_int'][0])
            ax.set_ylim(0,y_max)
            if ax.get_yticks()[0]==0:
                ax.set_yticks(ax.get_yticks()[1:])
            at=AnchoredText(label,loc=1,frameon=False)
            ax.add_artist(at)

    return(ax)

def getMS2fig(MS2spectra,fragments,title="",label=False,refrence="",error_tolerance=0.1):
    
    if refrence:
        fig=plt.figure(figsize=(9,6))
        ax1=plt.subplot(211)
        if title:
            ax1.set_title(title)
        ref=MS2spectra[MS2spectra['Location']==refrence]
        mz_values1=ref['raw_mz'][0]
        if len(mz_values1)==0:
            if label:
                ax_label=str(ref['condition'][0])
                ax1=getEmptyFig(ax1,flipped=False,chopped=True,label=ax_label)
            else:
                ax1=getEmptyFig(ax1,flipped=False,chopped=True,label=False)

        else:
            if label:
                ax_label=str(ref['condition'][0])
                ax1=getMS2subplot(ax1,ref,fragments,flipped=False,label=ax_label)
            else:    
                ax1=getMS2subplot(ax1,ref,fragments,flipped=False)
    
        ax2=plt.subplot(212,sharex=ax1)
        plt.setp(ax1.get_xticklabels(),visible=False)
        df=MS2spectra[MS2spectra['Location']!=refrence]
        mz_values2=df['raw_mz'][0]
        if  len(mz_values2)==0:
            if label:
                ax_label=str(df['condition'][0])
                ax2=getEmptyFig(ax2,flipped=True,chopped=True,label=ax_label)
            else:    
                ax2=getEmptyFig(ax2,flipped=True,chopped=True,label="")

        else:
            if label:
                ax_label=str(df['condition'][0])
                ax2=getMS2subplot(ax2,df,fragments,flipped=True,label=ax_label)
            else:
                ax2=getMS2subplot(ax2,df,fragments,flipped=True)
    
    else:
        fig=plt.figure(figsize=(9,3))
        ax1=plt.subplot(111)
        if title:
            ax1.set_title(title)
        df=MS2spectra
        mz_values1=df['raw_mz'].iloc[0]
        if len(mz_values1)==0:
            if label:
                ax_label=str(df['condition'][0])
                ax1=getEmptyFig(ax1,flipped=False,chopped=True,label=ax_label)
            else:
                ax1=getEmptyFig(ax1,flipped=False,chopped=True,label=False)

        else:
            if label:
                ax_label=str(df['condition'][0])
                ax1=getMS2subplot(ax1,df,fragments,flipped=False,label=ax_label)
            else:    
                ax1=getMS2subplot(ax1,df,fragments,flipped=False)

    fig.text(0.5,0.02,'m/z',ha='center')
    fig.text(0.02,0.5,'intensity',va='center',rotation='vertical')

    plt.tight_layout(pad=2, w_pad=0, h_pad=0.5)
    plt.close()
    return(fig)


def getEmptyFig(ax,flipped=False,chopped=False,label=""):
    if flipped:
        ax.set_ylim(-1,0)
        if chopped:
            ax.set_yticks([-0.2,-0.4,-0.6,-0.8,-1])
        #ax.annotate("NO SPECTRA",(0.5,-0.5),horizontalalignment='center',size=12,verticalalignment='center')
    else:
        ax.set_ylim(0,1)
        if chopped:
            ax.set_yticks(ax.get_yticks()[1:])
        #ax.annotate("NO SPECTRA",(0.5,0.5),horizontalalignment='center',size=12,verticalalignment='center')
    
    return(ax)
    
def createEmptyXML():
    XML="""<?xml version='1.0' encoding='UTF-8'?>
    <configuration></configuration>
    """
    return(XML)


def createConfigFile(param,path):
    emptyXML=createEmptyXML()
    dom= xml.parseString(emptyXML)
    root=dom.childNodes[0]
    for item in param._items:
        if item._name[0]=="_" or item._name[:4]=="file"  or item._name[-6:]=="button": continue
        newScript=dom.createElement("param")
        newScript.setAttribute("name",item._name)
        newScript.setAttribute("value",str(getattr(param, "_"+item._name)))
        root.appendChild(newScript)    

    file_handle = open(path,"w")
    dom.writexml(file_handle,addindent='\t',newl='\n')
    file_handle.close()

def loadConfiguration(fn):
    dict_bool={'False':False,'True':True}
    dom=xml.parse(fn)
    root=dom.childNodes[0]
    parameter=root.getElementsByTagName("param")
    dict_parameter={}
    for item in parameter:
        param_name=item.getAttribute("name")
        param_value=item.getAttribute("value")
        if param_name in ['XIC_fixedrange','XIC_fullrange','smooth','XIC_legend','MS_fullrange','MS_title','MS_legend','XIC_title','XIC_threshold','XIC_isotopicPattern','MS2_title','MS2_legend']:
            param_value=dict_bool[item.getAttribute("value")]
        elif re.match("\d+\.\d+",param_value):
            param_value=float(item.getAttribute("value"))
        elif re.match("\d+",param_value):
            param_value=int(item.getAttribute("value"))
        elif param_value[0]=='(' or param_value[0]=='[':
            param_value=re.findall("\d",param_value)
            param_value=[int(item) for item in param_value]
            
        dict_parameter[param_name]=param_value
        
    return(dict_parameter)
      
def createDefaultXML(path):
    import os
    fn="default_configuration.xml"
    fn=os.path.join(os.path.dirname(path),fn)
    fn=fn.replace('\\','/')
    file_handle = open(fn,"w")
    
    XML="""<?xml version="1.0" ?>
<configuration>
	<param name="XIC_fullrange" value="True"/>
	<param name="XIC_min_range" value="1.0"/>
	<param name="XIC_max_range" value="1.0"/>
	<param name="smooth" value="True"/>
	<param name="cycles" value="1"/>
	<param name="window" value="0.1"/>
	<param name="XIC_title" value="True"/>
	<param name="XIC_legend" value="True"/>
	<param name="XIC_linewidth" value="1.0"/>
	<param name="MS_fullrange" value="False"/>
	<param name="MS_min_range" value="0.5"/>
	<param name="MS_max_range" value="0.5"/>
	<param name="MS_title" value="True"/>
	<param name="MS_legend" value="True"/>
	<param name="MS_linewidth" value="1.0"/>
	<param name="ion_series" value="[1, 4, 6]"/>
	<param name="error_tolerance" value="0.5"/>
	<param name="error_unit" value="0"/>
	<param name="assign_to" value="0"/>
</configuration>
    """
    
    file_handle.write(XML)
    file_handle.close()
    return(fn)
  

def getUV(path):
    GetModule(r"C:\Program Files (x86)\Thermo\Foundation\XRawfile2.dll")
    obj = CreateObject('XRawfile.XRawfile')
    pnArraySize = c_long()
    scanFilter=''
    chroType1=0
    chroOperator=0
    chroType2=0
    delay=0.0
    smoothingType=0
    smoothingValue=0
    massRange2=""
    startTime=0
    endTime=0
    massrange=""
    obj.Open(path)
    obj.SetCurrentController(c_long(4), c_long(1))
    pvarChroData = comtypes.automation.VARIANT()
    pvarPeakFlags = comtypes.automation.VARIANT()
    
    obj.GetChroData(chroType1,
                      chroOperator,
                      chroType2,
                      scanFilter,
                      massrange,
                      massRange2,
                      c_double(delay),
                      byref(c_double(startTime)),
                      byref(c_double(endTime)),
                      smoothingType,
                      smoothingValue,
                      pvarChroData,
                      pvarPeakFlags,
                      byref(pnArraySize))
    UV_RT=[float(m) for m in pvarChroData[0][0]]
    UV_intens=[float(m) for m in pvarChroData[0][1]]
    obj.Close()
    
    return(UV_RT,UV_intens)
    

def getXIC(path,massrange,startTime=0,endTime=0,threshould=0,fullrange=False):
    GetModule(r"C:\Program Files (x86)\Thermo\Foundation\XRawfile2.dll")
    obj = CreateObject('XRawfile.XRawfile')
    pnArraySize = c_long()
    scanFilter="Full ms"
    chroType1=0
    chroOperator=0
    chroType2=0
    delay=0.0
    smoothingType=0
    smoothingValue=0
    massRange2=""
    massrange=massrange
    obj.Open(path)
    obj.SetCurrentController(c_long(0), c_long(1))
    pvarChroData = comtypes.automation.VARIANT()
    pvarPeakFlags = comtypes.automation.VARIANT()
    
    minStartTime=c_double()
    maxStartTime=c_double()
    
    
    
    if startTime:
        obj.GetStartTime(byref(minStartTime))
        if startTime<minStartTime.value:
            startTime=minStartTime.value
        if startTime>endTime:
            raise IOError("startTime must be lower than endTime!")
            
    
    if endTime:
        obj.GetEndTime(byref(maxStartTime))
        if endTime>maxStartTime.value:
            endTime=maxStartTime.value
    

    obj.GetChroData(chroType1,
                      chroOperator,
                      chroType2,
                      scanFilter,
                      massrange,
                      massRange2,
                      c_double(delay),
                      byref(c_double(startTime)),
                      byref(c_double(endTime)),
                      smoothingType,
                      smoothingValue,
                      pvarChroData,
                      pvarPeakFlags,
                      byref(pnArraySize))
    #XIC_RT=[float(m) for m in pvarChroData[0][0]]
    #XIC_intens=[float(m) for m in pvarChroData[0][1]]
    obj.Close()
    signal=np.dstack((pvarChroData[0][0],pvarChroData[0][1]))[0]
    signal=movaver(signal,0.1,2)
    return(signal[:,0],signal[:,1])          



def getMSspectra(path,RT,scanFilter="Full ms"):
    GetModule(r"C:\Program Files (x86)\Thermo\Foundation\XRawfile2.dll")
    obj = CreateObject('XRawfile.XRawfile')
    obj.Open(path)
    obj.SetCurrentController(c_long(0), c_long(1))
    RT=c_double(RT)
    numOfSpectra=c_long()
    pnScanNumber=c_long()
    obj.ScanNumFromRT(RT,pnScanNumber)
    pnMassOrder=c_long()
    obj.GetMSOrderForScanNum(pnScanNumber,pnMassOrder)
    obj.GetNumSpectra(numOfSpectra)
    while pnMassOrder.value>1 and pnScanNumber.value<=numOfSpectra.value:
        pnScanNumber=c_long(pnScanNumber.value+1)
        obj.GetMSOrderForScanNum(pnScanNumber,pnMassOrder)
    
    scanFilter=scanFilter
    intensityCufoffType=0
    intensityCufoffValue=0
    nMaxNumberOfPeaks=0
    CentroidResult=False
    centroidPeakwidth=c_double(0.0)
    pvarMassList=comtypes.automation.VARIANT()
    pvarPeakFlags = comtypes.automation.VARIANT()
    pnArraySize = c_long()
    
    obj.GetMassListFromScanNum(pnScanNumber,
                      scanFilter,
                      intensityCufoffType,
                      intensityCufoffValue,
                      nMaxNumberOfPeaks,
                      CentroidResult,
                      centroidPeakwidth,                      
                      pvarMassList,
                      pvarPeakFlags,
                      byref(pnArraySize))
    obj.Close()
    if not pvarMassList.value:
        return(0,0,0)
        
    return(np.array(pvarMassList[0][0]),np.array(pvarMassList[0][1]),pnScanNumber.value) 


def getMSLabel(path,RT=0,scanNumber=0):
    if not (RT or scanNumber): 
        raise IOError("RT or scannumber must be defined!")
     
    numOfSpectra=c_long()
    GetModule(r"C:\Program Files (x86)\Thermo\Foundation\XRawfile2.dll")
    obj = CreateObject('XRawfile.XRawfile')
    obj.Open(path)
    obj.SetCurrentController(c_long(0), c_long(1))
    obj.GetNumSpectra(numOfSpectra)
    if not scanNumber:
        retentionTime=c_double(RT)
        pnScanNumber=c_long()
        minStartTime=c_double()
        maxStartTime=c_double()    
        obj.GetStartTime(byref(minStartTime))
        obj.GetEndTime(byref(maxStartTime))
        if RT and minStartTime.value>RT:
            raise IOError("Retention time is too low!")
        if RT and maxStartTime.value<RT:
            raise IOError("Retention time is too high!")
        obj.ScanNumFromRT(retentionTime,pnScanNumber)
        pnMassOrder=c_long()
        obj.GetMSOrderForScanNum(pnScanNumber,pnMassOrder)
        while pnMassOrder.value>1:
            pnScanNumber=c_long(pnScanNumber.value+1)
            obj.GetMSOrderForScanNum(pnScanNumber,pnMassOrder)
    elif scanNumber> numOfSpectra.value:
        raise IOError("ScanNumber is too high!")           
    else: 
        pnScanNumber=c_long(scanNumber)

            
    pvarLabels=comtypes.automation.VARIANT()
    pvarFlags = comtypes.automation.VARIANT()

    
    obj.GetLabelData(pvarLabels,
                      pvarFlags,
                      pnScanNumber)

    
    obj.Close()
    return(np.array(pvarLabels[0][0]),np.array(pvarLabels[0][1]),np.array(pvarLabels[0][5]))

def getMS2Label(path,scanNumber):
    if scanNumber==-1:
        return([],[])
    numOfSpectra=c_long()
    GetModule(r"C:\Program Files (x86)\Thermo\Foundation\XRawfile2.dll")
    obj = CreateObject('XRawfile.XRawfile')
    obj.Open(path)
    obj.SetCurrentController(c_long(0), c_long(1))
    obj.GetNumSpectra(numOfSpectra)
    if scanNumber> numOfSpectra.value:
        raise IOError("ScanNumber is too high!")           
    pnScanNumber=c_long(scanNumber)            
    pvarMassList=comtypes.automation.VARIANT()
    pvarFlags = comtypes.automation.VARIANT()
    szFilter=""
    intensityCutoffType=0
    intensityCutoffValue=0
    nMaxNumberOfPeaks=0
    bCetroidResult=0
    centroidPeakWidth=0
    pnArraySize = c_long()
    obj.GetMassListFromScanNum(pnScanNumber,
                               szFilter, 
                               intensityCutoffType,
                               intensityCutoffValue,
                               nMaxNumberOfPeaks,
                               bCetroidResult,
                               c_double(centroidPeakWidth),
                               pvarMassList, 
                               pvarFlags,
                               byref(pnArraySize))

    
    obj.Close()
    
    return(np.array(pvarMassList[0][0]),np.array(pvarMassList[0][1]))

def defineIsotopeMZCluster(mz,charge):
    ISOTOPE_DISTANCE=1.00287
    difference=ISOTOPE_DISTANCE/charge
    cluster=[mz]
    for i in range(1,9):
        cluster.append(mz+i*difference)
    return(cluster,difference)         



def searchIsotopeCluster(MSLabel,mz,intens=0,charge=0,monoisotopicMass=0,mz_error=0.1,int_error=0.9):
    
    mz_values,int_values,charge_values=MSLabel
    lower_index=bisect.bisect_left(mz_values,mz-mz*20e-6)
    upper_index=bisect.bisect_right(mz_values,mz+mz*20e-6)
    if lower_index==upper_index:
            print "no mass match"
            return([])

    buff_mz=np.array(mz_values[lower_index:upper_index])
    buff_mz_error=abs(buff_mz-mz)
    mz=buff_mz[buff_mz_error.argmin()]
    
    if not intens:
        intens=int_values[mz_values==mz][0]
    if not charge:
        charge=charge_values[mz_values==mz][0]
    if not monoisotopicMass:
        monoisotopicMass=mz*charge-1.0078250321*charge
    isotopeMZcandidates,difference=defineIsotopeMZCluster(mz,charge)

    
    cluster=np.array([[mz,intens]])
    int_monoisotopic=intens
    
    mass = int(monoisotopicMass/200)
    pattern = patternLookupTable[mass]
    limit = min(len(pattern), len(isotopeMZcandidates))

    for n,cand in enumerate(isotopeMZcandidates[1:limit]):
        isotope=n+1
        lower_index = bisect.bisect_left(mz_values, (float(cand - mz_error)))
        upper_index = bisect.bisect_right(mz_values, (float(cand + mz_error)))        
        if lower_index==upper_index:
            break
        cand_int=int_values[lower_index:upper_index]
        cand_mz=mz_values[lower_index:upper_index]
        cand_charge=charge_values[lower_index:upper_index]

        valid_candidate=[]
        for m,mz in enumerate(cand_mz):            
            if cand_charge[m]!= charge:continue        
            intTheoretical = (int_monoisotopic / pattern[0]) * pattern[isotope]
            intError = cand_int[m] - intTheoretical
            if abs(intError) <= (intTheoretical * int_error):
                valid_candidate.append([mz,cand_int[m],int_error]) 
        if valid_candidate:
            valid_candidate=np.array(valid_candidate)
            valid_candidate=valid_candidate[valid_candidate[:,2].argmin()][:2]           
            cluster=np.append(cluster, [valid_candidate], axis=0)
        else:
            break

    return(cluster)
                

    
def getMSspectraLimits(cluster):
    cluster_xmax=max(cluster[:,0])
    cluster_xmin=min(cluster[:,0])
    cluster_xdiffrence=cluster_xmax-cluster_xmin
    if len(cluster)==1:
        return(cluster_xmin-0.3,cluster_xmin+3)        
    elif len(cluster)<=3:
        min_diff=min(np.diff(cluster[:,0]))
        return(cluster_xmin-0.3,cluster_xmin+3.25*min_diff)
    return(cluster_xmin-0.1*cluster_xdiffrence,cluster_xmax+0.1*cluster_xdiffrence)
    
def getMSFigAxes(ax,x_data,y_data,MS1Labels,flipped=False,title="",label="",color='black'):    
    if MS1Labels==[]:
        ax=getEmptyFig(ax,flipped,True,label)
        return(ax,0,0,False)
        
    if not isinstance(MS1Labels, np.ndarray):
        MS1Labels=np.array(MS1Labels)
        
    if flipped: 
        y_data=y_data*-1
    ax.plot(x_data,y_data,linewidth=1.0,color=color,label=label)
    if title: ax.set_title(title,size=10)
    cluster=MS1Labels
    cluster_ymax=max(cluster[:,1])
    if flipped:
        y_max=cluster_ymax*(-1)-0.1*cluster_ymax
        ax.set_ylim(y_max,0)
    else:
        y_max=0.1*cluster_ymax+cluster_ymax
        ax.set_ylim(0,y_max)
    for n,peak in enumerate(cluster):
        text='{:.4f}'.format(round(peak[0],4))
        if flipped:
            ax.annotate(text,(peak[0],peak[1]*-1-5),horizontalalignment='right',size=9,verticalalignment='top')
        else:
            ax.annotate(text,(peak[0],peak[1]+5),horizontalalignment='right',size=9,verticalalignment='bottom')
    yfmt = ScalarFormatterForceFormat()
    yfmt.set_powerlimits((0,0))
    ax.yaxis.set_major_formatter(yfmt)         
    cluster_xmin,cluster_xmax=getMSspectraLimits(cluster)
    return(ax,cluster_xmin,cluster_xmax,True)

def getMSFig(df_spectra,RT,mz,charge,monoisotopicMass,title="",label=True, refrence=""):
    if refrence:
        ref=df_spectra[df_spectra['Location']==refrence]
        data=df_spectra[df_spectra['Location']!=refrence]
        fig=plt.figure(figsize=(9,6))
        ax1=plt.subplot(211)
        ax2=plt.subplot(212,sharex=ax1)
        try:ref_label=ref['condition'][0]
        except:
            raise IOError("There is no refrence defined!")
        data_label=data['condition'][0]
        ref_color=ref['color'][0]
        data_color=data['color'][0]
        ref_loc=str(ref['Location'][0])
        data_loc=str(data['Location'][0])
        ref_mz_values,ref_intens_values,ref_scanNumber=getMSspectra(ref_loc,RT)
        data_mz_values,data_intens_values,data_scanNumber=getMSspectra(data_loc,RT)
        if ref_scanNumber==0 or data_scanNumber==0: raise IOError('Retention time might be wrong')
        ref_MSLabel=getMSLabel(ref_loc,RT=0,scanNumber=ref_scanNumber)
        data_MSLabel=getMSLabel(data_loc,RT=0,scanNumber=data_scanNumber)
        ref_cluster=searchIsotopeCluster(ref_MSLabel,mz,charge=charge,monoisotopicMass=monoisotopicMass)
        data_cluster=searchIsotopeCluster(data_MSLabel,mz,charge=charge,monoisotopicMass=monoisotopicMass)
        ax1,ref_min,ref_max,bool_ax1=getMSFigAxes(ax1,ref_mz_values,ref_intens_values,ref_cluster,label=ref_label,color=ref_color)
        ax2,data_min,data_max,bool_ax2=getMSFigAxes(ax2,data_mz_values,data_intens_values,data_cluster,label=data_label,color=data_color,flipped=True)
        
        if bool_ax1 and bool_ax2:
            ax_min=min(ref_min,data_min)
            ax_max=max(ref_max,data_max)
        elif bool_ax1:
            ax_min=ref_min
            ax_max=ref_max
        elif bool_ax2:
            ax_min=data_min
            ax_max=data_max
        else:
            ax_min=0
            ax_max=1
        ax1.set_xlim(ax_min,ax_max)
        ax2.set_xlim(ax_min,ax_max)
        
        if not bool_ax1:
            ax1.annotate("NO SPECTRA",((ax_max-ax_min)/2,0.5),horizontalalignment='center',size=12,verticalalignment='center')
        if not bool_ax2:
            ax2.annotate("NO SPECTRA",((ax_max-ax_min)/2,-0.5),horizontalalignment='center',size=12,verticalalignment='center')
        ax2.set_xlabel('mz')
        ax1.set_yticks(ax1.get_yticks()[1:])
        plt.setp(ax1.get_xticklabels(),visible=False)
        ax2.set_yticks(ax2.get_yticks()[:-1])
        if title:
            ax1.set_title(title,size=10)
        if label:
            ax1.legend(loc='upper right')
            ax2.legend(loc='lower right')                                  
    else:
        fig=plt.figure(figsize=(9,3))
        ax1=plt.subplot(111)
    
        label=df_spectra['condition'][0]
        color=df_spectra['color'][0]
        loc=df_spectra['Location'][0]
        mz_values,intens_values,scanNumber=getMSspectra(loc,RT)
        if scanNumber==0: raise IOError('Retention time might be wrong')
        MSLabel=getMSLabel(loc,RT=0,scanNumber=scanNumber)
        cluster=searchIsotopeCluster(MSLabel,mz,charge=charge,monoisotopicMass=monoisotopicMass)
        print cluster
        ax1,buff_min,buff_max,ax1_bool=getMSFigAxes(ax1,mz_values,intens_values,cluster,label=label,color=color)
        ax1.set_xlim(buff_min,buff_max)
        if not ax1_bool:
            ax1.annotate("NO SPECTRA",((buff_max-buff_min)/2,0.5),horizontalalignment='center',size=12,verticalalignment='center')
        ax1.set_xlabel('mz')
        if title:
            ax1.set_title(title,size=10)
        if label:
            ax1.legend(loc='upper right')
    
    fig.text(0.0,0.5,'intensity',va='center',rotation='vertical')
    plt.tight_layout(pad=1, w_pad=0, h_pad=0.5)
    plt.close()
    return(fig)
    
def getXICDataFrame(data_frame,mz,charge,mono_mass,startTime,endTime,mz_error,isotope_pattern,RT,threshould=0):
    rows=list()
    for n,path in enumerate(data_frame['Location']):
        if isotope_pattern:
            MSLabel=getMSLabel(str(path),RT)
            cluster=searchIsotopeCluster(MSLabel,mz,charge=charge,monoisotopicMass=mono_mass)
            if cluster==[]:
                massrange="%f-%f"%(mz-mz_error*mz,mz+mz_error*mz)
            elif not isinstance(cluster, np.ndarray):
                cluster=np.ndarray(cluster)
                massrange=""
                for mass in cluster[:,0]:
                    massrange+="%f-%f,"%(mass-mz_error*mass,mass+mz_error*mass)
                massrange=massrange[:-1]
            else:
                massrange=""
                for mass in cluster[:,0]:
                    massrange+="%f-%f,"%(mass-mz_error*mass,mass+mz_error*mass)
                massrange=massrange[:-1]
                
            
        else:
            massrange="%f-%f"%(mz-mz_error*mz,mz+mz_error*mz)
        XIC_RT,XIC_intens=getXIC(path,massrange,startTime,endTime,threshould)       
        area=np.trapz(XIC_intens,XIC_RT)                  
        rows.append([path,XIC_RT,XIC_intens,area])
    return (pd.DataFrame(rows, columns=['Location', 'XIC_RT','XIC_intens','area']))     

def getXICDataFrame_1(locations,mz,refrence="",mz_error=10e-5,RT=0,isotope_pattern=False,
                      charge=0,mono_mass=0,minRT=0,maxRT=0,threshold=0,
                      highest_intense=False):
    from mMass.mspy import mod_peakpicking as peak_picking
    rows=list()
    massrange="%f-%f"%(mz-mz_error*mz,mz+mz_error*mz)
    if refrence:
        XIC_RT,XIC_intens=getXIC(refrence,massrange)
        signal=np.dstack((XIC_RT,XIC_intens))[0]
        peaklist=peak_picking.labelscan(signal,relThreshold=threshold)
        ls_peakRT=np.array([])
        ls_peakInt=np.array([])
        for peak in peaklist:
            ls_peakRT=np.append(ls_peakRT,peak.mz)
            ls_peakInt=np.append(ls_peakInt,peak.ai)
        if RT==0:
            realRT=ls_peakRT[ls_peakInt.argmax()]
        else:
            realRT=ls_peakRT[(abs(ls_peakRT-RT)).argmin()]
        MSLabel=getMSLabel(str(refrence),realRT)
        ref_cluster=searchIsotopeCluster(MSLabel,mz,charge=charge,monoisotopicMass=mono_mass)
        
        if isotope_pattern:
            if ref_cluster==[]:
                pass
            elif not isinstance(ref_cluster, np.ndarray):
                ref_cluster=np.ndarray(ref_cluster)
                massrange=""
                for mass in ref_cluster[:,0]:
                    massrange+="%f-%f,"%(mass-mz_error*mass,mass+mz_error*mass)
                massrange=massrange[:-1]

            else:
                massrange=""
                for mass in ref_cluster[:,0]:
                    massrange+="%f-%f,"%(mass-mz_error*mass,mass+mz_error*mass)
                massrange=massrange[:-1]

        
        elif highest_intense:
            if ref_cluster==[]:
                pass
            elif not isinstance(ref_cluster, np.ndarray):
                ref_cluster=np.ndarray(ref_cluster)
                ref_cluster=ref_cluster[ref_cluster[:,1].argmax()]
                massrange=""
                mass=ref_cluster[0]
                massrange+="%f-%f"%(mass-mz_error*mass,mass+mz_error*mass)

            else:
                massrange=""
                ref_cluster=ref_cluster[ref_cluster[:,1].argmax()]
                mass=ref_cluster[0]
                massrange+="%f-%f"%(mass-mz_error*mass,mass+mz_error*mass)

            
        
    for n,path in enumerate(locations):        
        XIC_RT,XIC_intens=getXIC(path,massrange)
        signal=np.dstack((XIC_RT,XIC_intens))[0]
        peaklist=peak_picking.labelscan(signal,relThreshold=threshold)
        ls_peakRT=np.array([])
        ls_peakInt=np.array([])
        for peak in peaklist:
            ls_peakRT=np.append(ls_peakRT,peak.mz)
            ls_peakInt=np.append(ls_peakInt,peak.ai)
        if RT==0:
            realRT=ls_peakRT[ls_peakInt.argmax()]
        else:
            realRT=ls_peakRT[(abs(ls_peakRT-RT)).argmin()]
        
        if minRT==0:
            startTime=0
            endTime=0
            #XIC_RT,XIC_intens=getXIC(path,massrange,startTime,endTime)
        else:
            startTime=realRT-minRT
            endTime=realRT+minRT
            XIC_RT,XIC_intens=getXIC(path,massrange,startTime,endTime)
        area=np.trapz(XIC_intens,XIC_RT)                  
        rows.append([path,XIC_RT,XIC_intens,area,realRT])
    return (pd.DataFrame(rows, columns=['Location', 'XIC_RT','XIC_intens','area','RT']))

def getFilters(path):
    header = OrderedDict()
    header['numPackets'] = c_long()
    header['StartTime'] = c_double()
    header['LowMass'] = c_double()
    header['HighMass'] = c_double()
    header['TIC'] = c_double()
    header['BasePeakMass'] = c_double()
    header['BasePeakIntensity'] = c_double()
    header['numChannels'] = c_long()
    header['uniformTime'] = c_long()
    header['Frequency'] = c_double()
    GetModule(r"C:\Program Files (x86)\Thermo\Foundation\XRawfile2.dll")
    obj = CreateObject('XRawfile.XRawfile')
    pnScanNumber=c_long()
    numOfSpectra=c_long()
    obj.Open(path)
    obj.SetCurrentController(c_long(0), c_long(1))    

    obj.GetNumSpectra(numOfSpectra)
    ls_MS2Scans=list()
    for scan in range(1,100):#numOfSpectra.value+1):
        pnMassOrder=c_long()
        pnScanNumber=c_long(scan)
        obj.GetMSOrderForScanNum(pnScanNumber,pnMassOrder)
        pnNumberOfMSOrders = c_long(1)
        if pnMassOrder.value>1:
            ls_MS2Scans.append(pnScanNumber.value)
            print(pnScanNumber.value)
            print(pnMassOrder.value)   
            pvarLabels = comtypes.automation.VARIANT()
            pvarFlags = comtypes.automation.VARIANT(4)
            obj.GetAllMSOrderData(pnScanNumber,pvarLabels,pvarFlags,byref(pnNumberOfMSOrders))
            print(pvarLabels)
            #break                            
    
    obj.Close() 
    return(pvarFlags)

def getPrecursorInformation(path):
    GetModule(r"C:\Program Files (x86)\Thermo\Foundation\XRawfile2.dll")
    obj = CreateObject('XRawfile.XRawfile')
    obj.Open(path)
    obj.SetCurrentController(c_long(0), c_long(1))
    pnMassOrder=c_long()
    numOfSpectra=c_long()
    rows=list()
    obj.GetNumSpectra(numOfSpectra)
    for scan in range(1,numOfSpectra.value+1):
        pnScanNumber=c_long(scan)
        obj.GetMSOrderForScanNum(pnScanNumber,pnMassOrder)
        if pnMassOrder.value==1:
            parentScanNumber=pnScanNumber
            parentMSLabel=comtypes.automation.VARIANT()
            pvarFlags = comtypes.automation.VARIANT()
            obj.GetLabelData(parentMSLabel,
                               pvarFlags,
                               parentScanNumber)
            parent_mzValues=np.array(parentMSLabel[0][0])
            parent_intensValues=np.array(parentMSLabel[0][1])
            parent_charge_values=np.array(parentMSLabel[0][5])
        if pnMassOrder.value>1:
            precursorMass = c_double()
            numPackets = c_long()
            StartTime = c_double()
            LowMass = c_double()
            HighMass = c_double()
            TIC = c_double()
            BasePeakMass = c_double()
            BasePeakIntensity = c_double()
            numChannels = c_long()
            uniformTime = c_long()
            Frequency = c_double()
            obj.GetPrecursorMassForScanNum(pnScanNumber, pnMassOrder,byref(precursorMass))
            obj.GetScanHeaderInfoForScanNum(pnScanNumber, numPackets, StartTime,
                                                        LowMass, HighMass,
                                                        TIC, BasePeakMass,
                                                        BasePeakIntensity, numChannels,
                                                        uniformTime, Frequency)
            precursor_mz=precursorMass.value
            lower_index=0
            upper_index=0
            n=1
            while lower_index==upper_index:              
                lower_index=bisect.bisect_left(parent_mzValues,precursor_mz-precursor_mz*5e-6*n)
                upper_index=bisect.bisect_right(parent_mzValues,precursor_mz+precursor_mz*5e-6*n)
                n+=1
            if upper_index-lower_index>1:
                buff_mz=np.array(parent_mzValues[lower_index:upper_index])
                buff_mz_error=abs(buff_mz-precursor_mz)
                buff_mz=buff_mz[buff_mz_error.argmin()]
                intens=parent_intensValues[parent_mzValues==buff_mz][0]
                charge=parent_charge_values[parent_mzValues==buff_mz][0]
            else: 
                buff_mz=parent_mzValues[lower_index:upper_index][0]
                intens=parent_intensValues[lower_index:upper_index][0]
                charge=parent_charge_values[lower_index:upper_index][0]
            rows.append([scan,precursor_mz,StartTime.value,TIC.value,intens,charge,buff_mz])
    df=pd.DataFrame(rows,columns=['scanNumber','PrecursorMass','StartTime','TIC','precursorIntens','precursorCharge','observed_mz'])
    obj.Close()
    return df
    
def getMS2ScanNumber(path,ls_mz,ls_RT,ls_charge=[],ppm=10e-6):
    df_MS2=getPrecursorInformation(path)
    ls_scanNumber=list()
    for n,mz in enumerate(ls_mz):
        RT=ls_RT[n]
        temp_df=df_MS2[(df_MS2['StartTime']>RT-3)&(df_MS2['StartTime']<RT+3)]
        if ls_charge:
            charge=ls_charge[n]
            result_df=temp_df[df_MS2['precursorCharge']==charge]
        result_df=result_df[(result_df['PrecursorMass']>(mz-mz*ppm))&(result_df['PrecursorMass']<(mz+mz*ppm))]
        if len(result_df)==0:
            result_df=temp_df[(result_df['observed_mz']>(mz-mz*ppm))&(temp_df['observed_mz']<(mz+mz*ppm))]
        if len(result_df)==0:
            ls_scanNumber.append(-1)
        elif len(result_df)==1:
            ls_scanNumber.append(int(result_df['scanNumber'].iloc[0]))
        else:
            result_df=result_df[result_df['precursorIntens']==max(result_df['precursorIntens'])]
            ls_scanNumber.append(int(result_df['scanNumber'].iloc[0]))
    return(ls_scanNumber)
            
            
            

patternLookupTable = (
    (1.000, 0.059, 0.003), #0
    (1.000, 0.122, 0.013), #200
    (1.000, 0.241, 0.040, 0.005), #400
    (1.000, 0.303, 0.059, 0.008), #600
    (1.000, 0.426, 0.109, 0.020, 0.003), #800
    (1.000, 0.533, 0.166, 0.038, 0.006), #1000
    (1.000, 0.655, 0.244, 0.066, 0.014, 0.002), #1200
    (1.000, 0.786, 0.388, 0.143, 0.042, 0.009, 0.001), #1400
    (1.000, 0.845, 0.441, 0.171, 0.053, 0.013, 0.002), #1600
    (1.000, 0.967, 0.557, 0.236, 0.080, 0.021, 0.005), #1800
    (0.921, 1.000, 0.630, 0.291, 0.107, 0.032, 0.007, 0.001), #2000
    (0.828, 1.000, 0.687, 0.343, 0.136, 0.044, 0.011, 0.002), #2200
    (0.752, 1.000, 0.744, 0.400, 0.171, 0.060, 0.017, 0.004), #2400
    (0.720, 1.000, 0.772, 0.428, 0.188, 0.068, 0.020, 0.005), #2600
    (0.667, 1.000, 0.825, 0.487, 0.228, 0.088, 0.028, 0.007), #2800
    (0.616, 1.000, 0.884, 0.556, 0.276, 0.113, 0.039, 0.010, 0.002), #3000
    (0.574, 1.000, 0.941, 0.628, 0.330, 0.143, 0.052, 0.015, 0.003), #3200
    (0.536, 0.999, 1.000, 0.706, 0.392, 0.179, 0.069, 0.022, 0.005), #3400
    (0.506, 0.972, 1.000, 0.725, 0.412, 0.193, 0.077, 0.025, 0.006), #3600
    (0.449, 0.919, 1.000, 0.764, 0.457, 0.226, 0.094, 0.033, 0.009, 0.001), #3800
    (0.392, 0.853, 1.000, 0.831, 0.543, 0.295, 0.136, 0.053, 0.017, 0.004), #4000
    (0.353, 0.812, 1.000, 0.869, 0.593, 0.336, 0.162, 0.067, 0.023, 0.006), #4200
    (0.321, 0.776, 1.000, 0.907, 0.644, 0.379, 0.190, 0.082, 0.030, 0.009), #4400
    (0.308, 0.760, 1.000, 0.924, 0.669, 0.401, 0.205, 0.090, 0.033, 0.011, 0.001), #4600
    (0.282, 0.729, 1.000, 0.962, 0.723, 0.451, 0.239, 0.110, 0.042, 0.014, 0.003), #4800
    (0.258, 0.699, 1.000, 1.000, 0.780, 0.504, 0.277, 0.132, 0.053, 0.018, 0.004), #5000
    (0.228, 0.645, 0.962, 1.000, 0.809, 0.542, 0.308, 0.153, 0.065, 0.023, 0.007), #5200
    (0.203, 0.598, 0.927, 1.000, 0.839, 0.581, 0.343, 0.176, 0.078, 0.029, 0.010), #5400
    (0.192, 0.577, 0.911, 1.000, 0.854, 0.602, 0.361, 0.189, 0.086, 0.033, 0.011), #5600
    (0.171, 0.536, 0.880, 1.000, 0.884, 0.644, 0.399, 0.216, 0.102, 0.040, 0.014, 0.003), #5800
    (0.154, 0.501, 0.851, 1.000, 0.912, 0.686, 0.439, 0.244, 0.120, 0.050, 0.018, 0.004), #6000
    (0.139, 0.468, 0.823, 1.000, 0.942, 0.730, 0.482, 0.278, 0.141, 0.062, 0.023, 0.007), #6200
    (0.126, 0.441, 0.799, 1.000, 0.969, 0.772, 0.524, 0.310, 0.162, 0.073, 0.028, 0.009), #6400
    (0.121, 0.427, 0.787, 1.000, 0.983, 0.794, 0.547, 0.328, 0.174, 0.080, 0.031, 0.011), #6600
    (0.104, 0.381, 0.732, 0.971, 1.000, 0.848, 0.614, 0.390, 0.219, 0.109, 0.045, 0.016, 0.004), #6800
    (0.092, 0.349, 0.691, 0.944, 1.000, 0.872, 0.648, 0.422, 0.244, 0.125, 0.054, 0.020, 0.006), #7000
    (0.082, 0.321, 0.654, 0.919, 1.000, 0.894, 0.682, 0.456, 0.270, 0.143, 0.063, 0.024, 0.008), #7200
    (0.073, 0.296, 0.620, 0.895, 1.000, 0.917, 0.718, 0.492, 0.299, 0.162, 0.077, 0.030, 0.011), #7400
    (0.069, 0.284, 0.604, 0.884, 1.000, 0.929, 0.735, 0.509, 0.313, 0.172, 0.084, 0.033, 0.012), #7600
    (0.062, 0.262, 0.573, 0.861, 1.000, 0.952, 0.772, 0.548, 0.345, 0.195, 0.098, 0.040, 0.015, 0.003), #7800
    (0.056, 0.243, 0.544, 0.839, 1.000, 0.976, 0.811, 0.589, 0.380, 0.220, 0.114, 0.049, 0.019, 0.005), #8000
    (0.051, 0.227, 0.521, 0.821, 1.000, 0.997, 0.846, 0.628, 0.413, 0.244, 0.130, 0.058, 0.022, 0.007), #8200
    (0.045, 0.206, 0.486, 0.786, 0.980, 1.000, 0.869, 0.660, 0.444, 0.268, 0.147, 0.070, 0.027, 0.010), #8400
    (0.042, 0.196, 0.468, 0.767, 0.968, 1.000, 0.879, 0.676, 0.460, 0.281, 0.156, 0.075, 0.030, 0.011), #8600
    (0.038, 0.179, 0.437, 0.733, 0.947, 1.000, 0.899, 0.705, 0.491, 0.307, 0.173, 0.086, 0.036, 0.013, 0.002), #8800
    (0.033, 0.163, 0.408, 0.701, 0.926, 1.000, 0.919, 0.736, 0.524, 0.335, 0.193, 0.099, 0.043, 0.016, 0.004), #9000
    (0.030, 0.149, 0.382, 0.670, 0.906, 1.000, 0.938, 0.768, 0.558, 0.364, 0.215, 0.113, 0.051, 0.020, 0.006), #9200
    (0.026, 0.132, 0.348, 0.629, 0.877, 1.000, 0.971, 0.823, 0.620, 0.420, 0.258, 0.143, 0.069, 0.028, 0.010), #9400
    (0.024, 0.126, 0.337, 0.616, 0.868, 1.000, 0.981, 0.839, 0.638, 0.437, 0.271, 0.153, 0.074, 0.031, 0.011), #9600
    (0.022, 0.116, 0.317, 0.592, 0.851, 1.000, 1.000, 0.872, 0.676, 0.472, 0.298, 0.172, 0.087, 0.037, 0.014, 0.002), #9800
    (0.020, 0.106, 0.294, 0.561, 0.822, 0.983, 1.000, 0.888, 0.700, 0.498, 0.320, 0.188, 0.099, 0.043, 0.017, 0.004), #10000
    (0.017, 0.096, 0.272, 0.529, 0.790, 0.965, 1.000, 0.905, 0.727, 0.526, 0.346, 0.207, 0.113, 0.050, 0.020, 0.006), #10200
    (0.015, 0.087, 0.251, 0.499, 0.761, 0.946, 1.000, 0.922, 0.755, 0.556, 0.373, 0.227, 0.126, 0.061, 0.024, 0.008), #10400
    (0.014, 0.083, 0.242, 0.486, 0.747, 0.937, 1.000, 0.930, 0.768, 0.570, 0.385, 0.237, 0.134, 0.065, 0.026, 0.009), #10600
    (0.013, 0.075, 0.225, 0.459, 0.720, 0.920, 1.000, 0.947, 0.796, 0.602, 0.415, 0.260, 0.149, 0.075, 0.032, 0.012, 0.001), #10800
    (0.012, 0.069, 0.208, 0.435, 0.695, 0.904, 1.000, 0.963, 0.824, 0.633, 0.443, 0.284, 0.165, 0.085, 0.037, 0.015, 0.002), #11000
    (0.010, 0.063, 0.194, 0.412, 0.669, 0.888, 1.000, 0.980, 0.852, 0.667, 0.475, 0.309, 0.184, 0.098, 0.044, 0.018, 0.005), #11200
    (0.009, 0.057, 0.180, 0.391, 0.646, 0.872, 1.000, 0.997, 0.882, 0.702, 0.509, 0.336, 0.204, 0.113, 0.052, 0.021, 0.006), #11400
    (0.009, 0.054, 0.173, 0.379, 0.631, 0.861, 0.995, 1.000, 0.892, 0.717, 0.523, 0.350, 0.214, 0.119, 0.057, 0.023, 0.008), #11600
    (0.008, 0.049, 0.160, 0.355, 0.602, 0.834, 0.980, 1.000, 0.906, 0.739, 0.548, 0.373, 0.231, 0.132, 0.066, 0.026, 0.010), #11800
    (0.007, 0.042, 0.141, 0.321, 0.557, 0.791, 0.953, 1.000, 0.931, 0.781, 0.596, 0.417, 0.268, 0.158, 0.082, 0.037, 0.014, 0.002), #12000
    (0.006, 0.038, 0.130, 0.301, 0.531, 0.767, 0.939, 1.000, 0.945, 0.805, 0.624, 0.443, 0.289, 0.174, 0.093, 0.043, 0.017, 0.004), #12200
    (0.005, 0.035, 0.120, 0.283, 0.507, 0.744, 0.925, 1.000, 0.960, 0.830, 0.653, 0.470, 0.312, 0.191, 0.106, 0.051, 0.020, 0.006), #12400
    (0.005, 0.033, 0.115, 0.274, 0.495, 0.732, 0.918, 1.000, 0.967, 0.842, 0.668, 0.485, 0.324, 0.200, 0.112, 0.054, 0.023, 0.007), #12600
    (0.004, 0.030, 0.107, 0.257, 0.472, 0.710, 0.904, 1.000, 0.982, 0.868, 0.699, 0.515, 0.351, 0.219, 0.126, 0.063, 0.027, 0.010), #12800
    (0.004, 0.027, 0.098, 0.242, 0.450, 0.689, 0.890, 1.000, 0.997, 0.894, 0.731, 0.547, 0.378, 0.241, 0.141, 0.072, 0.032, 0.012, 0.002), #13000
    (0.003, 0.025, 0.090, 0.224, 0.426, 0.661, 0.867, 0.989, 1.000, 0.911, 0.756, 0.574, 0.402, 0.260, 0.155, 0.082, 0.037, 0.014, 0.003), #13200
    (0.003, 0.022, 0.082, 0.208, 0.402, 0.633, 0.843, 0.975, 1.000, 0.925, 0.777, 0.598, 0.425, 0.279, 0.169, 0.092, 0.043, 0.017, 0.005), #13400
    (0.003, 0.021, 0.079, 0.202, 0.392, 0.621, 0.833, 0.969, 1.000, 0.930, 0.786, 0.609, 0.435, 0.288, 0.176, 0.097, 0.046, 0.018, 0.006), #13600
    (0.003, 0.019, 0.073, 0.188, 0.370, 0.595, 0.810, 0.955, 1.000, 0.943, 0.808, 0.634, 0.460, 0.309, 0.191, 0.108, 0.053, 0.022, 0.007), #13800
    (0.002, 0.017, 0.067, 0.175, 0.350, 0.570, 0.787, 0.942, 1.000, 0.956, 0.831, 0.662, 0.487, 0.331, 0.209, 0.121, 0.062, 0.026, 0.010), #14000
    (0.002, 0.016, 0.061, 0.163, 0.330, 0.547, 0.765, 0.929, 1.000, 0.968, 0.855, 0.690, 0.515, 0.356, 0.227, 0.135, 0.070, 0.031, 0.012, 0.002), #14200
    (0.002, 0.014, 0.056, 0.151, 0.312, 0.524, 0.743, 0.916, 1.000, 0.982, 0.878, 0.718, 0.544, 0.382, 0.247, 0.149, 0.079, 0.037, 0.014, 0.003), #14400
    (0.002, 0.013, 0.054, 0.146, 0.304, 0.514, 0.733, 0.909, 1.000, 0.989, 0.890, 0.733, 0.559, 0.395, 0.257, 0.156, 0.084, 0.039, 0.016, 0.004), #14600
    (0.001, 0.012, 0.047, 0.131, 0.276, 0.478, 0.697, 0.881, 0.989, 1.000, 0.920, 0.777, 0.605, 0.437, 0.292, 0.182, 0.102, 0.051, 0.022, 0.007), #14800
    (0.001, 0.010, 0.043, 0.121, 0.259, 0.454, 0.671, 0.859, 0.977, 1.000, 0.932, 0.797, 0.629, 0.460, 0.312, 0.197, 0.114, 0.058, 0.025, 0.008, 0.001), #15000
)


"""from mMass.mspy import mod_peakpicking as peak_picking
import numpy as np
import pandas as pd
import os
rows=list()
ref_loc=input_table['Location'][0]
peaknumber=0
for n,loc in enumerate(input_table['Location']):
    x_data=np.array(input_table['XIC_RT'][n])
    y_data=np.array(input_table['XIC_intens'][n])
    signal=np.dstack((x_data,y_data))[0]
    peaklist=peak_picking.labelscan(signal,relThreshold=0.05)
    if loc==ref_loc:peaknumber+=1
    else:
        ref_loc=loc
        peaknumber=1
    for peak in peaklist:
        rows.append([loc,peak.mz,peak.ai,round(peak.fwhm*60,0),input_table['mz(theoretical)'][n]])
    
output_table = pd.DataFrame(rows, columns=['Location','RT','Intens','FWHM','mz(t)'])
"""
    