# -*- coding: utf-8 -*-
"""
Created on Sat May  4 10:51:27 2019

@author: KTMS458
"""
import numpy as np
from mMass.mspy import mod_peakpicking as peak_picking
#from mMass.mspy import mod_signal 
from pyqms.isotopologue_library import IsotopologueLibrary
from pyqms.params import params
import bisect
from mMass.mspy import obj_sequence
from  mMass.mspy.mod_proteo import fragment as fragment
from  mMass.mspy.mod_proteo import fragmentlosses as fragmentlosses
import pandas as pd
import re



       
class Component():
    def __init__(self,sequence,mz,charge,ThermoRawFile=None, modification=None,RT=0,name=None):
        self.rawfile_dict=dict()
        self.name=str(name)
        self.mz=mz
        self.RT=RT
        self.sequence=sequence
        self.charge=charge
        self.modification=modification
        # [[name, position=[#|symbol], state=[f|v]], ] (f->fixed, v->variable)
        self.PyqmsSequence=self.sequence
        self.mMassSequence=obj_sequence.sequence(chain=self.sequence,agentFormula="H")        
        if self.modification:
            from mMass.mspy import blocks
            blocks.loadModifications()
            for mod in self.modification:
                if re.match("[Na|K]\+",mod[0]):
                    self.mMassSequence.agentFormula=mod
                    continue                
                self.mMassSequence.modify(mod[0],mod[1],mod[2])
                mod=mod[0]
                mod=blocks.modifications[mod]
                if mod.gainFormula:self.PyqmsSequence+='+'+mod.gainFormula
                if mod.lossFormula:self.PyqmsSequence+='-'+mod.lossFormula
                    
        self.mMassSequence=obj_sequence.sequence(chain=self.sequence,)
        self.compInRawfile=False
        self.cc=self.getPyqms()
        if ThermoRawFile:
            self.setRawFile(ThermoRawFile)

    def getPyqms(self):
        cc=IsotopologueLibrary(molecules=[self.PyqmsSequence],charges=[self.charge],params=params)
        key=cc.keys()[0]
        self.chemicalComposition=cc[key]['cc']
        self.STRchemicalComposition=key
        self.theoIsotopologues=np.dstack((cc[key]["env"][(('N', '0.000'),)][self.charge]['mz'],cc[key]["env"][(('N', '0.000'),)]['relabun']))[0]
        self.mzIsotopologueHighestIntensity=self.theoIsotopologues[:,0][self.theoIsotopologues[:,1].argmax()]
        self.monoisotopicMass=cc[key]["env"][(('N', '0.000'),)]['mass'][0]
        return(cc)    
    
    def setRawFile(self,rawfile,RTwindow=2,mz_error=10e-6,flag="highest"):
        self.rawfile=rawfile
        self.rawfileBaseName=self.rawfile.baseName
        color=self.rawfile.color
        condition=self.rawfile.condition
        compInRawfile,realRT,matchedIsotopologues,highestIntensityIsotopologue,scanNumber,score=self.getRawFileAttributes(RTwindow,mz_error,flag)
        self.rawfile_dict[self.rawfileBaseName]={'compInRawfile':compInRawfile,
                                                 'realRT':realRT,
                                                 'matchedIsotopologues':matchedIsotopologues,
                                                 'highestIntensityIsotopologue':highestIntensityIsotopologue,
                                                 'scanNumber':scanNumber,
                                                 'MScore':score,
                                                 'location':rawfile.path,
                                                 'color':color,
                                                 'condition':condition,
                                                 'quant':{},
                                                 'MS2scanNumber':0}

        self.getQuantValues()



    
    def getRawFileAttributes(self,RTwindow=2,mz_error=10e-6,flag="highest"):
        massrange="%f-%f"%(self.mz-mz_error*self.mz,self.mz+mz_error*self.mz)
        startTime=self.RT-RTwindow
        endTime=self.RT+RTwindow
        XIC=self.rawfile.getXIC(massrange,startTime,endTime,smooth=False,scanFilter="Full ms")
        
        signal=np.dstack((XIC['RT'],XIC['intensity']))[0]
        peaklist=peak_picking.labelscan(signal,relThreshold=0.01)
        emptyPeaklist=np.array([[0,0],[0,0]])
                
        if len(peaklist)==0:
            return(0,0,emptyPeaklist,0,0,0)
        ls_peakRT=np.array([])
        ls_peakInt=np.array([])
        for peak in peaklist:
            ls_peakRT=np.append(ls_peakRT,peak.mz)
            ls_peakInt=np.append(ls_peakInt,peak.ai)
        
        peaklist=np.dstack((ls_peakRT,ls_peakInt))[0]
        if flag=="highest":
            realRT=peaklist[:,0][peaklist[:,1].argmax()]
        elif flag=="closest":
            realRT=peaklist[(abs(peaklist[:,0]-self.RT)).argmin()]
        else:
            return(0,0,emptyPeaklist,0,0,0)
        scanNumber=self.rawfile.getScanNumberFromRT(realRT)
        MSLabel=self.rawfile.getMSLabelByScanNumber(scanNumber)
        MSPeaklist=np.dstack((MSLabel['mz'],MSLabel['intensity']))[0]
        results=self.cc.match_all(MSPeaklist)
        if results:
            matchedIsotopologues,highestIntensityIsotopologue,score=self.getMatchedPeaklist(results)
            if score>=0.5:
                realRT=self.rawfile.getRTFromScanNumber(scanNumber)
                compInRawfile=True
                return(compInRawfile,realRT,matchedIsotopologues,highestIntensityIsotopologue,scanNumber,score)
            
        #sort peaklist descending
        peaklist=peaklist[peaklist[:,1].argsort()[::-1]]
        
        for RT in peaklist[:,0]:
            MSLabel,scanNumber=self.rawfile.getMSLabelByRT(RT)
            MSPeaklist=np.dstack((MSLabel['mz'],MSLabel['intensity']))[0]
            results=self.cc.match_all(MSPeaklist)
            if results:
                matchedIsotopologues,highestIntensityIsotopologue,score=self.getMatchedPeaklist(results)
                if score>=0.5:
                    realRT=self.rawfile.getRTFromScanNumber(scanNumber)
                    compInRawfile=True
                    return(compInRawfile,realRT,matchedIsotopologues,highestIntensityIsotopologue,scanNumber,score)


        return(0,0,emptyPeaklist,0,0,0)
            
    def getMatchedPeaklist(self,results):
        matchedPeaks=results[results.keys()[0]]['data'][0].peaks
        score=results[results.keys()[0]]['data'][0].score
        rows=list()
        highestIntensity=0
        for mmz, mi, rel_i, cmz, ci in matchedPeaks:
            if mi>highestIntensity:
                highestIntensity=mi
            row=[mmz,mi]
            rows.append(row)
        
        matchedIsotopologues=np.array(rows)
        highestIntensityIsotopologues=highestIntensity
        return(matchedIsotopologues,highestIntensityIsotopologues,score)          
  
    def assignFragments(self,MS2spectra,fragments):
        ls_ionMatch=list()
        ls_ionMatchIntens=list()
        for ion in fragments['mz']:
            minMZ=ion-0.1
            maxMZ=ion+0.1
            minpos=bisect.bisect(list(MS2spectra['mz']),minMZ)
            maxpos=bisect.bisect(list(MS2spectra['mz']),maxMZ)
            candidates=MS2spectra['mz'][minpos:maxpos]
            if len(candidates)==0:
                ls_ionMatch.append(-1)
                ls_ionMatchIntens.append(-1)
                continue
            maxIonIntens=max(MS2spectra['intensity'][minpos:maxpos])
            ionMatch=float(candidates[MS2spectra['intensity'][minpos:maxpos]==max(MS2spectra['intensity'][minpos:maxpos])])
            ls_ionMatch.append(ionMatch)
            ls_ionMatchIntens.append(maxIonIntens)
        fragments['ionMatches']=ls_ionMatch
        fragments['ionMatchIntens']=ls_ionMatchIntens
        return(fragments)    

    
    def getFragments(self,charge=1,series=['y','b'],fragment_losses=True,agentFormula='H'):
        sequence=self.mMassSequence
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
        
    def getQuantValues(self,RTwindow=1): 
        key=self.rawfileBaseName               
        if not self.rawfile_dict[key]['compInRawfile']:
            RT=self.RT
            startTime=RT-RTwindow
            endTime=RT+RTwindow                                    
            self.rawfile_dict[key]['quant']['PYQMS']={'value':0,
                                                     'intensity':[0,0],
                                                     'RT':[startTime,endTime],
                                                     'limit':[startTime,endTime]}
            

            quantValue,intensity,RT=self.getQuantMXIC(startTime,endTime)              
            self.rawfile_dict[key]['quant']['mXIC']={'value':quantValue,
                                                 
                                                 'intensity':intensity,
                                                 'RT':RT,
                                                 'limit':[startTime,endTime]}
            
            quantValue,intensity,RT=self.getQuantHXIC(startTime,endTime)  
            self.rawfile_dict[key]['quant']['hXIC']={'value':quantValue,
                                                     'intensity':intensity,
                                                     'RT':RT,
                                                     'limit':[startTime,endTime]}
            return
            
        
        RT=self.rawfile_dict[key]['realRT']
        if not RT:
            RT=self.RT
        startTime=RT-RTwindow
        endTime=RT+RTwindow
        
        quantValue,intensity,RT=self.getQuantMXIC(startTime,endTime)              
        self.rawfile_dict[key]['quant']['mXIC']={'value':quantValue,
                                                 'intensity':intensity,
                                                 'RT':RT,
                                                 'limit':[startTime,endTime]}
        
        quantValue,intensity,RT=self.getQuantHXIC(startTime,endTime)  
        self.rawfile_dict[key]['quant']['hXIC']={'value':quantValue,
                                                 'intensity':intensity,
                                                 'RT':RT,
                                                 'limit':[startTime,endTime]}
            

        quantValue,intensity,RT=self.getQuantHXIC(startTime,endTime)  
        self.rawfile_dict[key]['quant']['PYQMS']={'value':quantValue,
                                                 'intensity':intensity,
                                                 'RT':RT,
                                                 'limit':[startTime,endTime]}
        




    def getQuantMXIC(self,startTime,endTime,mz_error=10e-6):
        key=self.rawfile.baseName
        try:
            mz=self.rawfile_dict[key]['matchedIsotopologues'][:,0][0]
            if not mz:
                mz=self.theoIsotopologues[:,0][0] 
        except:
            mz=self.theoIsotopologues[:,0][0]           
        print startTime
        print endTime        
        massrange="%f-%f"%(mz-mz_error*mz,mz+mz_error*mz)
        XIC=self.rawfile.getXIC(massrange,startTime,endTime,smooth=False,scanFilter="Full ms")                
        quantValue=np.trapz(XIC['intensity'],XIC['RT'])
        return(quantValue,list(XIC['intensity']),list(XIC['RT']))


    def getQuantHXIC(self,startTime,endTime,mz_error=10e-6):
        key=self.rawfile.baseName
        try:
            mz=self.rawfile_dict[key]['matchedIsotopologues'][:,0][self.rawfile_dict[key]['matchedIsotopologues'][:,1].argmax()]
        except:
            mz=self.mzIsotopologueHighestIntensity     
        massrange="%f-%f"%(mz-mz_error*mz,mz+mz_error*mz)

        XIC=self.rawfile.getXIC(massrange,startTime,endTime,smooth=False,scanFilter="Full ms")                
        quantValue=np.trapz(XIC['intensity'],XIC['RT'])
        return(quantValue,list(XIC['intensity']),list(XIC['RT']))
        
    def getQuantPYQMS(self,startTime,endTime,mz_error=5e-6):
        scan_df=self.rawfile.getScanInformation()
        scan_df=scan_df[(scan_df['MassOrder']==1)&(scan_df['StartTime']>=startTime)&(scan_df['StartTime']<=endTime)]
        ls_RT=list()
        ls_intensity=list()
        for n,row in scan_df.iterrows():
            scanNumber=int(row['scanNumber'])
            RT=row['StartTime']
            MSLabel=self.rawfile.getMSLabelByScanNumber(scanNumber)
            MSPeaklist=np.dstack((MSLabel['mz'],MSLabel['intensity']))[0]
            results=self.cc.match_all(MSPeaklist)
            if results:
                matchedIsotopologues,highestIntensityIsotopologues,score=self.getMatchedPeaklist(results)
                if score >= 0.7:
                    ls_RT.append(RT)
                    ls_intensity.append(highestIntensityIsotopologues)
                else:
                    ls_RT.append(RT)
                    ls_intensity.append(0)
            else:
                ls_RT.append(RT)
                ls_intensity.append(0)
        
        quantValue=np.trapz(ls_intensity,ls_RT)
        return(quantValue,ls_intensity,ls_RT)
                
        
    def prepareForPickle(self,cc=True,rawfile=True):
        if cc: 
            del self.cc
            self.cc=None
        if rawfile: 
            self.rawfile.close()
            self.rawfile=None
        
        