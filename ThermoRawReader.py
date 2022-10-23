# -*- coding: utf-8 -*-
"""
Created on Wed Oct 24 14:25:11 2018

@author: KSLC065r
"""
import comtypes
from comtypes.client import GetModule, CreateObject
from ctypes import c_long,c_double,byref,cast, POINTER
import bisect
import numpy as np
import pandas as pd
from SignalProcessing import moveAverage
import os.path as Path

class ThermoRawReader():
    """
    This class creates an object which has functions to parse information of 
    Thermo RAW files.
    An instance is created by providing the path the Thermo RAW file is saved   
    """
    
    def __init__(self,path,color=None,condition=None):
        GetModule("./XRawfile2.dll")
        self.path=path
        self.obj = CreateObject('XRawfile.XRawfile')
        self.obj.Open(path)
        self.obj.SetCurrentController(c_long(0), c_long(1))
        self.minStartTime=c_double()
        self.maxStartTime=c_double()
        self.obj.GetStartTime(byref(self.minStartTime))
        self.obj.GetEndTime(byref(self.maxStartTime))
        self.minStartTime=self.minStartTime.value
        self.maxStartTime=self.maxStartTime.value
        self.color=color
        self.condition=condition
        self.baseName=Path.basename(path)

    
    def __del__(self):
        self.obj.Close()
        
    def close(self):
        """
        Use this function to close a raw file.  
        """
        
        self.obj.Close()


    def getMetaData(self):
        """
        This function returns a tuple containing information about:
            Instrument method
            Injection volume
            Date     
            Number of scans
            Number of MS1 scans
            Number of MS2 scans
            Number of MS3 scans
        """
        import datetime
        pbstrInstrumentMethod=comtypes.BSTR()
        pdInjVol=c_double()
        pCreationDate=c_double()
        self.obj.GetSeqRowInstrumentMethod(byref(pbstrInstrumentMethod))
        self.obj.GetSeqRowInjectionVolume(pdInjVol)
        self.obj.GetCreationDate(byref(pCreationDate))
        serial=pCreationDate.value
        seconds=(serial-25569)*86400
        try:
            date=datetime.datetime.utcfromtimestamp(seconds)
            date=date.isoformat()
        except:
            date=""
        scanInformation=self.getScanInformation()
        scanNumbers=len(scanInformation)
        a=scanInformation.groupby("MassOrder").apply(len)
        numberMS1=0
        numberMS2=0
        numberMS3=0
        for n,item in enumerate(a):
            if a.index[n]==1:
                numberMS1=item
            elif a.index[n]==2:
                numberMS2=item
            elif a.index[n]==3:
                numberMS3=item
        
        return (pbstrInstrumentMethod.value,pdInjVol.value,date,
                scanNumbers,numberMS1,numberMS2,numberMS3)

    def getScanInformation(self):
        """
        This function returns a pandas dataframe containing following columns:
            scanNumber
            Massorder (MS1,M2,...)
            Starttime
        """
        pnMassOrder=c_long()
        numOfSpectra=c_long()
        rows=list()
        self.obj.GetNumSpectra(numOfSpectra)
        for scan in range(1,numOfSpectra.value+1):
            pnScanNumber=c_long(scan)
            self.obj.GetMSOrderForScanNum(pnScanNumber,pnMassOrder)
            massOrder=pnMassOrder.value
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
            self.obj.GetScanHeaderInfoForScanNum(pnScanNumber, numPackets, StartTime,
                                                        LowMass, HighMass,
                                                        TIC, BasePeakMass,
                                                        BasePeakIntensity, numChannels,
                                                        uniformTime, Frequency)
            rows.append([scan,massOrder,StartTime.value])
        df=pd.DataFrame(rows,columns=['scanNumber','MassOrder','StartTime'])
        df.name="ScanInformation@"+self.baseName
        return df

    def getPrecursorInformation(self):
        """
        This function returns a pandas dataframe containing information about all tandem mass
        spectra in one run. The dataframe includes following columns:
            scanNumber
            PrecursorMass
            StartTime
            TIC            
        """
        pnMassOrder=c_long()
        numOfSpectra=c_long()
        rows=list()
        self.obj.GetNumSpectra(numOfSpectra)
        for scan in range(1,numOfSpectra.value+1):
            pnScanNumber=c_long(scan)
            self.obj.GetMSOrderForScanNum(pnScanNumber,pnMassOrder)
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
                self.obj.GetPrecursorMassForScanNum(pnScanNumber, pnMassOrder,byref(precursorMass))
                self.obj.GetScanHeaderInfoForScanNum(pnScanNumber, numPackets, StartTime,
                                                            LowMass, HighMass,
                                                            TIC, BasePeakMass,
                                                            BasePeakIntensity, numChannels,
                                                            uniformTime, Frequency)
                precursor_mz=precursorMass.value

                rows.append([scan,precursor_mz,StartTime.value,TIC.value])
        df=pd.DataFrame(rows,columns=['scanNumber','PrecursorMass','StartTime','TIC'])
        return df
      
         
    def getXIC(self,massrange,startTime=0,endTime=0,smooth=True,scanFilter="Full ms"):
        """
        Description of the parameters:
               name        type                description
               massrange   string              for example "201.01-203.4" or "201.01-203.4,502.1-503.2"
               startTime   float               
               endTime     float               l
               threshould  float               Error in ppm
        
        This function returns a pandas dataframe containing following columns:
            RT   retention time
            intensity
        """
        if startTime>endTime:
            raise IOError("startTime must be lower than endTime!")
        
        pnArraySize = c_long()
        scanFilter=scanFilter
        chroType1=0
        chroOperator=0
        chroType2=0
        delay=0.0
        smoothingType=0
        smoothingValue=0
        massRange2=""
        massrange=massrange
        pvarChroData = comtypes.automation.VARIANT()
        pvarPeakFlags = comtypes.automation.VARIANT()
        

        
        if startTime:
            if startTime<self.minStartTime:
                startTime=self.minStartTime

        if endTime:
            if endTime>self.maxStartTime:
                endTime=self.maxStartTime
                
        self.obj.GetChroData(chroType1,
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

        try:
            signal=np.dstack((pvarChroData[0][0],pvarChroData[0][1]))[0]
        except:
            df=pd.DataFrame()
            df['RT']=[startTime,endTime]
            df['intensity']=[0,0]
            return(df)
            
        if smooth:       
            signal=moveAverage(signal,0.1,2)
        
        df=pd.DataFrame()
        df['RT']=signal[:,0]
        df['intensity']=signal[:,1]
        return(df) 
 


        
    def getTIC(self,startTime=0,endTime=0,smooth=True):
        """
        Description of the parameters:
               name        type                description
               massrange   string              for example "201.01-203.4" or "201.01-203.4,502.1-503.2"
               startTime   float               
               endTime     float               l
               threshould  float               Error in ppm        
        
        This function returns a pandas dataframe containing following columns:
            RT          retention time
            intensity
        
        """
        if startTime>endTime:
            raise IOError("startTime must be lower than endTime!")
        
        pnArraySize = c_long()
        scanFilter="Full ms"
        chroType1=1
        chroOperator=0
        chroType2=1
        delay=0.0
        smoothingType=0
        smoothingValue=0
        massRange2=""
        massrange=""
        pvarChroData = comtypes.automation.VARIANT()
        pvarPeakFlags = comtypes.automation.VARIANT()
        
        
        if startTime:
            if startTime<self.minStartTime:
                startTime=self.minStartTime

        if endTime:
            if endTime>self.maxStartTime:
                endTime=self.maxStartTime     
                
        self.obj.GetChroData(chroType1,
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

        signal=np.dstack((pvarChroData[0][0],pvarChroData[0][1]))[0]
        if smooth:       
            signal=moveAverage(signal,0.1,2)
        df=pd.DataFrame()
        df['RT']=signal[:,0]
        df['intensity']=signal[:,1]
        return(df) 


    def getBPC(self,startTime=0,endTime=0,smooth=True):
        """
        Description of the parameters:
               name        type                description
               startTime   float               
               endTime     float 
               smooth      bool                if this is equal to True, a move average algorthm is applied
                                               with window=0.1 and cycles=2

        This function returns a pandas dataframe containing following columns:
            RT          retention time
            intensity       
        """
        if startTime>endTime:
            raise IOError("startTime must be lower than endTime!")
        
        pnArraySize = c_long()
        scanFilter="Full ms"
        chroType1=2
        chroOperator=0
        chroType2=0
        delay=0.0
        smoothingType=0
        smoothingValue=0
        massRange2=""
        massrange=""
        pvarChroData = comtypes.automation.VARIANT()
        pvarPeakFlags = comtypes.automation.VARIANT()
        
        if startTime:
            if startTime<self.minStartTime:
                startTime=self.minStartTime

        if endTime:
            if endTime>self.maxStartTime:
                endTime=self.maxStartTime  
                
        self.obj.GetChroData(chroType1,
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

        signal=np.dstack((pvarChroData[0][0],pvarChroData[0][1]))[0]
        if smooth:       
            signal=moveAverage(signal,0.1,2)
        df=pd.DataFrame()
        df['RT']=signal[:,0]
        df['intensity']=signal[:,1]
        return(df) 

    def getUV(self):
        """
        This function returns a pandas dataframe containing following columns:
            RT          retention time
            intensity  
       
        """              
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
        self.obj.SetCurrentController(c_long(4), c_long(1))
        pvarChroData = comtypes.automation.VARIANT()
        pvarPeakFlags = comtypes.automation.VARIANT()
        
        self.obj.GetChroData(chroType1,
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
        self.obj.SetCurrentController(c_long(0), c_long(1))
        df=pd.DataFrame()
        df['RT']=UV_RT
        df['intensity']=UV_intens
        return(df)

    def getMS1spectraByScanNumber(self,scanNumber):
        """
        Define the scan number of the MS1 spectra
        This function returns a pandas dataframe containing following columns:
            mz          
            intensity  
        
        """
        
        scanFilter="Full ms"
        intensityCufoffType=0
        intensityCufoffValue=0
        nMaxNumberOfPeaks=0
        CentroidResult=False
        centroidPeakwidth=c_double(0.0)
        pvarMassList=comtypes.automation.VARIANT()
        pvarPeakFlags = comtypes.automation.VARIANT()
        pnArraySize = c_long()
        pnScanNumber=c_long(scanNumber)
        self.obj.GetMassListFromScanNum(pnScanNumber,
                          scanFilter,
                          intensityCufoffType,
                          intensityCufoffValue,
                          nMaxNumberOfPeaks,
                          CentroidResult,
                          centroidPeakwidth,                      
                          pvarMassList,
                          pvarPeakFlags,
                          byref(pnArraySize))
        if not pvarMassList.value:
            df=pd.DataFrame(columns=['mz','intensity'])
            return(df)
        df=pd.DataFrame()
        df['mz']=np.array(pvarMassList[0][0])
        df['intensity']=np.array(pvarMassList[0][1])
        return(df) 


    def getMS1spectraByRT(self,RT):
        """
        Define the retention time in minutes
        This function returns a pandas dataframe containing following columns:
            mz          
            intensity  
        
        """        
        if RT<self.minStartTime:
            RT=self.minStartTime

        if RT>self.maxStartTime:
            RT=self.maxStartTime  
        scanFilter="Full ms"
        intensityCufoffType=0
        intensityCufoffValue=0
        nMaxNumberOfPeaks=0
        CentroidResult=False
        centroidPeakwidth=c_double(0.0)
        pvarMassList=comtypes.automation.VARIANT()
        pvarPeakFlags = comtypes.automation.VARIANT()
        pnArraySize = c_long()
        pdRT=c_double(RT)
        self.obj.GetMassListFromRT(pdRT,
                          scanFilter,
                          intensityCufoffType,
                          intensityCufoffValue,
                          nMaxNumberOfPeaks,
                          CentroidResult,
                          centroidPeakwidth,                      
                          pvarMassList,
                          pvarPeakFlags,
                          byref(pnArraySize))
        if not pvarMassList.value:
            return(0,0,0)
        df=pd.DataFrame()
        df['mz']=np.array(pvarMassList[0][0])
        df['intensity']=np.array(pvarMassList[0][1])
        return(df) 

    
    def getScanNumberFromRT(self,RT,MSorder=1):
        """
        This function returns the scan number of a defined retention time
        in minutes and MS order
        
        """    
        
        if RT>self.maxStartTime:
            raise IOError("Retention time is too high!")
        if RT<self.minStartTime:
            raise IOError("Retention time is too high!")
        pnScanNumber=c_long()
        pnMassOrder=c_long()
        self.obj.ScanNumFromRT(RT,pnScanNumber)
        self.obj.GetMSOrderForScanNum(pnScanNumber,pnMassOrder)
        while not pnMassOrder.value==MSorder:
            pnScanNumber=c_long(pnScanNumber.value+1)
            self.obj.GetMSOrderForScanNum(pnScanNumber,pnMassOrder)
        
        return(pnScanNumber.value)


    def getRTFromScanNumber(self,scanNumber):
        """
        This function returns the scan number of a defined retention time
        in minutes and MS order
        
        """    
        pnScanNumber=c_long(scanNumber)
        pdRT=c_double()
        self.obj.RTFromScanNum(pnScanNumber,pdRT)       
        return(pdRT.value)        

    def getMSLabelByScanNumber(self,scanNumber):
        """
        This function returns a pandas DataFrame containing following MS labels:
            mz
            intensity
            resolution
            Base
            noise
            charge
        
        Define the scanNumber of the spectra
        
        """  
        pnScanNumber=c_long(scanNumber)
        pvarLabels=comtypes.automation.VARIANT()
        pvarFlags = comtypes.automation.VARIANT()
  
        self.obj.GetLabelData(pvarLabels,
                          pvarFlags,
                          pnScanNumber)
        df=pd.DataFrame()
        df['mz']=np.array(pvarLabels[0][0])
        df['intensity']=np.array(pvarLabels[0][1])
        df['resolution']=np.array(pvarLabels[0][2])
        df['Base']=np.array(pvarLabels[0][3])
        df['noise']=np.array(pvarLabels[0][4])
        df['charge']=np.array(pvarLabels[0][5])
        
        return(df)
        
    def getMS2LabelByScanNumber(self,scanNumber):
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
        self.obj.GetMassListFromScanNum(pnScanNumber,
                                   szFilter, 
                                   intensityCutoffType,
                                   intensityCutoffValue,
                                   nMaxNumberOfPeaks,
                                   bCetroidResult,
                                   c_double(centroidPeakWidth),
                                   pvarMassList, 
                                   pvarFlags,
                                   byref(pnArraySize))
    
        

        df=pd.DataFrame()
        try:
            df['mz']=np.array(pvarMassList[0][0])
            df['intensity']=np.array(pvarMassList[0][1])
            return(df)
        except:                      
            df['mz']=[]
            df['intensity']=[]
            return(df)
            
        
       


    
    def getMSLabelByRT(self,RT,MSorder=1):
        """
        This function returns a pandas DataFrame containing following MS labels:
            mz
            intensity
            resolution
            Base
            noise
            charge
        
        Define the retention time (RT) in minutes and the MSorder        
        """  
        scanNumber=self.getScanNumberFromRT(RT,MSorder)
        df=self.getMSLabelByScanNumber(scanNumber)
        return(df,scanNumber)
        
    


    def getFilters(self):
        """
        This function returns a list of filters       
        """

        pvarFilterArray = comtypes.automation.VARIANT()
        pnArraySize = c_long()
        self.obj.GetFilters(byref(pvarFilterArray), byref(pnArraySize))
        return(pvarFilterArray.value)
        
    def getSummedMassSpectra(self,listOfScanNumbers,centroidResult=False):
        """
        This function returns a list of filters       
        """
        
        x = (c_long * len(listOfScanNumbers))()
        cast(x, POINTER(c_long))
        for n, scannumber in enumerate(listOfScanNumbers):
            x[n]=scannumber
        
        peakList = comtypes.automation.VARIANT()
        peakFlags = comtypes.automation.VARIANT()
        pnArraySize = c_long()
        self.obj.GetSummedMassSpectrum(x,
                                       len(x),
                                       centroidResult,
                                       byref(peakList),
                                       byref(peakFlags),
                                       byref(pnArraySize)) 
        
        df=pd.DataFrame()
        df['mz']=np.array(peakList[0][0])
        df['intensity']=np.array(peakList[0][1])
        return(df) 
        
        
        
    def getMS2ScanNumber(self,ls_mz,ls_RT,ls_charge=[],ppm=5e-6):
        """This function returns a list of potential tandem mass spectra for a list
           of mz- and RT-values.
           Input variable:
               name        type                description
               ls_mz       list                list of mz-values
               ls_RT       list                list of retention time values
               ls_charge   list(optional)      list of charge values
               ppm         float               Error in ppm        
        """
        df_MS2=self.getPrecursorInformation()
        ls_scanNumber=list()
        for n,mz in enumerate(ls_mz):
            RT=ls_RT[n]
            temp_df=df_MS2[(df_MS2['StartTime']>RT-3)&(df_MS2['StartTime']<RT+3)]
            if ls_charge:
                charge=ls_charge[n]
                result_df=temp_df[df_MS2['precursorCharge']==charge]
            else:
                result_df=temp_df
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



    def getDetailedPrecursorInformation(self):
        """
        This function returns a pandas dataframe containing information about all tandem mass
        spectra in one run. The dataframe includes following columns:
            scanNumber
            PrecursorMass
            StartTime
            TIC
            precursorIntensity
            precursorCharge
            observedMZ               
        """
        pnMassOrder=c_long()
        numOfSpectra=c_long()
        rows=list()
        self.obj.GetNumSpectra(numOfSpectra)
        for scan in range(1,numOfSpectra.value+1):
            pnScanNumber=c_long(scan)
            self.obj.GetMSOrderForScanNum(pnScanNumber,pnMassOrder)
            if pnMassOrder.value==1:
                parentScanNumber=pnScanNumber
                parentMSLabel=comtypes.automation.VARIANT()
                pvarFlags = comtypes.automation.VARIANT()
                self.obj.GetLabelData(parentMSLabel,
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
                self.obj.GetPrecursorMassForScanNum(pnScanNumber, pnMassOrder,byref(precursorMass))
                self.obj.GetScanHeaderInfoForScanNum(pnScanNumber, numPackets, StartTime,
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
        df.name="PrecursorList@"+self.baseName
        return df