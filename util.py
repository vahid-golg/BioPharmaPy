# -*- coding: utf-8 -*-
"""
Created on Wed Apr 17 14:08:15 2019

@author: KTMS458
"""
import pandas as pd
import re
from . import ThermoRawReader.ThermoRawReader
from . import Graphs4.Graphs as graph
from io import BytesIO
from PIL import Image

def getPeptideNames(input_table_1,input_table_2):
    """"
    input_table_1 should be a pandas DataFrame containing two columns acc and sequence
    input_table_2 should be a pandas DataFrame containing one column Peptide Sequence
    
    """
    from mMass.mspy import obj_sequence as seq
    from mMass.mspy.mod_proteo import digest
    import re
    d={}
    for n,item in input_table_1.iterrows():
        sequence=item['sequence']
        sequence=seq.sequence(sequence)
        acc=item['acc']
        peptides=digest(sequence,"Trypsin")
        d[acc]=peptides
    
    
    rows=list()
    for acc in d.keys():
        end=0
        numerator=1
        for pep in d[acc]:
            result="".join(pep.chain)
            begin=end+1
            end=begin+len(result)-1
            name="T%i"%(numerator)
            numerator+=1
            rows.append([acc,result,begin,end,name])
            
    reference_table = pd.DataFrame(data=rows,columns=['acc',
    				'Peptide Sequence','begin','end','name'])
    results=list()
    for pep in input_table_2['Peptide Sequence']:
        output_name=""
        ls_output=list()
        for n,row in input_table_1.iterrows():
            prot_seq=row['sequence']
            prot_name=row['acc']
            a=re.finditer(pep,prot_seq)
            output_name=""
            for item in a:
                output_name=""
                begin,end=item.span()
                buff=reference_table[reference_table['acc']==prot_name]
                buff=buff[(buff['begin']>=begin)&(buff['end']<=end)]
                if len(buff)==0:
                    continue
                for name in buff['name']:
                    output_name+=name+"+"
                output_name=output_name[:-1]
                ls_output.append(prot_name+output_name)
                
        if ls_output:
           output_string=",".join(ls_output)
        else:
           output_string='None'
        results.append(output_string)
                
                
    output_table=input_table_2.copy()
    output_table['peptide_name']=results
    return(output_table)

def generatePepModQuant(input_table):
    from mMass.mspy import blocks
    from mMass.mspy import obj_sequence
    blocks.loadModifications("../mMass/configs/modifications.xml")
    ls_begin=list()
    ls_end=list()
    ls_pep_modification=list()
    
    for n,row in input_table.iterrows():
        begin_end=row['begin_end']
        pepideSequence=row['Peptide Sequence']
        #To check if it is a valide sequence, do not delete this line
        mMassSequence=obj_sequence.sequence(chain=pepideSequence,agentFormula="H")
        del mMassSequence
        modification=row['Modification']
        site=str(row['Site'])
        if not begin_end:
            ls_begin.append(0)
            ls_end.append(0)
            ls_pep_modification.append('None')
            continue
        begin_end=begin_end.split(';')[0].split(",")
        begin,end=int(begin_end[1]),int(begin_end[2])
        ls_begin.append(begin)
        ls_end.append(end)
    	
        ls_modName=list()
    	 #Carbamidomethylation
        indexC=[i+1 for i,residue in enumerate(pepideSequence) if residue=='C']
    
        for item in indexC:
    	    mod_name="Carbamidomethyl:C:%i"%(item)
    	    ls_modName.append(mod_name)
            
    
        if modification=="None":
            if ls_modName:
                mod_name=','.join(ls_modName)  
                ls_pep_modification.append(mod_name)
            else: 
                ls_pep_modification.append('None')
            continue
        modifications=modification.split(',')
        residue_sites=str(site.split(','))
        sites=re.findall("[\d+]|cTerm|nTerm",residue_sites)
        ls_sites=list()
        for site in sites:
            if site=="nTerm":
                ls_sites.append("nTerm")
            elif site=="cTerm":
                ls_sites.append("cTerm")
            else:
                ls_sites.append(str(int(site)-begin+1))
                       
        for n,mod in enumerate(modifications):
            if re.match("[Na|K]\+",mod):
                mod_name="%s:%s"%(mod,"0")
                ls_modName.append(mod_name)
                continue                
            blocks.modifications[mod]                
            mod_name="%s:%s"%(mod,sites[n])
            ls_modName.append(mod_name)    
        mod_name=','.join(ls_modName)  
        ls_pep_modification.append(mod_name)      
        


    output_table=input_table.copy()
    output_table['begin']=ls_begin
    output_table['end']=ls_end
    output_table['pep_modification']=ls_pep_modification
    return(output_table)

def getModList(ModificationString):
    ls_modification=list()
    modifications=ModificationString.split(',')
    if modifications[0]=='None': return(ls_modification)
    else:
        for mod in modifications:
            modification=mod.split(':')
            mod=modification[0]
            site=modification[1]
            if site=="nTerm" or site=="cTerm":
                pass
            else:
                site=int(site)
            modification=[mod,site,'f']
            ls_modification.append(modification)
    return(ls_modification)

def modTable(input_table):
    from mMass.mspy import blocks
    blocks.loadModifications("../mMass/configs/modifications.xml")
    grouped=input_table.groupby(['Site','Protein','Modification']).apply(lambda x: 0).reset_index(name='remove')
    rows=list()
    buff_dict=dict()
    for n,item in grouped.iterrows():
        site=item['Site']
        protein=item['Protein']
        mod=item['Modification']
        mods=mod.split(",")
        sites=site.split(",")
        for n,site in enumerate(sites):        
            mod=mods[n]
            if re.match("[Na|K]\+",mod) or mod=="None":
                continue
            if mod=="Gln->Pyro-Glu":
                mod="Gln2Pyro-Glu"
            buff=blocks.modifications[mod]
            if buff.description=="Glycan":
                key=site+'_'+protein
                if not key in buff_dict:
                    buff_dict[key]=[]
                buff_dict[key].append(mod)
                continue
            else:
                row=[protein,site,mod]
            rows.append(row)
    
    for key in buff_dict.keys():
        site,protein=key.split("_")
        output_string=""
        for item in buff_dict[key]:
            output_string+=item+","
        output_string=output_string[:-1]
        row=[protein,site,output_string]
        rows.append(row)
                
    output_table=pd.DataFrame(data=rows, columns=['protein','site','mod'])
    return output_table   
   

"""

Fragmentation

"""
def getMS2ScanNumber(DataFramePrecursors,component,mzError=10e-6,RTwindow=3):
    try:
        RT=component.rawfile_dict[DataFramePrecursors.name]['realRT']
    except:
        RT=component.RT
    charge=component.charge
    mz=component.mz
    result_df=DataFramePrecursors[(DataFramePrecursors['StartTime']>RT-3)
                                 &(DataFramePrecursors['StartTime']<RT+3)
                                 &(DataFramePrecursors['precursorCharge']==charge)
                                 &(DataFramePrecursors['PrecursorMass']>(mz-mz*mzError))
                                 &(DataFramePrecursors['PrecursorMass']<(mz+mz*mzError))].copy()
    if len(result_df)==0:
        result_df=DataFramePrecursors[(DataFramePrecursors['StartTime']>RT-3)
                                     &(DataFramePrecursors['StartTime']<RT+3)
                                     &(DataFramePrecursors['precursorCharge']==charge)
                                     &(DataFramePrecursors['observed_mz']>(mz-mz*mzError))
                                     &(DataFramePrecursors['observed_mz']<(mz+mz*mzError))].copy()
        
    
    if len(result_df)==0:
        return(0)
    elif len(result_df)==1:
        return(result_df['scanNumber'].iloc[0])
    else:
        scanNumber=result_df['scanNumber'][result_df['precursorIntens'].idxmax()]
    return(scanNumber)
            

def generateMS2DF(dfComponents,referenceFile=""):
    from scipy import spatial
    
    comp=dfComponents['component'][0]
    ls_rawfile=list()
    for item in comp.rawfile_dict.keys():
        loc=comp.rawfile_dict[item]['location']
        rawfile=ThermoRawReader(loc)
        ls_rawfile.append(rawfile)
    
    for rawfile in ls_rawfile:
         df_precursors=rawfile.getDetailedPrecursorInformation()
         baseName=rawfile.baseName
         for comp in dfComponents['component']:
             scanNumberMS2=getMS2ScanNumber(df_precursors,comp)
             comp.rawfile_dict[baseName]['MS2scanNumber']=scanNumberMS2
             if scanNumberMS2:
                  MS2Spectrum=rawfile.getMS2LabelByScanNumber(scanNumberMS2)
                  comp.rawfile_dict[baseName]['MS2Spectrum']=MS2Spectrum
         rawfile.close()
    
    rowsImage=list()
    rowsData=list()
    for n,row in dfComponents.iterrows():
        component=row['component']
        componentNumber=row['No.']
        rawfileDict=component.rawfile_dict
        fragments=component.getFragments()
        buffDict={}
        maxAssignment=(0,referenceFile)
        
        for item in rawfileDict.keys():
            MS2ScanNumber=rawfileDict[item]['MS2scanNumber']
            if not MS2ScanNumber:
                buffDict[item]={'MS2Spectrum':pd.DataFrame(columns=['mz','intensity']),
                                'fragmentsAssigned':0,
                                'numberAssigned':0,
                                'basename':baseName,
                                'scanNumber':0}
                continue

            MS2Spectrum=rawfileDict[item]['MS2Spectrum']           
            fragmentsAssigned=component.assignFragments(MS2Spectrum,fragments)
            numberAssigned=len(fragmentsAssigned[fragmentsAssigned['ionMatches']>0])            
            if numberAssigned>maxAssignment[0]:
                maxAssignment=(numberAssigned,item)                
            
            buffDict[item]={'MS2Spectrum':MS2Spectrum,
                            'fragmentsAssigned':fragmentsAssigned,
                            'numberAssigned':numberAssigned,
                            'basename':item,
                            'scanNumber':MS2ScanNumber}
        
        referenceFile=maxAssignment[1]
        if not referenceFile:
            referenceFile=rawfileDict.keys()[0]
        referenceMS2Spectrum=buffDict[referenceFile]['MS2Spectrum']
        referenceFragmentsAssigned=buffDict[referenceFile]['fragmentsAssigned']
        title="file:%s scan number:%s"%(buffDict[referenceFile]['basename'],buffDict[referenceFile]['scanNumber'])
        Graph=graph()
        Graph.getMS2fig(referenceMS2Spectrum,referenceFragmentsAssigned,title=title)                        
        fig=Graph.fig
        fig_buffer=BytesIO()
        fig.savefig(fig_buffer)
        img=Image.open(fig_buffer)
        rowsImage.append([componentNumber,img])
        
        for item in buffDict.keys():
            MS2spectrum=buffDict[item]['MS2Spectrum']
            ls1,ls2=alignSpectra(referenceMS2Spectrum,MS2spectrum)
            if not ls1:
                cos=0
            else:
                cos=1-spatial.distance.cosine(ls1,ls2)
            numberAssigned=buffDict[item]['numberAssigned']
            scanNumber=buffDict[item]['scanNumber']
            rowsData.append([componentNumber,item,numberAssigned,scanNumber,cos])
            


    dfImage=pd.DataFrame(data=rowsImage,columns=['No.','MS2Image'])
    dfData=pd.DataFrame(data=rowsData,columns=['No.','rawfile','scanNumber','NumberFragmentsAssigned','cos'])
    return (dfImage,dfData)
                           


def alignSpectra(df1,df2):
    if len(df1)==0 or len(df2)==0:
        return ([],[])
    import numpy as np
    mz1=np.array(df1['mz'])
    mz2=np.array(df2['mz'])
    int1=np.array(df1['intensity'])
    int2=np.array(df2['intensity'])
    peaklist1=np.dstack((mz1,int1))[0]
    peaklist2=np.dstack((mz2,int2))[0]
    peaklist_combined=np.concatenate((peaklist1,peaklist2))
    peaklist_combined=peaklist_combined[peaklist_combined[:,0].argsort()]
    ls_Spectrum1=list()
    ls_Spectrum2=list()
    xmax=0
    for item in peaklist_combined:
        if item[0]<=xmax: continue
        buff_min=item[0]-0.2
        buff_max=item[0]+0.2
        xmax=buff_max
        ls_Spectrum1.append(peaklist1[:,1][(peaklist1[:,0]>buff_min)&(peaklist1[:,0]<=buff_max)].sum())
        ls_Spectrum2.append(peaklist2[:,1][(peaklist2[:,0]>buff_min)&(peaklist2[:,0]<=buff_max)].sum())    
    return(ls_Spectrum1,ls_Spectrum2) 


def generateMS1DF(compDF,referenceFile=""):
    rows=list()
    for n,row in compDF.iterrows():
        comp=row['component']
        compNR=row['No.']
        rawfile_dict=comp.rawfile_dict
        if referenceFile:
            referenceFile=referenceFile
        else:
            referenceFile=rawfile_dict.keys()[0]
        
        MScore=rawfile_dict[referenceFile]['MScore']
        if MScore<0.5:
            max_MScore=MScore
            for item in rawfile_dict.keys():
                buff_MScore=rawfile_dict[item]['MScore']
                if buff_MScore>max_MScore:
                    referenceFile=item
        location=rawfile_dict[referenceFile]['location']
        rawfile=ThermoRawReader(location)
        scanNumber=rawfile_dict[referenceFile]['scanNumber']
        MSSpectrum=rawfile.getMS1spectraByScanNumber(scanNumber)
        title="file:%s scan number:%s"%(referenceFile,str(scanNumber))
        Graph=graph(comp)
        Graph.createMS1Fig(MSSpectrum,referenceFile,title=title)
        fig=Graph.fig
        fig_buffer=BytesIO()
        fig.savefig(fig_buffer)
        img=Image.open(fig_buffer)
        rows.append([compNR,img])
    
    imgDF=pd.DataFrame(data=rows,columns=['No.','img'])
  
    rows=list()
    for n,row in compDF.iterrows():
        comp=row['component']
        compNR=row['No.']
        for item in comp.rawfile_dict.keys():
            realRT=comp.rawfile_dict[item]['realRT']
            scanNumber=comp.rawfile_dict[item]['scanNumber']
            MScore=comp.rawfile_dict[item]['MScore']
            rows.append([compNR,item,realRT,scanNumber,MScore])
    dataDF=pd.DataFrame(data=rows,columns=['No.','rawfile','observedRT','scanNumber','MScore'])
    
    return(imgDF,dataDF)
                

def generateXIC_DF(compDF):
    comp=compDF['component'][0]
    ls_rawfile=list()
    for item in comp.rawfile_dict.keys():
        loc=comp.rawfile_dict[item]['location']
        rawfile=ThermoRawReader(loc)
        ls_rawfile.append(rawfile)
    
    ls_img=list()
    for n,row1 in compDF.iterrows():
        comp=row1['component']
        componentNumber=row1['No.']
        RT=comp.RT
        start=RT-2
        end=RT+2
        mz=comp.mzIsotopologueHighestIntensity
        mz_start=mz-(mz*20e-6)
        mz_end=mz+(mz*20e-6)
        massrange="%f-%f"%(mz_start,mz_end)
        ls_rows=list()
        for rawfile in ls_rawfile:
            baseName=rawfile.baseName
            color=comp.rawfile_dict[baseName]['color']
            label=comp.rawfile_dict[baseName]['condition']
            XIC=rawfile.getXIC(massrange,start,end,False)
            row=[list(XIC['RT']),list(XIC['intensity']),color,label]
            ls_rows.append(row)

        XIC_df=pd.DataFrame(data=ls_rows, columns=['XIC_RT','XIC_intens','color','label'])
        Graph=graph()
        Graph.getXICFig(XIC_df)
        fig=Graph.fig
        fig_buffer=BytesIO()
        fig.savefig(fig_buffer)
        img=Image.open(fig_buffer)
        ls_img.append([componentNumber,img])
    
    for rawfile in ls_rawfile:
        rawfile.close()
    
    imgDF=pd.DataFrame(data=ls_img, columns=['No.','XIC_image'])
    return(imgDF)
    
def loadConfigFile(path):
    import xml.dom.minidom as xml
    doc = xml.parse(path)
    sections=doc.getElementsByTagName("section")
    rows=list()
    for section in sections:
        if section.getAttribute("name")=="project":  
            elements=section.getElementsByTagName("param")
            for element in elements:                        
                if element.getAttribute("name")=="project name":
                    projectName=element.getAttribute("value")
                    rows.append(['projectName',projectName])
                
                elif element.getAttribute("name")=="analyst":
                    analyst=element.getAttribute("value")
                    rows.append(['analyst',analyst])

                elif element.getAttribute("name")=="project directory":
                    projectDirectory=element.getAttribute("value")
                    rows.append(['projectDirectory',projectDirectory]) 

                elif element.getAttribute("name")=="project code":
                    projectCode=element.getAttribute("value")
                    rows.append(['projectCode',projectCode])   
                
                elif element.getAttribute("name")=="DLIMS project":
                    DLIMSCode=element.getAttribute("value")
                    rows.append(['DLIMSCode',DLIMSCode])  
        
        elif section.getAttribute("name")=="XIC":
            elements=section.getElementsByTagName("param")
            for element in elements:
                if element.getAttribute("name")=="based on":
                    XIC_basedOn=element.getAttribute("value")
                    rows.append(["XIC_basedOn",XIC_basedOn])
                
                elif element.getAttribute("name")=="time window":
                    XIC_timeWindow=element.getAttribute("value")
                    rows.append(['XIC_timeWindow',XIC_timeWindow])
                    
                elif element.getAttribute("name")=="apply smooth":
                    XIC_applySmooth=element.getAttribute("value")
                    rows.append(['XIC_applySmooth',XIC_applySmooth])


                elif element.getAttribute("name")=="apply type":
                    XIC_smoothType=element.getAttribute("value")
                    rows.append(['XIC_smoothType',XIC_smoothType])
                    
                elif element.getAttribute("name")=="smooth cycles":
                    XIC_cycles=element.getAttribute("value")
                    rows.append(['XIC_cycles',XIC_cycles])                     
                
                elif element.getAttribute("name")=="smooth window":
                    XIC_smomthWindow=element.getAttribute("value")
                    rows.append(['XIC_smomthWindow',XIC_smomthWindow])
            
        elif section.getAttribute("name")=="Full MS":  
            elements=section.getElementsByTagName("param")
            for element in elements:                        
                if element.getAttribute("name")=="PYQMSppm":
                    PYQMS_massErrorPPM=element.getAttribute("value")
                    rows.append(['PYQMS_massErrorPPM',PYQMS_massErrorPPM])

        elif section.getAttribute("name")=="MSMS":  
            elements=section.getElementsByTagName("param")
            for element in elements:                        
                if element.getAttribute("name")=="error tolerance":
                   MSMSerrorTolerance=element.getAttribute("value")
                   rows.append(['MSMSerrorTolerance',MSMSerrorTolerance])
                
                elif element.getAttribute("name")=="error unit":
                    MSMSerrorUnit=element.getAttribute("value")
                    rows.append(['MSMSerrorUnit',MSMSerrorUnit])
                
                elif element.getAttribute("name")=="assign to":
                    MSMSassignTo=element.getAttribute("value")
                    rows.append(['MSMSassignTo',MSMSassignTo])

                
                elif element.getAttribute("name")=="a-ion":
                    MSMS_aIon=element.getAttribute("value")
                    rows.append(['MSMS_aIon',MSMS_aIon])


                elif element.getAttribute("name")=="b-ion":
                    MSMS_bIon=element.getAttribute("value")
                    rows.append(['MSMS_bIon',MSMS_bIon])

                elif element.getAttribute("name")=="c-ion":
                    MSMS_cIon=element.getAttribute("value")
                    rows.append(['MSMS_cIon',MSMS_cIon])       

                elif element.getAttribute("name")=="x-ion":
                    MSMS_xIon=element.getAttribute("value")
                    rows.append(['MSMS_xIon',MSMS_xIon])

                elif element.getAttribute("name")=="y-ion":
                    MSMS_yIon=element.getAttribute("value")
                    rows.append(['MSMS_yIon',MSMS_yIon])

                elif element.getAttribute("name")=="z-ion":
                    MSMS_zIon=element.getAttribute("value")
                    rows.append(['MSMS_zIon',MSMS_zIon])

                elif element.getAttribute("name")=="fragment losses":
                    MSMS_fragmentLosses=element.getAttribute("value")
                    rows.append(['MSMS_fragmentLosses',MSMS_fragmentLosses])
    
    output_table=pd.DataFrame(data=rows,columns=['Attribute','Value'])
    return(output_table)
    
