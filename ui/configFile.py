# -*- coding: utf-8 -*-
"""
Created on Sun Jun  9 10:40:47 2019

@author: KTMS458
"""
import xml.dom.minidom as xml

def createDefaultXML(path):
    import os
    fn="default_configuration.xml"
    fn=os.path.join(os.path.dirname(path),fn)
    fn=fn.replace('\\','/')
    file_handle = open(fn,"w")
    
    XML="""<?xml version="1.0" ?>
<configuration>
    <section 
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