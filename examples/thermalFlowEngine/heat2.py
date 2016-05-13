# -*- coding: utf-8 -*-
"""
Created on Wed Jan 13 10:55:38 2016

@author: amminex
"""


from yade.pack import *
import numpy as np

#global heatBodies
heatBodies=[] # bodies concerned by the heat exchange
#
incTemp=[] # temperature increment

#global bodyTemp
bodyTemp=[] # temperature of the body
#
listCp=[] # heat capacity in 
listHeatCond=[] # heat conductivity

def initHeatProps(bodies, iniTemp, Cp, heatCond, flowThermalExchange=False):
  
  global bodyTemp
  global heatBodies
  global listCp
  global listHeatCond
  
  heatBodies=[]
  #
  bodyTemp = []
  incTemp =[]
  #
  listCp =[]
  listHeatCond=[]
  
  #
  for b in bodies:
    heatBodies.append(b.id)
    #
    bodyTemp.append(iniTemp) 
    incTemp.append(0) 
    #
    listCp.append(Cp)
    listHeatCond.append(heatCond)

def initHeatPropsList(bodyList, iniTemp, Cp, heatCond, flowThermalExchange=False):
  
  global bodyTemp
  global heatBodies
  global listCp
  global listHeatCond
  
  heatBodies=[]
  #
  bodyTemp = []
  incTemp =[]
  #
  listCp =[]
  listHeatCond=[]
  
  #
  for b in bodyList:
    heatBodies.append(b)
    #
    bodyTemp.append(iniTemp) 
    incTemp.append(0) 
    #
    listCp.append(Cp)
    listHeatCond.append(heatCond)

def heatConductivity(bodyTemp):
  
  #global bodyTemp
  global heatBodies
  global listCp
  global listHeatCond
  global incTemp
  
  #
  if len(heatBodies) >1:
    #
    incTemp = [0.]*len(heatBodies) # reinitialize the temperature increment
    #
    for i in O.interactions:
      if i.isActive:
	
	if i.id1 in heatBodies and i.id2 in heatBodies:# check if the bodies are registered for heat exchange
	  #
	  index1=heatBodies.index(i.id1)
	  index2=heatBodies.index(i.id2)
	  #
	  try:
	    #check if body 1 is a sphere
	    r1=O.bodies[i.id1].shape.radius
	  except :
	    # if not (e.g. wall or facet), the radius for the conductivity is 0
	    r1=0.0
	  #
	  try:
	    #check if body 2 is a sphere
	    r2=O.bodies[i.id2].shape.radius
	  except :
	    r2=0.0
	  
	  #
	  intVect = O.bodies[i.id1].state.pos-O.bodies[i.id2].state.pos
	  pene= r1+r2 - np.sqrt(intVect[0]**2 + intVect[1]**2 + intVect[2]**2)
	  #
	  #
	  if r1==0.0:
	    surf=( 2*r2*pene - pene**2) * np.pi
	  elif r2==0.0:
	     surf=( 2*r1*pene - pene**2) * np.pi
	  else:
	     surf=( 2*min(r1,r2)*pene - pene**2) * np.pi
	  #
	  dQ= -surf * listHeatCond[index1]*listHeatCond[index2]/((r1 -pene/2.)*listHeatCond[index2] + listHeatCond[index1]*(r2 -pene/2.)) * (bodyTemp[index1]-bodyTemp[index2])
	  	  #
	  incTemp[index1]+=O.dt *dQ/(O.bodies[i.id1].state.mass*listCp[index1])
	  incTemp[index2]-=O.dt *dQ/(O.bodies[i.id2].state.mass*listCp[index2])
    #if flowThermalExchange==True:
      
    #
    #print range(len(bodyTemp))
    for bIndex in range(len(bodyTemp)):
      bodyTemp[bIndex] += incTemp[bIndex]
    #bodyTemp = [sum(x) for x in zip(bodyTemp, incTemp)]
    return bodyTemp
