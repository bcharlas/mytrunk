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
incBodyTemp=[] # temperature increment

#global bodyTemp
bodyTemp=[] # temperature of the body
cellTemp =[]
#
listCp=[] # heat capacity in 
listHeatCond=[] # heat conductivity
#
isoThermalList=[]
#
g_pi=3.14

surfaceRatio=1.

def initHeatProps(bodies, iniTemp, Cp, heatCond, flowThermalExchange=False):
  
  global bodyTemp
  global cellTemp
  global heatBodies
  global listCp
  global listHeatCond
  global incBodyTemp
  
  #heatBodies=[]
  #
  bodyTemp = []
  cellTemp =[]
  incBodyTemp =[]
  #
  listCp =[]
  listHeatCond=[]
  
  #
  for b in bodies:
    #heatBodies.append(b.id)
    #
    bodyTemp.append(iniTemp)
    cellTemp.append(iniTemp) 
    incBodyTemp.append(0.) 
    #
    listCp.append(Cp)
    listHeatCond.append(heatCond)

def initHeatPropsList(bodyList, iniTemp, Cp, heatCond, flowThermalExchange=False):
  
  global bodyTemp
  global heatBodies
  global listCp
  global listHeatCond
  global incBodyTemp
  
  heatBodies=[]
  #
  bodyTemp = []
  incBodyTemp =[]
  #
  listCp =[]
  listHeatCond=[]
  
  #
  for b in bodyList:
    heatBodies.append(b)
    #
    bodyTemp.append(iniTemp) 
    incBodyTemp.append(0.) 
    #
    listCp.append(Cp)
    listHeatCond.append(heatCond)

def heatConductivity():
    global bodyTemp
    global heatBodies
    global listCp
    global listHeatCond
    global incBodyTemp
    #
    global surfaceRatio
    #
    if len(heatBodies) >1:
        #incBodyTemp = [0.]*len(heatBodies) # reinitialize the temperature increment
        #
        for i in O.interactions:
            #
            if i.id1 in heatBodies and i.id2 in heatBodies:
                # check if the bodies are registered for heat exchange
                #
                index1=heatBodies.index(i.id1)
                index2=heatBodies.index(i.id2)
                #
                if bodyTemp[index1] != bodyTemp[index2]: # no thermal exchange if they have the same temperature ...
                    #
                    try:
                        #check if body 1 is a sphere
                        r1=O.bodies[i.id1].shape.radius
                    except :
                        # if not (e.g. wall or facet), the radius for the conductivity is 0
                        r1=0.0
                    try:
                        #check if body 2 is a sphere
                        r2=O.bodies[i.id2].shape.radius
                    except :
                        # if not (e.g. wall or facet), the radius for the conductivity is 0
                        r2=0.0
                    #
                    #intVect = O.bodies[i.id1].state.pos-O.bodies[i.id2].state.pos
                    #pene= r1+r2 - np.sqrt(intVect[0]**2 + intVect[1]**2 + intVect[2]**2)
                    pene=i.geom.penetrationDepth
                    #
                    #print('r1'+ str(r1) + '  r2: ' + str(r2))
                    if  r1!=0. or r2!=0.:
                        #
                        if pene >0. :
                            #
                            k1=listHeatCond[index2]
                            k2=listHeatCond[index2]
                            #
                            T1=bodyTemp[index1]
                            T2=bodyTemp[index2]
                            #
                            m1=O.bodies[i.id1].state.mass
                            m2=O.bodies[i.id2].state.mass
                            #
                            if r1==0.0:
                                dij=2*np.sqrt(r2*pene-pene**2/4.)
                                #dQ = -L2 * g_pi * r2 * pene * (T1-T2)/r2
                            elif r2== 0.:
                                dij=2*np.sqrt(r1*pene-pene**2/4.)
                                #dQ = -L1* g_pi * r1 * pene * (T1-T2)/r1
                            else:
                                dij = 2*np.sqrt(0.5*(r1+r2)*pene-pene**2/4)
                                #dQ = L1*L2 / ( (r1-pene/2.)*L2 + (r2-pene/2.)*L1) * g_pi * min(r1,r2) * pene * (T1-T2)
                            #
                            dQ = surfaceRatio* 1./4. * dij**2 * np.pi / ((r1/k1+r2/k2)) * (T1-T2) 
                            #
                            if (i.id1 not in isoThermalList) and (abs(O.dt *dQ/(m1*listCp[index1])) < abs(T1-T2)):
                                incBodyTemp[index1]-=O.dt *dQ/(O.bodies[i.id1].state.mass*listCp[index1])
                            if (i.id2 not in isoThermalList) and (abs(O.dt *dQ/(m2*listCp[index2])) < abs(T1-T2)):
                                incBodyTemp[index2]+=O.dt *dQ/(O.bodies[i.id2].state.mass*listCp[index2])
        '''    
        for bIndex in range(len(bodyTemp)):
            bodyTemp[bIndex] += incBodyTemp[bIndex]
        '''
