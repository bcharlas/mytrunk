# -*- coding: utf-8 -*-
"""
Created on Wed Jan 13 10:55:38 2016

@author: amminex
"""


from yade.pack import *
import numpy as np

import bodyHeat

gasCp=1.
gasHeatCond=1.
gasMolarMass=1.
gasR=8.31
#
g_pi=3.14
#
recordVTK=0

#heatCells=[] # bodies concerned by the heat exchange
#
incTemp=[] # temperature increment
incP=[] # pressure increment

#global cellTemp
cellTemp=[] # temperature of the ce3ll
cellMass=[] # mass of gas in the cell
#
#listCp=[] # heat capacity in 
#listHeatCond=[] # heat conductivity
#
isoThermalList=[]

'''
def initHeatProps(bodies, iniTemp, Cp, heatCond, flowThermalExchange=False):
  
  global cellTemp
  global heatCells
  global listCp
  global listHeatCond
  global incTemp
  
  #heatCells=[]
  #
  cellTemp = []
  incTemp =[]
  #
  listCp =[]
  listHeatCond=[]
  
  #
  for b in bodies:
    #heatCells.append(b.id)
    #
    cellTemp.append(iniTemp)

    incTemp.append(0.) 
    #
    listCp.append(Cp)
    listHeatCond.append(heatCond)
'''

def initHeatPropsList(iniTemp,flow):
  
  global cellTemp
  global cellMass
  #global heatCells
  global incTemp
  global incP
  global gasCp
  global gasHeatCond
  global gasMolarMass
  
  
  #heatCells=[]
  #
  cellTemp = []
  incTemp =[]
  incP=[]
  #
  #gasCp=Cp
  #gasHeatCond=heatCond
  #gasMolarMass=molarMass
  #
  #
  for c in range(flow.nCells()):
    #heatCells.append(b)
    #
    cellTemp.append(iniTemp) 
    m= flow.getCellPressure(c) * flow.volumeVoidPore(c) *gasMolarMass /(gasR * iniTemp)
    if m>0:
        cellMass.append(m)
    else:
        cellMass.append(0.)
    incTemp.append(0.)
    incP.append(0.)




def heatConductivity(exceptions, flow):
  
  global cellTemp
  #
  global gasCp
  global gasHeatCond
  global gasMolarMass
  #
  global incTemp
  global incP
  
  #
  if flow.nCells() >0:
    #
    #incTemp = [0.]*len(heatCells) # reinitialize the temperature increment
    #
    for i in range(flow.nCells()):
        Tc=cellTemp[i]
        #print Tc
        pc=flow.getCellPressure(i)
        Vc=flow.volumeVoidPore(i)
        #
        if Vc >0 and pc>0:
            cellMass[i] = pc * Vc *gasMolarMass /(gasR * Tc)
            #
            # Temperature increment due to the pressure
            # incTemp[i] += Tc -pc * Vc *gasMolarMass /(gasR * cellMass[i])
            #
            
            # Exchange with neighboring bodies
            for b in flow.getVertices(i):
                if b not in exceptions:                    
                    Tb = bodyHeat.bodyTemp[bodyHeat.heatBodies.index(b)]
                    Rb = O.bodies[b].shape.radius
                    #
                    cellRho= pc*gasMolarMass/(Tc*gasR)
                    #
                    dQ=(1.14+0.85/Rb)/O.bodies[b].coordNum()*4*g_pi*pow(Rb,2)*(Tc-Tb)*pow(abs(Tc-Tb),0.25)*O.dt
                    #
                    DT=dQ/cellRho/gasCp
                    incTemp[i]-=DT
                    #
                    bodyHeat.incTemp[bodyHeat.heatBodies.index(b)]+=dQ/bodyHeat.listCp[bodyHeat.heatBodies.index(b)] /O.bodies[b].state.mass # to be checked ...
                    #
                    incP[i] += -DT* cellMass[i] / gasMolarMass * gasR /Vc # inc of Pressure due to T change in cell
                    
            # Thermal exchange with neighboring cells: 
            for ngb in range(4):#flow.getCellNeighbor(i):
                # WARNING : I had to put this one because the algorithm found a ngb = flow.nCells() --> Why?
                Tn = cellTemp[flow.getCellNeighbor(i)[ngb]]
                DT = flow.kNorm(b,ngb) * gasHeatCond * (Tc-Tn) * O.dt # missing distance
                incTemp[i] -= DT
                #incTemp[ngb] += DT # WARNING: if this is "on" it will add 2 times the temperature increment to each cell !!!
                #
                incP[i] += -DT* cellMass[i] / gasMolarMass * gasR /Vc
                #incP[ngb] += DT* cellMass[i] / gasMolarMass * gasR /Vc # WARNING: if this is "on" it will add 2 times the pressure increment to each cell !!!
            #
            # MISSING: heat exchange from mass exchange: dm / Cp
            #
    '''
    # Pressure variation due to temperature variation (the pressure variation due to volume variation is already taken into account)
    for i in range(flow.nCells):
        pc=flow.getCellPressure(i)
        flow.setCellPressure(i) = pc*(1 + incTemp[i]/cellTemp[i])
    '''
    #
    '''
    # Incrementing the temperature of cells - here only cells updated
    for i in range(flow.nCells):
        cellTemp[i]+=incTemp[i]
        flow.setCellPressure(i, flow.getCellPressure(i)*(1-incTemp[i]/cellTemp[i]))
        incTemp[i]=0.
    '''
	
