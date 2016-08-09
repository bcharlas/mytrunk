# -*- coding: utf-8 -*-
"""
Created on Mon Oct 12 14:21:58 2015

@author: amminex
"""

from yade.pack import *
import numpy as np

import bodyHeat
import cellHeat


mass_scaling=1.



# Enthalpies & Entropies of reactions
#-------------------------------------
# 1-8
DH=41.4*1e6# mJ/mol
DH_surf=39.2*1e6# mJ/mol
DS=228.1*1e3 # mJ/mol/K
#
p0=1e-6# MPa
#
E_a=32.9*1e6 # mJ/mol
k_a=9.5e3 #mol/s
M_a=1. # check unit
N_a=1. # check unit
#
E_d=32.9*1e6 #mJ/mol
k_d=9.5e3 #mol/s
M_d=1. # check unit
N_d=1. # check unit

M_0=158.*1e-6*mass_scaling #t/mol
M_NH3=17.05*1e-6*mass_scaling #t/mol


# Volumetric mass 
Rho_0=2.6*1e-9*mass_scaling #kg/m3
Rho_1=1.9*1e-9*mass_scaling #kg/m3
Rho_8=1.3*1e-9*mass_scaling #kg/m3

R=8.314*1e3 #mJ/mol/K

# Thermal
Cp=75.55*1e3 #mJ/mol/K
k_1= 0.1 #mJ/mm/K/s
k_8= 0.48 #mJ/mm/K/s
DL=0.025# fictive distance between salt and wall in m

s_max=8

# Initial parameters
#---------------
'''
s_ini=0.# initial saturation degree 0 < s_0 <1
s=s_ini


s_max=8.
#
n_s_0= M_Ad_0/M_0*8  # initial number of moles of NH3 which can be absorbed in the salt
#
V_Ad_0= M_Ad_0 / Rho_0#cm3
#
M_mol=M_0 + s_0 * M_NH3 * 8# molar mass of salt, kg/mol
'''
# Simulation props
# ------------------
dt=O.dt

sorpBodies=[]
#
incBodySorp=[]
#
bodySat=[]
#
bodyMass0=[] # in moles
bodyAbsMax=[] # in moles
#
v0=[]

def initSorpPropsList(bodyList, s_ini, s_max, Cp):
  
  global sorpBodies
  global incBodySorp
  global bodyMass0
  global bodySat
  global V0
  
  sorpBodies=[]
  #
  incBodySorp = []
  #
  V0=[]

  #
  for b in bodyList:
    sorpBodies.append(b)
    #
    incBodySorp.append(0.) 
    #
    bodyMass0.append(O.bodies[b].state.mass/ (0.8594*s_ini + 1.)) # mass of the body without nh3 absorbed
    bodySat.append(s_ini)
    #
    bodyHeat.listCp[bodyHeat.heatBodies.index(b)] = Cp / (M_0 + s_ini * M_NH3)
    bodyHeat.listHeatCond[bodyHeat.heatBodies.index(b)] = (0.5343-0.0543*s_ini)
    #
    V0.append(4./3. * np.pi * np.power(O.bodies[b].shape.radius,3.) / (2.65*s_ini + 1.08))


# SIMULATION
#--------------
# WARNING no temperature of the cell yet!
def sorption(flow):
    global sorpBodies
    global incBodySorp
    global bodySat
    global bodyMass0
    global s_max
    
    # ABSORPTION / DESORPTION
    for i in range(flow.nCells()):
        Tc=cellHeat.cellTemp[i]
        pc=flow.getCellPressure(i)
        Vc=flow.volumeVoidPore(i)
        #
        if Vc>0:
            for b in flow.getVertices(i):
                if b in sorpBodies:  
                    s= bodySat[sorpBodies.index(b)] # saturation degree
                    #
                    Tb = bodyHeat.bodyTemp[bodyHeat.heatBodies.index(b)]
                    Rb = O.bodies[b].shape.radius
                    #
                    peq=np.exp(-DH_surf/R/Tb + DS/R) * p0
                    p_rel=(pc-peq)/peq
                    #
                    #
                    r_a=M_NH3*k_a*np.exp(-E_a/(R*Tb))*(1-s)**M_a*p_rel**N_a # absorption reaction rate
                    r_d=M_NH3*k_d*np.exp(-E_d/(R*Tb))*(s)**M_d*p_rel**N_d # desorption reaction rate
                    #
                    r=r_a-r_d # overal reaction rate
                    #
                    if not( ((r>0) and (s>= s_max))   or   ((r<0) and (s <= 0))  ):
                        incBodySorp[sorpBodies.index(b)] += M_NH3*r*O.dt # MASS increment of the body
                        cellHeat.incP[i] -= r*O.dt * R * Tc / Vc # pressure increment of the cell
                        bodyHeat.incBodyTemp[bodyHeat.heatBodies.index(b)] +=DH_surf*r/M_NH3/bodyHeat.listCp[bodyHeat.heatBodies.index(b)] * O.dt # temperature increment of the body due to reaction exo/endo-thermy
    
    # DIFFUSION IN THE BODIES: transfer as soon as the saturation degree is higher
    for i in O.interactions:
        if i.id1 in sorpBodies and i.id2 in sorpBodies:
            #
            index1=sorpBodies.index(i.id1)
            index2=sorpBodies.index(i.id2)
            
            s1=bodySat[sorpBodies.index(i.id1)]
            s2=bodySat[sorpBodies.index(i.id2)]
            
            m01=bodyMass0[index1]
            m02=bodyMass0[index2]
            
            if s1 != s2:
                #
                r1=O.bodies[i.id1].shape.radius
                r2=O.bodies[i.id2].shape.radius
                #
                pene=i.geom.penetrationDepth
                if pene >0. :
                    rij = np.sqrt(0.5*(r1+r2)*pene-pene**2/4)# radius of the contact interface
                    #
                    T1=bodyHeat.bodyTemp[index1]
                    T2=bodyHeat.bodyTemp[index2]
                    
                    diffCoef1=M_NH3*k_a*np.exp(-E_a/(R*T1))# t / mm ARBITRARY NOW
                    diffCoef2=M_NH3*k_a*np.exp(-E_a/(R*T2))# t / mm ARBITRARY NOW
                    
                    # dSat = surface * diff coef/distance  * Delta_s * time
                    dSat =  rij**2 * np.pi / ((r1/diffCoef1+r2/diffCoef2)) * (s1-s2) * O.dt
                    #
                    # mass increment
                    incBodySorp[index1] += dSat*m02/(m01+m02)*0.8594 # scaling difference to the mass - 
                    incBodySorp[index2] -= dSat*m01/(m01+m02)*0.8594 # scaling difference to the mass - 
                    #
                    # temperature increments of the bodies
                    bodyHeat.incBodyTemp[bodyHeat.heatBodies.index(i.id1)] +=DH*dSat*m02/(m01+m02)*0.8594/M_NH3/bodyHeat.listCp[bodyHeat.heatBodies.index(i.id1)] * O.dt # temperature increment of the body due to reaction exo/endo-thermy
                    bodyHeat.incBodyTemp[bodyHeat.heatBodies.index(i.id2)] +=DH*dSat*m01/(m01+m02)*0.8594/M_NH3/bodyHeat.listCp[bodyHeat.heatBodies.index(i.id2)] * O.dt
                    
                

'''
def update(flow):
    # updating the sorption value, the temperature increment and co
    for b in bodyList:
'''







