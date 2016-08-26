# -*- coding: utf-8 -*-
"""
Created on Mon Oct 12 14:21:58 2015

@author: amminex
"""

from yade.pack import *
import numpy as np

import bodyHeat
import cellHeat


mass_scaling=1.0e14



# Enthalpies & Entropies of reactions
#-------------------------------------
# 1-8
DH=41.4*1e6# mJ/mol
DH_surf=39.2*1e6# mJ/mol
DS=228.1*1e3 # mJ/mol/K
#
p0=1e-6# MPa
#
#reference parameters in which the model was calculated
refWeight=5.1*1e-6*mass_scaling# weight in t corresponding to paper "Surface adsorption in strontium chloride ammines"
refRad=70e-3/2#radius in mm of the average grain size in the reference mesurement
refVol=90.*1e3#volume of gas used in the reference mesurement
#
k_m=0.8594

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

s_max=1

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
    bodyMass0.append(O.bodies[b].state.mass/ (k_m*s_ini + 1.)) # mass of the body without nh3 absorbed
    bodySat.append(s_ini)
    #
    bodyHeat.listCp[bodyHeat.heatBodies.index(b)] = Cp #/ (M_0 + s_ini * M_NH3)
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
        if Vc>0. and Tc>0. and pc>0.:
            # HERE INSERT A LIST OF 4 "r" so that the r will be adjusted if it is too big
            r_List=[]
            dr=0.
            dP = 0
            #
            for b in flow.getVertices(i):
                if b in sorpBodies:
                    s= bodySat[sorpBodies.index(b)] # saturation degree
                    #
                    Tb = bodyHeat.bodyTemp[bodyHeat.heatBodies.index(b)]
                    #
                    peq=pc# to be checked
                    if Tb> DH_surf/DS:
                        Rb = O.bodies[b].shape.radius
                        #
                        Ksurf=-DH_surf/R/Tb + DS/R 
                        peq=(np.exp(Ksurf) * p0)
                        p_rel=(pc-peq)/peq
                        
                        if Ksurf<0:
                            print('error with Ksurf: ' + str(Ksurf))
                            break
                        #
                        ###
                        '''
                        mu = 30 # Hz
                        D_E_in=10.0*1e6# mJ/mol
                        D_E_tot=41.4*1e6#mJ/mol
                        D_E_surf=DH_surf
                        D_E_out=D_E_tot-D_E_surf+D_E_in
                        gamma=100.
                        #
                        kp=mu*np.exp(-D_E_in/(R*T))
                        theta = Ksurf*pc/(1+Ksurf*pc)#adsorption reaction from Langmuir adsorption
                        nS=bodyMass0[sorpBodies.index(b)](1- 0.8594*s)/M_0 * 8 # number of free sites in moles
                        '''
                        ###
                        #
                        r_a=k_a*np.exp(-E_a/(R*Tb))*(1-s)*p_rel # absorption reaction rate
                        r_d=k_d*np.exp(-E_d/(R*Tb))*(s)*p_rel # desorption reaction rate
                        #
                        #r=(r_a-r_d) * bodyMass0[sorpBodies.index(b)] /refWeight * pow(abs(refRad/Rb),0.25)/len(O.bodies[b].intrs()) * Vc/refVol * 10. #overal reaction rate

                        if p_rel >0.:
                            r= r_a 
                        else:
                            r= r_d
                        
                        r *= bodyMass0[sorpBodies.index(b)] /refWeight * pow(abs(refRad/Rb),0.25)/len(O.bodies[b].intrs()) * Vc/refVol * 100.
                        
                        if r>0 and (((incBodySorp[sorpBodies.index(b)]+M_NH3*r*O.dt)/(k_m*bodyMass0[sorpBodies.index(b)]) +s < s_max)) and ( ((incBodySorp[sorpBodies.index(b)]+M_NH3*r*O.dt)/k_m/bodyMass0[sorpBodies.index(b)] +s>0.)):
                            #print (incBodySorp[sorpBodies.index(b)]+M_NH3*r*O.dt)/(k_m*bodyMass0[sorpBodies.index(b)])
                            r_List.append((b,r))
                            dP -= r*O.dt * R * Tc / Vc # incP
                        

            
            ratioP = 1.
            if (pc + cellHeat.incP[i] + dP < peq) and dP != 0:# correction if pressure increment is too big
                ratioP = abs((pc + cellHeat.incP[i])/dP)
            #
            cellHeat.incP[i] += ratioP * dP
            #
            #print('rList' + str(r_List))
            for i in range(len(r_List)):
                if r_List[i][0] in sorpBodies:
                    incBodySorp[sorpBodies.index(r_List[i][0])] += ratioP *r_List[i][1] * M_NH3 * O.dt # MASS increment of the body
                    dQ = ratioP * DH_surf *r_List[i][1] * O.dt
                    bodyHeat.incBodyTemp[bodyHeat.heatBodies.index(r_List[i][0])] += dQ / (O.bodies[r_List[i][0]].state.mass * bodyHeat.listCp[bodyHeat.heatBodies.index(r_List[i][0])]) #inc Body Temperature
            
                    
'''
                        if ((incBodySorp[sorpBodies.index(b)]+M_NH3*r*O.dt)/0.8594/bodyMass0[sorpBodies.index(b)] + s < s_max):
                            
                            if (pc + cellHeat.incP[i] - r*O.dt * R * Tc / Vc /len(O.bodies[b].intrs())) > 0:# safety to avoid having negative pressure
                                incBodySorp[sorpBodies.index(b)] += M_NH3*r*O.dt # MASS increment of the body
                                # pressure increment of the cell
                                cellHeat.incP[i] -= r*O.dt * R * Tc / Vc /len(O.bodies[b].intrs())
                                # temperature increment of the body due to reaction exo/endo-thermy WARNING it should be scaled to the mass !!!
                                bodyHeat.incBodyTemp[bodyHeat.heatBodies.index(b)] +=DH_surf*r/M_NH3/bodyHeat.listCp[bodyHeat.heatBodies.index(b)] * O.dt 
                            else:
                                print("Pressure increment too big !!!  dP:" + str(r*O.dt * R * Tc / Vc))
                        else:
                            print('s :' + str(s) + '    (incBodySorp[sorpBodies.index(b)]+M_NH3*r*O.dt)/0.8594/bodyMass0[sorpBodies.index(b)] :' + str((incBodySorp[sorpBodies.index(b)]+M_NH3*r*O.dt)/0.8594/bodyMass0[sorpBodies.index(b)]))
                    elif r<0 and (s>0):
                        if ((incBodySorp[sorpBodies.index(b)]+M_NH3*r*O.dt)/0.8594/bodyMass0[sorpBodies.index(b)] + s > 0):
                            incBodySorp[sorpBodies.index(b)] += M_NH3*r*O.dt # MASS increment of the body
                            # pressure increment of the cell
                            cellHeat.incP[i] -= r*O.dt * R * Tc / Vc /len(O.bodies[b].intrs())
                            # temperature increment of the body due to reaction exo/endo-thermy WARNING it should be scaled to the mass !!!
                            bodyHeat.incBodyTemp[bodyHeat.heatBodies.index(b)] +=DH_surf*r/M_NH3/bodyHeat.listCp[bodyHeat.heatBodies.index(b)] * O.dt 
'''
    
    # DIFFUSION IN THE BODIES: transfer as soon as the saturation degree is higher
'''
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
                    # mass increment - checking if mass remains >0 
                    if not( ((dSat>0) and (s1>= s_max))   and   ((dSat<0) and (s1 <= 0))  ):
                        if not( ((dSat<0) and (s2>= s_max))   and   ((dSat>0) and (s2 <= 0))  ):
                            if ((incBodySorp[index1] + dSat*m02/(m01+m02)*0.8594)>0) or (abs(incBodySorp[index1] + dSat*m02/(m01+m02)*0.8594)<O.bodies[i.id1].state.mass):
                                if ((incBodySorp[index2] + dSat*m01/(m01+m02)*0.8594)>0) or (abs(incBodySorp[index2] + dSat*m01/(m01+m02)*0.8594)<O.bodies[i.id2].state.mass):
                                    incBodySorp[index1] += dSat*m02/(m01+m02)*0.8594 # scaling difference to the mass - 
                                    incBodySorp[index2] -= dSat*m01/(m01+m02)*0.8594 # scaling difference to the mass - 
                    #
                    # temperature increments of the bodies
                    bodyHeat.incBodyTemp[bodyHeat.heatBodies.index(i.id1)] +=DH*dSat*m02/(m01+m02)*0.8594/M_NH3/bodyHeat.listCp[bodyHeat.heatBodies.index(i.id1)] * O.dt # temperature increment of the body due to reaction exo/endo-thermy
                    bodyHeat.incBodyTemp[bodyHeat.heatBodies.index(i.id2)] +=DH*dSat*m01/(m01+m02)*0.8594/M_NH3/bodyHeat.listCp[bodyHeat.heatBodies.index(i.id2)] * O.dt
'''
                

'''
def update(flow):
    # updating the sorption value, the temperature increment and co
    for b in bodyList:
'''







