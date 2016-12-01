#!/usr/bin/python
# -*- coding: utf-8 -*-


from yade import plot, pack, qt, export, ymport, bodyHeat, cellHeat, sorption, utils
from yade.pack import *
import numpy as np

import os

from multiprocessing import Process

#import heat

#import reaction2 as reaction

fric_Salt = 1.5
fric_Wall= 1.

massScaling =1.0
stiffScaling=1.0e6
#
dtFactorFlow=10
dtFactorHeat=100
dtFactorSorp=dtFactorHeat

ini_DL=10

mn,mx=Vector3(-2.5,-2.5,0.),Vector3(2.5,2.5,7.+ini_DL) # corners of the initial packing

O.materials.append(FrictMat(young=70e3/stiffScaling,poisson=0.3,frictionAngle=fric_Wall,density=7.85*1e-9 * massScaling,label='walls'))
walls=aabbWalls([mn,mx],thickness=1.,material='walls')
wallIds=O.bodies.append(walls)
#

for i in wallIds:
    O.bodies[i].state.blockedDOFs='xyzXYZ' #Blocked DOF

# SPHERES
E_Salt = 150.0*1e3/stiffScaling #GPa
nu_Salt = 0.3

rho_Salt = 3.0*1e-9 * massScaling# t/mm3
#
O.materials.append(CohFrictMat(young=E_Salt,poisson=nu_Salt,\
                               frictionAngle=fric_Salt,density=rho_Salt,\
                               momentRotationLaw=False,alphaKr=2., alphaKtw=2.,\
                               normalCohesion=200./stiffScaling, shearCohesion=200./stiffScaling,\
                               label='spheres'))


sp=pack.SpherePack()
sp.makeCloud(mn,mx,0.6,0.,100,seed=1) #"seed" make the "random" generation always the same
sphBodies=sp.toSimulation(material='spheres')
'''

b=ymport.text('base_geom_100')

sphBodies=O.bodies.append(b)
#print sphBodies
'''

#---
s_ini=0. # initial saturation degree
Cp=(476.58 -s_ini*27.5)*1e6 / massScaling #mJ/mol/K
#---
# Initialization of thermal properties of the bodies
bodyHeat.heatBodies= sphBodies # All spheres exchange temperature
#
bodyHeat.heatBodies.extend([0,1,2,3]) # adding lower wall
bodyHeat.initHeatPropsList(bodyHeat.heatBodies, 293., Cp, 0.480) # heat properties of the spheres
#
bodyHeat.isoThermalList.extend([0,1,2,3]) # lower wall is isothermal
# thermal properties of lower wall
bodyHeat.listCp[bodyHeat.heatBodies.index(0)]=900.0*1e6 / massScaling
bodyHeat.listHeatCond[bodyHeat.heatBodies.index(0)]=50.
bodyHeat.listCp[bodyHeat.heatBodies.index(1)]=900.0*1e6 / massScaling
bodyHeat.listHeatCond[bodyHeat.heatBodies.index(1)]=50.
bodyHeat.listCp[bodyHeat.heatBodies.index(1)]=900.0*1e6 / massScaling
bodyHeat.listHeatCond[bodyHeat.heatBodies.index(0)]=50.
bodyHeat.listCp[bodyHeat.heatBodies.index(1)]=900.0*1e6 / massScaling
bodyHeat.listHeatCond[bodyHeat.heatBodies.index(1)]=50.
#---
#Initialization of sorption properties
sorption.massScaling = massScaling# doesnot really work, it just updates the mass scaling but not other values : TO BE CORRECTED
sorption.initSorpPropsList(bodyHeat.heatBodies[:-4], s_ini, 1., Cp)# only spherical bodies are included in the sorption
#---


# Simulation parameters
#----------------------
#
dispSize=0.4
#
#dtSafety=0.0000000000001# factor applied to get stable flow calculations too
#dtSafety=0.0000000000001
dtSafety=0.1
#
# Colors
sphereColor=(0,0,1)
brokenSphereColor=(0,0,1)
pinColor=(1,1,1)
#

O.dt=dtSafety*PWaveTimeStep()#min(dtSafety*PWaveTimeStep(),1.0e-6)

multG=100

O.engines=[
    	ForceResetter(),
    	InsertionSortCollider([Bo1_Sphere_Aabb(),Bo1_Box_Aabb(),Bo1_Facet_Aabb()]),
    	InteractionLoop(
		[Ig2_Sphere_Sphere_ScGeom6D(),Ig2_Facet_Sphere_ScGeom(),Ig2_Box_Sphere_ScGeom()],
		[Ip2_FrictMat_FrictMat_FrictPhys(),
                 Ip2_CohFrictMat_CohFrictMat_CohFrictPhys(setCohesionOnNewContacts=True,setCohesionNow=True, label='cohPhysProps')],
		[Law2_ScGeom6D_CohFrictPhys_CohesionMoment(),
                 Law2_ScGeom_FrictPhys_CundallStrack()]
	)
        ,NewtonIntegrator(damping=0.2,gravity=(0,0,-9.81e3*multG),label='Nuton')
        #
        ,PyRunner(command='sorpAndHeat()',iterPeriod=dtFactorHeat, dead=True, label='thermExc')
        #
        ,TranslationEngine(ids=[5],translationAxis=(0,0,1),velocity=-ini_DL/(O.dt*50000), dead=False, label='trans')
        ,FlowEngine(dead=1,useSolver=3, label="flow")#Flow engine inserted "dead" does not work
        ,PyRunner(command='vtkRec()',iterPeriod=dtFactorHeat*10, dead=True, label='vtkRecord')
        ,PyRunner(command='addPlotData()',iterPeriod=dtFactorHeat, dead=True, label='plotTemp')
	]

O.run(50000)
O.wait()

trans.dead=True
O.bodies[5].state.blockedDOFs='xyzXYZ'
O.bodies[5].state.vel=Vector3(0.,0.,0.)
Nuton.gravity=(0,0,-9.81e3)
O.wait()
#

#
# FLOW Parameters
flow.dead=0

flow.isCompressible=True
flow.fluidBulkModulus=1.#1.0 #case d
flow.defTolerance=0.3
flow.meshUpdateInterval=10
flow.useSolver=3
flow.permeabilityFactor=1
flow.viscosity=9.19*1e-5 / 1e6 # viscosity of NH3 from Pa.s to MPa.s
flow.bndCondIsPressure=[0,0,0,0,1,0]
flow.bndCondValue=[0,0,0,0,0.01,0]
flow.boundaryUseMaxMin=[1,1,1,1,1,1]
flow.dt=O.dt*dtFactorFlow

O.step()
cohPhysProps.setCohesionOnNewContacts=False
O.wait()




cellHeat.gasCp = 2.17*1e9 / massScaling # in mJ/t/K
cellHeat.gasMolarMass = 17.03*1e-6 * massScaling # in t/mol
cellHeat.gasHeatCond = 22.916*1e-3 # in mJ/s/mm/K or mW/mm/K
cellHeat.gasR = 8.314*1e3 #mJ/mol/K

cellHeat.initHeatPropsList(311.,flow)
#


O.wait()
# Functions
#----------------------
#
def sorpEx():
    #
    #initialization
    flow.clearImposedFlux()
    sorption.incBodyMass = [0.]*len(sorption.sorpBodies)
    sorption.incBodyTemp = [0.]*len(bodyHeat.heatBodies)
    sorption.incP = [0.]*flow.nCells()
    #
    sorption.sorption(flow, dtFactorSorp)

def heatEx():
    # initialization
    bodyHeat.incBodyTemp = [0.]*len(bodyHeat.heatBodies)
    cellHeat.incTemp = [0.]*flow.nCells()
    cellHeat.incP = [0.]*flow.nCells()
    
    #cellHeat.heatConductivity(wallIds,flow)
    bodyHeat.heatConductivity(dtFactorHeat)# always putbodies bodyHeat after cellHeat
    #
    
def updateCells():
    print sorption.incP
    for c in range(flow.nCells()-1):
        #print flow.getCellPressure(c)
        #print cellHeat.incP[c]
        #print sorption.incP[c]
        flow.imposeFluxFromId(c, - (cellHeat.incP[c] +sorption.incP[c])/flow.getCellPressure(c) * flow.volumeVoidPore(c))
        #flow.setCellPressure(c, flow.getCellPressure(c)+cellHeat.incP[c])
        #cellHeat.cellTemp[c] += cellHeat.incTemp[c]

def updateHeatBodies():
    for bIndex in range(len(bodyHeat.bodyTemp)):
        if bodyHeat.heatBodies[bIndex] not in bodyHeat.isoThermalList:
            bodyHeat.bodyTemp[bIndex] += bodyHeat.incBodyTemp[bIndex]
            if bodyHeat.heatBodies[bIndex] in sorption.sorpBodies:
                bodyHeat.bodyTemp[bIndex] += sorption.incBodyTemp[bIndex]
        #
def updateSorpBodies():
    for bIndex in range(len(sorption.sorpBodies)):
        if sorption.incBodyMass[bIndex] != 0:# if there is a mass increment
            b= sorption.sorpBodies[bIndex]
            s_previous = sorption.bodySat[bIndex]
            O.bodies[b].state.mass += sorption.incBodyMass[bIndex] # mass increment
            #print(sorption.incBodyMass[bIndex])
            sorption.bodySat[bIndex] = (O.bodies[b].state.mass - sorption.bodyMass0[bIndex]) /(sorption.k_m) # increment of saturation degree
            s= sorption.bodySat[bIndex]
            #
            mult=(sorption.k_v*s + 1.)/(sorption.k_v*s_previous + 1.) # Volume multiplier
            utils.growParticle(b, multiplier=mult, updateMass=False)
            #
            bodyHeat.listCp[bodyHeat.heatBodies.index(b)] = (476.58 -s*27.5)*1e6 / massScaling # updating the Cp as a function of the saturation degree

def runInParallel(*fns):
  proc = []
  for fn in fns:
    p = Process(target=fn)
    p.start()
    proc.append(p)
  for p in proc:
    p.join()
    
def sorpAndHeat():
    runInParallel(sorpEx, heatEx)
    #O.wait()
    runInParallel(updateCells,updateHeatBodies,updateSorpBodies)
    #O.wait()

# saving results - missing temperature in the gas
vtkExporter = export.VTKExporter('vtkExporterTesting')

vtkExporter.exportSpheres(what=[('pos','b.state.pos')
                                    ,('radius','b.shape.radius')
                                    ,('mass','b.state.mass')
                                    ,('temp','bodyHeat.bodyTemp[bodyHeat.heatBodies.index(b.id)]')
                                    ,('sat','sorption.bodySat[sorption.sorpBodies.index(b.id)]')
                                    ])
flow.saveVtk()
#


def vtkRec():
    global fileNr
    
    vtkExporter.exportSpheres(what=[('pos','b.state.pos')
                                    ,('radius','b.shape.radius')
                                    ,('mass','b.state.mass')
                                    ,('temp','bodyHeat.bodyTemp[bodyHeat.heatBodies.index(b.id)]')
                                    ,('sat','sorption.bodySat[sorption.sorpBodies.index(b.id)]')
                                    ])
                      
    flow.saveVtk()

# ploting temp and sat degree
def addPlotData():
    #
    plot.addData(temp6=bodyHeat.bodyTemp[bodyHeat.heatBodies.index(0+6)],
                 temp55=bodyHeat.bodyTemp[bodyHeat.heatBodies.index(49+6)],
                 temp105=bodyHeat.bodyTemp[bodyHeat.heatBodies.index(100)],
                 tempSide=bodyHeat.bodyTemp[bodyHeat.heatBodies.index(0)],
                 t=O.time,
                 sat6=sorption.bodySat[sorption.sorpBodies.index(6)]
                 ,sat55=sorption.bodySat[sorption.sorpBodies.index(55)]
                 ,sat105=sorption.bodySat[sorption.sorpBodies.index(100)]
                 ,avg_P=flow.averagePressure()
                 )
    
# SIMULATION
O.wait()

#

plot.plots={'t':('temp6','temp55','temp105','tempSide'),
            't ':('sat6','sat55','sat105')
            ,'t  ':('avg_P')}

#plot.plots={'t':('tempCell1',)}
plot.plot()

plot.addData(temp6=bodyHeat.bodyTemp[bodyHeat.heatBodies.index(0+6)],
             temp55=bodyHeat.bodyTemp[bodyHeat.heatBodies.index(49+6)],
             temp105=bodyHeat.bodyTemp[bodyHeat.heatBodies.index(100)],
             tempSide=bodyHeat.bodyTemp[bodyHeat.heatBodies.index(0)],
             t=O.time,
             sat6=sorption.bodySat[sorption.sorpBodies.index(6)]
             ,sat55=sorption.bodySat[sorption.sorpBodies.index(55)]
             ,sat105=sorption.bodySat[sorption.sorpBodies.index(100)]
             ,avg_P=flow.averagePressure()
             )
#vtkRecord.dead=False

plotTemp.dead=False
thermExc.dead=False
vtkRecord.dead=False

O.run(100000)


'''
# cycling
for i in range(5):
    print('cycle ' + str(i)) 
    # desorption
    flow.bndCondIsPressure=[0,0,0,0,1,0]
    flow.bndCondValue=[0,0,0,0,1e-3,0]
    flow.boundaryUseMaxMin=[1,1,1,1,1,1]
    #
    cohPhysProps.setCohesionOnNewContacts=False
    #
    bodyHeat.bodyTemp[bodyHeat.heatBodies.index(0)]=293.
    bodyHeat.bodyTemp[bodyHeat.heatBodies.index(1)]=293.
    bodyHeat.bodyTemp[bodyHeat.heatBodies.index(2)]=293.
    bodyHeat.bodyTemp[bodyHeat.heatBodies.index(3)]=293.
    O.wait()
    O.run(50000)
    O.wait()
    #
    #Absorption
    flow.bndCondIsPressure=[0,0,0,0,1,0]
    flow.bndCondValue=[0,0,0,0,1.4,0]
    flow.boundaryUseMaxMin=[1,1,1,1,1,1]
    #
    cohPhysProps.setCohesionOnNewContacts=True
    #
    bodyHeat.bodyTemp[bodyHeat.heatBodies.index(0)]=373.
    bodyHeat.bodyTemp[bodyHeat.heatBodies.index(1)]=373.
    bodyHeat.bodyTemp[bodyHeat.heatBodies.index(2)]=373.
    bodyHeat.bodyTemp[bodyHeat.heatBodies.index(3)]=373.
    O.wait()
    #
    O.run(50000)
    O.wait()
'''

