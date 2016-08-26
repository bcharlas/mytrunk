#!/usr/bin/python
# -*- coding: utf-8 -*-


from yade import plot, pack, qt, export, ymport, bodyHeat, cellHeat, sorption, utils
from yade.pack import *
import numpy as np

import os

#import heat

#import reaction2 as reaction

fric_Salt = 1.
fric_Wall= 1.

massScaling = 1.0e14

mn,mx=Vector3(-2.5,-2.5,0.),Vector3(2.5,2.5,7.6935) # corners of the initial packing

O.materials.append(FrictMat(young=70e3,poisson=0.3,frictionAngle=fric_Wall,density=7.85*1e-9 * massScaling,label='walls'))
walls=aabbWalls([mn,mx],thickness=0.5,material='walls')
wallIds=O.bodies.append(walls)
#
'''
for i in wallIds:
    O.bodies[i].state.blockedDOFs='xyzXYZ' #Blocked DOF
'''
# SPHERES
E_Salt = 80.0*1e3 #GPa
nu_Salt = 0.3

rho_Salt = 3.0*1e-9 * massScaling# t/mm3
#
O.materials.append(CohFrictMat(young=E_Salt,poisson=nu_Salt,\
                               frictionAngle=fric_Salt,density=rho_Salt,\
                               momentRotationLaw=False,alphaKr=2., alphaKtw=2.,\
                               normalCohesion=1000000., shearCohesion=100000000.,\
                               label='spheres'))

'''
sp=pack.SpherePack()
sp.makeCloud(mn,mx,-1,0.3333,num_spheres,False, 0.3,seed=1) #"seed" make the "random" generation always the same
sphBodies=sp.toSimulation(material='spheres')
'''

b=ymport.text('base_geom_100')

sphBodies=O.bodies.append(b)
#print sphBodies

bodyHeat.heatBodies= sphBodies
bodyHeat.heatBodies.extend([4,5]) # adding lower wall
bodyHeat.isoThermalList.extend([4,5])

bodyHeat.initHeatPropsList(bodyHeat.heatBodies, 303., 430.0*1e6 / massScaling, 0.480)

# changing properties of lower and upper walls
bodyHeat.listCp[-1] = 900.0*1e6 / massScaling
bodyHeat.listHeatCond[-1] = 220.
#
bodyHeat.listCp[-2] = 900.0*1e6 / massScaling
bodyHeat.listHeatCond[-2] = 220.
#print heat.incTemp
#
# Simulation parameters
#----------------------
#
dispSize=0.4
#
#dtSafety=0.0000000000001# factor applied to get stable flow calculations too
#dtSafety=0.0000000000001
dtSafety=0.9
#
# Colors
sphereColor=(0,0,1)
brokenSphereColor=(0,0,1)
pinColor=(1,1,1)
#
vWalls=0.01
startPorosity=0.52




O.dt=dtSafety*PWaveTimeStep()

O.engines=[
	ForceResetter(),
	InsertionSortCollider([Bo1_Sphere_Aabb(),Bo1_Box_Aabb(),Bo1_Facet_Aabb()]),
	InteractionLoop(
		[Ig2_Sphere_Sphere_ScGeom6D(),Ig2_Facet_Sphere_ScGeom(),Ig2_Box_Sphere_ScGeom()],
		[Ip2_FrictMat_FrictMat_FrictPhys(),
                 Ip2_CohFrictMat_CohFrictMat_CohFrictPhys(setCohesionOnNewContacts=False,setCohesionNow=True, label='cohPhysProps')],
		[Law2_ScGeom6D_CohFrictPhys_CohesionMoment(),
                 Law2_ScGeom_FrictPhys_CundallStrack()]
	)
                #
        ,NewtonIntegrator(damping=.2,gravity=(0,0,0))#-9.81e-3),label='Nuton')
        #
        ,PyRunner(command='heatEx()',iterPeriod=1, dead=False, label='thermExc')
        #
        ,FlowEngine(dead=0,bndCondIsPressure=[0,0,0,0,1,0],bndCondValue=[0,0,0,0,1.4,0],
                           defTolerance=0.3, useSolver=3, label="flow")#Flow engine inserted "dead" does not work
        ,PyRunner(command='vtkRec()',iterPeriod=1000, dead=False, label='vtkRec')
        ,PyRunner(command='addPlotData()',iterPeriod=100, dead=False, label='plotTemp')
	]

flow.meshUpdateInterval=1000000000
flow.permeabilityFactor=1
flow.viscosity=9.19*1e-5 / 1e6 # viscosity of NH3 from Pa.s to MPa.s
flow.boundaryUseMaxMin=[0,0,0,0,1,0]


O.step()

for i in range(flow.nCells()):
    flow.setCellPressure(i,1.4) #14 bar so 1.4 MPa

cellHeat.gasCp = 2.17*1e9 / massScaling # in mJ/t/K
cellHeat.gasMolarMass = 17.03*1e-6 * massScaling # in t/mol
cellHeat.gasHeatCond = 22.916*1e-3 # in mJ/s/mm/K or mW/mm/K
cellHeat.gasR = 8.314*1e3 #mJ/mol/K

cellHeat.initHeatPropsList(311.,flow)
#
sorption.massScaling = massScaling# doesnot really work, it just updates the mass scaling but not other values : TO BE CORRECTED
#
Cp=476.58*1e6 / massScaling #mJ/mol/K
sorption.initSorpPropsList(bodyHeat.heatBodies[:-2], 0., 1., Cp)# if one wall is included in the thermal exchange


# Functions
#----------------------
#
def heatEx():
    sorption.incBodySorp = [0.]*len(sorption.sorpBodies)
    bodyHeat.incBodyTemp = [0.]*len(bodyHeat.heatBodies)
    cellHeat.incTemp = [0.]*flow.nCells()
    cellHeat.incP = [0.]*flow.nCells()
    #
    sorption.sorption(flow)
    #cellHeat.heatConductivity(wallIds,flow)
    bodyHeat.heatConductivity()# always putbodies bodyHeat after cellHeat
    #
    #for c in range(flow.nCells()):
        #flow.setCellPressure(c, flow.getCellPressure(i)+cellHeat.incP[c])# WARNING: the setter does not actually work :(
        #cellHeat.cellTemp[c] += cellHeat.incTemp[c]
    for bIndex in range(len(bodyHeat.bodyTemp)):
        bodyHeat.bodyTemp[bIndex] += bodyHeat.incBodyTemp[bIndex]
        #
    
    for bIndex in range(len(sorption.sorpBodies)):
        if sorption.incBodySorp[bIndex] != 0:# if there is a mass increment
            b= sorption.sorpBodies[bIndex]
            s_previous = sorption.bodySat[bIndex]
            O.bodies[b].state.mass += sorption.incBodySorp[bIndex] # mass increment
            #print(sorption.incBodySorp[bIndex])
            sorption.bodySat[bIndex] += sorption.incBodySorp[bIndex]/(sorption.k_m*sorption.bodyMass0[bIndex]) # increment of saturation degree
            s= sorption.bodySat[bIndex]
            #
            mult= (1.*s + 1.)/(1.*s_previous + 1.) # Volume multiplier
            utils.growParticle(b, multiplier=mult, updateMass=False)
            #
            #bodyHeat.listCp[bodyHeat.heatBodies.index(bIndex)] = (476.58 -s*27.5)*1e6 / massScaling # updating the Cp as a function of the saturation degree
    

       



#bodyHeat.bodyTemp[-2]=150.

#

# saving results - missing temperature in the gas
vtkExporter = export.VTKExporter('vtkExporterTesting')
'''
vtkExporter.exportSpheres(what=[('mass','O.bodies[b].state.mass')
                                    ,('temp','bodyHeat.bodyTemp[bodyHeat.heatBodies.index(b.id)]')
                                    ,('sat','sorption.bodySat[sorption.sorpBodies.index(b.id)]')
                                    ])
'''
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
                 temp105=bodyHeat.bodyTemp[bodyHeat.heatBodies.index(99+6)],
                 tempPlateBot=bodyHeat.bodyTemp[bodyHeat.heatBodies.index(4)],
                 t=O.time,
                 sat6=sorption.bodySat[sorption.sorpBodies.index(6)]
                 ,sat55=sorption.bodySat[sorption.sorpBodies.index(55)]
                 ,sat105=sorption.bodySat[sorption.sorpBodies.index(105)]
                 )


plot.plots={'t':('temp6','temp55','temp105','tempPlateBot'),
            't ':('sat6','sat55','sat105')}
#plot.plots={'t':('tempCell1',)}
plot.plot()

plot.addData(temp6=bodyHeat.bodyTemp[bodyHeat.heatBodies.index(0+6)],
             temp55=bodyHeat.bodyTemp[bodyHeat.heatBodies.index(49+6)],
             temp105=bodyHeat.bodyTemp[bodyHeat.heatBodies.index(99+6)],
             tempPlateBot=bodyHeat.bodyTemp[bodyHeat.heatBodies.index(4)],
             t=O.time,
             sat6=sorption.bodySat[sorption.sorpBodies.index(6)]
             ,sat55=sorption.bodySat[sorption.sorpBodies.index(55)]
             ,sat105=sorption.bodySat[sorption.sorpBodies.index(105)]
             )

    
# SIMULATION
#O.run(720000000)# ~1 min

O.wait()
O.run(112500)
#O.run(720000000)# ~1 min
