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

massScaling = 1e12

mn,mx=Vector3(-2.5,-2.5,0.),Vector3(2.5,2.5,7.6935) # corners of the initial packing

O.materials.append(FrictMat(young=70e3,poisson=0.3,frictionAngle=fric_Wall,density=7.85*1e-9 * massScaling,label='walls'))
walls=aabbWalls([mn,mx],thickness=0.5,material='walls')
wallIds=O.bodies.append(walls)
#
for i in wallIds:
    O.bodies[i].state.blockedDOFs='xyzXYZ' #Blocked DOF

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
#bodyHeat.heatBodies.extend([4,5]) # adding upper and lower wall

bodyHeat.initHeatPropsList(bodyHeat.heatBodies, 500., 430.0*1e6 / massScaling, 0.480)

#print heat.incTemp
#
# Simulation parameters
#----------------------
#
dispSize=0.4
#
dtSafety=0.00000001# factor applied to get stable flow calculations too
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
        ,FlowEngine(dead=0,bndCondIsPressure=[0,0,0,0,0,0],bndCondValue=[0,0,0,0,0,0],
                           defTolerance=0.3, useSolver=3, label="flow")#Flow engine inserted "dead" does not work
        ,PyRunner(command='vtkRec()',iterPeriod=1000, dead=False, label='vtkRec')
	]

O.step()

# Initialize pressure

for i in range(flow.nCells()):
    flow.setCellPressure(i,0.1) #1 bar so 0.1 MPa


cellHeat.gasCp = 2.17*1e9 / massScaling # in mJ/t/K
cellHeat.gasMolarMass = 17.03*1e-6 * massScaling # in t/mol
cellHeat.gasHeatCond = 22.916*1e-3 # in mJ/s/mm/K or mW/mm/K
cellHeat.gasR = 8.314*1e3 #mJ/mol/K

cellHeat.initHeatPropsList(273.,flow)

O.step()

#O.step()

for i in range(flow.nCells()):
    print  flow.getCellPressure(i)

#cellHeat.recordVTK=1
'''
print flow.nCells()
print flow.getCellNeighbor(0)
'''
# Functions
#----------------------
#
def heatEx():
    # INITIALIZE increments
    #sorption.incBodySorp = [0.]*len(sorption.sorpBodies)
    #
    bodyHeat.incTemp = [0.]*len(bodyHeat.heatBodies)
    cellHeat.incTemp = [0.]*flow.nCells()
    cellHeat.incP = [0.]*flow.nCells()
    #
    # SORPTION of gas by the solid + diffusion in the solid
    #sorption.sorption(flow)
    #
    # TEMPERATURE EXCHANGE
    cellHeat.heatConductivity(wallIds,flow)
    bodyHeat.heatConductivity()# always putbodies bodyHeat after cellHeat
    #
    # UPDATE cell temp, cell Pressure, body temp, body mass, body volume
    for c in range(flow.nCells()):
        #flow.setCellPressure(c, flow.getCellPressure(i)+cellHeat.incP[c])
        cellHeat.cellTemp[c] += cellHeat.incTemp[c]
        
    for bIndex in range(len(bodyHeat.bodyTemp)):
        bodyHeat.bodyTemp[bIndex] += bodyHeat.incBodyTemp[bIndex]
        #
    '''
    for bIndex in range(len(sorption.sorpBodies)):
        b= sorption.sorpBodies[index]
        s_previous = sorption.bodySat[bIndex]
        O.bodies[b].state.mass += sorption.incBodySorp[bIndex]
        sorption.bodySat[bIndex] += sorption.incBodySorp[bIndex]/(0.8594*sorption.bodyMass0[bIndex])
        s= sorption.bodySat[bIndex]
        #
        mult= (2.65*s + 1.08)/(2.65*s_previous + 1.08)
        utils.growParticle(b, multiplier=mult, updateMass=False)
    '''
        

flow.meshUpdateInterval=100000


#bodyHeat.bodyTemp[-2]=150.
#bodyHeat.isoThermalList.append(4)
#

# saving results - missing temperature in the gas
vtkExporter = export.VTKExporter('vtkExporterTesting')
vtkExporter.exportSpheres(what=[('pos','b.state.pos')
                                    ,('radius','b.shape.radius')
                                    ,('temp','bodyHeat.bodyTemp[bodyHeat.heatBodies.index(b.id)]')
                                    ])
flow.saveVtk()
#
fileNr=0
with open('VTK/out_'+str(fileNr) +'.vtk', 'a') as fTemp:
    #fTemp.write('CELL_DATA 286')
    fTemp.write('SCALARS Temperature float 1\n')
    fTemp.write('LOOKUP_TABLE default\n')
    vtkID=0
    for i in range(len(cellHeat.cellTemp)):
        if flow.isDrawable(i):
            fTemp.write(str(float(cellHeat.cellTemp[i])) + '\n')
            vtkID +=1
    fTemp.write('\n')

def vtkRec():
    global fileNr
    vtkExporter.exportSpheres(what=[('pos','b.state.pos')
                                    ,('radius','b.shape.radius')
                                    ,('temp','bodyHeat.bodyTemp[bodyHeat.heatBodies.index(b.id)]')
                                    ])
    flow.saveVtk()
    #
    fileNr +=1
    with open('VTK/out_'+str(fileNr) +'.vtk', 'a') as fTemp:
        #fTemp.write('CELL_DATA 286') 
        fTemp.write('SCALARS Temperature float 1\n')
        fTemp.write('LOOKUP_TABLE default\n')
        vtkID=0
        for i in range(len(cellHeat.cellTemp)):
            if flow.isDrawable(i):
                fTemp.write(str(float(cellHeat.cellTemp[i])) + '\n')
                vtkID +=1
        fTemp.write('\n')
    

O.wait()
