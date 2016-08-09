#!/usr/bin/python
# -*- coding: utf-8 -*-


from yade import plot, pack, qt, export, ymport, heat#, heat,reaction#, vtk #,utils
from yade.pack import *
import numpy as np

import os

#from yade import heat_test as heat

#import heat

#import reaction2 as reaction

fric_Salt = 1.
fric_Wall= 1.

massScaling = 1e12

mn,mx=Vector3(-2.5,-2.5,0.),Vector3(2.5,2.5,7.6935) # corners of the initial packing

O.materials.append(FrictMat(young=70e3,poisson=0.3,frictionAngle=fric_Wall,density=7.85*1e-9*massScaling,label='walls'))
walls=aabbWalls([mn,mx],thickness=0.5,material='walls')
wallIds=O.bodies.append(walls)
#
for i in wallIds:
    O.bodies[i].state.blockedDOFs='xyzXYZ' #Blocked DOF

# SPHERES
E_Salt = 80.0*1e3 #GPa
nu_Salt = 0.3

rho_Salt = 3.0*1e-9*massScaling# g/cm3
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

for i in sphBodies:
    O.bodies[i].state.blockedDOFs='xyzXYZ' #Blocked DOF

#print sphBodies

heat.heatBodies= sphBodies
heat.heatBodies.extend([4]) # adding upper and lower wall

heat.initHeatPropsList(heat.heatBodies, 273., 430.0*1e6 / massScaling, 0.480)# units mm, N and consequences: mJ
#heat.initHeatPropsList(heat.heatBodies, 273., 897.0*1e6, 237.0)# units mm, N and consequences: mJ

#print heat.incTemp
#
# Simulation parameters
#----------------------
#
dispSize=0.4
#
dtSafety=0.8
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
                 Ip2_CohFrictMat_CohFrictMat_CohFrictPhys(setCohesionOnNewContacts=False,setCohesionNow=False, label='cohPhysProps')],
		[Law2_ScGeom6D_CohFrictPhys_CohesionMoment(),
                 Law2_ScGeom_FrictPhys_CundallStrack()]
	)
                #
        ,NewtonIntegrator(damping=.2,gravity=(0,0,-9.81*1e3 / massScaling))#-9.81e-3),label='Nuton')
        #
        ,PyRunner(command='heatEx()',iterPeriod=1, dead=False, label='thermExc')
        #
        #,FlowEngine(dead=0,bndCondIsPressure=[0,0,0,0,1,0],bndCondValue=[0,0,0,0,10.,0],
        #                   defTolerance=0.3, useSolver=3, label="flow")#Flow engine inserted "dead"
        ,PyRunner(command='vtkRec()',iterPeriod=10, dead=False, label='vtkRec')
        ,PyRunner(command='addPlotData()',iterPeriod=1, dead=False, label='plotTemp')
	]
heat.bodyTemp[heat.heatBodies.index(4)]=500.
heat.isoThermalList.append(4)

heat.surfaceRatio=1.
# Functions
#----------------------
#
def heatEx():
    heat.incBodyTemp = [0.]*len(heat.heatBodies)
    heat.heatConductivity()
    

vtkExporter = export.VTKExporter('vtkExporterTesting')
vtkExporter.exportSpheres(what=[('pos','b.state.pos'),
                                    ('radius','b.shape.radius'),
                                    ('temp','heat.bodyTemp[heat.heatBodies.index(b.id)]')])
def vtkRec():
    vtkExporter.exportSpheres(what=[('pos','b.state.pos'),
                                    ('radius','b.shape.radius'),
                                    ('temp','heat.bodyTemp[heat.heatBodies.index(b.id)]')])

O.wait()


def addPlotData():
    #
    plot.addData(temp0=heat.bodyTemp[heat.heatBodies.index(0+6)],
                 temp24=heat.bodyTemp[heat.heatBodies.index(24+6)],
                 temp49=heat.bodyTemp[heat.heatBodies.index(49+6)],
                 temp74=heat.bodyTemp[heat.heatBodies.index(74+6)],
                 temp99=heat.bodyTemp[heat.heatBodies.index(99+6)],
                 tempPlate4=heat.bodyTemp[heat.heatBodies.index(4)],
                 #tempPlate5=heat.bodyTemp[heat.heatBodies.index(5)],
                 t=O.iter)




plot.plots={'t':('temp0','temp24','temp49','temp74','temp99','tempPlate4','tempPlate5',)}
#plot.plots={'t':('tempCell1',)}
plot.plot()
plot.addData(temp0=heat.bodyTemp[heat.heatBodies.index(0+6)],
             temp24=heat.bodyTemp[heat.heatBodies.index(24+6)],
             temp49=heat.bodyTemp[heat.heatBodies.index(49+6)],
             temp74=heat.bodyTemp[heat.heatBodies.index(74+6)],
             temp99=heat.bodyTemp[heat.heatBodies.index(99+6)],
             tempPlate4=heat.bodyTemp[heat.heatBodies.index(4)],
             #tempPlate5=heat.bodyTemp[heat.heatBodies.index(5)],
             t=O.iter)