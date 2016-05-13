#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  9 10:44:56 2015

@author: amminex
"""

from yade import plot, pack, qt, export, ymport#, vtk #,utils
from yade.pack import *
import numpy as np
#from yade import heat_ as heat

from yade import heat2 as heat

import os

O.materials.append(FrictMat(young=70e3,poisson=0.3,frictionAngle=0.3,density=1.9,label='material'))

sList=[]
r=1
ratio=0.85
N=10

#create a column of N spheres
b=ymport.text('pack')

sphBodies=O.bodies.append(b)

heat.heatBodies= sphBodies
heat.initHeatProps(O.bodies, 273., 0.01, 1000.)



for i in sphBodies:
    O.bodies[i].state.blockedDOFs='xyzXYZ'


O.engines=[
	ForceResetter(),
	InsertionSortCollider([Bo1_Sphere_Aabb()]),
	InteractionLoop(
		[Ig2_Sphere_Sphere_ScGeom6D()],
		[Ip2_FrictMat_FrictMat_FrictPhys()],
		[Law2_ScGeom_FrictPhys_CundallStrack()]
	)
                #
        # Newton integrator and scaled gravity
        ,NewtonIntegrator(damping=.2,gravity=(0,0,-9.81e-3),label='Nuton')
        ,ThermalFlowEngine(dead=0,bndCondIsPressure=[0,0,0,0,1,0],bndCondValue=[0,0,0,0,0,0],
                           defTolerance=0.3, useSolver=3, label="flow")
        #
        ,PyRunner(command='heatEx()',iterPeriod=1, dead=True, label='thermExc')
        ,PyRunner(command='vtkRec()',iterPeriod=2000, dead=True, label='vtkRec')
	]

flow.meshUpdateInterval=100000
print 't'
flow.initCellsProperties(273.,0.01,1000.)
print 'o'
flow.initBodiesProperties(sphBodies, 273., 0.01, 1000.)
print 't'
#O.step()
#thermExc.dead=False

#change initial temperature of one sphere
for i in range(0,len(b),50):
    flow.setBodyTemperature(i,500.)
print 'o'

def heatEx():
    bodiesTemp=[]
    for i in sphBodies:
        bodiesTemp.append(flow.getBodyTemperature(i))
        
    bodiesTemp=heat.heatConductivity(bodiesTemp)
    
    for i in sphBodies:
        flow.setBodyTemperature(i,bodiesTemp[i])
 


vtkExporter = export.VTKExporter('vtkExporterTesting')
def vtkRec():
    print 'what'
    vtkExporter.exportSpheres(what=[('pos','b.state.pos'),
                                    ('radius','b.shape.radius'),
                                    ('temp','flow.getBodyTemperature(b.id)')])



#O.run(10000000)
'''
vtkExporter.exportFacets(what=[('pos','b.state.pos')])
vtkExporter.exportInteractions(what=[('kn','i.phys.kn')])
vtkExporter.exportContactPoints(what=[('nn','i.geom.normal')])
vtkExporter.exportSomething(what=[('n','bodyTemp[b.id]')])

'''
#export([heat.bodyTemp[i] for i in b])

"""
# to use on flow
defTolerance(=0.05)¶

    Cumulated deformation threshold for which retriangulation of pore space is performed. If negative, the triangulation update will occure with a fixed frequency on the basis of FlowEngine::meshUpdateInterval


doInterpolate(=false)

    Force the interpolation of cell’s info while remeshing. By default, interpolation would be done only for compressible fluids. It can be forced with this flag.

"""
