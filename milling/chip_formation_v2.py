import sys
sys.path.append("C:\Program Files\Python36\Lib\site-packages")
from abaqus import *
from abaqusConstants import *
from caeModules import *
from driverUtils import executeOnCaeStartup
import math
import tempfile
import sys
import os
import argparse
from math import *
import numpy as np
from sympy.geometry import *

MAX_TIME = 0.001

def mm(value):
    return value * 1e-3

def m(value):
    return value

def cm(value):
    return value * 1e-2

def rad(value):
    return value

def deg(value):
    return (value * pi)/180

import inspect
filename = inspect.getframeinfo(inspect.currentframe()).filename
cd = os.path.dirname(os.path.abspath(filename))

executeOnCaeStartup()
Mdb()
model = mdb.models['Model-1']

class Ti6AlV:

    def __init__(self):
        self.name = 'Ti6AlV'
        model = mdb.models['Model-1']
        model.Material(name=self.name)
        model.materials[self.name].Conductivity(table=((7.2, ), ))
        model.materials[self.name].JohnsonCookDamageInitiation(table=((-0.09, 0.25, 0.5, 0.014, 3.87, 1560.0, 25.0, 0.0011), ))
        model.materials[self.name].Density(table=((4470.0, ), ))
        model.materials[self.name].Depvar(deleteVar=1, n=1)
        model.materials[self.name].Elastic(table=((115000000000.0, 0.32), ))
        model.materials[self.name].Expansion(table=((8.85e-06, ), ))
        model.materials[self.name].InelasticHeatFraction()
        model.materials[self.name].Plastic(hardening=JOHNSON_COOK, table=((1098000000.0, 1092000000.0, 0.93, 1.1, 1560.0, 25.0), ))
        model.materials[self.name].plastic.RateDependent(type=JOHNSON_COOK, table=((0.014, 0.011), ))
        model.materials[self.name].SpecificHeat(table=((560.0, ), ))


class Carbide:

    def __init__(self):
        self.name = "Carbide"
        model = mdb.models['Model-1']
        model.Material(name=self.name)
        model.Material(name=self.name)
        model.materials[self.name].Conductivity(table=((55.0, ), ))
        model.materials[self.name].Density(table=((15000.0, ), ))
        model.materials[self.name].Elastic(table=((640000000000.0, 0.22), ))
        model.materials[self.name].Expansion(table=((6e-06, ), ))
        model.materials[self.name].SpecificHeat(table=((240.0, ), ))


class Tool:

    def __init__(self, name, width, height, thickness):
        self.name = name
        self.height = height
        self.width = width
        self.thickness = thickness
        self.radius = mm(5)

    def create(self):
        sketch = model.ConstrainedSketch(name='__profile__', sheetSize=0.02)
        g, v, d, c = sketch.geometry, sketch.vertices, sketch.dimensions, sketch.constraints
        corner1 = (0, 0)
        corner2 = (0, self.height)
        corner3 = (self.width, 0)
        corner4 = (self.width, self.height)
        sketch.rectangle(point1=(0.0, 0.0), point2=corner4)
        sketch.FilletByRadius(radius=self.radius, curve1=g[5], nearPoint1=corner1, curve2=g[2], nearPoint2=corner1)
        sketch.FilletByRadius(radius=self.radius, curve1=g[2], nearPoint1=corner2, curve2=g[3], nearPoint2=corner2)
        sketch.FilletByRadius(radius=self.radius, curve1=g[3], nearPoint1=corner4, curve2=g[4], nearPoint2=corner4)
        sketch.FilletByRadius(radius=self.radius, curve1=g[4], nearPoint1=corner3, curve2=g[5], nearPoint2=corner3)
        p = model.Part(name=self.name, dimensionality=THREE_D, type=DEFORMABLE_BODY)
        p.BaseSolidExtrude(sketch=sketch, depth=self.thickness)
        self.part = model.parts[self.name]
        del model.sketches['__profile__']

    def set_material(self, material):
        model.HomogeneousSolidSection(name=material.name, material=material.name, thickness=None)
        region = self.part.Set(cells=self.part.cells, name='Material region')
        self.part.SectionAssignment(region=region, sectionName=material.name, offset=0.0, 
            offsetType=MIDDLE_SURFACE, offsetField='', thicknessAssignment=FROM_SECTION)        

    def mesh(self):
        elemType1 = mesh.ElemType(elemCode=C3D8RT, elemLibrary=EXPLICIT, kinematicSplit=AVERAGE_STRAIN, secondOrderAccuracy=OFF, 
                                    hourglassControl=DEFAULT, distortionControl=DEFAULT)
        elemType2 = mesh.ElemType(elemCode=C3D6T, elemLibrary=EXPLICIT)
        elemType3 = mesh.ElemType(elemCode=C3D4T, elemLibrary=EXPLICIT)
        pickedRegions =(self.part.cells, )
        self.part.setElementType(regions=pickedRegions, elemTypes=(elemType1, elemType2, elemType3))
        self.part.generateMesh()

        p = mdb.models['Model-1'].parts['Tool']
        p.generateMesh()



class Workpiece(object):
    
    def __init__(self, name, length, width, height, radius, base_thikness):
        self.name = name
        self.length = length
        self.width = width
        self.height = height
        self.base_thikness = base_thikness
        self.radius = radius

    def create(self):
        sketch = model.ConstrainedSketch(name='__profile__', sheetSize=0.02)
        g, v, d, c = sketch.geometry, sketch.vertices, sketch.dimensions, sketch.constraints
        corner1 = (0, 0)
        corner2 = (0, self.width)
        corner3 = (self.length, 0)
        corner4 = (self.length, self.width)
        sketch.rectangle(point1=(0.0, 0.0), point2=corner4)


        sketch.CircleByCenterPerimeter(center=corner3, point1=(self.length, self.radius))
        sketch.CoincidentConstraint(entity1=v[4], entity2=g[5], addUndoState=False)
        sketch.autoTrimCurve(curve1=g[4], point1=corner3)
        sketch.autoTrimCurve(curve1=g[5], point1=corner3)
        sketch.autoTrimCurve(curve1=g[6], point1=(self.length + 0.01, -0.01))
        sketch.autoTrimCurve(curve1=g[9], point1=(0.110796421766281, 0.045832198113203))
        
        p = model.Part(name=self.name, dimensionality=THREE_D, type=DEFORMABLE_BODY)
        p.BaseSolidExtrude(sketch=sketch, depth=self.height)
        self.part = model.parts[self.name]
        del model.sketches['__profile__']

        p = mdb.models['Model-1'].parts['Workpiece']
        f, e = p.faces, p.edges
        t = p.MakeSketchTransform(sketchPlane=f[6], sketchUpEdge=e[6], sketchPlaneSide=SIDE1, sketchOrientation=RIGHT, origin=(0.0, 0.0, 0.0))
        s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=0.228, gridSpacing=0.005, transform=t)
        g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints

        p.projectReferencesOntoSketch(sketch=s, filter=COPLANAR_EDGES)
        s.rectangle(point1=corner1, point2=(self.width, self.length))
        
        f, e = self.part.faces, self.part.edges
        p.SolidExtrude(sketchPlane=f[6], sketchUpEdge=e[6], sketchPlaneSide=SIDE1, sketchOrientation=RIGHT, sketch=s, depth=self.base_thikness, flipExtrudeDirection=OFF)
        s.unsetPrimaryObject()
        del mdb.models['Model-1'].sketches['__profile__']

    def set_material(self, material):
        model.HomogeneousSolidSection(name=material.name, material=material.name, thickness=None)
        region = self.part.Set(cells=self.part.cells, name='Material region')
        self.part.SectionAssignment(region=region, sectionName=material.name, offset=0.0, 
            offsetType=MIDDLE_SURFACE, offsetField='', thicknessAssignment=FROM_SECTION)


class Assembly(object):

    def __init__(self, workpiece, tool):
        self.workpiece = workpiece
        self.tool = tool
        self.a = model.rootAssembly

    def create(self):
        self.a.DatumCsysByDefault(CARTESIAN)
        self.a.Instance(name=workpiece.name, part=workpiece.part, dependent=ON)
        self.a.Instance(name=tool.name, part=tool.part, dependent=ON)       
        self.a.rotate(instanceList=('Tool', ), axisPoint=(0.0, 0.0, 0.0), axisDirection=(1, 0.0, 0.0), angle=90.0)
        self.a.translate(instanceList=('Tool', ), vector=((self.workpiece.length - self.workpiece.radius)/2.0, 0.0, 0.0))


if __name__ == "__main__":
    titan = Ti6AlV()
    carbide = Carbide()

    tool = Tool("Tool", mm(30), mm(50), mm(5))
    tool.create()
    tool.set_material(carbide)
    tool.mesh()

    workpiece = Workpiece("Workpiece", mm(70), mm(90), mm(40), mm(60), mm(15))
    workpiece.create()    
    workpiece.set_material(titan)

    assembly = Assembly(workpiece, tool)
    assembly.create()