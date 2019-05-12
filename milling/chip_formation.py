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


class Material:

    def __init__(self, name, density, young, poisson, A, B, n, d1, d2, d3, ref_strain_rate, disp_at_failure):
        self.name = name
        mdb.models['Model-1'].Material(name=name)
        self.material = mdb.models['Model-1'].materials[name]
        self.material.Density(table=((density, ), ))
        self.material.Elastic(table=((young, poisson), ))
        self.material.Plastic(hardening=JOHNSON_COOK, table=((A, B, n, 0.0, 0.0, 0.0), ))
        self.material.JohnsonCookDamageInitiation(table=((d1, d2, d3, 0.0, 0.0, 0.0, 0.0, ref_strain_rate), ))
        self.material.johnsonCookDamageInitiation.DamageEvolution(type=DISPLACEMENT, table=((disp_at_failure, ), ))


class Step:

    def __init__(self, name, previous='Initial'):
        self.name = name
        mdb.models['Model-1'].ExplicitDynamicsStep(name=name, previous=previous, description='', timePeriod=MAX_TIME)


class Interaction:

    def __init__(self, step):
        # create interaction propery
        model = mdb.models['Model-1']
        model.ContactProperty('IntProp-1')
        model.interactionProperties['IntProp-1'].TangentialBehavior(
            formulation=PENALTY, directionality=ISOTROPIC, slipRateDependency=OFF, 
            pressureDependency=OFF, temperatureDependency=OFF, dependencies=0, table=((
            0.308, ), ), shearStressLimit=None, maximumElasticSlip=FRACTION, 
            fraction=0.005, elasticSlipStiffness=None)
        model.interactionProperties['IntProp-1'].NormalBehavior(
            pressureOverclosure=HARD, allowSeparation=ON, 
            constraintEnforcementMethod=DEFAULT)
        model.interactionProperties['IntProp-1'].HeatGeneration(
            conversionFraction=1.0, slaveFraction=0.5)
        
        model.ContactExp(name='Interaction', createStepName=step.name)
        mdb.models['Model-1'].interactions['Interaction'].includedPairs.setValuesInStep(
            addPairs=(
                (ALLSTAR, mdb.models['Model-1'].rootAssembly.surfaces['Tool Surface']), 
                (ALLSTAR, mdb.models['Model-1'].rootAssembly.surfaces['Workpiece Surface']), 
                (mdb.models['Model-1'].rootAssembly.surfaces['Tool Surface'], mdb.models['Model-1'].rootAssembly.surfaces['Workpiece Surface'])), 
            stepName=step.name, useAllstar=OFF)
        model.interactions['Interaction'].contactPropertyAssignments.appendInStep(stepName=step.name, assignments=((GLOBAL, SELF, 'IntProp-1'), ))




class Tool:

    def __init__(self, name, scale):
        self.name = name
        self.step_file = cd + './step_models/milling_cutter.stp'
        self.diameter = mm(10)  # from file
        self.length = mm(100)
        self.cutter_length = mm(22)
        self.scale = scale
        self.part = None

    def create(self):
        step = mdb.openStep(self.step_file)
        mdb.models['Model-1'].PartFromGeometryFile(name=self.name, geometryFile=step, 
                combine=False, dimensionality=THREE_D, type=DEFORMABLE_BODY, scale=self.scale)
        self.part = mdb.models['Model-1'].parts[self.name]


    def set_material(self, material):
        model = mdb.models['Model-1']
        model.HomogeneousSolidSection(name=material.name, material=material.name, thickness=None)  
        region = self.part.Set(cells=self.part.cells, name='Material region')
        self.part.SectionAssignment(region=region, sectionName=material.name, offset=0.0, 
            offsetType=MIDDLE_SURFACE, offsetField='', thicknessAssignment=FROM_SECTION)

    def mesh(self):
        self.part.seedPart(size=0.0025, deviationFactor=0.1, minSizeFactor=0.1)
        c = self.part.cells
        pickedRegions = c.getSequenceFromMask(mask=('[#1 ]', ), )
        self.part.setMeshControls(regions=pickedRegions, elemShape=TET, technique=FREE)
        elemType1 = mesh.ElemType(elemCode=C3D8T, elemLibrary=EXPLICIT)
        elemType2 = mesh.ElemType(elemCode=C3D6T, elemLibrary=EXPLICIT)
        elemType3 = mesh.ElemType(elemCode=C3D4T, elemLibrary=EXPLICIT, 
            secondOrderAccuracy=OFF, distortionControl=DEFAULT)
        cells = c.getSequenceFromMask(mask=('[#1 ]', ), )
        pickedRegions =(cells, )
        self.part.setElementType(regions=pickedRegions, elemTypes=(elemType1, elemType2, elemType3))
        self.part.generateMesh()

    def __str__(self):
        return "Tool: name={}, scale={}".format(self.name, self.scale)

    def __repr__(self):
        return str(self)


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


class Workpiece:
    def __init__(self, name, length, w_height, w_width, b_height, b_width):
        self.name = name
        self.length = length
        self.w_height = w_height
        self.w_width = w_width
        self.b_height = b_height
        self.b_width = b_width
        self.t = mm(5)
        self.part = None
        self.c_point = (0.5 * self.b_width, self.b_height + 0.5 * self.w_height, 0.5 * self.length)
        self.c1_point = (0.25 * (self.b_width - self.w_width), 0.5 * self.b_height, 0.5 * self.length)
        self.c2_point = (0.5 * self.b_width, 0.5 * self.b_height, 0.5 * self.length)
        self.c3_point = (self.b_width - 0.25 * (self.b_width - self.w_width), 0.5 * self.b_height, 0.5 * self.length)

    def create(self):
        s = mdb.models['Model-1'].ConstrainedSketch(name=self.name + '_sketch', sheetSize=max(self.length, self.w_height, self.w_width, 
            self.b_height, self.b_width))
        g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints    
        s.sketchOptions.setValues(decimalPlaces=4)
        s.setPrimaryObject(option=STANDALONE)

        p0 = (0, 0)
        p1 = (self.b_width, 0)
        p2 = (self.b_width, self.b_height)
        p3 = (0.5 * self.b_width + 0.5 * self.w_width, self.b_height)
        p4 = (0.5 * self.b_width + 0.5 * self.w_width, self.b_height + self.w_height)
        p5 = (0.5 * self.b_width - 0.5 * self.w_width, self.b_height + self.w_height)
        p6 = (0.5 * self.b_width - 0.5 * self.w_width, self.b_height)
        p7 = (0, self.b_height)

        s.Line(point1=p0, point2=p1)
        s.Line(point1=p1, point2=p2)
        s.Line(point1=p2, point2=p3)
        s.Line(point1=p3, point2=p4)
        s.Line(point1=p4, point2=p5)
        s.Line(point1=p5, point2=p6)
        s.Line(point1=p6, point2=p7)
        s.Line(point1=p7, point2=p0)


        p = mdb.models['Model-1'].Part(name=self.name, dimensionality=THREE_D, type=DEFORMABLE_BODY)
        self.part = mdb.models['Model-1'].parts[self.name]
        self.part.BaseSolidExtrude(sketch=s, depth=self.length)
        s.unsetPrimaryObject()
        # create 
        c = self.part.cells
        pickedCells = c.getSequenceFromMask(mask=('[#1 ]', ), )
        x = 0.5*(self.b_width - self.w_width)
        workpiece.part.PartitionCellByPlaneThreePoints(cells=pickedCells, point1=(x, mm(0), mm(0)), point2=(x, mm(1), mm(0)), point3=(x, mm(0), mm(1)))
        x = 0.5*(self.b_width + self.w_width)

        c = self.part.cells
        pickedCells = c.getSequenceFromMask(mask=('[#2 ]', ), )  # TODO replace
        workpiece.part.PartitionCellByPlaneThreePoints(cells=pickedCells, point1=(x, mm(0), mm(0)), point2=(x, mm(1), mm(0)), point3=(x, mm(0), mm(1)))
        pickedCells = c.getSequenceFromMask(mask=('[#4 ]', ), )  # TODO replace
        workpiece.part.PartitionCellByPlaneThreePoints(cells=pickedCells, point1=(mm(0), self.b_height, mm(1)), point2=(mm(1), self.b_height, mm(0)), point3=(mm(1), self.b_height, mm(1)))

        c = self.part.cells
        pickedCells = c.getSequenceFromMask(mask=('[#9 ]', ), ) # TODO replace
        e1, v1, d1 = p.edges, p.vertices, p.datums
        x,y,z = 0.5*(self.b_width + self.w_width) - self.t, 0, 0
        p.PartitionCellByPlanePointNormal(normal=e1[5], cells=pickedCells, 
            point=(x,y,z))


    def set_material(self, material):
        model = mdb.models['Model-1']
        model.HomogeneousSolidSection(name=material.name, material=material.name, thickness=None)
        c = self.part.cells
        region = self.part.Set(cells=c, name='Material region')
        self.part.SectionAssignment(region=region, sectionName=material.name, offset=0.0, 
            offsetType=MIDDLE_SURFACE, offsetField='', thicknessAssignment=FROM_SECTION)
        

    def mesh(self):
        e = self.part.edges
        e1_point = (0.5* self.b_width, self.b_height, 0)
        e2_point = (0.5* self.b_width, self.b_height+ self.w_height, 0)
        e3_point = (0.5* self.b_width, self.b_height, self.length)
        e4_point = (0.5* self.b_width, self.b_height+ self.w_height, self.length)
        pickedEdges = e.findAt((e1_point, ), (e2_point, ), (e3_point, ), (e4_point, ))
        self.part.seedEdgeByNumber(edges=pickedEdges, number=int(0.3*self.w_width * 1e3), constraint=FINER)
        
        e1_point = (0.5* (self.b_width - self.w_width), self.b_height + 0.5 * self.w_height, 0)
        e2_point = (0.5* (self.b_width + self.w_width), self.b_height + 0.5 * self.w_height, 0)
        e3_point = (0.5* (self.b_width - self.w_width), self.b_height + 0.5 * self.w_height, self.length)
        e4_point = (0.5* (self.b_width + self.w_width), self.b_height + 0.5 * self.w_height, self.length)
        pickedEdges = e.findAt((e1_point, ), (e2_point, ), (e3_point, ), (e4_point, ))
        self.part.seedEdgeByNumber(edges=pickedEdges, number=int(0.5*self.w_height * 1e3), constraint=FINER)

        e1_point = (0.5* (self.b_width - self.w_width), self.b_height, 0.5 * self.length)
        e2_point = (0.5* (self.b_width + self.w_width), self.b_height, 0.5 * self.length)
        e3_point = (0.5* (self.b_width - self.w_width), self.b_height + self.w_height, 0.5 * self.length)
        e4_point = (0.5* (self.b_width + self.w_width), self.b_height + self.w_height, 0.5 * self.length)
        pickedEdges = e.findAt((e1_point, ), (e2_point, ), (e3_point, ), (e4_point, ))
        self.part.seedEdgeByNumber(edges=pickedEdges, number=int(1*self.length * 1e3), constraint=FINER)
        
        e1_point = (0.5* (self.b_width + self.w_width) - 0.5*self.t , self.b_height, 0)
        print("Len:"+str(len(e1_point)))
        self.part.DatumPointByCoordinate(coords=e1_point)
        pickedEdges = e.findAt((e1_point, ), )
        self.part.seedEdgeByNumber(edges=pickedEdges, number=int(6*self.t * 1e3), constraint=FINER)
        self.part.seedPart(size=0.01, deviationFactor=0.051, minSizeFactor=0.1)
    
        c = self.part.cells
        pickedRegions = c.findAt((self.c1_point, ), (self.c2_point, ), (self.c3_point, ))
        self.part.setMeshControls(regions=pickedRegions, technique=SWEEP, algorithm=ADVANCING_FRONT)


        elemType1 = mesh.ElemType(elemCode=C3D8RT, elemLibrary=EXPLICIT, 
            secondOrderAccuracy=OFF, distortionControl=ON, 
            lengthRatio=0.100000001490116, elemDeletion=ON)
        elemType2 = mesh.ElemType(elemCode=C3D6T, elemLibrary=EXPLICIT)
        elemType3 = mesh.ElemType(elemCode=C3D4T, elemLibrary=EXPLICIT)
        
        pickedRegions =(self.part.cells, )
        self.part.setElementType(regions=pickedRegions, elemTypes=(elemType1, elemType2, elemType3))
        self.part.generateMesh()



class Assembly:

    def __init__(self, workpiece, tool):        
        model = mdb.models['Model-1']
        self.workpiece = workpiece
        self.tool = tool
        self.a = model.rootAssembly
        self.a.DatumCsysByDefault(CARTESIAN)
        self.a.Instance(name=workpiece.name, part=workpiece.part, dependent=ON)
        self.a.Instance(name=tool.name, part=tool.part, dependent=ON)
        self.a.rotate(instanceList=('Tool', ), axisPoint=(0.0, 0.0, 0.0), axisDirection=(0.00625, 0.0, 0.0), angle=90.0)
        self.a.rotate(instanceList=('Tool', ), axisPoint=(0.0, 0.0, 0.0), axisDirection=(0.0, 1.0, 0.0), angle=45.0)
        self.workpieceInstance = self.a.instances[workpiece.name]
        self.toolInstance = self.a.instances[tool.name]
        self.a.translate(instanceList=('Tool', ), vector=(0.0645, 0.07, -0.0037))
        self.a.translate(instanceList=('Tool', ), vector=(0, 0.005, 0))
        self.a.translate(instanceList=('Tool', ), vector=(0.0, 0.0, -0.0125))
        self.a.translate(instanceList=(tool.name, ), vector=(mm(-18), 0.0, 0))
        self.a.translate(instanceList=('Tool', ), vector=(-0.003, -0.02, 0.0))
        self.surface_workpiece()
        self.surface_tool()
        session.viewports['Viewport: 1'].setValues(displayedObject=self.a)
        session.viewports['Viewport: 1'].view.fitView()    

    def surface_workpiece(self):
        mdb.models['Model-1'].rootAssembly.Surface(face1Elements=
            mdb.models['Model-1'].rootAssembly.instances['Workpiece'].elements.getSequenceFromMask(
            mask=('[#ffffffff:890 #fffff #0:48 #fffffffc #ffffffff:30 #3ff ]', ), ), 
            face2Elements=
            mdb.models['Model-1'].rootAssembly.instances['Workpiece'].elements.getSequenceFromMask(
            mask=('[#0:46 #f0000000 #ffffffff:890 #ffff ]', ), ), face3Elements=
            mdb.models['Model-1'].rootAssembly.instances['Workpiece'].elements.getSequenceFromMask(
            mask=('[#ffffffff:937 #ffff ]', ), ), face4Elements=
            mdb.models['Model-1'].rootAssembly.instances['Workpiece'].elements.getSequenceFromMask(
            mask=(
            '[#ffffffff #fffdffff #ffffffff #fffffff7 #ffdfffff #ffffffff #ffffff7f', 
            ' #fdffffff #ffffffff #fffff7ff #dfffffff #ffffffff #ffff7fff #ffffffff', 
            ' #fffffffd #fff7ffff #ffffffff #ffffffdf #ff7fffff #ffffffff #fffffdff', 
            ' #f7ffffff #ffffffff #ffffdfff #7fffffff #ffffffff #fffdffff #ffffffff', 
            ' #fffffff7 #ffdfffff #ffffffff #ffffff7f #fdffffff #ffffffff #fffff7ff', 
            ' #dfffffff #ffffffff #ffff7fff #ffffffff #fffffffd #fff7ffff #ffffffff', 
            ' #ffffffdf #ff7fffff #ffffffff #fffffdff #f7ffffff #ffffffff #ffffdfff', 
            ' #7fffffff #ffffffff #fffdffff #ffffffff #fffffff7 #ffdfffff #ffffffff', 
            ' #ffffff7f #fdffffff #ffffffff #fffff7ff #dfffffff #ffffffff #ffff7fff', 
            ' #ffffffff #fffffffd #fff7ffff #ffffffff #ffffffdf #ff7fffff #ffffffff', 
            ' #fffffdff #f7ffffff #ffffffff #ffffdfff #7fffffff #ffffffff #fffdffff', 
            ' #ffffffff #fffffff7 #ffdfffff #ffffffff #ffffff7f #fdffffff #ffffffff', 
            ' #fffff7ff #dfffffff #ffffffff #ffff7fff #ffffffff #fffffffd #fff7ffff', 
            ' #ffffffff #ffffffdf #ff7fffff #ffffffff #fffffdff #f7ffffff #ffffffff', 
            ' #ffffdfff #7fffffff #ffffffff #fffdffff #ffffffff #fffffff7 #ffdfffff', 
            ' #ffffffff #ffffff7f #fdffffff #ffffffff #fffff7ff #dfffffff #ffffffff', 
            ' #ffff7fff #ffffffff #fffffffd #fff7ffff #ffffffff #ffffffdf #ff7fffff', 
            ' #ffffffff #fffffdff #f7ffffff #ffffffff #ffffdfff #7fffffff #ffffffff', 
            ' #fffdffff #ffffffff #fffffff7 #ffdfffff #ffffffff #ffffff7f #fdffffff', 
            ' #ffffffff #fffff7ff #dfffffff #ffffffff #ffff7fff #ffffffff #fffffffd', 
            ' #fff7ffff #ffffffff #ffffffdf #ff7fffff #ffffffff #fffffdff #f7ffffff', 
            ' #ffffffff #ffffdfff #7fffffff #ffffffff #fffdffff #ffffffff #fffffff7', 
            ' #ffdfffff #ffffffff #ffffff7f #fdffffff #ffffffff #fffff7ff #dfffffff', 
            ' #ffffffff #ffff7fff #ffffffff #fffffffd #fff7ffff #ffffffff #ffffffdf', 
            ' #ff7fffff #ffffffff #fffffdff #f7ffffff #ffffffff #ffffdfff #7fffffff', 
            ' #ffffffff #fffdffff #ffffffff #fffffff7 #ffdfffff #ffffffff #ffffff7f', 
            ' #fdffffff #ffffffff #fffff7ff #dfffffff #ffffffff #ffff7fff #ffffffff', 
            ' #fffffffd #fff7ffff #ffffffff #ffffffdf #ff7fffff #ffffffff #fffffdff', 
            ' #f7ffffff #ffffffff #ffffdfff #7fffffff #ffffffff #fffdffff #ffffffff', 
            ' #fffffff7 #ffdfffff #ffffffff #ffffff7f #fdffffff #ffffffff #fffff7ff', 
            ' #dfffffff #ffffffff #ffff7fff #ffffffff #fffffffd #fff7ffff #ffffffff', 
            ' #ffffffdf #ff7fffff #ffffffff #fffffdff #f7ffffff #ffffffff #ffffdfff', 
            ' #7fffffff #ffffffff #fffdffff #ffffffff #fffffff7 #ffdfffff #ffffffff', 
            ' #ffffff7f #fdffffff #ffffffff #fffff7ff #dfffffff #ffffffff #ffff7fff', 
            ' #ffffffff #fffffffd #fff7ffff #ffffffff #ffffffdf #ff7fffff #ffffffff', 
            ' #fffffdff #f7ffffff #ffffffff #ffffdfff #7fffffff #ffffffff #fffdffff', 
            ' #ffffffff #fffffff7 #ffdfffff #ffffffff #ffffff7f #fdffffff #ffffffff', 
            ' #fffff7ff #dfffffff #ffffffff #ffff7fff #ffffffff #fffffffd #fff7ffff', 
            ' #ffffffff #ffffffdf #ff7fffff #ffffffff #fffffdff #f7ffffff #ffffffff', 
            ' #ffffdfff #7fffffff #ffffffff #fffdffff #ffffffff #fffffff7 #ffdfffff', 
            ' #ffffffff #ffffff7f #fdffffff #ffffffff #fffff7ff #dfffffff #ffffffff', 
            ' #ffff7fff #ffffffff #fffffffd #fff7ffff #ffffffff #ffffffdf #ff7fffff', 
            ' #ffffffff #fffffdff #f7ffffff #ffffffff #ffffdfff #7fffffff #ffffffff', 
            ' #fffdffff #ffffffff #fffffff7 #ffdfffff #ffffffff #ffffff7f #fdffffff', 
            ' #ffffffff #fffff7ff #dfffffff #ffffffff #ffff7fff #ffffffff #fffffffd', 
            ' #fff7ffff #ffffffff #ffffffdf #ff7fffff #ffffffff #fffffdff #f7ffffff', 
            ' #ffffffff #ffffdfff #7fffffff #ffffffff #fffdffff #ffffffff #fffffff7', 
            ' #ffdfffff #ffffffff #ffffff7f #fdffffff #ffffffff #fffff7ff #dfffffff', 
            ' #ffffffff #ffff7fff #ffffffff #fffffffd #fff7ffff #ffffffff #ffffffdf', 
            ' #ff7fffff #ffffffff #fffffdff #f7ffffff #ffffffff #ffffdfff #7fffffff', 
            ' #ffffffff #fffdffff #ffffffff #fffffff7 #ffdfffff #ffffffff #ffffff7f', 
            ' #fdffffff #ffffffff #fffff7ff #dfffffff #ffffffff #ffff7fff #ffffffff', 
            ' #fffffffd #fff7ffff #ffffffff #ffffffdf #ff7fffff #ffffffff #fffffdff', 
            ' #f7ffffff #ffffffff #ffffdfff #7fffffff #ffffffff #fffdffff #ffffffff', 
            ' #fffffff7 #ffdfffff #ffffffff #ffffff7f #fdffffff #ffffffff #fffff7ff', 
            ' #dfffffff #ffffffff #ffff7fff #ffffffff #fffffffd #fff7ffff #ffffffff', 
            ' #ffffffdf #ff7fffff #ffffffff #fffffdff #f7ffffff #ffffffff #ffffdfff', 
            ' #7fffffff #ffffffff #fffdffff #ffffffff #fffffff7 #ffdfffff #ffffffff', 
            ' #ffffff7f #fdffffff #ffffffff #fffff7ff #dfffffff #ffffffff #ffff7fff', 
            ' #ffffffff #fffffffd #fff7ffff #ffffffff #ffffffdf #ff7fffff #ffffffff', 
            ' #fffffdff #f7ffffff #ffffffff #ffffdfff #7fffffff #ffffffff #fffdffff', 
            ' #ffffffff #fffffff7 #ffdfffff #ffffffff #ffffff7f #fdffffff #ffffffff', 
            ' #fffff7ff #dfffffff #ffffffff #ffff7fff #ffffffff #fffffffd #fff7ffff', 
            ' #ffffffff #ffffffdf #ff7fffff #ffffffff #fffffdff #f7ffffff #ffffffff', 
            ' #ffffdfff #7fffffff #ffffffff #fffdffff #ffffffff #fffffff7 #ffdfffff', 
            ' #ffffffff #ffffff7f #fdffffff #ffffffff #fffff7ff #dfffffff #ffffffff', 
            ' #ffff7fff #ffffffff #fffffffd #fff7ffff #ffffffff #ffffffdf #ff7fffff', 
            ' #ffffffff #fffffdff #f7ffffff #ffffffff #ffffdfff #7fffffff #ffffffff', 
            ' #fffdffff #ffffffff #fffffff7 #ffdfffff #ffffffff #ffffff7f #fdffffff', 
            ' #ffffffff #fffff7ff #dfffffff #ffffffff #ffff7fff #ffffffff #fffffffd', 
            ' #fff7ffff #ffffffff #ffffffdf #ff7fffff #ffffffff #fffffdff #f7ffffff', 
            ' #ffffffff #ffffdfff #7fffffff #ffffffff #fffdffff #ffffffff #fffffff7', 
            ' #ffdfffff #ffffffff #ffffff7f #fdffffff #ffffffff #fffff7ff #dfffffff', 
            ' #ffffffff #ffff7fff #ffffffff #fffffffd #fff7ffff #ffffffff #ffffffdf', 
            ' #ff7fffff #ffffffff #fffffdff #f7ffffff #ffffffff #ffffdfff #7fffffff', 
            ' #ffffffff #fffdffff #ffffffff #fffffff7 #ffdfffff #ffffffff #ffffff7f', 
            ' #fdffffff #ffffffff #fffff7ff #dfffffff #ffffffff #ffff7fff #ffffffff', 
            ' #fffffffd #fff7ffff #ffffffff #ffffffdf #ff7fffff #ffffffff #fffffdff', 
            ' #f7ffffff #ffffffff #ffffdfff #7fffffff #ffffffff #fffdffff #ffffffff', 
            ' #fffffff7 #ffdfffff #ffffffff #ffffff7f #fdffffff #ffffffff #fffff7ff', 
            ' #dfffffff #ffffffff #ffff7fff #ffffffff #fffffffd #fff7ffff #ffffffff', 
            ' #ffffffdf #ff7fffff #ffffffff #fffffdff #f7ffffff #ffffffff #ffffdfff', 
            ' #7fffffff #ffffffff #fffdffff #ffffffff #fffffff7 #ffdfffff #ffffffff', 
            ' #ffffff7f #fdffffff #ffffffff #fffff7ff #dfffffff #ffffffff #ffff7fff', 
            ' #ffffffff #fffffffd #fff7ffff #ffffffff #ffffffdf #ff7fffff #ffffffff', 
            ' #fffffdff #f7ffffff #ffffffff #ffffdfff #7fffffff #ffffffff #fffdffff', 
            ' #ffffffff #fffffff7 #ffdfffff #ffffffff #ffffff7f #fdffffff #ffffffff', 
            ' #fffff7ff #dfffffff #ffffffff #ffff7fff #ffffffff #fffffffd #fff7ffff', 
            ' #ffffffff #ffffffdf #ff7fffff #ffffffff #fffffdff #f7ffffff #ffffffff', 
            ' #ffffdfff #7fffffff #ffffffff #fffdffff #ffffffff #fffffff7 #ffdfffff', 
            ' #ffffffff #ffffff7f #fdffffff #ffffffff #fffff7ff #dfffffff #ffffffff', 
            ' #ffff7fff #ffffffff #fffffffd #fff7ffff #ffffffff #ffffffdf #ff7fffff', 
            ' #ffffffff #fffffdff #f7ffffff #ffffffff #ffffdfff #7fffffff #ffffffff', 
            ' #fffdffff #ffffffff #fffffff7 #ffdfffff #ffffffff #ffffff7f #fdffffff', 
            ' #ffffffff #fffff7ff #dfffffff #ffffffff #ffff7fff #ffffffff #fffffffd', 
            ' #fff7ffff #ffffffff #ffffffdf #ff7fffff #ffffffff #fffffdff #f7ffffff', 
            ' #ffffffff #ffffdfff #7fffffff #ffffffff #fffdffff #ffffffff #fffffff7', 
            ' #ffdfffff #ffffffff #ffffff7f #fdffffff #ffffffff #fffff7ff #dfffffff', 
            ' #ffffffff #ffff7fff #ffffffff #fffffffd #fff7ffff #ffffffff #ffffffdf', 
            ' #ff7fffff #ffffffff #fffffdff #f7ffffff #ffffffff #ffffdfff #7fffffff', 
            ' #ffffffff #fffdffff #ffffffff #fffffff7 #ffdfffff #ffffffff #ffffff7f', 
            ' #fdffffff #ffffffff #fffff7ff #dfffffff #ffffffff #ffff7fff #ffffffff', 
            ' #fffffffd #fff7ffff #ffffffff #ffffffdf #ff7fffff #ffffffff #fffffdff', 
            ' #f7ffffff #ffffffff #ffffdfff #7fffffff #ffffffff #fffdffff #ffffffff', 
            ' #fffffff7 #ffdfffff #ffffffff #ffffff7f #fdffffff #ffffffff #fffff7ff', 
            ' #dfffffff #ffffffff #ffff7fff #ffffffff #fffffffd #fff7ffff #ffffffff', 
            ' #ffffffdf #ff7fffff #ffffffff #fffffdff #f7ffffff #ffffffff #ffffdfff', 
            ' #7fffffff #ffffffff #fffdffff #ffffffff #fffffff7 #ffdfffff #ffffffff', 
            ' #ffffff7f #fdffffff #ffffffff #fffff7ff #dfffffff #ffffffff #ffff7fff', 
            ' #ffffffff #fffffffd #fff7ffff #ffffffff #ffffffdf #ff7fffff #ffffffff', 
            ' #fffffdff #f7ffffff #ffffffff #ffffdfff #7fffffff #ffffffff #fffdffff', 
            ' #ffffffff #fffffff7 #ffdfffff #ffffffff #ffffff7f #fdffffff #ffffffff', 
            ' #fffff7ff #dfffffff #ffffffff #ffff7fff #ffffffff #fffffffd #fff7ffff', 
            ' #ffffffff #ffffffdf #ff7fffff #ffffffff #fffffdff #f7ffffff #ffffffff', 
            ' #ffffdfff #7fffffff #ffffffff #fffdffff #ffffffff #fffffff7 #ffdfffff', 
            ' #ffffffff #ffffff7f #fdffffff #ffffffff #fffff7ff #dfffffff #ffffffff', 
            ' #ffff7fff #ffffffff #fffffffd #fff7ffff #ffffffff #ffffffdf #ff7fffff', 
            ' #ffffffff #fffffdff #f7ffffff #ffffffff #ffffdfff #7fffffff #ffffffff', 
            ' #fffdffff #ffffffff #fffffff7 #ffdfffff #ffffffff #ffffff7f #fdffffff', 
            ' #ffffffff #fffff7ff #dfffffff #ffffffff #ffff7fff #ffffffff #fffffffd', 
            ' #fff7ffff #ffffffff #ffffffdf #ff7fffff #ffffffff #fffffdff #f7ffffff', 
            ' #ffffffff #ffffdfff #7fffffff #ffffffff #fffdffff #ffffffff #fffffff7', 
            ' #ffdfffff #ffffffff #ffffff7f #fdffffff #ffffffff #fffff7ff #dfffffff', 
            ' #ffffffff #ffff7fff #ffffffff #fffffffd #fff7ffff #ffffffff #ffffffdf', 
            ' #ff7fffff #ffffffff #fffffdff #f7ffffff #ffffffff #ffffdfff #7fffffff', 
            ' #ffffffff #fffdffff #ffffffff #fffffff7 #ffdfffff #ffffffff #ffffff7f', 
            ' #fdffffff #ffffffff #fffff7ff #dfffffff #ffffffff #ffff7fff #ffffffff', 
            ' #fffffffd #fff7ffff #ffffffff #ffffffdf #ff7fffff #ffffffff #fffffdff', 
            ' #f7ffffff #ffffffff #ffffdfff #7fffffff #ffffffff #fffdffff #ffffffff', 
            ' #fffffff7 #ffdfffff #ffffffff #ffffff7f #fdffffff #ffffffff #fffff7ff', 
            ' #dfffffff #ffffffff #ffff7fff #ffffffff #fffffffd #fff7ffff #ffffffff', 
            ' #ffffffdf #ff7fffff #ffffffff #fffffdff #f7ffffff #ffffffff #ffffdfff', 
            ' #7fffffff #ffffffff #fffdffff #ffffffff #fffffff7 #ffdfffff #ffffffff', 
            ' #ffffff7f #fdffffff #ffffffff #fffff7ff #dfffffff #ffffffff #7fff ]', ), 
            ), face5Elements=
            mdb.models['Model-1'].rootAssembly.instances['Workpiece'].elements.getSequenceFromMask(
            mask=(
            '[#ffffffff:45 #3ff #f0000000 #ffffffff:45 #3f #ff000000 #ffffffff:45', 
            ' #3 #fff00000 #ffffffff:44 #3fffffff #0 #ffff0000 #ffffffff:44', 
            ' #3ffffff #0 #fffff000 #ffffffff:44 #3fffff #0 #ffffff00', 
            ' #ffffffff:44 #3ffff #0 #fffffff0 #ffffffff:44 #3fff #0', 
            ' #ffffffff:45 #3ff #f0000000 #ffffffff:45 #3f #ff000000 #ffffffff:45', 
            ' #3 #fff00000 #ffffffff:44 #3fffffff #0 #ffff0000 #ffffffff:44', 
            ' #3ffffff #0 #fffff000 #ffffffff:44 #3fffff #0 #ffffff00', 
            ' #ffffffff:44 #3ffff #0 #fffffff0 #ffffffff:44 #3fff #0', 
            ' #ffffffff:45 #3ff #f0000000 #ffffffff:45 #3f #ff000000 #ffffffff:45', 
            ' #3 #fff00000 #ffffffff:44 #3fffffff ]', ), ), face6Elements=
            mdb.models['Model-1'].rootAssembly.instances['Workpiece'].elements.getSequenceFromMask(
            mask=(
            '[#fffffffe #fffbffff #ffffffff #ffffffef #ffbfffff #ffffffff #fffffeff', 
            ' #fbffffff #ffffffff #ffffefff #bfffffff #ffffffff #fffeffff #ffffffff', 
            ' #fffffffb #ffefffff #ffffffff #ffffffbf #feffffff #ffffffff #fffffbff', 
            ' #efffffff #ffffffff #ffffbfff #ffffffff #fffffffe #fffbffff #ffffffff', 
            ' #ffffffef #ffbfffff #ffffffff #fffffeff #fbffffff #ffffffff #ffffefff', 
            ' #bfffffff #ffffffff #fffeffff #ffffffff #fffffffb #ffefffff #ffffffff', 
            ' #ffffffbf #feffffff #ffffffff #fffffbff #efffffff #ffffffff #ffffbfff', 
            ' #ffffffff #fffffffe #fffbffff #ffffffff #ffffffef #ffbfffff #ffffffff', 
            ' #fffffeff #fbffffff #ffffffff #ffffefff #bfffffff #ffffffff #fffeffff', 
            ' #ffffffff #fffffffb #ffefffff #ffffffff #ffffffbf #feffffff #ffffffff', 
            ' #fffffbff #efffffff #ffffffff #ffffbfff #ffffffff #fffffffe #fffbffff', 
            ' #ffffffff #ffffffef #ffbfffff #ffffffff #fffffeff #fbffffff #ffffffff', 
            ' #ffffefff #bfffffff #ffffffff #fffeffff #ffffffff #fffffffb #ffefffff', 
            ' #ffffffff #ffffffbf #feffffff #ffffffff #fffffbff #efffffff #ffffffff', 
            ' #ffffbfff #ffffffff #fffffffe #fffbffff #ffffffff #ffffffef #ffbfffff', 
            ' #ffffffff #fffffeff #fbffffff #ffffffff #ffffefff #bfffffff #ffffffff', 
            ' #fffeffff #ffffffff #fffffffb #ffefffff #ffffffff #ffffffbf #feffffff', 
            ' #ffffffff #fffffbff #efffffff #ffffffff #ffffbfff #ffffffff #fffffffe', 
            ' #fffbffff #ffffffff #ffffffef #ffbfffff #ffffffff #fffffeff #fbffffff', 
            ' #ffffffff #ffffefff #bfffffff #ffffffff #fffeffff #ffffffff #fffffffb', 
            ' #ffefffff #ffffffff #ffffffbf #feffffff #ffffffff #fffffbff #efffffff', 
            ' #ffffffff #ffffbfff #ffffffff #fffffffe #fffbffff #ffffffff #ffffffef', 
            ' #ffbfffff #ffffffff #fffffeff #fbffffff #ffffffff #ffffefff #bfffffff', 
            ' #ffffffff #fffeffff #ffffffff #fffffffb #ffefffff #ffffffff #ffffffbf', 
            ' #feffffff #ffffffff #fffffbff #efffffff #ffffffff #ffffbfff #ffffffff', 
            ' #fffffffe #fffbffff #ffffffff #ffffffef #ffbfffff #ffffffff #fffffeff', 
            ' #fbffffff #ffffffff #ffffefff #bfffffff #ffffffff #fffeffff #ffffffff', 
            ' #fffffffb #ffefffff #ffffffff #ffffffbf #feffffff #ffffffff #fffffbff', 
            ' #efffffff #ffffffff #ffffbfff #ffffffff #fffffffe #fffbffff #ffffffff', 
            ' #ffffffef #ffbfffff #ffffffff #fffffeff #fbffffff #ffffffff #ffffefff', 
            ' #bfffffff #ffffffff #fffeffff #ffffffff #fffffffb #ffefffff #ffffffff', 
            ' #ffffffbf #feffffff #ffffffff #fffffbff #efffffff #ffffffff #ffffbfff', 
            ' #ffffffff #fffffffe #fffbffff #ffffffff #ffffffef #ffbfffff #ffffffff', 
            ' #fffffeff #fbffffff #ffffffff #ffffefff #bfffffff #ffffffff #fffeffff', 
            ' #ffffffff #fffffffb #ffefffff #ffffffff #ffffffbf #feffffff #ffffffff', 
            ' #fffffbff #efffffff #ffffffff #ffffbfff #ffffffff #fffffffe #fffbffff', 
            ' #ffffffff #ffffffef #ffbfffff #ffffffff #fffffeff #fbffffff #ffffffff', 
            ' #ffffefff #bfffffff #ffffffff #fffeffff #ffffffff #fffffffb #ffefffff', 
            ' #ffffffff #ffffffbf #feffffff #ffffffff #fffffbff #efffffff #ffffffff', 
            ' #ffffbfff #ffffffff #fffffffe #fffbffff #ffffffff #ffffffef #ffbfffff', 
            ' #ffffffff #fffffeff #fbffffff #ffffffff #ffffefff #bfffffff #ffffffff', 
            ' #fffeffff #ffffffff #fffffffb #ffefffff #ffffffff #ffffffbf #feffffff', 
            ' #ffffffff #fffffbff #efffffff #ffffffff #ffffbfff #ffffffff #fffffffe', 
            ' #fffbffff #ffffffff #ffffffef #ffbfffff #ffffffff #fffffeff #fbffffff', 
            ' #ffffffff #ffffefff #bfffffff #ffffffff #fffeffff #ffffffff #fffffffb', 
            ' #ffefffff #ffffffff #ffffffbf #feffffff #ffffffff #fffffbff #efffffff', 
            ' #ffffffff #ffffbfff #ffffffff #fffffffe #fffbffff #ffffffff #ffffffef', 
            ' #ffbfffff #ffffffff #fffffeff #fbffffff #ffffffff #ffffefff #bfffffff', 
            ' #ffffffff #fffeffff #ffffffff #fffffffb #ffefffff #ffffffff #ffffffbf', 
            ' #feffffff #ffffffff #fffffbff #efffffff #ffffffff #ffffbfff #ffffffff', 
            ' #fffffffe #fffbffff #ffffffff #ffffffef #ffbfffff #ffffffff #fffffeff', 
            ' #fbffffff #ffffffff #ffffefff #bfffffff #ffffffff #fffeffff #ffffffff', 
            ' #fffffffb #ffefffff #ffffffff #ffffffbf #feffffff #ffffffff #fffffbff', 
            ' #efffffff #ffffffff #ffffbfff #ffffffff #fffffffe #fffbffff #ffffffff', 
            ' #ffffffef #ffbfffff #ffffffff #fffffeff #fbffffff #ffffffff #ffffefff', 
            ' #bfffffff #ffffffff #fffeffff #ffffffff #fffffffb #ffefffff #ffffffff', 
            ' #ffffffbf #feffffff #ffffffff #fffffbff #efffffff #ffffffff #ffffbfff', 
            ' #ffffffff #fffffffe #fffbffff #ffffffff #ffffffef #ffbfffff #ffffffff', 
            ' #fffffeff #fbffffff #ffffffff #ffffefff #bfffffff #ffffffff #fffeffff', 
            ' #ffffffff #fffffffb #ffefffff #ffffffff #ffffffbf #feffffff #ffffffff', 
            ' #fffffbff #efffffff #ffffffff #ffffbfff #ffffffff #fffffffe #fffbffff', 
            ' #ffffffff #ffffffef #ffbfffff #ffffffff #fffffeff #fbffffff #ffffffff', 
            ' #ffffefff #bfffffff #ffffffff #fffeffff #ffffffff #fffffffb #ffefffff', 
            ' #ffffffff #ffffffbf #feffffff #ffffffff #fffffbff #efffffff #ffffffff', 
            ' #ffffbfff #ffffffff #fffffffe #fffbffff #ffffffff #ffffffef #ffbfffff', 
            ' #ffffffff #fffffeff #fbffffff #ffffffff #ffffefff #bfffffff #ffffffff', 
            ' #fffeffff #ffffffff #fffffffb #ffefffff #ffffffff #ffffffbf #feffffff', 
            ' #ffffffff #fffffbff #efffffff #ffffffff #ffffbfff #ffffffff #fffffffe', 
            ' #fffbffff #ffffffff #ffffffef #ffbfffff #ffffffff #fffffeff #fbffffff', 
            ' #ffffffff #ffffefff #bfffffff #ffffffff #fffeffff #ffffffff #fffffffb', 
            ' #ffefffff #ffffffff #ffffffbf #feffffff #ffffffff #fffffbff #efffffff', 
            ' #ffffffff #ffffbfff #ffffffff #fffffffe #fffbffff #ffffffff #ffffffef', 
            ' #ffbfffff #ffffffff #fffffeff #fbffffff #ffffffff #ffffefff #bfffffff', 
            ' #ffffffff #fffeffff #ffffffff #fffffffb #ffefffff #ffffffff #ffffffbf', 
            ' #feffffff #ffffffff #fffffbff #efffffff #ffffffff #ffffbfff #ffffffff', 
            ' #fffffffe #fffbffff #ffffffff #ffffffef #ffbfffff #ffffffff #fffffeff', 
            ' #fbffffff #ffffffff #ffffefff #bfffffff #ffffffff #fffeffff #ffffffff', 
            ' #fffffffb #ffefffff #ffffffff #ffffffbf #feffffff #ffffffff #fffffbff', 
            ' #efffffff #ffffffff #ffffbfff #ffffffff #fffffffe #fffbffff #ffffffff', 
            ' #ffffffef #ffbfffff #ffffffff #fffffeff #fbffffff #ffffffff #ffffefff', 
            ' #bfffffff #ffffffff #fffeffff #ffffffff #fffffffb #ffefffff #ffffffff', 
            ' #ffffffbf #feffffff #ffffffff #fffffbff #efffffff #ffffffff #ffffbfff', 
            ' #ffffffff #fffffffe #fffbffff #ffffffff #ffffffef #ffbfffff #ffffffff', 
            ' #fffffeff #fbffffff #ffffffff #ffffefff #bfffffff #ffffffff #fffeffff', 
            ' #ffffffff #fffffffb #ffefffff #ffffffff #ffffffbf #feffffff #ffffffff', 
            ' #fffffbff #efffffff #ffffffff #ffffbfff #ffffffff #fffffffe #fffbffff', 
            ' #ffffffff #ffffffef #ffbfffff #ffffffff #fffffeff #fbffffff #ffffffff', 
            ' #ffffefff #bfffffff #ffffffff #fffeffff #ffffffff #fffffffb #ffefffff', 
            ' #ffffffff #ffffffbf #feffffff #ffffffff #fffffbff #efffffff #ffffffff', 
            ' #ffffbfff #ffffffff #fffffffe #fffbffff #ffffffff #ffffffef #ffbfffff', 
            ' #ffffffff #fffffeff #fbffffff #ffffffff #ffffefff #bfffffff #ffffffff', 
            ' #fffeffff #ffffffff #fffffffb #ffefffff #ffffffff #ffffffbf #feffffff', 
            ' #ffffffff #fffffbff #efffffff #ffffffff #ffffbfff #ffffffff #fffffffe', 
            ' #fffbffff #ffffffff #ffffffef #ffbfffff #ffffffff #fffffeff #fbffffff', 
            ' #ffffffff #ffffefff #bfffffff #ffffffff #fffeffff #ffffffff #fffffffb', 
            ' #ffefffff #ffffffff #ffffffbf #feffffff #ffffffff #fffffbff #efffffff', 
            ' #ffffffff #ffffbfff #ffffffff #fffffffe #fffbffff #ffffffff #ffffffef', 
            ' #ffbfffff #ffffffff #fffffeff #fbffffff #ffffffff #ffffefff #bfffffff', 
            ' #ffffffff #fffeffff #ffffffff #fffffffb #ffefffff #ffffffff #ffffffbf', 
            ' #feffffff #ffffffff #fffffbff #efffffff #ffffffff #ffffbfff #ffffffff', 
            ' #fffffffe #fffbffff #ffffffff #ffffffef #ffbfffff #ffffffff #fffffeff', 
            ' #fbffffff #ffffffff #ffffefff #bfffffff #ffffffff #fffeffff #ffffffff', 
            ' #fffffffb #ffefffff #ffffffff #ffffffbf #feffffff #ffffffff #fffffbff', 
            ' #efffffff #ffffffff #ffffbfff #ffffffff #fffffffe #fffbffff #ffffffff', 
            ' #ffffffef #ffbfffff #ffffffff #fffffeff #fbffffff #ffffffff #ffffefff', 
            ' #bfffffff #ffffffff #fffeffff #ffffffff #fffffffb #ffefffff #ffffffff', 
            ' #ffffffbf #feffffff #ffffffff #fffffbff #efffffff #ffffffff #ffffbfff', 
            ' #ffffffff #fffffffe #fffbffff #ffffffff #ffffffef #ffbfffff #ffffffff', 
            ' #fffffeff #fbffffff #ffffffff #ffffefff #bfffffff #ffffffff #fffeffff', 
            ' #ffffffff #fffffffb #ffefffff #ffffffff #ffffffbf #feffffff #ffffffff', 
            ' #fffffbff #efffffff #ffffffff #ffffbfff #ffffffff #fffffffe #fffbffff', 
            ' #ffffffff #ffffffef #ffbfffff #ffffffff #fffffeff #fbffffff #ffffffff', 
            ' #ffffefff #bfffffff #ffffffff #fffeffff #ffffffff #fffffffb #ffefffff', 
            ' #ffffffff #ffffffbf #feffffff #ffffffff #fffffbff #efffffff #ffffffff', 
            ' #ffffbfff #ffffffff #fffffffe #fffbffff #ffffffff #ffffffef #ffbfffff', 
            ' #ffffffff #fffffeff #fbffffff #ffffffff #ffffefff #bfffffff #ffffffff', 
            ' #fffeffff #ffffffff #fffffffb #ffefffff #ffffffff #ffffffbf #feffffff', 
            ' #ffffffff #fffffbff #efffffff #ffffffff #ffffbfff #ffffffff #fffffffe', 
            ' #fffbffff #ffffffff #ffffffef #ffbfffff #ffffffff #fffffeff #fbffffff', 
            ' #ffffffff #ffffefff #bfffffff #ffffffff #fffeffff #ffffffff #fffffffb', 
            ' #ffefffff #ffffffff #ffffffbf #feffffff #ffffffff #fffffbff #efffffff', 
            ' #ffffffff #ffffbfff #ffffffff #fffffffe #fffbffff #ffffffff #ffffffef', 
            ' #ffbfffff #ffffffff #fffffeff #fbffffff #ffffffff #ffffefff #bfffffff', 
            ' #ffffffff #fffeffff #ffffffff #fffffffb #ffefffff #ffffffff #ffffffbf', 
            ' #feffffff #ffffffff #fffffbff #efffffff #ffffffff #ffffbfff #ffffffff', 
            ' #fffffffe #fffbffff #ffffffff #ffffffef #ffbfffff #ffffffff #fffffeff', 
            ' #fbffffff #ffffffff #ffffefff #bfffffff #ffffffff #fffeffff #ffffffff', 
            ' #fffffffb #ffefffff #ffffffff #ffffffbf #feffffff #ffffffff #fffffbff', 
            ' #efffffff #ffffffff #ffffbfff #ffffffff #fffffffe #fffbffff #ffffffff', 
            ' #ffffffef #ffbfffff #ffffffff #fffffeff #fbffffff #ffffffff #ffffefff', 
            ' #bfffffff #ffffffff #fffeffff #ffffffff #fffffffb #ffefffff #ffffffff', 
            ' #ffffffbf #feffffff #ffffffff #fffffbff #efffffff #ffffffff #ffffbfff', 
            ' #ffffffff #fffffffe #fffbffff #ffffffff #ffffffef #ffbfffff #ffffffff', 
            ' #fffffeff #fbffffff #ffffffff #ffffefff #bfffffff #ffffffff #ffff ]', ), 
            ), name='Workpiece Surface')

        
    def surface_tool(self):
        mdb.models['Model-1'].rootAssembly.Surface(face1Elements=
            mdb.models['Model-1'].rootAssembly.instances['Tool'].elements.getSequenceFromMask(
            mask=(
            '[#ffffffff:2 #fffffffd #ffffffff:2 #fffffffb #ffffffff:2 #fffbffff #ffffffff:3', 
            ' #bfffffff #ffffffff:3 #fffdffbf #ffffffff:6 #ffffff7f #ffffffff:2 #fdffffff', 
            ' #ffffffff:11 #ffffffef #ffffffff:5 #fffffeff #ffffffff:11 #dfff7fff #ffffefff', 
            ' #ffffffff:64 #fdffffff #ffffffff:14 #ffffbfff #ffffffff:2 #f7ffffff #ffffffff:7', 
            ' #ffff7fff #ffffffff:4 #ffffffbf #ffffffff:11 #ffff7fff #fffffff7 #ffffffff:19', 
            ' #fdbfffff #ffffffff:15 #bfffffff #ffffffff:5 #ff7fffdf #dfffffff #ffffffff:2', 
            ' #fffffbff #ffffffff:7 #fffeffff #fdffffff #ffffffff #fff7ffff #ffffffff:16', 
            ' #ffbfffff #fffdefff #ffffffff #fffbf9ff #fdfffbff #ffffffff #ff7fffff', 
            ' #ffffffff:2 #fbffffff #ffffffff #fbffffff #fffdffff #ffffffff:3 #feffffff', 
            ' #ffffffff #ffffffdf #ffbfffff #ffffffff:2 #fffffbff #7fffffff #ffffffff:4', 
            ' #ffffefff #ffffffff:2 #fffffeff #ffffffff #fffffeff #ffffffff:4 #ffffff7f', 
            ' #ffffffff:18 #fffffff7 #ffffffff:2 #ffdfffff #ff7ff7ff #ffffffff #fffffff3', 
            ' #fffffffd #fffffffe #feffffff #ffffffff #fffbffff #ffffffff:2 #fffffff7', 
            ' #ffffffff #feffffff #ffffffff:2 #fffdffff #ffffffff #fff7efff #fefffffe', 
            ' #fffdfbff #fd3ffff7 #ffffffff #75ffffff #feffffff #fffffffe #ffff7fff', 
            ' #ffffffff:3 #dfffffff #ffffffff #fffffffd #ffffffff #ffefffff #ffffffff', 
            ' #fff7ffff #ffffffff:3 #f7ffffff #ffffffff:4 #effeffff #fffbffff #ffffffff:4', 
            ' #ffffbfff #ffffffff #bbffffff #ffbff6ff #fbefffff #ffffffff #ffffff7e', 
            ' #ffffffff #ffdfffff #fffeffff #ffffffff:6 #7bff7fff #ffdfbfff #ffffefff', 
            ' #fff7ffff #fffbffff #fffffff7 #ffffff77 #ffffbf7f #feffff7f #ffffffff', 
            ' #fffffffb #ffffffff:3 #ff3fffff #ffffffff #ffffff7f #ffffffff:5 #fffffeff', 
            ' #ffffffff #fbffffff #ffffffff #fffffeff #ffffffff:2 #ffdfffff #ffffffff:2', 
            ' #75ffffbf #dfffdffe #fffffffe #ffffffff:4 #fffbfffe #ffffffff:6 #f7ffffff', 
            ' #fffff7ff #ffffffff #ffbff7ff #fffffdff #fffeffff #ffff7fff #defffffe', 
            ' #fdfaffff #ffffffff:3 #fffffeff #ffffffff:2 #fdffffff #ffffffff #f7ffffff', 
            ' #efbfffff #ffffffff:7 #fffffeff #ffffffff:6 #dfbfffff #7fefffdf #ffffdfff', 
            ' #feffff5e #ffffffff #fffdffff #ffffffff:2 #ffffffdd #ffffffff:2 #feffffff', 
            ' #ffffffff:8 #fffffdff #fff7feff #ffffffff:9 #ffffefff #ffffffff:8 #fffff7ff', 
            ' #ffffffff #fffffdff #ffffffff #ffdfffff #ffffbfff #bffffeff #ffbf7fff', 
            ' #fffffbff #ff7dffff #ffffffff:2 #fb7fffff #fffffdef #ffffdbff #fffffffb', 
            ' #ffffffff #fffdffff #ffffffff:4 #eeffffff #ffdfffff #ffffffff:3 #fffffff7', 
            ' #dffffffc #fffdffff #fffffeff #bff6feff #fdffffff #ffffffff:2 #fffbf7ff', 
            ' #ffffffff #ffff6fff #ffffffff:7 #fffffffb #ffffffff:4 #ffffff7f #ffffffff:2', 
            ' #bfffffff #ffffffff #ffefffff #ffbfffff #ffffffff #ee7fffff #ffffff7f', 
            ' #ffffffff #fffffffd #ffffffff:3 #ffffffdf #ffffffff:2 #fffffffe #ffffffff', 
            ' #fbdfffff #fdffffbf #deffffff #ffffbfff #ffffeff7 #ffffbfdf #ffffffff:7', 
            ' #3fffffff #ffffffff:2 #ffbebfff #ffffffff:4 #fffffffe #ffffffff:4 #fffffffd', 
            ' #fffffeff #ffffffff #fffeffff #ffffffff #ffbfffff #ffffffff:4 #ffffedfb', 
            ' #ffffffff:2 #ffefbfff #fff7f7ff #fffffffa #fffffeff #ffbfffff #fff7ffff', 
            ' #f7ebfb5f #ffffffff:3 #fbffafff #f7ffffef #ffffffff #fffdffff #ffffffbf', 
            ' #ffffffff #fffffff7 #ffffffff:3 #ffbfffff #fbffffff #ffdfffff #fbffffff', 
            ' #fff7fffe #ffffffff:2 #fff5ffff #ffffffff #ffffdfff #ffffffff #fffffbff', 
            ' #ffffffef #fffffff7 #ffffffff #fdfff7ff #ffffffff:10 #fff7ffff #ffbfedff', 
            ' #fffdffff #ffffffff #fffdf7ff #d7fffdff #fffffffe #f7ffffff #ffefffff', 
            ' #ffffbfff #fffffffd #fffffff7 #ffffffff:3 #f7ffffff #ffffffff:2 #ffffefff', 
            ' #ffff7fff #ffffffff #ffffbfff #ffffffff:2 #bfffffff #ffffffff:10 #effbffff', 
            ' #ffffffff:3 #efffffff #bfdfffbf #fffffdff #ffffffff:3 #ffefffff #ffffffff:2', 
            ' #fffffffb #ffffffff #f7beffff #dfffffff #fffffff7 #fffffffd #fffffff7', 
            ' #ffffffff:5 #fff7ffff #ffffffff:2 #fdffffff #fffff7ff #feffffff #ffffffff:3', 
            ' #ffffefff #ffffffdf #ffffffff:14 #ffefffff #ffffffff #fbffffff #ffffffff', 
            ' #bfffffff #ffffefff #ffffffff:12 #dfffffff #ffffffff:17 #fffffbff #ffffffff:11', 
            ' #fffffffd #ffffffff:12 #fffffeff #ffffffff:3 #ffffffbe #ffffffff:4 #ffffdfff', 
            ' #ffffffff:15 #efffffff #ffffffff:5 #f7ffffff #ffffffff:2 #fffffb7f #ffffffff:16', 
            ' #ffffffdf #ffffffff:23 #fffff ]', ), ), face2Elements=
            mdb.models['Model-1'].rootAssembly.instances['Tool'].elements.getSequenceFromMask(
            mask=(
            '[#fffeefff #ffffffff:2 #bfffbfff #ffffffff:3 #fffffffb #ffffffff #fdffffff', 
            ' #ffffffff #ffffefff #7ffffbff #fffbffff #ffffffff:2 #fffffdff #ffffffff:2', 
            ' #7fffefff #ffffffff:5 #ffffffef #ffffffff:15 #ffffbfff #ffffffff:30 #feffffff', 
            ' #ffffffff:21 #fbffffff #ffffffff:26 #feffffff #ffffffff #fffbffff #ffffffff:3', 
            ' #ffffffbf #ffffffdf #ffffffff #fffffff7 #ffdfffff #ffffffff #ffffffdf', 
            ' #ffffffff #ffefffff #fffeffff #ffffffff:2 #ffffffdf #ffffffff:3 #ffffefff', 
            ' #ffffffff #fffbffff #ffffffff:2 #ffff7fff #ffffffff #fff7ffff:2 #ffffffff', 
            ' #ffffffef #ffffffff:3 #ffffffef #ffffffff:2 #ffefffff #ffffffff #ffefffff', 
            ' #7fbffffb #ffffffff:10 #ffffdfff #ffffffff:4 #7fffffff #ffffffff:7 #fffeffff', 
            ' #ffffffff:15 #ffffff7f #ffffffff #ffffffdf #f7ffffff #fffffff5 #fffbffff', 
            ' #ffffffff:12 #fdffffff #ffffdfff #fffffdff #ffffffff:2 #fff7ffff #ffffffff:9', 
            ' #fffffbff #ffffffff:3 #ffbfffff #ffffffff #ffff7fff #ffffffff #fffffffb', 
            ' #ffffffff:11 #fffbffff #fffdffff #ffffffff:3 #fdffffff #ffffffff #ffeffffe', 
            ' #ffffffed #ffffffff #fffbffff #7fffffff #ffffffff:2 #ffffffdf #ffffffff:4', 
            ' #fffffeff #fffbffff #ffffffef #ffffffff:2 #fffffffd #ffffffff #fffefffb', 
            ' #fffffffe #ffffffff:5 #fffffffb #ffffffff #ff7fffff #fdffffff #f7fffffb', 
            ' #9fffff7f #ffffffff #ffefffff #effffffe #ffffffff #f7bfffff #efffffff', 
            ' #ffffbfff #fffffff7:2 #ffffffff #7fffffff #ffffffff #fffffffd #fffbffff', 
            ' #ffffffff:4 #fffffff7 #ffffffff:3 #bffffeff #ffffffff:2 #fffffbff #ffe7f7ff', 
            ' #f65effff #7fffffff #ffffffff:2 #fffebfff #fffffdff #ffffffff:3 #ffffeebf', 
            ' #fffffffb #fffbffff #feffffff #ff7fffff #ffffffff:3 #fffffdff #ffffffff', 
            ' #fffffffe #bdffffdf #fffff7ff #ff7fffdf #ff7bffff #fdffffff #ffffffff:2', 
            ' #fffffffe #fdffffff #f7ffff7f #ffffffff #f7ffffff #fffdfeff #fffffffe', 
            ' #fdffffff #ffffffff #ffefffff #ffffdfff #fffffffe #ffffffff:3 #fffffdff', 
            ' #7fffffff #ffffffff #fdfffdff #fffffeff #ffffffff #efffffff #f77fffdf', 
            ' #fffffffd #efffffff #ffffff7f #fffffc7f #ffffffff:5 #dffffffb #ffffffff', 
            ' #fffffffd #ffffffff:2 #fffeffff #ffffffff #feffffff #ffffffff #fffffffb', 
            ' #7fffffff #ffffffff:2 #ffffffdf #ffffffff #fffff7ff #fbfff7ff #fffffbff', 
            ' #f5ffff7f #efffffff #f7ffffff #ffffffff:2 #ffefffff #ffffffff #ffbeffff', 
            ' #7ffdfdff #fefffffb #ffffffff:3 #ffdfffff #ffff7edf #fdffffff #ffbfffff', 
            ' #f3bfffff #ffffffff #7ffffbff #ffbfcfff #ffffffff #ffff7fff #bfffffff', 
            ' #7fff7fff #ffffffff:2 #fffeffff #ff7fffff #ffffffff #fff7bfff #ffffffff', 
            ' #fffbffff #fffff7ff #fdffff7f #ffffffff:7 #ffdfffff #fdffffff #fffffffe', 
            ' #feffffff #ffffffff #ffdfffdf #fffffeff #dfffe7ff #ffffffff #7ffdffff', 
            ' #fffffbff #ffffffff:3 #feffffff #fffffff3 #fffffebf #ffffffff:4 #feffffff', 
            ' #ffffffef #ffffdfff #fdffff7f #ffffffff #fbffffff #ffffefff #ffffff7f', 
            ' #ffffffdf #ffffffff:2 #fffffbff #feffffff #ffffffff:3 #ffffbfff #ffffffff', 
            ' #ffffff7f #ffffffff #ffffedff #7fffffff #ffffffff:4 #ffffff7f #ffffffff:4', 
            ' #fffffeff #ffffffff:2 #ffbfffff #feffffff #fdff1fff #fffeffff #fffffffb', 
            ' #7f7ffffe #fffffbff #ffddefbf #f7ffffff #ffffffff:2 #ffffefff #ffffffff:2', 
            ' #ffbfffff #fffffffb #ffffffff #ffffffdf #ffffffff #7fbfffef #ffffffff:3', 
            ' #ffbfffff #dfffffcf #ffffffef #ffffffff:3 #fffffdff #ffffff77 #ffffffff', 
            ' #f7fff7ff #ffffffff:2 #ffffffbf #fdffefff #ffffffdf #feffffff #fdbfffdf', 
            ' #ffffffff #ffff7fff #ffffffff #ffdffddf #ffffffff #fffffffb #deffffff', 
            ' #fbffdffb #bfffffff #ffffff7f #ffffffff #ffdff7ef #fdffdfff #ffffffff', 
            ' #dfffffff #dffffff7 #ffffffff #fffff7ff #ffffffff:2 #fefffeff #ffffffff:2', 
            ' #f3ffffff #ffffffff #ff7fffff #fffff7ff #3fffffff #ffffffff #ffbdffff', 
            ' #feffffff #ffffffff:5 #fffffff7 #ffffffff #ffff7fff #ff7fffff #ffffffff', 
            ' #fffffaff #ffffffff #ffdfffff #ffffffff:2 #ffdf7fff #fffffe7b #ffffffff', 
            ' #ffbffbff #fffff7ff #ffffffff:5 #ffffffbf #fffffffb #ffffffff #bdff9fff', 
            ' #ffffffff #bfff1bfe #fffffeff #dffeffff #ffffeeff #ffdfffdf #ffffffff:3', 
            ' #ebffffff #ffffffef #ffffffff #fffdffff #ffffffff #fffffffe #ffffffff:2', 
            ' #dfffffff #ffffffff:2 #ffffbfff #fdf7ffff #ffffffbd #ffffffff #fffffff7', 
            ' #ffffffff:3 #ffffffdf #ffffffff #ffefffff #7bffdffe #ffffffe7 #fffbffff', 
            ' #ffffffff #fff7ffff #ffffff7f #ffefffff #ffffffff #fffffbff #ffffffff:4', 
            ' #7fffffff #ffdfefff #fffbffff #ffffffff #7ffdffff #fbfffeff #5ff7f7ff', 
            ' #ffdfffff #ffffffff #ffefffff #ffffffff #fffefeff #ffffdff7 #f7feffff', 
            ' #ffbfffff #fffffffe #ffffffff:4 #ff7fffef #ffffffff #fdffffff #ffffffff:3', 
            ' #f7ffefff #fffffeef #ff777fff #fff7ffff #ffefffff #fdffffff #fffdffff', 
            ' #ffdfffff #ffffffff #fbffffff #ffffeeef #ffffffff #7fffffff #fff7fdeb', 
            ' #dfffffff #ffdfbfff #ffffffff:2 #ffffffbf #ffffffff #fdfdfefe #ffffffff', 
            ' #dffffffb #ffffffdf #fffffeff #ff7fffff #bfbffff7 #ff7fffff #ffffffff:9', 
            ' #dffeffff #ffffffff #fffbffdf #feffffff #ffffff7f #ffffffff #fffff7ff', 
            ' #ffffffff:49 #ffff7fff #ffffffff:30 #ffff7fff #ffffffff:24 #ffffff7f #ffffffff:40', 
            ' #efffffff #ffffffff:5 #fffff ]', ), ), face3Elements=
            mdb.models['Model-1'].rootAssembly.instances['Tool'].elements.getSequenceFromMask(
            mask=(
            '[#7ffbffff #ffff7fef #ffffffff #fffeffff #fffffbff #efffffff #ffffffff', 
            ' #ffbffffb #ffffdff7 #ff7ffffd #ffdfffff #bfffffef #f7ffffdf #ffffffff', 
            ' #7fffffff #ffffffff #ff7ffdff #ffffffff:3 #fffdffff #ffffffff:2 #feffffff', 
            ' #ffdfffff #ffffefef #fdffffff #ffffffff:2 #fffefeff #ffffffef #ffffffdf', 
            ' #ffffffff #f7ffffff #feffffff #ffffffff:2 #fffffdff #fff7ffff #ffffffff:2', 
            ' #ffdfb7ff #bfbfffbf #ffefefff #bfdfff7f #fdffdfff #dfffffbf #dfeff7fd', 
            ' #ffffffff #fffcffff #ffff7fff #ffffffff #fffebfef #ffffffbf #ffffffff', 
            ' #f7fff7fe #fffeff6d #fbefefff #fffffff6 #ffffffbf #ffffeffd #fffffffe', 
            ' #e5ffffff #7fffbbef #ffffffff:2 #fff7ffff #ffffffff #f7fbffef #efffffff', 
            ' #ffffeffd #fff7ffff #bfffefff #ffffdfff #ffffffff #fffffdf7 #f7ffffff', 
            ' #fffffff7 #bffffff7 #ffdfffff #ffffefff #ffafffff #ffffffff #fdfffbfd', 
            ' #ffefffef #ffffffff:3 #ff7fffff #ffffffff:4 #ffefffff #ffffffff:2 #fffffbff', 
            ' #ffffffff:4 #7fbfffff #fbffffff #ffffffff:3 #feffffff #f7ffffff #ffffffff:7', 
            ' #fffffbff #ffffffff:2 #dfffffff #ffffeeaf #ffe7ff77 #f7fbbbfb #fffefbff', 
            ' #abbfdfab #efffffdd #7bf5fa7f #fffebfff #ff6d73ff #f3f7ffdf #dfff8bfd', 
            ' #ebb2fde7 #7dfeeefd #7bffebfe #fff5fbf9 #f5f4ffe7 #ff6ef3bf #fdbbefff', 
            ' #fc37ff7f #f7ff7fe3 #fcfe7dff #f6f7ff3d #3fafd967 #faefafee #f5ebeb47', 
            ' #fefffff4 #bcd5fffe #f7f97b7d #ffdd87f7 #aee7ee5a #ed3fd7fb #a9ffaddf', 
            ' #df3ffb8f #faff3fff #773dd77f #bfdbe667 #fddf1ffb #fff79e7b #ffbffdeb', 
            ' #fbbbe7e9 #6fdef7fd #ffedbffd #fdb6f7f7 #dffff3ff #7bfafb7d #f3ef5fb7', 
            ' #fefeffcd #fcffe7fa #edeffdda #f9f99fff #fe5ffcbf #9beff2fb #ffbdef83', 
            ' #fef3f887 #cfe71bba #dfdde9da #edd25bc7 #f1ff3bf7 #fcff6dff #bffb2d7d', 
            ' #fff1ff7f #bb7eebd3 #3ffeeff7 #fbf2ffdf #1fafff75 #bf7e7dfd #7ff77af6', 
            ' #df3fadef #ffaf57d7 #b63b97f6 #fffbb73d #ffffae77 #97f5fd6e #efff6577', 
            ' #fdffb5bf #fde7fa5e #e7269fbd #277ecf6e #3b9c19f3 #f6e7ff2f #fffafc7e', 
            ' #5ffd9fcb #d3ebffef #ef6ad3e5 #dffeff2b #79e7d7df #73fcff92 #fef435af', 
            ' #f4ff516f #bf7ff65b #dfef76ff #dfeef9ff #ab7dff4e #f73b7f5a #7cd3dffe', 
            ' #7dd5e0af #efcfd7f8 #9feddcbf #fbff79f7 #e6fabfdf #99fcbf1f #eff7b7cb', 
            ' #f7bbfbdf #7f7df3df #7db7fb7d #f79b9755 #93ede77e #effbcfbf #edf9c7f7', 
            ' #6fbebded #ffaebf6e #3f5f9df5 #cf7fbffe #bfc7df2d #5fff5fd7 #9ce5cdfd', 
            ' #6f6e3cfd #fdf77ddf #b9dffdb7 #f6c7777f #ffd6766e #f5fb6fec #f7f9efcb', 
            ' #fbffbf67 #dffbf9ff #7eeefb73 #fe77ddc3 #7f2ffffc #6fbddaee #fffa7ec6', 
            ' #9befdfff #ffff6bef #ffbffffc #efdf6b7d #bfdebcdd #f38e9dde #fffffbdf', 
            ' #7a9fdffa #efddddbb #4df9bc7f #ff76fc96 #bee1b2dd #fbbfeffa #cfffdfbd', 
            ' #6defbfdf #dbf8fdad #bbf5bfff #ff7f6bfb #7cfdfffd #8efebbbd #237dffa7', 
            ' #f77f5bfe #fffcfeff #fdffff29 #6c7e5eff #67abdf5b #fef5fdbf #ff70ffff', 
            ' #e7ebbf71 #f6ffc6f7 #7bdbffbb #7ffdffed #f7b6be7f #7ffcefbf #d7f4f9c2', 
            ' #7f3dfd7e #effdf5dc #faa9bfff #fe6f7e27 #effff3f8 #fefdfef3 #befedf7f', 
            ' #ff05ffeb #fcfeffe7 #fffb77df #9dffb66f #7e3f9cfd #ffff4efd #f7fffefe', 
            ' #ff7fd7ff #ff99ffff #3ebfbe5f #ff773bf7 #ffddffb3 #ffffff79 #ff7fbbef', 
            ' #f7f3fffe #f6cafeff #dffffb7f #7f7fffec #2ffbbbdd #6557fff6 #ef3b6ddd', 
            ' #bf7e46bd #21ff75d9 #d9a2f9ff #efb77c27 #776df97f #ffd3feac #affffeb6', 
            ' #8b1f7efd #f6fdd6ff #bddfffbe #7f7ffa78 #f0ebfe7f #ffb7fbf4 #d773dfb8', 
            ' #ffedc927 #fbfbffbf #7f7ffb3e #ff7ffffd #fd7dbff3 #fff7ff5f #676f59f7', 
            ' #bf5feefb #ffaee592 #f71ffdfa #af9fbedf #7fffb4ef #bfefdefe #fddbfddf', 
            ' #f77ddbfa #c7ffcfcc #bdf7fde5 #eff673df #de298b8f #77fff78e #bfbdeff9', 
            ' #fbfff5ff #fde7ffee #efbfe5fe #f5e7fbf7 #f7e506f7 #ffd6f7ae #db7fbdff', 
            ' #dfff1b3d #ffbfee7e #659bbb53 #22342cff #ed00ffab #bfb9ffff #7efb676f', 
            ' #ebe9ff9d #96e79947 #7edf7ef6 #ffdaf2f5 #93ff74f3 #d3fce525 #ffb2fff7', 
            ' #fabbffbe #7b7fd37a #fbe4fbd7 #e428b87f #8ffd16f2 #6793bfef #ffff48ff', 
            ' #7ffd77ad #a7575bc7 #dbbe97f7 #cfdd69bf #f078cc7f #9976475e #7ffda95b', 
            ' #ffbe75ae #ffc9fcfb #caedf677 #c77ff7ee #cefeda7f #fffdfbfc #ffdee7bf', 
            ' #fbe9effb #f7fb7ad8 #7fabab9d #edd0fbfc #ffb7a68e #f87fffcf #ffbf7fef', 
            ' #d6dddf7f #9c7edacc #7fdfc7ff #3bbfceff #f7ffefdf #e7ddd15f #aff9f6dd', 
            ' #bff5ff7f #78c7ffff #fdfecece #ef83cb7b #5defb0e9 #af6b7feb #6aaab7cd', 
            ' #ee9f77b7 #7fdd576f #fbb9ddd7 #b99cbfdf #f59fb9e8 #fefedfe5 #ffdadb7f', 
            ' #f7f5fbc3 #afdff6e9 #2efe7f5d #94dfcdad #f7a27ba6 #39edfff9 #2fbcfb7e', 
            ' #fdff6d7e #d3bf94ff #cfcfddff #f5eff1ba #edb6bcdd #d3fbf7fd #beff7b7f', 
            ' #b6ecfffd #77a7fefb #f545d77f #fbeddd1f #defffcb9 #6ffffdf6 #eef9dd5f', 
            ' #77dffb5d #fdf2feee #f7defbff #47b7ffdb #d4b6137e #ffef9fbd #9ffdff5e', 
            ' #bffbbffb #ecbfffdb #d7e3ddd7 #de5bdabb #dfff3d5e #395dffef #f7f9f9fc', 
            ' #ffd9ffb7 #ff75f9af #6b7c7bdd #8bf6d787 #f9bdf3d9 #8a4ffbfa #eb5a75ff', 
            ' #a3fff299 #effbffaf #7b9dfd7d #7bcebfdf #7efeffce #fe747cfe #6fcc891f', 
            ' #ff8f56ff #bffffff0 #afdbcf1f #dffbdf7a #dbefebcd #affffbb4 #ef7dfebf', 
            ' #7dfbd8ce #bebcfd3c #5afde9fe #79affaf5 #f15fedeb #6ffe6fb5 #f3fbd1f8', 
            ' #bdf73fbd #9ffdcff2 #fff4ff7b #3cbf8eef #eb5fefdb #fded7ce7 #df5fc79e', 
            ' #5bffebea #d7f7fd7f #df677ddf #7f5b7efe #7ffeffed #ffffffdf #b70f7f72', 
            ' #337fbbef #fabfbffd #3bbbffff #cfbfdfbb #fd7ed9fe #da75e9c6 #fe7edbff', 
            ' #bf7bf7f2 #edfdea7c #3d76b09f #f32ed46e #e6bffd5d #fd5e7efa #5f7e57aa', 
            ' #ebdff6a6 #f7ff7cfd #a84fd762 #c5cda6f6 #5fff72bf #ffbecfb9 #737f6d57', 
            ' #dbff77fd #bb7df7cb #befd7eeb #afcfff0f #fbffe729 #beeef97e #eff7f1fd', 
            ' #f6e6fbf9 #f77ffc93 #a78effbf #ff7dfed7 #fffa66ff #f9637fdf #ff9fdff7', 
            ' #fffdfdeb #ae6fbb5f #b559e9f6 #f6b8e77f #ffffefdb #3fffbfef #ffeffbcf', 
            ' #efb3dfdf #8ffdfafb #fdffffff #ebbffefb #fddf7fef #ff7dffff #ffffffff', 
            ' #ffefffef #ffff7fff #fefff7ff #eff3eb26 #fbffb6ed #dd7eee7b #236f67f7', 
            ' #ebdfbe2f #6fbdfdfb #db75becf #cdfff9f9 #fd73b3ff #f3ec9d64 #77ffe6ff', 
            ' #dbe7ffff #ff3ffdbc #5fefffbf #3eed7f7b #9d3501ff #c79b6dd #ef9ebf7a', 
            ' #f7f55fd7 #df9fffff #bfaddfef #ff7acb27 #1f7ffff7 #e7fb3713 #d1d7beea', 
            ' #fafff5fd #ffaff7ff #dfffd5ff #adb7dfbf #ff5b713f #277bbc5f #f7fafa7a', 
            ' #7f7eff77 #ffabbbfe #ff62ffb3 #ff7dff7f #ff7f7f7f #ee7fffe7 #7fbea4b3', 
            ' #e6bf7bfd #ffdbf7ff #337bf2f9 #6d7ffe6e #fa7bffff #fb777fff #dc6df7bb', 
            ' #fae8719e #ff9fff7b #ff7fcf4d #baeefbff #f7fffbff #fff7efcf #c3d7dff9', 
            ' #c7ff7fed #dd6dcffb #4fbcffff #33bef45f #df9ffff7 #cdffcffb #f7ef6fff', 
            ' #fffffbff #af3fb777 #cbf787ea #6ff7d5bd #fefebbfe #ff17b7ce #1ffefebe', 
            ' #cfacfddf #8f7efffc #ba7efbfe #f9efbf9f #e2d3ffff #eb7fdeff #8f57f6ff', 
            ' #cefffd47 #5ba6fdff #3bb4bfde #7feff715 #efcff93f #f4fff7ff #ffefefff', 
            ' #f5afffef #4b8cffdf #cfa65fd5 #7fbed973 #77f71aaa #1ba496ef #b375ffe1', 
            ' #fd7bd9bd #cddf1b8e #dd61c4fb #ffeff7df #5bdfdffa #bbadfbfd #55ffdfeb', 
            ' #e5fdfbff #777fdfe6 #fdfdfffe #7fbdd7d7 #fdf4e67e #fffaef7f #69fafeff', 
            ' #fbeefff6 #d9bdbd6e #fc7eff69 #cfffffff #bd5fefff #5af7f9f7 #7f7d7a3e', 
            ' #2dbdff8d #eb4ddcce #efba5c4b #757ec747 #e8f87b7d #fff5fff7 #9ab67ace', 
            ' #5f7edfff #739315c7 #3effb740 #cffff6ff #dffbfddf #fdafbed7 #fefd78bd', 
            ' #dffdf7fc #fe2fead7 #ddfd7b7d #fefe97ad #7dde7f7b #ab9fafff #75d35d5f', 
            ' #bced5fee #fffd7f7f #7efddd6e #e5b6b8cf #ffbfffde #d2bffffe #eff5ccff', 
            ' #ffcdeda6 #6fb7bebf #cbf5f7d5 #73f1fdbc #77fbdbea #7fa5daf9 #7a3fb77b', 
            ' #8dcffd2f #7fffd7f2 #ef7ffb5f #ffc7e8fb #ebff737b #c6fdfbaf #aefb3eff', 
            ' #bbfeeff6 #ebffe7e5 #ea664fd8 #d773fe7 #fffcddf3 #e94cffad #3f7ef6ff', 
            ' #ce7fffdf #6bfbf5ef #ddefdff9 #2e96fbff #ad9a1a76 #937dfb3b #7f6ffc37', 
            ' #d81359ef #fbbbffbd #ef562fdb #afdc5ff9 #e76dfbff #dfbaffbf #7bdfd9fd', 
            ' #a77fff7f #edb9e577 #fbeffbf5 #fd7fff5f #dfbd0fef #fff9d3ff #7fffff7f', 
            ' #bfff7b5f #57deefdf #ffffffcf #edfffe79 #7feefdff #feffdffb #f7fb7eff', 
            ' #b74fffff #ffffffbf #7fffdfef #fdffffff #1ffe7fff #b9effffd #9fffdff7', 
            ' #3bffd1fd #ffffd7ff #b7b79fff #ffbcfdff #dffdffff #ffff7fff #efefffff', 
            ' #fff5f8ff #ffef7fcf #9ebfffff #8ffbfffb #57fffbff #fffdd7fb #ff7f7fff', 
            ' #ffffffb6 #fffdffff #f3fffdff #9fdfffff #ffffffff #f7ffffff #f7fdafff', 
            ' #fffffecf #dfdfffdf #ffffffff #efffffef #ffffffff #bfff7fff #ffffffdf:2', 
            ' #fffffff4 #ffffffdf #fff7fbbf #fffffffb #3efffff1 #fffffeff #fffeffff', 
            ' #ffffffff #fffbffff #ffffbfff #ffffffff #ffffefff #ffffff7f #f7ffffff', 
            ' #ffff6dff #fdffffff #ff7fdfff #ff7ffffc #afffeffb #ffffffff:2 #fffbffef', 
            ' #7ffbffff #ffffbfff #fffff7fd #e5bbffff #7fffffff #ffff9ebd #f3ffffff', 
            ' #fffffeef #efffffff #ffff7a7f #ffffffff:2 #fff58dff #ffffd7ff #fddfffff', 
            ' #ffffffff #ffffff7d #95ffbfff #ffffffef #f77fffff #ffffffff #bfffdeef', 
            ' #ff6fbffb #ffffffbf #bfbfffff #ffdf3ef7 #ffffffff #7fff7fff #ffffffff:2', 
            ' #7fffff7f #ffffffe7 #ffffffff #ffffffdb #ffffd7ff #3ff7ffff #ffffffff', 
            ' #fdffefff #efffffff #fffdffef #ffffffff #feffffff #ffffffff:2 #ffff7fff', 
            ' #fffff7bf #ffffffff:2 #ffdfffef #bfffffff #ffdfffff #fffefdff #ffffffff:2', 
            ' #f7ffffff #ffffffff:4 #ffdf7fff #27ffffff #efffefff #efffffff #ffffffff', 
            ' #ffffebfb #ffffffff:2 #fdfeb7ff #ffffffff #cffeffff #ffffffff #deffff7f', 
            ' #fffff ]', ), ), face4Elements=
            mdb.models['Model-1'].rootAssembly.instances['Tool'].elements.getSequenceFromMask(
            mask=(
            '[#ffffffff #ffbfffbf #ffffffff:3 #fffffdff #ffffffff:5 #bfffffff #3fffffff', 
            ' #ffffffff #ffffffbf #ffffffff:4 #ffffffbf #ffffffff:5 #bfffeffd #ffffdfff', 
            ' #ffffffff #fffffffe #ffffffff:6 #efffffff #ffffffff:2 #ffffffef #ffffffff:2', 
            ' #7effffff #ffffffbf #fffeffff #ffffffff:11 #fffffeff #ffffffff:10 #ffff7fff', 
            ' #ffffffff:40 #feffffff #ffffffff:5 #fffffdff #ffffffff:9 #ffbfefff #ffffffff:12', 
            ' #ffffdfff #ffffffff:4 #ffffffbf #7fffffbf #ffffffff:5 #fffffffb #ffefffff', 
            ' #ffffffff:7 #fbffffff #ffffffff:3 #fff7ffff #ffffffff #fffeffff #ffffff7f', 
            ' #ffffffff #efffbfff #ffffffff:6 #ffffffef #ffffffff:5 #ffff7fff #ffffffff:4', 
            ' #ffffffdf #ffffffff #fdffffbf #ffffffff #dfffffff #ffffdfff #ffffffff:7', 
            ' #ffffdfff #ffffffff:4 #bfffffff #ffffffff:2 #efffffff #ffffffff #ffeffff5', 
            ' #fffbffff:2 #ffffffff #fffffffd #ffffffff:3 #fff7ffff #ffffffff:3 #fefeffff', 
            ' #fffeffff #fefffdff #ffffffff #f7ffffff #ffffffff:4 #ffff7fff #ffffffff:2', 
            ' #ffefffff #ffffffff:4 #ffffffbf #ffffefff #ffffffff:2 #ffdfffff #ffffffff:6', 
            ' #f7fffffb #ffffffff #fffeffff #fffffdff #ffffffff:5 #ffffbfdf #f7ffffff', 
            ' #ffffffff:2 #ffefffff #fbffffff #ffffdbff #ffffffff #dfffbfff #fffffffd', 
            ' #ffffffff #fffdffff #fffeffff #ffdfffff #ffffffff #dfffffff #ffffffff', 
            ' #fbfffff7 #ffffffff #effffeff #ffefffff #ffffffff #ffffff7f #ffffffff:3', 
            ' #fffbffff #ffffffff:3 #fefbffff #ffffffbf #fffffbff #ffffffff:4 #f7fff7ff', 
            ' #ffffffff #fffffdff #fffffff6 #efffffff #fbfbffff #ffffddff #ffffffef', 
            ' #fdfff7f7 #ffffffff #ffffffdf #fffffffd #fffeffff #ffdfffff #ffff7bfd', 
            ' #ffffffff #fffffef7 #ffffffff #ffffffdf #ffffbfff #fffffff7 #7fffffff', 
            ' #ffffffff:2 #fffffffd #fffffffb #f7fbffcf #ffffffff #fffffffd #fdffffbf', 
            ' #ffffffff #ffffffbf #ffffffff:3 #ffffffbf #ffffffff:4 #ff7fffff #ffffffff:2', 
            ' #ffffefff #ffffffff:2 #fffbefff #ffffffff #fffffbff #ffffffff:5 #f7fffff7', 
            ' #1fffffff #ffffffff #ffffdfff #feefffff #ffffffff:2 #ffbffef7 #ffefffff', 
            ' #ffffffff #fffeffff #ffffffff #ffffffbf #ffffffff:6 #fffffffe #7fffffff', 
            ' #feefffff #ffffffff #f5ffffff #fffffff7 #ffffffff:2 #fffffdff:2 #ffffffff', 
            ' #fbffffff #ffffffff #ffffffdf #feffffff #fffffeff #ffffffff:6 #feffffff', 
            ' #ffffffff:3 #fffffffb #ffffffff #fffffffb #ffffffff:5 #efffdff7 #efffffff', 
            ' #ffbfffff #ffffffff #fffffbff #ff7fffff #ffffffff:2 #fbffffff #fffffeff', 
            ' #ffffffff:2 #fff7ffbf #ffffefff #ffffffff:2 #fffffffe #ffffffff #ffbfffff', 
            ' #ffffffff #ff5fffff #ffffffff:2 #ffffbfff #ffffffff:3 #ff7fffff #ffff7fff', 
            ' #ffffffff:2 #ffdfffff #fffffbff #ffffffff:3 #ffdfffff:2 #ffffffff:4 #fffffeff', 
            ' #ffffffff:2 #efffffff #ffffffff:2 #fffbf7ff #ffffffff:4 #ff7fffff #ffffffff:3', 
            ' #fbffffff #ffffffff:3 #fffffff7 #fdffffff #f7ffffff #ffffffbf #fffefffb', 
            ' #fffffff7 #ffffe7ff #ffffffff:2 #efffffff #ffffffff #ffffffdf #ffffffff:2', 
            ' #ffffffef #ffffffff:2 #fdfdffff #ffffffff:5 #ffffffef #ffffffff:7 #fffffeff', 
            ' #ffffffff #dfffffff #fffdffff #fd7ffeff #fffffffb #ffffffff #7fffffff', 
            ' #ffffffff:2 #fbffffff #ffeffeff #ffffffff #bfffffff #ff7fffff #ffffffff', 
            ' #ffff7fff #ff7ffff7 #ffffffff #fffffffe #bfbfffff #ffffbfbf #ffffffff:4', 
            ' #fdffffff #ffffffff:2 #fdffffef #fffffbff #ffffffff #fffdffff #ffffffff:3', 
            ' #ffffffcf #ffffffff:2 #fffbffef #ffffffff #ffffefff #ffffffff #fffffffb', 
            ' #efffffff #fdffffff #ffffffff:2 #dfffffff #effff7bf #fbffffff #ffffffff', 
            ' #ffff3fff #ffffff77 #ffffffff:3 #ddfffeff #fffffffb #ffffffff:3 #fffffffe', 
            ' #effbffff #fffffdff #ffffffff #fffffffb #ffff3fff #ffffffff #fffb7fff', 
            ' #ffffffff #ffffffef #ffffffff #fffdffff #ffffffff:8 #ffffffef #ffffffff:4', 
            ' #f7ffffff #ffffffff #ffbfffff #ffffffff:2 #fdffffff #ffffffff #fffffffd', 
            ' #ffffffff:2 #fffff7ff #fffdffff #ffffffdf #ffffefff #ffffffff #fffff7ff', 
            ' #ffffffff #fffffbff #ffffffff:2 #ffbfffff #ffffffff #f7febfff #ffffffff:4', 
            ' #f7ffbfff #ffffffff:4 #fffefffb #ffffffff:9 #feffffff #ffffffff:15 #f7fbffff', 
            ' #ffffffff #ffffdfff #ffffffff #7bfff7ff #efbfffff #ffffffff:3 #e7dfffff', 
            ' #ffdfffff #ffffffff:2 #fbffffff #ffffefff #ffffffff:5 #ffdfffff #ffffffff', 
            ' #fffffffb #ffffffff #ffffff7f #ffff7fe7 #fffffeff #fff7ffff #bfffffff', 
            ' #ffdfffff #ffefffff #ff7fffff #ffffffff #dfffffff #ffffffff #fffdffff', 
            ' #ffffffff #fffbffff #ffffffff #fbffffff #ffffffff #fffbffff #ffffffff', 
            ' #ffbfffff #ffffffff:7 #cfffffff #fffffbff #ffffdfff #fff7ffff #fffeffff', 
            ' #ffffffff:3 #ffffefff #ffffdfff #ffffffff:3 #ffffffdf #ffffffff #ffefffff', 
            ' #ffffffff:4 #fffdffff #ffffffff:11 #ffffff7f #feffffff #ffffffff:3 #fbffffff', 
            ' #ffffffff #fffbffff #ffffffff:7 #fffffbff #ffffffff:15 #ffffffdf #ffffffff', 
            ' #ffff7fff #ffffffff:4 #efffffff #ffffffff:2 #ffefffff #ffffffff:22 #ffff7fff', 
            ' #ffffffff:8 #fffeffff #ffffffff:4 #ffffff7f #ffffffff:7 #fffffdff #fffffbff', 
            ' #ffffffff #fffffffe #fffff7ff #ffffffff:3 #ffffdfff #ffffffff:22 #ffff7fff', 
            ' #ffffffff:2 #ffefffff #ffffffff:15 #efffffff #ffffffff:2 #feffffef #ffffffff:2', 
            ' #fffff ]', ), ), name='Tool Surface')
        

    def workpiece_bc(self):
        name = 'Workpiece Bottom Encastre'
        a = mdb.models['Model-1'].rootAssembly
        f1 = a.instances['Workpiece'].faces
        faces1 = f1.getSequenceFromMask(mask=('[#4082010 ]', ), )
        region = a.Set(faces=faces1, name='Set-1')
        mdb.models['Model-1'].EncastreBC(name=name, createStepName='Initial', 
            region=region, localCsys=None)

    def tool_bc(self, step):
        name = 'Tool RP Velocity'
        model = mdb.models['Model-1']
        # create reference point
        refPoint = self.a.ReferencePoint(point=self.toolInstance.InterestingPoint(edge=self.toolInstance.edges[170], rule=CENTER))
        refPointRegion=regionToolset.Region(referencePoints=(self.a.referencePoints[refPoint.id],))
        toolBodyRegion = self.a.Set(cells=self.toolInstance.cells, name=name)
        model.RigidBody(name='Tool RP RigidBody', refPointRegion=refPointRegion, bodyRegion=toolBodyRegion)
        model.rootAssembly.features.changeKey(fromName='RP-1', toName='Tool Top RP')
        # create velocity BC
        model.VelocityBC(name=name, createStepName='Initial', 
            region=refPointRegion, v1=0.0, v2=0.0, v3=0.0, vr1=0.0, vr2=0.0, vr3=0.0, 
            amplitude=UNSET, localCsys=None, distributionType=UNIFORM, fieldName='')
        model.boundaryConditions[name].setValuesInStep(stepName=step.name, v3=30)
        model.boundaryConditions[name].setValuesInStep(stepName=step.name, vr2=-6000)        


if __name__ == "__main__":
    executeOnCaeStartup()
    Mdb()  # clear all
    model = mdb.models['Model-1']
    

    tool = Tool("Tool",  scale=1e-3)
    tool.create()
    carbide = Carbide()
    tool.set_material(carbide)
    tool.mesh()

    workpiece = Workpiece("Workpiece", length=mm(50), w_height=mm(40), w_width=mm(8), b_height=mm(10), b_width=mm(25))
    workpiece.create()
    ti6alv = Ti6AlV()
    mdb.models['Model-1'].materials['Ti6AlV'].johnsonCookDamageInitiation.DamageEvolution(type=DISPLACEMENT, table=((0.002, ), ))

    alu = Material("Alu", density= 2700, young=70e9, poisson=0.33, 
        A=3.241e8, B=1.138e8, n=0.42, 
        d1=-0.77, d2=1.45, d3=-0.47, ref_strain_rate=1.0, disp_at_failure=1e-4)

    mdb.models['Model-1'].materials["Alu"].Conductivity(table=((7.2, ), ))
    mdb.models['Model-1'].materials["Alu"].SpecificHeat(table=((560.0, ), ))

    workpiece.set_material(alu)
    workpiece.mesh()

    step1 = Step("Step-1")
    
    assembly = Assembly(workpiece, tool)
    assembly.workpiece_bc()
    assembly.tool_bc(step1)
    Interaction(step1)    
    
    model.fieldOutputRequests['F-Output-1'].setValues(variables=('A', 'CSTRESS', 'EVF', 'LE', 'PE', 'PEEQ', 'PEEQVAVG', 'PEVAVG', 'RF', 'S', 
        'STATUS', # important!
        'SVAVG', 'U', 'V'))
    model.fieldOutputRequests['F-Output-1'].setValues(
    numIntervals=1000)
    model.steps['Step-1'].setValues(timePeriod=0.002)



    a = mdb.models['Model-1'].rootAssembly
    session.viewports['Viewport: 1'].setValues(displayedObject=a)

    mdb.Job(name='Job-1', model='Model-1', description='', type=ANALYSIS, 
    atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, 
    memoryUnits=PERCENTAGE, explicitPrecision=SINGLE, 
    nodalOutputPrecision=SINGLE, echoPrint=OFF, modelPrint=OFF, 
    contactPrint=OFF, historyPrint=OFF, userSubroutine='', scratch='', 
    resultsFormat=ODB, parallelizationMethodExplicit=DOMAIN, numDomains=4, 
    activateLoadBalancing=False, multiprocessingMode=DEFAULT, numCpus=4)
    # mdb.jobs['Job-1'].submit(consistencyChecking=OFF)