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
# -*- coding: mbcs -*-
from part import *
from material import *
from section import *
from assembly import *
from step import *
from interaction import *
from load import *
from mesh import *
from optimization import *
from job import *
from sketch import *
from visualization import *
from connectorBehavior import *


import inspect
filename = inspect.getframeinfo(inspect.currentframe()).filename
cd = os.path.dirname(os.path.abspath(filename))

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


class Workpiece:
    def __init__(self, name, length, width, height, t=mm(1)):
        self.name = name
        self.width = width
        self.length = length
        self.height = height
        self.part = None
        self.t = t

    def create(self):
        s = mdb.models['Model-1'].ConstrainedSketch(name=self.name + '_sketch', sheetSize=max(self.width, self.height, self.length))
        g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints    
        s.sketchOptions.setValues(decimalPlaces=4)
        s.setPrimaryObject(option=STANDALONE)
        s.rectangle(point1=(0, 0), point2=(self.length, self.height))
        p = mdb.models['Model-1'].Part(name=self.name, dimensionality=THREE_D, type=DEFORMABLE_BODY)
        self.part = mdb.models['Model-1'].parts[self.name]
        self.part.BaseSolidExtrude(sketch=s, depth=self.width)
        s.unsetPrimaryObject()
        # create 
        c = self.part.cells
        pickedCells = c.getSequenceFromMask(mask=('[#1 ]', ), )
        h = self.height - self.t
        workpiece.part.PartitionCellByPlaneThreePoints(cells=pickedCells, point1=(mm(0), h, mm(0)), point2=(mm(1), h, mm(0)), point3=(mm(1), h, mm(1)))

    def set_section(self, section):
        c = self.part.cells
        cells = c.getSequenceFromMask(mask=('[#1 ]', ), )
        region = self.part.Set(cells=cells, name='Material region top')
        self.part.SectionAssignment(region=region, sectionName=section.name, offset=0.0, 
            offsetType=MIDDLE_SURFACE, offsetField='', thicknessAssignment=FROM_SECTION)
        cells = c.getSequenceFromMask(mask=('[#2 ]', ), )
        region = self.part.Set(cells=cells, name='Material region bottom')
        self.part.SectionAssignment(region=region, sectionName=section.name, offset=0.0, 
            offsetType=MIDDLE_SURFACE, offsetField='', thicknessAssignment=FROM_SECTION)

    def mesh(self):
        e = self.part.edges
        pickedEdges = e.getSequenceFromMask(mask=('[#9050 ]', ), )
        self.part.seedEdgeByNumber(edges=pickedEdges, number=20, constraint=FINER)
        pickedEdges = e.getSequenceFromMask(mask=('[#1002a ]', ), )
        self.part.seedEdgeByNumber(edges=pickedEdges, number=200, constraint=FINER)
        self.part.seedPart(size=0.002, deviationFactor=0.1, minSizeFactor=0.1)

        self.part.seedPart(size=0.002, deviationFactor=0.1, minSizeFactor=0.1)
        pickedEdges = e.getSequenceFromMask(mask=('[#2000 ]', ), )
        
        c = self.part.cells
        pickedRegions = c.getSequenceFromMask(mask=('[#2 ]', ), )
        self.part.setMeshControls(regions=pickedRegions, technique=SWEEP, algorithm=ADVANCING_FRONT)
        elemType1 = mesh.ElemType(elemCode=C3D8R, elemLibrary=EXPLICIT, 
            kinematicSplit=AVERAGE_STRAIN, secondOrderAccuracy=OFF, 
            hourglassControl=DEFAULT, distortionControl=DEFAULT, elemDeletion=ON)
        elemType2 = mesh.ElemType(elemCode=C3D6, elemLibrary=EXPLICIT)
        elemType3 = mesh.ElemType(elemCode=C3D4, elemLibrary=EXPLICIT)
        cells = c.getSequenceFromMask(mask=('[#1 ]', ), )
        pickedRegions =(cells, )
        self.part.setElementType(regions=pickedRegions, elemTypes=(elemType1, elemType2, elemType3))
        self.part.generateMesh()

        

class Cutter:
    def __init__(self, name, length, width, height, angle):
        self.name = name
        self.width = width
        self.length = length
        self.height = height
        self.part = None
        self.angle = angle

    def create(self):
        s = mdb.models['Model-1'].ConstrainedSketch(name=self.name + '_sketch', sheetSize=max(self.width, self.height, self.length))
        g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints    
        s.sketchOptions.setValues(decimalPlaces=4)
        s.setPrimaryObject(option=STANDALONE)
        s.rectangle(point1=(0, 0), point2=(self.length, self.height))
        p = mdb.models['Model-1'].Part(name=self.name, dimensionality=THREE_D, type=DEFORMABLE_BODY)
        p = mdb.models['Model-1'].parts[self.name]
        p.BaseSolidExtrude(sketch=s, depth=self.width)
        s.unsetPrimaryObject()
        self.part = mdb.models['Model-1'].parts[self.name]

    def extrude_cut(self):
        f, e = self.part.faces[3], self.part.edges[8]
        t = self.part.MakeSketchTransform(sketchPlane=f, sketchUpEdge=e, sketchPlaneSide=SIDE1, sketchOrientation=RIGHT, origin=(0, 0, 0))
        # draw sketch
        s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=0.01, gridSpacing=0.0005, transform=t)
        max_l = max(self.width, self.length, self.height)
        p1 = (0, 0)
        p2 = (max_l, max_l * tan(self.angle))
        s.Line(point1=p1, point2=p2)
        p3 = (p2[0], p2[1] - max_l)
        s.Line(point1=p2, point2=p3)
        p4 = (0, p3[1])
        s.Line(point1=p3, point2=p4)
        s.Line(point1=p4, point2=p1)
        self.part.CutExtrude(sketchPlane=f, sketchUpEdge=e, sketchPlaneSide=SIDE1, 
            sketchOrientation=RIGHT, sketch=s, flipExtrudeDirection=OFF, draftAngle=37.0)
        s.unsetPrimaryObject()
        del mdb.models['Model-1'].sketches['__profile__']

    def set_section(self, section):
        c = self.part.cells
        cells = c.getSequenceFromMask(mask=('[#1 ]', ), )
        region = self.part.Set(cells=cells, name='Material region')
        self.part.SectionAssignment(region=region, sectionName=section.name, offset=0.0, 
            offsetType=MIDDLE_SURFACE, offsetField='', thicknessAssignment=FROM_SECTION)

    def mesh(self):
        self.part.seedPart(size=0.0005, deviationFactor=0.1, minSizeFactor=0.1)
        self.part.generateMesh()


    def __str__(self):
        return "Cutter: name={}, length={}, height={}, width={}, angle={}".format(self.name,self.length,self.height,self.width, self.angle)

    def __repr__(self):
        return str(self)

class Assembly:

    def __init__(self, workpiece, cutter):
        self.workpiece = workpiece
        self.cutter = cutter
        self.a = mdb.models['Model-1'].rootAssembly
        self.a.DatumCsysByDefault(CARTESIAN)
        self.a.Instance(name=workpiece.name, part=workpiece.part, dependent=ON)
        self.a.Instance(name=cutter.name, part=cutter.part, dependent=ON)
        self.a.translate(instanceList=('Tool', ), vector=(0.001754, 0.015977, 0.0))
        self.a.translate(instanceList=('Tool', ), vector=(-0.01, 0.0, 0.0))
        self.a.instances['Tool'].translateTo(clearance=0.0, direction=(1.0, 0.0, 0.0), fixedList=(mdb.models['Model-1'].rootAssembly.instances['Brick'].faces[5], ), movableList=(mdb.models['Model-1'].rootAssembly.instances['Tool'].faces[5], ))
        
        self.a.ReferencePoint(point=
            mdb.models['Model-1'].rootAssembly.instances['Tool'].InterestingPoint(
            mdb.models['Model-1'].rootAssembly.instances['Tool'].edges[170], CENTER))

        self.a.Set(cells=
            self.a.instances['Tool'].cells.getSequenceFromMask(('[#1 ]', ), ), name='b_Set-2')
        mdb.models['Model-1'].RigidBody(bodyRegion=
            self.a.sets['b_Set-2'], name='Constraint-1', 
            refPointRegion=Region(referencePoints=(self.a.referencePoints[6], )))

        mdb.models['Model-1'].rootAssembly.Set(name='m_Set-3', referencePoints=(
            mdb.models['Model-1'].rootAssembly.referencePoints[6], ))
        mdb.models['Model-1'].rootAssembly.Surface(name='s_Surf-1', side1Faces=
            mdb.models['Model-1'].rootAssembly.instances['Tool'].faces.getSequenceFromMask(
            ('[#ffffffff:2 #3ff ]', ), ))
        mdb.models['Model-1'].Coupling(controlPoint=
            mdb.models['Model-1'].rootAssembly.sets['m_Set-3'], couplingType=KINEMATIC, 
            influenceRadius=WHOLE_SURFACE, localCsys=None, name='Constraint-2', 
            surface=mdb.models['Model-1'].rootAssembly.surfaces['s_Surf-1'], u1=ON, u2=
            ON, u3=ON, ur1=ON, ur2=ON, ur3=ON)
        session.viewports['Viewport: 1'].setValues(displayedObject=self.a)
        session.viewports['Viewport: 1'].view.fitView()
    
    def workpiece_bc(self):
        f1 = self.a.instances[self.workpiece.name].faces
        faces1 = f1.getSequenceFromMask(mask=('[#100 ]', ), )
        region = self.a.Set(faces=faces1, name='Workpiece bottom')
        mdb.models['Model-1'].EncastreBC(name='BC workpiece bottom', createStepName='Initial', region=region, localCsys=None)

    def cutter_bc(self, step):
        mdb.models['Model-1'].rootAssembly.Set(name='Set-5', referencePoints=(
            mdb.models['Model-1'].rootAssembly.referencePoints[6], ))
        mdb.models['Model-1'].VelocityBC(amplitude=UNSET, createStepName='Initial', 
            distributionType=UNIFORM, fieldName='', localCsys=None, name='BC-2', 
            region=mdb.models['Model-1'].rootAssembly.sets['Set-5'], v1=0.0, v2=UNSET, 
            v3=UNSET, vr1=UNSET, vr2=UNSET, vr3=0.0)
        mdb.models['Model-1'].boundaryConditions['BC-2'].setValuesInStep(stepName=
            'Step-1', v1=40.0, vr3=6000.0)

class Step:

    def __init__(self, name, previous='Initial'):
        self.name = name
        mdb.models['Model-1'].ExplicitDynamicsStep(name=name, previous=previous, description='', timePeriod=MAX_TIME)

class Interaction:

    def __init__(self, step, friction):
        # create interaction propery
        model = mdb.models['Model-1']
        model.ContactProperty('InteractionProperty-1')
        model.interactionProperties['InteractionProperty-1'].TangentialBehavior(
            formulation=PENALTY, directionality=ISOTROPIC, slipRateDependency=OFF, 
            pressureDependency=OFF, temperatureDependency=OFF, dependencies=0, table=((friction, ), ), 
            shearStressLimit=None, maximumElasticSlip=FRACTION, 
            fraction=0.005, elasticSlipStiffness=None)
        # create interaction
        model.ContactExp(name='Interaction', createStepName=step.name)
        model.interactions['Interaction'].includedPairs.setValuesInStep(
            stepName=step.name, useAllstar=ON)
        model.interactions['Interaction'].contactPropertyAssignments.appendInStep(
            stepName=step.name, assignments=((GLOBAL, SELF, 'InteractionProperty-1'), ))


class Tool:

    def __init__(self, name, scale):
        self.name = name
        self.step_file = cd + '\milling_cutter.stp'
        self.diameter = mm(10)  # from file
        self.length = mm(100)
        self.cutter_length = mm(22)
        self.scale = scale
        self.part = None

    def create(self):
        print(self.step_file)
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
        self.part.seedPart(size=0.0015, deviationFactor=0.1, minSizeFactor=0.1)
        c = self.part.cells
        pickedRegions = c.getSequenceFromMask(mask=('[#1 ]', ), )
        self.part.setMeshControls(regions=pickedRegions, elemShape=TET, technique=FREE)
        elemType1 = mesh.ElemType(elemCode=C3D8, elemLibrary=EXPLICIT)
        elemType2 = mesh.ElemType(elemCode=C3D6, elemLibrary=EXPLICIT)
        elemType3 = mesh.ElemType(elemCode=C3D4, elemLibrary=EXPLICIT, 
            secondOrderAccuracy=OFF, distortionControl=DEFAULT)
        cells = c.getSequenceFromMask(mask=('[#1 ]', ), )
        pickedRegions =(cells, )
        self.part.setElementType(regions=pickedRegions, elemTypes=(elemType1, elemType2, elemType3))
        self.part.generateMesh()

    def __str__(self):
        return "Tool: name={}, scale={}".format(self.name, self.scale)

    def __repr__(self):
        return str(self)


if __name__ == "__main__":
    executeOnCaeStartup()
    Mdb()  # clear all
    workpiece = Workpiece("Brick", length=mm(20), width=mm(2), height=mm(4))
    workpiece.create()
    workpiece.mesh()

    tool = Tool("Tool",  scale=0.4e-3)
    tool.create()
    carbide = Carbide()
    tool.set_material(carbide)
    tool.mesh()

    # Create materials
    steel = Material("Steel", density= 7870, young=2e11, poisson=0.29, 
        A=375e6, B=552e6, n=0.457, 
        d1=0.25, d2=4.38, d3=2.68, ref_strain_rate=1.0, disp_at_failure=0.1)
    steel_section = mdb.models['Model-1'].HomogeneousSolidSection(name='Steel section', material=steel.name, thickness=None)
    alu = Material("Alu", density= 2700, young=70e9, poisson=0.33, 
        A=3.241e8, B=1.138e8, n=0.42, 
        d1=-0.77, d2=1.45, d3=-0.47, ref_strain_rate=1.0, disp_at_failure=1e-4)
    alu_section = mdb.models['Model-1'].HomogeneousSolidSection(name='Alu section', material=alu.name, thickness=None)
    # assign materials
    workpiece.set_section(alu_section)
    tool.set_material(steel)

    assembly = Assembly(workpiece, tool)
    step1 = Step("Step-1")
    Interaction(step1, friction=0.15)
    assembly.workpiece_bc()
    assembly.cutter_bc(step1)
    mdb.Job(name='Chipformation', model='Model-1', description='', type=ANALYSIS, 
        atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, 
        memoryUnits=PERCENTAGE, explicitPrecision=SINGLE, 
        nodalOutputPrecision=SINGLE, echoPrint=OFF, modelPrint=OFF, 
        contactPrint=OFF, historyPrint=OFF, userSubroutine='', scratch='', 
        resultsFormat=ODB, parallelizationMethodExplicit=DOMAIN, numDomains=4, 
        activateLoadBalancing=False, multiprocessingMode=DEFAULT, numCpus=4)
    mdb.models['Model-1'].fieldOutputRequests['F-Output-1'].setValues(variables=(
    'S', 'SVAVG', 'PE', 'PEVAVG', 'PEEQ', 'PEEQVAVG', 'LE', 'U', 'V', 'A', 
    'RF', 'CSTRESS', 'EVF', 'STATUS'))
    mdb.models['Model-1'].fieldOutputRequests['F-Output-1'].setValues(numIntervals=400)

    # mdb.jobs['Chipformation'].submit(consistencyChecking=OFF)

