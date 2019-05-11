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
        pickedEdges = e.getSequenceFromMask(mask=('[#2000 ]', ), )
        self.part.seedEdgeByNumber(edges=pickedEdges, number=3, constraint=FINER)
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
        self.a.rotate(instanceList=(cutter.name, ), axisPoint=(0, 0, 0), axisDirection=(0.0, 1, 0.0), angle=90.0)
        self.a.translate(instanceList=(cutter.name, ), vector=(workpiece.length, workpiece.height - 0.5 * workpiece.t, workpiece.width))
        self.a.translate(instanceList=('Cutter', ), vector=(0.0, 0.0, 0.001))
        session.viewports['Viewport: 1'].setValues(displayedObject=self.a)
        session.viewports['Viewport: 1'].view.fitView()
    
    def workpiece_bc(self):
        f1 = self.a.instances[self.workpiece.name].faces
        faces1 = f1.getSequenceFromMask(mask=('[#100 ]', ), )
        region = self.a.Set(faces=faces1, name='Workpiece bottom')
        mdb.models['Model-1'].EncastreBC(name='BC workpiece bottom', createStepName='Initial', region=region, localCsys=None)

    def cutter_bc(self, step):
        name = "BC cutter"
        cutter_inst = self.a.instances['Cutter']
        model = mdb.models['Model-1']
        c1 = cutter_inst.cells
        cells1 = c1.getSequenceFromMask(mask=('[#1 ]', ), )
        f1 = cutter_inst.faces
        faces1 = f1.getSequenceFromMask(mask=('[#3b ]', ), )
        e1 = cutter_inst.edges
        edges1 = e1.getSequenceFromMask(mask=('[#a65 ]', ), )
        v1 = cutter_inst.vertices
        verts1 = v1.getSequenceFromMask(mask=('[#20 ]', ), )
        region = self.a.Set(vertices=verts1, edges=edges1, faces=faces1, cells=cells1, name='Cutter set')
        
        model.DisplacementBC(name=name, createStepName='Initial', 
            region=region, u1=SET, u2=SET, u3=SET, ur1=SET, ur2=SET, ur3=SET, 
            amplitude=UNSET, distributionType=UNIFORM, fieldName='', localCsys=None)
        model.TabularAmplitude(name='Amplitude-1', timeSpan=STEP, 
            smooth=SOLVER_DEFAULT, data=((0.0, 0.0), (MAX_TIME, 1.0)))
        model.boundaryConditions[name].setValuesInStep(stepName=step.name, u1=-self.workpiece.length, amplitude='Amplitude-1')


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


if __name__ == "__main__":
    executeOnCaeStartup()
    Mdb()  # clear all
    workpiece = Workpiece("Brick", length=mm(20), width=mm(2), height=mm(4))
    workpiece.create()
    workpiece.mesh()
    cutter = Cutter("Cutter",  length=mm(4),width=mm(4), height=mm(4), angle=deg(10))
    cutter.create()
    cutter.extrude_cut()
    cutter.mesh()
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
    cutter.set_section(steel_section)

    assembly = Assembly(workpiece, cutter)
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
    # mdb.jobs['Chipformation'].submit(consistencyChecking=OFF)

