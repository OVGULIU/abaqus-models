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

    def __init__(self, name, inner_d, outer_d, length):
        self.name = name
        self.length = length
        self.inner_d = inner_d
        self.outer_d = outer_d
        self.part = None

    def create(self):
        s = mdb.models['Model-1'].ConstrainedSketch(name=self.name + '_sketch', sheetSize=max(self.length, self.outer_d, self.outer_d))
        g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints    
        s.sketchOptions.setValues(decimalPlaces=4)
        s.setPrimaryObject(option=STANDALONE)
        s.CircleByCenterPerimeter(center=(0.0, 0.0), point1=(0, self.inner_d/2))
        s.CircleByCenterPerimeter(center=(0.0, 0.0), point1=(0, self.outer_d/2))
        self.part = mdb.models['Model-1'].Part(name=self.name, dimensionality=THREE_D, type=DEFORMABLE_BODY)
        self.part.BaseSolidExtrude(sketch=s, depth=self.length)

    def set_section(self, section):
        c = self.part.cells
        region = self.part.Set(cells=c, name='Workpiece material region')
        self.part.SectionAssignment(region=region, sectionName=section.name, offset=0.0, 
            offsetType=MIDDLE_SURFACE, offsetField='', thicknessAssignment=FROM_SECTION)        

    def mesh(self):
        # seed edges
        e = self.part.edges
        e1_point = (self.outer_d/2, 0, 0)
        e2_point = (self.inner_d/2, 0, 0)
        e3_point = (self.outer_d/2, 0, self.length)
        e4_point = (self.inner_d/2, 0, self.length)
        pickedEdges = e.findAt((e1_point, ), (e2_point, ), (e3_point, ), (e4_point, ))
        self.part.seedEdgeByNumber(edges=pickedEdges, number=90, constraint=FINER)
        # set mesh controls
        c = self.part.cells
        pickedRegions = self.part.Set(cells=c, name='Mesh controls region')
        self.part.setMeshControls(regions=c, technique=SWEEP, algorithm=ADVANCING_FRONT)
        # pick element type and generate mesh
        elemType1 = mesh.ElemType(elemCode=C3D8R, elemLibrary=EXPLICIT, 
            kinematicSplit=AVERAGE_STRAIN, secondOrderAccuracy=OFF, 
            hourglassControl=DEFAULT, distortionControl=DEFAULT, elemDeletion=ON)
        elemType2 = mesh.ElemType(elemCode=C3D6, elemLibrary=EXPLICIT)
        elemType3 = mesh.ElemType(elemCode=C3D4, elemLibrary=EXPLICIT)

        self.part.setElementType(regions=pickedRegions, elemTypes=(elemType1, elemType2, elemType3))
        self.part.generateMesh()

class Tool:

    def __init__(self, name, axis1, axis2, thickness):
        self.name=name
        self.axis1=axis1
        self.axis2=axis2
        self.thickness=thickness

    def create(self):
        s = mdb.models['Model-1'].ConstrainedSketch(name=self.name + '_sketch', sheetSize=max(self.axis1, self.axis2, self.thickness))
        g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints    
        s.sketchOptions.setValues(decimalPlaces=4)
        s.setPrimaryObject(option=STANDALONE)
        axis1, axis2 = self.axis1, self.axis2
        s.Line(point1=(axis1,0), point2=(0,axis2))
        s.Line(point1=(0,axis2), point2=(-axis1,0))
        s.Line(point1=(-axis1,0), point2=(0,-axis2))
        s.Line(point1=(0,-axis2), point2=(axis1,0))
        self.part = mdb.models['Model-1'].Part(name=self.name, dimensionality=THREE_D, type=DEFORMABLE_BODY)
        self.part.BaseSolidExtrude(sketch=s, depth=self.thickness)

    def set_section(self, section):
        c = self.part.cells
        region = self.part.Set(cells=c, name='Tool material region')
        self.part.SectionAssignment(region=region, sectionName=section.name, offset=0.0, 
            offsetType=MIDDLE_SURFACE, offsetField='', thicknessAssignment=FROM_SECTION)        

    def mesh(self):
        self.part.seedPart(size=0.0025, deviationFactor=0.1, minSizeFactor=0.1)
        c = self.part.cells
        pickedRegions = c.getSequenceFromMask(mask=('[#1 ]', ), )
        # self.part.setMeshControls(regions=c, technique=SWEEP, algorithm=ADVANCING_FRONT)
        self.part.setMeshControls(regions=pickedRegions, elemShape=TET, technique=FREE)
        elemType1 = mesh.ElemType(elemCode=C3D20R)
        elemType2 = mesh.ElemType(elemCode=C3D15)
        elemType3 = mesh.ElemType(elemCode=C3D10M)
        cells = c.getSequenceFromMask(mask=('[#1 ]', ), )
        pickedRegions =(cells, )
        self.part.setElementType(regions=pickedRegions, elemTypes=(elemType1, elemType1, elemType1))
        self.part.generateMesh()


class Assembly:

    def __init__(self, workpiece, tool):
        self.workpiece = workpiece
        self.tool = tool
        self.a = mdb.models['Model-1'].rootAssembly
        self.a.DatumCsysByDefault(CARTESIAN)
        self.a.Instance(name=workpiece.name, part=workpiece.part, dependent=ON)
        self.a.Instance(name=tool.name, part=tool.part, dependent=ON)
        self.a.translate(instanceList=(tool.name, ), vector=(0,(workpiece.outer_d+workpiece.inner_d)/4 + tool.axis2,0))
        self.a.translate(instanceList=(tool.name, ), vector=(0,0,-tool.thickness))
        self.a.rotate(instanceList=(tool.name, ), axisPoint=(0, 0, 0), axisDirection=(0.0, 1.0, 0), angle=-90.0)
        self.a.translate(instanceList=(tool.name, ), vector=(0,0,-tool.axis1/4))
        # self.a.rotate(instanceList=(tool.name, ), axisPoint=(0, 0, 0), axisDirection=(1, 0.0, 0.0), angle=-90.0)
        # self.a.translate(instanceList=(tool.name, ), vector=(0.5*(workpiece.w_width + workpiece.b_width), workpiece.w_height + workpiece.b_height, workpiece.length))
        # self.a.translate(instanceList=(tool.name, ), vector=(0, -0.75 * tool.cutter_length, 0))
        # self.a.translate(instanceList=(tool.name, ), vector=(-0.9 * tool.diameter, 0, 0))
        # self.a.translate(instanceList=(tool.name, ), vector=(0, 0, 0.5 * tool.diameter))
        session.viewports['Viewport: 1'].setValues(displayedObject=self.a)
        session.viewports['Viewport: 1'].view.fitView()
    
    def tool_bc(self):
        f1 = self.a.instances[self.tool.name].faces
        # faces1 = f1.findAt((0,0,0))
        region = self.a.Set(faces=f1, name='Whole tool')
        mdb.models['Model-1'].EncastreBC(name='BC whole tool', createStepName='Initial', region=region, localCsys=None)

    def workpiece_bc(self, step):
        name = "BC workpiece"
        workpiece_inst = self.a.instances[self.workpiece.name]
        model = mdb.models['Model-1']
        c1 = workpiece_inst.cells
        f1 = workpiece_inst.faces
        e1 = workpiece_inst.edges
        v1 = workpiece_inst.vertices
        region = self.a.Set(vertices=v1, edges=e1, faces=f1, cells=c1, name='Tool set')
        
        model.DisplacementBC(name=name, createStepName='Initial', 
            region=region, u1=SET, u2=SET, u3=SET, ur1=SET, ur2=SET, ur3=SET, 
            amplitude=UNSET, distributionType=UNIFORM, fieldName='', localCsys=None)
        model.TabularAmplitude(name='Amplitude-1', timeSpan=STEP, 
            smooth=SOLVER_DEFAULT, data=((0.0, 0.0), (MAX_TIME, 1.0)))
        model.boundaryConditions[name].setValuesInStep(stepName=step.name, u3=-self.workpiece.length, amplitude='Amplitude-1')


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
    # create material and sections
    steel = Material("Steel", density= 7870, young=2e11, poisson=0.29, 
        A=375e6, B=552e6, n=0.457, 
        d1=0.25, d2=4.38, d3=2.68, ref_strain_rate=1.0, disp_at_failure=0.1)
    steel_section = mdb.models['Model-1'].HomogeneousSolidSection(name='Steel section', material=steel.name, thickness=None)
    alu = Material("Alu", density= 2700, young=70e9, poisson=0.33, 
        A=3.241e8, B=1.138e8, n=0.42, 
        d1=-0.77, d2=1.45, d3=-0.47, ref_strain_rate=1.0, disp_at_failure=1e-4)
    alu_section = mdb.models['Model-1'].HomogeneousSolidSection(name='Alu section', material=alu.name, thickness=None)
    # create workpiece
    workpiece = Workpiece("Tube", inner_d=mm(40), outer_d=mm(52),length=mm(80))
    workpiece.create()
    workpiece.set_section(alu_section)
    workpiece.mesh()
    tool = Tool("Tool", axis1=mm(10), axis2=mm(20), thickness=mm(3))
    tool.create()
    tool.set_section(steel_section)
    tool.mesh()
    assembly = Assembly(workpiece,tool)
    step1 = Step("Step-1")
    Interaction(step1, friction=0.15)
    assembly.workpiece_bc(step1)
    assembly.tool_bc()
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



