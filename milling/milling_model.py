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

executeOnCaeStartup()
Mdb()
MAX_TIME = 0.001
model = mdb.models['Model-1']

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
        model.Material(name=name)
        self.material = model.materials[name]
        self.material.Density(table=((density, ), ))
        self.material.Elastic(table=((young, poisson), ))
        self.material.Plastic(hardening=JOHNSON_COOK, table=((A, B, n, 0.0, 0.0, 0.0), ))
        self.material.JohnsonCookDamageInitiation(table=((d1, d2, d3, 0.0, 0.0, 0.0, 0.0, ref_strain_rate), ))
        self.material.johnsonCookDamageInitiation.DamageEvolution(type=DISPLACEMENT, table=((disp_at_failure, ), ))


class Workpiece:
    def __init__(self, name, length, w_height, w_width, b_height, b_width):
        self.name = name
        self.length = length
        self.w_height = w_height
        self.w_width = w_width
        self.b_height = b_height
        self.b_width = b_width
        self.part = None
        self.c_point = (0.5 * self.b_width, self.b_height + 0.5 * self.w_height, 0.5 * self.length)
        self.c1_point = (0.25 * (self.b_width - self.w_width), 0.5 * self.b_height, 0.5 * self.length)
        self.c2_point = (0.5 * self.b_width, 0.5 * self.b_height, 0.5 * self.length)
        self.c3_point = (self.b_width - 0.25 * (self.b_width - self.w_width), 0.5 * self.b_height, 0.5 * self.length)

    def create(self):
        s = model.ConstrainedSketch(name=self.name + '_sketch', sheetSize=max(self.length, self.w_height, self.w_width, 
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


        p = model.Part(name=self.name, dimensionality=THREE_D, type=DEFORMABLE_BODY)
        self.part = model.parts[self.name]
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


    def set_section(self, section):
        c = self.part.cells
        pickedRegions = c.findAt((self.c1_point, ), (self.c2_point, ), (self.c3_point, ), (self.c_point, ))
        region = self.part.Set(cells=pickedRegions, name='Material')
        self.part.SectionAssignment(region=region, sectionName=section.name, offset=0.0, 
            offsetType=MIDDLE_SURFACE, offsetField='', thicknessAssignment=FROM_SECTION)
        

    def mesh(self):
        e = self.part.edges
        e1_point = (0.5* self.b_width, self.b_height, 0)
        e2_point = (0.5* self.b_width, self.b_height+ self.w_height, 0)
        e3_point = (0.5* self.b_width, self.b_height, self.length)
        e4_point = (0.5* self.b_width, self.b_height+ self.w_height, self.length)
        pickedEdges = e.findAt((e1_point, ), (e2_point, ), (e3_point, ), (e4_point, ))
        self.part.seedEdgeByNumber(edges=pickedEdges, number=int(2*self.w_width * 1e3), constraint=FINER)
        
        e1_point = (0.5* (self.b_width - self.w_width), self.b_height + 0.5 * self.w_height, 0)
        e2_point = (0.5* (self.b_width + self.w_width), self.b_height + 0.5 * self.w_height, 0)
        e3_point = (0.5* (self.b_width - self.w_width), self.b_height + 0.5 * self.w_height, self.length)
        e4_point = (0.5* (self.b_width + self.w_width), self.b_height + 0.5 * self.w_height, self.length)
        pickedEdges = e.findAt((e1_point, ), (e2_point, ), (e3_point, ), (e4_point, ))
        self.part.seedEdgeByNumber(edges=pickedEdges, number=int(1*self.w_height * 1e3), constraint=FINER)

        e1_point = (0.5* (self.b_width - self.w_width), self.b_height, 0.5 * self.length)
        e2_point = (0.5* (self.b_width + self.w_width), self.b_height, 0.5 * self.length)
        e3_point = (0.5* (self.b_width - self.w_width), self.b_height + self.w_height, 0.5 * self.length)
        e4_point = (0.5* (self.b_width + self.w_width), self.b_height + self.w_height, 0.5 * self.length)
        pickedEdges = e.findAt((e1_point, ), (e2_point, ), (e3_point, ), (e4_point, ))
        self.part.seedEdgeByNumber(edges=pickedEdges, number=int(2*self.length * 1e3), constraint=FINER)
        self.part.seedPart(size=0.005, deviationFactor=0.1, minSizeFactor=0.1)
        
        c = self.part.cells
        pickedRegions = c.findAt((self.c1_point, ), (self.c2_point, ), (self.c3_point, ))
        self.part.setMeshControls(regions=pickedRegions, technique=SWEEP, algorithm=ADVANCING_FRONT)

        elemType1 = mesh.ElemType(elemCode=C3D8R, elemLibrary=EXPLICIT, 
            kinematicSplit=AVERAGE_STRAIN, secondOrderAccuracy=OFF, 
            hourglassControl=DEFAULT, distortionControl=DEFAULT, elemDeletion=ON)
        elemType2 = mesh.ElemType(elemCode=C3D6, elemLibrary=EXPLICIT)
        elemType3 = mesh.ElemType(elemCode=C3D4, elemLibrary=EXPLICIT)
        cells = c.findAt((self.c_point, ))
        pickedRegions =(cells, )
        self.part.setElementType(regions=pickedRegions, elemTypes=(elemType1, elemType2, elemType3))
        self.part.generateMesh()

        
class Tool:
    def __init__(self, name, scale):
        self.name = name
        self.step_file = 'D:/ereme/GoogleDrive/PostGraduate/AbaqusModels/mill3Dmodel/4-flute-flat-endmill-1.snapshot.4/endmill.stp'
        self.diameter = mm(10)  # from file
        self.length = mm(100)
        self.cutter_length = mm(22)
        self.scale = scale
        self.part = None

    def create(self):
        step = mdb.openStep(self.step_file)
        model.PartFromGeometryFile(name=self.name, geometryFile=step, 
                combine=False, dimensionality=THREE_D, type=DEFORMABLE_BODY, scale=self.scale)
        self.part = model.parts[self.name]

    def set_section(self, section):
        c = self.part.cells
        cells = c.getSequenceFromMask(mask=('[#1 ]', ), )
        region = self.part.Set(cells=cells, name='Material')
        self.part.SectionAssignment(region=region, sectionName=section.name, offset=0.0, 
            offsetType=MIDDLE_SURFACE, offsetField='', thicknessAssignment=FROM_SECTION)

    def mesh(self):
        self.part.seedPart(size=0.0025, deviationFactor=0.1, minSizeFactor=0.1)
        c = self.part.cells
        pickedRegions = c.getSequenceFromMask(mask=('[#1 ]', ), )
        self.part.setMeshControls(regions=pickedRegions, elemShape=TET, technique=FREE)
        elemType1 = mesh.ElemType(elemCode=C3D20R)
        elemType2 = mesh.ElemType(elemCode=C3D15)
        elemType3 = mesh.ElemType(elemCode=C3D10M)
        cells = c.getSequenceFromMask(mask=('[#1 ]', ), )
        pickedRegions =(cells, )
        self.part.setElementType(regions=pickedRegions, elemTypes=(elemType1, elemType2, elemType3))
        self.part.generateMesh()


    def __str__(self):
        return "Tool: name={}, scale={}".format(self.name, self.scale)

    def __repr__(self):
        return str(self)


class Assembly:

    def __init__(self, workpiece, tool):
        x_delta_mill = mm(0)
        y_delta_mill = mm(-5)
        z_delta_mill = mm(0)
        self.workpiece = workpiece
        self.tool = tool
        self.a = model.rootAssembly
        self.a.DatumCsysByDefault(CARTESIAN)
        self.a.Instance(name=workpiece.name, part=workpiece.part, dependent=ON)
        self.a.Instance(name=tool.name, part=tool.part, dependent=ON)
        self.a.rotate(instanceList=(tool.name, ), axisPoint=(0, 0, 0), axisDirection=(1, 0.0, 0.0), angle=-90.0)
        self.a.translate(instanceList=(tool.name, ), vector=(0.5*(workpiece.w_width + workpiece.b_width), workpiece.w_height + workpiece.b_height, workpiece.length))
        self.a.translate(instanceList=(tool.name, ), vector=(0, -0.75 * tool.cutter_length, 0))
        self.a.translate(instanceList=(tool.name, ), vector=(-0.9 * tool.diameter, 0, 0))
        self.a.translate(instanceList=(tool.name, ), vector=(0, 0, 0.5 * tool.diameter))
        self.a.translate(instanceList=(tool.name, ), vector=(x_delta_mill, y_delta_mill, z_delta_mill))
        self.workpieceInstance = self.a.instances[workpiece.name]
        self.toolInstance = self.a.instances[tool.name]
        
        self.millCenterPoint=(0.5*(workpiece.w_width + workpiece.b_width) -0.9 * tool.diameter + x_delta_mill,
                tool.length + workpiece.w_height + workpiece.b_height -0.75 * tool.cutter_length + y_delta_mill,
                workpiece.length + 0.5 * tool.diameter + z_delta_mill)

        
        session.viewports['Viewport: 1'].setValues(displayedObject=self.a)
        session.viewports['Viewport: 1'].view.fitView()    


    def workpiece_bc(self):
        name = 'Workpiece Bottom Encastre'
        p1 = (0.25 * (workpiece.b_width - workpiece.w_width), 0, 0.5 * workpiece.length)
        p2 = (0.5 * workpiece.b_width, 0, 0.5 * workpiece.length)
        p3 = (workpiece.b_width - 0.25 * (workpiece.b_width - workpiece.w_width), 0, 0.5 * workpiece.length)
        f1 = self.workpieceInstance.faces
        faces1 = f1.findAt((p1, ), (p3, ), (p2, ) )
        region = self.a.Set(faces=faces1, name=name)
        model.EncastreBC(name=name, createStepName='Initial', region=region, localCsys=None)

    def tool_bc(self, step):
        name = 'Tool RP Velocity'
        # create reference point
        refPoint = self.a.ReferencePoint(point=self.millCenterPoint)
        refPointRegion=regionToolset.Region(referencePoints=(self.a.referencePoints[refPoint.id],))
        toolBodyRegion = self.a.Set(cells=self.toolInstance.cells, name=name)
        model.RigidBody(name='Tool RP RigidBody', refPointRegion=refPointRegion, bodyRegion=toolBodyRegion)
        model.rootAssembly.features.changeKey(fromName='RP-1', toName='Tool Top RP')
        # create velocity BC
        model.VelocityBC(name=name, createStepName='Initial', 
            region=refPointRegion, v1=0.0, v2=0.0, v3=0.0, vr1=0.0, vr2=0.0, vr3=0.0, 
            amplitude=UNSET, localCsys=None, distributionType=UNIFORM, fieldName='')
        model.boundaryConditions[name].setValuesInStep(stepName=step.name, v3=-120)
        model.boundaryConditions[name].setValuesInStep(stepName=step.name, vr2=6000)


class Step:

    def __init__(self, name, previous='Initial'):
        self.name = name
        model.ExplicitDynamicsStep(name=name, previous=previous, description='', timePeriod=MAX_TIME)

class Interaction:

    def __init__(self, step, friction):
        # create interaction propery
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
    workpiece = Workpiece("Workpiece", length=mm(100), w_height=mm(60), w_width=mm(5), b_height=mm(10), b_width=mm(50))
    workpiece.create()
    workpiece.mesh()
    
    tool = Tool("Tool",  scale=1e-3)
    tool.create()
    tool.mesh()
    # Create materials
    steel = Material("Steel", density= 7870, young=2e11, poisson=0.29, 
        A=375e6, B=552e6, n=0.457, 
        d1=0.25, d2=4.38, d3=2.68, ref_strain_rate=1.0, disp_at_failure=0.1)
    steel_section = model.HomogeneousSolidSection(name='Steel section', material=steel.name, thickness=None)
    alu = Material("Alu", density= 2700, young=70e9, poisson=0.33, 
        A=3.241e8, B=1.138e8, n=0.42, 
        d1=-0.77, d2=1.45, d3=-0.47, ref_strain_rate=1.0, disp_at_failure=1e-4)
    alu_section = model.HomogeneousSolidSection(name='Alu section', material=alu.name, thickness=None)
    # assign materials 
    workpiece.set_section(alu_section)
    tool.set_section(steel_section)

    assembly = Assembly(workpiece, tool)
    assembly.workpiece_bc()
    
    step1 = Step("Step-1")
    Interaction(step1, friction=0.15)    
    assembly.tool_bc(step1)
    model.fieldOutputRequests['F-Output-1'].setValues(variables=('A', 'CSTRESS', 'EVF', 'LE', 'PE', 'PEEQ', 'PEEQVAVG', 'PEVAVG', 'RF', 'S', 
        'STATUS', # important!
        'SVAVG', 'U', 'V'))

    mdb.Job(name='Chipformation-v2', model='Model-1', description='', type=ANALYSIS, 
        atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, 
        memoryUnits=PERCENTAGE, explicitPrecision=SINGLE, 
        nodalOutputPrecision=SINGLE, echoPrint=OFF, modelPrint=OFF, 
        contactPrint=OFF, historyPrint=OFF, userSubroutine='', scratch='', 
        resultsFormat=ODB, parallelizationMethodExplicit=DOMAIN, numDomains=4, 
        activateLoadBalancing=False, multiprocessingMode=DEFAULT, numCpus=4)

    mdb.jobs['Chipformation-v2'].submit(consistencyChecking=OFF)