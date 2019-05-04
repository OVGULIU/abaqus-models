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
        model.interactions['Interaction'].includedPairs.setValuesInStep(stepName=step.name, useAllstar=ON)
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
        self.part.seedEdgeByNumber(edges=pickedEdges, number=int(1*self.w_height * 1e3), constraint=FINER)

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
        self.part.seedEdgeByNumber(edges=pickedEdges, number=int(2*self.t * 1e3), constraint=FINER)
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
        # self.a.rotate(instanceList=(tool.name, ), axisPoint=(0, 0, 0), axisDirection=(1, 0.0, 0.0), angle=-90.0)
        # a.rotate(instanceList=('Tool', ), axisPoint=(0.0375, 0.07, 0.0), axisDirection=(0.0, 0.0, 0.05), angle=180.0)
        # self.a.translate(instanceList=(tool.name, ), vector=(0.5*(workpiece.w_width + workpiece.b_width), workpiece.w_height + workpiece.b_height, workpiece.length))
        self.a.rotate(instanceList=('Tool', ), axisPoint=(0.0, 0.0, 0.0), axisDirection=(0.00625, 0.0, 0.0), angle=90.0)
        self.a.rotate(instanceList=('Tool', ), axisPoint=(0.0, 0.0, 0.0), axisDirection=(0.0, 1.0, 0.0), angle=45.0)
        self.workpieceInstance = self.a.instances[workpiece.name]
        self.toolInstance = self.a.instances[tool.name]
        
        self.a.translate(instanceList=('Tool', ), vector=(0.0645, 0.07, -0.0037))
        self.a.translate(instanceList=('Tool', ), vector=(0, 0.005, 0))
        self.a.translate(instanceList=('Tool', ), vector=(0.0, 0.0, -0.0125))
        session.viewports['Viewport: 1'].setValues(displayedObject=self.a)
        session.viewports['Viewport: 1'].view.fitView()    

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

    workpiece = Workpiece("Workpiece", length=mm(50), w_height=mm(60), w_width=mm(25), b_height=mm(10), b_width=mm(50))
    workpiece.create()
    ti6alv = Ti6AlV()
    mdb.models['Model-1'].materials['Ti6AlV'].johnsonCookDamageInitiation.DamageEvolution(type=DISPLACEMENT, table=((0.002, ), ))

    alu = Material("Alu", density= 2700, young=70e9, poisson=0.33, 
        A=3.241e8, B=1.138e8, n=0.42, 
        d1=-0.77, d2=1.45, d3=-0.47, ref_strain_rate=1.0, disp_at_failure=1e-4)

    mdb.models['Model-1'].materials["Alu"].Conductivity(table=((7.2, ), ))

    workpiece.set_material(alu)
    workpiece.mesh()

    step1 = Step("Step-1")
    Interaction(step1)    
    
    assembly = Assembly(workpiece, tool)
    assembly.workpiece_bc()
    assembly.tool_bc(step1)
    model.fieldOutputRequests['F-Output-1'].setValues(variables=('A', 'CSTRESS', 'EVF', 'LE', 'PE', 'PEEQ', 'PEEQVAVG', 'PEVAVG', 'RF', 'S', 
        'STATUS', # important!
        'SVAVG', 'U', 'V'))
    model.fieldOutputRequests['F-Output-1'].setValues(
    numIntervals=400)
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