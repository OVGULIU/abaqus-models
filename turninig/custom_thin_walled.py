from abaqus import *
from abaqusConstants import *
from caeModules import *
from driverUtils import executeOnCaeStartup
import math
import tempfile
import sys
import os
import argparse
import math
import numpy as np
executeOnCaeStartup()

MAX_TIME = 0.001



def cart2pol(x, y):
    rho = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y, x)
    return(rho, math.degrees(phi))

def pol2cart(rho, phi):
    phi = math.radians(phi)
    x = rho * np.cos(phi)
    y = rho * np.sin(phi)
    return(x, y)   

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

class Material_Explicit:

    def __init__(self, name, density, young, poisson, A, B, n, d1, d2, d3, ref_strain_rate, disp_at_failure):
        self.name = name
        mdb.models['Model-1'].Material(name=name)
        self.material = mdb.models['Model-1'].materials[name]
        self.material.Density(table=((density, ), ))
        self.material.Elastic(table=((young, poisson), ))
        self.material.Plastic(hardening=JOHNSON_COOK, table=((A, B, n, 0.0, 0.0, 0.0), ))
        self.material.JohnsonCookDamageInitiation(table=((d1, d2, d3, 0.0, 0.0, 0.0, 0.0, ref_strain_rate), ))
        self.material.johnsonCookDamageInitiation.DamageEvolution(type=DISPLACEMENT, table=((disp_at_failure, ), ))


class Material:

    def __init__(self, name, young, poisson):
        self.name = name
        mdb.models['Model-1'].Material(name=name)
        mdb.models['Model-1'].materials[name].Elastic(table=((young, poisson), ))
        self.manterial =  mdb.models['Model-1'].materials[name]


class Workpiece:

    def __init__(self, path, outer_radius):
        self.name = 'Workpiece'
        self.outer_radius = outer_radius
        step = mdb.openStep(path, scaleFromFile=OFF)
        mdb.models['Model-1'].PartFromGeometryFile(name=self.name, geometryFile=step, 
            combine=False, dimensionality=THREE_D, type=DEFORMABLE_BODY, scale=0.001)
        self.part = mdb.models['Model-1'].parts[self.name]


    def partition(self, number):
        for angle in np.arange(0, 360, 360/number):
            self.part.PartitionCellByPlaneThreePoints(cells=self.part.cells, point1=(0,0,0), point2=(1,0,0), point3=(0,) +pol2cart(1, 30+angle))
        

    def mesh(self):
        self.part.seedPart(size=0.005, deviationFactor=0.01, minSizeFactor=0.01)
        self.part.generateMesh()

    def set_section(self, section):
        region = self.part.Set(cells=self.part.cells, name='Material region')
        self.part.SectionAssignment(region=region, sectionName=section.name, offset=0.0, 
            offsetType=MIDDLE_SURFACE, offsetField='', thicknessAssignment=FROM_SECTION)
        


class Jaw:

    
    def __init__(self, length, width, height):
        self.name = "Jaw"
        self.length = length
        self.width = width
        self.height = height
        s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=0.1) # 
        s.sketchOptions.setValues(decimalPlaces=3)
        s.setPrimaryObject(option=STANDALONE)
        s.rectangle(point1=(0.0, 0.0), point2=(length, height)) #
        self.part = mdb.models['Model-1'].Part(name=self.name, dimensionality=THREE_D,  type=DEFORMABLE_BODY) # 
        self.part.BaseSolidExtrude(sketch=s, depth=width) #
        s.unsetPrimaryObject()

    def partition(self):
        p = mdb.models['Model-1'].parts[self.name]
        length, width, height = self.length, self.width , self.height 
        p1 = p.DatumPointByCoordinate(coords=(length/2, 0.0, 0.0)).id
        p2 = p.DatumPointByCoordinate(coords=(length/2, height, width)).id
        p3 = p.DatumPointByCoordinate(coords=(length/2, 0, width)).id
        c, v1, e1, d1 = p.cells, p.vertices, p.edges, p.datums
        p.PartitionCellByPlaneThreePoints(cells=c, point1=d1[p1], point2=d1[p2], point3=d1[p3])
            
        p1 = p.DatumPointByCoordinate(coords=(0.0, 0.0, width/2)).id
        p2 = p.DatumPointByCoordinate(coords=(length, 0.0, width/2)).id
        p3 = p.DatumPointByCoordinate(coords=(length, height, width/2)).id
        c, v1, e1, d1 = p.cells, p.vertices, p.edges, p.datums
        p.PartitionCellByPlaneThreePoints(cells=c, point1=d1[p1], point2=d1[p2], point3=d1[p3])
        
    def mesh(self):
        self.part.seedPart(size=0.01, deviationFactor=0.1, minSizeFactor=0.1)
        self.part.generateMesh()

    def set_section(self, section):
        region = self.part.Set(cells=self.part.cells, name='Material region')
        self.part.SectionAssignment(region=region, sectionName=section.name, offset=0.0, 
            offsetType=MIDDLE_SURFACE, offsetField='', thicknessAssignment=FROM_SECTION)

class InteractionProperty:

    def __init__(self):
        self.name = "Interaction-Property"
        mdb.models['Model-1'].ContactProperty(self.name)
        mdb.models['Model-1'].interactionProperties[self.name].TangentialBehavior(formulation=ROUGH)
        mdb.models['Model-1'].interactionProperties[self.name].NormalBehavior(
            pressureOverclosure=HARD, allowSeparation=OFF, contactStiffness=DEFAULT, 
            contactStiffnessScaleFactor=1.0, clearanceAtZeroContactPressure=0.0, 
            stiffnessBehavior=LINEAR, constraintEnforcementMethod=PENALTY)

class Step:

    def __init__(self, name, previous='Initial'):
        self.name = name
        mdb.models['Model-1'].StaticStep(name=name, previous=previous)


class Assembly:

    class Jaw:

        def __init__(self, index, name):
            self.index = index
            self.name = name
            self.angle = 360/3 * (index-1)

        def __hash__(self):
            return self.index

    def __init__(self, workpiece, jaw, jaw_force):
        self.workpiece = workpiece
        self.jaw = jaw
        self.jaw_force = jaw_force
        self.a = mdb.models['Model-1'].rootAssembly
        self.a.DatumCsysByDefault(CARTESIAN)
        self.a.Instance(name=workpiece.name, part=workpiece.part, dependent=ON)
        self.a.Instance(name=jaw.name, part=jaw.part, dependent=ON)
        self.a.translate(instanceList=('Workpiece', ), vector=(-0.1, 0.0, 0.0))
        self.a.rotate(instanceList=(jaw.name, ), axisPoint=(0, 0, 0), axisDirection=(1, 0.0, 0.0), angle=-90.0)
        self.a.translate(instanceList=(jaw.name, ), vector=(0,-jaw.width/2, 0))
        self.a.translate(instanceList=('Jaw', ), vector=(-jaw.length, 0.0, workpiece.outer_radius + jaw.height))
        self.a.RadialInstancePattern(instanceList=('Jaw', ), point=(0.0, 0.0, 0.0), 
            axis=(1.0, 0.0, 0.0), number=3, totalAngle=360.0)
        self.a.features.changeKey(fromName='Jaw', toName='Jaw-rad-1')
        self.interactionProperty = InteractionProperty()
        self.jaws = [Assembly.Jaw(i, 'Jaw-rad-'+str(i)) for i in range(1,4)]
        self.CSYSs = dict()
        for a_jaw in self.jaws:
            self._create_interaction(a_jaw, self.interactionProperty)
            self.CSYSs[a_jaw] = self._create_CSYS(a_jaw)
            self._create_jaw_BSs(a_jaw)
        self.step = Step("Step-1")
        for a_jaw in self.jaws:
            self._apply_jaw_force(a_jaw, jaw_force)


    def _create_CSYS(self, jaw):
        rho, phi = cart2pol(self.workpiece.outer_radius, 0)
        print(jaw.angle)
        z, y = pol2cart(rho, phi - jaw.angle)
        origin = (0, y, z)

        rho, phi = cart2pol(self.workpiece.outer_radius, self.jaw.width)
        z, y = pol2cart(rho, phi - jaw.angle)
        point2 = (0, y, z)

        return self.a.DatumCsysByThreePoints(origin=origin, point1=(0,0,0), point2=point2, name = jaw.name + '_CSYS', coordSysType=CARTESIAN)
        

    def _create_interaction(self, jaw, property):

        # create jaw slave region
        rho1, phi1 = cart2pol(self.workpiece.outer_radius, 0.25 * self.jaw.width)
        z1, y1 = pol2cart(rho1, phi1 - jaw.angle)

        rho2, phi2 = cart2pol(self.workpiece.outer_radius, -0.25 * self.jaw.width)
        z2, y2 = pol2cart(rho2, phi2 - jaw.angle)

        jaw_instance = self.a.instances[jaw.name]
        workpiece = self.a.instances['Workpiece']

        jaw_faces = jaw_instance.faces.findAt(
            ((-0.25 * self.jaw.length, y1, z1), ), 
            ((-0.25 * self.jaw.length, y2, z2), ),
            ((-0.75 * self.jaw.length, y1, z1), ), 
            ((-0.75 * self.jaw.length, y2, z2), )
            )
        region_jaw=self.a.Surface(side1Faces=jaw_faces, name=jaw.name + '_master_surf')       


        # create workpiece master region   
        rho1_w, phi1_w = cart2pol(z1, y1)
        z1_w, y1_w = pol2cart(self.workpiece.outer_radius, phi1_w)

        rho2_w, phi2_w = cart2pol(z2, y2)
        z2_w, y2_w = pol2cart(self.workpiece.outer_radius, phi2_w)
        
        workpiece_faces = workpiece.faces.findAt(
            ((-0.25 * self.jaw.length, y1_w, z1_w), ), 
            ((-0.25 * self.jaw.length, y2_w, z2_w), ),
            ((-0.75 * self.jaw.length, y1_w, z1_w), ), 
            ((-0.75 * self.jaw.length, y2_w, z2_w), )
            )
        region_workpiece=self.a.Surface(side1Faces=workpiece_faces, name=jaw.name + '_slave_surf')
        # return
        mdb.models['Model-1'].SurfaceToSurfaceContactStd(name='Interaction-'+jaw.name, 
            createStepName='Initial', master=region_jaw, slave=region_workpiece, sliding=FINITE, 
            thickness=ON, interactionProperty=property.name, adjustMethod=NONE, 
            initialClearance=OMIT, datumAxis=None, clearanceRegion=None)


    def _create_jaw_BSs(self, jaw):
        jaw_instance = self.a.instances[jaw.name]
        jaw_faces = jaw_instance.faces.getSequenceFromMask(mask=('[#d5428 ]', ), ) # hack, replace it with findAt
        region = self.a.Set(faces = jaw_faces, name='Jaw_BS_set-'+jaw.name)
        datum = self.a.datums[self.CSYSs[jaw].id]
        mdb.models['Model-1'].DisplacementBC(name='BC-'+jaw.name, createStepName='Initial', 
            region=region, u1=UNSET, u2=SET, u3=SET, ur1=UNSET, ur2=UNSET, ur3=UNSET, 
            amplitude=UNSET, distributionType=UNIFORM, fieldName='', localCsys=datum)

    def _apply_jaw_force(self, jaw, value):
        vertices = self.a.instances[jaw.name].vertices
        verts1 = vertices.getSequenceFromMask(mask=('[#1 ]', ), ) # hack, replace it with findAt
        region = self.a.Set(vertices=verts1, name='Jaw-force-region-'+jaw.name)        
        datum = self.a.datums[self.CSYSs[jaw].id]
        mdb.models['Model-1'].ConcentratedForce(name='Load-'+jaw.name, createStepName=self.step.name, 
            region=region, cf1=value, distributionType=UNIFORM, field='', localCsys=datum)
   


if __name__ =="__main__":
    Mdb()
    # create materials
    steel = Material_Explicit("Steel", density= 7870, young=2e11, poisson=0.29, 
        A=375e6, B=552e6, n=0.457, 
        d1=0.25, d2=4.38, d3=2.68, ref_strain_rate=1.0, disp_at_failure=0.1)
    steel = Material('Steel', 210e9, 0.29)
    steel_section = mdb.models['Model-1'].HomogeneousSolidSection(name='Steel section', material=steel.name, thickness=None)
    
    alu = Material_Explicit("Alu", density= 2700, young=70e9, poisson=0.33, 
        A=3.241e8, B=1.138e8, n=0.42, 
        d1=-0.77, d2=1.45, d3=-0.47, ref_strain_rate=1.0, disp_at_failure=1e-4)
    alu = Material('Alu', 70e9, 0.28)
    # alu = Material('Alu', 100e15, 0.28)
    alu_section = mdb.models['Model-1'].HomogeneousSolidSection(name='Alu section', material=alu.name, thickness=None)

    # create workpiece
    workpiece = Workpiece('D:/ereme/GoogleDrive/PostGraduate/ZhAD/parts/meshed/part1.STEP', outer_radius=0.046)
    workpiece.partition(3)
    workpiece.mesh()
    workpiece.set_section(alu_section)
    session.viewports['Viewport: 1'].setValues(displayedObject=workpiece.part)

    #  create jaw
    jaw = Jaw(length=0.02, width=0.02, height=0.02)
    jaw.partition()
    jaw.mesh()
    jaw.set_section(steel_section)
    session.viewports['Viewport: 1'].setValues(displayedObject=jaw.part)

    # assembly
    assembly = Assembly(workpiece, jaw, jaw_force=400)
    session.viewports['Viewport: 1'].setValues(displayedObject=assembly.a)
    session.viewports['Viewport: 1'].view.setValues(session.views['Iso'])