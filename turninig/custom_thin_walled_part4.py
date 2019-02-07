"""This model doesn't work yet"""

"""If an error occures on import nearestNodeModule go manually
in Abaqus main menu to Plugins->Tools->Find Nearest Node. 
It will initialize the module"""
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
Mdb()
model = mdb.models['Model-1']
OX = (1, 0, 0)
OY = (0, 1, 0)
OZ = (0, 0, 1)
MAX_TIME = 0.001

import inspect
filename = inspect.getframeinfo(inspect.currentframe()).filename
cd = os.path.dirname(os.path.abspath(filename))

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

def rotate(point, axis, theta):
    """
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians.
    """
    def rotation_matrix(axis, theta):
        axis = np.asarray(axis)
        axis = axis / math.sqrt(np.dot(axis, axis))
        a = math.cos(theta / 2.0)
        b, c, d = -axis * math.sin(theta / 2.0)
        aa, bb, cc, dd = a * a, b * b, c * c, d * d
        bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
        return np.array([[aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)], [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)], [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc]])
    return tuple(np.dot(rotation_matrix(axis, theta), point))


class Material_Explicit:

    def __init__(self, name, density, young, poisson, A, B, n, d1, d2, d3, ref_strain_rate, disp_at_failure):
        self.name = name
        model.Material(name=name)
        self.material = model.materials[name]
        self.material.Density(table=((density, ), ))
        self.material.Elastic(table=((young, poisson), ))
        self.material.Plastic(hardening=JOHNSON_COOK, table=((A, B, n, 0.0, 0.0, 0.0), ))
        self.material.JohnsonCookDamageInitiation(table=((d1, d2, d3, 0.0, 0.0, 0.0, 0.0, ref_strain_rate), ))
        self.material.johnsonCookDamageInitiation.DamageEvolution(type=DISPLACEMENT, table=((disp_at_failure, ), ))


class Material:

    def __init__(self, name, young, poisson):
        self.name = name
        model.Material(name=name)
        model.materials[name].Elastic(table=((young, poisson), ))
        self.manterial =  model.materials[name]


class Workpiece:

    def __init__(self, path, outer_radius, length):
        self.name = 'Workpiece'
        self.outer_radius = outer_radius
        self.length = length
        step = mdb.openStep(path, scaleFromFile=OFF)
        model.PartFromGeometryFile(name=self.name, geometryFile=step, 
            combine=False, dimensionality=THREE_D, type=DEFORMABLE_BODY, scale=0.001)
        self.part = model.parts[self.name]


    def partition(self, number):
        p = mdb.models['Model-1'].parts['Workpiece']
        c = p.cells
        pickedCells = c.getSequenceFromMask(mask=('[#1 ]', ), )
        e, v, d = p.edges, p.vertices, p.datums
        p.PartitionCellByPlanePointNormal(point=v[62], normal=e[94], cells=pickedCells)
        p = mdb.models['Model-1'].parts['Workpiece']
        c = p.cells
        pickedCells = c.getSequenceFromMask(mask=('[#2 ]', ), )
        e1, v1, d1 = p.edges, p.vertices, p.datums
        p.PartitionCellByPlanePointNormal(point=v1[93], normal=e1[101], cells=pickedCells)
        p = mdb.models['Model-1'].parts['Workpiece']
        c = p.cells
        pickedCells = c.getSequenceFromMask(mask=('[#4 ]', ), )
        e1, v1, d1 = p.edges, p.vertices, p.datums
        p.PartitionCellByPlanePointNormal(point=v1[86], normal=e1[95], cells=pickedCells)
        p = mdb.models['Model-1'].parts['Workpiece']
        c = p.cells
        pickedCells = c.getSequenceFromMask(mask=('[#f ]', ), )
        v1, e1, d1 = p.vertices, p.edges, p.datums
        p.PartitionCellByPlaneThreePoints(cells=pickedCells, point1=p.InterestingPoint(
            edge=e1[123], rule=MIDDLE), point2=p.InterestingPoint(edge=e1[39], 
            rule=MIDDLE), point3=p.InterestingPoint(edge=e1[126], rule=MIDDLE))
        p = mdb.models['Model-1'].parts['Workpiece']
        c = p.cells
        pickedCells = c.getSequenceFromMask(mask=('[#10 ]', ), )
        e, d = p.edges, p.datums
        pickedEdges =(e[53], )
        p.PartitionCellByExtrudeEdge(line=e[17], cells=pickedCells, edges=pickedEdges, 
            sense=REVERSE)
        p = mdb.models['Model-1'].parts['Workpiece']
        c = p.cells
        pickedCells = c.getSequenceFromMask(mask=('[#10 ]', ), )
        e1, d1 = p.edges, p.datums
        pickedEdges =(e1[74], )
        p.PartitionCellByExtrudeEdge(line=e1[23], cells=pickedCells, edges=pickedEdges, 
            sense=REVERSE)
        p = mdb.models['Model-1'].parts['Workpiece']
        c = p.cells
        pickedCells = c.getSequenceFromMask(mask=('[#200 ]', ), )
        e, v, d = p.edges, p.vertices, p.datums
        p.PartitionCellByPlanePointNormal(point=v[17], normal=e[24], cells=pickedCells)

        p = mdb.models['Model-1'].parts['Workpiece']
        c = p.cells
        pickedCells = c.getSequenceFromMask(mask=('[#8 ]', ), )
        e1, v1, d1 = p.edges, p.vertices, p.datums
        p.PartitionCellByPlanePointNormal(normal=e1[35], cells=pickedCells, 
            point=p.InterestingPoint(edge=e1[101], rule=MIDDLE))
        p = mdb.models['Model-1'].parts['Workpiece']
        c = p.cells
        pickedCells = c.getSequenceFromMask(mask=('[#810 ]', ), )
        e1, d1 = p.edges, p.datums
        pickedEdges =(e1[64], e1[157])
        p.PartitionCellByExtrudeEdge(line=e1[38], cells=pickedCells, edges=pickedEdges, 
            sense=REVERSE)
        p = mdb.models['Model-1'].parts['Workpiece']
        c = p.cells
        pickedCells = c.getSequenceFromMask(mask=('[#3fff ]', ), )
        v, e, d = p.vertices, p.edges, p.datums
        p.PartitionCellByPlaneThreePoints(cells=pickedCells, point1=p.InterestingPoint(
            edge=e[201], rule=MIDDLE), point2=p.InterestingPoint(edge=e[172], 
            rule=MIDDLE), point3=p.InterestingPoint(edge=e[199], rule=MIDDLE))
        p = mdb.models['Model-1'].parts['Workpiece']
        c = p.cells
        pickedCells = c.getSequenceFromMask(mask=('[#dd05fc5 ]', ), )
        v1, e1, d1 = p.vertices, p.edges, p.datums
        p.PartitionCellByPlaneThreePoints(cells=pickedCells, point1=p.InterestingPoint(
            edge=e1[171], rule=MIDDLE), point2=p.InterestingPoint(edge=e1[310], 
            rule=MIDDLE), point3=p.InterestingPoint(edge=e1[169], rule=MIDDLE))

        
        p = mdb.models['Model-1'].parts['Workpiece']
        c = p.cells
        pickedCells = c.getSequenceFromMask(mask=('[#16207f1f #204 ]', ), )
        v1, e1, d1 = p.vertices, p.edges, p.datums
        p.PartitionCellByPlaneThreePoints(cells=pickedCells, point1=p.InterestingPoint(
            edge=e1[89], rule=MIDDLE), point2=p.InterestingPoint(edge=e1[260], 
            rule=MIDDLE), point3=p.InterestingPoint(edge=e1[85], rule=MIDDLE))
        p = mdb.models['Model-1'].parts['Workpiece']
        c = p.cells
        pickedCells = c.getSequenceFromMask(mask=('[#b2004011 #a3fa0b ]', ), )
        v1, e1, d1 = p.vertices, p.edges, p.datums
        p.PartitionCellByPlaneThreePoints(cells=pickedCells, point1=p.InterestingPoint(
            edge=e1[396], rule=MIDDLE), point2=p.InterestingPoint(edge=e1[389], 
            rule=MIDDLE), point3=p.InterestingPoint(edge=e1[395], rule=MIDDLE))
        p = mdb.models['Model-1'].parts['Workpiece']
        c = p.cells
        pickedCells = c.getSequenceFromMask(mask=('[#10046004 #401f148e #37 ]', ), )
        v, e, d = p.vertices, p.edges, p.datums
        p.PartitionCellByPlaneThreePoints(cells=pickedCells, point1=p.InterestingPoint(
            edge=e[404], rule=MIDDLE), point2=p.InterestingPoint(edge=e[523], 
            rule=MIDDLE), point3=p.InterestingPoint(edge=e[402], rule=MIDDLE))
        # exit()

    def mesh(self, size, dev_factor, min_size_factor):
        p = mdb.models['Model-1'].parts['Workpiece']
        p.seedPart(size=0.0031, deviationFactor=0.1, minSizeFactor=0.1)
        p = mdb.models['Model-1'].parts['Workpiece']
        p.generateMesh()

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
        sketch = model.ConstrainedSketch(name='-profile', sheetSize=0.1) 
        sketch.sketchOptions.setValues(decimalPlaces=3)
        sketch.setPrimaryObject(option=STANDALONE)
        sketch.rectangle(point1=(-0.5 * length, -0.5 * height), point2=(0.5 * length, 0.5 * height))
        self.part = model.Part(name=self.name, dimensionality=THREE_D,  type=DEFORMABLE_BODY) # 
        self.part.BaseSolidExtrude(sketch=sketch, depth=width)
        sketch.unsetPrimaryObject()
    
    def set_section(self, section):
        region = self.part.Set(cells=self.part.cells, name='Material-region')
        self.part.SectionAssignment(region=region, sectionName=section.name, offset=0.0, 
            offsetType=MIDDLE_SURFACE, offsetField='', thicknessAssignment=FROM_SECTION)

    def partition(self):
        point1 = (0, 0, 0)
        point2 = (1, 0, 0)
        point3 = (0, 0, 1)
        self.part.PartitionCellByPlaneThreePoints(cells=self.part.cells, point1=point1, point2=point2, point3=point3)
            
        point1 = (0, 0, 0)
        point2 = (0, 1, 0)
        point3 = (0, 0, 1)
        self.part.PartitionCellByPlaneThreePoints(cells=self.part.cells, point1=point1, point2=point2, point3=point3)

    def mesh(self, size, dev_factor, min_size_factor):
        self.part.setMeshControls(regions=self.part.cells, technique=SWEEP, algorithm=MEDIAL_AXIS)
        self.part.seedPart(size=size, deviationFactor=dev_factor, minSizeFactor=min_size_factor)
        self.part.generateMesh()


class InteractionProperty:

    def __init__(self):
        self.name = "Interaction-Property"
        model.ContactProperty(self.name)
        model.interactionProperties[self.name].TangentialBehavior(formulation=ROUGH)
        model.interactionProperties[self.name].NormalBehavior(
            pressureOverclosure=HARD, allowSeparation=OFF, contactStiffness=DEFAULT, 
            contactStiffnessScaleFactor=1.0, clearanceAtZeroContactPressure=0.0, 
            stiffnessBehavior=LINEAR, constraintEnforcementMethod=PENALTY)

class Step:

    def __init__(self, name, previous='Initial'):
        self.name = name
        model.StaticStep(name=name, previous=previous)


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
        self.a = model.rootAssembly
        self.a.DatumCsysByDefault(CARTESIAN)
        self.a.Instance(name=workpiece.name, part=workpiece.part, dependent=ON)
        self.a.translate(instanceList=('Workpiece', ), vector=(0.188909, 0.0, 0.0))
        self.a.translate(instanceList=('Workpiece', ), vector=(-0.015, 0.0, 0.0))
        self.a.Instance(name=jaw.name, part=jaw.part, dependent=ON)
        self.a.translate(instanceList=('Workpiece', ), vector=(-0.1, 0.0, 0.0))
        # self.a.rotate(instanceList=(jaw.name, ), axisPoint=(0, 0, 0), axisDirection=(1, 0.0, 0.0), angle=-90.0)
        # self.a.translate(instanceList=(jaw.name, ), vector=(0,-jaw.width/2, 0))
        self.a.translate(instanceList=(jaw.name, ), vector=(-jaw.length/2, 0.0, workpiece.outer_radius))
        self.a.RadialInstancePattern(instanceList=(jaw.name, ), point=(0.0, 0.0, 0.0), 
            axis=(1.0, 0.0, 0.0), number=3, totalAngle=360.0)
        self.a.features.changeKey(fromName=jaw.name, toName=jaw.name+'-rad-1')
        self.interactionProperty = InteractionProperty()
        self.jaws = [Assembly.Jaw(i, 'Jaw-rad-'+str(i)) for i in range(1,4)]
        self.CSYSs = dict()
        for a_jaw in self.jaws:
            self._create_interaction(a_jaw, self.interactionProperty)
            self.CSYSs[a_jaw] = self._create_CSYS(a_jaw)
            self._create_jaw_BSs(a_jaw)
        self.step = Step("Step-1")
        for a_jaw in self.jaws:
            self._create_tie(a_jaw)
            model.boundaryConditions['BC-'+a_jaw.name].setValuesInStep(stepName=self.step.name, u1=0.0001)
        self.step_2 = Step("Step-2")
        # self._apply_cutting_force(force_rad=300, force_tan=100, force_axial=-100, dz=0.05, phi=deg(-45))
        # for a_jaw in self.jaws:
        #     self._apply_jaw_force(a_jaw, jaw_force)

    def _create_tie(self, jaw):
        region1=self.a.surfaces[jaw.name+'_master_surf']
        region2=self.a.surfaces[jaw.name+'_slave_surf']
        model.Tie(name=jaw.name + 'tie_constr', master=region1, slave=region2, 
            positionToleranceMethod=COMPUTED, adjust=ON, tieRotations=ON, thickness=ON)


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
        # no need for interaction with tie constraint!
        return
        model.SurfaceToSurfaceContactStd(name='Interaction-'+jaw.name, 
            createStepName='Initial', master=region_jaw, slave=region_workpiece, sliding=FINITE, 
            thickness=ON, interactionProperty=property.name, adjustMethod=NONE, 
            initialClearance=OMIT, datumAxis=None, clearanceRegion=None)


    def _create_jaw_BSs(self, jaw):
        jaw_instance = self.a.instances[jaw.name]
        jaw_faces = jaw_instance.faces.getSequenceFromMask(mask=('[#d5428 ]', ), ) # hack, replace it with findAt
        region = self.a.Set(faces = jaw_faces, name='Jaw_BS_set-'+jaw.name)
        datum = self.a.datums[self.CSYSs[jaw].id]
        model.DisplacementBC(name='BC-'+jaw.name, createStepName='Initial', 
            region=region, u1=UNSET, u2=SET, u3=SET, ur1=UNSET, ur2=UNSET, ur3=UNSET, 
            amplitude=UNSET, distributionType=UNIFORM, fieldName='', localCsys=datum)

    def _apply_jaw_force(self, jaw, value):
        vertices = self.a.instances[jaw.name].vertices
        verts1 = vertices.getSequenceFromMask(mask=('[#1 ]', ), ) # hack, replace it with findAt
        region = self.a.Set(vertices=verts1, name='Jaw-force-region-'+jaw.name)        
        datum = self.a.datums[self.CSYSs[jaw].id]
        model.ConcentratedForce(name='Load-'+jaw.name, createStepName=self.step.name, 
            region=region, cf1=value, distributionType=UNIFORM, field='', localCsys=datum)

    def _apply_cutting_force(self, force_rad, force_tan, force_axial, dz, phi):
        dct = dict((k,v) for (k,v) in [('cf1', force_rad), ('cf2', force_tan), ('cf3', force_axial)] if v !=0)

        x = dz-self.workpiece.length
        x, y, z = rotate(point=(x, 0, self.workpiece.outer_radius), axis=OX, theta=phi)

        import nearestNodeModule
        session.viewports['Viewport: 1'].assemblyDisplay.setValues(mesh=ON)
        session.viewports['Viewport: 1'].assemblyDisplay.meshOptions.setValues(
            meshTechnique=ON)
        pickedSelectedNodes = self.a.instances[workpiece.name].nodes
        force_node_info = nearestNodeModule.findNearestNode(xcoord=x, ycoord=y, zcoord=z, name='', 
            selectedNodes=pickedSelectedNodes, instanceName="'Workpiece'")
        nearestNodeModule.hideTextAndArrow()
        node_index = force_node_info[0]
        force_node = self.a.instances[workpiece.name].nodes.sequenceFromLabels((node_index,))

        nf_x, nf_y, nf_z = force_node[0].coordinates
        
        csys = self.a.DatumCsysByThreePoints(origin=(nf_x, nf_y, nf_z), 
            point1=(x, 0, 0), 
            point2=(nf_x-0.01, nf_y, nf_z), 
            name = 'force_CSYS', coordSysType=CARTESIAN)
        
        region = self.a.Set(nodes= force_node,
         name='Force-nodes')
        
        datum = self.a.datums[csys.id]
        model.ConcentratedForce(name='cutting_force', createStepName=self.step_2.name, 
            region=region, distributionType=UNIFORM, field='', localCsys=datum, **dct)
        

def run_job():
    Job1 = mdb.Job(name='Job-1', model='Model-1', description='', type=ANALYSIS, 
        atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, 
        memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
        explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
        modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', 
        scratch='', resultsFormat=ODB, multiprocessingMode=DEFAULT, numCpus=3, 
        numDomains=3, numGPUs=0)
    Job1.submit(consistencyChecking=OFF)
    Job1.waitForCompletion()


if __name__ =="__main__":
    # create materials
    steel = Material_Explicit("Steel", density= 7870, young=2e11, poisson=0.29, 
        A=375e6, B=552e6, n=0.457, 
        d1=0.25, d2=4.38, d3=2.68, ref_strain_rate=1.0, disp_at_failure=0.1)
    steel = Material('Steel', 210e9, 0.29)
    steel_section = model.HomogeneousSolidSection(name='Steel section', material=steel.name, thickness=None)
    
    alu = Material_Explicit("Alu", density= 2700, young=70e9, poisson=0.33, 
        A=3.241e8, B=1.138e8, n=0.42, 
        d1=-0.77, d2=1.45, d3=-0.47, ref_strain_rate=1.0, disp_at_failure=1e-4)
    alu = Material('Alu', 70e9, 0.28)
    # alu = Material('Alu', 100e15, 0.28)
    alu_section = model.HomogeneousSolidSection(name='Alu section', material=alu.name, thickness=None)

    # create workpiece
    workpiece = Workpiece(cd+'\step_models\part4.STEP', outer_radius=99.746E-03/2, length=0.200)
    workpiece.partition(6)
    workpiece.mesh(size=0.006, dev_factor=0.1, min_size_factor = 0.1)
    workpiece.set_section(alu_section)
    session.viewports['Viewport: 1'].setValues(displayedObject=workpiece.part)

    #  create jaw
    jaw = Jaw(length=0.02, width=0.02, height=0.02)
    jaw.partition()
    jaw.mesh(size=0.003, dev_factor=0.1, min_size_factor=0.1)
    jaw.set_section(steel_section)
    session.viewports['Viewport: 1'].setValues(displayedObject=jaw.part)

    # assembly
    assembly = Assembly(workpiece, jaw, jaw_force=400)
    session.viewports['Viewport: 1'].setValues(displayedObject=assembly.a)
    session.viewports['Viewport: 1'].view.setValues(session.views['Iso'])

    run_job()