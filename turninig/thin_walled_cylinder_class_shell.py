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

def N(value):
    return value;

def m(value):
    return value

def mm(value): # gets mm, returns m
    return value * 1e-3

def cm(value):
    return value * 1e-2

def rad(rad_value): # gets radians, returns radians
    return rad_value

def deg(value): # gets degrees, return radians
    return (value * pi)/180

def cart2pol(x, y):
    rho = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y, x)
    return(rho, math.degrees(phi))

def pol2cart(rho, phi):
    phi = math.radians(phi)
    x = rho * np.cos(phi)
    y = rho * np.sin(phi)
    return(x, y)   

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
    


class MaterialExplicit:

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
        self.material =  model.materials[name]

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



class Workpiece:

    def __init__(self, length, inner, outer, p_num):
        self.name = "Workpiece"
        self.length = length
        self.inner = inner
        self.outer = outer
        self.p_num = p_num

        sketch = model.ConstrainedSketch(name=self.name + '-profile', sheetSize=0.1)
        sketch.sketchOptions.setValues(decimalPlaces=3)
        sketch.setPrimaryObject(option=STANDALONE)
        sketch.CircleByCenterPerimeter(center=(0.0, 0.0), point1=(0.0, self.inner))
        sketch.CircleByCenterPerimeter(center=(0.0, 0.0), point1=(0.0, self.outer))
        model.Part(name=self.name, dimensionality=THREE_D, type=DEFORMABLE_BODY)
        self.part = model.parts[self.name]
        self.part.BaseSolidExtrude(sketch=sketch, depth=self.length)
        sketch.unsetPrimaryObject()

    def set_section(self, section):
        region = self.part.Set(cells=self.part.cells, name='Material-region')
        self.part.SectionAssignment(region=region, sectionName=section.name, offset=0.0, 
            offsetType=MIDDLE_SURFACE, offsetField='', thicknessAssignment=FROM_SECTION)

    def mesh(self, size, dev_factor, min_size_factor):
        self.part.setMeshControls(regions=self.part.cells, technique=SWEEP, algorithm=MEDIAL_AXIS)
        self.part.seedPart(size=size, deviationFactor=dev_factor, minSizeFactor=min_size_factor)
        self.part.generateMesh()

    def partition(self, p_num):
        for p in range(0, p_num):
            try:
                self.part.PartitionCellByPlaneThreePoints(cells=self.part.cells, 
                    point1=(0, 0, 0), 
                    point2=(0, 0, 1),
                    point3=rotate(point=(0, 1, 0), axis=OZ, theta = deg(p * 360/p_num)))
            except:
                pass


class SketchWorkpiece():

    def __init__(self, length, inner, outer, p_num):
        self.name = "Workpiece"
        self.length = length
        self.inner = inner
        self.outer = outer
        self.p_num = p_num

        s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=0.1)
        g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
        s.sketchOptions.setValues(decimalPlaces=3)
        s.setPrimaryObject(option=STANDALONE)
        s.ConstructionLine(point1=(0.0, -0.05), point2=(0.0, 0.05))
        s.FixedConstraint(entity=g[2])
        session.viewports['Viewport: 1'].view.setValues(nearPlane=0.082073, 
            farPlane=0.106489, cameraPosition=(0.0278429, 0.00694274, 0.0942809), 
            cameraTarget=(0.0278429, 0.00694274, 0))
        
        points=[(0.030, 0), (0.033, 0), (0.033, 0.025), (0.040, 0.025), (0.040, 0.028),(0.030, 0.028), (0.030, 0)]
        for i,point in enumerate(points[:-1]):
            s.Line(point1=points[i], point2=points[i+1])
            
        p = mdb.models['Model-1'].Part(name='Part-3', dimensionality=THREE_D, 
            type=DEFORMABLE_BODY)
        p = mdb.models['Model-1'].parts['Part-3']
        p.BaseSolidRevolve(sketch=s, angle=360.0, flipRevolveDirection=OFF)
        s.unsetPrimaryObject()
        p = mdb.models['Model-1'].parts['Part-3']
        self.part = p
        session.viewports['Viewport: 1'].setValues(displayedObject=p)
        del mdb.models['Model-1'].sketches['__profile__']


    def set_section(self, section):
        region = self.part.Set(cells=self.part.cells, name='Material-region')
        self.part.SectionAssignment(region=region, sectionName=section.name, offset=0.0, 
            offsetType=MIDDLE_SURFACE, offsetField='', thicknessAssignment=FROM_SECTION)

    def mesh(self, size, dev_factor, min_size_factor):
        self.part.setMeshControls(regions=self.part.cells, technique=SWEEP, algorithm=MEDIAL_AXIS)
        self.part.seedPart(size=size, deviationFactor=dev_factor, minSizeFactor=min_size_factor)
        self.part.generateMesh()

    def partition(self, p_num):
        for p in range(0, p_num):
            try:
                self.part.PartitionCellByPlaneThreePoints(cells=self.part.cells, 
                    point1=(0, 0, 0), 
                    point2=(0, 1, 0),
                    point3=rotate(point=(1, 0, 0), axis=OY, theta = deg(30+p * 360/p_num)))
            except:
                pass



class CustomWorkpiece:

    def __init__(self, path, outer_radius):
        self.name = 'Workpiece'
        self.outer = outer_radius
        step = mdb.openStep(path, scaleFromFile=OFF)
        model.PartFromGeometryFile(name=self.name, geometryFile=step, 
            combine=False, dimensionality=THREE_D, type=DEFORMABLE_BODY, scale=0.001)
        self.part = model.parts[self.name]
        self.partition(3)

    def partition(self, number):
        for angle in np.arange(-30, -30+360, 360/number):
            try:
                self.part.PartitionCellByPlaneThreePoints(cells=self.part.cells, point1=(0,0,0), point2=(1,0,0), point3=(0,) +pol2cart(1, 30+angle))
            except:
                print("One partition failed")
        
    def mesh(self, size, dev_factor, min_size_factor):
        self.part.setMeshControls(regions=self.part.cells, technique=SWEEP, algorithm=MEDIAL_AXIS)
        self.part.seedPart(size=size, deviationFactor=dev_factor, minSizeFactor=min_size_factor)
        self.part.generateMesh()

    def set_section(self, section):
        region = self.part.Set(cells=self.part.cells, name='Material-region')
        self.part.SectionAssignment(region=region, sectionName=section.name, offset=0.0, 
            offsetType=MIDDLE_SURFACE, offsetField='', thicknessAssignment=FROM_SECTION)


class ShellWorkpiece:

    def __init__(self, length, inner, outer, p_num):
        self.name = "Workpiece"
        self.length = length
        self.inner = inner
        self.outer = outer
        self.p_num = p_num
        self.thickness = (outer-inner)/2

        sketch = model.ConstrainedSketch(name=self.name + '-profile', sheetSize=0.1)
        sketch.sketchOptions.setValues(decimalPlaces=3)
        sketch.setPrimaryObject(option=STANDALONE)
        sketch.CircleByCenterPerimeter(center=(0.0, 0.0), point1=(0.0, self.outer))
        model.Part(name=self.name, dimensionality=THREE_D, type=DEFORMABLE_BODY)
        self.part = model.parts[self.name]
        self.part.BaseShellExtrude(sketch=sketch, depth=self.length)
        sketch.unsetPrimaryObject()


    def set_section(self, section):
        region = self.part.Set(faces=self.part.faces, name='Material-region')
        self.part.SectionAssignment(region=region, sectionName='Aluminum-shell', offset=0.0, 
            offsetType=MIDDLE_SURFACE, offsetField='', 
            thicknessAssignment=FROM_SECTION)

        # region = self.part.Set(faces=self.part.faces, name='Material-region')
        # self.part.SectionAssignment(region=region, sectionName='Aluminum-shell', offset=0.0, 
        #     offsetType=TOP_SURFACE, offsetField='', thicknessAssignment=FROM_SECTION)


    def mesh(self, size, dev_factor, min_size_factor):
        self.part.setMeshControls(regions=self.part.cells, technique=SWEEP, algorithm=MEDIAL_AXIS)
        self.part.seedPart(size=size, deviationFactor=dev_factor, minSizeFactor=min_size_factor)
        self.part.generateMesh()

    def partition(self, p_num):
        for p in range(0, p_num):
            try:
                self.part.PartitionCellByPlaneThreePoints(cells=self.part.cells, 
                    point1=(0, 0, 0), 
                    point2=(0, 1, 0),
                    point3=rotate(point=(1, 0, 0), axis=OY, theta = deg(30+p * 360/p_num)))
            except:
                pass


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



class Assembly:
    
    class Jaw:

        def __init__(self, index, name):
            self.index = index
            self.name = name
            self.angle = 360/3 * (1-index)
            self.CSYS = None

        def __hash__(self):
            return self.index

    def __init__(self, workpiece, jaw, jaw_num, jaw_force):
        self.workpiece = workpiece
        self.jaw = jaw
        self.jaw_force = jaw_force
        self.a = model.rootAssembly    
        self.a.DatumCsysByDefault(CARTESIAN)
        self.a.Instance(name=workpiece.name, part=workpiece.part, dependent=ON)
        
        # self.a.rotate(instanceList=(workpiece.name, ), axisPoint=(0, 0, 0), axisDirection=OX, angle=rad(+90.0))
        # self.a.translate(instanceList=('Workpiece', ), vector=(0.0, 0.0, 0.015-0.1))

        self.a.Instance(name=jaw.name, part=jaw.part, dependent=ON)
        self.a.rotate(instanceList=(jaw.name, ), axisPoint=(0, 0, 0), axisDirection=OX, angle=rad(-90.0))
        self.a.translate(instanceList=(jaw.name, ), vector=(0, workpiece.outer, 0.5 * jaw.height))
        self.a.RadialInstancePattern(instanceList=(jaw.name, ), point=(0, 0, 0), axis=OZ, number=jaw_num, totalAngle=-360)
        self.a.features.changeKey(fromName=jaw.name, toName=jaw.name+'-rad-1')
        self.assembly_jaws = [Assembly.Jaw(i, 'Jaw-rad-'+str(i)) for i in range(1, 1+jaw_num)]
        self.interactionProperty = InteractionProperty()
        
        for a_jaw in self.assembly_jaws:
            self._create_CSYS(a_jaw)
            self._create_interaction(a_jaw)
            self._create_jaw_BSs(a_jaw)
        self.step = Step("Step-1")
        for a_jaw in self.assembly_jaws:
            self._apply_jaw_force(a_jaw, jaw_force)

        self.workpiece.partition(jaw_num)
        # self.workpiece.mesh(size=0.002, dev_factor=0.1, min_size_factor=0.1)
        self.workpiece.part.generateMesh()
        self.a.regenerate()


    def _create_CSYS(self, jaw):    
        jaw.CSYS = self.a.DatumCsysByThreePoints(origin=rotate((0, workpiece.outer, 0), OZ, deg(jaw.angle)), 
            point1=(0,0,0), 
            point2=rotate((1, workpiece.outer, 0), OZ, deg(jaw.angle)), 
            name=jaw.name + '_CSYS', coordSysType=CARTESIAN)

    def _create_interaction(self, jaw):
        jaw_instance = self.a.instances[jaw.name]
        workpiece = self.a.instances[self.workpiece.name]
        p1 = rotate((0.25 * self.jaw.length, self.workpiece.outer, 0.25 * self.jaw.height), OZ, deg(jaw.angle))
        p2 = rotate((-0.25 * self.jaw.length, self.workpiece.outer, 0.25 * self.jaw.height), OZ, deg(jaw.angle))
        p3 = rotate((0.25 * self.jaw.length, self.workpiece.outer, 0.75 * self.jaw.height), OZ, deg(jaw.angle))
        p4 = rotate((-0.25 * self.jaw.length, self.workpiece.outer, 0.75 * self.jaw.height), OZ, deg(jaw.angle))

        jaw_faces = jaw_instance.faces.findAt( (p1,),(p2,), (p3,), (p4,), )
        
        region_jaw=self.a.Surface(side1Faces=jaw_faces, name=jaw.name + '_master_surf')       
        
        # create workpiece master region   
        def translate_to_workpiece(point):
            x, y, z = point
            rho, phi = cart2pol(x,y)
            x_w, y_w = pol2cart(self.workpiece.outer, phi)
            return x_w, y_w, z

        workpiece_faces = workpiece.faces.findAt( (translate_to_workpiece(p1),),(translate_to_workpiece(p2),), (translate_to_workpiece(p3),), (translate_to_workpiece(p4),), )
        region_workpiece=self.a.Surface(side1Faces=workpiece_faces, name=jaw.name + '_slave_surf')

        model.SurfaceToSurfaceContactStd(name='Interaction-'+jaw.name, 
            createStepName='Initial', master=region_jaw, slave=region_workpiece, sliding=FINITE, 
            thickness=ON, interactionProperty=self.interactionProperty.name, adjustMethod=NONE, 
            initialClearance=OMIT, datumAxis=None, clearanceRegion=None)

    def _create_jaw_BSs(self, jaw):
        jaw_instance = self.a.instances[jaw.name]
        
        jaw_faces = jaw_instance.faces.findAt(
            (rotate((0.25 * self.jaw.length, self.workpiece.outer + 0.5* self.jaw.width, 0), OZ, deg(jaw.angle)),),
            (rotate((-0.25 * self.jaw.length, self.workpiece.outer + 0.5* self.jaw.width, 0), OZ, deg(jaw.angle)),),
            )
        
        region = self.a.Set(faces = jaw_faces, name='Jaw_BS_set-'+jaw.name)
        
        datum = self.a.datums[jaw.CSYS.id]
        model.DisplacementBC(name='BC-'+jaw.name, createStepName='Initial', 
            region=region, u1=UNSET, u2=SET, u3=SET, ur1=UNSET, ur2=UNSET, ur3=UNSET, 
            amplitude=UNSET, distributionType=UNIFORM, fieldName='', localCsys=datum)

    def _apply_jaw_force(self, jaw, value):
        vertices = self.a.instances[jaw.name].vertices
        verts1 = vertices.getSequenceFromMask(mask=('[#1 ]', ), ) # hack, replace it with findAt
        region = self.a.Set(vertices=verts1, name='Jaw-force-region-'+jaw.name)        
        datum = self.a.datums[jaw.CSYS.id]
        model.ConcentratedForce(name='Load-'+jaw.name, createStepName=self.step.name, 
            region=region, cf1=value, distributionType=UNIFORM, field='', localCsys=datum)
   

    def create_force(self, tanf, radf, axlf, angle, fparameter):
        if tanf**2+ radf**2 + axlf**2 == 0:
            return
        d = self.a.datums
        n = self.workpiece.part.nodes
        n = self.a.instances[workpiece.name].nodes
        force_sys = self.a.DatumCsysByTwoLines(CYLINDRICAL, line1=d[1].axis1, line2=d[1].axis2, name='Cutting_force_csys')
        v = self.a.instances[workpiece.name].vertices

        datum = model.rootAssembly.datums[force_sys.id]
        dct = dict((k,v) for (k,v) in [('cf1', radf), ('cf2', tanf), ('cf3', axlf)] if v !=0)
        
        delta = 0.5e-4
        nodes1=[]
        x,y = pol2cart(self.workpiece.outer,angle)
        z = self.workpiece.length * (1-fparameter)
        while not nodes1:
            xmin, ymin, zmin = x-delta, y-delta, z-delta
            xmax, ymax, zmax = x+delta, y+delta, z+delta 
            nodes1 = n.getByBoundingBox(xmin, ymin, zmin, xmax, ymax, zmax)
            delta+=delta

        region = self.a.Set(nodes=nodes1, name='Force-nodes')
        mdb.models['Model-1'].StaticStep(name='Step-2', previous='Step-1')
        model.ConcentratedForce(name='Cutting_force', createStepName='Step-2', 
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

jcount = 3
JAWS = range(1,jcount+1)

if __name__ == "__main__":
    jaw_num = 3
    steel = Material('Steel', 210e15, 0.29)
    steel_section = model.HomogeneousSolidSection(name='Steel-section', material=steel.name, thickness=None)

    aluminum = Material('Aluminum', 0.7e9, 0.28)
    aluminum_section = model.HomogeneousSolidSection(name='Aluminum-section', material=aluminum.name, thickness=None)
    
    workpiece = ShellWorkpiece(length=mm(60), inner=mm(59/2), outer=mm(68/2), p_num=jaw_num)
    # workpiece = Workpiece(length=mm(60), inner=mm(59/2), outer=mm(68/2), p_num=jaw_num)
    # workpiece = CustomWorkpiece('D:/ereme/GoogleDrive/PostGraduate/ZhAD/parts/meshed/part1.STEP', outer_radius=0.046)
    # workpiece = SketchWorkpiece(length=mm(60), inner=mm(59/2), outer=0.033, p_num=jaw_num)
    
    aluminum_shell_section = mdb.models['Model-1'].HomogeneousShellSection(name='Aluminum-shell', 
    preIntegrate=OFF, material='Aluminum', thicknessType=UNIFORM, 
    thickness=workpiece.thickness, thicknessField='', idealization=NO_IDEALIZATION, 
    poissonDefinition=DEFAULT, thicknessModulus=None, temperature=GRADIENT, 
    useDensity=OFF, integrationRule=SIMPSON, numIntPts=5)

    workpiece.set_section(aluminum_shell_section)
    workpiece.partition(3)
    workpiece.mesh(size=0.002, dev_factor=0.1, min_size_factor = 0.1 )

    jaw = Jaw(length=mm(15), width=mm(15), height=mm(15))
    jaw.set_section(steel_section)
    jaw.partition()
    jaw.mesh(size=0.0015, dev_factor=0.1, min_size_factor=0.1)

    assembly = Assembly(workpiece=workpiece, jaw=jaw, jaw_num=jaw_num, jaw_force=N(1000))
    assembly.create_force(300, -400, -500, 30, 0.2)

    a = mdb.models['Model-1'].rootAssembly
    session.viewports['Viewport: 1'].setValues(displayedObject=a)
    session.viewports['Viewport: 1'].view.setValues(session.views['Iso'])

    run_job()

    session.graphicsOptions.setValues(backgroundStyle=SOLID, backgroundColor='#FFFFFF')
    session.viewports['Viewport: 1'].viewportAnnotationOptions.setValues(triad=ON, 
    legend=ON, title=OFF, state=OFF, annotations=OFF, compass=OFF)

    o3 = session.openOdb(name='C:/Program Files/SIMULIA/Temp/Job-1.odb')
    session.viewports['Viewport: 1'].setValues(displayedObject=o3)
    o3 = session.openOdb(name='C:/Program Files/SIMULIA/Temp/Job-1.odb')
    session.viewports['Viewport: 1'].setValues(displayedObject=o3)


    session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(CONTOURS_ON_DEF, ))
    session.viewports['Viewport: 1'].view.setValues(session.views['Iso'])
    session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(variableLabel='U', outputPosition=NODAL, refinement=(INVARIANT, 
        'Magnitude'), )
    session.viewports['Viewport: 1'].odbDisplay.basicOptions.setValues(renderShellThickness=ON)

