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

    def mesh(self, size=0.0015, deviationFactor=0.1, minSizeFactor=0.1):
        self.part.setMeshControls(regions=self.part.cells, technique=SWEEP, algorithm=MEDIAL_AXIS)
        self.part.seedPart(size=size, deviationFactor=deviationFactor, minSizeFactor=minSizeFactor)
        self.part.generateMesh()

    def partition(self):
        for p in range(0, self.p_num):
            try:
                self.part.PartitionCellByPlaneThreePoints(cells=self.part.cells, 
                    point1=(0, 0, 0), 
                    point2=(0, 0, 1),
                    point3=rotate(point=(0, 1, 0), axis=OZ, theta = deg(p * 360/self.p_num)))
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

    def mesh(self, size=0.0015, deviationFactor=0.1, minSizeFactor=0.1):
        self.part.setMeshControls(regions=self.part.cells, technique=SWEEP, algorithm=MEDIAL_AXIS)
        self.part.seedPart(size=size, deviationFactor=deviationFactor, minSizeFactor=minSizeFactor)
        self.part.generateMesh()


def __create_workpiece_partion():
    p = mdb.models['Model-1'].parts['Workpiece']
    d1 = p.DatumAxisByPrincipalAxis(principalAxis=ZAXIS)
    v, e1 = p.vertices, p.edges
    d2 = p.DatumPlaneByThreePoints(point1=v[0], point3=v[1], point2=p.InterestingPoint(edge=e1[0], rule=CENTER))
    # d2 = p.DatumPlaneByThreePoints(point2=(0,0,0), point3=(0,0,0.1), point1=(0,0.1,0.1))
    d = p.datums
    for j in JAWS:
        jaw = 'Jaw-1-rad-'+str(j)
        # JAWS=1,2,3; n = 0,1,2
        n = j - 1
        new_d = p.DatumPlaneByRotation(plane=d[d2.id], axis=d[d1.id], angle=n*360/len(JAWS))
        try:
            p.PartitionCellByDatumPlane(datumPlane=d[new_d.id], cells=p.cells)
        except:
            pass
    
def __create_force_point(parameter):    
    p = mdb.models['Model-1'].parts['Workpiece']
    pickedEdges = mdb.models['Model-1'].parts['Workpiece'].edges.findAt((0,outer,length/2))
    p.PartitionEdgeByParam(edges=pickedEdges, parameter=parameter)
    
def create_jaw_part(length, height, width):
    """Creates a jaw part"""
    s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=0.1) # 
    g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints #
    s.sketchOptions.setValues(decimalPlaces=3)
    s.setPrimaryObject(option=STANDALONE)
    s.rectangle(point1=(0.0, 0.0), point2=(length, height)) #
    p = mdb.models['Model-1'].Part(name='Jaw', dimensionality=THREE_D,  type=DEFORMABLE_BODY) # 
    p = mdb.models['Model-1'].parts['Jaw']
    p.BaseSolidExtrude(sketch=s, depth=width) #
    s.unsetPrimaryObject()
    p = mdb.models['Model-1'].parts['Jaw']
    del mdb.models['Model-1'].sketches['__profile__'] #
    __create_jaw_partion(p.name, length, height, width)
    sys.__stdout__.write("Create jaw part "+" Length:"+str(length)+" Height:"+str(height)+" Width:"+str(width)+"\n")
    return p.name
    
def __create_jaw_partion(part_name, length, height, width):
    p = mdb.models['Model-1'].parts[part_name]
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
    

def create_section(section_name,material_name, type, thickness=None):
    if type == 'solid':
        mdb.models['Model-1'].HomogeneousSolidSection(name=section_name, material=material_name, thickness=None)
        return section_name


def create_material(material_name,young,poisson):
    """Creates material with specified name, Young's modulus and Poisson ratio"""
    mdb.models['Model-1'].Material(name=material_name)
    mdb.models['Model-1'].materials[material_name].Elastic(table=((young, poisson), ))
    return material_name

def mesh_part(part_name, size, dev_factor, min_size_factor):
    p = mdb.models['Model-1'].parts[part_name]
    p.setMeshControls(regions=p.cells, technique=SWEEP, algorithm=MEDIAL_AXIS)
    p.seedPart(size=size, deviationFactor=dev_factor, minSizeFactor=min_size_factor)
    p.generateMesh()

def create_property(friction_coeff):
    mdb.models['Model-1'].ContactProperty('IntProp-1')
    mdb.models['Model-1'].interactionProperties['IntProp-1'].TangentialBehavior(formulation=ROUGH)
    mdb.models['Model-1'].interactionProperties['IntProp-1'].NormalBehavior(
        pressureOverclosure=HARD, allowSeparation=OFF, contactStiffness=DEFAULT, 
        contactStiffnessScaleFactor=1.0, clearanceAtZeroContactPressure=0.0, 
        stiffnessBehavior=LINEAR, constraintEnforcementMethod=PENALTY)

def create_assembly(jawf, outer, length, tanf, radf, axlf, angle):
    create_property(f_coeff)# pass friction coeff
    a = mdb.models['Model-1'].rootAssembly
    a.DatumCsysByDefault(CARTESIAN)
    p = mdb.models['Model-1'].parts['Workpiece']
    a.Instance(name='Part-1', part=p, dependent=ON)
    p = mdb.models['Model-1'].parts['Jaw']
    a.Instance(name='Jaw-1', part=p, dependent=ON)
    a.rotate(instanceList=('Jaw-1', ), axisPoint=(0, 0, 0), axisDirection=OX, angle=rad(-90.0))
    a.translate(instanceList=('Jaw-1', ), vector=(0, outer, 0.5 * jaw_height))
    a.RadialInstancePattern(instanceList=('Jaw-1', ), point=(0, 0, 0), axis=OZ, number=3, totalAngle=-360)
    mdb.models['Model-1'].rootAssembly.features.changeKey(fromName='Jaw-1', toName='Jaw-1-rad-1')
    sys.__stdout__.write("Create assembly"+"\n") 
    # a.rotate(instanceList=('Part-1', ), axisPoint=(0.0, 0.0, 0.0), axisDirection=(0.0, 0.0, 1.0), angle = angle)
    create_interaction()
    c_systems = create_CSYS()
    create_step()
    create_jaw_BSs(c_systems)
    apply_jaw_force(jawf, c_systems)
    __create_workpiece_partion()
    mdb.models['Model-1'].parts['Workpiece'].generateMesh()
    a.regenerate()
    
def reset_force_node():
    x,y = pol2cart(outer,anglef)
    z = length * fparameter
    print "z="+str(z)
    a = mdb.models['Model-1'].rootAssembly
    
    node = min(a.instances['Part-1'].nodes, key = lambda n: (n.coordinates[0]-x)**2 + (n.coordinates[1]-y)**2+(n.coordinates[2]-z)**2)
    nodes1 = a.instances['Part-1'].nodes.sequenceFromLabels((node.label,))
    a.Set(nodes=(nodes1,), name='Force-nodes')


def apply_jaw_force(value, c_systems):
    for j in JAWS:
        jaw = 'Jaw-1-rad-'+str(j)
        # JAWS=1,2,3; csys_n = 0,1,2
        csys_n = (-j+2) %3
        a = mdb.models['Model-1'].rootAssembly
        v1 = a.instances[jaw].vertices
        verts1 = v1.getSequenceFromMask(mask=('[#1 ]', ), )
        region = a.Set(vertices=verts1, name=jaw+'force_set')
        datum = mdb.models['Model-1'].rootAssembly.datums[c_systems[csys_n].id]
        mdb.models['Model-1'].ConcentratedForce(name=jaw + 'load', createStepName='Step-1', 
            region=region, cf1=value, distributionType=UNIFORM, field='', localCsys=datum)
    

def create_interaction():
    z = jaw_length/100
    # mid_d = float(inner-outer)/2+inner
    for j in JAWS:
        jaw = 'Jaw-1-rad-'+str(j)
        angle = (j-1)*(360.0/len(JAWS))+90
        
        # select areas on Jaw
        a = mdb.models['Model-1'].rootAssembly
        s1 = a.instances[jaw].faces
        jaw_instance = a.instances[jaw]
        workpiece = a.instances['Part-1']
        jaw_angle = 360/3 * (1-j)
        p1 = rotate((0.25 * jaw_length, outer, 0.25 * jaw_height), OZ, deg(jaw_angle))
        p2 = rotate((-0.25 * jaw_length, outer, 0.25 * jaw_height), OZ, deg(jaw_angle))
        p3 = rotate((0.25 * jaw_length, outer, 0.75 * jaw_height), OZ, deg(jaw_angle))
        p4 = rotate((-0.25 * jaw_length, outer, 0.75 * jaw_height), OZ, deg(jaw_angle))

        jaw_faces = jaw_instance.faces.findAt( (p1,),(p2,), (p3,), (p4,), )
        region1=a.Surface(side1Faces=jaw_faces, name=jaw + '_master_surf')
        
        # select areas on Part
        a = mdb.models['Model-1'].rootAssembly
        s1 = a.instances['Part-1'].faces
        
        def translate_to_workpiece(point):
            x, y, z = point
            rho, phi = cart2pol(x,y)
            x_w, y_w = pol2cart(outer, phi)
            return x_w, y_w, z

        workpiece_faces = workpiece.faces.findAt( (translate_to_workpiece(p1),),(translate_to_workpiece(p2),), (translate_to_workpiece(p3),), (translate_to_workpiece(p4),), )
        region2=a.Surface(side1Faces=workpiece_faces, name='Part_slave_surf' + str(j))
        
        # apply interaction to Jaw and Part
        mdb.models['Model-1'].SurfaceToSurfaceContactStd(name='Int-'+str(j), 
            createStepName='Initial', master=region1, slave=region2, sliding=FINITE, 
            thickness=ON, interactionProperty='IntProp-1', adjustMethod=NONE, 
            initialClearance=OMIT, datumAxis=None, clearanceRegion=None)

def create_CSYS():
    res=[]
    for j in JAWS:
        jaw = 'Jaw-1-rad-'+str(j)
        csys_n = j+1
        a = mdb.models['Model-1'].rootAssembly
        d1 = a.instances[jaw].datums
        v1 = a.instances[jaw].vertices
        jaw_angle = csys_n *360/3
        res.append(a.DatumCsysByThreePoints(origin=rotate((0, outer, 0), OZ, deg(jaw_angle)), 
            point1=(0,0,0), 
            point2=rotate((1, outer, 0), OZ, deg(jaw_angle)), 
            name=jaw + '_csys', coordSysType=CARTESIAN))
    return res

def create_jaw_BSs(c_systems):
    for j in JAWS:
        jaw = 'Jaw-1-rad-'+str(j)
        #starts from 0 cause c_systems is passed to function
        # JAWS=1,2,3; csys_n = 0,1,2
        csys_n = (-j+2) %3
        a = mdb.models['Model-1'].rootAssembly
        f1 = a.instances['Jaw-1-rad-'+str(j)].faces
        #faces1 = f1.getSequenceFromMask(mask=('[#402 ]', ), )
        x1, y1 = pol2cart(outer+jaw_width/2, 360/len(JAWS)*csys_n+0.5)
        x2, y2 = pol2cart(outer+jaw_width/2, 360/len(JAWS)*csys_n-0.5)
        
        # faces1 = f1.findAt((x1,y1, 0),(x2,y2, 0))
        # sys.__stdout__.write( str(faces1))
        faces1 = f1.getSequenceFromMask(mask=('[#402 ]', ), )
        # faces1 = f1.getSequenceFromMask(mask=('[#8200 ]', ), )
        region = a.Set(faces = faces1, name='Jaw_BS_set-'+str(j))
        datum = mdb.models['Model-1'].rootAssembly.datums[c_systems[csys_n].id]
        mdb.models['Model-1'].DisplacementBC(name='BC-'+str(j), createStepName='Initial', 
            region=region, u1=UNSET, u2=SET, u3=SET, ur1=UNSET, ur2=UNSET, ur3=UNSET, 
            amplitude=UNSET, distributionType=UNIFORM, fieldName='', localCsys=datum)

            
def create_step():
    mdb.models['Model-1'].StaticStep(name='Step-1', previous='Initial')


def assign_section(part_name, section_name):
    p = mdb.models['Model-1'].parts[part_name]
    if p.cells:  # if solid
        region = p.Set(cells=p.cells, name=part_name+'_material_region')
    else:        # if shell
        region = p.Set(faces=p.faces, name=part_name+'_material_region')
    p.SectionAssignment(region=region, sectionName=section_name, offset=0.0, 
                offsetType=TOP_SURFACE, offsetField='', 
                thicknessAssignment=FROM_SECTION)          

def run_job(home):
    sys.__stdout__.write("Run job"+"\n")
    Job1 = mdb.Job(name='Job-1', model='Model-1', description='', type=ANALYSIS, 
    atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, 
    memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
    explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
    modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', 
    scratch=home, resultsFormat=ODB, multiprocessingMode=DEFAULT, numCpus=1, 
    numGPUs=0)
    Job1.submit(consistencyChecking=OFF)
    Job1.waitForCompletion()





length, outer_diameter, inner_diameter = 0.060, 0.068, 0.059
outer, inner = outer_diameter/2, inner_diameter/2
anglef = 45
fparameter = 0.6 # axial position of force (?)
jcount = 3
JAWS = range(1,jcount+1)
jaw_length, jaw_height, jaw_width= 0.015, 0.015, 0.015
jaw_young, jaw_poisson = 210e15, 0.29
jaw_mesh_size = 0.0015
jawf = 1000
refresh = False
f_coeff = 0.3

workpiece_young, workpiece_poisson = 0.7e9, 0.28
workpiece_mesh_size = 0.002

tanf,radf,axlf = 200, 300, 100

rplane = 0
tplane = 1

if __name__ == "__main__":
    
    home = 'C:\\Program Files\\SIMULIA\\Workspace'
    os.chdir(home)

    # workpiece_part=create_workpiece_part(length, inner, outer, anglef, fparameter)
    workpiece_part=Workpiece(length, inner, outer, p_num=3)
    jaw_part = Jaw(jaw_length, jaw_width, jaw_height)
    jaw_part.partition()
    # jaw_part = create_jaw_part(jaw_width, jaw_height, jaw_length)

    steel = Material('Steel', jaw_young, jaw_poisson)
    steel_section = model.HomogeneousSolidSection(name='Steel-section', material=steel.name, thickness=None)

    aluminum = Material('Aluminum', workpiece_young, workpiece_poisson)
    aluminum_section = model.HomogeneousSolidSection(name='Aluminum-section', material=aluminum.name, thickness=None)
    
    mesh_part(part_name = jaw_part.name, size = jaw_mesh_size, dev_factor=0.1, min_size_factor = 0.1 )
    mesh_part(part_name = workpiece_part.name, size = workpiece_mesh_size, dev_factor=0.1, min_size_factor = 0.1 )

    assign_section(part_name = workpiece_part.name, section_name = aluminum_section.name)
    assign_section(part_name = jaw_part.name, section_name = steel_section.name)
    
    create_assembly(jawf, outer, length, tanf, radf, axlf, anglef)
    a = mdb.models['Model-1'].rootAssembly
    session.viewports['Viewport: 1'].setValues(displayedObject=a)
    session.viewports['Viewport: 1'].assemblyDisplay.setValues(optimizationTasks=OFF, geometricRestrictions=OFF, stopConditions=OFF)
    session.viewports['Viewport: 1'].odbDisplay.basicOptions.setValues(renderShellThickness=ON)
    session.writeOBJFile(fileName=home+"/tmp_model.obj",canvasObjects= (session.viewports['Viewport: 1'], ))
    # COMMENT_THIS_LINE
    run_job(home)
