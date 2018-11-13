""""
replace force with display
"""
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

def cart2pol(x, y):
    rho = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y, x)
    return(rho, math.degrees(phi))

def pol2cart(rho, phi):
    phi = math.radians(phi)
    x = rho * np.cos(phi)
    y = rho * np.sin(phi)
    return(x, y)   

def create_workpiece_part(length, inner, outer, angle, parameter):   
    s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=0.1)
    g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
    s.sketchOptions.setValues(decimalPlaces=3)
    s.setPrimaryObject(option=STANDALONE)
    s.CircleByCenterPerimeter(center=(0.0, 0.0), point1=(0.0, inner))
    s.CircleByCenterPerimeter(center=(0.0, 0.0), point1=(0.0, outer))
    p = mdb.models['Model-1'].Part(name='Part', dimensionality=THREE_D, type=DEFORMABLE_BODY)
    p = mdb.models['Model-1'].parts['Part']
    p.BaseSolidExtrude(sketch=s, depth=length)
    s.unsetPrimaryObject()
    p = mdb.models['Model-1'].parts['Part']
    session.viewports['Viewport: 1'].setValues(displayedObject=p)
    del mdb.models['Model-1'].sketches['__profile__']
    
    # __create_force_point(parameter)
    sys.__stdout__.write("Create workpiece model"+" Length:"+str(length)+" Inner:"+str(inner)+" Outer:"+str(outer)+"\n")
    return p.name

def __create_workpiece_partion():
    p = mdb.models['Model-1'].parts['Part']
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
    p = mdb.models['Model-1'].parts['Part']
    pickedEdges = mdb.models['Model-1'].parts['Part'].edges.findAt((0,outer,length/2))
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
    p = mdb.models['Model-1'].parts['Part']
    a.Instance(name='Part-1', part=p, dependent=ON)
    p = mdb.models['Model-1'].parts['Jaw']
    a.Instance(name='Jaw-1', part=p, dependent=ON)
    a.translate(instanceList=('Jaw-1', ), vector=(-jaw_width/2, 0.0, 0.0))
    a.rotate(instanceList=('Jaw-1', ), axisPoint=(0.0, 0.0, 0.0), axisDirection=(0.0, 0.0, 1), angle=-90.0)
    mid_d = float(inner-outer)/2+inner
    a.translate(instanceList=('Jaw-1', ), vector=(outer, 0.0, 0.0))
    a.rotate(instanceList=('Jaw-1', ), axisPoint=(0.0, 0.0, 0.0), axisDirection=(0.0, 0.0, 1), angle=90.0)
    a.RadialInstancePattern(instanceList=('Jaw-1', ), point=(0.0, 0.0, 0.0), axis=(0.0, 0.0, 1.0), number=jcount, totalAngle=360.0)
    mdb.models['Model-1'].rootAssembly.features.changeKey(fromName='Jaw-1', toName='Jaw-1-rad-1')
    sys.__stdout__.write("Create assembly"+"\n") 
    # a.rotate(instanceList=('Part-1', ), axisPoint=(0.0, 0.0, 0.0), axisDirection=(0.0, 0.0, 1.0), angle = angle)
    create_interaction()
    c_systems = create_CSYS()
    create_step()
    create_jaw_BSs(c_systems)
    apply_jaw_force(jawf, c_systems)
    __create_workpiece_partion()
    mdb.models['Model-1'].parts['Part'].generateMesh()
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
        csys_n = j - 1
        a = mdb.models['Model-1'].rootAssembly
        v1 = a.instances[jaw].vertices
        verts1 = v1.getSequenceFromMask(mask=('[#2 ]', ), )
        region = a.Set(vertices=verts1, name=jaw+'force_set')
        datum = mdb.models['Model-1'].rootAssembly.datums[c_systems[csys_n].id]
        mdb.models['Model-1'].ConcentratedForce(name=jaw + 'load', createStepName='Step-1', 
            region=region, cf1=value, distributionType=UNIFORM, field='', localCsys=datum)
    

def create_force(tanf, radf, axlf, angle):
    a = mdb.models['Model-1'].rootAssembly
    d = a.datums
    n = a.instances['Part-1'].nodes
    force_sys = a.DatumCsysByTwoLines(CYLINDRICAL, line1=d[1].axis1, line2=d[1].axis2, name='Cutting_force_csys')
    v = a.instances['Part-1'].vertices
    verts = v.findAt((pol2cart(outer,90+angle)[0],pol2cart(outer,90+angle)[1], length*(1-fparameter)))
    a.instances['Part-1'].vertices.getSequenceFromMask(mask=('[#2 ]', ), )
    verts1 = v.getSequenceFromMask(mask=('[#2 ]', ), )
    region = a.Set(vertices=verts1, name='Set-7')
    datum = mdb.models['Model-1'].rootAssembly.datums[force_sys.id]
    dct = dict((k,v) for (k,v) in [('cf1', radf), ('cf2', tanf), ('cf3', axlf)] if v !=0)
        
    delta = 1.0e-4
    nodes1=[]
    x,y = pol2cart(outer,angle)
    z = length * (fparameter-1)
    while not nodes1:
        xmin, ymin, zmin = x-delta, y-delta, z-delta
        xmax, ymax, zmax = x+delta, y+delta, z+delta 
        nodes1 = n.getByBoundingBox(xmin, ymin, zmin, xmax, ymax, zmax)
        delta+=delta
    region = a.Set(nodes=nodes1, name='Force-nodes')
    mdb.models['Model-1'].StaticStep(name='Step-2', previous='Step-1')
    if tanf**2+ radf**2 + axlf**2 !=0:
        mdb.models['Model-1'].ConcentratedForce(name='Cutting_force', createStepName='Step-2', 
            region=region, distributionType=UNIFORM, field='', localCsys=datum, **dct)
    

def __get_xy(radius, ref_angle):
    x1 = radius * math.cos(math.radians(ref_angle+0.05))
    x2 = radius * math.cos(math.radians(ref_angle-0.05))
    y1 = radius * math.sin(math.radians(ref_angle+0.05))
    y2 = radius * math.sin(math.radians(ref_angle-0.05))
    return x1, y1, x2,y2

def create_interaction():
    z = jaw_length/100
    # mid_d = float(inner-outer)/2+inner
    for j in JAWS:
        jaw = 'Jaw-1-rad-'+str(j)
        angle = (j-1)*(360.0/len(JAWS))+90
        x1, y1, x2, y2 = __get_xy(outer, angle)
        
        # select areas on Jaw
        a = mdb.models['Model-1'].rootAssembly
        s1 = a.instances[jaw].faces
        side1Faces1 = s1.findAt(((x1, y1, z), ), ((x2, y2, z), ),((x1, y1, z+jaw_length/2), ),((x2, y2, z+jaw_length/2), ))
        region1=a.Surface(side1Faces=side1Faces1, name=jaw + '_master_surf')
        
        # select areas on Part
        a = mdb.models['Model-1'].rootAssembly
        s1 = a.instances['Part-1'].faces
        side1Faces1 = s1.findAt(((x1, y1, z), ), ((x2, y2, z), ),((x1, y1, z+jaw_length/2), ),((x2, y2, z+jaw_length/2), ))
        region2=a.Surface(side1Faces=side1Faces1, name='Part_slave_surf' + str(j))
        
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
        res.append(a.DatumCsysByThreePoints(origin=d1[2], point1=(0,0,0), point2=v1[7], name = jaw + '_csys-'+str(csys_n), coordSysType=CARTESIAN))
    return res

def create_jaw_BSs(c_systems):
    for j in JAWS:
        jaw = 'Jaw-1-rad-'+str(j)
        #starts from 0 cause c_systems is passed to function
        # JAWS=1,2,3; csys_n = 0,1,2
        csys_n = j - 1
        a = mdb.models['Model-1'].rootAssembly
        f1 = a.instances['Jaw-1-rad-'+str(j)].faces
        #faces1 = f1.getSequenceFromMask(mask=('[#402 ]', ), )
        x1, y1 = pol2cart(outer+jaw_width/2, 360/len(JAWS)*csys_n+0.5)
        x2, y2 = pol2cart(outer+jaw_width/2, 360/len(JAWS)*csys_n-0.5)
        
        # faces1 = f1.findAt((x1,y1, 0),(x2,y2, 0))
        # sys.__stdout__.write( str(faces1))
        faces1 = f1.getSequenceFromMask(mask=('[#402 ]', ), )
        region = a.Set(faces = faces1, name='Jaw_BS_set-'+str(j))
        datum = mdb.models['Model-1'].rootAssembly.datums[c_systems[csys_n].id]
        mdb.models['Model-1'].DisplacementBC(name='BC-'+str(j), createStepName='Initial', 
            region=region, u1=UNSET, u2=SET, u3=SET, ur1=UNSET, ur2=UNSET, ur3=UNSET, 
            amplitude=UNSET, distributionType=UNIFORM, fieldName='', localCsys=datum)

            
def create_step():
    mdb.models['Model-1'].StaticStep(name='Step-1', previous='Initial')

def apply_jaw_displ(value, c_systems):
    z = jaw_length/100
    for j in JAWS:
        jaw = 'Jaw-1-rad-'+str(j)
        # JAWS=1,2,3; csys_n = 0,1,2
        csys_n = j - 1
        
        angle = (j-1)*(360.0/len(JAWS))+90
        x1, y1, x2, y2 = __get_xy(outer+jaw_height, angle)
        
        # select areas on Jaw
        a = mdb.models['Model-1'].rootAssembly
        s1 = a.instances[jaw].faces
        side1Faces1 = s1.findAt(((x1, y1, z), ), ((x2, y2, z), ),((x1, y1, z+jaw_length/2), ),((x2, y2, z+jaw_length/2), ))
        region=a.Set(faces=side1Faces1, name=jaw + '_displ_surf')
        
        datum = mdb.models['Model-1'].rootAssembly.datums[c_systems[csys_n].id]
        mdb.models['Model-1'].DisplacementBC(name='DISPL-BC-'+str(j), createStepName='Step-1', region=region, u1=value, u2=UNSET, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET, amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='', localCsys=datum)


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
    
Mdb()

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

    workpiece_part=create_workpiece_part(length, inner, outer, anglef, fparameter)
    jaw_part = create_jaw_part(jaw_width, jaw_height, jaw_length)
    j_material = create_material('Jaw_material', jaw_young, jaw_poisson)
    w_material = create_material('Workpiece_material', workpiece_young, workpiece_poisson)
    j_section = create_section(section_name='Jaw_section', material_name=j_material,type='solid')
    w_section = create_section(section_name='Workpiece_section', material_name=w_material,type='solid')#,thickness=0.002)
    
    mesh_part(part_name = jaw_part, size = jaw_mesh_size, dev_factor=0.1, min_size_factor = 0.1 )
    mesh_part(part_name = workpiece_part, size = workpiece_mesh_size, dev_factor=0.1, min_size_factor = 0.1 )

    assign_section(part_name = workpiece_part, section_name = w_section)
    assign_section(part_name = jaw_part, section_name = j_section)
    
    create_assembly(jawf, outer, length, tanf, radf, axlf, anglef)
    a = mdb.models['Model-1'].rootAssembly
    session.viewports['Viewport: 1'].setValues(displayedObject=a)
    session.viewports['Viewport: 1'].assemblyDisplay.setValues(optimizationTasks=OFF, geometricRestrictions=OFF, stopConditions=OFF)
    session.viewports['Viewport: 1'].odbDisplay.basicOptions.setValues(renderShellThickness=ON)
    session.writeOBJFile(fileName=home+"/tmp_model.obj",canvasObjects= (session.viewports['Viewport: 1'], ))
    # COMMENT_THIS_LINE
    run_job(home)
