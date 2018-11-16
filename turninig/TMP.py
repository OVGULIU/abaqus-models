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

    def mesh(self, size, dev_factor, min_size_factor):
        self.part.setMeshControls(regions=self.part.cells, technique=SWEEP, algorithm=MEDIAL_AXIS)
        self.part.seedPart(size=size, deviationFactor=dev_factor, minSizeFactor=min_size_factor)
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

    def mesh(self, size, dev_factor, min_size_factor):
        self.part.setMeshControls(regions=self.part.cells, technique=SWEEP, algorithm=MEDIAL_AXIS)
        self.part.seedPart(size=size, deviationFactor=dev_factor, minSizeFactor=min_size_factor)
        self.part.generateMesh()



class Assembly:

    def __init__(self, workpiece, jaw):
        self.workpiece = workpiece
        self.jaw = jaw

    def __create_workpiece_partion(self):
        p = workpiece.part
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
        

    def create_assembly(self, jawf, outer, length):
        self.interaction_property = InteractionProperty()
        a = mdb.models['Model-1'].rootAssembly
        a.DatumCsysByDefault(CARTESIAN)
        p = mdb.models['Model-1'].parts['Workpiece']
        a.Instance(name='Part-1', part=p, dependent=ON)
        p = mdb.models['Model-1'].parts['Jaw']
        a.Instance(name='Jaw-1', part=p, dependent=ON)
        a.rotate(instanceList=('Jaw-1', ), axisPoint=(0, 0, 0), axisDirection=OX, angle=rad(-90.0))
        a.translate(instanceList=('Jaw-1', ), vector=(0, outer, 0.5 * self.jaw.height))
        a.RadialInstancePattern(instanceList=('Jaw-1', ), point=(0, 0, 0), axis=OZ, number=3, totalAngle=-360)
        mdb.models['Model-1'].rootAssembly.features.changeKey(fromName='Jaw-1', toName='Jaw-1-rad-1')
        sys.__stdout__.write("Create assembly"+"\n") 
        self.create_interaction()
        c_systems = self.create_CSYS()
        self.create_step()
        self.create_jaw_BSs(c_systems)
        self.apply_jaw_force(jawf, c_systems)
        self.__create_workpiece_partion()
        mdb.models['Model-1'].parts['Workpiece'].generateMesh()
        a.regenerate()
        

    def apply_jaw_force(self,value, c_systems):
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
        

    def create_interaction(self):
        z = self.jaw.length/100
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
            p1 = rotate((0.25 * self.jaw.length, outer, 0.25 * self.jaw.height), OZ, deg(jaw_angle))
            p2 = rotate((-0.25 * self.jaw.length, outer, 0.25 * self.jaw.height), OZ, deg(jaw_angle))
            p3 = rotate((0.25 * self.jaw.length, outer, 0.75 * self.jaw.height), OZ, deg(jaw_angle))
            p4 = rotate((-0.25 * self.jaw.length, outer, 0.75 * self.jaw.height), OZ, deg(jaw_angle))

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
                thickness=ON, interactionProperty=self.interaction_property.name, adjustMethod=NONE, 
                initialClearance=OMIT, datumAxis=None, clearanceRegion=None)

    def create_CSYS(self):
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

    def create_jaw_BSs(self,c_systems):
        for j in JAWS:
            jaw = 'Jaw-1-rad-'+str(j)
            #starts from 0 cause c_systems is passed to function
            # JAWS=1,2,3; csys_n = 0,1,2
            csys_n = (-j+2) %3
            a = mdb.models['Model-1'].rootAssembly
            f1 = a.instances['Jaw-1-rad-'+str(j)].faces
            #faces1 = f1.getSequenceFromMask(mask=('[#402 ]', ), )
            x1, y1 = pol2cart(outer+self.jaw.width/2, 360/len(JAWS)*csys_n+0.5)
            x2, y2 = pol2cart(outer+self.jaw.width/2, 360/len(JAWS)*csys_n-0.5)
            
            # faces1 = f1.findAt((x1,y1, 0),(x2,y2, 0))
            # sys.__stdout__.write( str(faces1))
            faces1 = f1.getSequenceFromMask(mask=('[#402 ]', ), )
            # faces1 = f1.getSequenceFromMask(mask=('[#8200 ]', ), )
            region = a.Set(faces = faces1, name='Jaw_BS_set-'+str(j))
            datum = mdb.models['Model-1'].rootAssembly.datums[c_systems[csys_n].id]
            mdb.models['Model-1'].DisplacementBC(name='BC-'+str(j), createStepName='Initial', 
                region=region, u1=UNSET, u2=SET, u3=SET, ur1=UNSET, ur2=UNSET, ur3=UNSET, 
                amplitude=UNSET, distributionType=UNIFORM, fieldName='', localCsys=datum)

                
    def create_step(self):
        mdb.models['Model-1'].StaticStep(name='Step-1', previous='Initial')
        

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
fparameter = 0.6 # axial position of force (?)
jcount = 3
JAWS = range(1,jcount+1)

jawf = 1000
f_coeff = 0.3

workpiece_mesh_size = 0.002


if __name__ == "__main__":
    jaw_num = 3
    steel = Material('Steel', 210e15, 0.29)
    steel_section = model.HomogeneousSolidSection(name='Steel-section', material=steel.name, thickness=None)

    aluminum = Material('Aluminum', 0.7e9, 0.28)
    aluminum_section = model.HomogeneousSolidSection(name='Aluminum-section', material=aluminum.name, thickness=None)
    
    workpiece = Workpiece(length=mm(60), inner=mm(59/2), outer=mm(68/2), p_num=jaw_num)
    workpiece.set_section(aluminum_section)
    workpiece.mesh(size=0.002, dev_factor=0.1, min_size_factor = 0.1 )


    jaw = Jaw(length=mm(15), width=mm(15), height=mm(15))
    jaw.set_section(steel_section)
    jaw.partition()
    jaw.mesh(size = 0.0015, dev_factor=0.1, min_size_factor = 0.1 )

    assembly = Assembly(workpiece=workpiece, jaw=jaw)
    assembly.create_assembly(jawf, outer, length)
    a = mdb.models['Model-1'].rootAssembly
    session.viewports['Viewport: 1'].setValues(displayedObject=a)
    session.viewports['Viewport: 1'].view.setValues(session.views['Iso'])

    run_job(home)




