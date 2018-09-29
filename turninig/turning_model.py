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
        region = self.part.Set(cells=c, name='Material region')
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