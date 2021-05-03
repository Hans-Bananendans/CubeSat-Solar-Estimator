# -*- coding: utf-8 -*-
"""
Created on Mon May  3 08:57:58 2021

@author: Johan Monster

Frame experiment 3
"""

import numpy as np
from copy import deepcopy
from numpy import sin, cos
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import mpl_toolkits.mplot3d as mp3d

# from cp_vertex import Vertex
# from cp_face import Face
# from cp_geometry import Geometry
# from cp_vector import Vector
# from cp_frame import Frame
from cp_utilities import d2r#, r2d
from cp_plotting import plot_global_tripod, plot_frame


class Frame:
    def __init__(self, x=0., y=0., z=0.):
        self.x = x
        self.y = y
        self.z = z
        self.xdir = np.array([1,0,0])
        self.ydir = np.array([0,1,0])
        self.zdir = np.array([0,0,1])
        
        self.dcm = None
        self.geometries = []
        
        print("[DEBUG] {} constructed.".format(self))
    
    def __del__(self):
        print("[DEBUG] {} destructed.".format(self))
        
    def add_geometry(self, geometry):
        if type(geometry) == Geometry:
            self.geometries.append(geometry)
        elif type(geometry) == list:
            for geo in geometry:
                if type(geo) == Geometry:
                    self.geometries.append(geometry)
                else:
                    raise TypeError("Frame.add_geometry() argument must be a \
                                    list of geometries!")
        else:
            raise TypeError("Invalid argument provided to \
                            Frame.add_geometry()!")

class Vertex:
    def __init__(self, x=0., y=0., z=0.):
        self.x = x
        self.y = y
        self.z = z
        
        print("[DEBUG] {} constructed.".format(self))
    
    def __del__(self):
        print("[DEBUG] {} destructed.".format(self))

class Face:
    def __init__(self, p1: Vertex, p2: Vertex, p3: Vertex, p4: Vertex):
        self.p1 = p1
        self.p2 = p2
        self.p3 = p3
        self.p4 = p4
        
        print("[DEBUG] {} constructed.".format(self))
    
    def __del__(self):
        print("[DEBUG] {} destructed.".format(self))
    
class Geometry:
    def __init__(self, faces: list = []):
        
        self.faces = faces
            
        print("[DEBUG] {} constructed.".format(self))
    
    def __del__(self):
        print("[DEBUG] {} destructed.".format(self))
        
        
p1 = Vertex([-0.05, -0.05, -0.1])
p2 = Vertex([-0.05, -0.05,  0.1])
p3 = Vertex([ 0.05, -0.05,  0.1])
p4 = Vertex([ 0.05, -0.05, -0.1])
p5 = Vertex([-0.05,  0.05, -0.1])
p6 = Vertex([-0.05,  0.05,  0.1])
p7 = Vertex([ 0.05,  0.05,  0.1])
p8 = Vertex([ 0.05,  0.05, -0.1])

fA = Face(p4, p3, p2, p1)
fB = Face(p2, p3, p7, p6)
fC = Face(p3, p4, p8, p7)
fD = Face(p4, p1, p5, p8)
fE = Face(p1, p2, p6, p5)
fF = Face(p5, p6, p7, p8)

cubesat = Geometry([fA, fB, fC, fD, fE, fF])

frame1 = Frame(1,1,1)
frame1.add_geometry(cubesat)

# f1 = Face(Vertex(0,0,0), Vertex(1,0,0), Vertex(1,1,0), Vertex(0,1,0))
# f2 = Face(Vertex(0,0,0), Vertex(1,0,0), Vertex(1,1,0), Vertex(0,1,0))
# f3 = Face(Vertex(0,0,0), Vertex(1,0,0), Vertex(1,1,0), Vertex(0,1,0))
# f4 = Face(Vertex(0,0,0), Vertex(1,0,0), Vertex(1,1,0), Vertex(0,1,0))
# f5 = Face(Vertex(0,0,0), Vertex(1,0,0), Vertex(1,1,0), Vertex(0,1,0))
# f6 = Face(Vertex(0,0,0), Vertex(1,0,0), Vertex(1,1,0), Vertex(0,1,0))