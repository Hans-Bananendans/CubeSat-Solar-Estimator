# -*- coding: utf-8 -*-
"""
Created on Sat May 01 5:43:09 2021   

@author: Johan Monster

Cube Projector - Frame class
Version: 2.0

"""

import numpy as np
from numpy import sin, cos

from cp_vertex import Vertex
from cp_face import Face
from cp_geometry import Geometry
from cp_vector import Vector

class Frame:
    """Custom implentation of a 3D point expressed in cartesians."""

    def __init__(self, x=0., y=0., z=0.):
        self.x = x
        self.y = y
        self.z = z
        self.xdir = np.array([1,0,0])
        self.ydir = np.array([0,1,0])
        self.zdir = np.array([0,0,1])
        
        self.dcm = None
        self.recalculate_dcm()
        
        self.vertices = []
        self.faces = []
        self.geometries = []
        self.vectors = []
    
    def __del__(self):
        print("[DEBUG] {} destructed.".format(self))

    def add_vertex(self, vertex: Vertex):
        self.vertices.append(vertex)

    def remove_vertex(self, vertex: Vertex):
        self.vertices.remove(vertex)
        
    def add_face(self, face: Face):
        """TODO: merge into one add_child method."""
        self.faces.append(face)

    def remove_face(self, face: Face):
        """TODO: merge into one add_child method."""
        self.faces.remove(face)

    def add_geometry(self, geometry: Geometry):
        """TODO: merge into one add_child method."""
        self.geometries.append(geometry)

    def remove_geometry(self, geometry: Geometry):
        """TODO: merge into one add_child method."""
        self.geometries.remove(geometry)
    
    def add_vector(self, vector: Vector):
        """TODO: merge into one add_child method."""
        self.vectors.append(vector)

    def remove_vector(self, vector: Vector):
        """TODO: merge into one add_child method."""
        self.vectors.remove(vector)
    
    def origin(self):
        return np.array([self.x, self.y, self.z])
    
    def recalculate_dcm(self):
        gx = np.array([1,0,0])
        gy = np.array([0,1,0])
        gz = np.array([0,0,1])
        self.dcm = np.array([
            [np.dot(gx, self.xdir), 
             np.dot(gx, self.ydir), 
             np.dot(gx, self.zdir)],
            
            [np.dot(gy, self.xdir), 
             np.dot(gy, self.ydir), 
             np.dot(gy, self.zdir)],
            
            [np.dot(gz, self.xdir), 
             np.dot(gz, self.ydir), 
             np.dot(gz, self.zdir)]])
        
    def translate(self, dx=0., dy=0., dz=0.):
        self.x += dx 
        self.y += dy
        self.z += dz

    def rotate(self, a, b, c, seq='321'):
        """ 
        """
        Rx = np.array([[1,      0,       0],
                        [0, cos(a), -sin(a)],
                        [0, sin(a),  cos(a)]])
        
        Ry = np.array([[ cos(b), 0, sin(b)],
                        [      0, 1,      0],
                        [-sin(b), 0, cos(b)]])
        
        Rz = np.array([[cos(c), -sin(c), 0],
                        [sin(c),  cos(c), 0],
                        [     0,       0, 1]])
        
        for axis in reversed(seq):
            temp = (self.x, self.y, self.z)
            self.translate(-1*temp[0], -1*temp[1], -1*temp[2])
            
            if axis == '1':
                self.xdir = np.dot(Rx, self.xdir) # Redundant?
                self.ydir = np.dot(Rx, self.ydir)
                self.zdir = np.dot(Rx, self.zdir)
            if axis == '2':
                self.xdir = np.dot(Ry, self.xdir)
                self.ydir = np.dot(Ry, self.ydir) # Redundant?
                self.zdir = np.dot(Ry, self.zdir)
            if axis == '3':
                self.xdir = np.dot(Rz, self.xdir)
                self.ydir = np.dot(Rz, self.ydir)
                self.zdir = np.dot(Rz, self.zdir) # Redundant?
                
            self.translate(1*temp[0], 1*temp[1], 1*temp[2])
        
        self.recalculate_dcm()
    
    def readout(self, dec=4):
        """More detailed information about Vertex.
           Use dec to define the number of decimals before rounding.
           Set dec to -1 to disable rounding.
           """
        if dec == -1:  # Set dec to -1 to disable rounding
            origin = self.origin()
        else:
            origin = np.round(self.origin(), dec)
        
        children = len(self.vertices)   + \
                   len(self.faces)      + \
                   len(self.geometries) + \
                   len(self.vectors)
        
        print("Carthesian coordinate frame.")
        print("Frame origin: {}\n".format(origin))
        print("Children: {}".format(children))
        print("Associated Vertices:   {}".format(len(self.vertices)))
        print("Associated Faces:      {}".format(len(self.faces)))
        print("Associated Geometries: {}".format(len(self.geometries)))
        print("Associated Vectors:    {}".format(len(self.vectors)))