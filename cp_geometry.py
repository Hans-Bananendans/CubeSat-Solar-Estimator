# -*- coding: utf-8 -*-
"""
Created on Sat May 01 5:43:09 2021

@author: Johan Monster

Cube Projector - Geometry class
Version: 2.0

"""

import numpy as np

from cp_vertex import Vertex
from cp_face import Face

class Geometry:
    """In this code, a Geometry is defined as a loose collection of Faces.
       The constructor merely creates an empty Geometry object, and can
       subsequently be filled with Faces using the add_face() and 
       add_faces methods.
       Some methods of this class interact with the Face objects, whilst
       others interact instead with the individual Vertices out of which
       the Faces are made, to avoid e.g. applying transformations more than
       once to the same Vertex object."""
       
    def __init__(self, faces: list = [], parenttype="global"):

        self.faces = []        
        self.add_faces(faces)
        
        self.parenttype= "global"
        
        # Re-assign self.parenttype to check for validity of given parenttype
        if parenttype != "global":
            self.set_parenttype(parenttype)
        
        
    def set_parenttype(self, new_parenttype):
        """Sets own parenttype to a specified parenttype. First verifies 
        specified parenttype.
        """
        if new_parenttype in ["global", "frame"]:
            self.parenttype=new_parenttype
        else:
            raise ValueError("The parenttype cannot be anything other than: \
                             'global', 'frame'!")


    def add_face(self, face: Face):
        """Add a singular face to the geometry."""
        self.faces.append(face)
        face.set_parenttype("geometry")


    def add_faces(self, faces: list):
        """Add a list of faces to the geometry."""
        for face in faces:
            if isinstance(face, Face):
                self.faces.append(face)
                face.set_parenttype("geometry")
            else:
                raise TypeError("This function can only add Face objects, \
                                but object of type {} was given!"
                                .format(type(face)))
    
    
    def vertices(self):
        vertices = []
        for face in self.faces:
            vertices.extend(face.vertices())
        return list(set(vertices))


    def translate(self, dx=0., dy=0., dz=0.):
        """Translate all the vertices in the geometry."""
        vertices = self.vertices()
        for vertex in vertices:
            vertex.translate(dx=dx, dy=dy, dz=dz)
       
            
    def rotate(self, a, b, c, cor, seq='321'):
        """Rotate all the vertices in the geometry around a point in space.
           a is Euler angle of rotation around x, etc...
               expressed in radians
           seq is a string with the rotation sequence, e.g. '321' for:
               Rz(c).Ry(b).Rx(a).vertex
           cor is the centre of rotation, which must be specified as
               a vertex, coordinate list, or coordinate numpy.ndarray
           """
        vertices = self.vertices()
        for vertex in vertices:
            vertex.rotate(a, b, c, cor=cor, seq=seq)
    
    
    def area(self):
        """Total area of all faces in the geometry.
        Note that faces have a single side, not two sides. """
        A = 0
        for face in self.faces:
            A += face.area
        return A
    
    
    def find_fcentroids(self):
        """Returns a list of the centroid coordinates of all faces currently 
            in the geometry. Technically, these are vertex centroids.
            Optionally, centroids can be returned as a Numpy array by 
            setting returnarray to True."""
        fcentroids = []
        for face in self.faces:
            fcentroids.append(list(face.find_centroid()))
        return fcentroids
    
    
    def find_cuboid_centroid(self):
        """Locate the centroid of a cuboid geometry.
           Cuboids have eight corners, which this method checks first.
           """
        vertices = self.vertices()
        if len(vertices) != 8:
            raise ValueError("Tried to use method find_cuboid_centroid() \
                             with a geometry containing {} vertices. \
                             All cuboids must have exactly 8 vertices. " \
                             .format(len(vertices)))
        xc = 0
        yc = 0
        zc = 0
        for vertex in vertices:
            xc += vertex.x/len(vertices)
            yc += vertex.y/len(vertices)
            zc += vertex.z/len(vertices)
        return np.array([xc, yc, zc])
    
    
    def make_cuboid_centroid(self):
        """Locates the centroid of a cuboid geometry and creates a
        Vertex at this point."""
        (xc, yc, zc) = self.find_cuboid_centroid()
        return Vertex([xc, yc, zc])

