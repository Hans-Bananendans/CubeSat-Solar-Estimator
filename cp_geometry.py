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
        
        self.parenttype=parenttype # Default: "global"
        
        print("[DEBUG] {} constructed.".format(self))
    
    def __del__(self):
        print("[DEBUG] {} destructed.".format(self))
    
    def set_parenttype(self, new_parenttype):
        if new_parenttype in ["global", "frame"]:
            self.parent_type=new_parenttype
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
            if type(face) == Face:
                self.faces.append(face)
                face.set_parenttype("geometry")
            else:
                raise TypeError("Geometries can only be made from faces!")
    
    def vertices(self):
        """This method loops through all the Face objects in the geometry,
           lists the vertices making up each face, and removes the duplicates
           in the list. It then returns a list of unique vertices."""
        vertices = []
        for face in self.faces:
            face_vertices = face.vertices()
            for vertex in face_vertices:
                if vertex not in vertices:
                    vertices.append(vertex)
        return vertices
    
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
        return Vertex(xc, yc, zc, parent=self.parent)

    def translate(self, dx=0., dy=0., dz=0.):
        """Translate all the vertices in the geometry."""
        vertices = self.vertices()
        for vertex in vertices:
            vertex.translate(dx=dx, dy=dy, dz=dz)
    
    def unique_vertices(self):
        unique_vertices = []
        for face in self.faces:
            unique_vertices.extend(face.vertices())
        return list(set(unique_vertices))
            
    def rotate(self, a, b, c, cor, seq='321'):
        """Rotate all the vertices in the geometry around a point in space.
           a is Euler angle of rotation around x, etc...
               expressed in radians
           seq is a string with the rotation sequence, e.g. '321' for:
               Rz(c).Ry(b).Rx(a).vertex
           cor is the centre of rotation, which MUST be specified as
               a vertex.
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
    
    # def find_fcentroids(self, returnarray=False):
    #     """Returns a list of the centroid coordinates of all faces currently 
    #        in the geometry. Technically, these are vertex centroids.
    #        Optionally, centroids can be returned as a Numpy array by 
    #        setting returnarray to True."""
    #     fcentroids = []
    #     for i in range(len(self.faces)):
    #         if returnarray == True:
    #             fcentroids.append(self.faces[i].find_centroid().xyz())
    #         else:
    #             fcentroids.append(self.faces[i].find_centroid())
    #     return fcentroids

    # def find_perpendiculars(self, scale=1.):
    #     """Returns a list with a vector perpendicular to each face in the
    #        geometry. If scale=1, these vectors will be unit vectors.
            
    #        TODO: Guarantee that all perpendiculars point away from the 
    #        centre of the geometry, so that the output is not a mishmash of
    #        inward and outward pointing vertices. Can achieve this by taking
    #        inner product of:
    #            - vector pointing from geometry centroid to Face centroid.
    #            - computed perpendicular vector for given face.
    #        If the result is positive, the perpendicular is pointing away from 
    #        geometry, and if result is negative, flip the perpendicular.
    #        Downside: Does not work well for concave geometries whose centroid
    #        is outside the geometry boundaries.
    #        """
    #     perpendiculars = []
    #     for i in range(len(self.faces)):
    #         perp = self.faces[i].find_perpendicular()
    #         perpendiculars.append(perp * scale)

    #     # Code that ensures arrows point outwards (TODO):
    #     # ...

    #     return perpendiculars

    # def perpendiculars_plotlist(self):
    #     """Returns a list of lists with coordinates used for plotting.
    #        Each item in the list is a list of six coordinates:
    #        - The first three coordinates indicate the xyz of the tail point.
    #        - The second three coordinates are point a vector from this tail
    #          point.
           
    #        TODO: Pass scaling as argument.
    #        """
    #     c = self.find_fcentroids(returnarray=True)
    #     p = self.find_perpendiculars(scale=0.05)  # Vary 0.05 if needed.
    #     if len(c) != len(p):
    #         raise ValueError("Number of centroids and perpendiculars is "
    #                          "somehow unequal! This should not happen.")
    #     plotarray = np.hstack([np.vstack(c), np.vstack(p)]).transpose()
    #     return plotarray.tolist()
