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

    def __init__(self, xyz: list = [0., 0., 0.], geometry=[]):
        self.x = xyz[0]
        self.y = xyz[1]
        self.z = xyz[2]
        self.xdir = np.array([1,0,0])
        self.ydir = np.array([0,1,0])
        self.zdir = np.array([0,0,1])
        
        self.parenttype = "global"
        
        self.dcm = None
        self.recalculate_dcm()
        
        # self.vectors = []
        # self.vertices = []
        # self.faces = []
        self.geometries = []
        self.add_geometry(geometry)
        
        # print("[DEBUG] {} constructed.".format(self))
    
    # def __del__(self):
    #     print("[DEBUG] {} destructed.".format(self))


    # def add_vertex(self, vertex: Vertex):
    #     self.vertices.append(vertex)
    #     vertex.set_parenttype="frame"

    # def remove_vertex(self, vertex: Vertex):
    #     self.vertices.remove(vertex)
    #     vertex.set_parenttype="global"
        
    # def add_face(self, face: Face):
    #     """TODO: merge into one add_child method."""
    #     self.faces.append(face)
    #     face.set_parenttype="frame"

    # def remove_face(self, face: Face):
    #     """TODO: merge into one add_child method."""
    #     self.faces.remove(face)
    #     face.set_parenttype="global"

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

    def remove_geometry(self, geometry: Geometry):
        """TODO: merge into one add_child method."""
        self.geometries.remove(geometry)
        geometry.set_parenttype="global"
    
    # def add_vector(self, vector: Vector):
    #     """TODO: merge into one add_child method."""
    #     self.vectors.append(vector)
    #     vector.set_parenttype="frame"

    # def remove_vector(self, vector: Vector):
    #     """TODO: merge into one add_child method."""
    #     self.vectors.remove(vector)
    #     vector.set_parenttype="global"
    
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

    def rotate(self, a, b, c, cor=None, seq='321'):
        """TODO: Implement Centre of Rotation functionality."""
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
        
    def vertex_xyz_global(self, vertex):
        """Vertex coordinates in terms of parent frame."""
        
        # If vertex object was given:
        if type(vertex) == Vertex:
            xyz_local = vertex.xyz()
        elif type(vertex) == list or type(vertex) == np.ndarray:
            if len(vertex) == 3:
                xyz_local = np.array(vertex)
        else:
            raise TypeError("Function vertex_xyz_global does not accept \
                            arguments of this type!")
        
        # Coordinates of local frame origin in terms of global frame
        o_local = np.array([self.x, self.y, self.z])
        
        # Vertex coordinates in terms of global frame
        xyz_global = np.dot(self.dcm, xyz_local) + o_local
        
        # return xyz_global[0], xyz_global[1], xyz_global[2]
        return xyz_global
    
    def geometry_vertices(self):
        geometry_vertices = []
        for geometry in self.geometries:
            geometry_vertices.extend(geometry.vertices())
        return list(set(geometry_vertices))
    
    def plotlist(self, face: Face):
        """Return three lists with the x, y, and z-components of all four
           vertices in the face."""
        
        xlist = []
        ylist = []
        zlist = []
        
        for vertex in face.vertices():
            xg, yg, zg = self.vertex_xyz_global(vertex)
            xlist.append(xg)
            ylist.append(yg)
            zlist.append(zg)
        
        return xlist, ylist, zlist
    
    def plotlist_xyz(self, face: Face):
        """Return a list of lists with the xyz coordinates of each vertex."""
        # If object has no parent, return plot list in terms of global frame:
        
        plotlist_xyz = []
        
        for vertex in face.vertices():
            plotlist_xyz.append(list(self.vertex_xyz_global(vertex)))

        return plotlist_xyz
    
    def face_perpendicular(self, face: Face, scale=1.):
        
        p1 = self.vertex_xyz_global(face.p1)
        p2 = self.vertex_xyz_global(face.p2)
        p4 = self.vertex_xyz_global(face.p4)
        
        perpendicular = np.cross(p2-p1, p4-p1)
        return perpendicular / np.linalg.norm(perpendicular)
    
    def geometry_perpendiculars(self, geometry: Geometry, scale=1.):
        """Returns a list with a vector perpendicular to each face in the
            geometry. If scale=1, these vectors will be unit vectors.
            
            TODO: Guarantee that all perpendiculars point away from the 
            centre of the geometry, so that the output is not a mishmash of
            inward and outward pointing vertices. Can achieve this by taking
            inner product of:
                - vector pointing from geometry centroid to Face centroid.
                - computed perpendicular vector for given face.
            If the result is positive, the perpendicular is pointing away from 
            geometry, and if result is negative, flip the perpendicular.
            Downside: Does not work well for concave geometries whose centroid
            is outside the geometry boundaries.
            """
        perpendiculars = []
        for face in geometry.faces:
            perp = self.face_perpendicular(face)
            perpendiculars.append(perp * scale)

        # Code that ensures arrows point outwards (TODO):
        # ...

        return perpendiculars
    
    def face_centroid(self, face: Face):
        """Find the vertex centroid of the face."""
        
        xc = 0
        yc = 0
        zc = 0
        
        for vertex in face.vertices():
            (dx, dy, dz) = self.vertex_xyz_global(vertex)
            xc += 0.25*dx
            yc += 0.25*dy
            zc += 0.25*dz
    
        return np.array([xc, yc, zc])
    
    def geometry_face_centroids(self, geometry: Geometry):
        """Find the vertex centroid of the face."""
        face_centroids = []
        for face in geometry.faces:
            face_centroids.append(self.face_centroid(face))
        return np.array(face_centroids)
    # def make_centroid(self):
    #     """Find the vertex centroid of the face, and turn it into a Vertex."""
    #     (xc, yc, zc) = self.find_centroid()
    #     return Vertex(xc, yc, zc, parenttype="face")
    
    def perpendiculars_plotlist(self, geometry: Geometry, scale=.1):
        """Returns a list of lists with coordinates used for plotting.
            Each item in the list is a list of six coordinates:
            - The first three coordinates indicate the xyz of the tail point.
            - The second three coordinates are point a vector from this tail
              point.
              
            TODO: Pass scaling as argument.
            """
        c = self.geometry_face_centroids(geometry)
        p = self.geometry_perpendiculars(geometry, scale)  # Note: 0.05

        if len(c) != len(p):
            raise ValueError("Number of centroids and perpendiculars is "
                              "somehow unequal! This should not happen.")

        # Generate plotarray, which is structured as [xyzuvw1, xyzuvw2, ...]
        plotarray = np.hstack([np.vstack(c), np.vstack(p)]).transpose()
        return plotarray.tolist()
    
    def readout(self, dec=4):
        """More detailed information about Vertex.
           Use dec to define the number of decimals before rounding.
           Set dec to -1 to disable rounding.
           """
        if dec == -1:  # Set dec to -1 to disable rounding
            origin = self.origin()
        else:
            origin = np.round(self.origin(), dec)
        
        # children = len(self.vertices)   + \
        #            len(self.faces)      + \
        #            len(self.geometries) + \
        #            len(self.vectors)
        
        print("Carthesian coordinate frame.")
        print("Frame origin: {}\n".format(origin))
        # print("Children:              {}".format(children))
        # print("Associated Vertices:   {}".format(len(self.vertices)))
        # print("Associated Faces:      {}".format(len(self.faces)))
        print("Associated Geometries: {}".format(len(self.geometries)))
        # print("Associated Vectors:    {}".format(len(self.vectors)))