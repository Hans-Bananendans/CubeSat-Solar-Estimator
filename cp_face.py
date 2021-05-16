# -*- coding: utf-8 -*-
"""
Created on Sat May 01 5:43:09 2021 

@author: Johan Monster

Cube Projector - Face class
Version: 2.0

"""

import numpy as np
from copy import deepcopy
from cp_vertex import Vertex

class Face:
    """Class for a quadrilateral face (4-point polygon).
    
       A face consists of four vertices that make up the corners of a plane.
       All four vertices must be coplanar, and the class constructor will
       verify this. The order in which the vertices are supplied DOES matter, 
       as the faces have one side only (as opposed to a front and a backside)
       
       For best results the corners must be supplied sequentially and 
       counterclockwise, i.e.:
           
          4  o---o 3
             |   |  
          1  o---o 2
          
       This will create a face that has a front area, and the normal of the 
       face will point upward (towards the reader). If supplied in clockwise
       direction instead, this face would have identical points, but its
       normal would point down (away from the reader, into the screen).
          
       Construction example:
              
           f1 = Face(Vertex([ 1, 0, 0]), 
                     Vertex([ 1, 1, 0]),
                     Vertex([-1, 1, 0]), 
                     Vertex([-1,-1, 0])
                     )
               Creates an instance of the Face class named 'f1', which forms
               a square face in the XY-plane, with sides of length 2, centred
               on the local origin.             
       """

    def __init__(self, p1: Vertex, p2: Vertex, p3: Vertex, p4: Vertex,
                 parenttype="global"):
        
        self.p1 = p1
        self.p2 = p2
        self.p3 = p3
        self.p4 = p4
                
        # Ensure provided vertices are coplanar
        if self.check_coplanarity() == False:
            raise ValueError("Error! Set of four points is not coplanar!")
            
        self.parenttype = "global"
        
        # Re-assign self.parenttype to check for validity of given parenttype
        if parenttype != "global":
            self.set_parenttype(parenttype)
        
        # Let child vertices know their parenttype is 'face'
        for vertex in self.vertices():
            vertex.set_parenttype("face")


    def set_parenttype(self, new_parenttype):
        """Sets own parenttype to a specified parenttype. First verifies 
        specified parenttype.
        """
        if new_parenttype in ["global", "frame", "geometry"]:
            self.parenttype=new_parenttype
            
        else:
            raise ValueError("The parenttype cannot be anything other than: \
                             'global', 'frame', 'geometry'!")
                             
    def check_coplanarity(self):
        """Check if all points are coplanar, by investigating determinant
               of [p12, p13, p14]. Points are coplanar only if it is zero."""
        plane_matrix = np.vstack((self.p2.xyz() - self.p1.xyz(),
                                  self.p3.xyz() - self.p1.xyz(),
                                  self.p4.xyz() - self.p1.xyz()))
        if round(np.linalg.det(plane_matrix).transpose(), 8) == 0:
            # Rounding here happens to prevent floating point errors from
            # interfering with the function. 8 decimals chosen arbitrarily.
            return True
        else:
            return False
    
    def area(self):
        """Calculate area of the face, by computing the diagonals.
        
           Returns: face area (float)"""
        diag13 = self.p3.xyz() - self.p1.xyz()
        diag24 = self.p4.xyz() - self.p2.xyz()
        cpdiag = np.cross(diag13, diag24)
        facearea = 0.5 * np.sqrt(np.einsum('...i,...i', cpdiag, cpdiag))
        return facearea
    
    def vertices(self):
        """Returns a list of the four vertices spanning the face."""
        return [self.p1, self.p2, self.p3, self.p4]
    
    def corners(self):
        """Synonym function of self.vertices()."""
        return self.vertices()
    
    def plotlist_local(self):
        """Return three lists with the x, y, and z-components of all four
            vertices in the face."""
        
        xlist = []
        ylist = []
        zlist = []
        
        for vertex in self.vertices():
            xlist.append(vertex.x)
            ylist.append(vertex.y)
            zlist.append(vertex.z)
  
        return xlist, ylist, zlist
    
    def project(self, plane='xy', new_frame=None):
        """Project a copy of the face onto a plane that is spanned by two
            axes. The projection is orthographic.
            
            If this face is not attached to the global coordinate frame, then
            the face returned by this function will be in terms of the local
            frame.
            
            TODO: A "new_frame" can be specified, which will make the 
            projection elements children of this frame. Otherwise, the 
            elements will be children of the global coordinate frame.
            
            TODO: Generalize this function for arbitrary projection plane. 
            """
        
        # Eliminate minus sign in front of the plane:    
        if plane[0] == '-':
            plane = plane[1:]
        
        pj_xyz = []
        
        # Flatten the coordinates of each vertex onto the projection plane:
        for vertex in [self.p1, self.p2, self.p3, self.p4]:
            if plane == 'xy':
                pj_x = vertex.x
                pj_y = vertex.y
                pj_z = 0
            elif plane == 'xz':
                pj_x = vertex.x
                pj_y = 0
                pj_z = vertex.z
            elif plane == 'yz':
                pj_x = 0
                pj_y = vertex.y
                pj_z = vertex.z
            else:
                raise ValueError("No valid projection plane given to "
                                  "projection method! "
                                  "Valid options: 'xy', 'xz', 'yz'")
            pj_xyz.append(Vertex([pj_x, pj_y, pj_z], parenttype="face"))
        
        # Create the projected face using the projected vertices.
        if not new_frame: 
            pj = Face(pj_xyz[0], pj_xyz[1], pj_xyz[2], pj_xyz[3])
        else:
            pj = Face(pj_xyz[0], pj_xyz[1], pj_xyz[2], pj_xyz[3],
                      parenttype="frame")
            new_frame.add_face(pj)
        return pj
            
    def find_centroid(self):
        """Find the vertex centroid of the face, and then return its
           local coordinates in an numpy.ndarray.
           
           Returns: numpy.ndarray([xc, yc, zc])
           """
        xc = 0.25 * (self.p1.x + self.p2.x + self.p3.x + self.p4.x)
        yc = 0.25 * (self.p1.y + self.p2.y + self.p3.y + self.p4.y)
        zc = 0.25 * (self.p1.z + self.p2.z + self.p3.z + self.p4.z)
        return np.array([xc, yc, zc])

    def make_centroid(self):
        """Find the vertex centroid of the face, and turn it into a Vertex.
        
           Returns: Vertex([xc, yc, zc])
           """
        (xc, yc, zc) = self.find_centroid()
        return Vertex([xc, yc, zc])
    
    def find_perpendicular(self):
        """Find a vector perpendicular to a face (direction is ambiguous).
           TODO: Generalize the direction of the perpendicular vector.
           Currently the direction depends on how the face is defined. """
           
        # Find perpendicular vector to plane spanned by p12, p14.
        p12 = (self.p2.xyz() - self.p1.xyz())
        p14 = (self.p4.xyz() - self.p1.xyz())
        perpendicular = np.cross(p12, p14)
        
        # Then normalize it to a unit vector and return
        return perpendicular / np.linalg.norm(perpendicular)
    
    def readout(self, dec=4):
        """Print the vertex coordinates of all corners."""
        for vertex in [self.p1, self.p2, self.p3, self.p4]:
            print(vertex.xyz())