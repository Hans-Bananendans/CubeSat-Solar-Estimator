# -*- coding: utf-8 -*-
"""
Created on Sat May 01 5:43:09 2021 

@author: Johan Monster

Cube Projector - Vertex class
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
       and for best results they must be supplied sequentially and 
       counterclockwise, i.e.:
          4  o---o 3
             |   |  
          1  o---o 2
          
       TODO: Add vertex cleanup?
       TODO: Add projection
       """

    def __init__(self, p1: Vertex, p2: Vertex, p3: Vertex, p4: Vertex,
                 parenttype="global"):
        
        self.p1 = p1
        self.p2 = p2
        self.p3 = p3
        self.p4 = p4
        
        for vertex in self.vertices():
            vertex.set_parenttype("face")
            
        # self.special_vertices = [] # This is for special vertices only!
        
        self.parenttype=parenttype # Default: "global"
        
        if self.check_coplanarity() == False:
            raise ValueError("Error! Set of four points is not coplanar!")

        
        # print("[DEBUG] {} constructed.".format(self))
        
    # def __del__(self):
    #     """Custom deconstructor to clean up child-parent relationships."""
        
        # for vertex in [self.p1, self.p2, self.p3, self.p4]:
        #     del(vertex)
        # print("[DEBUG] Erasing c-p relation of {}".format(self))
    
        # print("[DEBUG] {} deconstructed..".format(self))

    def set_parenttype(self, new_parenttype):
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
        """Calculate area of the face, by computing the diagonals."""
        diag13 = self.p3.xyz() - self.p1.xyz()
        diag24 = self.p4.xyz() - self.p2.xyz()
        cpdiag = np.cross(diag13, diag24)
        facearea = 0.5 * np.sqrt(np.einsum('...i,...i', cpdiag, cpdiag))
        return facearea
    
    def vertices(self):
        """Returns a list of the four vertices spanning the face."""
        return [self.p1, self.p2, self.p3, self.p4]
    
    def corners(self):
        """Synonym function of self.vertices():"""
        return self.vertices()

    # def vertices_special(self):
    #     return self.special_vertices
    
    # def vertices_all(self):
    #     return self.vertices() + self.vertices_special()
    
    def readout(self, dec=4):
        """Print the vertex coordinates."""
        for vertex in [self.p1, self.p2, self.p3, self.p4]:
            print(vertex.xyz())
    
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

    """MOVE TO FRAME!"""
    # def plotlist_xyz_local(self, uselocal=False):
    #     """Return a list of lists with the xyz coordinates of each vertex."""
    #     # If object has no parent, return plot list in terms of global frame:
    #     if self.parent == None:
    #         uselocal = True
        
    #     if uselocal:
    #         # Fetch local vertex coordinates, and apply them to global frame
    #         grid1 = np.array([self.p1.xyz(), self.p2.xyz(),
    #                           self.p3.xyz(), self.p4.xyz()])
    #     else:
    #         # Fetch global vertex coordinates, and apply them to global frame
    #         grid1 = np.array([self.p1.xyz_global(), self.p2.xyz_global(),
    #                           self.p3.xyz_global(), self.p4.xyz_global()])
    #     return grid1.tolist()
    
    def project(self, new_frame=None, plane='xy'):
        """Project a copy of the frame onto a plane that is spanned by two
            axes. The projection is orthographic.
            
            TODO: A "new_frame" can be specified, which will make the 
            projection elements children of this frame. Otherwise, the 
            elements will be children of the global coordinate frame.
            
            TODO: Generalize this function for arbitrary projection frame. """
        
        # Eliminate minus sign in front of the plane:    
        if plane[0] == '-':
            plane = plane[1:]
        
        pj_xyz = []
        
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
            pj = Face(pj_xyz[0], pj_xyz[1], pj_xyz[2], pj_xyz[3],
                      parenttype="global")
        return pj
            
    def find_centroid(self):
        """Find the vertex centroid of the face."""
        xc = 0.25 * (self.p1.x + self.p2.x + self.p3.x + self.p4.x)
        yc = 0.25 * (self.p1.y + self.p2.y + self.p3.y + self.p4.y)
        zc = 0.25 * (self.p1.z + self.p2.z + self.p3.z + self.p4.z)
        return np.array([xc, yc, zc])

    def make_centroid(self):
        """Find the vertex centroid of the face, and turn it into a Vertex."""
        (xc, yc, zc) = self.find_centroid()
        return Vertex([xc, yc, zc])
    
    def find_perpendicular(self):
        """Find a vector perpendicular to a face (direction is ambiguous).
           TODO: Generalize the direction of the perpendicular vector.
           Currently the direction depends on how the face is defined. """
        # Find perpendicular vector to plane spanned by p12, p14.
        # Then normalize it to a unit vector
        p12 = (self.p2.xyz() - self.p1.xyz())
        p14 = (self.p4.xyz() - self.p1.xyz())
        perpendicular = np.cross(p12, p14)
        return perpendicular / np.linalg.norm(perpendicular)