# -*- coding: utf-8 -*-
"""
Created on Sat May 01 5:43:09 2021 

@author: Johan Monster

Cube Projector - Vertex class
Version: 2.0

"""

import numpy as np

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
                 parent=None):
        self.p1 = p1
        self.p2 = p2
        self.p3 = p3
        self.p4 = p4
        if self.check_coplanarity() == False:
            raise ValueError("Error! Set of four points is not coplanar!")
            
        self.area = self.calc_area()
        
        self.parent = parent
        # If parent is not None, also add vertex to the parent as a child.
        if parent is not None: 
            self.parent.add_face(self)

    def remove_parent(self):
        """If Vertex has a parent frame, removes it as a parent, and ensures 
           it is removed from the internal list in the frame.
           
           TODO: Check if vertices aren't used in other faces.
                 If not, remove that vertex as well.
           """
        if self.parent != None:
            self.parent.remove_face(self)
            # Fetch list of faces from parent.faces
            self.parent = None

        
    def change_parent(self, new_parent):
        """Connects vertex to another frame. If vertex was already attached 
           to a frame, it undoes this first.
           
           TODO: Add points to parent too.
           """
        # Remove old parent first (if applicable):
        if self.parent is not None:
            self.remove_parent()
        
        # Update parent in child
        self.parent = new_parent
        
        # Edit new parent to add new child (unless new parent is None):
        if new_parent is not None:
            self.parent.add_face(self)
    
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
    
    def calc_area(self):
        """Calculate area of the face, by computing the diagonals."""
        diag13 = self.p3.xyz() - self.p1.xyz()
        diag24 = self.p4.xyz() - self.p2.xyz()
        cpdiag = np.cross(diag13, diag24)
        facearea = 0.5 * np.sqrt(np.einsum('...i,...i', cpdiag, cpdiag))
        return facearea
    
    def vertices(self):
        return [self.p1, self.p2, self.p3, self.p4]

    def readout(self, dec=4):
        """Print the vertex coordinates."""
        for vertex in [self.p1, self.p2, self.p3, self.p4]:
            print(vertex)
            
    def plotlist(self, uselocal=False):
        """Return three lists with the x, y, and z-components of all four
           vertices in the face."""
        if self.parent == None:
            uselocal = True
        
        if uselocal:
            # Fetch local vertex coordinates, and apply them to global frame
            xlist = [self.p1.x, self.p2.x, self.p3.x, self.p4.x]
            ylist = [self.p1.y, self.p2.y, self.p3.y, self.p4.y]
            zlist = [self.p1.z, self.p2.z, self.p3.z, self.p4.z]
        else:
            # Fetch global vertex coordinates, and apply them to global frame
            xlist = [self.p1.global_coordinates()[0], 
                     self.p2.global_coordinates()[0],
                     self.p3.global_coordinates()[0],
                     self.p4.global_coordinates()[0]]
            ylist = [self.p1.global_coordinates()[1], 
                     self.p2.global_coordinates()[1],
                     self.p3.global_coordinates()[1],
                     self.p4.global_coordinates()[1]]
            zlist = [self.p1.global_coordinates()[2], 
                     self.p2.global_coordinates()[2],
                     self.p3.global_coordinates()[2],
                     self.p4.global_coordinates()[2]]        
        return xlist, ylist, zlist

    def plotlist2(self, uselocal=False):
        """Return a list of lists with the xyz coordinates of each vertex."""
        # If object has no parent, return plot list in terms of global frame:
        if self.parent == None:
            uselocal = True
        
        if uselocal:
            # Fetch local vertex coordinates, and apply them to global frame
            grid1 = np.array([self.p1.xyz(), self.p2.xyz(),
                              self.p3.xyz(), self.p4.xyz()])
        else:
            # Fetch global vertex coordinates, and apply them to global frame
            grid1 = np.array([self.p1.xyz_global(), self.p2.xyz_global(),
                              self.p3.xyz_global(), self.p4.xyz_global()])
        return grid1.tolist()
    
    # def project(self, plane='xy'):
    #     """Project a copy of the frame onto a plane that is spanned by two
    #        axes. The projection is orthographic.
    #        TODO: Generalize this function for arbitrary projection frame. """
    #     pj = deepcopy(self)
    #     for vertex in [pj.p1, pj.p2, pj.p3, pj.p4]:
    #         vertex.remove_parent()
    #         if plane == 'xy':
    #             vertex.z = 0
    #         elif plane == 'xz':
    #             vertex.y = 0
    #         elif plane == 'yz':
    #             vertex.x = 0
    #         else:
    #             raise ValueError("No valid projection plane given to "
    #                              "projection method! "
    #                              "Valid options: 'xy', 'xz', 'yz'")
    #     # Reinitialize the face with the projected vertices.
    #     pj = Face(pj.p1, pj.p2, pj.p3, pj.p4, pj.parent)
    #     pj.remove_parent()
    #     return pj
            
    def find_centroid(self):
        """Find the vertex centroid of the face."""
        xc = 0.25 * (self.p1.x + self.p2.x + self.p3.x + self.p4.x)
        yc = 0.25 * (self.p1.y + self.p2.y + self.p3.y + self.p4.y)
        zc = 0.25 * (self.p1.z + self.p2.z + self.p3.z + self.p4.z)
        return np.array([xc, yc, zc])

    def make_centroid(self):
        """Find the vertex centroid of the face, and turn it into a Vertex."""
        (xc, yc, zc) = self.find_centroid()
        return Vertex(xc, yc, zc, parent=self.parent)
    
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