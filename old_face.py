# -*- coding: utf-8 -*-
"""
Created on Thu April 27 15:02:11 2021    

@author: Johan Monster

Cube Projector - Face class
Version: 1.2

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
    """

    def __init__(self, p1: Vertex, p2: Vertex, p3: Vertex, p4: Vertex):
        self.p1 = p1
        self.p2 = p2
        self.p3 = p3
        self.p4 = p4
        if self.check_coplanarity() == False:
            raise ValueError("Error! Set of four points is not coplanar!")
        self.area = self.area()

    def check_coplanarity(self):
        """Check if all points are coplanar, by investigating determinant
               of [p12, p13, p14]. Points are coplanar only if it is zero."""
        plane_matrix = np.vstack((self.p2.array() - self.p1.array(),
                                  self.p3.array() - self.p1.array(),
                                  self.p4.array() - self.p1.array()
                                  )).transpose()
        if round(np.linalg.det(plane_matrix), 8) == 0:
            # Rounding here happens to prevent floating point errors from
            # interfering with the function. 8 decimals chosen arbitrarily.
            return True
        else:
            return False

    def area(self):
        """Calculate area of the face, by computing the diagonals."""
        diag13 = self.p3.array() - self.p1.array()
        diag24 = self.p4.array() - self.p2.array()
        cpdiag = np.cross(diag13, diag24)
        facearea = 0.5 * np.sqrt(np.einsum('...i,...i', cpdiag, cpdiag))
        return facearea
    
    def vertices(self):
        return [self.p1, self.p2, self.p3, self.p4]

    def project(self, plane='xy'):
        """Project a copy of the frame onto a plane that is spanned by two
           axes. The projection is orthographic.
           TODO: Generalize this function for arbitrary projection frame. """
        pj = deepcopy(self)
        for vertex in [pj.p1, pj.p2, pj.p3, pj.p4]:
            if plane == 'xy':
                vertex.z = 0
            elif plane == 'xz':
                vertex.y = 0
            elif plane == 'yz':
                vertex.x = 0
            else:
                raise ValueError("No valid projection plane given to "
                                 "projection method! "
                                 "Valid options: 'xy', 'xz', 'yz'")
        # Reinitialize the face with the projected vertices.
        pj = Face(pj.p1, pj.p2, pj.p3, pj.p4)
        return pj

    def find_centroid(self):
        """Find the vertex centroid of the face."""
        xc = 0.25 * (self.p1.x + self.p2.x + self.p3.x + self.p4.x)
        yc = 0.25 * (self.p1.y + self.p2.y + self.p3.y + self.p4.y)
        zc = 0.25 * (self.p1.z + self.p2.z + self.p3.z + self.p4.z)
        return Vertex(xc, yc, zc)

    def find_perpendicular(self):
        """Find a vector perpendicular to a face (direction is ambiguous).
           TODO: Generalize the direction of the perpendicular vector.
           Currently the direction depends on how the face is defined. """
        # Find perpendicular vector to plane spanned by p12, p14.
        # Then normalize it to a unit vector
        p12 = (self.p2.array() - self.p1.array())
        p14 = (self.p4.array() - self.p1.array())
        perpendicular = np.cross(p12, p14)
        return perpendicular / np.linalg.norm(perpendicular)

    def readout(self):
        """Print vertices of face to console."""
        for point in [self.p1, self.p2, self.p3, self.p4]:
            point.readout()

    def plotlist(self):
        """Return three lists with the x, y, and z-components of all four
           vertices in the face."""
        xlist = [self.p1.x, self.p2.x, self.p3.x, self.p4.x]
        ylist = [self.p1.y, self.p2.y, self.p3.y, self.p4.y]
        zlist = [self.p1.z, self.p2.z, self.p3.z, self.p4.z]
        return xlist, ylist, zlist

    def plotlist2(self):
        """Return a list of lists with the xyz coordinates of each vertex."""
        grid1 = np.array([self.p1.array(), self.p2.array(),
                          self.p3.array(), self.p4.array()])
        return grid1.tolist()
