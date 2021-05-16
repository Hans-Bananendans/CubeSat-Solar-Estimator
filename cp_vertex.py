# -*- coding: utf-8 -*-
"""
Created on Sat May 01 5:43:09 2021

@author: Johan Monster

Cube Projector - Vertex class
Version: 2.0

"""

import numpy as np
from numpy import sin, cos


class Vertex:
    """Custom implentation of a 3D point expressed in cartesians.
    
    Usage: Vertex has 3D carthesian coordinates, and some knowledge of its
    parent. The parent can be used to guide certain behaviour.
    
    Construction examples:
        p1 = Vertex()
            Creates a Vertex instance 'p1' at XYZ(0, 0, 0).
        
        p1 = Vertex([2, 4, 1])
            Creates a Vertex instance 'p1' at XYZ (2, 4, 1).
        
        p1 = Vertex([1, 1, 1], parenttype='frame')
            Creates a Vertex instance 'p1' at XYZ (1, 1, 1) and manually
            set the vertex parenttype to 'frame' (generally not necessary).
        
    """

    def __init__(self, xyz=[0, 0, 0], parenttype="global"):
        
        if not (isinstance(xyz,list) or isinstance(xyz, np.ndarray)):
            raise TypeError("Vertex constructor was not supplied with correct\
                            Vertex coordinates! Use a list or numpy.ndarray!")

        self.x = xyz[0]
        self.y = xyz[1]
        self.z = xyz[2]
        self.parenttype = "global"
        
        # Re-assign self.parenttype to check for validity of given parenttype
        if parenttype != "global":
            self.set_parenttype(parenttype)
        
        
    def set_parenttype(self, new_parenttype):
        """Sets own parenttype to a specified parenttype. First verifies 
        specified parenttype.
        """
        if new_parenttype in ["global","frame","geometry","face","vector"]:
            self.parenttype=new_parenttype
        else:
            raise ValueError("The parenttype cannot be anything other than: \
                             'global','frame','geometry','face','vector'!")
                             
    
    def translate(self, dx=0., dy=0., dz=0.):
        """Move the vertex' carthesian components by a relative amount.
        """
        self.x += dx
        self.y += dy
        self.z += dz
        
    
    def rotate(self, a, b, c, cor=None, seq='321'):
        """ Rotates the vertex around a point in space.
            a, b, c is Euler angle of rotation around local x, y, z axes
                expressed in radians
                
            seq is a string with the rotation sequence, e.g. '321' for:
                Rz(c).Ry(b).Rx(a).vertex
                Default rotation sequence: '321' (zyx)
                
            cor is the centre of rotation, which can be specified as
                a Vertex, or a list or numpy.ndarray of local coordinates.
        """
        
        # ==== COR argument handling ====
        # Case: COR is given as a list/ndarray of coordinates  
        if (isinstance(cor, list) or isinstance(cor, np.ndarray)):
            if len(cor) == 3:
                cor = Vertex(list(cor))
                
        # Case: COR is given as a Vertex
        elif isinstance(cor, Vertex):
            pass
            
        # All other cases
        else:
            raise TypeError("Error when rotating vertex. Centre of rotation \
                            argument not provided correctly!")
        
        # ==== Define relevant rotation matrices ====
        Rx = np.array([[1,      0,       0],
                       [0, cos(a), -sin(a)],
                       [0, sin(a),  cos(a)]])
        
        Ry = np.array([[ cos(b), 0, sin(b)],
                       [      0, 1,      0],
                       [-sin(b), 0, cos(b)]])
        
        Rz = np.array([[cos(c), -sin(c), 0],
                       [sin(c),  cos(c), 0],
                       [     0,       0, 1]])

        # ==== Perform rotation ====
        for axis in reversed(seq):
            
            # If no Centre of Rotation (cor) is given, take frame origin.
            # If not, subtract coordinates of cor first before applying 
            # the rotation, and then reverse this afterwards.
            if cor:
                self.translate(-1*cor.x, -1*cor.y, -1*cor.z)
            
            if axis == '1':
                (self.x, self.y, self.z) = \
                np.dot(Rx, np.array([self.x, self.y, self.z]))
                
            if axis == '2':
                (self.x, self.y, self.z) = \
                np.dot(Ry, np.array([self.x, self.y, self.z]))
                
            if axis == '3':
                (self.x, self.y, self.z) = \
                np.dot(Rz, np.array([self.x, self.y, self.z]))
            
            # Reverse temporary translation.
            if cor:
                self.translate(1*cor.x, 1*cor.y, 1*cor.z)
            
    
    def xyz(self):
        """Returns a numpy array with the local x, y, z coordinates.
        
            Output: np.array([x, y, z])
        """
        return np.array([self.x, self.y, self.z])
    
            
    def readout(self, dec=4):
        """More detailed information about Vertex.
           Use dec to define the number of decimals before rounding.
           Set dec to -1 to disable rounding.
           """
        if dec == -1:  # Set dec to -1 to disable rounding
            localcoords = self.xyz()
        else:
            localcoords = np.round(self.xyz(), dec)
        print("Local coordinates:   {}".format(localcoords))