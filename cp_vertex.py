# -*- coding: utf-8 -*-
"""
Created on Sat May 01 5:43:09 2021

@author: Johan Monster

Cube Projector - Vertex class
Version: 2.0

"""

import numpy as np
# from copy import deepcopy
from numpy import sin, cos
# import matplotlib.pyplot as plt
# import matplotlib.animation as animation
# import mpl_toolkits.mplot3d as mp3d

class Vertex:
    """Custom implentation of a 3D point expressed in cartesians.
    
    TODO: Implement ability to project points
    """

    def __init__(self, x=0., y=0., z=0., parent=None):
        self.x = x
        self.y = y
        self.z = z
        self.connect_parent(parent)
    
    def remove_parent(self):
        """If Vertex has a parent frame, removes it as a parent, and ensures 
           it is removed from the internal list in the frame.
           """
        if self.parent != None:
            self.parent.remove_vertex(self)
            self.parent = None

    def connect_parent(self, new_parent):
        """Connects vertex to a frame, unless:
            - the given new_parent is None
            - the vertex is already in the frame (duplicate control)
           """
        # Set own parent to new_parent:
        self.parent = new_parent
        
        # Ensure new_parent is not None:
        if new_parent != None:
            # Ensure vertex is not already connected to parent
            if self not in self.parent.vertices:
                self.parent.add_vertex(self)
        
    def change_parent(self, new_parent):
        """Connects vertex to another frame. If vertex was already attached 
           to a frame, it undoes this first.
           """
        # Remove old parent first (if applicable):
        if self.parent != None:
            self.remove_parent()
        
        # Connect new parent:
        self.parent = new_parent
        self.parent.add_vertex(self)
    
    def global_coordinates(self):
        # Vertex coordinates in terms of parent frame
        xyz_local = np.array([self.x, self.y, self.z])
        
        # Coordinates of local frame origin in terms of global frame
        o_local = np.array([self.parent.x, self.parent.y, self.parent.z])
        
        # Vertex coordinates in terms of global frame
        xyz_global = np.dot(self.parent.dcm, xyz_local) + o_local
        
        return xyz_global[0], xyz_global[1], xyz_global[2]
    
    def translate(self, dx=0., dy=0., dz=0.):
        self.x += dx
        self.y += dy
        self.z += dz
    
    def rotate(self, a, b, c, cor=None, seq='321'):
        """ Rotates the vertex around a point in space.
            a is Euler angle of rotation around x, etc...
                expressed in radians
            seq is a string with the rotation sequence, e.g. '321' for:
                Rz(c).Ry(b).Rx(a).vertex
            cor is the centre of rotation, which MUST be specified as
                a vertex.
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
            
            # If no Centre of Rotation (cor) is given, take frame origin.
            # If not, subtract coordinates of cor first before applying 
            # the rotation, and then reverse this afterwards.
            if cor != None:
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
            if cor != None:
                self.translate(1*cor.x, 1*cor.y, 1*cor.z)
            
    def __str__(self):
        """Print the vertex coordinates."""
        return str([round(self.x, 5),
                    round(self.y, 5),
                    round(self.z, 5)])
    
    def xyz(self):
        """Returns a numpy array with the local x, y, z coordinates."""
        return np.array([self.x, self.y, self.z])
    
    def xyz_global(self):
        """Returns a numpy array with the global x, y, z coordinates."""
        return np.array(self.global_coordinates())
        
    def readout(self, dec=4):
        """More detailed information about Vertex.
           Use dec to define the number of decimals before rounding.
           Set dec to -1 to disable rounding.
           """
        if dec == -1:  # Set dec to -1 to disable rounding
            localcoords = self.xyz()
            globalcoords = self.xyz_global()
            parentorigin = self.parent.origin()
        else:
            localcoords = np.round(self.xyz(), dec)
            globalcoords = np.round(self.xyz_global(), dec)
            parentorigin = np.round(self.parent.origin(), dec)
        
        
        print("Parent frame: {}".format(self.parent))
        print("Parent frame origin: {}".format(parentorigin))
        print("Local coordinates:   {}".format(localcoords))
        print("Global coordinates:  {}".format(globalcoords))