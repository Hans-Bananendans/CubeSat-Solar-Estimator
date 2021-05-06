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

    def __init__(self, xyz: list, parenttype="global"):
        
        self.x = xyz[0]
        self.y = xyz[1]
        self.z = xyz[2]
        self.parenttype=parenttype
        
        # print("[DEBUG] {} constructed.".format(self))
        
    # def __del__(self):
        # print("[DEBUG] Deleting {}.".format(self))
        
    def set_parenttype(self, new_parenttype):
        if new_parenttype in ["global", "frame", "geometry", "face"]:
            self.parenttype=new_parenttype
        else:
            raise ValueError("The parenttype cannot be anything other than: \
                             'global', 'frame', 'geometry', 'face'!")
    
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
        
        # COR argument handling:
        if type(cor) == list or type(cor) == np.array:
            if len(cor) == 3:
                cor = Vertex(list(cor))
        elif type(cor) == Vertex:
            pass
        else:
            raise TypeError("Error when rotating vertex. Centre of rotation \
                            argument not provided correctly!")
        
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
            
    # def __str__(self):
    #     """Print the vertex coordinates."""
    #     return str([round(self.x, 5),
    #                 round(self.y, 5),
    #                 round(self.z, 5)])
    
    def xyz(self):
        """Returns a numpy array with the local x, y, z coordinates."""
        return np.array([self.x, self.y, self.z])
    
    """MOVE TO FRAME!"""
    # def xyz_global(self):
    #     """Returns a numpy array with the global x, y, z coordinates."""
    #     return np.array(self.global_coordinates())
        
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