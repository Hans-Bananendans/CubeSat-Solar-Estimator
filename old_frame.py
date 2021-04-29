# -*- coding: utf-8 -*-
"""
Created on Thu April 27 15:02:11 2021    

@author: Johan Monster

Cube Projector - Vertex class
Version: 1.2

"""

import numpy as np
from numpy import sin, cos

class Frame:
    """Custom implentation of a 3D point expressed in cartesians."""

    def __init__(self, x=0., y=0., z=0.):
        self.origin = np.array([0,0,0])
        self.x = x
        self.y = y
        self.z = z

    def translate(self, dx=0., dy=0., dz=0.):
        self.x += dx
        self.y += dy
        self.z += dz

    def rotate_origin(self, a, b, c, seq='321'):
        """ Rotates the vertex around the origin.
            a is Euler angle of rotation around x, etc...
                expressed in radians
            seq is a string with the rotation sequence, e.g. '321' for:
                Rz(c).Ry(b).Rx(a).vertex
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
            if axis == '1':
                newcoor = np.dot(Rx, np.array([self.x, self.y, self.z]))
            if axis == '2':
                newcoor = np.dot(Ry, np.array([self.x, self.y, self.z]))
            if axis == '3':
                newcoor = np.dot(Rz, np.array([self.x, self.y, self.z]))
            self.__init__(newcoor[0], newcoor[1], newcoor[2])

    def rotate(self, a, b, c, cor, seq='321'):
        """ Rotates the vertex around a point in space.
            a is Euler angle of rotation around x, etc...
                expressed in radians
            seq is a string with the rotation sequence, e.g. '321' for:
                Rz(c).Ry(b).Rx(a).vertex
            cor is the centre of rotation, which MUST be specified as
                a vertex.
        """

        Rx = np.array([[1, 0, 0],
                       [0, cos(a), -sin(a)],
                       [0, sin(a), cos(a)]])

        Ry = np.array([[cos(b), 0, sin(b)],
                       [0, 1, 0],
                       [-sin(b), 0, cos(b)]])

        Rz = np.array([[cos(c), -sin(c), 0],
                       [sin(c), cos(c), 0],
                       [0, 0, 1]])

        for axis in reversed(seq):
            self.translate(-1 * cor.x, -1 * cor.y, -1 * cor.z)

            if axis == '1':
                newcoor = np.dot(Rx, np.array([self.x, self.y, self.z]))
            if axis == '2':
                newcoor = np.dot(Ry, np.array([self.x, self.y, self.z]))
            if axis == '3':
                newcoor = np.dot(Rz, np.array([self.x, self.y, self.z]))
            self.__init__(newcoor[0], newcoor[1], newcoor[2])

            self.translate(1 * cor.x, 1 * cor.y, 1 * cor.z)

    def array(self):
        """Returns vertex coordinates as a Numpy array."""
        return np.array([self.x, self.y, self.z])

    def readout(self, dec=4):
        """Print the vertex coordinates."""
        if dec == -1:  # Set dec to -1 to disable rounding
            print([self.x, self.y, self.z])
        else:
            print([round(self.x, dec),
                   round(self.y, dec),
                   round(self.z, dec)])