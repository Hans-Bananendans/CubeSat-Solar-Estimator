# -*- coding: utf-8 -*-
"""
Created on Thu April 27 15:02:11 2021    

@author: Johan Monster

Cube Projector - Utility functions
Version: 1.2

"""

import numpy as np
from numpy import sin, cos, sqrt


def r2d(a):
    """Convert radians to degrees."""
    return a * 180 / np.pi


def d2r(a):
    """Convert degrees to radians."""
    return a * np.pi / 180


def e2q(a, b, c, seq='321'):
    """Converts Euler angles to quaternions."""
    q = np.array([
        [sin(.5*c)*cos(.5*b)*cos(.5*a) - cos(.5*c)*sin(.5*b)*sin(.5*a)],
        [cos(.5*c)*sin(.5*b)*cos(.5*a) + sin(.5*c)*cos(.5*b)*sin(.5*a)],
        [cos(.5*c)*cos(.5*b)*sin(.5*a) - sin(.5*c)*sin(.5*b)*cos(.5*a)],
        [cos(.5*c)*cos(.5*b)*cos(.5*a) + sin(.5*c)*sin(.5*b)*sin(.5*a)]
        ])
    return q

def q2e(quaternion, seq='321', returnarray=False):
    """Converts a unit quaternion to a set of euler angles."""
    q1, q2, q3, q4 = quaternion
    q1 = float(q1)
    q2 = float(q2)
    q3 = float(q3)
    q4 = float(q4)
    length = round(sqrt(q1*q1 + q2*q2 + q3*q3 + q4*q4), 8)
    if length != 1:
        raise ValueError("Quaternion given is not a unit quaternion. "
                         "|q|^2 = {} != 1".format(length))
    
    a = np.arctan2(2*(q4*q3 + q1*q2), 1-2*(q2*q2 + q3*q3))
    b = np.arcsin( 2*(q4*q2 - q3*q1))
    c = np.arctan2(2*(q4*q1 + q2*q3), 1-2*(q1*q1 + q2*q2))
    
    if returnarray == True:
        return np.array([round(a,10), round(b,10), round(c,10)])
    else:
        return round(a,10), round(b,10), round(c,10)
    
def euler_rotation(self, a, b, c, vector, seq='321'):
    """ Rotates a vector using sequential Euler rotations.
        a is Euler angle of rotation around x, etc...
            expressed in radians
        vector is the vector to be rotated
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
            vector = np.dot(Rx, vector)
        if axis == '2':
            vector = np.dot(Ry, vector)
        if axis == '3':
            vector = np.dot(Rz, vector)
    
    return vector