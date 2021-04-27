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