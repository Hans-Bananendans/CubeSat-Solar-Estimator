# -*- coding: utf-8 -*-
"""
Created on Sat May 01 5:43:09 2021

@author: Johan Monster

Cube Projector - Utility functions
Version: 2.0

"""

import numpy as np
from numpy import sin, cos
    
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
    length = round(np.sqrt(q1*q1 + q2*q2 + q3*q3 + q4*q4), 8)
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
    
def hex2nRGB(hex_colour):
    """Converts 3, 6, or 9-length hexidecimal codes to normalized RGB values.
        The normalization is around 1, i.e. 00 -> 0.0 and FF -> 1.0 
        
        Returns: tuple(normalized R, normalized G, normalized B)
        """
    # Strip hash off colour code:
    hex_colour = hex_colour.lstrip('#')
    
    # Verifying formatting (check if all hexidecimal numbers):
    if not all(c.lower() in "0123456789abcdef" for c in hex_colour):
        raise ValueError("Hex colour code '#"+hex_colour+\
                         "' contains illegal characters!")
    else:
        # 1-digit numbers for r/g/b:
        if len(hex_colour) == 3:
            return tuple(int(hex_colour[i:i+1], 16)/15 for i in (0, 1, 2))
        # 2-digit numbers for r/g/b:
        elif len(hex_colour) == 6 :
            return tuple(int(hex_colour[i:i+2], 16)/255 for i in (0, 2, 4))
        # 3-digit numbers for r/g/b:
        elif len(hex_colour) == 9:
            return tuple(int(hex_colour[i:i+3], 16)/4095 for i in (0, 3, 6))
        else:
            raise ValueError("Hex colour code '#"+hex_colour+\
                             "' is of incorrect length!")