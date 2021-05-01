# -*- coding: utf-8 -*-
"""
Created on Sat May 01 5:43:09 2021

@author: Johan Monster

Cube Projector - Vector class
Version: 2.0

"""

import numpy as np
from numpy import sin, cos


class Vector:
    """Custom vector implementation, mainly for display purposes.
    
    TODO: Implement vector class.
    """

    def __init__(self, head, tail=None, parent=None):
        self.head = None # Vertex or None
        self.tail = None # Vertex or None
        self.parent = None # Frame or None