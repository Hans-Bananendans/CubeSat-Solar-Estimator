# -*- coding: utf-8 -*-
"""
Created on Sun May 16 19:40:08 2021

@author: Johan Monster

Cube Projector - Geogens
Version: 2.0


This file will contain utility functions for automatic creation of common/basic shapes.
Instead of creating geometries by individually specifying each vertex, these functions will generate the
geometries automatically from a set of much simpler input parameters.
"""

from cp_vertex import Vertex
from cp_face import Face
from cp_geometry import Geometry
from cp_vector import Vector
from cp_frame import Frame
from cp_utilities import r2d, d2r