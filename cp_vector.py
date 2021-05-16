# -*- coding: utf-8 -*-
"""
Created on Sat May 01 5:43:09 2021

@author: Johan Monster

Cube Projector - Vector class
Version: 2.0

"""

# import numpy as np
# from numpy import sin, cos
from cp_vertex import Vertex


class Vector:
    """Custom vector implementation. Unfinished.
    
    TODO: Implement translation.
    TODO: Implement rotation.
    TODO: Implement length/norm.
    TODO: Implement various numpy linalg compatibility.
    TODO: ...
    """

    def __init__(self, head: Vertex, tail:Vertex=Vertex([0,0,0]), 
                 parenttype="global"):
        
        self.head = head # Vertex
        self.tail = tail # Vertex
        
        self.parenttype = "global"
        
        # Re-assign self.parenttype to check for validity of given parenttype
        if parenttype != "global":
            self.set_parenttype(parenttype)
        
        # Let child vertices know their parenttype is 'vector'
        for vertex in self.vertices():
            vertex.set_parenttype("vector")
     
        
    def set_parenttype(self, new_parenttype):
        if new_parenttype in ["global", "frame"]:
            self.parenttype=new_parenttype
        else:
            raise ValueError("The parenttype cannot be anything other than: \
                             'global', 'frame'!")
    
    def vertices(self):
        """Returns a list of the two vertices defining the vector."""
        return [self.head, self.tail]
    