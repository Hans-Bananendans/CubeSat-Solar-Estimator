# -*- coding: utf-8 -*-
"""
Created on Sat May 01 6:39:27 2021   

@author: Johan Monster

Empty sandbox for testing Cube Projector functionality

"""

#%% ==== Import dependencies ====

import numpy as np
from copy import deepcopy
from numpy import sin, cos
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import mpl_toolkits.mplot3d as mp3d

from cp_vertex import Vertex
from cp_face import Face
from cp_geometry import Geometry
from cp_vector import Vector
from cp_frame import Frame
from cp_utilities import d2r#, r2d
from cp_plotting import plot_global_tripod, plot_frame, \
    plot_vertex, plot_face


#%% ==== Manipulate objects ====

frame1 = Frame()
frame1.translate(1,1,1)

frame2 = Frame()

p1 = Vertex(0,0,0)
# p2 = Vertex(1,0,0)
# p3 = Vertex(1,1,0)
# p4 = Vertex(0,1,0)
# p5 = Vertex(0,1,1)
# p6 = Vertex(0,0,1)

# f1 = Face(p1, p2, p3, p4, frame1)
# f2 = Face(p1, p4, p5, p6, frame1)

f1 = Face(Vertex(0,0,0), Vertex(1,0,0), Vertex(1,1,0), Vertex(0,1,0))








#%% ==== Plotting ====

# Toggle plotting functionality:
if False:
 
    # Setting up the plot:
    fig = plt.figure(figsize=(8, 8))
    ax = mp3d.Axes3D(fig)

    """ TO CHANGE THE DEFAULT CAMERA VIEW, CHANGE THESE: """
    ax.view_init(elev=20, azim=-60)
    
    # Setting up the axes object        
    ax.set_title("Wireframe visualization.")
    
    ax.set_xlim(-2, 2)
    ax.set_ylim(-2, 2)
    ax.set_zlim(-2, 2)

    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    
    # Plotting the global XYZ tripod
    plot_global_tripod(ax)
    
    # ==== Plot custom objects ====

    for vertex in [p1, p2, p3, p4, p5, p6]:
        plot_vertex(ax, vertex)
    
    for face in [f1, f2]:
        plot_face(ax, face)
        
    plot_frame(ax, frame1)

    # =============================
    
    plt.show()
