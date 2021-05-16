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
from cp_plotting import plot_global_tripod, plot_frame


#%% ==== Manipulate objects ====
           













#%% ==== Plotting ====

# Toggle plotting functionality:
if True:
 
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






    # =============================
    
    plt.show()
