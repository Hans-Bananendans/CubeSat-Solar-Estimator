# -*- coding: utf-8 -*-
"""
Created on Tue April 27 15:02:11 2021    

@author: Johan Monster

Frame experiment 1
"""

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


# Toggle plotting functionality:
if True:
    
    # Setting up the plot:
    fig = plt.figure(figsize=(10, 7))
    ax = mp3d.Axes3D(fig)

    """ TO CHANGE THE DEFAULT CAMERA VIEW, CHANGE THESE: """
    ax.view_init(elev=20, azim=-60)
    
    steps = 40
    angle_step = d2r(360/steps)
    
    frame1 = Frame()
    
    vertex1 = Vertex(0.5,0,0,parent=frame1)
    vertex2 = Vertex(0.5,0.5,0,parent=frame1)
    
    def update(i):
        
        # Transforming frame1
        frame1.translate(2/steps,2/steps,2/steps)
        # frame1.rotate(0,0,-2*np.pi/(steps))
        
        # Local transformation of vertex1 - rotate around vertex2
        vertex1.rotate(2*np.pi/steps,0,0, cor=vertex2)
        
        # Setting up the axes object
        ax.clear()
        
        ax.set_title("Wireframe visualization. Frame: {}".format(str(i)))
        
        ax.set_xlim(0, 0.4*10)
        ax.set_ylim(0, 0.4*10)
        ax.set_zlim(0, 0.3*10)
    
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
        
        # Plotting the global XYZ tripod
        plot_global_tripod(ax)
        
        # Plot tripod of frame1:
        plot_frame(ax, frame1)

    ani = animation.FuncAnimation(fig, update, np.arange(1,steps), 
                                  interval = 75, repeat = False)

    plt.show()
