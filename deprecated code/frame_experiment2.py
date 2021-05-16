# -*- coding: utf-8 -*-
"""
Created on Tue April 27 15:02:11 2021

@author: Johan Monster

Frame experiment 2
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
    
    steps = 128
    angle_step = d2r(360/steps)
    
    frame1 = Frame()

    projection_frame = Frame()
    
    p1 = Vertex(-0.05, -0.05, -0.1, frame1)
    p2 = Vertex(-0.05, -0.05, 0.1, frame1)
    p3 = Vertex(0.05, -0.05, 0.1, frame1)
    p4 = Vertex(0.05, -0.05, -0.1, frame1)
    p5 = Vertex(-0.05, 0.05, -0.1, frame1)
    p6 = Vertex(-0.05, 0.05, 0.1, frame1)
    p7 = Vertex(0.05, 0.05, 0.1, frame1)
    p8 = Vertex(0.05, 0.05, -0.1, frame1)
    
    fA = Face(p4, p3, p2, p1, frame1)
    fB = Face(p2, p3, p7, p6, frame1)
    fC = Face(p3, p4, p8, p7, frame1)
    fD = Face(p4, p1, p5, p8, frame1)
    fE = Face(p1, p2, p6, p5, frame1)
    fF = Face(p5, p6, p7, p8, frame1)
    
    
    def update(i):
        
        # Transforming frame1
        if i < steps/2:
            frame1.translate(1/steps,2/steps,2/steps)
        else:
            frame1.translate(1/steps,-2/steps,-2/steps)
        frame1.rotate(0,-2*np.pi/steps,-2*np.pi/steps)
        # frame1.rotate(0,0,0)
        
        # vertex1.rotate(2*np.pi/steps,0,0, cor=vertex2)
        
        # Setting up the axes object
        ax.clear()
        
        ax.set_title("Wireframe visualization. Frame: {}".format(str(i)))
        
        plotscale = 1
        ax.set_xlim(0, 1.25*plotscale)
        ax.set_ylim(0, 1.25*plotscale)
        ax.set_zlim(0, plotscale)
        # ax.set_xlim(-1, 1)
        # ax.set_ylim(-1, 1)
        # ax.set_zlim(-1, 1)
    
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
        
        # Plotting the global XYZ tripod
        plot_global_tripod(ax, scaling=plotscale/2)
        
        # Plot tripod of frame1:
        plot_frame(ax, frame1, tripod_scale=plotscale/8, facealpha=0.4)

        
        # Plot vertex1
        # plot_vertex(ax, vertex1)

    ani = animation.FuncAnimation(fig, update, np.arange(1,steps), 
                                  interval = 75, repeat = False)

    plt.show()
