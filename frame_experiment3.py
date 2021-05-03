# -*- coding: utf-8 -*-
"""
Created on Mon May  3 08:57:58 2021

@author: Johan Monster

Frame experiment 3
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
    frame1.translate(0.5,0.5,0.5)
    frame1.rotate(0,d2r(7*360/12),d2r(1*360/12))


    
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
    
    cubesat = Geometry(frame1)
    cubesat.add_faces([fA, fB, fC, fD, fE, fF])
    
    def update(i):
        
        # Transforming frame1
        cubesat.rotate(angle_step,0,0,cor=cubesat.make_cuboid_centroid())
        
        projection_frame = Frame()
        
        for face in cubesat.faces:
            face.project(new_frame=projection_frame, plane='xz') 
        
        
        print("[DEBUG] Plotting frame {}".format(i))
        
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
        plot_frame(ax, frame1, tripod_scale=plotscale/8,facealpha=0.3)
        plot_frame(ax, projection_frame, tripod_scale=plotscale/8,
                   facefill=False, linecolour="#AA2")

        
        # Plot vertex1
        # plot_vertex(ax, vertex1)

    ani = animation.FuncAnimation(fig, update, np.arange(1,steps), 
                                  interval = 75, repeat = True)

    plt.show()
