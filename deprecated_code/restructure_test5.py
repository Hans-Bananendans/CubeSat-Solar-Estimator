# -*- coding: utf-8 -*-
"""
Created on Mon May  3 08:57:58 2021

@author: Johan Monster

"""

import numpy as np
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
from cp_plotting import plot_global_tripod, plot_frame,\
    plot_vertex, plot_face, plot_geometry_perpendiculars, plot_illumination,\
        plot_A_ill


#%% ==== Plotting ====

# Toggle plotting functionality:
if True:
    
    # Setting up the plot:
    fig = plt.figure(figsize=(10, 7))
    ax = mp3d.Axes3D(fig)
    
    """ TO CHANGE THE DEFAULT CAMERA VIEW, CHANGE THESE: """
    ax.view_init(elev=20, azim=-50)
    
    steps = 64
    angle_step = d2r(360/steps)
    
    p0 = Vertex([  0.5,   0.5,  0.5])
    p1 = Vertex([-0.05, -0.05, -0.1])
    p2 = Vertex([-0.05, -0.05,  0.1])
    p3 = Vertex([ 0.05, -0.05,  0.1])
    p4 = Vertex([ 0.05, -0.05, -0.1])
    p5 = Vertex([-0.05,  0.05, -0.1])
    p6 = Vertex([-0.05,  0.05,  0.1])
    p7 = Vertex([ 0.05,  0.05,  0.1])
    p8 = Vertex([ 0.05,  0.05, -0.1])
    
    fA = Face(p4, p3, p2, p1)
    fB = Face(p2, p3, p7, p6)
    fC = Face(p3, p4, p8, p7)
    fD = Face(p4, p1, p5, p8)
    fE = Face(p1, p2, p6, p5)
    fF = Face(p5, p6, p7, p8)
    
    cubesat = Geometry([fA, fB, fC, fD, fE, fF])
    
    frame1 = Frame()
    frame1.add_geometry(cubesat)
    frame1.translate(0.3, 0.3, 0.25)
    
    frame1.rotate(0, 0, d2r(1.5*360/12))
    
    # DO A BARN DOOR!
    # cubesat.rotate(0,d2r(-90),0,cor=list(cubesat.find_cuboid_centroid()))
    
    
    # Model the fake Earth
    Ep = Face(Vertex([0.001,0.001,0]), Vertex([0.001,-0.001,0]),
              Vertex([-0.001,-0.001,0]), Vertex([-0.001,0.001,0]))
    E = Geometry([Ep])
    frameE = Frame()
    frameE.add_geometry(E)
    frameE.translate(0.3, 0.3, 0.25)
    # frameE.rotate(0, 0, d2r(1.5*360/12))
    
    E.translate(0,.15,0)
    
    
    # Pre-allocate the illumination area vector
    A_ill = np.zeros(steps-1)
    
    def update(i):
        
        # Setting up the axes object
        ax.clear()
        
        
        
        plotscale = 0.5
        ax.set_xlim(0, 1.25*plotscale)
        ax.set_ylim(0, 1.25*plotscale)
        ax.set_zlim(0, plotscale)
        
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
        
        
        # Transforming cubesat
        cubesat.rotate(angle_step,0,0,cor=list(cubesat.find_cuboid_centroid()))
        
        # Transform fake Earth
        E.rotate(angle_step,0,0,cor=list(cubesat.find_cuboid_centroid()))

        # Shitty shadow simulation:
        if True: # Toggle
            shadow = False  
            if 0.5-0.361/2 <= i/(steps-1) <= 0.5+0.361/2:
                shadow = True
            else:
                shadow = False
        
        global A_ill
        if not shadow:
            A_ill[i-1] = round(frame1.illuminated_area(cubesat, plane='xz'), 4)
        else:
            A_ill[i-1] = 0
        
        print("[DEBUG] Frame {}, i/(steps-1) = {}, shadow {}"\
              .format(i, round(i/(steps-1),2), shadow))
        
        # Plotting the global XYZ tripod
        plot_global_tripod(ax, scaling=plotscale/2)
        
        # Plot tripod of frame1:
        if not shadow:
            plot_frame(ax, frame1, tripod_scale=plotscale/8,
                       perpfill=False, perpscale=0.1, 
                       illumination=True, ill_plane='xz',
                       facecolour="#0F0F0F", facealpha=1
                       )
            
            plot_illumination(ax, frame1, plane='xz')
            
        else:
            plot_frame(ax, frame1, tripod_scale=plotscale/8,
                       perpfill=False, perpscale=0.1, 
                       illumination=False, ill_plane='xz',
                       facecolour="#0F0F0F", facealpha=1
                       )
            
            plot_illumination(ax, frame1, plane='xz',linecolour="#000")
        

        
        plot_frame(ax, frameE, tripod_scale=plotscale/8, show_tripod=False,
                   linefill=False, facefill=False,
                   illumination=False,
                   vertexcolour="#0000DD", vertexsize=300
                   )

        ax.set_title("3D Visualization.    Frame: {}.    A_ill = {} m^2.    Shadow is {}.\
                     ".format(str(i), str(A_ill[i-1]), shadow))
        
    ani = animation.FuncAnimation(fig, update, np.arange(1,steps), 
                                  interval = 100, repeat = False)

    plt.show()
