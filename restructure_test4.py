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
from cp_plotting import plot_global_tripod, plot_frame, \
    plot_vertex, plot_face, plot_geometry_perpendiculars, plot_illumination

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
frame1.translate(0.5, 0.5, 0.5)
frame1.rotate(0, 0, d2r(1*360/12))

cubesat.rotate(d2r(-45),0,0,cor=list(cubesat.find_cuboid_centroid()))

f = Face(Vertex([0,0,0]),Vertex([1,0,0]),Vertex([1,1,0]),Vertex([0,1,0]))

for face in frame1.illuminated_faces(cubesat, 'xz'):
    print(frame1.illuminated_strength(face, 'xz'))


#%% ==== Plotting ====

# Toggle plotting functionality:
if True:

    # Setting up the plot:
    fig = plt.figure(figsize=(8, 8))
    ax = mp3d.Axes3D(fig)
    
    """ TO CHANGE THE DEFAULT CAMERA VIEW, CHANGE THESE: """
    ax.view_init(elev=20, azim=-60)
    
    # Setting up the axes object        
    ax.set_title("Wireframe visualization. Illuminated area: {} m^2."
                 .format(round(frame1.illuminated_area(cubesat, plane='xz'),4)))
    
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_zlim(0, 1)
    
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    
    # Plotting the global XYZ tripod
    plot_global_tripod(ax)
    
    # ==== Plot custom objects ====
    
    # plot_vertex(ax, p0)
    # plot_face(ax, f)    
    
    # plot_face(ax, fA, frame1.plotlist(fA))
    plot_frame(ax, frame1, tripod_scale=0.2, \
               illumination=True, ill_plane='xz', \
               facecolour="#468", facealpha=0.95, \
               perpfill=True, perpscale=0.1, perpalpha=0.25)
        
    plot_illumination(ax, frame1, plane='xz')

    
    # =============================
    
    plt.show()
