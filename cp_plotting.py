# -*- coding: utf-8 -*-
"""
Created on Sat May 01 5:43:09 2021 

@author: Johan Monster

Cube Projector - Plotting functions
Version: 2.0

"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import mpl_toolkits.mplot3d as mp3d

from cp_vertex import Vertex
from cp_face import Face
from cp_geometry import Geometry
from cp_vector import Vector
from cp_frame import Frame
     
def plot_vertex(axes: plt.matplotlib.axes, vertex: Vertex, 
                colour="#000", size=10):
    
    x_global, y_global, z_global = vertex.global_coordinates()
    
    axes.scatter(x_global, 
                 y_global,
                 z_global,
                 c=colour, s=size)

def plot_face(axes: plt.matplotlib.axes, face: Face, 
              fill=1, alpha=0.9,
              linecolour='black', linealpha=1,
              illumination=False, ill_value=0):
    """ Plots an individual Face object. Check source code for kwargs.
    """
    (xlist, ylist, zlist) = face.plotlist()

    # Plot individual vertices:
    # axes.scatter(xlist, ylist, zlist, c=linecolour, s=10)

    # Plot edges that connect vertices:
    for i in range(-1, len(xlist) - 1):
        axes.plot3D([xlist[i], xlist[i + 1]],
                    [ylist[i], ylist[i + 1]],
                    [zlist[i], zlist[i + 1]],
                    linecolour, alpha=linealpha, lw=2)

    # Plot the face surface:
    if fill == 1:

        fs = mp3d.art3d.Poly3DCollection([face.plotlist2()],
                                         linewidth=0)
        fs.set_facecolor((0.1, 0.1, 0.1, alpha))
        axes.add_collection3d(fs)
        
# def plot_geometry(axes: plt.matplotlib.axes, geometry: Geometry, 
#                   fill=0, alpha=0.2,
#                   linecolour='black', linealpha=1,
#                   illumination=False, illumination_plane='xy'):
#     """Plot geometry by individually plotting its faces. If illumination 
#        is turned on, uses illuminated_faces() method to fetch which faces 
#        are illuminated, and passes this to the plot_face() function.
       
#        TODO: Allow passing of 'plane'
#        """
#     if illumination == True:
#         ill_faces = geometry.illuminated_faces(plane=illumination_plane)
#         for idx, face in enumerate(geometry.faces):
#             plot_face(axes, face, fill=fill, alpha=alpha,
#                       linecolour=linecolour, linealpha=linealpha,
#                       illumination=True, ill_value=ill_faces[idx])
#     else:
#         for idx, face in enumerate(geometry.faces):
#             plot_face(axes, face, fill=fill, alpha=alpha,
#                       linecolour=linecolour, linealpha=linealpha)

# def plot_geometry_perpendiculars(axes, geometry: Geometry, \
#                                  colour='gray', alpha=.5):
#     """Plots the normals/perpendiculars of each face in a geometry, and 
#        displays them as little gray arrows.
#        TODO: Implement scaling."""
#     # First plot centroids of each plane:
#     xc = []
#     yc = []
#     zc = []

#     centroids = geometry.find_fcentroids()

#     for vertex in centroids:
#         xc.append(vertex.x)
#         yc.append(vertex.y)
#         zc.append(vertex.z)
#     axes.scatter(xc, yc, zc, c=colour, s=5, alpha=alpha)

#     # Then, attach plane-perpendicular arrows to the centroids:
#     xyzuvw = geometry.perpendiculars_plotlist()
#     axes.quiver(xyzuvw[0][:], xyzuvw[1][:], xyzuvw[2][:],
#                 xyzuvw[3][:], xyzuvw[4][:], xyzuvw[5][:],
#                 color=colour, alpha=alpha)

def plot_frame_tripod(axes: plt.matplotlib.axes, frame: Frame, 
                      alpha=1, scaling=1.):
    """Plots an rgb xyz tripod at the global origin.
        - alpha changes the opacity of the arrows. (default: 1)
        - scaling changes the length of the arrows (default: 1)
        """
    sc = scaling
    axes.quiver(frame.x, frame.y, frame.z, 
                sc*frame.xdir[0], sc*frame.xdir[1], sc*frame.xdir[2],
              arrow_length_ratio=0.15, color='red', alpha=alpha)
    axes.quiver(frame.x, frame.y, frame.z, 
                sc*frame.ydir[0], sc*frame.ydir[1], sc*frame.ydir[2],
              arrow_length_ratio=0.15, color='green', alpha=alpha)
    axes.quiver(frame.x, frame.y, frame.z, 
                sc*frame.zdir[0], sc*frame.zdir[1], sc*frame.zdir[2],
              arrow_length_ratio=0.15, color='blue', alpha=alpha)

def plot_frame(axes: plt.matplotlib.axes, frame: Frame, show_tripod=True):
    
    if show_tripod:
        plot_frame_tripod(axes, frame, scaling=0.25)
    
    for vertex in frame.vertices:
        plot_vertex(axes, vertex)
    
    for face in frame.faces:
        plot_face(axes, face)

def plot_global_tripod(axes: plt.matplotlib.axes, alpha=1, scaling=1.):
    """Plots an rgb xyz tripod at the global origin.
        - alpha changes the opacity of the arrows. (default: 1)
        - scaling changes the length of the arrows (default: 1)
        """
    axes.quiver(0, 0, 0, scaling, 0, 0,
              arrow_length_ratio=0.15, color='red', alpha=alpha)
    axes.quiver(0, 0, 0, 0, scaling, 0,
              arrow_length_ratio=0.15, color='green', alpha=alpha)
    axes.quiver(0, 0, 0, 0, 0, scaling,
              arrow_length_ratio=0.15, color='blue', alpha=alpha)