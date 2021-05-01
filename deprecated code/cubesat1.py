# -*- coding: utf-8 -*-
"""
Created on Thu April 27 15:02:11 2021     

@author: Johan Monster

Version: 1.2

TODO: Build tools to model dynamic behaviour
 - For Earth-facing spacecraft: 
     
     J . d/dt(omega) + omega x J . omega = 3n^2 * a3 xJ.a3 
     
     [ tdot1 ]    1  [ct2 st1*st2  ct1*st2] [omega1]    n  [   st3   ]
     [ tdot2 ] = --- [ 0  ct1*ct2 -st1*ct2] [omega2] + --- [ ct2*ct3 ]
     [ tdot3 ]   ct2 [ 0    st1      ct1  ] [omega3]   ct2 [ st2*st3 ]
     
TODO: Add ways to compute power using solar flux and incidence angle
TODO: Add second point of light as albedo source
TODO: Add orbital propagation of point mass

"""

import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d as mp3d

from cp_utilities import r2d, d2r, q2e, e2q
from cp_vertex import Vertex
from cp_face import Face
from cp_geometry import Geometry
from cp_plotting import plot_xyz_tripod, plot_geometry, \
    plot_geometry_perpendiculars

# %% Defining the body

""" === CHANGE ROTATION OF BODY HERE === """

# Note: choosing bdeg=90 will result in gimbal lock. 
adeg = 30 # Rotation around x axis in degrees
bdeg = 40  # Rotation around y axis in degrees
cdeg = 0   # Rotation around z axis in degrees

""" ==================================== """

a = d2r(adeg)
b = d2r(bdeg)
c = d2r(cdeg)


"""Basic Cubesat model:"""
p1 = Vertex(0.0, 0.0, 0.0)
p2 = Vertex(0.0, 0.0, 0.2)
p3 = Vertex(0.1, 0.0, 0.2)
p4 = Vertex(0.1, 0.0, 0.0)
p5 = Vertex(0.0, 0.1, 0.0)
p6 = Vertex(0.0, 0.1, 0.2)
p7 = Vertex(0.1, 0.1, 0.2)
p8 = Vertex(0.1, 0.1, 0.0)


# Manually calculate centroid of Cubesat geometry
# COR = Vertex(0.05, 0.05, 0.1)
# COR.translate(0.1, 0.1, 0.1)

# Define faces of Cubesat
fA = Face(p4, p3, p2, p1)
fB = Face(p2, p3, p7, p6)
fC = Face(p3, p4, p8, p7)
fD = Face(p4, p1, p5, p8)
fE = Face(p1, p2, p6, p5)
fF = Face(p5, p6, p7, p8)

# Initialize Cubesat model as a Geometry object.
geometry1 = Geometry()
geometry1.add_faces([fA, fB, fC, fD, fE, fF])

# Translate geometry outward, so it is not positioned on the global axes.
geometry1.translate(0.1, 0.1, 0.1)

# Apply rotations as specified earlier.
geometry1.rotate_cuboid_centroid(a, b, c)

# Create a projection object that can be plotted later.
geometry1xy = geometry1.project_illuminated_faces('xy')


# %% Plotting code

# Check value for illuminated area. Because Cubesat is a cuboid, illuminated
# area and the unilluminated area should be the same UNLESS only one side is 
# facing the light. Used for debugging.
A_shadow_check = round(geometry1xy.area(), 4)

# Properly compute the illuminated area using Geometry method:
A_shadow = round(geometry1.illuminated_area(plane='xy'), 4)

print('Coarse method: A_xy =', A_shadow_check, 'm^2')
print('Projected A_xy =', A_shadow, "m^2")

# Toggle plotting functionality:
if True:
    # Setting up the plot:

    fig = plt.figure(figsize=(10, 7))
    ax = mp3d.Axes3D(fig)

    """ TO CHANGE THE DEFAULT CAMERA VIEW, CHANGE THESE: """
    ax.view_init(elev=20, azim=-90)

    # # Forcing aspect ratio for Axes3D object is still fucked up in
    # # matplotlib source, even after many years. :(
    # ax.set_aspect('equal')
    # ax.set_box_aspect((1,1,1))
    
    ax.set_title("Wireframe visualization, A_xy = {} m^2".format(A_shadow))

    ax.set_xlim(0, 0.4)
    ax.set_ylim(0, 0.4)
    ax.set_zlim(0, 0.3)

    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')

    # Plotting the XYZ tripod
    plot_xyz_tripod(ax, scaling=0.25)

    # Plotting the centroid of the geometry (manually)
    COR = geometry1.find_cuboid_centroid()
    ax.scatter(COR.x, COR.y, COR.z, c='cyan', s=20, alpha=0.5)

    # Plot Cubesat model
    plot_geometry(ax, geometry1, linecolour='black', fill=1, alpha=1,
                  illumination=True)
    # Plot projection of Cubesat
    plot_geometry(ax, geometry1xy, linecolour='orange')

    # Highlight one of the faces in bright yellow (useful for debugging)
    # plot_face(ax, geometry1.faces[5], colour='yellow')    

    # Plotting the perpendiculars
    plot_geometry_perpendiculars(ax, geometry1)

    fig.show()
