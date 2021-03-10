# -*- coding: utf-8 -*-
"""
Created on Thu Feb 25 12:55:55 2021

Sources:
https://math.stackexchange.com/questions/2725780/orthogonal-projection-area-of-a-3-d-cube
    

@author: Johan Monster


Tests:
Minimal area: 
    a = 0, b = 0        -> A = 1
    
Max area with one rotation:
    a = 45, b = 0       -> A = sqrt(2)
    
Max area with any rotations:
    a = 45, b = 35.26   -> A = sqrt(3)
"""

import numpy as np
from numpy import sin, cos
import matplotlib.pyplot as plt

adeg = 45  # Rotation around x axis
bdeg = 35.3  # Rotation around y axis
cdeg = 0

a = np.pi/180*adeg
b = np.pi/180*bdeg
c = np.pi/180*cdeg

AV = np.array([[0.0, 0.0, 0.0],
               [0.1, 0.0, 0.0],
               [0.1, 0.1, 0.0],
               [0.0, 0.1, 0.0],
               [0.0, 0.0, 0.2],
               [0.1, 0.0, 0.2],
               [0.1, 0.1, 0.2],
               [0.0, 0.1, 0.2]])
choice = [0, 4, 5, 6, 2, 3]

# Defining six prominent vertices of the shape.
# Note that these vertices should be chosen such that the final rotated
# set of vertices produces a convex hull.
# V = np.array([[0, 0, 1],
#               [1, 0, 1],
#               [1, 0, 0],
#               [1, 1, 0],
#               [0, 1, 0],
#               [0, 1, 1]])

V = AV[choice[0],:]
for i in range(1, len(choice)):
    V = np.vstack([V, AV[choice[i],:]])

# Transpose V
V = np.transpose(V)


# Defining rotation matrices:
Rx = np.array([[1,      0,       0], 
               [0, cos(a), -sin(a)],
               [0, sin(a),  cos(a)]])

Ry = np.array([[ cos(b), 0, sin(b)], 
               [      0, 1,       0],
               [-sin(b), 0, cos(b)]])

Rz = np.array([[cos(c), -sin(c), 0], 
               [sin(c),  cos(c), 0],
               [     0,       0, 1]])


# Rotate object by computing Rz.Ry.Rx.V
V_r = np.dot(Rz, np.dot(Ry, np.dot(Rx, V)))



#%% Using the shoelace algorithm, compute area of the projection:

n = len(V_r[0,:])
lace1 = V_r[0,n-1]*V_r[1,0]
lace2 = V_r[0,0]*V_r[1,n-1]

for i in range(n-1):
    lace1 += V_r[0,i]*V_r[1,i+1]
    lace2 += V_r[0,i+1]*V_r[1,i]
    
A = 0.5 * np.abs(lace1 - lace2)
    
# Quick dirty vertex plot
if True:
    xpoints = list(V_r[0,:])
    xpoints.append(V_r[0,0])
    ypoints = list(V_r[1,:])
    ypoints.append(V_r[1,0])
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_aspect(1)
    ax.plot(xpoints, ypoints)
    fig.show

