# -*- coding: utf-8 -*-
"""
Created on Thu April 27 15:02:11 2021    

@author: Johan Monster

Cube Projector - Vertex class
Version: 2.0

"""

import numpy as np
from numpy import sin, cos
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import mpl_toolkits.mplot3d as mp3d

class Frame:
    """Custom implentation of a 3D point expressed in cartesians."""

    def __init__(self, x=0., y=0., z=0.):
        self.x = x
        self.y = y
        self.z = z
        self.xdir = np.array([1,0,0])
        self.ydir = np.array([0,1,0])
        self.zdir = np.array([0,0,1])

    def translate(self, dx=0., dy=0., dz=0.):
        self.x += dx 
        self.y += dy
        self.z += dz
        
    def rotate(self, a, b, c, seq='321'):
        """ 
        """

        Rx = np.array([[1,      0,       0],
                       [0, cos(a), -sin(a)],
                       [0, sin(a),  cos(a)]])

        Ry = np.array([[ cos(b), 0, sin(b)],
                       [      0, 1,      0],
                       [-sin(b), 0, cos(b)]])

        Rz = np.array([[cos(c), -sin(c), 0],
                       [sin(c),  cos(c), 0],
                       [     0,       0, 1]])

        for axis in reversed(seq):
            temp = (self.x, self.y, self.z)
            self.translate(-1*temp[0], -1*temp[1], -1*temp[2])

            if axis == '1':
                self.xdir = np.dot(Rx, self.xdir) # Redundant?
                self.ydir = np.dot(Rx, self.ydir)
                self.zdir = np.dot(Rx, self.zdir)
            if axis == '2':
                self.xdir = np.dot(Ry, self.xdir)
                self.ydir = np.dot(Ry, self.ydir) # Redundant?
                self.zdir = np.dot(Ry, self.zdir)
            if axis == '3':
                self.xdir = np.dot(Rz, self.xdir)
                self.ydir = np.dot(Rz, self.ydir)
                self.zdir = np.dot(Rz, self.zdir) # Redundant?

            self.translate(1*temp[0], 1*temp[1], 1*temp[2])


class Vertex:
    """Custom implentation of a 3D point expressed in cartesians."""

    def __init__(self, x=0., y=0., z=0., parent:Frame=None):
        self.x = x
        self.y = y
        self.z = z
        self.parent = parent

    def translate(self, dx=0., dy=0., dz=0.):
        self.x += dx
        self.y += dy
        self.z += dz

def plot_frame(axes, frame, alpha=1, scaling=1.):
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

def plot_global_tripod(axes, alpha=1, scaling=1.):
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
    
def r2d(a):
    """Convert radians to degrees."""
    return a * 180 / np.pi

def d2r(a):
    """Convert degrees to radians."""
    return a * np.pi / 180

# Toggle plotting functionality:
if True:
    
    # Setting up the plot:
    fig = plt.figure(figsize=(10, 7))
    ax = mp3d.Axes3D(fig)

    """ TO CHANGE THE DEFAULT CAMERA VIEW, CHANGE THESE: """
    ax.view_init(elev=20, azim=-60)
    
    steps = 180
    angle_step = d2r(360/steps)
    
    frame1 = Frame()
    frame1.translate(2,2,2)
    
    # Pre-allocate illuminated area array
    A_illuminated = np.zeros(steps)
    
    def update(i):
        
        # Transforming frame1
        # frame1.translate(2/steps,2/steps,2/steps)
        frame1.rotate(np.pi/(2*steps),0,np.pi/(steps))

        # Setting up the axes object
        ax.clear()
        
        ax.set_title("Wireframe visualization.")
        
        ax.set_xlim(0, 0.4*10)
        ax.set_ylim(0, 0.4*10)
        ax.set_zlim(0, 0.3*10)
    
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
        
        # Plotting the global XYZ tripod
        plot_global_tripod(ax)
        
        # Plot tripod of frame1:
        plot_frame(ax, frame1, scaling=0.25)


    ani = animation.FuncAnimation(fig, update, np.arange(steps), interval = 75, repeat = False)

    plt.show()
