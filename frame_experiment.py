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

class Vertex:
    """Custom implentation of a 3D point expressed in cartesians."""

    def __init__(self, x=0., y=0., z=0., parent=None):
        self.x = x
        self.y = y
        self.z = z
        self.parent = parent
        self.connect_parent(parent)
    
    def remove_parent(self):
        """If Vertex has a parent frame, removes it as a parent, and ensures 
           it is removed from the internal list in the frame.
           """
        if self.parent != None:
            self.parent.remove_vertex(self)
            self.parent = None

    def connect_parent(self, new_parent):
        """Connects vertex to a frame.
           """
        # Connect new parent:
        self.parent = new_parent
        self.parent.add_vertex(self)
        
    def change_parent(self, new_parent):
        """Connects vertex to another frame. If vertex was already attached 
           to a frame, it undoes this first.
           """
        # Remove old parent first (if applicable):
        if self.parent != None:
            self.remove_parent()
        
        # Connect new parent:
        self.parent = new_parent
        self.parent.add_vertex(self)
    
    def translate(self, dx=0., dy=0., dz=0.):
        self.x += dx
        self.y += dy
        self.z += dz
    
    def global_coordinates(self):
        xyz_local = np.array([self.x, self.y, self.z])
        o_local = np.array([self.parent.x, self.parent.y, self.parent.z])
        
        # xyz_global = np.dot(np.transpose(self.parent.dcm), xyz_local)
        xyz_global = np.dot(self.parent.dcm, xyz_local)
        xyz_global += o_local
        
        return xyz_global[0], xyz_global[1], xyz_global[2]
        
    def update_translation(self):
        pass
    
    def rotate(self, a, b, c, cor=None, seq='321'):
        """ Rotates the vertex around a point in space.
            a is Euler angle of rotation around x, etc...
                expressed in radians
            seq is a string with the rotation sequence, e.g. '321' for:
                Rz(c).Ry(b).Rx(a).vertex
            cor is the centre of rotation, which MUST be specified as
                a vertex.
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
            
            # If no Centre of Rotation (cor) is given, take frame origin.
            # If not, subtract coordinates of cor first before applying 
            # the rotation, and then reverse this afterwards.
            if cor != None:
                self.translate(-1*cor.x, -1*cor.y, -1*cor.z)
            
            if axis == '1':
                (self.x, self.y, self.z) = \
                np.dot(Rx, np.array([self.x, self.y, self.z]))
                
            if axis == '2':
                (self.x, self.y, self.z) = \
                np.dot(Ry, np.array([self.x, self.y, self.z]))
            if axis == '3':
                (self.x, self.y, self.z) = \
                np.dot(Rz, np.array([self.x, self.y, self.z]))
            
            # Reverse temporary translation.
            if cor != None:
                self.translate(1*cor.x, 1*cor.y, 1*cor.z)
            
    def __str__(self):
        """Print the vertex coordinates."""
        return str([round(self.x, 5),
                    round(self.y, 5),
                    round(self.z, 5)])
    
    def xyz(self):
        """Returns a numpy array with the local x, y, z coordinates."""
        return np.array([self.x, self.y, self.z])
    
    def xyz_global(self):
        """Returns a numpy array with the global x, y, z coordinates."""
        return np.array(self.global_coordinates())
        
    def readout(self, dec=4):
        """More detailed information about Vertex.
           Use dec to define the number of decimals before rounding.
           Set dec to -1 to disable rounding.
           """
        if dec == -1:  # Set dec to -1 to disable rounding
            localcoords = self.xyz()
            globalcoords = self.xyz_global()
            parentorigin = self.parent.origin()
        else:
            localcoords = np.round(self.xyz(), dec)
            globalcoords = np.round(self.xyz_global(), dec)
            parentorigin = np.round(self.parent.origin(), dec)
        
        
        print("Parent frame: {}".format(self.parent))
        print("Parent frame origin: {}".format(parentorigin))
        print("Local coordinates:   {}".format(localcoords))
        print("Global coordinates:  {}".format(globalcoords))


class Frame:
    """Custom implentation of a 3D point expressed in cartesians."""

    def __init__(self, x=0., y=0., z=0.):
        self.x = x
        self.y = y
        self.z = z
        self.xdir = np.array([1,0,0])
        self.ydir = np.array([0,1,0])
        self.zdir = np.array([0,0,1])
        
        self.dcm = None
        self.recalculate_dcm()
        
        self.vertices = []
        self.faces = []
        self.geometries = []
        self.vectors = []

    def add_vertex(self, vertex: Vertex):
        self.vertices.append(vertex)

    def remove_vertex(self, vertex: Vertex):
        self.vertices.remove(vertex)
    
    def origin(self):
        return np.array([self.x, self.y, self.z])
    
    def recalculate_dcm(self):
        gx = np.array([1,0,0])
        gy = np.array([0,1,0])
        gz = np.array([0,0,1])
        self.dcm = np.array([
            [np.dot(gx, self.xdir), np.dot(gx, self.ydir), np.dot(gx, self.zdir)],
            [np.dot(gy, self.xdir), np.dot(gy, self.ydir), np.dot(gy, self.zdir)],
            [np.dot(gz, self.xdir), np.dot(gz, self.ydir), np.dot(gz, self.zdir)]])
        
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
        
        self.recalculate_dcm()
    
    def readout(self, dec=4):
        """More detailed information about Vertex.
           Use dec to define the number of decimals before rounding.
           Set dec to -1 to disable rounding.
           """
        if dec == -1:  # Set dec to -1 to disable rounding
            origin = self.origin()
        else:
            origin = np.round(self.origin(), dec)
        
        children = len(self.vertices)   + \
                   len(self.faces)      + \
                   len(self.geometries) + \
                   len(self.vectors)
        
        print("Carthesian coordinate frame.")
        print("Frame origin: {}\n".format(origin))
        print("Children: {}".format(children))
        print("Associated Vertices:   {}".format(len(self.vertices)))
        print("Associated Faces:      {}".format(len(self.faces)))
        print("Associated Geometries: {}".format(len(self.geometries)))
        print("Associated Vectors:    {}".format(len(self.vectors)))
             
     
def plot_vertex(axes, vertex, colour="#000", size=10):
    
    x_global, y_global, z_global = vertex.global_coordinates()
    
    axes.scatter(x_global, 
                 y_global,
                 z_global,
                 c=colour, s=size)

def plot_frame_tripod(axes, frame, alpha=1, scaling=1.):
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

def plot_frame(axes, frame, show_tripod=True):
    
    if show_tripod:
        plot_frame_tripod(axes, frame, scaling=0.25)
    
    for vertex in frame.vertices:
        plot_vertex(axes, vertex)

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
    
    steps = 40
    angle_step = d2r(360/steps)
    
    frame1 = Frame()
    # frame1.translate(2,2,2)
    
    vertex1 = Vertex(0.5,0,0,parent=frame1)
    vertex2 = Vertex(0.5,0.5,0,parent=frame1)
    
    # Pre-allocate illuminated area array
    # A_illuminated = np.zeros(steps)
    
    def update(i):
        
        # Transforming frame1
        frame1.translate(2/steps,2/steps,2/steps)
        # frame1.rotate(0,0,-2*np.pi/(steps))
        # frame1.rotate(0,0,0)
        
        vertex1.rotate(2*np.pi/steps,0,0, cor=vertex2)
        
        # Setting up the axes object
        ax.clear()
        
        ax.set_title("Wireframe visualization. Frame: {}".format(str(i)))
        
        ax.set_xlim(0, 0.4*10)
        ax.set_ylim(0, 0.4*10)
        ax.set_zlim(0, 0.3*10)
        # ax.set_xlim(-1, 1)
        # ax.set_ylim(-1, 1)
        # ax.set_zlim(-1, 1)
    
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
        
        # Plotting the global XYZ tripod
        plot_global_tripod(ax)
        
        # Plot tripod of frame1:
        plot_frame(ax, frame1)

        
        # Plot vertex1
        # plot_vertex(ax, vertex1)

    ani = animation.FuncAnimation(fig, update, np.arange(1,steps), interval = 75, repeat = False)

    plt.show()
