# -*- coding: utf-8 -*-
"""
Created on Thu April 27 15:02:11 2021    

@author: Johan Monster

Cube Projector - Vertex class
Version: 2.0

"""

import numpy as np
from copy import deepcopy
from numpy import sin, cos
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import mpl_toolkits.mplot3d as mp3d

class Vertex:
    """Custom implentation of a 3D point expressed in cartesians.
    
    TODO: Implement ability to project points
    """

    def __init__(self, x=0., y=0., z=0., parent=None):
        self.x = x
        self.y = y
        self.z = z
        self.connect_parent(parent)
    
    def remove_parent(self):
        """If Vertex has a parent frame, removes it as a parent, and ensures 
           it is removed from the internal list in the frame.
           """
        if self.parent != None:
            self.parent.remove_vertex(self)
            self.parent = None

    def connect_parent(self, new_parent):
        """Connects vertex to a frame, unless:
            - the given new_parent is None
            - the vertex is already in the frame (duplicate control)
           """
        # Set own parent to new_parent:
        self.parent = new_parent
        
        # Ensure new_parent is not None:
        if new_parent != None:
            # Ensure vertex is not already connected to parent
            if self not in self.parent.vertices:
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
    
    def global_coordinates(self):
        xyz_local = np.array([self.x, self.y, self.z])
        o_local = np.array([self.parent.x, self.parent.y, self.parent.z])
        
        # xyz_global = np.dot(np.transpose(self.parent.dcm), xyz_local)
        xyz_global = np.dot(self.parent.dcm, xyz_local) + o_local
        # xyz_global = xyz_global + o_local
        
        return xyz_global[0], xyz_global[1], xyz_global[2]
    
    def translate(self, dx=0., dy=0., dz=0.):
        self.x += dx
        self.y += dy
        self.z += dz
    
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

class Face:
    """Class for a quadrilateral face (4-point polygon).
    
       A face consists of four vertices that make up the corners of a plane.
       All four vertices must be coplanar, and the class constructor will
       verify this. The order in which the vertices are supplied DOES matter,
       and for best results they must be supplied sequentially and 
       counterclockwise, i.e.:
          4  o---o 3
             |   |  
          1  o---o 2
          
       TODO: Add vertex cleanup?
       TODO: Add projection
       """

    def __init__(self, p1: Vertex, p2: Vertex, p3: Vertex, p4: Vertex, 
                 parent=None):
        self.p1 = p1
        self.p2 = p2
        self.p3 = p3
        self.p4 = p4
        if self.check_coplanarity() == False:
            raise ValueError("Error! Set of four points is not coplanar!")
            
        self.area = self.calc_area()
        
        self.parent = parent
        self.connect_parent(parent)
    
    def remove_parent(self):
        """If Vertex has a parent frame, removes it as a parent, and ensures 
           it is removed from the internal list in the frame.
           
           TODO: Check if vertices aren't used in other faces.
                 If not, remove that vertex as well.
           """
        if self.parent != None:
            self.parent.remove_face(self)
            self.parent = None

    def connect_parent(self, new_parent):
        """Connects vertex to a frame.
        
           TODO: Add points to parent too.
           """
        # Connect new parent:
        if self.parent != None:
            self.parent = new_parent
            self.parent.add_face(self)
        
    def change_parent(self, new_parent):
        """Connects vertex to another frame. If vertex was already attached 
           to a frame, it undoes this first.
           
           TODO: Add points to parent too.
           """
        # Remove old parent first (if applicable):
        if self.parent != None:
            self.remove_parent()
        
        # Connect new parent:
        self.parent = new_parent
        self.parent.add_face(self)
    
    def check_coplanarity(self):
        """Check if all points are coplanar, by investigating determinant
               of [p12, p13, p14]. Points are coplanar only if it is zero."""
        plane_matrix = np.vstack((self.p2.xyz() - self.p1.xyz(),
                                  self.p3.xyz() - self.p1.xyz(),
                                  self.p4.xyz() - self.p1.xyz()))
        if round(np.linalg.det(plane_matrix).transpose(), 8) == 0:
            # Rounding here happens to prevent floating point errors from
            # interfering with the function. 8 decimals chosen arbitrarily.
            return True
        else:
            return False
    
    def calc_area(self):
        """Calculate area of the face, by computing the diagonals."""
        diag13 = self.p3.xyz() - self.p1.xyz()
        diag24 = self.p4.xyz() - self.p2.xyz()
        cpdiag = np.cross(diag13, diag24)
        facearea = 0.5 * np.sqrt(np.einsum('...i,...i', cpdiag, cpdiag))
        return facearea
    
    def vertices(self):
        return [self.p1, self.p2, self.p3, self.p4]

    def readout(self, dec=4):
        """Print the vertex coordinates."""
        for vertex in [self.p1, self.p2, self.p3, self.p4]:
            print(vertex)
            
    def plotlist(self):
        """Return three lists with the x, y, and z-components of all four
           vertices in the face."""
        xlist = [self.p1.global_coordinates()[0], 
                 self.p2.global_coordinates()[0],
                 self.p3.global_coordinates()[0],
                 self.p4.global_coordinates()[0]]
        ylist = [self.p1.global_coordinates()[1], 
                 self.p2.global_coordinates()[1],
                 self.p3.global_coordinates()[1],
                 self.p4.global_coordinates()[1]]
        zlist = [self.p1.global_coordinates()[2], 
                 self.p2.global_coordinates()[2],
                 self.p3.global_coordinates()[2],
                 self.p4.global_coordinates()[2]]                
                 
        # ylist = [self.p1.y, self.p2.y, self.p3.y, self.p4.y]
        # zlist = [self.p1.z, self.p2.z, self.p3.z, self.p4.z]
        return xlist, ylist, zlist

    def plotlist2(self):
        """Return a list of lists with the xyz coordinates of each vertex."""
        grid1 = np.array([self.p1.xyz_global(), self.p2.xyz_global(),
                          self.p3.xyz_global(), self.p4.xyz_global()])
        return grid1.tolist()
    
    # def project(self, plane='xy'):
    #     """Project a copy of the frame onto a plane that is spanned by two
    #        axes. The projection is orthographic.
    #        TODO: Generalize this function for arbitrary projection frame. """
    #     pj = deepcopy(self)
    #     for vertex in [pj.p1, pj.p2, pj.p3, pj.p4]:
    #         vertex.remove_parent()
    #         if plane == 'xy':
    #             vertex.z = 0
    #         elif plane == 'xz':
    #             vertex.y = 0
    #         elif plane == 'yz':
    #             vertex.x = 0
    #         else:
    #             raise ValueError("No valid projection plane given to "
    #                              "projection method! "
    #                              "Valid options: 'xy', 'xz', 'yz'")
    #     # Reinitialize the face with the projected vertices.
    #     pj = Face(pj.p1, pj.p2, pj.p3, pj.p4, pj.parent)
    #     pj.remove_parent()
    #     return pj
            
    def find_centroid(self):
        """Find the vertex centroid of the face."""
        xc = 0.25 * (self.p1.x + self.p2.x + self.p3.x + self.p4.x)
        yc = 0.25 * (self.p1.y + self.p2.y + self.p3.y + self.p4.y)
        zc = 0.25 * (self.p1.z + self.p2.z + self.p3.z + self.p4.z)
        return np.array([xc, yc, zc])

    def make_centroid(self):
        """Find the vertex centroid of the face, and turn it into a Vertex."""
        (xc, yc, zc) = self.find_centroid()
        return Vertex(xc, yc, zc, parent=self.parent)
    
    def find_perpendicular(self):
        """Find a vector perpendicular to a face (direction is ambiguous).
           TODO: Generalize the direction of the perpendicular vector.
           Currently the direction depends on how the face is defined. """
        # Find perpendicular vector to plane spanned by p12, p14.
        # Then normalize it to a unit vector
        p12 = (self.p2.xyz() - self.p1.xyz())
        p14 = (self.p4.xyz() - self.p1.xyz())
        perpendicular = np.cross(p12, p14)
        return perpendicular / np.linalg.norm(perpendicular)
    


class Geometry:
    """In this code, a Geometry is defined as a loose collection of Faces.
       The constructor merely creates an empty Geometry object, and can
       subsequently be filled with Faces using the add_face() and 
       add_faces methods.
       Some methods of this class interact with the Face objects, whilst
       others interact instead with the individual Vertices out of which
       the Faces are made, to avoid e.g. applying transformations more than
       once to the same Vertex object."""

    def __init__(self, parent=None):
        
        self.faces = []
        self.parent = parent
        self.connect_parent(parent)
    
    def remove_parent(self):
        """If geometry has a parent frame, removes it as a parent, and 
           ensures it is removed from the internal list in the frame.
           
           TODO: Garbage removal.
           """
        if self.parent != None:
            self.parent.remove_geometry(self)
            self.parent = None

    def connect_parent(self, new_parent):
        """Connects geometry to a frame.
        
           TODO: Garbage removal.
           """
        # Connect new parent:
        self.parent = new_parent
        self.parent.add_geometry(self)
        
    def change_parent(self, new_parent):
        """Connects geometry to another frame. If vertex was already attached 
           to a frame, it undoes this first.
           
           TODO: Garbage removal
           """
        # Remove old parent first (if applicable):
        if self.parent != None:
            self.remove_parent()
        
        # Connect new parent:
        self.parent = new_parent
        self.parent.add_geometry(self)
    
    def add_face(self, face: Face):
        """Add a singular face to the geometry
        
        TODO: Cascade management. """
        self.faces.append(face)

    def add_faces(self, faces: list):
        """Add a list of faces to the geometry
        
        TODO: Cascade management. """
        for face in faces:
            self.faces.append(face)
    
    def vertices(self):
        """This method loops through all the Face objects in the geometry,
           lists the vertices making up each face, and removes the duplicates
           in the list. It then returns a list of unique vertices."""
        vertices = []
        for face in self.faces:
            face_vertices = face.vertices()
            for vertex in face_vertices:
                if vertex not in vertices:
                    vertices.append(vertex)
        return vertices
    
    def find_cuboid_centroid(self):
        """Locate the centroid of a cuboid geometry.
           Cuboids have eight corners, which this method checks first.
           """
        vertices = self.vertices()
        if len(vertices) != 8:
            raise ValueError("Tried to use method find_cuboid_centroid() "
                             "with a geometry containing {} vertices. "
                             "All cuboids must have exactly 8 vertices."
                             "".format(len(vertices)))
        xc = 0
        yc = 0
        zc = 0
        for vertex in vertices:
            xc += vertex.x/len(vertices)
            yc += vertex.y/len(vertices)
            zc += vertex.z/len(vertices)
        return np.array([xc, yc, zc])
    
    def make_cuboid_centroid(self):
        """Locates the centroid of a cuboid geometry and creates a
        Vertex at this point."""
        (xc, yc, zc) = self.find_cuboid_centroid()
        return Vertex(xc, yc, zc, parent=self.parent)

    def translate(self, dx=0., dy=0., dz=0.):
        """Translate all the vertices in the geometry."""
        vertices = self.vertices()
        for vertex in vertices:
            vertex.translate(dx=dx, dy=dy, dz=dz)
            
    def rotate(self, a, b, c, cor, seq='321'):
        """Rotate all the vertices in the geometry around a point in space.
           a is Euler angle of rotation around x, etc...
               expressed in radians
           seq is a string with the rotation sequence, e.g. '321' for:
               Rz(c).Ry(b).Rx(a).vertex
           cor is the centre of rotation, which MUST be specified as
               a vertex.
           """
        vertices = self.vertices()
        for vertex in vertices:
            vertex.rotate(a, b, c, cor=cor, seq=seq)
            
    def area(self):
        """Total area of all faces in the geometry.
        Note that faces have a single side, not two sides. """
        A = 0
        for face in self.faces:
            A += face.area
        return A
    
    # def find_fcentroids(self, returnarray=False):
    #     """Returns a list of the centroid coordinates of all faces currently 
    #        in the geometry. Technically, these are vertex centroids.
    #        Optionally, centroids can be returned as a Numpy array by 
    #        setting returnarray to True."""
    #     fcentroids = []
    #     for i in range(len(self.faces)):
    #         if returnarray == True:
    #             fcentroids.append(self.faces[i].find_centroid().xyz())
    #         else:
    #             fcentroids.append(self.faces[i].find_centroid())
    #     return fcentroids

    # def find_perpendiculars(self, scale=1.):
    #     """Returns a list with a vector perpendicular to each face in the
    #        geometry. If scale=1, these vectors will be unit vectors.
            
    #        TODO: Guarantee that all perpendiculars point away from the 
    #        centre of the geometry, so that the output is not a mishmash of
    #        inward and outward pointing vertices. Can achieve this by taking
    #        inner product of:
    #            - vector pointing from geometry centroid to Face centroid.
    #            - computed perpendicular vector for given face.
    #        If the result is positive, the perpendicular is pointing away from 
    #        geometry, and if result is negative, flip the perpendicular.
    #        Downside: Does not work well for concave geometries whose centroid
    #        is outside the geometry boundaries.
    #        """
    #     perpendiculars = []
    #     for i in range(len(self.faces)):
    #         perp = self.faces[i].find_perpendicular()
    #         perpendiculars.append(perp * scale)

    #     # Code that ensures arrows point outwards (TODO):
    #     # ...

    #     return perpendiculars

    # def perpendiculars_plotlist(self):
    #     """Returns a list of lists with coordinates used for plotting.
    #        Each item in the list is a list of six coordinates:
    #        - The first three coordinates indicate the xyz of the tail point.
    #        - The second three coordinates are point a vector from this tail
    #          point.
           
    #        TODO: Pass scaling as argument.
    #        """
    #     c = self.find_fcentroids(returnarray=True)
    #     p = self.find_perpendiculars(scale=0.05)  # Vary 0.05 if needed.
    #     if len(c) != len(p):
    #         raise ValueError("Number of centroids and perpendiculars is "
    #                          "somehow unequal! This should not happen.")
    #     plotarray = np.hstack([np.vstack(c), np.vstack(p)]).transpose()
    #     return plotarray.tolist()
        
    
class Vector:
    pass

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
        
    def add_face(self, face: Face):
        """TODO: merge into one add_child method."""
        self.faces.append(face)

    def remove_face(self, face: Face):
        """TODO: merge into one add_child method."""
        self.faces.remove(face)

    def add_geometry(self, geometry: Geometry):
        """TODO: merge into one add_child method."""
        self.geometries.append(geometry)

    def remove_geometry(self, geometry: Geometry):
        """TODO: merge into one add_child method."""
        self.geometries.remove(geometry)
    
    def add_vector(self, vector: Vector):
        """TODO: merge into one add_child method."""
        self.vectors.append(vector)

    def remove_vector(self, vector: Vector):
        """TODO: merge into one add_child method."""
        self.vectors.remove(vector)
    
    def origin(self):
        return np.array([self.x, self.y, self.z])
    
    def recalculate_dcm(self):
        gx = np.array([1,0,0])
        gy = np.array([0,1,0])
        gz = np.array([0,0,1])
        self.dcm = np.array([
            [np.dot(gx, self.xdir), 
             np.dot(gx, self.ydir), 
             np.dot(gx, self.zdir)],
            
            [np.dot(gy, self.xdir), 
             np.dot(gy, self.ydir), 
             np.dot(gy, self.zdir)],
            
            [np.dot(gz, self.xdir), 
             np.dot(gz, self.ydir), 
             np.dot(gz, self.zdir)]])
        
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

def plot_face(axes, face: Face, fill=1, alpha=0.9,
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
    
    for face in frame.faces:
        plot_face(axes, face)

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
            frame1.translate(2/steps,4/steps,4/steps)
        else:
            frame1.translate(2/steps,-4/steps,-4/steps)
        frame1.rotate(0,-2*np.pi/steps,-2*np.pi/steps)
        # frame1.rotate(0,0,0)
        
        # vertex1.rotate(2*np.pi/steps,0,0, cor=vertex2)
        
        # Setting up the axes object
        ax.clear()
        
        ax.set_title("Wireframe visualization. Frame: {}".format(str(i)))
        
        ax.set_xlim(0, 2.5)
        ax.set_ylim(0, 2.5)
        ax.set_zlim(0, 2)
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

    ani = animation.FuncAnimation(fig, update, np.arange(1,steps), 
                                  interval = 75, repeat = False)

    plt.show()
