# -*- coding: utf-8 -*-
"""
Created on Thu Mar 4 15:02:11 2021    

@author: Johan Monster

Version: 1.0

"""

from copy import deepcopy
import numpy as np
from numpy import sin, cos
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d as mp3d


def r2d(a):
    """Convert radians to degrees."""
    return a * 180 / np.pi


def d2r(a):
    """Convert degrees to radians."""
    return a * np.pi / 180


""" CHANGE ROTATION OF BODY HERE: """
# Note: choosing bdeg=90 will result in gimbal lock. 
adeg = 30  # Rotation around x axis in degrees
bdeg = 45  # Rotation around y axis in degrees
cdeg = 0   # Rotation around z axis in degrees

a = d2r(adeg)
b = d2r(bdeg)
c = d2r(cdeg)


class Vertex:
    """Custom implentation of a 3D point expressed in cartesians."""

    def __init__(self, x=0., y=0., z=0.):
        self.x = x
        self.y = y
        self.z = z

    def translate(self, dx=0., dy=0., dz=0.):
        self.x += dx
        self.y += dy
        self.z += dz

    def rotate_origin(self, a, b, c, seq='321'):
        """ Rotates the vertex around the origin.
            a is Euler angle of rotation around x, etc...
                expressed in radians
            seq is a string with the rotation sequence, e.g. '321' for:
                Rz(c).Ry(b).Rx(a).vertex
        """

        Rx = np.array([[1, 0, 0],
                       [0, cos(a), -sin(a)],
                       [0, sin(a), cos(a)]])

        Ry = np.array([[cos(b), 0, sin(b)],
                       [0, 1, 0],
                       [-sin(b), 0, cos(b)]])

        Rz = np.array([[cos(c), -sin(c), 0],
                       [sin(c), cos(c), 0],
                       [0, 0, 1]])

        for axis in reversed(seq):
            if axis == '1':
                newcoor = np.dot(Rx, np.array([self.x, self.y, self.z]))
            if axis == '2':
                newcoor = np.dot(Ry, np.array([self.x, self.y, self.z]))
            if axis == '3':
                newcoor = np.dot(Rz, np.array([self.x, self.y, self.z]))
            self.__init__(newcoor[0], newcoor[1], newcoor[2])

    def rotate(self, a, b, c, cor, seq='321'):
        """ Rotates the vertex around a point in space.
            a is Euler angle of rotation around x, etc...
                expressed in radians
            seq is a string with the rotation sequence, e.g. '321' for:
                Rz(c).Ry(b).Rx(a).vertex
            cor is the centre of rotation, which MUST be specified as
                a vertex.
        """

        Rx = np.array([[1, 0, 0],
                       [0, cos(a), -sin(a)],
                       [0, sin(a), cos(a)]])

        Ry = np.array([[cos(b), 0, sin(b)],
                       [0, 1, 0],
                       [-sin(b), 0, cos(b)]])

        Rz = np.array([[cos(c), -sin(c), 0],
                       [sin(c), cos(c), 0],
                       [0, 0, 1]])

        for axis in reversed(seq):
            self.translate(-1 * cor.x, -1 * cor.y, -1 * cor.z)

            if axis == '1':
                newcoor = np.dot(Rx, np.array([self.x, self.y, self.z]))
            if axis == '2':
                newcoor = np.dot(Ry, np.array([self.x, self.y, self.z]))
            if axis == '3':
                newcoor = np.dot(Rz, np.array([self.x, self.y, self.z]))
            self.__init__(newcoor[0], newcoor[1], newcoor[2])

            self.translate(1 * cor.x, 1 * cor.y, 1 * cor.z)

    def array(self):
        """Returns vertex coordinates as a Numpy array."""
        return np.array([self.x, self.y, self.z])

    def display(self, dec=4):
        """Print the vertex coordinates."""
        if dec == -1:  # Set dec to -1 to disable rounding
            print([self.x, self.y, self.z])
        else:
            print([round(self.x, dec),
                   round(self.y, dec),
                   round(self.z, dec)])


class Line:
    """Unused line class."""

    def __init__(self, point1=None, point2=None):
        self.point1 = point1 or Vertex()
        self.point2 = point2 or Vertex()

        self.i = self.point1.x - self.point2.x
        self.j = self.point1.y - self.point2.y
        self.k = self.point1.z - self.point2.z

        self.length = np.sqrt(self.i * self.i +
                              self.j * self.j +
                              self.k * self.k)

    def magnitude(self):
        return self.length


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
    """

    def __init__(self, p1: Vertex, p2: Vertex, p3: Vertex, p4: Vertex):
        self.p1 = p1
        self.p2 = p2
        self.p3 = p3
        self.p4 = p4
        if self.check_coplanarity() == False:
            raise ValueError("Error! Set of four points is not coplanar!")
        self.area = self.area()

    def check_coplanarity(self):
        """Check if all points are coplanar, by investigating determinant
               of [p12, p13, p14]. Points are coplanar only if it is zero."""
        plane_matrix = np.vstack((self.p2.array() - self.p1.array(),
                                  self.p3.array() - self.p1.array(),
                                  self.p4.array() - self.p1.array()
                                  )).transpose()
        if round(np.linalg.det(plane_matrix), 8) == 0:
            # Rounding here happens to prevent floating point errors from
            # interfering with the function. 8 decimals chosen arbitrarily.
            return True
        else:
            return False

    def area(self):
        """Calculate area of the face, by computing the diagonals."""
        diag13 = self.p3.array() - self.p1.array()
        diag24 = self.p4.array() - self.p2.array()
        cpdiag = np.cross(diag13, diag24)
        facearea = 0.5 * np.sqrt(np.einsum('...i,...i', cpdiag, cpdiag))
        return facearea

    def project(self, plane='xy'):
        """Project a copy of the frame onto a plane that is spanned by two
           axes. The projection is orthographic.
           TODO: Generalize this function for arbitrary projection frame. """
        pj = deepcopy(self)
        for vertex in [pj.p1, pj.p2, pj.p3, pj.p4]:
            if plane == 'xy':
                vertex.z = 0
            elif plane == 'xz':
                vertex.y = 0
            elif plane == 'yz':
                vertex.x = 0
            else:
                raise ValueError("No valid projection plane given to "
                                 "projection method! "
                                 "Valid options: 'xy', 'xz', 'yz'")
        # Reinitialize the face with the projected vertices.
        pj = Face(pj.p1, pj.p2, pj.p3, pj.p4)
        return pj

    def find_centroid(self):
        """Find the vertex centroid of the face."""
        xc = 0.25 * (self.p1.x + self.p2.x + self.p3.x + self.p4.x)
        yc = 0.25 * (self.p1.y + self.p2.y + self.p3.y + self.p4.y)
        zc = 0.25 * (self.p1.z + self.p2.z + self.p3.z + self.p4.z)
        return Vertex(xc, yc, zc)

    def find_perpendicular(self):
        """Find a vector perpendicular to a face (direction is ambiguous).
           TODO: Generalize the direction of the perpendicular vector.
           Currently the direction depends on how the face is defined. """
        # Find perpendicular vector to plane spanned by p12, p14.
        # Then normalize it to a unit vector
        p12 = (self.p2.array() - self.p1.array())
        p14 = (self.p4.array() - self.p1.array())
        perpendicular = np.cross(p12, p14)
        return perpendicular / np.linalg.norm(perpendicular)

    def display(self):
        """Print vertices of face to console."""
        for point in [self.p1, self.p2, self.p3, self.p4]:
            point.display()

    def plotlist(self):
        """Return three lists with the x, y, and z-components of all four
           vertices in the face."""
        xlist = [self.p1.x, self.p2.x, self.p3.x, self.p4.x]
        ylist = [self.p1.y, self.p2.y, self.p3.y, self.p4.y]
        zlist = [self.p1.z, self.p2.z, self.p3.z, self.p4.z]
        return xlist, ylist, zlist

    def plotlist2(self):
        """Return a list of lists with the xyz coordinates of each vertex."""
        grid1 = np.array([self.p1.array(), self.p2.array(),
                          self.p3.array(), self.p4.array()])
        return grid1.tolist()


class Geometry():
    """In this code, a Geometry is defined as a loose collection of Faces.
       The constructor merely creates an empty Geometry object, and can
       subsequently be filled with Faces using the add_face() and 
       add_faces methods.
       Some methods of this class interact with the Face objects, whilst
       others interact instead with the individual Vertices out of which
       the Faces are made, to avoid e.g. applying transformations more than
       once to the same Vertex object."""

    def __init__(self):
        """Creates an empty list that can be filled with faces."""
        self.faces = []

    def add_face(self, face):
        """Add a singular face to the geometry"""
        self.faces.append(face)

    def add_faces(self, faces):
        """Add a list of faces to the geometry"""
        for face in faces:
            self.faces.append(face)

    def find_centroid(self):
        """Locate the centroid of the geometry.
        
           TODO: Write the function."""
        pass

    def find_centroids(self, returnarray=False):
        """Returns a list of the centroid coordinates of all faces currently 
           in the geometry. Technically, these are vertex centroids.
           Optionally, centroids can be returned as a Numpy array by 
           setting returnarray to True."""
        centroids = []
        for i in range(len(self.faces)):
            if returnarray == True:
                centroids.append(self.faces[i].find_centroid().array())
            else:
                centroids.append(self.faces[i].find_centroid())
        return centroids

    def find_perpendiculars(self, scale=1.):
        """Returns a list with a vector perpendicular to each face in the
           geometry. If scale=1, these vectors will be unit vectors.
            
           TODO: Guarantee that all perpendiculars point away from the 
           centre of the geometry, so that the output is not a mishmash of
           inward and outward pointing vertices. Can achieve this by taking
           inner product of:
               - vector pointing from geometry centroid to Face centroid.
               - computed perpendicular vector for given face.
           If the result is positive, the perpendicular is pointing away from 
           geometry, and if result is negative, flip the perpendicular.
           Downside: Does not work well for concave geometries whose centroid
           is outside the geometry boundaries.
           """
        perpendiculars = []
        for i in range(len(self.faces)):
            perp = self.faces[i].find_perpendicular()
            perpendiculars.append(perp * scale)

        # Code that ensures arrows point outwards (TODO):
        # ...

        return perpendiculars

    def perpendiculars_plotlist(self):
        """Returns a list of lists with coordinates used for plotting.
           Each item in the list is a list of six coordinates:
           - The first three coordinates indicate the xyz of the tail point.
           - The second three coordinates are point a vector from this tail
             point.
           
           TODO: Pass scaling as argument.
           """
        c = self.find_centroids(returnarray=True)
        p = self.find_perpendiculars(scale=0.05)  # Vary 0.05 if needed.
        if len(c) != len(p):
            raise ValueError("Number of centroids and perpendiculars is "
                             "somehow unequal! This should not happen.")
        plotarray = np.hstack([np.vstack(c), np.vstack(p)]).transpose()
        return plotarray.tolist()

    def area(self):
        """Total area of all faces in the geometry.
        Note that faces have a single side, not two sides. """
        A = 0
        for face in self.faces:
            A += face.area
        return A

    def illumination_vector(self, plane='xy'):
        """Generate an illumination vector based on a given plane.
           This vector essentially "simulates" where the light is coming from.
           Currently limited to directions parallel to global XYZ system.
           """
        if plane == 'xy':
            iv = np.array([0, 0, 1])
        elif plane == '-xy':
            iv = np.array([0, 0, -1])
        elif plane == 'xz':
            iv = np.array([0, 1, 0])
        elif plane == '-xz':
            iv = np.array([0, -1, 0])
        elif plane == 'yz':
            iv = np.array([1, 0, 0])
        elif plane == '-yz':
            iv = np.array([-1, 0, 0])
        else:
            raise ValueError("Invalid plane given. Options: 'xy', '-xy', "
                             "'xz', '-xz', 'yz', '-yz'")
        return iv

    def illuminated_faces(self, plane='xy'):
        """Computes which faces of the geometry are illuminated given
        a certain illumination vector. It returns a vector representing
        each face in the geometry. If the value is 0, the face is not
        illuminated. If it is non-zero, it represents the angle in radians
        that the respective face makes with the illumination vector."""

        # Fetch direction of illumination and face perpendiculars.
        # Both MUST be unit vectors (and should be).
        iv = self.illumination_vector(plane=plane)
        perps = self.find_perpendiculars(scale=1)

        # Pre-allocate illuminated_faces vector
        illuminated_faces = np.zeros(len(self.faces))
        for i in range(len(perps)):
            # Calculate dot product between iv and face i
            ill_status = np.dot(iv, perps[i])
            # Face is illuminated only if ill_status < 0:
            if ill_status < 0:
                illuminated_faces[i] = np.arccos(ill_status) - np.pi / 2
            else:
                pass
        return illuminated_faces

    def illuminated_area(self, plane='xy'):
        """Computes the projected area of the illuminated area.
        Requires that projection plane and illumination vector are 
        exactly perpendicular, for example:
            if projection plane is xy, then iv must be [0, 0, +1].
        Units of the area are in terms of the same units as the vectors."""

        # Fetch projected shape and illuminated faces information:
        ill_faces = self.illuminated_faces(plane=plane)
        projection = self.project_geometry(plane=plane)

        # Loop over illuminated faces and sum area:
        illuminated_area = 0
        for i in range(len(ill_faces)):
            # If face is illuminated, do something:
            if ill_faces[i] != 0:
                # Add the area of the current projected face:
                illuminated_area += projection.faces[i].area
        return illuminated_area

    def project_geometry(self, plane='xy'):
        """Project a whole geometry by individually projecting each Face
           object in self.faces.
           """
        projected_faces = []
        for i in range(len(self.faces)):
            projected_faces.append(self.faces[i].project(plane=plane))

        projected_geometry = Geometry()
        projected_geometry.add_faces(projected_faces)
        return projected_geometry


# %% Define body

# # Old junk:
# point1 = Vertex(0, 0, 0)
# point2 = Vertex(1, 0, 0)
# point3 = Vertex(1, 0, 2)
# point4 = Vertex(0, 0, 2)

# point1 = Vertex(0, 0, 0)
# point2 = Vertex(1, 0, 1)
# point3 = Vertex(2, 2, 4)
# point4 = Vertex(0, 1, 1)


# face1 = Face(point1, point2, point3, point4)
# face1.display()
# print("A_face1 = ", face1.area)

# face1xy = face1.project('xy')
# face1xy.display()
# print("A_pj1_xy =", face1xy.area)

# face1xz = face1.project('xz')
# face1yz = face1.project('yz')

# # Some more junk:
# tp1 = Vertex(0.0, 0.0, 0.0)
# tp2 = Vertex(0.1, 0.0, 0.0)
# tp3 = Vertex(0.1, 0.1, 0.0)
# tp4 = Vertex(0.0, 0.1, 0.0)

# tf1 = Face(tp1, tp2, tp3, tp4)

# geometryt1 = Geometry()
# geometryt1.add_face(tf1)


"""Basic Cubesat model:"""
p1 = Vertex(0.0, 0.0, 0.0)
p2 = Vertex(0.0, 0.0, 0.2)
p3 = Vertex(0.1, 0.0, 0.2)
p4 = Vertex(0.1, 0.0, 0.0)
p5 = Vertex(0.0, 0.1, 0.0)
p6 = Vertex(0.0, 0.1, 0.2)
p7 = Vertex(0.1, 0.1, 0.2)
p8 = Vertex(0.1, 0.1, 0.0)

pointcollection = [p1, p2, p3, p4, p5, p6, p7, p8]

# Manually calculate centroid of Cubesat geometry
COR = Vertex(0.05, 0.05, 0.1)
COR.translate(0.1, 0.1, 0.1)

# Transform individual vertices. This is a really lame way of doing this.
# TODO: move the circuitry for transforming vertices inside class structure.
# Problem: How to prevent points from being transformed more than once?
for vertex in pointcollection:
    vertex.translate(0.1, 0.1, 0.1)
    vertex.rotate(a, b, c, COR)

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

# Create a projection object that can be plotted later.
geometry1xy = geometry1.project_geometry('xy')

# %% Plotting code

# Check value for illuminated area. Because Cubesat is a cuboid, illuminated
# area and the unilluminated area should be the same UNLESS only one side is 
# facing the light. Used for debugging.
A_shadow_check = round(geometry1xy.area() / 2, 4)

# Properly compute the illuminated area using Geometry method:
A_shadow = round(geometry1.illuminated_area(plane='xy'), 4)

print('Coarse method: A_xy =', A_shadow_check, 'm^2')
print('Projected A_xy =', A_shadow, "m^2")

# Toggle plotting functionality:
if True:

    fig = plt.figure(figsize=(10, 7))
    ax = mp3d.Axes3D(fig)

    ax.set_title("Wireframe visualization, A_xy = {} m^2".format(A_shadow))

    """ TO CHANGE THE DEFAULT CAMERA VIEW, CHANGE THESE: """
    ax.view_init(elev=20, azim=-90)

    # # Forcing aspect ratio for Axes3D object is still fucked up in
    # # matplotlib source, even after many years. :(
    # ax.set_aspect('equal')
    # ax.set_box_aspect((1,1,1))

    ax.set_xlim(0, 0.4)
    ax.set_ylim(0, 0.4)
    ax.set_zlim(0, 0.3)

    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')


    def plot_xyz_tripod(axes, alpha=1, scaling=1.):
        """Plots an rgb xyz tripod at the global origin.
            - alpha changes the opacity of the arrows. (default: 1)
            - scaling changes the length of the arrows (default: 1)
            """
        ax.quiver(0, 0, 0, scaling, 0, 0,
                  arrow_length_ratio=0.15, color='red', alpha=alpha)
        ax.quiver(0, 0, 0, 0, scaling, 0,
                  arrow_length_ratio=0.15, color='green', alpha=alpha)
        ax.quiver(0, 0, 0, 0, 0, scaling,
                  arrow_length_ratio=0.15, color='blue', alpha=alpha)


    # Plotting the XYZ tripod
    plot_xyz_tripod(ax, scaling=0.25)

    # Plotting the centroid of the geometry (manually)
    ax.scatter(COR.x, COR.y, COR.z, c='cyan', s=20, alpha=0.5)


    def plot_face(axes, face, fill=0, alpha=0.2,
                  linecolour='black', linealpha=1,
                  illumination=False, ill_value=0):
        """ Plots an individual Face object. Check source code for kwargs.
        """
        (xlist, ylist, zlist) = face.plotlist()

        # Plot individual vertices:
        axes.scatter(xlist, ylist, zlist, c=linecolour, s=10)

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
            # If illumination functionality is turned on, face is coloured
            # based on its illumination value.
            if illumination == True:
                # No illumination: rgb = [0.1, 0.1, 0.1] (dark gray)
                # Total illumination: [1, 1, 1] (white)
                value = round(0.1 + 0.9 / 1.6 * ill_value, 2)
                fs.set_facecolor((value, value, value, alpha))
            else:
                fs.set_facecolor((0.1, 0.1, 0.1, alpha))
            axes.add_collection3d(fs)


    def plot_geometry(axes, geometry, fill=0, alpha=0.2,
                      linecolour='black', linealpha=1,
                      illumination=False):
        """Plot geometry by individually plotting its faces. If illumination 
           is turned on, uses illuminated_faces() method to fetch which faces 
           are illuminated, and passes this to the plot_face() function.
           
           TODO: Allow passing of 'plane'
           """
        if illumination == True:
            ill_faces = geometry.illuminated_faces(plane='xy')
            for idx, face in enumerate(geometry.faces):
                plot_face(axes, face, fill=fill, alpha=alpha,
                          linecolour=linecolour, linealpha=linealpha,
                          illumination=True, ill_value=ill_faces[idx])
        else:
            for idx, face in enumerate(geometry.faces):
                plot_face(axes, face, fill=fill, alpha=alpha,
                          linecolour=linecolour, linealpha=linealpha)


    # plot_face(ax, geometry1.faces[0], colour='black')

    # Plot Cubesat model
    plot_geometry(ax, geometry1, linecolour='black', fill=1, alpha=1,
                  illumination=True)
    # Plot projection of Cubesat
    plot_geometry(ax, geometry1xy, linecolour='orange')


    def plot_geometry_perpendiculars(axes, geometry, colour='gray', alpha=0.5):
        """Plots the normals/perpendiculars of each face in a geometry, and 
           displays them as little gray arrows.
           TODO: Implement scaling."""
        # First plot centroids of each plane:
        xc = []
        yc = []
        zc = []

        centroids = geometry.find_centroids()

        for vertex in centroids:
            xc.append(vertex.x)
            yc.append(vertex.y)
            zc.append(vertex.z)
        axes.scatter(xc, yc, zc, c=colour, s=5, alpha=alpha)

        # Then, attach plane-perpendicular arrows to the centroids:
        xyzuvw = geometry.perpendiculars_plotlist()
        ax.quiver(xyzuvw[0][:], xyzuvw[1][:], xyzuvw[2][:],
                  xyzuvw[3][:], xyzuvw[4][:], xyzuvw[5][:],
                  color=colour, alpha=alpha)


    # Highlight one of the faces in bright yellow (useful for debugging)
    # plot_face(ax, geometry1.faces[5], colour='yellow')    

    # Plotting the perpendiculars
    plot_geometry_perpendiculars(ax, geometry1)

    fig.show()
