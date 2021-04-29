# -*- coding: utf-8 -*-
"""
Created on Thu April 27 15:02:11 2021    

@author: Johan Monster

Cube Projector - Geometry class
Version: 1.2

"""

import numpy as np
from cp_vertex import Vertex
from cp_face import Face

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
        self.frame = np.array([[1, 0, 0],
                               [0, 1, 0],
                               [0, 0, 1]])

    def add_face(self, face: Face):
        """Add a singular face to the geometry"""
        self.faces.append(face)

    def add_faces(self, faces: list):
        """Add a list of faces to the geometry"""
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

    def find_centroid(self):
        """Locate the centroid of the geometry.
        
           TODO: Write the function.
           Almost impossible for arbitrary geometry.
           """
        pass

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
        return Vertex(xc, yc, zc)

    def translate(self, dx=0., dy=0., dz=0.):
        """Translate all the vertices in the geometry."""
        vertices = self.vertices()
        for vertex in vertices:
            vertex.translate(dx=dx, dy=dy, dz=dz)
    
    def rotate_origin(self, a, b, c, seq='321'):
        """Rotate all the vertices in the geometry around the origin.
           a is Euler angle of rotation around x, etc...
               expressed in radians
           seq is a string with the rotation sequence, e.g. '321' for:
               Rz(c).Ry(b).Rx(a).vertex
           """
        vertices = self.vertices()
        for vertex in vertices:
            vertex.rotate_origin(a=a, b=b, c=c, seq=seq)
    
    # TODO: Find way to unify rotation functions, as this is a bit messy.
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
            vertex.rotate(a=a, b=b, c=c, cor=cor, seq=seq)
        
    def rotate_cuboid_centroid(self, a, b, c, seq='321'):
        """Rotate all the vertices in the geometry around its cuboid centroid.
           WARNING: Only works if the geometry is a valid cuboid.
           
           a is Euler angle of rotation around x, etc...
               expressed in radians
           seq is a string with the rotation sequence, e.g. '321' for:
               Rz(c).Ry(b).Rx(a).vertex
           cor is the centre of rotation, which MUST be specified as
               a vertex.
           """
        centroid = self.find_cuboid_centroid()
        vertices = self.vertices()
        for vertex in vertices:
            vertex.rotate(a=a, b=b, c=c, cor=centroid, seq=seq)

    def find_fcentroids(self, returnarray=False):
        """Returns a list of the centroid coordinates of all faces currently 
           in the geometry. Technically, these are vertex centroids.
           Optionally, centroids can be returned as a Numpy array by 
           setting returnarray to True."""
        fcentroids = []
        for i in range(len(self.faces)):
            if returnarray == True:
                fcentroids.append(self.faces[i].find_centroid().array())
            else:
                fcentroids.append(self.faces[i].find_centroid())
        return fcentroids

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
        c = self.find_fcentroids(returnarray=True)
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
    
    def illumination_vector_new(self, vector):
        """Converts a given vector into a unit vector for use as an
           illumination vector.
           """
        return vector / np.linalg.norm(vector)
 
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

    def illuminated_faces_new(self, illumination_vector):
        """Computes which faces of the geometry are illuminated given
        a certain illumination vector. It returns a vector representing
        each face in the geometry. If the value is 0, the face is not
        illuminated. If it is non-zero, it represents the angle in radians
        that the respective face makes with the illumination vector."""

        # Fetch direction of illumination and face perpendiculars.
        # Both MUST be unit vectors (and should be).
        iv = self.illumination_vector(illumination_vector)
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
        Units of the area are in terms of the same units as the vectors.
        TODO: Possibly merge with project_illuminated_area() method."""

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

    def project_illuminated_faces(self, plane='xy'):
        """Identical to the project_geometry() method, except it only
           projects faces that are illuminated.
           """
        projected_faces = []
        ill_faces = self.illuminated_faces(plane=plane)
        
        for i in range(len(self.faces)):
            # Filter illuminated faces
            if ill_faces[i] != 0:
                projected_faces.append(self.faces[i].project(plane=plane))

        projected_geometry = Geometry()
        projected_geometry.add_faces(projected_faces)
        return projected_geometry

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
    
    def readout(self):
        """Print the vertices in in the geometry to console."""
        for vertex in self.vertices():
            vertex.readout()
