# -*- coding: utf-8 -*-
"""
Created on Sat May 01 5:43:09 2021   

@author: Johan Monster

Cube Projector - Frame class
Version: 2.0

"""

import numpy as np
from numpy import sin, cos

from cp_vertex import Vertex
from cp_face import Face
from cp_geometry import Geometry
from cp_vector import Vector

class Frame:
    """Custom implentation of a 3D point expressed in cartesians.
    
       TODO: Clean up projection/illumination functions in a more logical
               manner.
       TODO: Add support for parenttype being "frame", in which case frames
               can be daisy-chained right up to the global frame.
    """

    def __init__(self, xyz: list = [0., 0., 0.], parenttype="global"):
        self.x = xyz[0]
        self.y = xyz[1]
        self.z = xyz[2]
        self.xdir = np.array([1,0,0])
        self.ydir = np.array([0,1,0])
        self.zdir = np.array([0,0,1])
        
        self.dcm = None
        self.recalculate_dcm()
        
        self.vertices = []
        self.faces = []
        self.geometries = []
        self.vectors = []
        
        self.parenttype = "global"
        # Re-assign self.parenttype to check for validity of given parenttype
        if parenttype != "global":
            self.set_parenttype(parenttype)


    def add_vertex(self, vertex: Vertex):
        self.vertices.append(vertex)
        vertex.set_parenttype="frame"


    def remove_vertex(self, vertex: Vertex):
        self.vertices.remove(vertex)
        vertex.set_parenttype="global"
        
        
    def add_face(self, face: Face):
        self.faces.append(face)
        face.set_parenttype="frame"


    def remove_face(self, face: Face):
        self.faces.remove(face)
        face.set_parenttype="global"


    def add_geometry(self, geometry):
        if isinstance(geometry, Geometry):
            self.geometries.append(geometry)
        elif isinstance(geometry, list):
            for geo in geometry:
                if isinstance(geo, Geometry):
                    self.geometries.append(geometry)
                else:
                    raise TypeError("Frame.add_geometry() argument must be a \
                                    list of geometries!")
        else:
            raise TypeError("Invalid argument provided to \
                            Frame.add_geometry()!")
    
    
    def remove_geometry(self, geometry: Geometry):
        """TODO: merge into one add_child method."""
        self.geometries.remove(geometry)
        geometry.set_parenttype="global"
    
    
    def add_vector(self, vector: Vector):
        """TODO: merge into one add_child method."""
        self.vectors.append(vector)
        vector.set_parenttype="frame"


    def remove_vector(self, vector: Vector):
        """TODO: merge into one add_child method."""
        self.vectors.remove(vector)
        vector.set_parenttype="global"
    
    
    def origin(self):
        """Returns the origin of the frame.
        
           Returns: origin (numpy.ndarray)
        """
        return np.array([self.x, self.y, self.z])
    
    
    def recalculate_dcm(self):
        """Computes the direction cosine matrix (DCM) for conversion between
           the frame's local coordinates and the global frame.
           
           Returns: dcm (3-by-3 numpy.ndarray)
           """
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
        """Move the frame's origin by a relative amount."""
        self.x += dx 
        self.y += dy
        self.z += dz
    
    
    def rotate(self, a, b, c, cor=None, seq='321'):
        """ Rotates the frame around a point in the global frame
        .
            a, b, c is Euler angle of rotation around local x, y, z axes
                expressed in radians
                
            seq is a string with the rotation sequence, e.g. '321' for:
                Rz(c).Ry(b).Rx(a).vertex
                Default rotation sequence: '321' (zyx)
                
            cor is the centre of rotation, which can be specified as
                a Vertex, or a list or numpy.ndarray of local coordinates.
        """
        
        # ==== COR argument handling ====
        # Case: COR is given as a list/ndarray of coordinates  
        if (isinstance(cor, list) or isinstance(cor, np.ndarray)):
            if len(cor) == 3:
                cor = Vertex(list(cor))
                
        # Case: COR is given as a Vertex
        elif isinstance(cor, Vertex):
            pass
        
        # If COR is None, ignore now and handle later
        elif not cor:
            pass
            
        # All other cases
        else:
            raise TypeError("Error when rotating frame. Centre of rotation \
                            argument not provided correctly!")
        
        # ==== Define relevant rotation matrices ====
        Rx = np.array([[1,      0,       0],
                       [0, cos(a), -sin(a)],
                       [0, sin(a),  cos(a)]])
        
        Ry = np.array([[ cos(b), 0, sin(b)],
                       [      0, 1,      0],
                       [-sin(b), 0, cos(b)]])
        
        Rz = np.array([[cos(c), -sin(c), 0],
                       [sin(c),  cos(c), 0],
                       [     0,       0, 1]])
        
        # ==== Perform rotation ====
        for axis in reversed(seq):
            
            # If no Centre of Rotation (cor) is given, take frame origin.
            # If not, subtract coordinates of cor first before applying 
            # the rotation, and then reverse this afterwards.
            temp = (self.x, self.y, self.z)
            if cor:
                self.translate(-1*temp[0], -1*temp[1], -1*temp[2])
            
            if axis == '1':
                self.xdir = np.dot(Rx, self.xdir)
                self.ydir = np.dot(Rx, self.ydir)
                self.zdir = np.dot(Rx, self.zdir)
            if axis == '2':
                self.xdir = np.dot(Ry, self.xdir)
                self.ydir = np.dot(Ry, self.ydir)
                self.zdir = np.dot(Ry, self.zdir)
            if axis == '3':
                self.xdir = np.dot(Rz, self.xdir)
                self.ydir = np.dot(Rz, self.ydir)
                self.zdir = np.dot(Rz, self.zdir)
              
            # Reverse temporary translation.
            if cor:
                self.translate(1*temp[0], 1*temp[1], 1*temp[2])
        
        # Update the frame's direction cosine matrix:
        self.recalculate_dcm()
    
    
    def vertex_xyz_global(self, vertex):
        """Vertex coordinates in terms of parent frame."""
        
        # If vertex object was given:
        if isinstance(vertex, Vertex):
            xyz_local = vertex.xyz()
        elif (isinstance(vertex, list) or isinstance(vertex, np.ndarray)):
            if len(vertex) == 3:
                xyz_local = np.array(vertex)
            else:
                raise ValueError("This function can only convert 3D vertices,\
                                  but {} coordinates were given!"
                                  .format(len(vertex)))
        else:
            raise TypeError("Function vertex_xyz_global does not accept \
                            arguments of this type!")
        
        # Coordinates of local frame origin in terms of global frame
        o_local = np.array([self.x, self.y, self.z])
        
        # Vertex coordinates in terms of global frame
        xyz_global = np.dot(self.dcm, xyz_local) + o_local
        
        # return xyz_global[0], xyz_global[1], xyz_global[2]
        return xyz_global        
    
    
    def geometry_vertices(self):
        geometry_vertices = []
        for geometry in self.geometries:
            geometry_vertices.extend(geometry.vertices())
        return list(set(geometry_vertices))
    
    
    def plotlist(self, face: Face):
        """Return three lists with the x, y, and z-components of all four
           vertices in the face."""
        
        xlist = []
        ylist = []
        zlist = []
    
        for vertex in face.vertices():
            xg, yg, zg = self.vertex_xyz_global(vertex)
            xlist.append(xg)
            ylist.append(yg)
            zlist.append(zg)
        
        return xlist, ylist, zlist
    
    
    def plotlist_xyz(self, face: Face):
        """Return a list of lists with the xyz coordinates of each vertex."""
        # If object has no parent, return plot list in terms of global frame:
        
        plotlist_xyz = []
        
        for vertex in face.vertices():
            plotlist_xyz.append(list(self.vertex_xyz_global(vertex)))

        return plotlist_xyz
    
    
    def face_perpendicular(self, face: Face, scale=1.):
        """Returns the perpendicular of a specified face in terms of the
           global coordinate frame.
           
           Returns: perpendicular (3-by-1 numpy.ndarray)
           """
        p1 = self.vertex_xyz_global(face.p1)
        p2 = self.vertex_xyz_global(face.p2)
        p4 = self.vertex_xyz_global(face.p4)
        
        perpendicular = np.cross(p2-p1, p4-p1)
        return perpendicular / np.linalg.norm(perpendicular)
    
    
    def geometry_perpendiculars(self, geometry: Geometry, scale=1.):
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
        for face in geometry.faces:
            perp = self.face_perpendicular(face)
            perpendiculars.append(perp * scale)

        # Code that ensures arrows point outwards (TODO):
        # ...

        return perpendiculars
    
    
    def face_centroid(self, face: Face):
        """Find the vertex centroid of the face."""
        
        xc = 0
        yc = 0
        zc = 0
        
        for vertex in face.vertices():
            (dx, dy, dz) = self.vertex_xyz_global(vertex)
            xc += 0.25*dx
            yc += 0.25*dy
            zc += 0.25*dz
    
        return np.array([xc, yc, zc])
    
    
    def geometry_face_centroids(self, geometry: Geometry):
        """Find the vertex centroids of all the faces in specified geometry
           in terms of the global coordinate frame.
           
           Returns: np.ndarray
           """
        face_centroids = []
        for face in geometry.faces:
            face_centroids.append(self.face_centroid(face))
        return np.array(face_centroids)

    
    
    def perpendiculars_plotlist(self, geometry: Geometry, scale=.1):
        """Returns a list of lists with coordinates used for plotting.
            Each item in the list is a list of six coordinates:
            - The first three coordinates indicate the xyz of the tail point.
            - The second three coordinates are point a vector from this tail
              point.
              
            TODO: Pass scaling as argument.
            """
        c = self.geometry_face_centroids(geometry)
        p = self.geometry_perpendiculars(geometry, scale)  # Note: 0.05

        if len(c) != len(p):
            raise ValueError("Number of centroids and perpendiculars is "
                              "somehow unequal! This should not happen.")

        # Generate plotarray, which is structured as [xyzuvw1, xyzuvw2, ...]
        plotarray = np.hstack([np.vstack(c), np.vstack(p)]).transpose()
        return plotarray.tolist()
    
    
    def illumination_vector(self, plane='xy'):
        """Generate an illumination vector based on a given plane.
           This vector essentially "simulates" where the light is coming from.
           Currently limited to directions parallel to global XYZ system.
           
           Returns: 3-by-1 np.ndarray
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
    
    
    def illuminated_faces(self, geometry: Geometry, plane='xy'):
        """Input a geometry and a illumination plane, and this function will
            return the faces in the geometry that are illuminated.
        """
        iv = self.illumination_vector(plane=plane)
        
        illuminated_faces = []
        
        for face in geometry.faces:
            if np.dot(iv, self.face_perpendicular(face)) < 0:
                illuminated_faces.append(face)
            else:
                pass
        return illuminated_faces
    
    
    def illuminated_strength(self, face: Face, plane='xy'):
        iv = self.illumination_vector(plane=plane)
        return np.arccos(np.dot(iv, self.face_perpendicular(face))) - np.pi/2
    
    
    def cosi_face(self, face: Face, vector):
        """Calculates the cosine of the angle between the perpendicular of a
           specified face and a specified vector.
           """
        return np.dot(vector, self.face_perpendicular(face))
    
    
    def area_projected(self, geometry: Geometry, vector):
        """Takes a geometry and a vector in a certain direction, expressed
           in the global frame. It then only looks at the faces of this
           geometry whose perpendiculars are pointing in the opposite 
           direction as the specified vector. Then it computes the projected
           area of these faces with respect to the specified vector.
           
           In other words, if you were to look at a geometry, and the 
           specified vector indicates the direction from which you observe it, 
           this function returns the area of the geometry you would observe.
           
           returns: A_projected_total (float)
           """
        A_partial = []
        
        for face in geometry.faces:
            cosi = self.cosi_face(face, vector)
            if cosi < 0:
                A_partial.append(-cosi*face.area())
            else:
                A_partial.append(0)
        return sum(A_partial)
        
    
    def project_face(self, face: Face, plane='xy', new_frame=None):
        """Project a copy of the face onto a plane that is spanned by two
            axes. The projection is orthographic. Coordinates will be in terms
            of the GLOBAL frame.
            
            TODO: A "new_frame" can be specified, which will make the 
            projection elements children of this frame. Otherwise, the 
            elements will be children of the global coordinate frame.
            
            TODO: Generalize this function for arbitrary projection plane. 
            """
        # Eliminate minus sign in front of the plane:    
        if plane[0] == '-':
            plane = plane[1:]
        
        pj_xyz = []
        
        for vertex in [face.p1, face.p2, face.p3, face.p4]:
            
            global_coords = self.vertex_xyz_global(vertex)
            
            if plane == 'xy':
                pj_x = global_coords[0]
                pj_y = global_coords[1]
                pj_z = 0
            elif plane == 'xz':
                pj_x = global_coords[0]
                pj_y = 0
                pj_z = global_coords[2]
            elif plane == 'yz':
                pj_x = 0
                pj_y = global_coords[1]
                pj_z = global_coords[2]
            else:
                raise ValueError("No valid projection plane given to "
                                  "projection method! "
                                  "Valid options: 'xy', 'xz', 'yz'")
            pj_xyz.append(Vertex([pj_x, pj_y, pj_z], parenttype="face"))
        
        # Create the projected face using the projected vertices.
        if not new_frame: 
            pj = Face(pj_xyz[0], pj_xyz[1], pj_xyz[2], pj_xyz[3],
                      parenttype="global")
        else:
            raise TypeError("new_frame functionality is not yet implemented.")
        return pj
        
    
    def illuminated_area(self, geometry: Geometry, plane='xy'):
        area = 0
        for face in self.illuminated_faces(geometry, plane):
            area += self.project_face(face, plane).area()
        return area
    
    
    def readout(self, dec=4):
        """More detailed information about Vertex.
           Use dec to define the number of decimals before rounding.
           Set dec to -1 to disable rounding.
           """
        if dec == -1:  # Set dec to -1 to disable rounding
            origin = self.origin()
        else:
            origin = np.round(self.origin(), dec)
        
        # children = len(self.vertices)   + \
        #            len(self.faces)      + \
        #            len(self.geometries) + \
        #            len(self.vectors)
        
        print("Carthesian coordinate frame.")
        print("Frame origin: {}\n".format(origin))
        # print("Children:              {}".format(children))
        # print("Associated Vertices:   {}".format(len(self.vertices)))
        # print("Associated Faces:      {}".format(len(self.faces)))
        print("Associated Geometries: {}".format(len(self.geometries)))
        # print("Associated Vectors:    {}".format(len(self.vectors)))