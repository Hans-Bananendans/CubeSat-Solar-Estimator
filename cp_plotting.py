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
from cp_utilities import hex2nRGB


def plot_vertex(axes: plt.matplotlib.axes, vertex: Vertex, 
                xyz_global = None,
                
                vertexfill = True,   # If False, vertex will not be plotted
                vertexcolour="#000", # Specifies the vertex colour
                vertexsize=10,       # Size of the plotted vertex
                vertexalpha=1        # Opacity of the plotted vertex
                ):
    
    # print("[DEBUG] Plotting {}".format(vertex))
    
    # Check whether vertex should be plotted:
    if vertexfill:
        # If vertex is an orphan, use its local coordinates for plotting
        # if vertex.parenttype == "global":
        if not xyz_global.all():
            x, y, z = vertex.xyz()
        # If not, check and use xyz_global coordinates provided
        else:
            # Verifying type of xyz_global
            if (isinstance(xyz_global, list) or isinstance(xyz_global, np.ndarray)):
                # Verifying length of xyz_global
                if len(xyz_global) == 3:
                    [x, y, z] = xyz_global
                else:
                    raise ValueError("Error: Tried to plot vertex using {}\
                                     coordinates, but exactly 3 are needed!\
                                     ".format(len(xyz_global)))
            # elif xyz_global == None:
            #     raise TypeError("Something went wrong with vertex plot!")
            else:
                raise TypeError("Error: Tried to plot vertex, but argument \
                                'xyz_global' was invalid type ({})!\
                                ".format(type(xyz_global)))
     
        axes.scatter(x, y, z,
                     c=vertexcolour, s=vertexsize, alpha=vertexalpha)

def plot_frame(axes: plt.matplotlib.axes, frame: Frame, 
               # Tripod properties
               show_tripod=True,     # If False, does not plot the tripod
               tripod_scale=1,       # Sets the scale of the tripod
               plot_perpendiculars=True, # Plot perpendiculars
               
               # Vertex plotting properties:
               vertexfill = True,    # If False, vertex will not be plotted
               vertexcolour="#000",  # Specifies the vertex colour
               vertexsize=10,        # Size of the plotted vertex
               vertexalpha=1,        # Opacity of the plotted vertex
               
               # Face plotting properties:
               linefill=True,        # If False, does not plot face lines
               linecolour="#000",    # Colour of face lines
               linewidth=2,          # Thickness of face lines
               linealpha=1,          # Opacity of face lines
               facefill=True,        # If False, does not shade the face area
               facecolour="#555",    # Colour of the face area shading
               facealpha=1,          # Opacity of the face area shading

               # Face perpendicular arrow plotting properties:
               perpfill = False,     # If True, plot perpendiculars
               perpcolour="#888",    # Specifies the perp. arrow colour
               perpscale=1,          # Size of the plotted perp. arrow
               perpalpha=0.5,        # Opacity of the plotted perp. arrow             

               # Illumination:
               illumination=False,   # If True, plots illumination intensity
               ill_value=0,          # Used to plot illumination intensity
               ill_plane=None,       # If illumination is used, a plane
                                     #  MUST be satisfied.
               
               # Vector plotting properties:
               vectorfill=True,      # If False, does not plot vector arrow
               vectorcolour="#000",  # Colour of vector arrow
               vectoralpha=1,        # Opacity of vector arrow
               vectorscale=1,        # Scale the whole vector by a constant
               vectorratio=0.15      # Vector arrow length ratio
               ):
    
    if frame.geometries:
        """If frame contains geometries, plot these. Then find out if there
           are any leftover edges and vertices associated with the frame, plot
           these too.
           """
        
        for geometry in frame.geometries:
            
            # Plot unique vertices first:
            for vertex in geometry.vertices():
                plot_vertex(axes, vertex, frame.vertex_xyz_global(vertex),
                            vertexfill = vertexfill, 
                            vertexcolour=vertexcolour,
                            vertexsize=vertexsize,
                            vertexalpha=vertexalpha)
            
            # Then plot faces, without plotting the vertices
            for face in geometry.faces:
                if illumination:
                    if not ill_plane:
                        raise AssertionError("If illumination is plotted, an \
                                             illumination plane MUST be \
                                             specified, e.g. \
                                             (ill_plane = 'xy')")
                    else:
                        if face not in frame.illuminated_faces(geometry,
                                                               plane=ill_plane
                                                               ):
                            ill_value = 0
                        else:
                            ill_value = frame.illuminated_strength(face,
                                                                   plane=ill_plane)
                        
                plot_face(axes, face, frame.plotlist(face),
                          linefill=linefill, 
                          linecolour=linecolour,
                          linewidth=linewidth, 
                          linealpha=linealpha,
                          facefill=facefill,
                          facecolour=facecolour,
                          facealpha=facealpha,
                          vertexfill=False,
                          vertexcolour=vertexcolour,
                          vertexsize=vertexsize,
                          vertexalpha=vertexalpha, 
                          illumination=illumination,
                          ill_value=ill_value)
           
            if perpfill:
                # First plot centroids of each plane:
                plot_geometry_perpendiculars(axes, \
                            frame.perpendiculars_plotlist(geometry, 
                                                          perpscale),
                            colour = perpcolour,
                            alpha = perpalpha)
                
    if show_tripod:
        plot_frame_tripod(axes, frame, scaling=tripod_scale)

def plot_geometry_perpendiculars(axes, perpendiculars_plotlist, \
                                  colour='#888', alpha=.5, perp_scale=1.):
    """Plots the normals/perpendiculars of each face in a geometry, and 
        displays them as little gray arrows.
        """
    xyzuvw = perpendiculars_plotlist
    
    # First plot arrow bases as dots:
    axes.scatter(xyzuvw[0], xyzuvw[1], xyzuvw[2], c=colour, s=2, alpha=alpha)

    # Then, attach plane-perpendicular arrows to these points:
    axes.quiver(xyzuvw[0][:], xyzuvw[1][:], xyzuvw[2][:],
                xyzuvw[3][:], xyzuvw[4][:], xyzuvw[5][:],
                color=colour, alpha=alpha)


def plot_face(axes: plt.matplotlib.axes, face: Face, 
              plotlist: list = None,
              
              # Face plotting properties:
              linefill=True,        # If False, does not plot face lines
              linecolour="#000",    # Colour of face lines
              linewidth=2,          # Thickness of face lines
              linealpha=1,          # Opacity of face lines
              facefill=True,        # If False, does not shade the face area
              facecolour="#555",    # Colour of the face area shading
              facealpha=1,          # Opacity of the face area shading
              
              # Vertex plotting properties:
              vertexfill = True,    # If False, vertex will not be plotted
              vertexcolour="#000",  # Specifies the vertex colour
              vertexsize=10,        # Size of the plotted vertex
              vertexalpha=1 ,       # Opacity of the plotted vertex
              
              # Illumination:
              illumination=False,   # If True, plots illumination intensity
              ill_value=0           # Used to plot illumination intensity
              ):
    """ Plots an individual Face object. Check source code for kwargs.
    """
    
    def plotlist2plotlist_xyz(plotlist: tuple):
        """Return a list of lists with the xyz coordinates of each vertex."""
        # If object has no parent, return plot list in terms of global frame:
        
        p1_xyz = []
        p2_xyz = []
        p3_xyz = []
        p4_xyz = []
        
        for coordinate_row in plotlist:
            p1_xyz.append(coordinate_row[0])
            p2_xyz.append(coordinate_row[1])
            p3_xyz.append(coordinate_row[2])
            p4_xyz.append(coordinate_row[3])
            
        return [p1_xyz, p2_xyz, p3_xyz, p4_xyz]
    
    # Situations:
    # 1. Face is an orphan (parent == "global") and so it has to use
    #    local coordinates only
    # 2. Face is part of a geometry, which we will assume has a frame underlying
    #    it
        
    # print("[DEBUG] Plotting {}".format(face))
    
    # Case 1: Face is an orphan:
    if face.parenttype == "global":
        (xlist, ylist, zlist) = face.plotlist_local()

        # Plot individual vertices:
        if vertexfill: 
            for vertex in face.vertices():
                plot_vertex(axes, vertex)

        # Plot edges that connect vertices:
        for i in range(-1, len(xlist) - 1):
            axes.plot3D([xlist[i], xlist[i + 1]],
                        [ylist[i], ylist[i + 1]],
                        [zlist[i], zlist[i + 1]],
                        linecolour, alpha=linealpha, lw=linewidth)
    
        # Plot the face surface:
        if facefill:
            plotlist_xyz = plotlist2plotlist_xyz(face.plotlist_local())
            fs = mp3d.art3d.Poly3DCollection([plotlist_xyz],
                                              linewidth=0)
            
            (nR, nG, nB) = hex2nRGB(facecolour)
            fs.set_facecolor((nR, nG, nB, facealpha))
            axes.add_collection3d(fs)
    
    # Case 2: Face is in a frame
    else:
        (xlist, ylist, zlist) = plotlist

        # Plot individual vertices:
        # if vertexfill: 
        #     for vertex in face.vertices():
        #         plot_vertex(axes, vertex)

        # Plot edges that connect vertices:
        for i in range(-1, len(xlist) - 1):
            axes.plot3D([xlist[i], xlist[i + 1]],
                        [ylist[i], ylist[i + 1]],
                        [zlist[i], zlist[i + 1]],
                        linecolour, alpha=linealpha, lw=linewidth)
    
        # Plot the face surface:
        if facefill:
            plotlist_xyz = plotlist2plotlist_xyz(plotlist)
            fs = mp3d.art3d.Poly3DCollection([plotlist_xyz],
                                              linewidth=0)
            if not illumination:
                (nR, nG, nB) = hex2nRGB(facecolour)
                fs.set_facecolor((nR, nG, nB, facealpha))
                axes.add_collection3d(fs)
            else:
                value = round(0.1 + 0.9 / 1.6 * ill_value, 2)
                fs.set_facecolor((value, value, value, facealpha))
                axes.add_collection3d(fs)
                

     
def plot_illumination(axes: plt.matplotlib.axes, frame: Frame, 
                      plane='xy', geometry: Geometry=None, 
                      
                      # Tripod properties
                      show_tripod=False,    # If False, does not plot the tripod
                      tripod_scale=1,       # Sets the scale of the tripod
                      plot_perpendiculars=False, # Plot perpendiculars
                      
                      # Vertex plotting properties:
                      vertexfill = False,   # If False, vertex will not be plotted
                      vertexcolour="FFC125",# Specifies the vertex colour
                      vertexsize=10,        # Size of the plotted vertex
                      vertexalpha=1,        # Opacity of the plotted vertex
                      
                      # Face plotting properties:
                      linefill=True,        # If False, does not plot face lines
                      linecolour="#FFC125", # Colour of face lines
                      linewidth=2,          # Thickness of face lines
                      linealpha=1,          # Opacity of face lines
                      facefill=False,       # If False, does not shade the face area
                      facecolour="FFC125",  # Colour of the face area shading
                      facealpha=1,          # Opacity of the face area shading
                       
                      # Face perpendicular arrow plotting properties:
                      perpfill = False,     # If True, plot perpendiculars
                      perpcolour="#888",    # Specifies the perp. arrow colour
                      perpscale=1,          # Size of the plotted perp. arrow
                      perpalpha=0.5,        # Opacity of the plotted perp. arrow             
                       
                      # Illumination:
                      illumination=False,   # If True, plots illumination intensity
                      ill_value=0,          # Used to plot illumination intensity
                      
                      # Vector plotting properties:
                      vectorfill=False,     # If False, does not plot vector arrow
                      vectorcolour="#000",  # Colour of vector arrow
                      vectoralpha=1,        # Opacity of vector arrow
                      vectorscale=1,        # Scale the whole vector by a constant
                      vectorratio=0.15      # Vector arrow length ratio
                      ):
    # Gather the faces to be plotted in the geometry. If no geometry was
    #   provided (geometry=None), plot all the illuminated faces in the frame.
    
    if not geometry:
        geometries = frame.geometries
    else:
        geometries = list(geometry)
    
    illuminated_faces = []
    
    # Find out what faces are illuminated:
    for geo in geometries:
        faces = frame.illuminated_faces(geo, plane=plane)
        illuminated_faces.extend(faces)
    
    # Plot each face after projecting them in the global frame:
    for face in illuminated_faces:
        # projected_face = frame.project_face(face)
        plot_face(axes, frame.project_face(face, plane=plane), 
              linefill=linefill, linecolour=linecolour, linewidth=linewidth, 
              linealpha=linealpha, facefill=facefill, facecolour=facecolour,
              facealpha=facealpha, vertexfill=vertexfill, 
              vertexcolour=vertexcolour, vertexsize=vertexsize, 
              vertexalpha=vertexalpha, illumination=False)
        
        
    
    
    
            
# def plot_geometry(axes: plt.matplotlib.axes, geometry: Geometry,
#                   # Face plotting properties:
#                   linefill=True,        # If False, does not plot face lines
#                   linecolour="#000",    # Colour of face lines
#                   linewidth=2,          # Thickness of face lines
#                   linealpha=1,          # Opacity of face lines
#                   facefill=True,        # If False, does not shade the face area
#                   facecolour="#000",    # Colour of the face area shading
#                   facealpha=1,          # Opacity of the face area shading
                  
#                   # Vertex plotting properties:
#                   vertexfill = True,    # If False, vertex will not be plotted
#                   vertexcolour="#000",  # Specifies the vertex colour
#                   vertexsize=10,        # Size of the plotted vertex
#                   vertexalpha=1 ,       # Opacity of the plotted vertex
                  
#                   # Illumination:
#                   illumination=False,   # If True, plots illumination intensity
#                   ill_value=0           # Used to plot illumination intensity
#                   ):
#     """Plots an individual Face object. Check source code for kwargs."""
    
#     # print("[DEBUG] Plotting {}".format(geometry))
    
#     for face in geometry.faces:
#         plot_face(axes, face)

def plot_vector(axes: plt.matplotlib.axes, vector: Vector):
    """Plots an individual vector object.."""
    pass

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

def plot_arrow(axes, base, head, 
                      alpha=1, scaling=1.):
    """Plots an rgb xyz tripod at the global origin.
        - alpha changes the opacity of the arrows. (default: 1)
        - scaling changes the length of the arrows (default: 1)
        """
    sc = scaling
    axes.quiver(base[0], base[1], base[2], 
                sc*head[0], sc*head[1], sc*head[2],
              arrow_length_ratio=0.15, color='purple', alpha=alpha)


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
    
# def plot_frame(axes: plt.matplotlib.axes, frame: Frame, 
#                # Tripod properties
#                show_tripod=True,     # If False, does not plot the tripod
#                tripod_scale=1,       # Sets the scale of the tripod
               
#                # Face plotting properties:
#                linefill=True,        # If False, does not plot face lines
#                linecolour="#000",    # Colour of face lines
#                linewidth=2,          # Thickness of face lines
#                linealpha=1,          # Opacity of face lines
#                facefill=True,        # If False, does not shade the face area
#                facecolour="#000",    # Colour of the face area shading
#                facealpha=1,          # Opacity of the face area shading
              
#                # Vertex plotting properties:
#                vertexfill = True,    # If False, vertex will not be plotted
#                vertexcolour="#000",  # Specifies the vertex colour
#                vertexsize=10,        # Size of the plotted vertex
#                vertexalpha=1,        # Opacity of the plotted vertex
              
#                # Illumination:
#                illumination=False,   # If True, plots illumination intensity
#                ill_value=0,          # Used to plot illumination intensity
               
#                # Vector plotting properties:
#                vectorfill=True,      # If False, does not plot vector arrow
#                vectorcolour="#000",  # Colour of vector arrow
#                vectoralpha=1,        # Opacity of vector arrow
#                vectorscaling=1,      # Scale the whole vector by a constant
#                vectorratio=0.15      # Vector arrow length ratio
#                ):
    
#     if show_tripod:
#         plot_frame_tripod(axes, frame, scaling=tripod_scale)
    
#     if frame.geometries:
#         """If frame contains geometries, plot these. Then find out if there
#            are any leftover edges and vertices associated with the frame, plot
#            these too.
#            """
        
#         # Keep track of plotted faces and vertices
#         plotted_faces = []
#         plotted_vertices = []
        
#         # Plot the geometries first
#         for geometry in frame.geometries:
#             plot_geometry(axes, geometry,
#                           linefill=linefill, linecolour=linecolour,
#                           linewidth=linewidth, linealpha=linealpha,
#                           facefill=facefill,facecolour=facecolour,
#                           facealpha=facealpha,
#                           vertexfill = vertexfill, vertexcolour=vertexcolour,
#                           vertexsize=vertexsize,vertexalpha=vertexalpha, 
#                           illumination=illumination,ill_value=ill_value
#                           )
#             # Mark all faces in the geometry as plotted
#             for face in geometry.faces:
#                 plotted_faces.append(face)
#                 # Mark all vertices in these faces as plotted
#                 for vertex in face.vertices():
#                     plotted_vertices.append(vertex)
        
#         # Find faces in the frame that do not belong to any geometry
#         for face in set(frame.faces).difference(set(plotted_faces)):
#             # Plot these faces too
#             plot_face(axes, face,
#                       linefill=linefill, linecolour=linecolour,
#                       linewidth=linewidth, linealpha=linealpha,
#                       facefill=facefill,facecolour=facecolour,
#                       facealpha=facealpha,
#                       vertexfill = vertexfill, vertexcolour=vertexcolour,
#                       vertexsize=vertexsize,vertexalpha=vertexalpha, 
#                       illumination=illumination,ill_value=ill_value
#                       )
#             # Then mark all vertices in these faces as plotted too
#             for vertex in face.vertices():
#                     plotted_vertices.append(vertex)
        
#         # Find the vertices in the frame that do not belong to any faces
#         for vertex in set(frame.vertices).difference(set(plotted_vertices)):
#             # Then plot these too
#             plot_vertex(axes, vertex,
#                         vertexfill = vertexfill, vertexcolour=vertexcolour,
#                         vertexsize=vertexsize,vertexalpha=vertexalpha
#                         )
#         # By now, all geometries, faces, and vertices in the frame should
#         # have been plotted.
            
#     elif not frame.geometries and frame.faces:
#         """If frame contains no geometries, but does contain faces, plot
#            these faces. Then find out if there are any leftover vertices in
#            the frame, and plot these too.
#            """
           
#         # Keep track of plotted vertices
#         plotted_vertices = []
        
#         # Plot all faces first
#         for face in frame.faces:
#             plot_face(axes, face,
#                       linefill=linefill, linecolour=linecolour,
#                       linewidth=linewidth, linealpha=linealpha,
#                       facefill=facefill,facecolour=facecolour,
#                       facealpha=facealpha,
#                       vertexfill = vertexfill, vertexcolour=vertexcolour,
#                       vertexsize=vertexsize,vertexalpha=vertexalpha, 
#                       illumination=illumination,ill_value=ill_value
#                       )
#             # Then mark all vertices in these faces as plotted
#             for vertex in face.vertices():
#                     plotted_vertices.append(vertex)
        
#         # Find the vertices in the frame that do not belong to any faces
#         for vertex in set(frame.vertices).difference(set(plotted_vertices)):
#             # Then plot these too
#             plot_vertex(axes, vertex,
#                         vertexfill = vertexfill, vertexcolour=vertexcolour,
#                         vertexsize=vertexsize,vertexalpha=vertexalpha
#                         )
#         # By now, all faces and vertices in the frame should be plotted.
           
#     elif not frame.geometries and not frame.faces and frame.vertices:
#         """If frame contains no geometries nor faces, but does contain some
#         vertices, plot these.
#         """
#         for vertex in frame.vertices:
#             plot_vertex(axes, vertex,
#                         vertexfill = vertexfill, vertexcolour=vertexcolour,
#                         vertexsize=vertexsize,vertexalpha=vertexalpha
#                         )
#     else:
#         """In other cases, do nothing."""
#         pass
    
#     # Plot all vectors in frame too
#     for vector in frame.vectors:
#         plot_vector(axes, vector)

def plot_A_ill(A_ill: list):
    # Calculate average area
    A_avg = round(sum(A_ill)/len(A_ill), 4)
    
    # Set up plot
    fig_tmp = plt.figure(figsize=(10, 7))
    ax_tmp = fig_tmp.add_subplot(111)
    ax_tmp.set_ylim([0, 0.04])
    plt.title("Illuminated CubeSat area during simulation.")
    plt.xlabel("Simulation steps")
    plt.ylabel("Illuminated area [m^2]")    
    plt.text(0,0.035,"Average area: {} m^2".format(A_avg), color='r')
    # Area progression plot
    plt.plot(range(len(A_ill)), A_ill, 'k')
    
    # Average area plot
    plt.plot([0, len(A_ill)], [A_avg, A_avg], 'r:')