# -*- coding: utf-8 -*-
"""
Created on Mon May  3 08:57:58 2021

@author: Johan Monster

"""

import numpy as np
from numpy import sin, cos
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import mpl_toolkits.mplot3d as mp3d

from cp_vertex import Vertex
from cp_face import Face
from cp_geometry import Geometry
from cp_vector import Vector
from cp_frame import Frame
from cp_utilities import d2r#, r2d
from cp_plotting import plot_global_tripod, plot_frame,\
    plot_vertex, plot_face, plot_geometry_perpendiculars, plot_illumination,\
        plot_A_ill, plot_arrow


# Define some constant parameters:
ALB_earth_avg = 0.306  # [-]     Average geometric albedo of Earth
solar_flux = 1371      # [W/m2]  Solar flux in LEO
n_cell = 0.30          # [-]     Efficiency of a GaAs solar cell
n_packing = 0.604      # [-]     Cell packing efficiency on panels
degradation = 0.02     # [-]     Estimated yearly solar cell degradation factor
years_elapsed = 0      # [year]  How many years in orbit


#%% ==== Plotting ====

# Toggle plotting functionality:
if True:
    
    # Setting up the plot:
    fig = plt.figure(figsize=(10, 7))
    ax = mp3d.Axes3D(fig)
    
    """ TO CHANGE THE DEFAULT CAMERA VIEW, CHANGE THESE: """
    ax.view_init(elev=20, azim=-50)
    
    steps = 95
    angle_step = d2r(360/steps)
    
    p0 = Vertex([  0.5,   0.5,  0.5])
    p1 = Vertex([-0.05, -0.05, -0.1])
    p2 = Vertex([-0.05, -0.05,  0.1])
    p3 = Vertex([ 0.05, -0.05,  0.1])
    p4 = Vertex([ 0.05, -0.05, -0.1])
    p5 = Vertex([-0.05,  0.05, -0.1])
    p6 = Vertex([-0.05,  0.05,  0.1])
    p7 = Vertex([ 0.05,  0.05,  0.1])
    p8 = Vertex([ 0.05,  0.05, -0.1])
    
    fA = Face(p4, p3, p2, p1)
    fB = Face(p2, p3, p7, p6)
    fC = Face(p3, p4, p8, p7)
    fD = Face(p4, p1, p5, p8)
    fE = Face(p1, p2, p6, p5)
    fF = Face(p5, p6, p7, p8)
    
    cubesat = Geometry([fA, fB, fC, fD, fE, fF])
    
    frame1 = Frame()
    frame1.add_geometry(cubesat)
    frame1.translate(0.3, 0.3, 0.25)
    
    # frame1.rotate(0, 0, d2r(1.5*360/24))
    
    # Nadir pointing Z
    # cubesat.rotate(0,0,d2r(-1.5*360/24+45),cor=list(cubesat.find_cuboid_centroid()))
    # cubesat.rotate(0,0,d2r(-45),cor=list(cubesat.find_cuboid_centroid()))
    # cubesat.rotate(0,0,d2r(45),cor=list(cubesat.find_cuboid_centroid()))    
    # cubesat.rotate(d2r(25),0,0,cor=list(cubesat.find_cuboid_centroid()))    

    
    # DO A BARN DOOR!
    # cubesat.rotate(0,d2r(-90),0,cor=list(cubesat.find_cuboid_centroid()))
    
    # MOAR SURFACE!

    

    # Model the fake Earth
    E_anchor = Vertex([0,0,0])
    Ep = Face(E_anchor, Vertex([0,0.001,0]),
              Vertex([0.001,0.001,0]), Vertex([0,0.001,0]))
    E = Geometry([Ep])
    frameE = Frame()
    frameE.add_geometry(E)
    frameE.translate(0.3, 0.3, 0.25)
    frameE.rotate(0, 0, d2r(1.5*360/24))
    
    E.translate(0,.15,0)
    
    
    # Pre-allocate the illumination area vector
    A_ill = np.zeros(steps-1)
    A_ill2 = np.zeros(steps-1)
    multiplier = np.zeros(steps-1)

    P_sun = np.zeros(steps-1)
    P_alb = np.zeros(steps-1)
    
    
    
    
    def update(i):
        
        # Setting up the axes object
        ax.clear()
        
        
        
        plotscale = 0.5
        ax.set_xlim(0, 1.25*plotscale)
        ax.set_ylim(0, 1.25*plotscale)
        ax.set_zlim(0, plotscale)
        
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
        
        
        
        # Transforming cubesat
        cubesat.rotate(angle_step,0,0,cor=list(cubesat.find_cuboid_centroid()))
        
        # Transform fake Earth
        E.rotate(angle_step,0,0,cor=list(cubesat.find_cuboid_centroid()))



        # Find albedo vector:        
        albedo_vector = frame1.vertex_xyz_global(cubesat.make_cuboid_centroid()) \
            - frameE.vertex_xyz_global(E_anchor)
        # Normalize albedo vector to unit vector
        albedo_vector = albedo_vector/np.linalg.norm(albedo_vector)

        plot_arrow(ax, frameE.vertex_xyz_global(E_anchor), albedo_vector, scaling=0.05)



        # Finding albedo multiplier
        multiplier[i-1] = np.dot(albedo_vector, frame1.illumination_vector('xz'))
        if multiplier[i-1] >= 0:
            multiplier[i-1] = 0
        else:
            multiplier[i-1] = -multiplier[i-1]
        
        

        # Shitty shadow simulation:
        if True: # Toggle
            shadow = False  
            if 0.5-0.361/2 <= i/(steps-1) <= 0.5+0.361/2:
                shadow = True
            else:
                shadow = False
        
        global A_ill
        if not shadow:
            A_ill[i-1] = round(frame1.illuminated_area(cubesat, plane='xz'), 4)
            A_ill2[i-1] = round(frame1.area_projected(cubesat, albedo_vector),4)
        else:
            A_ill[i-1] = 0
            A_ill2[i-1] = 0
        
        
        
        # Calculate power
        
        # Power from sun:
        P_sun[i-1] = A_ill[i-1]\
            * solar_flux * n_cell * n_packing \
            * (1 - degradation * years_elapsed)
        
        # Power from albedo:
        P_alb[i-1] = A_ill2[i-1] * multiplier[i-1] * ALB_earth_avg \
            * solar_flux * n_cell * n_packing \
            * (1 - degradation * years_elapsed)    
        
        
        
        
        # print("[DEBUG] Frame {}, i/(steps-1) = {}, shadow {}"\
        #       .format(i, round(i/(steps-1),2), shadow))
        
        print("[DEBUG] Frame {}, A_ill = {}, A_ill2 = {}"\
              .format(i, A_ill[i-1], A_ill2[i-1]))        
        
            
        # Plotting the global XYZ tripod
        plot_global_tripod(ax, scaling=plotscale/2)
        
        # Plot tripod of frame1:
        if not shadow:
            plot_frame(ax, frame1, tripod_scale=plotscale/8,
                       perpfill=False, perpscale=0.1, 
                       illumination=True, ill_plane='xz',
                       facecolour="#0F0F0F", facealpha=0.7
                       )
            
            plot_illumination(ax, frame1, plane='xz')
            
        else:
            plot_frame(ax, frame1, tripod_scale=plotscale/8,
                       perpfill=False, perpscale=0.1, 
                       illumination=False, ill_plane='xz',
                       facecolour="#0F0F0F", facealpha=0.7
                       )
            
            plot_illumination(ax, frame1, plane='xz',linecolour="#000")
            
        
        plot_frame(ax, frameE, tripod_scale=plotscale/8, show_tripod=False,
                   linefill=False, facefill=False,
                   illumination=False,
                   vertexcolour="#0000DD", vertexsize=300
                   )

        ax.set_title("3D Visualization.    Frame: {}.    A_ill = {} m^2.    Shadow is {}.\
                     ".format(str(i), str(A_ill[i-1]), shadow))
        
    ani = animation.FuncAnimation(fig, update, np.arange(1,steps), 
                                  interval = 100, repeat = False)

    plt.show()



def plot_shit(A_ill: list, A_ill2: list, multiplier: list):
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
    plt.plot(range(len(A_ill2)), A_ill2, 'b')
    plt.plot(range(len(multiplier)), multiplier, 'gray')
    
    # Average area plot
    plt.plot([0, len(A_ill)], [A_avg, A_avg], 'r:')

def plot_power(P_sun: list, P_alb: list):
    # P_tot list:
    P_tot = [P_sun[i]+P_alb[i] for i in range(len(P_sun))]
    
    # Calculate average powers
    P_sun_avg = round(sum(P_sun)/len(P_sun), 4)
    P_alb_avg = round(sum(P_alb)/len(P_alb), 4)
    P_tot_avg = P_sun_avg + P_alb_avg
    
    print(" P_sun: ",round(P_sun_avg,2),"W \n",
          "P_alb: ",round(P_alb_avg,2),"W \n",
          "P_tot: ",round(P_tot_avg,2),"W \n",
          "P_alb/P_tot: ", round(P_alb_avg/P_tot_avg*100,1),"% ")
           
    # Set up plot
    fig_tmp = plt.figure(figsize=(10, 7))
    ax_tmp = fig_tmp.add_subplot(111)
    ax_tmp.set_ylim([0, 10])
    plt.title("Illuminated CubeSat power during simulation.")
    plt.xlabel("Time [min]")
    plt.ylabel("Power [W]")    
    # plt.text(0,0.035,"Average area: {} m^2".format(A_avg), color='r')
    # Area progression plot
    plt.plot(range(len(P_sun)), P_sun, 'r', label="P_sun")
    plt.plot(range(len(P_alb)), P_alb, 'b', label="P_albedo")
    plt.plot(range(len(P_tot)), P_tot, 'k', label="P_total")

    plt.text(33,9.6,"P_total:     {} W".format(round(P_tot_avg,2)), color='k')    
    plt.text(33,9.1,"P_sun:       {} W".format(round(P_sun_avg,2)), color='k')
    plt.text(33,8.8,"P_albedo:    {} W".format(round(P_alb_avg,2)), color='k')
    plt.text(33,8.5,"P_alb/P_tot: {} %".format(round(P_alb_avg/P_tot_avg*100,1)), 
             color='k')
    
    plt.legend()
    # Average area plot
    # plt.plot([0, len(A_ill)], [A_avg, A_avg], 'r:')