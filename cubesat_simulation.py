# -*- coding: utf-8 -*-
"""
Created on Mon May  3 08:57:58 2021

@author: Johan Monster

"""

import numpy as np
#from numpy import sin, cos
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import mpl_toolkits.mplot3d as mp3d

from cp_vertex import Vertex
from cp_face import Face
from cp_geometry import Geometry
#from cp_vector import Vector
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
    ax = mp3d.Axes3D(fig, auto_add_to_figure=False)
    fig.add_axes(ax)
    
    
    """========= TO CHANGE THE DEFAULT CAMERA VIEW, CHANGE THESE: ========="""
    ax.view_init(elev=20, azim=-50)
    
    # Define how many steps the simulation/animation will have:
    # High values give a more granular curve, but it takes longer to compute
    global steps 
    steps = 64
    # steps is a global variable, so we can access it outside the update()
    #   function. It should be noted that relying on global variables is
    #   bad practice, but I'm an adult and I'll do what I like!
    
    # Defining an angle step that ensures the CubeSat model rotates completely,
    #   but only once.
    angle_step = d2r(360/steps)
    
    
    """ ================= DEFINING THE GEOMETRIC OBJECTS ================= """
    # Defining the vertices that will be the corners of the cubesat. I 
    #   recommend that you keep the units in meters.
                                        #         ^ Z
                                        #         | 
    p1 = Vertex([-0.05, -0.05, -0.1])   #      2 ___ 6
    p2 = Vertex([-0.05, -0.05,  0.1])   #       |\3__\7 
    p3 = Vertex([ 0.05, -0.05,  0.1])   #       | |  |         Y
    p4 = Vertex([ 0.05, -0.05, -0.1])   #      1| | 5|  -------->
    p5 = Vertex([-0.05,  0.05, -0.1])   #        \|__|
    p6 = Vertex([-0.05,  0.05,  0.1])   #        4    8  
    p7 = Vertex([ 0.05,  0.05,  0.1])   #            \
    p8 = Vertex([ 0.05,  0.05, -0.1])   #             v X
    
    # Assembling the side panels of the Cubesat. Note that the order in which
    #   the points is specified DOES matter! (see definition of Face class)
    
                                #         ^ Z
                                #         | 
    fA = Face(p4, p3, p2, p1)   #     E  ___  
    fB = Face(p2, p3, p7, p6)   #       |\ B_\  
    fC = Face(p3, p4, p8, p7)   #     A | |  | F        Y
    fD = Face(p4, p1, p5, p8)   #       | | C|  -------->
    fE = Face(p1, p2, p6, p5)   #        \|__|
    fF = Face(p5, p6, p7, p8)   #          D       
                                #            \
                                #             v X
    
    # Assembling a Geometry instance named 'cubesat', and add all the faces
    #   we just defined to the geometry.
    cubesat = Geometry([fA, fB, fC, fD, fE, fF])
    
    # Define a reference frame which will be attached to the centre of gravity
    #   of the CubeSat. If you move/rotate this frame, the 'cubesat' 
    #   Geometry will move/rotate along with it (if it is attached of course).
    frame1 = Frame()
    # Attaching the cubesat Geometry to 'frame1'
    frame1.add_geometry(cubesat)
    # Translate 'frame1' away from the origin. The sole reasons for doing this
    #   is so the geometry will not overlap with the yellow projection of the
    #   illuminated area, and because the plot will generally look nicer.
    frame1.translate(0.3, 0.3, 0.25)
    
    
    """ == HERE WE START CHANGING THE INITIAL ATTITUDE OF THE SATELLITE == """
    # Change the LTAN by 1.5 hours
    frame1.rotate(0, 0, d2r(1.5/24*360))
    

    """EXAMPLE: Nadir pointing with +Z"""    
    # The way we've set it up, at the start of the simulation, the satellite
    #   will face Earth with its +Y face. We could change this, if we wish.
    # For example, to make the satellite nadir-pointing with the +Z face
    #   instead, we have to rotate by -90 degrees around the X-axis of frame1.
    # We have to also make sure that we rotate around the Centre of Gravity
    #   of the CubeSat, rather than some other point. So, we set the "cor"
    #   (=Centre Of Rotation) equal to the Centre of Gravity of 'cubesat':

#    cubesat.rotate(d2r(-90),0,0,cor=list(cubesat.find_cuboid_centroid()))  


    """EXAMPLE: Pointing the corner at p8 towards Earth"""
    # Its not a great representation, but this looks like what we want to
    #   achieve. We decide that the angle between line segment p7-p8 and the 
    #   dotted line should be 45 degrees.
    #
    #       45*:                  ||
    #      p7  :                 ||
    #    .-'\  :                ||
    #    \ \ \ :                ||  EARTH
    #     \ \ \: p8    ---->    ||
    #      \-\`                 ||
    #                            ||
    #                             ||
    # To achieve this, we need to make two rotations. We need to rotate the
    #   CubeSat model around the local Z-axis by +45 degrees, and also rotate
    #   around the local X-axis by 45 degrees. 
    # The Geometry.rotate() function uses a 3-2-1 rotation sequence by default.
    
#    cubesat.rotate(      0, 0, d2r(45), cor=list(cubesat.find_cuboid_centroid()))
#    cubesat.rotate(d2r(45), 0,       0, cor=list(cubesat.find_cuboid_centroid()))

    
    """EXAMPLE: Fly in 'barndoor' config instead of 'bus' config"""
    # In this case, we have to rotate either:
    #   +90 degrees around local Y-axis (-X face will be in direction of velocity)
    # or
    #   -90 degrees around local Y-axis (+X face will be in direction of velocity)
    
#    cubesat.rotate(0,d2r(90),0,cor=list(cubesat.find_cuboid_centroid()))
    
  
    """ =========================== FAKE EARTH ========================= """
    # In this code segment, we create a small Geometry which will represent 
    #   the Earth. This is both for visualization purposes, but also to
    #   calculate the light coming in from the albedo effect. In the plot, we
    #   will focus on the satellite, and have the Earth "rotate" around the 
    #   satellite, rather than the other way around. This may look a bit
    #   strange at first, but it will make the plot look neater (the projected
    #   areas of the satellite will remain in the same place).
    
    # We first create a small vertex at the origin
    E_anchor = Vertex([0,0,0])
    # We cannot add single vertices to a Geometry currently, so we do a little
    #   hack: We're going to make a really tiny face and just render the points
    #   really big later on, so it'll look like a blue dot regardless.
    Ep = Face(E_anchor, Vertex([0,0.001,0]),
              Vertex([0.001,0.001,0]), Vertex([0,0.001,0]))
    # Add the tiny face to a geometry 'E'
    E = Geometry([Ep])
    # Define a new frame 'frameE'
    frameE = Frame()
    # Attach geometry 'E' to 'frameE'
    frameE.add_geometry(E)
    # Translate it by the same amount as we translated 'frame1'
    frameE.translate(0.3, 0.3, 0.25)
    # Rotate it by the same amount as we rotated 'frame1'
    frameE.rotate(0, 0, d2r(1.5*360/24))
    # Then we translate geometry E by 0.15 in the local (frameE) Y-direction
    #   so that the fake earth does not collide with the CubeSat model.
    E.translate(0,.15,0)
    
    
    """Pre-allocating storage arrays. """
    # Pre-allocate the illumination area arrays.
    # A_ill  will store the total area of the CubeSat faces directly
    #   illuminated by the sun (primary illumination)
    # A_ill2 will store the total area of the CubeSat faces that are pointing
    #   at the Earth (for albedo)
    A_ill = np.zeros(steps-1)
    A_ill2 = np.zeros(steps-1)
    
    # Then we define a multiplier array. This store values of the albedo
    #   multiplier. Basically, the numerical value of this is just the dot 
    #   product between the vector pointing from the sun to the satellite and
    #   the vector pointing from the Earth to the satellite. We assume that 
    #   the albedo is proportional to this angle (its value will follow a
    #   sinusoidal as the satellite passes in front of the Earth.)
    multiplier = np.zeros(steps-1)

    # Pre-allocate some more empty arrays for the generated power due to direct
    #   illumination by the sun, indirect illumination due to albedo, and the
    #   sum of the two (P_tot).
    P_sun = np.zeros(steps-1)
    P_alb = np.zeros(steps-1)
    P_tot = np.zeros(steps-1)
    
    
    """ ================== DEFINING THE UPDATE FUNCTION ================== """
    # Now we define the update() function. It is necessary to define this for 
    #   matplotlib.animation.FuncAnimation() to work correctly.
    # This function will run every time a new frame is made by
    #   FuncAnimation(). Not only do we define here what happens to the plot,
    #   but we also use it to manipulate the satellite incrementally.
    
    def update(i):
        
        """ Defining some plotting parameters. """
        # Clearing the plot with each update, so matplotlib doesn't keep
        #   rendering overtop one another in the plot window.
        ax.clear()
        
        # Set the scale of the plot with a sizing factor.
        plotscale = 0.5
        # Set the limits of the X/Y/Z axes of the plot. You can change these,
        #   but the CubeSat may start looking oddly fat or skinny.
        ax.set_xlim(0, 1.25*plotscale)
        ax.set_ylim(0, 1.25*plotscale)
        ax.set_zlim(0, plotscale)
        
        # Labelling the axes
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
        
        
        """ ============== TRANSFORMING THE CUBESAT MODEL ================ """
        # We could rotate the whole CubeSate by rotating its parent frame
        #   'frame1'. You would do this like so:
        # frame1.rotate(angle_step,0,0,cor=frame1.origin())
        # This is bad, because it forces us to define the rotations of the 
        #   CubeSat in terms of the global frame of reference. This is hard.
        
        # However, instead we will merely rotate the 'cubesat' Geometry and
        #   keep 'frame1' stationary. This way, we can define the rototations
        #   of the CubeSat in terms of the local frame 'frame1'.
        # Note that we rotate by angle_step each time the function updates, 
        #   and that we rotate 'cubesat' around the CubeSat centroid.
        
        cubesat.rotate(angle_step,0,0,cor=list(cubesat.find_cuboid_centroid()))
        
        # We rotate the fake earth Geometry 'E' by the same amount, and also
        #   around the CubeSat centroid.
        E.rotate(angle_step,0,0,cor=list(cubesat.find_cuboid_centroid()))


        """ ================ SETTING UP ALBEDO SIMULATION ================= """
        
        # Calculate the albedo vector. This is needed to determine the albedo
        #   multiplier. We will also plot it in purple, so it will be easy
        #   to understand what is going on:        
        albedo_vector = frame1.vertex_xyz_global(cubesat.make_cuboid_centroid()) \
            - frameE.vertex_xyz_global(E_anchor)
            
        # Normalize albedo vector to unit vector
        albedo_vector = albedo_vector/np.linalg.norm(albedo_vector)

        # Add the albedo vector to the plot so it is easy to visualize.
        plot_arrow(ax, 
                   frameE.vertex_xyz_global(E_anchor), 
                   albedo_vector, 
                   scaling=0.05)


        # Finding albedo multiplier
        multiplier[i-1] = np.dot(albedo_vector, frame1.illumination_vector('xz'))
        # If the satellite is behind the Earth, no reflected albedo light
        #   reaches the satellite. If the satellite is behind the Earth, the 
        #   dot product between the illumination vector (from Sun->CubeSat)
        #   and the albedo vector will be positive. If this is the case, we 
        #   want to make the multiplier 0 (no albedo).
        if multiplier[i-1] >= 0:
            multiplier[i-1] = 0
        # Else, we just set the multiplier to positive and assign it.
        #   The value of multiplier goes from 0 to 1.
        else:
            multiplier[i-1] = -multiplier[i-1]
        
        
        """ ================ SETTING UP SHADOW SIMULATION ================ """        
        # Shitty shadow simulation: Here we use a simple if-statement to check
        #   whether the satellite is behind the Earth. If it is, none of the
        #   satellite sides are illuminated, and so we just set the illuminated
        #   area to zero.
        
        if True: # If you set this False, you turn the shadow simulation off
            shadow = False  
            if 0.5-0.361/2 <= i/(steps-1) <= 0.5+0.361/2:
                shadow = True
            else:
                shadow = False
        
        # Unnecessary but I'm keeping it here just in case
#        global A_ill
        
        # If the CubeSat is not in the Earth's shadow, calculate the 
        #   illuminated area (from sun) and the area subjected to albedo.
        if not shadow:
            A_ill[i-1] = round(frame1.illuminated_area(cubesat, plane='xz'), 4)
            A_ill2[i-1] = round(frame1.area_projected(cubesat, albedo_vector),4)
        else:
            A_ill[i-1] = 0
            A_ill2[i-1] = 0
        
        
        """ ================= CALCULATE THE POWER VALUES ================= """    
        
        # Power from sunlight:
        P_sun[i-1] = A_ill[i-1]\
            * solar_flux * n_cell * n_packing \
            * (1 - degradation * years_elapsed)
        
        # Power from albedo:
        P_alb[i-1] = A_ill2[i-1] * multiplier[i-1] * ALB_earth_avg \
            * solar_flux * n_cell * n_packing \
            * (1 - degradation * years_elapsed)    
        
        # Calculate the total power
        P_tot[i-1] = P_sun[i-1] + P_alb[i-1]
        
        # Debug line to monitor A_ill, A_ill2 for each frame.        
#        print("[DEBUG] Frame {}, A_ill = {}, A_ill2 = {}"\
#              .format(i, A_ill[i-1], A_ill2[i-1]))        
        
        
        """ ======================= MAKE THE 3D PLOT ====================== """  
        # This section uses several variables from cp_plotting.py
        
        # Plotting the global XYZ tripod
        plot_global_tripod(ax, scaling=plotscale/2)
        
        
        if not shadow:
            # Sunlight shining on satellite? Plot contents of 'frame1', and use
            #   illumination values to show which sides are illuminated
            plot_frame(ax, frame1, tripod_scale=plotscale/8,
                       perpfill=False, perpscale=0.1, 
                       illumination=True, ill_plane='xz',
                       facecolour="#0F0F0F", facealpha=0.7
                       )
            
            # Sunlight shining on satellite? Plot illuminated area in yellow
            plot_illumination(ax, frame1, plane='xz')
            
        else:
            # Satellite in shadow? Plot contents of 'frame1', but disable
            #   illumination.
            plot_frame(ax, frame1, tripod_scale=plotscale/8,
                       perpfill=False, perpscale=0.1, 
                       illumination=False, ill_plane='xz',
                       facecolour="#0F0F0F", facealpha=0.7
                       )
            
            # Satellite in shadow? Plot illuminated area in black        
            plot_illumination(ax, frame1, plane='xz',linecolour="#000")
            
        
        # Plot the fake earth by plotting the contents of 'frameE'
        plot_frame(ax, frameE, tripod_scale=plotscale/8, show_tripod=False,
                   linefill=False, facefill=False,
                   illumination=False,
                   vertexcolour="#0000DD", vertexsize=300
                   )

        # Set plot title
        ax.set_title("3D Visualization.    Frame: {}.    A_ill = {} m^2.    Shadow is {}.\
                     ".format(str(i), str(A_ill[i-1]), shadow))
    
        
        # If simulation is at the final step, plot the power curve:
        if i == steps-1:
            print("Function update() is done! Plotting power curves...")
            
            # Call power curve plotting function
            plot_power(P_sun, P_alb)
    
    """ =================== CALL THE ANIMATION FUNCTION =================== """  
    ani = animation.FuncAnimation(fig, update, np.arange(1,steps), 
                                  interval = 100, repeat = False)

    plt.show()


""" ===== Define plot_power(), to visualize the progression of power ===== """
# Call this after the simulation has gone through all its steps

def plot_power(P_sun: list, P_alb: list):
    # P_tot list (we calculate it again to minimize input arguments):
    P_tot = [P_sun[i]+P_alb[i] for i in range(len(P_sun))]
    
    # Calculate average powers
    P_sun_avg = round(sum(P_sun)/len(P_sun), 4)
    P_alb_avg = round(sum(P_alb)/len(P_alb), 4)
    P_tot_avg = P_sun_avg + P_alb_avg
    
    print("P_sun: ",round(P_sun_avg,2),"W \n",
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


"""Deprecated function"""
#def plot_shit(A_ill: list, A_ill2: list, multiplier: list):
#    # Calculate average area
#    A_avg = round(sum(A_ill)/len(A_ill), 4)
#    
#    # Set up plot
#    fig_tmp = plt.figure(figsize=(10, 7))
#    ax_tmp = fig_tmp.add_subplot(111)
#    ax_tmp.set_ylim([0, 0.04])
#    plt.title("Illuminated CubeSat area during simulation.")
#    plt.xlabel("Simulation steps")
#    plt.ylabel("Illuminated area [m^2]")    
#    plt.text(0,0.035,"Average area: {} m^2".format(A_avg), color='r')
#    # Area progression plot
#    plt.plot(range(len(A_ill)), A_ill, 'k')
#    plt.plot(range(len(A_ill2)), A_ill2, 'b')
#    plt.plot(range(len(multiplier)), multiplier, 'gray')
#    
#    # Average area plot
#    plt.plot([0, len(A_ill)], [A_avg, A_avg], 'r:')