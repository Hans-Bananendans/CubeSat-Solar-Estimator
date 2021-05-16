# 3D_Graphics
Collection of code fragments dealing with three-dimensional manipulation and rendering. Primary use for the author is propagation of a CubeSat model for estimating the effects of illumination over time. Matplotlib is used as the rendering engine.

## Example output in Matplotlib
![image](./demo_animation.gif?raw=true)

## Code structure: 
 * Adds a **Vertex** class that defines a vertex in 3D space. 
 * Four instances of the Vertex class can then be used to make a quadrilateral **Face**. A Face has only one side, and so its area is computed from one side only. The order in which the vertices are supplied does matter, and should be done with care.
 * A **Geometry** is a collection of Face objects. These can be translated and rotated by manipulating the underlying Vertex and Face objects. 
 * A **Frame** is a Cartesian coordinate frame to which Vertex, Face, and Geometry objects can be attached. Translations and rotations of this Frame will then be propagated to all the objects that are attached to the frame.
 * A number of utility functions, such as angle transformations, can be found in *cp_utilities.py*.
 * Various plotting functions for plotting geometries with matplotlib can be found in *cp_plotting.py*
 
![alt text](./layout2.0.png?raw=true)

## Python dependencies
 * **numpy** for vector math functions and ndarrays.
 * **matplotlib** for use as rendering engine. 
 * **mpl_toolkits.mplot3d** for plotting in 3D.