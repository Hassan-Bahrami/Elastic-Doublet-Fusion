# Elastic Doublet Fusion
This is a Python project for simulating an elastic doublet fusion using the SPH method.
The code is a part of the implementation of the "Physically-based simulation of elastic-plastic fusion of 3D bioprinted spheroids" article by Bahrami et al. (2023) published in the "Biofabrication" journal. The full article can be found in the below link:

[https://iopscience.iop.org/article/10.1088/1758-5090/acf2cb/meta](https://iopscience.iop.org/article/10.1088/1758-5090/acf2cb/meta)

## Code Description
The project has four main Python files:
- Particles.py
- Physics.py
- Initialization.py
- glWidget.py

### Particles.py
This file includes a class for particles we use in the SPH method. Each particle has multiple attributes, such as mass, density, etc.

### Physics.py
This file contains the main computation of the SPH methods. The main class is called "Environment", which refers to the physical environment of the particles and includes various functions for computing the physical quantities of the particles, such as density, pressure, momentum matrix, elastic force, etc.

### Initialization.py
This file initializes two spheres (doublet), samples particles inside them, and computes the initial values of some particle attributes such as SPH kernel radius, particle mass, initial densities, volumes, etc.

### glWidget.py
This file contains the class and related functions for visualizing the elastic fusion of the doublet. For visualization, we used OpenGL contents in the "pyQTGraph" package. This allows us to use QT GUI in our project. The GUI includes parameters such as Young modulus, Poison Ratio, volume constant, density constant, elastic limit and plastic limit. Changing this parameter will change the behaviour of the doublet in real time.

## How to use
Run glWidget.py and then click Run!

## Dependencies
All mathematical computations are done using "numpy". Also, we have used KDTree from the "scipy" package to find the nearest neighbour of each particle. Visualization depends on "OpenGL" and "pyQTGraph". So, the required packages for this project are:
- numpy
- scipy
- pyQTGraph
- pyOpenGL

## Simulation Preview
Here is a preview of the simulation

![Simulation Preview](./assets/shortened_example.gif)
 
