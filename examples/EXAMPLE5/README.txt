This example simulate a 1D pyroclastic avalanche with a Voellmy-Salm rheology over a constant slope topography. 

A Python script is provided to create the input file for this example. 
Please provide three arguments:

1) Number of cells

2) Topography Slope (in degrees)

3) Pile Slope (in degrees, relative to the topography)

4) Variables to reconstruct: phys or cons

5) Order of the RK scheme

Usage example of the script:

>> ./create_example5.py 400 13 20 phys 4

Once the input file (IMEX_SfloW2D.inp) is created create a simbolic link of the executable 
in this folder:

>> ln -s ../../../bin/IMEX_SfloW2D .

Finally, launch the solver:

>> ./IMEX_SfloW2D

A Python script to plot the results is provided. With this script you can choose the output and the variable to plot (h,hB,B,u,v)
Usage example:

>> ./plot_phys.py example5_0040.p_2d B hB

It is possible to plot the velocity u to check if the flow stops

>> ./plot_phys.py example5_0040.p_2d u


A Python script to create an animation (mp4) of the simulation is also provided. This script plot the topography and animate the flow over it. The interval between frames in milliseconds has to be given as input.
Usage example:

>> ./plot_animated.py example5 100



