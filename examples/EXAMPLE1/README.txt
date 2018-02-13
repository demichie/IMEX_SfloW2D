This example is a 1D test for a Riemann problem with a discontinuous topography, as presented in Andrianov, 2004 [1] and Kurganov and Petrova 2007 [2].
No friction is considered here. The exact solution consists of a left-going rarefaction wave and a right-going shock-wave.

A Python script is provided to create the input file for this example. 
Please provide three arguments:

1) Number of cells

2) Variables to reconstruct: phys or cons

3) Steps of the RK scheme

Usage example of the script:

>> ./create_example1.py 400 phys 3

Once the input file (IMEX_SfloW2D.inp) is created create a simbolic link of the executable 
in this folder:

>> ln -s ../../bin/IMEX_SfloW2D .

Finally, launch the solver:

>> ./IMEX_SfloW2D

A Python script to plot the results is provided. With this script you can choose the output and the variable to plot (h,hB,B,u,v)
Usage example:

>> ./plot_phys.py example1_0020.p_2d B hB

A Python script to create an animation (mp4,gif) of the simulation is also provided. This script plots the topography and animate the flow over it. The duration of the animation in seconds has to be given as input.
Usage example:

>> ./plot_animated.py example1 5


REFERENCES

[1] Andiranov, N. testing numerical schemes for the shallow water equations, preprint available at http://www.ian.math.uni-magdeburg.de/home/andrianov/CONSTRUCT/testing.ps.gz, 2004

[2] Kurganov, A. & Petrova, G. A second-order well-balanced positivity preserving central-upwind scheme for the Saint-Venant system, Communications in Mathematical Sciences, International Press of Boston, 2007, 5, 133-160


