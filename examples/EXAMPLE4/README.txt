Simulation example of an avalanche of finite granular mass sliding down an inclined plane and merging continuously into a horizontal plane is presented. The initial conditions and the topography of this tutorial are the same as in Example 4.1 from Wang et al., 2004 [1]. A hemispherical shell holding the material together is suddenly released so that the bulk material commences to slide on an inclined flat plane at 35° into a horizontal run-out plane connected by a smooth transition.

Please note that we do not use the same rheological model as in the original paper of the example, but a Voellmy-Salm rhology with mu=0.3 and xi=300.

A Python script (create_example4.py) is provided to create the input file for this example. 
Please provide three arguments:

1) Number of cells in the x-direction 

2) Variables to reconstruct: phys or cons

3) Steps of the IMEX-RK scheme

Usage example of the script:

>> ./create_example4.py 150 phys 2

Once the input file (IMEX_SfloW2D.inp) is created create a simbolic link of the executable in this folder:

>> ln -s ../../bin/IMEX_SfloW2D .

Finally, launch the solver:

>> ./IMEX_SfloW2D

A Python script to plot the results is provided. With this script you can choose the output to plot.
Usage example:

>> ./plot_p_2d.py example4_0075.p_2d

A Python script to create an animation of the simulation is also provided. This script plot the topography and animate the flow over it. An animated gif is saved by the script. The script requires two argument: run name (as given in IMEX_SfloW2D.inp) and duration of the animated gif in seconds

Usage example:

>> ./animate_asc.py example4 10


REFERENCES

[1] Wang, Y., K. Hutter and S. P. Pudasaini, The Savage-Hutter theory: A system of partial differential equations for avalanche flows of snow, debris, and mud, ZAMM · Z. Angew. Math. Mech. 84, No. 8, 507 – 527 / DOI 10.1002/zamm.200310123, 2004.

