# IMEX_SfloW2D
IMEX_SfloW is a FORTRAN90 code designed to model shallow water granular flow over digital elevation models (DEM) of natural terrain, with the friction forces described by the Voellmy-Salm rheology. The system is described by an hyperbolic system of partial differential equations with relaxation and source terms.

The model is discretized in time with an explicit-implicit Runge-Kutta method where the hyperbolic part is solved explicetely and the other terms (relaxation and surce) are treated implicitely to allow for larger time steps and to enforce the stopping condition.

The finite volume solver for the hyperbolic part of the system is based on the Kurganov and Petrova 2007 semidiscrete central scheme and it is not tied on the specific eigenstructure of the model.

The implicit part is solved with a Newton-Raphson method where the elements of the Jacobian of the nonlinear system are evaluated numerically with a complex step derivative technique. This automatic procedure allows for easy changes of the friction term.


### Compilation

To compile:

./configure

make

make install



The executable is copied in the bin folder.

To test the code copy the executable from the bin folder to the examples folder and:

./IMEX_SfloW2D

The code also need a DEM file in esri ascii format, named topography_dem.asc

### Webpage

https://demichie.github.io/IMEX_SfloW2D/


### Documentation

A wiki page describing the model is available at:

[https://github.com/demichie/IMEX_SfloW2D/wiki](https://github.com/demichie/IMEX_SfloW2D/wiki) 

Doxygen generated documentation of the code can be found at:

[http://demichie.github.io/IMEX_SfloW2D/html/](http://demichie.github.io/IMEX_SfloW2D/html/) 


### Authors

Mattia de' Michieli Vitturi, Giacomo Lari
