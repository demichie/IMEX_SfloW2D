# IMEX_SfloW2D
IMEX_SfloW is a FORTRAN90 code designed to solve an hyperbolic system of partial differential equations with relaxation and source terms. 

The model is discretized in time with an explicit-implicit Runge-Kutta method where the hyperbolic part is solved explicetely and the other terms (relaxation and surce) are treated implicitely.

The finite volume solver for the hyperbolic part of the system is based on a semidiscrete central scheme and it is not tied on the specific eigenstructure of the model.

The implicit part is solved with a Newton-Raphson method where the elements of the Jacobian of the nonlinear system are evaluated numerically with a complex step derivative technique.

### Compilation

To compile:

./configure
make
make install

The executable is copied in the bin folder.

To test the code copy the executable from the bin folder in the examples folder and:

./IMEX_SfloW2d

The code also need a DEM file in esri ascii format, named topography_dem.asc

### Webpage

https://demichie.github.io/IMEX_SfloW2D/

