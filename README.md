## Welcome to IMEX_SfloW2D


IMEX_SfloW is a FORTRAN90 code designed to solve the shallow water equations, written as a hyperbolic system of partial differential equations with relaxation and source terms. The model is discretized in time with an explicit-implicit Runge-Kutta method where the hyperbolic part is solved explicetely and the other terms are treated implicitely. The finite volume solver for the hyperbolic part of the system is based on a semidiscrete central scheme and it is not tied on the specific eigenstructure of the model. The implicit part is solved with a Newton-Raphson method where the elements of the Jacobian of the nonlinear system are evaluated numerically with a complex step derivative technique.

### Authors and Contributors

Mattia de' Michieli Vitturi (@demichie)

Giacomo Lari

### Installation and execution

Check first if you have the LAPACK library installed on your system.

Download the IMEX_SfloW package and create the executable with the following commands from a terminal:

>./configure
>
>make
>
>make install

This will create the executable and copy it in the bin folder. You can test the executable copying it in the EXAMPLES folder and running it.

### Documentation

A wiki page describing the model is available at:

[https://github.com/demichie/IMEX_SfloW2D/wiki](https://github.com/demichie/IMEX_SfloW2D/wiki) 

Doxygen generated documentation of the code can be found at:

[http://demichie.github.io/IMEX_SfloW2D/html/](http://demichie.github.io/IMEX_SfloW2D/html/) 


### References

M. de’ Michieli Vitturi, 2016.
