
AUTOMAKE_OPTIONS = foreign

###############################################################################
# THE PROGRAMS TO BUILD
###############################################################################

# the program to build (the names of the final binaries)
bin_PROGRAMS = IMEX_SfloW2D
bindir=@bindir@

# list of sources for the binary
IMEX_SfloW2D_SOURCES = \
                complexify.f90 \
                parameters_2d.f90 \
                geometry_2d.f90 \
                constitutive_2d.f90 \
                solver_2d.f90 \
                init_2d.f90 \
                inpout_2d.f90 \
                IMEX_Sflow_2d.f90

FCFLAGS         = -Wunused -O3 -funroll-all-loops
FCDEBUGFLAGS    = -Wall -g

LIBS            = -lm -llapack

FC              = gfortran
MPFC            = gfortran
LINKER          = gfortran


