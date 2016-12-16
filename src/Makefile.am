AUTOMAKE_OPTIONS = foreign

###############################################################################
# THE PROGRAMS TO BUILD
###############################################################################

# the program to build (the names of the final binaries)
bin_PROGRAMS = TEST
bindir=@bindir@

# list of sources for the binary
PLUMEMoM_SOURCES = \
		test.f90

FCFLAGS         = -g -O3 -ffast-math
FCDEBUGFLAGS    = -Wall -g

LIBS            = -lz -lm -llapack

FC              = gfortran
MPFC            = gfortran
LINKER          = gfortran

