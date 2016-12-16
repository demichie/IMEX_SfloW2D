AC_PREREQ([2.69])

# initial information about the project
AC_INIT(IMEX-SfloW2D, 0.9,mattia.demichielivitturi@ingv.it)
AM_INIT_AUTOMAKE(0.9)
# check if the source folder is available
AC_CONFIG_SRCDIR(src/IMEX_Sflow_2d.f90)

# check for fortran compiler
AC_PROG_FC

# files to generate via autotools (prepare .am or .in source files)
AC_CONFIG_FILES([Makefile])
AC_CONFIG_FILES([src/Makefile])


AC_PREFIX_DEFAULT([..])

# finally this generates the Makefiles etc. for the build
AC_OUTPUT
