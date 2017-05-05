!*********************************************************************
!> \brief Parameters
!
!> This module contains the parameters for numerical solution of the
!> model.
!*********************************************************************
MODULE parameters_2d

  IMPLICIT NONE

  REAL*8 :: eps_newton        !< threshold for the convergence of the
                              !< Newton's method 
  REAL*8 :: max_dt            !< Largest time step allowed
  REAL*8 :: cfl               !< Courant-Friedrichs-Lewy parameter 

  REAL*8 :: eps_sing               !< parameter for desingularization

  CHARACTER(LEN=20) :: reconstr_variables
  
  REAL*8 :: reconstr_coeff    !< Slope coefficient in the linear reconstruction

  !> Flag to add the relaxation terms after the linear reconstruction:\n
  !> - T      => evaluate the relaxation terms
  !> - F      => reconstruction without the relaxation 
  !> .
  LOGICAL :: interfaces_relaxation

  !> Flag to choose in which way we upload the topography
  !> - T      => through a function
  !> - F      => through points
  !> .
  LOGICAL :: topography_function_flag

  !> Flag for uploading topography from a different file (topography_dem.asc)
  !> - T      => from dem
  !> - F      => from points in file.inp
  !> .
  LOGICAL :: topography_demfile

  !> Flag to choose the sort of problem to solve
  !> - T      => riemann problem
  !> - F      => generic initial conditions (uploaded through functions, to be defined in inpout_2d.f90)
  !> .
  LOGICAL :: riemann_flag

  !> Flag to choose if we add the rheology
  !> - T      => rheology activated
  !> - F      => no rheology
  !> .
  LOGICAL :: rheology_flag
  
  !> choice of the rheology model
  !> - 1      => Voellmy-Salm rheology
  !> - 2      => plastic rheology
  !> .
  INTEGER :: rheology_model

  !> Flag to choose if we model temperature transport
  !> - T      => Solve for transport equation for the temperature
  !> - F      => Do no solve for transport equation for the temperature 
  !> .
  LOGICAL :: temperature_flag

  !> Flag to choose if there is a source of mass within the domain
  !> - T      => source term added to the system
  !> - F      => source term not added to the system 
  !> .
  LOGICAL :: source_flag
 
  REAL*8 :: x_source
  REAL*8 :: y_source
  REAL*8 :: r_source
  REAL*8 :: vfr_source
  REAL*8 :: vel_source
  REAL*8 :: T_source
  
  !> Initial volume of the flow
  REAL*8 :: released_volume

  !> Initial x-coordiante of the pile
  REAL*8 :: x_release

  !> Initial y-coordinate of the pile
  REAL*8 :: y_release
  
    !> Initial velocity module of the pile
  REAL*8 :: velocity_mod_release

  !> Initial velocity direction (angle in degree):\n
  !> - >=0    => departing from positive x-axis
  !> - <0     => departign from maximum slope direction
  !.
  
  REAL*8 :: velocity_ang_release

  !> Initial temperature of the pile of material
  REAL*8 :: T_init

  !> Ambient temperature
  REAL*8 :: T_ambient

  INTEGER :: n_vars        !< Number of conservative variables
  INTEGER :: n_eqns   !< Number of equations

  INTEGER :: n_nh     !< Number of non-hyperbolic terms

  INTEGER :: n_RK     !< Runge-Kutta order
  
  INTEGER, PARAMETER :: max_nl_iter = 100

  REAL*8, PARAMETER :: tol_abs = 1.D-5
  REAL*8, PARAMETER :: tol_rel = 1.D-5

  !> Limiter for the slope in the linear reconstruction:\n
  !> - 'none'     => no limiter (constant value);
  !> - 'minmod'   => minmod sloe;
  !> - 'superbee' => superbee limiter (Roe, 1985);
  !> - 'van_leer' => monotonized central-difference limiter (van Leer, 1977)
  !> .
  INTEGER, ALLOCATABLE :: limiter(:)

  !> Finite volume method:\n
  !> - 'LxF'       => lax-friedrichs scheme;
  !> - 'GFORCE '   => gforce scheme;
  !> - 'KT'        => Kurganov and Tadmor semidiscrete scheme;
  !> .
  CHARACTER(LEN=20) :: solver_scheme     

  REAL*8 :: theta             !< Van Leer limiter parameter
  REAL*8 :: t_start           !< initial time for the run
  REAL*8 :: t_end             !< end time for the run
  REAL*8 :: t_output          !< time of the next output
  REAL*8 :: dt_output         !< time interval for the output of the solution

  INTEGER :: verbose_level

  TYPE bc
     INTEGER :: flag
     REAL*8 :: value
  END TYPE bc

  ! -------boundary conditions variables

  !> bcW&flag defines the west boundary condition:\n
  !> - bcW%flag = 0     => Dirichlet boundary condition;
  !> - bcW%flag = 1     => Neumann boundary condition.
  !> .
  !> bcLWvalue is the value of the left boundary condition:\n
  !> - value of the variable for Dirichlet boundary condition (bcW%flag=0);
  !> - gradient of the variable for Neumann boundary condition (bcW%flag=1).
  !> .
  TYPE(bc), ALLOCATABLE :: bcW(:)

  !> bcE&flag defines the east boundary condition:\n
  !> - bcE%flag = 0     => Dirichlet boundary condition;
  !> - bcE%flag = 1     => Neumann boundary condition.
  !> .
  !> bcE%value is the value of the right boundary condition:\n
  !> - value of the variable for Dirichlet boundary condition (bcE%flag=0);
  !> - gradient of the variable for Neumann boundary condition (bcE%flag=1).
  !> .
  TYPE(bc), ALLOCATABLE :: bcE(:)

  !> bcS&flag defines the south boundary condition:\n
  !> - bcS%flag = 0     => Dirichlet boundary condition;
  !> - bcS%flag = 1     => Neumann boundary condition.
  !> .
  !> bcS%value is the value of the bottom boundary condition:\n
  !> - value of the variable for Dirichlet boundary condition (bcS%flag=0);
  !> - gradient of the variable for Neumann boundary condition (bcS%flag=1).
  !> .
  TYPE(bc), ALLOCATABLE :: bcS(:)

  !> bcN&flag defines the north boundary condition:\n
  !> - bcN%flag = 0     => Dirichlet boundary condition;
  !> - bcN%flag = 1     => Neumann boundary condition.
  !> .
  !> bcN%value is the value of the top boundary condition:\n
  !> - value of the variable for Dirichlet boundary condition (bcN%flag=0);
  !> - gradient of the variable for Neumann boundary condition (bcN%flag=1).
  !> .
  TYPE(bc), ALLOCATABLE :: bcN(:)

END MODULE parameters_2d
