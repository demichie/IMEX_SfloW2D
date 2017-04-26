!********************************************************************************
!> \brief Input/Output module
!
!> This module contains all the input/output subroutine and the 
!> realted variables.
!
!> \date 07/10/2016
!> @author 
!> Mattia de' Michieli Vitturi
!
!********************************************************************************

MODULE inpout_2d

  ! -- Variables for the namelist RUN_PARAMETERS
  USE parameters_2d, ONLY : t_start , t_end , dt_output ,                       &
       topography_function_flag, topography_demfile

  USE solver_2d, ONLY : verbose_level

  ! -- Variables for the namelist NEWRUN_PARAMETERS
  USE geometry_2d, ONLY : x0 , y0 , comp_cells_x , comp_cells_y , cell_size
  USE geometry_2d, ONLY : topography_profile , n_topography_profile_x ,         &
       n_topography_profile_y
  USE init_2d, ONLY : riemann_interface
  USE parameters_2d, ONLY : temperature_flag , riemann_flag , rheology_flag
  USE parameters_2d, ONLY : source_flag

  ! -- Variables for the namelist INITIAL_CONDITIONS
  USE parameters_2d, ONLY : released_volume , x_release , y_release
  USE parameters_2d, ONLY : velocity_mod_release , velocity_ang_release
  USE parameters_2d, ONLY : T_init , T_ambient

  ! -- Variables for the namelist LEFT_STATE
  USE init_2d, ONLY : hB_W , u_W , v_W , T_W

  ! -- Variables for the namelist RIGHT_STATE
  USE init_2d, ONLY : hB_E , u_E , v_E , T_E

  ! -- Variables for the namelists LEFT/RIGHT_BOUNDARY_CONDITIONS
  USE parameters_2d, ONLY : bc

  ! -- Variables for the namelist NUMERIC_PARAMETERS
  USE parameters_2d, ONLY : solver_scheme, max_dt , cfl, limiter , theta,       &
       reconstr_coeff , interfaces_relaxation , n_RK , min_dt   

  ! -- Variables for the namelist EXPL_TERMS_PARAMETERS
  USE constitutive_2d, ONLY : grav
  USE parameters_2d, ONLY : x_source , y_source , r_source , vel_source ,       &
       T_source


  
  ! -- Variables for the namelist TEMPERATURE_PARAMETERS
  USE constitutive_2d, ONLY : rho , emissivity , exp_area_fract , enne , emme , &
       atm_heat_transf_coeff , thermal_conductivity , T_env , T_ref , T_ground ,&
       visc_par , mu_ref , c_p
    
  ! -- Variables for the namelist RHEOLOGY_PARAMETERS
  USE parameters_2d, ONLY : rheology_model
  USE constitutive_2d, ONLY : mu , xi , tau , mu_ref , visc_par , T_ref , rho

  IMPLICIT NONE

  CHARACTER(LEN=40) :: run_name           !< Name of the run
  CHARACTER(LEN=40) :: bak_name           !< Backup file for the parameters
  CHARACTER(LEN=40) :: input_file         !< File with the run parameters
  CHARACTER(LEN=40) :: output_file        !< Name of the output files
  CHARACTER(LEN=40) :: restart_file       !< Name of the restart file 
  CHARACTER(LEN=40) :: dakota_file        !< Name of the dakota file 
  CHARACTER(LEN=40) :: output_file_2d     !< Name of the output files
  CHARACTER(LEN=40) :: output_esri_file   !< Name of the esri output files

  INTEGER, PARAMETER :: input_unit = 7    !< Input data unit
  INTEGER, PARAMETER :: backup_unit = 8   !< Backup input data unit
  INTEGER, PARAMETER :: output_unit = 9   !< Output data unit
  INTEGER, PARAMETER :: restart_unit = 10 !< Restart data unit
  INTEGER, PARAMETER :: dakota_unit = 11  !< Dakota data unit
  INTEGER, PARAMETER :: output_unit_2d = 12  !< Output data 2D unit
  INTEGER, PARAMETER :: output_esri_unit = 13  !< Esri Output unit
  INTEGER, PARAMETER :: dem_esri_unit = 13  !< Computational grid Esri fmt unit

  !> Counter for the output files
  INTEGER :: output_idx 

  !> Flag to start a run from a previous output:\n
  !> - T     => Restart from a previous output (.asc or .q_2d)
  !> - F     => Restart from initial condition read from two_phases.inp
  !> .
  LOGICAL :: restart

  !> Flag to save the output in esri ascii format *.asc
  !> - T     => write esri file
  !> - F     => do not write esri file
  !> .
  LOGICAL :: output_esri_flag

  !> Flag to save the physical variables on file *.p_2d
  !> - T     => write physical variables on file
  !> - F     => do not write the physical variables
  !> .
  LOGICAL :: output_phys_flag

  !> Flag to save the conservative variables on file *.q_2d
  !> - T     => write conservative variables on file
  !> - F     => do not write the conservative variables
  !> .
  LOGICAL :: output_cons_flag

  ! -- Variables for the namelists WEST_BOUNDARY_CONDITIONS
  TYPE(bc) :: hB_bcW , u_bcW , v_bcW , T_bcW

  ! -- Variables for the namelists EAST_BOUNDARY_CONDITIONS
  TYPE(bc) :: hB_bcE , u_bcE , v_bcE , T_bcE

  ! -- Variables for the namelists SOUTH_BOUNDARY_CONDITIONS
  TYPE(bc) :: hB_bcS , u_bcS , v_bcS , T_bcS

  ! -- Variables for the namelists NORTH_BOUNDARY_CONDITIONS
  TYPE(bc) :: hB_bcN , u_bcN , v_bcN , T_bcN


  ! parameters to read a dem file
  INTEGER :: ncols, nrows, nodata_value

  REAL*8 :: xllcorner, yllcorner, cellsize


  LOGICAL :: write_first_q

  NAMELIST / run_parameters / run_name , restart , topography_demfile ,         &
       t_start , t_end , dt_output , output_cons_flag , output_esri_flag ,      &
       output_phys_flag , verbose_level

  NAMELIST / restart_parameters / restart_file , T_init , T_ambient

  NAMELIST / newrun_parameters / x0 , y0 , comp_cells_x , comp_cells_y ,        &
       cell_size , temperature_flag , source_flag , rheology_flag , riemann_flag

  NAMELIST / initial_conditions /  released_volume , x_release , y_release ,    &
       velocity_mod_release , velocity_ang_release , T_init , T_ambient

  NAMELIST / left_state / hB_W , u_W , v_W , T_W

  NAMELIST / right_state / hB_E , u_E , v_E , T_E

  NAMELIST / west_boundary_conditions / hB_bcW , u_bcW , v_bcW , T_bcW

  NAMELIST / east_boundary_conditions / hB_bcE , u_bcE , v_bcE , T_bcE

  NAMELIST / south_boundary_conditions / hB_bcS , u_bcS , v_bcS , T_bcS

  NAMELIST / north_boundary_conditions / hB_bcN , u_bcN , v_bcN , T_bcN

  NAMELIST / numeric_parameters / solver_scheme, min_dt, max_dt , cfl, limiter, &
       theta, reconstr_coeff , interfaces_relaxation , n_RK   

  NAMELIST / expl_terms_parameters / grav , x_source , y_source , r_source ,    &
       vel_source , T_source

  NAMELIST / temperature_parameters / emissivity ,  atm_heat_transf_coeff ,     &
       thermal_conductivity , exp_area_fract , c_p , enne , emme , T_env ,      &
       T_ground

  NAMELIST / rheology_parameters / rheology_model , mu , xi , tau , mu_ref ,    &
       visc_par , T_ref , rho


CONTAINS

  !******************************************************************************
  !> \brief Initialization of the variables read from the input file
  !
  !> This subroutine initialize the input variables with default values
  !> that solve for a Riemann problem. If the input file does not exist
  !> one is created with the default values.
  !
  !> @author 
  !> Mattia de' Michieli Vitturi
  !> \date 07/10/2016
  !
  !******************************************************************************

  SUBROUTINE init_param

    USE parameters_2d , ONLY : n_vars

    IMPLICIT none

    LOGICAL :: lexist

    INTEGER :: j , k

    !-- Inizialization of the Variables for the namelist RUN_PARAMETERS
    run_name = 'default'
    restart = .FALSE.
    topography_function_flag=.FALSE.
    topography_demfile=.FALSE.
    t_start = 0.0
    t_end = 5.0d-4
    dt_output = 1.d-4
    output_cons_flag = .TRUE.
    output_esri_flag = .TRUE.
    output_phys_flag = .TRUE.
    verbose_level = 0

    !-- Inizialization of the Variables for the namelist restart parameters
    restart_file = ''
    T_init = 0.D0
    T_ambient = 0.D0

    !-- Inizialization of the Variables for the namelist newrun_parameters
    x0 = 0.D0
    comp_cells_x = 400
    y0 = 0.D0
    comp_cells_y = 1
    cell_size = 2.5D-3
    temperature_flag = .FALSE.
    source_flag = .FALSE.
    rheology_flag = .FALSE.
    riemann_flag=.TRUE.
    riemann_interface = 0.5D0

    !-- Inizialization of the Variables for the namelist left_state
    hB_W = 2.222D0
    u_W = 0.8246D0
    v_W = 0.D0
    T_W = -1.D0
    
    !-- Inizialization of the Variables for the namelist right_state
    hB_E = -1.0D0
    u_E = -1.6359D0
    u_E = 0.D0
    T_E = -1.D0

    !-- Inizialization of the Variables for the namelist west boundary conditions
    hB_bcW%flag = 1 
    hB_bcW%value = 0.d0 

    u_bcW%flag = 1 
    u_bcW%value = 0.d0 
    v_bcW%flag = 1 
    v_bcW%value = 0.d0 

    !-- Inizialization of the Variables for the namelist east boundary conditions
    hB_bcE%flag = 1 
    hB_bcE%value = 0.d0 

    u_bcE%flag = 1 
    u_bcE%value = 0.d0 
    v_bcE%flag = 1 
    v_bcE%value = 0.d0 

    !-- Inizialization of the Variables for the namelist south boundary conditions
    hB_bcS%flag = 1 
    hB_bcS%value = 0.d0 

    u_bcS%flag = 1 
    u_bcS%value = 0.d0 
    v_bcS%flag = 1 
    v_bcS%value = 0.d0 

    !-- Inizialization of the Variables for the namelist north boundary conditions
    hB_bcN%flag = 1 
    hB_bcN%value = 0.d0 

    u_bcN%flag = 1 
    u_bcN%value = 0.d0 
    v_bcN%flag = 1 
    v_bcN%value = 0.d0 

    !-- Inizialization of the Variables for the namelist NUMERIC_PARAMETERS
    min_dt = 1.d-6
    max_dt = 1.d-3
    solver_scheme = 'KT'
    n_RK = 1
    cfl = 0.24
    limiter(1:n_vars) = 0
    theta=1.0
    reconstr_coeff = 1.0

    !-- Inizialization of the Variables for the namelist EXPL_TERMS_PARAMETERS
    grav = 9.81D0
    x_source = 0.D0
    y_source = 0.D0
    r_source = 0.D0
    vel_source = 0.D0
    T_source = 0.D0

    !-- Inizialization of the Variables for the namelist TEMPERATURE_PARAMETERS
    exp_area_fract = 0.5D0
    emissivity = 0.0D0                 ! no radiation to atmosphere
    atm_heat_transf_coeff = 0.0D0      ! no convection to atmosphere
    thermal_conductivity = 0.0D0       ! no conduction to ground
    mu_ref = 0.0D0                     ! no viscous heating
    enne = 4.0D0
    emme = 12.D0
    T_env = 300.0D0
    T_ground = 1200.0D0
    c_p = 1200.D0
    
    !-- Inizialization of the Variables for the namelist RHEOLOGY_PARAMETERS
    rheology_model = 0
    mu = 0.D0
    xi = 0.D0
    tau = 0.D0
    rho = 0.D0
    T_ref = 0.D0
    visc_par = 0.0D0
    mu_ref = 0.0D0


    !-------------- Check if input file exists ----------------------------------
    input_file = 'IMEX_SfloW2D.inp'

    INQUIRE (FILE=input_file,exist=lexist)

    IF (lexist .EQV. .FALSE.) THEN

       OPEN(input_unit,FILE=input_file,STATUS='NEW')

       WRITE(input_unit, run_parameters )
       WRITE(input_unit, newrun_parameters )
       WRITE(input_unit, left_state )
       WRITE(input_unit, right_state )
       WRITE(input_unit, west_boundary_conditions )
       WRITE(input_unit, east_boundary_conditions )
       WRITE(input_unit, south_boundary_conditions )
       WRITE(input_unit, north_boundary_conditions )
       WRITE(input_unit, numeric_parameters )
       WRITE(input_unit, expl_terms_parameters )

       n_topography_profile_x = 2
       n_topography_profile_y = 2

       ALLOCATE( topography_profile( 2 , n_topography_profile_x ,               &
            n_topography_profile_y) )

       topography_profile(1,1,1) = 0.d0
       topography_profile(2,1,1) = 0.d0
       topography_profile(3,1,1) = 1.d0

       topography_profile(1,1,2) = 0.d0
       topography_profile(2,1,2) = 1.d0
       topography_profile(3,1,2) = 1.d0

       topography_profile(1,2,1) = 1.d0
       topography_profile(2,2,1) = 0.d0
       topography_profile(3,2,1) = 0.5d0

       topography_profile(1,2,2) = 1.d0
       topography_profile(2,2,2) = 1.d0
       topography_profile(3,2,2) = 0.5d0



       WRITE(input_unit,*) '''TOPOGRAPHY_PROFILE'''
       WRITE(input_unit,*) n_topography_profile_x
       WRITE(input_unit,*) n_topography_profile_y

       DO j = 1, n_topography_profile_x

          DO k = 1, n_topography_profile_y

            WRITE(input_unit,108) topography_profile(1:3,j,k)

108         FORMAT(3(1x,e14.7))

          ENDDO

       END DO

       CLOSE(input_unit)

       WRITE(*,*) 'Input file IMEX_SfloW2D.inp not found'
       WRITE(*,*) 'A new one with default values has been created'
       STOP

    ELSE

    END IF

    ! output file index
    output_idx = 0

  END SUBROUTINE init_param

  !******************************************************************************
  !> \brief Read the input file
  !
  !> This subroutine read the input parameters from the file 
  !> "two_phases.inp" and write a backup file of the input parameters 
  !> with name "run_name.bak", where run_name is read from the input
  !> file.
  !
  !> \date 07/10/2016
  !> @author 
  !> Mattia de' Michieli Vitturi
  !
  !******************************************************************************

  SUBROUTINE read_param

    USE parameters_2d, ONLY : n_vars , n_eqns
    USE parameters_2d, ONLY : limiter
    USE parameters_2d, ONLY : bcW , bcE , bcS , bcN

    IMPLICIT none

    REAL*8 :: max_cfl

    LOGICAL :: tend1 
    CHARACTER(LEN=80) :: card

    INTEGER :: j,k

    INTEGER :: dot_idx
    
    CHARACTER(LEN=3) :: check_file

    LOGICAL :: lexist

    CHARACTER(LEN=15) :: chara

    INTEGER :: ios

    OPEN(input_unit,FILE=input_file,STATUS='old')


    ! ------- READ run_parameters NAMELIST -----------------------------------
    READ(input_unit, run_parameters,IOSTAT=ios )

    IF ( ios .NE. 0 ) THEN

       WRITE(*,*) 'IOSTAT=',ios
       WRITE(*,*) 'ERROR: problem with namelist RUN_PARAMETERS'
       WRITE(*,*) 'Please check the input file'
       STOP

    ELSE

       REWIND(input_unit)

    END IF

    ! ------- READ newrun_parameters NAMELIST --------------------------------
    READ(input_unit,newrun_parameters,IOSTAT=ios)

    IF ( ios .NE. 0 ) THEN

       WRITE(*,*) 'IOSTAT=',ios
       WRITE(*,*) 'ERROR: problem with namelist NEWRUN_PARAMETERS'
       WRITE(*,*) 'Please check the input file'
       STOP

    ELSE

       REWIND(input_unit)

    END IF

    IF ( ( comp_cells_x .EQ. 1 ) .OR. ( comp_cells_y .EQ. 1 ) ) THEN

       WRITE(*,*) '----- 1D SIMULATION -----' 

    ELSE
       
       WRITE(*,*) '----- 2D SIMULATION -----' 

    END IF
       
       
    IF ( temperature_flag ) THEN

       n_vars = 4
       n_eqns = 4

       ELSE

       n_vars = 3
       n_eqns = 3

    END IF
   
    ALLOCATE( limiter(n_vars) )
    ALLOCATE( bcW(n_vars) , bcE(n_vars) , bcS(n_vars) , bcN(n_vars) , )

    IF ( restart ) THEN

       READ(input_unit,restart_parameters,IOSTAT=ios)

       IF ( ios .NE. 0 ) THEN
          
          WRITE(*,*) 'IOSTAT=',ios
          WRITE(*,*) 'ERROR: problem with namelist RESTART_PARAMETERS'
          WRITE(*,*) 'Please check the input file'
          STOP
          
       ELSE
          
          REWIND(input_unit)
          
       END IF

       IF ( temperature_flag .AND. restart ) THEN
          
          dot_idx = SCAN(restart_file, ".", .TRUE.)
          
          check_file = restart_file(dot_idx+1:dot_idx+3)
          
          IF ( ( check_file .EQ. 'asc' )  .AND.                                 &
               ( T_init*T_ambient .EQ. 0.D0 ) ) THEN
          
             WRITE(*,*) 'T_init=',T_init
             WRITE(*,*) 'T_ambient=',T_ambient
             WRITE(*,*) 'Add the variables to the namelist RESTART_PARAMETERS'
             STOP
             
          END IF

       END IF

    ELSE

       IF ( riemann_flag ) THEN

          READ(input_unit,left_state,IOSTAT=ios)

          IF ( ios .NE. 0 ) THEN
             
             WRITE(*,*) 'IOSTAT=',ios
             WRITE(*,*) 'ERROR: problem with namelist LEFT_PARAMETERS'
             WRITE(*,*) 'Please check the input file'
             STOP
             
          ELSE
             
             REWIND(input_unit)

             IF ( ( temperature_flag ) .AND. ( T_W .EQ. -1.D0 ) ) THEN

                WRITE(*,*) 'ERROR: problem with namelist LEFT_PARAMETERS'
                WRITE(*,*) 'Initial temperature T_W not defined'
                STOP
                
             END IF
             
          END IF
          
          READ(input_unit,right_state,IOSTAT=ios)

          IF ( ios .NE. 0 ) THEN
             
             WRITE(*,*) 'IOSTAT=',ios
             WRITE(*,*) 'ERROR: problem with namelist RIGHT_PARAMETERS'
             WRITE(*,*) 'Please check the input file'
             STOP
             
          ELSE
             
             REWIND(input_unit)

             IF ( ( temperature_flag ) .AND. ( T_E .EQ. -1.D0 ) ) THEN

                WRITE(*,*) 'ERROR: problem with namelist RIGHT_PARAMETERS'
                WRITE(*,*) 'Initial temperature T_E not defined'
                STOP
                
             END IF
             
             
          END IF
            
       ELSE

          READ(input_unit,initial_conditions,IOSTAT=ios)

          IF ( ios .NE. 0 ) THEN
             
             WRITE(*,*) 'IOSTAT=',ios
             WRITE(*,*) 'ERROR: problem with namelist INITIAL_CONDITIONS'
             WRITE(*,*) 'Please check the input file'
             STOP
             
          ELSE
             
             REWIND(input_unit)
 
          END IF
            
       IF ( ( temperature_flag ) .AND. ( T_init*T_ambient .EQ. 0.D0 ) ) THEN

          WRITE(*,*) 'T_init=',T_init
          WRITE(*,*) 'T_ambient=',T_ambient
          WRITE(*,*) 'Add the two variables to the namelist INITIAL_CONDITIONS'
          STOP

       END IF

     
       END IF

    END IF


    ! ------- READ boundary_conditions NAMELISTS --------------------------------

    IF ( COMP_CELLS_X .GT. 1 ) THEN
    
       ! West boundary conditions
       READ(input_unit,west_boundary_conditions,IOSTAT=ios)
       
       IF ( ios .NE. 0 ) THEN
          
          WRITE(*,*) 'IOSTAT=',ios
          WRITE(*,*) 'ERROR: problem with namelist WEST_BOUNDARY_CONDITIONS'
          WRITE(*,*) 'Please check the input file'
          STOP
          
       ELSE
          
          REWIND(input_unit)
          
       END IF
       
       bcW(1) = hB_bcW 
       bcW(2) = u_bcW 
       bcW(3) = v_bcW 
       
       ! East boundary conditions
       READ(input_unit,east_boundary_conditions,IOSTAT=ios)
       
       IF ( ios .NE. 0 ) THEN
          
          WRITE(*,*) 'IOSTAT=',ios
          WRITE(*,*) 'ERROR: problem with namelist EAST_BOUNDARY_CONDITIONS'
          WRITE(*,*) 'Please check the input file'
          STOP
          
       ELSE
          
          REWIND(input_unit)
          
       END IF
       
       bcE(1) = hB_bcE 
       bcE(2) = u_bcE 
       bcE(3) = v_bcE 

    END IF

    IF ( comp_cells_y .GT. 1 ) THEN
    
       ! South boundary conditions
       READ(input_unit,south_boundary_conditions,IOSTAT=ios)
       
       IF ( ios .NE. 0 ) THEN
          
          WRITE(*,*) 'IOSTAT=',ios
          WRITE(*,*) 'ERROR: problem with namelist SOUTH_BOUNDARY_CONDITIONS'
          WRITE(*,*) 'Please check the input file'
          STOP
          
       ELSE
          
          REWIND(input_unit)
          
       END IF
       
       bcS(1) = hB_bcS 
       bcS(2) = u_bcS 
       bcS(3) = v_bcS 
       
       ! North boundary conditions
       READ(input_unit,north_boundary_conditions,IOSTAT=ios)
       
       IF ( ios .NE. 0 ) THEN
          
          WRITE(*,*) 'IOSTAT=',ios
          WRITE(*,*) 'ERROR: problem with namelist NORTH_BOUNDARY_CONDITIONS'
          WRITE(*,*) 'Please check the input file'
          STOP
          
       ELSE
          
          REWIND(input_unit)
          
       END IF
       
       bcN(1) = hB_bcN 
       bcN(2) = u_bcN 
       bcN(3) = v_bcN 

    END IF
    
    IF ( temperature_flag ) THEN
       
       bcW(4) = T_bcW
       bcE(4) = T_bcE
       bcS(4) = T_bcS
       bcN(4) = T_bcN
       
    END IF
    
    ! ------- READ numeric_parameters NAMELIST ----------------------------------

    READ(input_unit,numeric_parameters)

    IF ( ios .NE. 0 ) THEN
       
       WRITE(*,*) 'IOSTAT=',ios
       WRITE(*,*) 'ERROR: problem with namelist NUMERIC_PARAMETERS'
       WRITE(*,*) 'Please check the input file'
       STOP
       
    ELSE
       
       REWIND(input_unit)
       
    END IF

    IF ( ( solver_scheme .NE. 'LxF' ) .AND. ( solver_scheme .NE. 'KT' ) .AND.   &
         ( solver_scheme .NE. 'GFORCE' ) ) THEN

       WRITE(*,*) 'WARNING: no correct solver scheme selected',solver_scheme
       WRITE(*,*) 'Chose between: LxF, GFORCE or KT'
       STOP

    END IF

    IF  ( ( solver_scheme.EQ.'LxF' ) .OR. ( solver_scheme.EQ.'GFORCE' ) ) THEN 

       max_cfl = 1.0

    ELSEIF ( solver_scheme .EQ. 'KT' ) THEN

       IF ( ( comp_cells_x .EQ. 1 ) .OR. ( comp_cells_y .EQ. 1 ) ) THEN

          max_cfl = 0.50D0

       ELSE

          max_cfl = 0.25D0

       END IF

    END IF


    IF ( ( cfl .GT. max_cfl ) .OR. ( cfl .LT. 0.D0 ) ) THEN

       WRITE(*,*) 'WARNING: wrong value of cfl ',cfl
       WRITE(*,*) 'Choose a value between 0.0 and ',max_cfl
       READ(*,*)

    END IF

    IF ( verbose_level .GE. 1 ) WRITE(*,*) 'Limiters',limiter

    IF ( ( MAXVAL(limiter) .GT. 3 ) .OR. ( MINVAL(limiter) .LT. 0 ) ) THEN

       WRITE(*,*) 'WARNING: wrong limiter ',limiter
       WRITE(*,*) 'Choose among: none, minmod,superbee,van_leer'
       STOP         

    END IF

    IF ( ( reconstr_coeff .GT. 1.0D0 ) .OR. ( reconstr_coeff .LT. 0.D0 ) ) THEN

       WRITE(*,*) 'WARNING: wrong value of reconstr_coeff ',reconstr_coeff
       WRITE(*,*) 'Change the value between 0.0 and 1.0 in the input file'
       READ(*,*)

    END IF

    ! ------- READ expl_terms_parameters NAMELIST -------------------------------

    READ(input_unit, expl_terms_parameters,IOSTAT=ios)

    IF ( ios .NE. 0 ) THEN
       
       WRITE(*,*) 'IOSTAT=',ios
       WRITE(*,*) 'ERROR: problem with namelist EXPL_TERMS_PARAMETERS'
       WRITE(*,*) 'Please check the input file'
       STOP
       
    ELSE
       
       REWIND(input_unit)
       
    END IF

    IF ( source_flag ) THEN

       IF ( r_source .EQ. 0.D0 ) THEN

          WRITE(*,*) 'ERROR: problem with namelist EXPL_TERMS_PARAMETERS'
          WRITE(*,*) 'R_SOURCE =',r_source
          WRITE(*,*) 'Please check the input file'
          STOP

       END IF

       IF ( vel_source .EQ. 0.D0 ) THEN

          WRITE(*,*) 'ERROR: problem with namelist EXPL_TERMS_PARAMETERS'
          WRITE(*,*) 'VEL_SOURCE =',vel_source
          WRITE(*,*) 'Please check the input file'
          STOP

       END IF

       IF ( temperature_flag .AND. ( r_source .EQ. 0.D0 ) ) THEN

          WRITE(*,*) 'ERROR: problem with namelist EXPL_TERMS_PARAMETERS'
          WRITE(*,*) 'TEMPERATURE_FLAG =',temperature_flag,' R_SOURCE =',r_source
          WRITE(*,*) 'Please check the input file'
          STOP

       END IF

    END IF
    
    ! ------- READ temperature_parameters NAMELIST ------------------------------

    IF ( temperature_flag ) THEN

       READ(input_unit, temperature_parameters,IOSTAT=ios)
       
       IF ( ios .NE. 0 ) THEN
          
          WRITE(*,*) 'IOSTAT=',ios
          WRITE(*,*) 'ERROR: problem with namelist TEMPERATURE_PARAMETERS'
          WRITE(*,*) 'Please check the input file'
          STOP
          
       ELSE
          
          REWIND(input_unit)
          
       END IF
       
       IF ( emissivity .EQ. 0.D0 ) THEN

          WRITE(*,*) 'No radiative term: emissivity =',emissivity

       END IF

       IF ( atm_heat_transf_coeff .EQ. 0.D0 ) THEN

          WRITE(*,*) 'No convective term: atm_heat_transf_coeff =',             &
               atm_heat_transf_coeff

       END IF
       
       IF ( thermal_conductivity .EQ. 0.D0 ) THEN

          WRITE(*,*) 'No conductive term: thermal_conductivity =',              &
               thermal_conductivity

       END IF

    END IF


    ! ------- READ rheology_parameters NAMELIST --------------------------------

    IF ( rheology_flag ) THEN

       READ(input_unit, rheology_parameters,IOSTAT=ios)
       
       IF ( ios .NE. 0 ) THEN
          
          WRITE(*,*) 'IOSTAT=',ios
          WRITE(*,*) 'ERROR: problem with namelist RHEOLOGY_PARAMETERS'
          WRITE(*,*) 'Please check the input file'
          STOP
          
       ELSE
          
          REWIND(input_unit)
          
       END IF
 
       IF ( rheology_model .EQ. 0 ) THEN

          WRITE(*,*) 'ERROR: problem with namelist RHEOLOGY_PARAMETERS'
          WRITE(*,*) 'RHEOLOGY_FLAG' , rheology_flag , 'RHEOLOGY_MODEL =' ,     &
               rheology_model
          WRITE(*,*) 'Please check the input file'
          STOP
          
       ELSEIF ( rheology_model .EQ. 1 ) THEN

          IF ( ( mu .EQ. 0.D0 ) .AND. ( xi .EQ. 0.D0 ) ) THEN

             WRITE(*,*) 'ERROR: problem with namelist RHEOLOGY_PARAMETERS'
             WRITE(*,*) 'RHEOLOGY_MODEL =' , rheology_model
             WRITE(*,*) 'MU =' , mu ,' XI =' , xi
             WRITE(*,*) 'Please check the input file'
             STOP
             
          END IF

          IF ( ( rho .NE. 0.D0 ) .OR. ( T_ref .NE. 0.D0 ) .OR.                  &
               ( mu_ref .NE. 0.D0 ) .OR. ( visc_par .NE. 0.D0 ) .OR.            &
               ( tau .NE. 0.D0 ) ) THEN

             WRITE(*,*) 'WARNING: parameters not used in RHEOLOGY_PARAMETERS'
             IF ( rho .NE. 0.D0 ) WRITE(*,*) 'rho =',rho 
             IF ( T_ref .NE. 0.D0 ) WRITE(*,*) 'T_ref =',T_ref 
             IF ( mu_ref .NE. 0.D0 ) WRITE(*,*) 'mu_ref =',mu_ref 
             IF ( visc_par .NE. 0.D0 ) WRITE(*,*) 'visc_par =',visc_par
             IF ( tau .NE. 0.D0 ) WRITE(*,*) 'tau =',tau 
             WRITE(*,*) 'Press ENTER to continue'
             READ(*,*)

          END IF

       ELSEIF ( rheology_model .EQ. 1 ) THEN

          IF ( tau .EQ. 0.D0 )  THEN

             WRITE(*,*) 'ERROR: problem with namelist RHEOLOGY_PARAMETERS'
             WRITE(*,*) 'RHEOLOGY_MODEL =' , rheology_model
             WRITE(*,*) 'TAU =' , tau
             WRITE(*,*) 'Please check the input file'
             STOP
             
          END IF
          
          IF ( ( rho .NE. 0.D0 ) .OR. ( T_ref .NE. 0.D0 ) .OR.                  &
               ( mu_ref .NE. 0.D0 ) .OR. ( visc_par .NE. 0.D0 ) .OR.            &
               ( mu .NE. 0.D0 ) .OR. ( xi .NE. 0.D0 ) ) THEN

             WRITE(*,*) 'WARNING: parameters not used in RHEOLOGY_PARAMETERS'
             IF ( rho .NE. 0.D0 ) WRITE(*,*) 'rho =',rho 
             IF ( T_ref .NE. 0.D0 ) WRITE(*,*) 'T_ref =',T_ref 
             IF ( mu_ref .NE. 0.D0 ) WRITE(*,*) 'mu_ref =',mu_ref 
             IF ( visc_par .NE. 0.D0 ) WRITE(*,*) 'visc_par =',visc_par
             IF ( mu .NE. 0.D0 ) WRITE(*,*) 'mu =',mu 
             IF ( xi .NE. 0.D0 ) WRITE(*,*) 'xi =',xi
             WRITE(*,*) 'Press ENTER to continue'
             READ(*,*)


          END IF

       ELSEIF ( rheology_model .EQ. 3 ) THEN

          IF ( mu_ref .LE. 0.D0 ) THEN
             
             WRITE(*,*) 'ERROR: problem with namelist RHEOLOGY_PARAMETERS'
             WRITE(*,*) 'MU_REF =' , mu_ref 
             WRITE(*,*) 'Please check the input file'
             STOP
             
          END IF

          IF ( T_ref .LE. 0.D0 ) THEN
             
             WRITE(*,*) 'ERROR: problem with namelist RHEOLOGY_PARAMETERS'
             WRITE(*,*) 'T_REF =' , T_ref 
             WRITE(*,*) 'Please check the input file'
             STOP
             
          END IF

          IF ( rho .LE. 0.D0 ) THEN
             
             WRITE(*,*) 'ERROR: problem with namelist RHEOLOGY_PARAMETERS'
             WRITE(*,*) 'RHO =' , rho
             WRITE(*,*) 'Please check the input file'
             STOP
             
          END IF

          IF ( visc_par .LT. 0.D0 ) THEN
             
             WRITE(*,*) 'ERROR: problem with namelist RHEOLOGY_PARAMETERS'
             WRITE(*,*) 'VISC_PAR =' , visc_par
             WRITE(*,*) 'Please check the input file'
             STOP
          
          ELSEIF ( visc_par .EQ. 0.D0 ) THEN
             
             WRITE(*,*) 'WARNING: temperature and momentum uncoupled'
             WRITE(*,*) 'VISC_PAR =' , visc_par
             WRITE(*,*) 'Press ENTER to continue'
             READ(*,*)
   
          END IF

          IF ( ( mu .NE. 0.D0 ) .OR. ( xi .NE. 0.D0 ) .OR. ( tau .NE. 0.D0 ) )  &
               THEN

             WRITE(*,*) 'WARNING: parameters not used in RHEOLOGY_PARAMETERS'
             IF ( mu .NE. 0.D0 ) WRITE(*,*) 'mu =',mu 
             IF ( xi .NE. 0.D0 ) WRITE(*,*) 'xi =',xi
             IF ( tau .NE. 0.D0 ) WRITE(*,*) 'tau =',tau 
             WRITE(*,*) 'Press ENTER to continue'
             READ(*,*)

          END IF

       ELSE

             WRITE(*,*) 'ERROR: problem with namelist RHEOLOGY_PARAMETERS'
             WRITE(*,*) 'RHEOLOGY_MODEL =' , rheology_model
             WRITE(*,*) 'Please check the input file'
             STOP

       END IF

    END IF

    ! ---------------------------------------------------------------------------

    IF ( verbose_level .GE. 1 ) WRITE(*,*) 
    
    ! read topography from .inp (recommended for simple batimetries) 
    IF ( .NOT.topography_demfile ) THEN

       tend1 = .FALSE.

       WRITE(*,*) 'Searching for topography_profile'

       topography_profile_search: DO

          READ(input_unit,*, END = 200 ) card

          IF( TRIM(card) == 'TOPOGRAPHY_PROFILE' ) THEN

             EXIT topography_profile_search

          END IF

       END DO topography_profile_search


       READ(input_unit,*) n_topography_profile_x

       IF ( verbose_level .GE. 1 ) WRITE(*,*) 'n_topography_profile_x' ,        &
            n_topography_profile_x

       READ(input_unit,*) n_topography_profile_y

       IF ( verbose_level .GE. 1 ) WRITE(*,*) 'n_topography_profile_y' ,        &
            n_topography_profile_y

       ALLOCATE( topography_profile( 3 , n_topography_profile_x ,               &
            n_topography_profile_y ) )

       DO j = 1, n_topography_profile_x

          DO k = 1, n_topography_profile_y

             READ(input_unit,*) topography_profile(1:3,j,k)

             IF ( verbose_level.GE.1 ) WRITE(*,*) j,k,topography_profile(1:3,j,k)

          ENDDO

       END DO

       GOTO 210
200    tend1 = .TRUE.
210    CONTINUE


       ! read topography from dem file (recommended for complex batimetries) 
    ELSE

       WRITE(*,*) 'Searching for DEM file'

       INQUIRE(FILE='topography_dem.asc',EXIST=lexist)

       IF(lexist)THEN

          OPEN(2001, file='topography_dem.asc', status='old', action='read')

       ELSE

          WRITE(*,*) 'no dem file'
          STOP

       ENDIF

       READ(2001,*) chara, ncols
       READ(2001,*) chara, nrows
       READ(2001,*) chara, xllcorner
       READ(2001,*) chara, yllcorner
       READ(2001,*) chara, cellsize
       READ(2001,*) chara, nodata_value

       IF ( x0 .LT. xllcorner ) THEN 

          WRITE(*,*) 'Computational domain problem'
          WRITE(*,*) 'x0 < xllcorner',x0,xllcorner
          STOP

       END IF

       IF ( x0+comp_cells_x*cell_size .GT. xllcorner+ncols*cellsize ) THEN 

          WRITE(*,*) 'Computational domain problem'
          WRITE(*,*) 'right edge > xllcorner+ncols*cellsize',                   &
               x0+comp_cells_x*cell_size , xllcorner+ncols*cellsize
          STOP

       END IF

       IF ( y0 .LT. yllcorner ) THEN 

          WRITE(*,*) 'Computational domain problem'
          WRITE(*,*) 'y0 < yllcorner',y0,yllcorner
          STOP

       END IF

       IF ( y0+comp_cells_y*cell_size .GT. yllcorner+nrows*cellsize ) THEN 

          WRITE(*,*) 'Computational domain problem'
          WRITE(*,*) 'top edge > yllcorner+nrows*cellsize',                     &
               y0+comp_cells_y*cell_size , yllcorner+nrows*cellsize
          STOP

       END IF


       WRITE(*,*) 'Reading DEM file' 
       WRITE(*,*) 'ncols',ncols
       WRITE(*,*) 'nrows',nrows

       n_topography_profile_x = ncols

       n_topography_profile_y = nrows

       ALLOCATE( topography_profile( 3 , n_topography_profile_x ,               &
            n_topography_profile_y) )

       topography_profile(1,1,:) = xllcorner

       DO j=2,n_topography_profile_x 

          topography_profile(1,j,:) = topography_profile(1,j-1,:) + cellsize

       ENDDO

       topography_profile(2,:,1) = yllcorner

       DO k=2,n_topography_profile_y

          topography_profile(2,:,k) = topography_profile(2,:,k-1) + cellsize

       ENDDO

       ! store in opposite sense

       DO k=1,n_topography_profile_y

          WRITE(*,FMT="(A1,A,t21,F6.2,A)",ADVANCE="NO") ACHAR(13),              &
               & " Percent Complete: " ,                                        &
               ( REAL(k) / REAL(n_topography_profile_y))*100.0, "%"

          READ(2001,*) topography_profile(3,:,n_topography_profile_y-k+1)

       ENDDO

       WRITE(*,*) ''


       CLOSE(2001)

    ENDIF

    CLOSE( input_unit )

    bak_name = TRIM(run_name)//'.bak'

    OPEN(backup_unit,file=bak_name,status='unknown')

    WRITE(backup_unit, run_parameters )

    IF ( restart ) THEN

       WRITE(backup_unit,restart_parameters)

    ELSE

       WRITE(backup_unit,newrun_parameters)

       IF ( riemann_flag) THEN

          WRITE(backup_unit,left_state)
          WRITE(backup_unit,right_state)

       ELSE

          WRITE(backup_unit,initial_conditions)

       END IF

    END IF

    IF ( comp_cells_x .GT. 1 ) THEN

       WRITE(backup_unit,west_boundary_conditions)
       WRITE(backup_unit,east_boundary_conditions)

    END IF

    IF ( comp_cells_y .GT. 1 ) THEN

       WRITE(backup_unit,north_boundary_conditions)
       WRITE(backup_unit,south_boundary_conditions)

    END IF

    WRITE(backup_unit, numeric_parameters )

    WRITE(backup_unit, expl_terms_parameters )

    IF ( temperature_flag ) WRITE(backup_unit,temperature_parameters)

    IF ( rheology_flag ) WRITE(backup_unit,rheology_parameters)

    IF ( .NOT.topography_demfile ) THEN

       WRITE(backup_unit,*) '''TOPOGRAPHY_PROFILE'''
       WRITE(backup_unit,*) n_topography_profile_x
       WRITE(backup_unit,*) n_topography_profile_y
       
       DO j = 1, n_topography_profile_x
          
          DO k = 1, n_topography_profile_y
             
             WRITE(backup_unit,107) topography_profile(1:3,j,k)
             
107          FORMAT(3(1x,e14.7))
             
          ENDDO
          
       END DO

    END IF

    CLOSE(backup_unit)


  END SUBROUTINE read_param


  !******************************************************************************
  !> \brief Read the solution from the restart unit
  !
  !> This subroutine is called when the parameter "restart" in the input 
  !> file is TRUE. Then the initial solution is read from a file. 
  !
  !> \date 07/10/2016
  !> @author 
  !> Mattia de' Michieli Vitturi
  !
  !******************************************************************************

  SUBROUTINE read_solution

    USE solver_2d, ONLY : allocate_solver_variables

    USE geometry_2d, ONLY : comp_cells_x , x0 , comp_cells_y , y0 , dx , dy
    USE geometry_2d, ONLY : B_cent
    USE init_2d, ONLY : thickness_init
    ! USE parameters_2d, ONLY : n_vars
    USE solver_2d, ONLY : q

    IMPLICIT none

    CHARACTER(LEN=15) :: chara

    INTEGER :: j,k

    INTEGER :: dot_idx

    LOGICAL :: lexist

    CHARACTER(LEN=30) :: string

    CHARACTER(LEN=3) :: check_file

    INTEGER :: ncols , nrows , nodata_value

    REAL*8 :: xllcorner , yllcorner , cellsize

    REAL*8 :: xj , yk

    INQUIRE (FILE=restart_file,exist=lexist)

    IF (lexist .EQV. .FALSE.) THEN

       WRITE(*,*) 'Restart: ',TRIM(restart_file),' not found'
       STOP

    ELSE

       OPEN(restart_unit,FILE=restart_file,STATUS='old')
       
       WRITE(*,*) 'Restart: ',TRIM(restart_file), ' found'

    END IF

    dot_idx = SCAN(restart_file, ".", .TRUE.)

    check_file = restart_file(dot_idx+1:dot_idx+3)

    IF ( check_file .EQ. 'asc' ) THEN
       
       READ(restart_unit,*) chara, ncols
       READ(restart_unit,*) chara, nrows
       READ(restart_unit,*) chara, xllcorner
       READ(restart_unit,*) chara, yllcorner
       READ(restart_unit,*) chara, cellsize
       READ(restart_unit,*) chara, nodata_value
       
       ALLOCATE( thickness_init(ncols,nrows) )

       IF ( ncols .NE. comp_cells_x ) THEN
          
          WRITE(*,*) 'ncols not equal to comp_cells_x',ncols,comp_cells_x
          STOP
          
       END IF
       
       IF ( nrows .NE. comp_cells_y ) THEN
          
          WRITE(*,*) 'nrows not equal to comp_cells_y',nrows,comp_cells_y
          STOP
          
       END IF
       
       IF ( DABS( xllcorner - (x0 + 0.5D0*cellsize) ) .GT. 1.D-5*cellsize ) THEN
          
          WRITE(*,*) 'xllcorner not equal to x0+0.5*cellsize', xllcorner ,      &
               x0+0.5*cellsize
          STOP
          
       END IF
       
       IF ( DABS( yllcorner - (y0 + 0.5D0*cellsize) ) .GT. 1.D-5*cellsize ) THEN
          
          WRITE(*,*) 'yllcorner not equal to y0+0.5*cellsize', yllcorner ,      &
               y0+0.5*cellsize
          STOP
          
       END IF
       
       IF ( cellsize .NE. cell_size ) THEN
          
          WRITE(*,*) 'cellsize not equal to cell_size', cellsize , cell_size
          STOP
          
       END IF
       
       DO k=1,nrows
          
          WRITE(*,FMT="(A1,A,t21,F6.2,A)",ADVANCE="NO") ACHAR(13),              &
               & " Percent Complete: ",( REAL(k) / REAL(comp_cells_y))*100.0, "%"
          
          READ(restart_unit,*) thickness_init(1:ncols,nrows-k+1)
          
       ENDDO

       WHERE ( thickness_init .EQ. nodata_value )

          thickness_init = 0.D0
          
       END WHERE
       
       WRITE(*,*)
       WRITE(*,*) 'Total volume =', cellsize**2 * SUM(thickness_init)

       q(1,:,:) = B_cent(:,:) + thickness_init(:,:)
       q(2,:,:) = 0.D0
       q(3,:,:) = 0.D0

       IF ( temperature_flag ) THEN

          q(4,:,:) = thickness_init(:,:) * T_init

          WHERE ( thickness_init .GT. 0.D0 )
          
             q(4,:,:) = thickness_init(:,:) * T_ambient
          
          END WHERE
    
       END IF

       output_idx = 0

       IF ( verbose_level .GE. 1 ) THEN

          WRITE(*,*) 'Min q(1,:,:) =',MINVAL(q(1,:,:))
          WRITE(*,*) 'Max q(1,:,:) =',MAXVAL(q(1,:,:))
          WRITE(*,*) 'SUM(q(1,:,:)) =',SUM(q(1,:,:))

          DO k=1,nrows

             WRITE(*,*) k,B_cent(:,k)
             READ(*,*)

          END DO

          WRITE(*,*) 'SUM(B_cent(:,:)) =',SUM(B_cent(:,:))
          READ(*,*)

       END IF

    ELSEIF ( check_file .EQ. 'q_2' ) THEN
    
       DO k=1,comp_cells_y
          
          DO j=1,comp_cells_x

             IF ( temperature_flag ) THEN
             
                READ(restart_unit,1004) xj , yk , q(:,j,k)

             ELSE
                
                READ(restart_unit,1003) xj , yk , q(:,j,k)

             END IF
                
             IF ( verbose_level .GE. 2 ) THEN
                
                WRITE(*,*) 'x,y,q,B',xj,yk,q(:,j,k),B_cent(j,k)
                WRITE(*,*) 'h,T',q(1,j,k)-B_cent(j,k), q(4,j,k)                 &
                     / ( q(1,j,k)-B_cent(j,k) )
                READ(*,*)
                
             END IF
             
             IF ( q(1,j,k) .LE. B_cent(j,k) ) THEN

                IF ( verbose_level .GE. 2 ) THEN

                   WRITE(*,*) 'q(1,j,k),B_cent(j,k)',j,k,q(1,j,k),B_cent(j,k)
                   READ(*,*)

                END IF

                q(1,j,k) = B_cent(j,k)

             END IF
             
          ENDDO
          
          READ(restart_unit,*)  
          
       END DO

       WRITE(*,*) 'Total volume =',dx*dy* SUM( q(1,:,:)-B_cent(:,:) )
          
1003   FORMAT(5e20.12)
1004   FORMAT(6e20.12)

       j = SCAN(restart_file, '.' , .TRUE. )
       
       string = TRIM(restart_file(j-4:j-1))
       READ( string,* ) output_idx
       
       WRITE(*,*) 
       WRITE(*,*) 'Starting from output index ',output_idx
       
       ! Set this flag to 0 to not overwrite the initial condition
    
    ELSE
   
       WRITE(*,*) 'boh'

    END IF
 
    CLOSE(restart_unit)

      
  END SUBROUTINE read_solution

  !******************************************************************************
  !> \brief Write the solution on the output unit
  !
  !> This subroutine write the parameters of the grid, the output time 
  !> and the solution to a file with the name "run_name.q****", where  
  !> run_name is the name of the run read from the input file and ****
  !> is the counter of the output.
  !
  !> \param[in]   t      output time
  !
  !> @author 
  !> Mattia de' Michieli Vitturi
  !> \date 07/10/2016
  !
  !******************************************************************************

  SUBROUTINE output_solution(time)

    ! external procedures
    USE constitutive_2d, ONLY : qc_to_qp

    ! external variables
    USE constitutive_2d, ONLY : h , u , v , T

    USE geometry_2d, ONLY : comp_cells_x , B_cent , comp_cells_y , x_comp, y_comp
    ! USE geometry_2d, ONLY : x0 , dx , B_ver , y0 , dy

    USE parameters_2d, ONLY : n_vars
    USE parameters_2d, ONLY : t_output , dt_output

    USE solver_2d, ONLY : q

    IMPLICIT none

    REAL*8, INTENT(IN) :: time

    CHARACTER(LEN=4) :: idx_string

    REAL*8 :: qp(n_vars)

    INTEGER :: j,k
    INTEGER :: i

    output_idx = output_idx + 1

    idx_string = lettera(output_idx-1)

    IF ( output_cons_flag ) THEN
       
       output_file_2d = TRIM(run_name)//'_'//idx_string//'.q_2d'
       
       WRITE(*,*) 'WRITING ',output_file_2d
       
       OPEN(output_unit_2d,FILE=output_file_2d,status='unknown',form='formatted')
       
       !WRITE(output_unit_2d,1002) x0,dx,comp_cells_x,y0,dy,comp_cells_y,t
       
       DO k = 1,comp_cells_y
          
          DO j = 1,comp_cells_x
             
             DO i = 1,n_vars
                
                ! Exponents with more than 2 digits cause problems reading
                ! into matlab... reset tiny values to zero:
                IF ( dabs(q(i,j,k)) .LT. 1d-99) q(i,j,k) = 0.d0
                
             ENDDO

             IF ( temperature_flag ) THEN
             
                WRITE(output_unit_2d,1007) x_comp(j), y_comp(k), q(:,j,k)
                
             ELSE
                
                WRITE(output_unit_2d,1008) x_comp(j), y_comp(k), q(:,j,k)
                
             END IF

          ENDDO
          
          WRITE(output_unit_2d,*) ' ' 
          
       END DO
       
       WRITE(output_unit_2d,*) ' '
       WRITE(output_unit_2d,*) ' '
       
       CLOSE(output_unit_2d)
       
    END IF
    
    IF ( output_phys_flag ) THEN
       
       output_file_2d = TRIM(run_name)//'_'//idx_string//'.p_2d'
       
       WRITE(*,*) 'WRITING ',output_file_2d
       
       OPEN(output_unit_2d,FILE=output_file_2d,status='unknown',form='formatted')
       
       DO k = 1,comp_cells_y
          
          DO j = 1,comp_cells_x
             
             DO i = 1,n_vars
                
                ! Exponents with more than 2 digits cause problems reading
                ! into matlab... reset tiny values to zero:
                IF ( dabs(q(i,j,k)) .LT. 1d-99) q(i,j,k) = 0.d0
                
             END DO
             
             CALL qc_to_qp(q(:,j,k),B_cent(j,k),qp(:))
             
             IF ( temperature_flag ) THEN

                WRITE(output_unit_2d,1010) x_comp(j), y_comp(k), REAL(h),       &
                     REAL(u), REAL(v) , B_cent(j,k) , REAL(h) + B_cent(j,k) ,   &
                     REAL(T)
                
             ELSE

                WRITE(output_unit_2d,1009) x_comp(j), y_comp(k), REAL(h),       &
                     REAL(u), REAL(v) , B_cent(j,k) , REAL(h) + B_cent(j,k)

             END IF

          END DO
          
          WRITE(output_unit_2d,*) ' ' 
          
       ENDDO
       
       WRITE(output_unit_2d,*) ' '
       WRITE(output_unit_2d,*) ' '
       
       CLOSE(output_unit_2d)

    END IF

1007 FORMAT(6e20.12)
1008 FORMAT(5e20.12)
1009 FORMAT(7e20.12)
1010 FORMAT(8e20.12)

    t_output = time + dt_output

    IF ( output_esri_flag ) CALL output_esri(output_idx)

  END SUBROUTINE output_solution

  !******************************************************************************
  !> \brief Write the thickness in ESRI format
  !
  !> This subroutine write the thickness in the ascii ESRI format. 
  !> A masking is applied to the region with thickness less than 1D-5.
  !
  !> \param[in]   output_idx      output index
  !
  !> @author 
  !> Mattia de' Michieli Vitturi
  !> \date 15/12/2016
  !
  !******************************************************************************

  SUBROUTINE output_esri(output_idx)

    USE geometry_2d, ONLY : dx , dy , B_cent , B_ver , grid_output
    USE geometry_2d, ONLY : comp_interfaces_x , comp_interfaces_y
    USE solver_2d, ONLY : q

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: output_idx

    CHARACTER(LEN=4) :: idx_string
    INTEGER :: j
    
    IF ( output_idx .EQ. 1 ) THEN
       
       OPEN(dem_esri_unit,FILE='dem_esri.asc',status='unknown',form='formatted')
       
       WRITE(dem_esri_unit,'(A,I5)') 'ncols ', comp_cells_x
       WRITE(dem_esri_unit,'(A,I5)') 'nrows ', comp_cells_y
       WRITE(dem_esri_unit,'(A,F15.3)') 'xllcorner ', x0 + 0.5D0 * dx
       WRITE(dem_esri_unit,'(A,F15.3)') 'yllcorner ', y0 + 0.5D0 * dy
       WRITE(dem_esri_unit,'(A,F15.3)') 'cellsize ', cell_size
       WRITE(dem_esri_unit,'(A,I5)') 'NODATA_value ', -9999
        
       DO j = comp_cells_y,1,-1
          
          WRITE(dem_esri_unit,*) B_cent(1:comp_cells_x,j)
          
       ENDDO
       
       CLOSE(dem_esri_unit)

       OPEN(dem_esri_unit,FILE='dem_interfaces_esri.asc',status='unknown',      &
            form='formatted')
       
       WRITE(dem_esri_unit,'(A,I5)') 'ncols ', comp_interfaces_x
       WRITE(dem_esri_unit,'(A,I5)') 'nrows ', comp_interfaces_y
       WRITE(dem_esri_unit,'(A,F15.3)') 'xllcorner ', x0 
       WRITE(dem_esri_unit,'(A,F15.3)') 'yllcorner ', y0 
       WRITE(dem_esri_unit,'(A,F15.3)') 'cellsize ', cell_size
       WRITE(dem_esri_unit,'(A,I5)') 'NODATA_value ', -9999
        
       DO j = comp_interfaces_y,1,-1
          
          WRITE(dem_esri_unit,*) B_ver(1:comp_interfaces_x,j)
          
       ENDDO
       
       CLOSE(dem_esri_unit)

       
    END IF
    
    idx_string = lettera(output_idx-1)
    
    output_esri_file = TRIM(run_name)//'_'//idx_string//'.asc'

    WRITE(*,*) 'WRITING ',output_esri_file

    OPEN(output_esri_unit,FILE=output_esri_file,status='unknown',form='formatted')

    grid_output = -9999 

    WHERE ( q(1,:,:) - B_cent(:,:) .GE. 1.D-5 )

       grid_output = q(1,:,:) - B_cent(:,:)

    END WHERE

    WRITE(output_esri_unit,'(A,I5)') 'ncols ', comp_cells_x
    WRITE(output_esri_unit,'(A,I5)') 'nrows ', comp_cells_y
    WRITE(output_esri_unit,'(A,F15.3)') 'xllcorner ', x0 + 0.5D0 * dx
    WRITE(output_esri_unit,'(A,F15.3)') 'yllcorner ', y0 + 0.5D0 * dy
    WRITE(output_esri_unit,'(A,F15.3)') 'cellsize ', cell_size
    WRITE(output_esri_unit,'(A,I5)') 'NODATA_value ', -9999
        
    DO j = comp_cells_y,1,-1

       WRITE(output_esri_unit,*) grid_output(1:comp_cells_x,j)

    ENDDO
    
    CLOSE(output_esri_unit)
    
  END SUBROUTINE output_esri

  SUBROUTINE close_units

    IMPLICIT NONE

  END SUBROUTINE close_units

  !******************************************************************************
  !> \brief Numeric to String conversion
  !
  !> This function convert the integer in input into a numeric string for
  !> the subfix of the output files.
  !
  !> \param[in]   k      integer to convert         
  !
  !> \date 07/10/2016
  !> @author 
  !> Mattia de' Michieli Vitturi
  !
  !******************************************************************************

  CHARACTER*4 FUNCTION lettera(k)
    IMPLICIT NONE
    CHARACTER ones,tens,hund,thou
    !
    INTEGER :: k
    !
    INTEGER :: iten, ione, ihund, ithou
    !
    ithou=INT(k/1000)
    ihund=INT((k-(ithou*1000))/100)
    iten=INT((k-(ithou*1000)-(ihund*100))/10)
    ione=k-ithou*1000-ihund*100-iten*10
    ones=CHAR(ione+48)
    tens=CHAR(iten+48)
    hund=CHAR(ihunD+48)
    thou=CHAR(ithou+48)
    lettera=thou//hund//tens//ones
    !
    RETURN
  END FUNCTION lettera

  SUBROUTINE output_dakota

    IMPLICIT NONE

    dakota_file = 'IMEX_SfloW2D.out'

    OPEN(dakota_unit,FILE=dakota_file,status='unknown',form='formatted')

    CLOSE(dakota_unit)

  END SUBROUTINE output_dakota

END MODULE inpout_2d

