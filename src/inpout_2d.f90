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
  USE parameters_2d, ONLY : solver_scheme, max_dt , t_start , t_end ,           &
       dt_output , cfl, limiter , theta, reconstr_coeff ,                       &
       interfaces_relaxation , n_RK, topography_function_flag,                  &
       topography_demfile, riemann_flag

  USE solver_2d, ONLY : verbose_level

  ! -- Variables for the namelist NEWRUN_PARAMETERS
  USE geometry_2d, ONLY : x0 , y0 , comp_cells_x , comp_cells_y , cell_size
  USE geometry_2d, ONLY : topography_profile , n_topography_profile_x ,         &
       n_topography_profile_y
  USE init_2d, ONLY : riemann_interface
  USE parameters_2d, ONLY : released_volume , x_release , y_release

  ! -- Variables for the namelist LEFT_STATE
  USE init_2d, ONLY : hB_L , u_L , v_L

  ! -- Variables for the namelist RIGHT_STATE
  USE init_2d, ONLY : hB_R , u_R , v_R

  ! -- Variables for the namelists LEFT/RIGHT_BOUNDARY_CONDITIONS
  USE parameters_2d, ONLY : bc

  ! -- Variables for the namelist SOURCE_PARAMETERS
  USE constitutive_2d, ONLY : grav, mu, xi

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
  !> - T     => Restart from a previous output
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

  ! -- Variables for the namelists LEFT_BOUNDARY_CONDITIONS
  TYPE(bc) :: hB_bcL , u_bcL , v_bcL

  ! -- Variables for the namelists RIGHT_BOUNDARY_CONDITIONS
  TYPE(bc) :: hB_bcR , u_bcR , v_bcR

  ! -- Variables for the namelists BOTTOM_BOUNDARY_CONDITIONS
  TYPE(bc) :: hB_bcD , u_bcD , v_bcD

  ! -- Variables for the namelists TOP_BOUNDARY_CONDITIONS
  TYPE(bc) :: hB_bcU , u_bcU , v_bcU


  ! parameters to read a dem file
  INTEGER :: ncols, nrows, nodata_value

  REAL*8 :: xllcorner, yllcorner, cellsize


  NAMELIST / run_parameters / run_name , restart , topography_demfile ,         &
       riemann_flag , max_dt , t_start , t_end , dt_output , solver_scheme ,    &
       cfl , limiter , theta , reconstr_coeff , n_RK , output_cons_flag ,       &
       output_esri_flag , output_phys_flag , verbose_level

  NAMELIST / restart_parameters / restart_file

  NAMELIST / newrun_parameters / x0 , y0 , comp_cells_x , comp_cells_y ,        &
       cell_size , released_volume , x_release , y_release

  NAMELIST / left_state / hB_L , u_L , v_L

  NAMELIST / right_state / hB_R , u_R , v_R

  NAMELIST / left_boundary_conditions / hB_bcL , u_bcL , v_bcL

  NAMELIST / right_boundary_conditions / hB_bcR , u_bcR , v_bcR

  NAMELIST / bottom_boundary_conditions / hB_bcD , u_bcD , v_bcD

  NAMELIST / top_boundary_conditions / hB_bcU , u_bcU , v_bcU

  NAMELIST / source_parameters / grav, mu, xi

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
    riemann_flag=.TRUE.
    max_dt = 1.d-3
    t_start = 0.0
    t_end = 5.0d-4
    dt_output = 1.d-4
    solver_scheme = 'KT'
    n_RK = 2
    output_cons_flag = .TRUE.
    output_esri_flag = .TRUE.
    output_phys_flag = .TRUE.
    verbose_level = 0


    cfl = 0.5
    limiter(1:n_vars) = 0
    theta=1.0
    reconstr_coeff = 1.0

    !-- Inizialization of the Variables for the namelist restart parameters
    restart_file = ''

    !-- Inizialization of the Variables for the namelist newrun_parameters
    x0 = 0.D0
    comp_cells_x = 500
    y0 = 0.D0
    comp_cells_y = 500
    cell_size = 0.1
    riemann_interface = 0.5D0

    !-- Inizialization of the Variables for the namelist left_state
    hB_L = 1.D0
    u_L = 0.D0
    v_L = 0.D0

    !-- Inizialization of the Variables for the namelist right_state
    hB_R = 0.5D0
    u_R = 0.D0
    u_R = 0.D0

    !-- Inizialization of the Variables for the namelist left boundary conditions

    hB_bcL%flag = 1 
    hB_bcL%value = 0.d0 

    u_bcL%flag = 1 
    u_bcL%value = 0.d0 

    !-- Inizialization of the Variables for the namelist right boundary conditions

    hB_bcR%flag = 1 
    hB_bcR%value = 0.d0 

    u_bcR%flag = 1 
    u_bcR%value = 0.d0 

    !-- Inizialization of the Variables for the namelist bottom boundary conditions

    hB_bcD%flag = 1 
    hB_bcD%value = 0.d0 

    u_bcD%flag = 1 
    u_bcD%value = 0.d0 

    !-- Inizialization of the Variables for the namelist top boundary conditions

    hB_bcU%flag = 1 
    hB_bcU%value = 0.d0 

    u_bcU%flag = 1 
    u_bcU%value = 0.d0 


    !-- Inizialization of the Variables for the namelist source_parameters
    grav = -9.81D0
    mu = 0.225
    xi = 130

    input_file = 'shallow_water_2d.inp'

    INQUIRE (FILE=input_file,exist=lexist)

    IF (lexist .EQV. .FALSE.) THEN

       OPEN(input_unit,FILE=input_file,STATUS='NEW')

       WRITE(input_unit, run_parameters )
       WRITE(input_unit, restart_parameters )
       WRITE(input_unit, newrun_parameters )
       WRITE(input_unit, left_state )
       WRITE(input_unit, right_state )
       WRITE(input_unit, left_boundary_conditions )
       WRITE(input_unit, right_boundary_conditions )
       WRITE(input_unit, source_parameters )

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

       WRITE(*,*) 'Input file shallow_water_2d.inp not found'
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

    USE parameters_2d, ONLY : bcL , bcR , bcD , bcU

    IMPLICIT none

    REAL*8 :: max_cfl

    LOGICAL :: tend1 
    CHARACTER(LEN=80) :: card

    INTEGER :: j,k

    LOGICAL :: lexist

    CHARACTER(LEN=15) :: chara

    OPEN(input_unit,FILE=input_file,STATUS='old')


    ! ------- READ run_parameters NAMELIST -----------------------------------
    READ(input_unit, run_parameters )

    IF ( ( solver_scheme .NE. 'LxF' ) .AND. ( solver_scheme .NE. 'KT' ) .AND. &
         ( solver_scheme .NE. 'GFORCE' ) ) THEN

       WRITE(*,*) 'WARNING: no correct solver scheme selected',solver_scheme
       WRITE(*,*) 'Chose between: LxF, GFORCE or KT'
       STOP

    END IF

    IF  ( ( solver_scheme.EQ.'LxF' ) .OR. ( solver_scheme.EQ.'GFORCE' ) ) THEN 

       max_cfl = 1.0

    ELSEIF ( solver_scheme .EQ. 'KT' ) THEN

       max_cfl = 0.5

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

    IF ( restart ) THEN

       READ(input_unit,restart_parameters)

    ELSE

       READ(input_unit,newrun_parameters)

       IF ( riemann_flag ) THEN

          READ(input_unit,left_state)
          READ(input_unit,right_state)
          
       END IF

    END IF

    READ(input_unit,left_boundary_conditions)

    bcL(1) = hB_bcL 
    bcL(2) = u_bcL 
    bcL(3) = v_bcL 

    READ(input_unit,right_boundary_conditions)

    bcR(1) = hB_bcR 
    bcR(2) = u_bcR 
    bcR(3) = v_bcR 

    READ(input_unit,bottom_boundary_conditions)

    bcD(1) = hB_bcD 
    bcD(2) = u_bcD 
    bcD(3) = v_bcD 

    READ(input_unit,top_boundary_conditions)

    bcU(1) = hB_bcU 
    bcU(2) = u_bcU 
    bcU(3) = v_bcU 

    READ(input_unit, source_parameters )

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

       END IF

    END IF

    WRITE(backup_unit,left_boundary_conditions)
    WRITE(backup_unit,right_boundary_conditions)

    WRITE(backup_unit, source_parameters )


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

    USE geometry_2d, ONLY : init_grid
    USE solver_2d, ONLY : allocate_solver_variables

    USE geometry_2d, ONLY : comp_cells_x , x0 , dx, comp_cells_y , y0 , dy
    USE parameters_2d, ONLY : n_vars
    USE parameters_2d, ONLY : t_start
    USE solver_2d, ONLY : q

    IMPLICIT none

    INTEGER :: j,k
    INTEGER :: i

    CHARACTER(LEN=30) :: string

    OPEN(restart_unit,FILE=restart_file,STATUS='old')

    READ(restart_unit,1001) x0,dx,comp_cells_x,y0,dy,comp_cells_y,t_start
1001 FORMAT(e18.8,'    x0', /, e18.8,'    dx', /, i18,'    cells_x', /,         &
          e18.8,'    y0', /, e18.8,'    dy', /, i18,'    cells_y', /, e18.8,    &
          '    t', /)


    CALL init_grid

    CALL allocate_solver_variables


    DO i = 1,n_vars

       DO j = 1,comp_cells_x

          DO k = 1,comp_cells_y

             ! Exponents with more than 2 digits cause problems reading
             ! into matlab... reset tiny values to zero:
             IF ( dabs(q(i,j,k)) .LT. 1d-99) q(i,j,k) = 0.d0

          ENDDO

       END DO

       READ(restart_unit,1003) q(:,j,k)
1003   FORMAT(3e20.12)

    END DO

    j = SCAN(restart_file, '.' , .TRUE. )

    string = TRIM(restart_file(j+2:j+5))

    READ( string,* ) output_idx

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

  SUBROUTINE output_solution(t)

    ! external procedures
    USE constitutive_2d, ONLY : qc_to_qp

    ! external variables
    USE constitutive_2d, ONLY : h , u , v

    USE geometry_2d, ONLY : comp_cells_x , B_cent , comp_cells_y , x_comp, y_comp
    ! USE geometry_2d, ONLY : x0 , dx , B_ver , y0 , dy

    USE parameters_2d, ONLY : n_vars
    USE parameters_2d, ONLY : t_output , dt_output

    USE solver_2d, ONLY : q

    IMPLICIT none

    REAL*8, INTENT(IN) :: t

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
             
             WRITE(output_unit_2d,1008) x_comp(j), y_comp(k), q(:,j,k)
             
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
             
             WRITE(output_unit_2d,1009) x_comp(j), y_comp(k), REAL(h), REAL(u), &
                  REAL(v) , B_cent(j,k) , REAL(h) + B_cent(j,k)
             
          END DO
          
          WRITE(output_unit_2d,*) ' ' 
          
       ENDDO
       
       WRITE(output_unit_2d,*) ' '
       WRITE(output_unit_2d,*) ' '
       
       CLOSE(output_unit_2d)

    END IF

1008 FORMAT(5e20.12)
1009 FORMAT(7e20.12)

    t_output = t + dt_output

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

    USE geometry_2d, ONLY : dx , dy , B_cent , grid_output
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

    dakota_file = 'shallow_water_2d.out'

    OPEN(dakota_unit,FILE=dakota_file,status='unknown',form='formatted')

    CLOSE(dakota_unit)

  END SUBROUTINE output_dakota

END MODULE inpout_2d

