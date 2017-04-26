!********************************************************************************
!> \brief Initial solution
!
!> This module contains the variables and the subroutine for the
!> initialization of the solution for a Riemann problem.
!********************************************************************************

MODULE init_2d

  USE parameters_2d, ONLY : temperature_flag , verbose_level

  IMPLICIT none

  REAL*8, ALLOCATABLE :: q_init(:,:,:)

  REAL*8, ALLOCATABLE :: thickness_init(:,:)
  
  !> Riemann problem interface relative position. It is a value
  !> between 0 and 1
  REAL*8 :: riemann_interface  


  REAL*8 :: hB_W         !< Left height
  REAL*8 :: u_W          !< Left velocity x
  REAL*8 :: v_W          !< Left velocity y
  REAL*8 :: T_W          !< Left temperature

  REAL*8 :: hB_E         !< Right height
  REAL*8 :: u_E          !< Right velocity x
  REAL*8 :: v_E          !< Right velocity y
  REAL*8 :: T_E          !< Right temperature


CONTAINS


  !******************************************************************************
  !> \brief Riemann problem initialization
  !
  !> This subroutine initialize the solution for a Riemann problem. The 
  !> values for the left and right states and the interface location 
  !> are read from the input file.\
  !> \date 26/08/2011
  !******************************************************************************

  SUBROUTINE riemann_problem

    USE constitutive_2d, ONLY : qp_to_qc

    USE geometry_2d, ONLY : x_comp , comp_cells_x , comp_cells_y , B_cent

    ! USE geometry_2d, ONLY : x0 , xN , y0 , yN

    USE parameters_2d, ONLY : n_vars , verbose_level

    USE solver_2d, ONLY : q

    IMPLICIT none

    ! REAL*8 :: hB            !< height + topography
    ! REAL*8 :: u             !< velocity
    ! REAL*8 :: v             !< velocity

    REAL*8 :: qp(n_vars,comp_cells_x,comp_cells_y) , qj(n_vars)

    INTEGER :: j,k          !< loop counter
    INTEGER :: i1           !< last index with left state

    REAL*8 :: eps

    IF ( verbose_level .GE. 1 ) WRITE(*,*) 'Riemann problem initialization'

    i1 = 0
    
    riemann_int_search:DO j = 1,comp_cells_x
       
       IF ( x_comp(j) .LT. riemann_interface ) THEN

          i1 = j

       ELSE

          EXIT riemann_int_search

       END IF

    END DO riemann_int_search
    
    eps = 1.D-10

    ! Left initial state
    qp(1,1:i1,:) = hB_W
    qp(2,1:i1,:) = u_W
    qp(3,1:i1,:) = v_W
    IF ( temperature_flag ) qp(4,1:i1,:) = T_W


    IF ( verbose_level .GE. 1 ) WRITE(*,*) 'Left state'

    DO j = 1,i1

       DO k = 1,comp_cells_y

         ! evaluate the vector of conservative variables
         CALL qp_to_qc( qp(:,j,k) , B_cent(j,k) , qj )

         q(1:n_vars,j,k) = qj

         IF ( verbose_level .GE. 1 ) THEN 
            
            WRITE(*,*) j,k,B_cent(j,k)
            WRITE(*,*) qp(:,j,k)
            WRITE(*,*) q(1:n_vars,j,k)

         END IF

       ENDDO

    END DO

    IF ( verbose_level .GE. 1 ) READ(*,*)

    ! Right initial state
    qp(1,i1+1:comp_cells_x,:) = hB_E
    qp(2,i1+1:comp_cells_x,:) = u_E
    qp(3,i1+1:comp_cells_x,:) = v_E
    IF ( temperature_flag ) qp(4,i1+1:comp_cells_x,:) = T_E


    IF ( verbose_level .GE. 1 ) WRITE(*,*) 'Right state'

    DO j = i1+1,comp_cells_x

       DO k = 1,comp_cells_y

         ! evaluate the vector of conservative variables
         CALL qp_to_qc( qp(:,j,k) , B_cent(j,k) , qj )

         q(1:n_vars,j,k) = qj

         IF ( verbose_level .GE. 1 ) THEN 
            
            WRITE(*,*) j,k,B_cent(j,k)
            WRITE(*,*) qp(:,j,k)
            WRITE(*,*) q(1:n_vars,j,k)

         END IF
    
      END DO

    ENDDO

    IF ( verbose_level .GE. 1 ) READ(*,*)

    RETURN

  END SUBROUTINE riemann_problem

  !******************************************************************************
  !> \brief Problem initialization
  !
  !> This subroutine initialize the solution for a a generic problem. The 
  !> values are read from thickness and velocity functions
  !> \date OCTOBER 2016
  !******************************************************************************

  SUBROUTINE initial_conditions

    USE constitutive_2d, ONLY : qp_to_qc

    USE geometry_2d, ONLY : x_comp , comp_cells_x, y_comp , comp_cells_y ,      &
         B_cent

    USE geometry_2d, ONLY : dx , dy

    USE parameters_2d, ONLY : n_vars , released_volume , verbose_level ,        &
         source_flag

    USE solver_2d, ONLY : q

    IMPLICIT none

    REAL*8 :: qp(n_vars,comp_cells_x,comp_cells_y) , qj(n_vars)

    INTEGER :: j,k          !< loop counter

    DO j=1,comp_cells_x

      DO k=1,comp_cells_y

         qp(1,j,k) = thickness_function(x_comp(j),y_comp(k),B_cent(j,k))

         qp(2,j,k) = velocity_u_function(x_comp(j),y_comp(k),B_cent(j,k))

         qp(3,j,k) = velocity_v_function(x_comp(j),y_comp(k),B_cent(j,k))

         IF ( temperature_flag ) THEN

            qp(4,j,k) = temperature_function(x_comp(j),y_comp(k))

         END IF

      ENDDO

    ENDDO

    IF ( released_volume .GT. ( dx * dy * SUM( qp(1,:,:)-B_cent(:,:) ) ) ) THEN

       ! Correction for the released volume
       qp(1,:,:) = B_cent(:,:) + ( qp(1,:,:)-B_cent(:,:) ) * released_volume    &
            / ( dx * dy * SUM( qp(1,:,:)-B_cent(:,:) ) )
       
    END IF

    WRITE(*,*) 'Initial volume =',dx * dy * SUM( qp(1,:,:)-B_cent(:,:) ) 


    DO j = 1,comp_cells_x

       DO k = 1,comp_cells_y

         ! evaluate the vector of conservative variables
         CALL qp_to_qc( qp(:,j,k) , B_cent(j,k) , qj )

         q(1:n_vars,j,k) = qj

         IF ( verbose_level .GE. 1 ) WRITE(*,*) j,k,B_cent(j,k),qp(:,j,k)
    
      END DO

    ENDDO

    IF ( verbose_level .GE. 1 ) READ(*,*)

    IF ( source_flag ) CALL init_source

    RETURN

  END SUBROUTINE initial_conditions

  !******************************************************************************
  !> \brief Source initialization
  !
  !> This subroutine initialize the source terms. 
  !> \date 25/04/2017
  !******************************************************************************

  SUBROUTINE init_source

    USE geometry_2d, ONLY : x_comp , comp_cells_x, y_comp , comp_cells_y

    USE parameters_2d, ONLY : x_source , y_source , r_source

    USE solver_2d, ONLY : source_xy
    
    IMPLICIT NONE

    INTEGER :: j,k          !< loop counter

    source_xy(:,:) = 0.D0

    DO j=1,comp_cells_x

      DO k=1,comp_cells_y

         IF ( ( x_comp(j) - x_source )**2 + ( y_comp(k) - y_source ) **2 .LE.   &
              r_source**2 ) THEN

            !WRITE(*,*) x_comp(j) , x_source , y_comp(k) , y_source , r_source
            !WRITE(*,*) ( x_comp(j) - x_source )**2 + ( y_comp(k) - y_source ) **2
            !WRITE(*,*) r_source**2
            !READ(*,*) 

            source_xy(j,k) = 1.D0 
         
         END IF

      ENDDO

    ENDDO
    

  END SUBROUTINE init_source
  
  !******************************************************************************
  !> Thickness function
  !
  !> This subroutine defines thickness height hB=h+B
  !> in the input (x,y) grid point
  !> \date OCTOBER 2016
  !> \param    x           original grid                (\b input)
  !> \param    y           original grid                (\b input)
  !> \param    Bj          original grid                (\b input)
  !******************************************************************************

  REAL*8 FUNCTION thickness_function(x,y,Bj)

    USE parameters_2d, ONLY : released_volume , x_release , y_release

    IMPLICIT NONE
    
    REAL*8, INTENT(IN) :: x,y
    REAL*8, INTENT(IN) :: Bj
    
    REAL*8, PARAMETER :: pig = 4.0*ATAN(1.0)
    REAL*8 :: R

    ! example 1D from Kurganov and Petrova 2007    
    !IF(y.LE.0.d0)THEN
    !   thickness_function=0.4d0+Bj
    !ELSE
    !   thickness_function=Bj
    !ENDIF

    ! example 2D from Kurganov and Petrova 2007    
    ! thickness_function=MAX(0.25,Bj)
    
    R = ( released_volume / pig )**(1.0/3.0)
    
    IF ( DSQRT( (x-x_release)**2 + (y-y_release)**2 ) .LE. R ) THEN

      thickness_function = Bj + R

    ELSE

      thickness_function = Bj

    ENDIF

  END FUNCTION thickness_function

  !******************************************************************************
  !> Velocity u function
  !
  !> This subroutine defines x component of the velocity
  !> in the input (x,y) grid point
  !> \date OCTOBER 2016
  !> \param    x           original grid                (\b input)
  !> \param    y           original grid                (\b input)
  !******************************************************************************

  REAL*8 FUNCTION velocity_u_function(x,y,Bj)
  
    USE parameters_2d, ONLY : released_volume , x_release , y_release
    USE parameters_2d, ONLY : velocity_mod_release, velocity_ang_release

    USE geometry_2d, ONLY : x0 , y0 , dx , dy
    USE geometry_2d, ONLY : x_stag , y_stag , B_ver
    
    IMPLICIT NONE
    
    REAL*8, INTENT(IN) :: x,y
    REAL*8, INTENT(IN) :: Bj
    
    REAL*8, PARAMETER :: pig = 4.0*ATAN(1.0)
    REAL*8 :: R

    INTEGER :: idx1 , idy1
    REAL*8 :: coeff_x , coeff_y

    REAL*8 :: dBdx1 , dBdx2 , dBdx
    REAL*8 :: dBdy1 , dBdy2 , dBdy

    REAL*8 :: max_slope_angle_rad
    REAL*8 :: velocity_ang_release_rad

    
    R = ( released_volume / pig )**(1.d0/3.d0)

    ! indexes of the lower-left vertex of the cell in the computational grid
    idx1 = CEILING( ( x_release - x0 ) / dx )
    idy1 = CEILING( ( y_release - y0 ) / dy )

    ! relative coordinates of (x_release,y_release) within the cell
    coeff_x = ( x_release - x_stag(idx1) ) / dx
    coeff_y = ( y_release - y_stag(idy1) ) / dy


    ! x-derivatives at bottom and top edges of the cell
    dBdx1 = ( B_ver(idx1+1,idy1) - B_ver(idx1,idy1) ) / dx
    dBdx2 = ( B_ver(idx1+1,idy1+1) - B_ver(idx1,idy1+1) ) / dx

    ! x-derivative of terrain at (x_release,y_release)
    dBdx = coeff_y * dBdx2 + ( 1.d0 - coeff_y ) * dBdx1
    
    ! y-derivatives at left and right edges of the cell
    dBdy1 = ( B_ver(idx1,idy1+1) - B_ver(idx1,idy1) ) / dy
    dBdy2 = ( B_ver(idx1+1,idy1+1) - B_ver(idx1+1,idy1) ) / dx

    ! y-derivative of terrain at (x_release,y_release)
    dBdy = coeff_x * dBdy2 + ( 1.d0 - coeff_x ) * dBdy1


    ! direction of maximum slope in radians
    max_slope_angle_rad = datan2(dbdy,dbdx)

    IF ( verbose_level .GE. 1 ) THEN

       WRITE(*,*) x0,x_stag(1),y0,y_stag(1),dx,dy
       WRITE(*,*) '---',x_stag(idx1),x_release,x_stag(idx1+1),coeff_x
       WRITE(*,*) '---',y_stag(idy1),y_release,y_stag(idy1+1),coeff_y
       
       WRITE(*,*) B_ver(idx1,idy1),B_ver(idx1+1,idy1)
       WRITE(*,*) B_ver(idx1,idy1+1),B_ver(idx1+1,idy1+1)
       WRITE(*,*) 'dbdx,dbdy,angle',dbdx,dbdy,max_slope_angle_rad*180.d0/pig
    
       READ(*,*)

    END IF
    
    IF ( DSQRT( (x-x_release)**2 + (y-y_release)**2 ) .LE. R ) THEN

       velocity_ang_release_rad = velocity_ang_release * pig / 180.d0
       
       IF ( velocity_ang_release .GE. 0.D0 ) THEN 

          ! angle departing from positive x-axis
          velocity_u_function = velocity_mod_release                            &
               * COS( velocity_ang_release * ( 2.D0 * pig / 360.d0 ) )

       ELSE

          ! angle departing from maximum slope direction
          velocity_u_function = velocity_mod_release                            &
               * COS( max_slope_angle_rad + velocity_ang_release_rad )

       END IF
          
    ELSE

      velocity_u_function = 0.d0

    ENDIF

  END FUNCTION velocity_u_function

  !******************************************************************************
  !> Velocity v function
  !
  !> This subroutine defines y component of the velocity
  !> in the input (x,y) grid point
  !> \date OCTOBER 2016
  !> \param    x           original grid                (\b input)
  !> \param    y           original grid                (\b input)
  !******************************************************************************
  
  REAL*8 FUNCTION velocity_v_function(x,y,Bj)
  
    USE parameters_2d, ONLY : released_volume , x_release , y_release
    USE parameters_2d, ONLY : velocity_mod_release, velocity_ang_release

    USE geometry_2d, ONLY : x0 , y0 , dx , dy
    USE geometry_2d, ONLY : x_stag , y_stag , B_ver

    IMPLICIT NONE
    
    REAL*8, INTENT(IN) :: x,y
    REAL*8, INTENT(IN) :: Bj
   
    REAL*8, PARAMETER :: pig = 4.0*ATAN(1.0)
    REAL*8 :: R
    
    INTEGER :: idx1 , idy1
    REAL*8 :: coeff_x , coeff_y
    
    REAL*8 :: dBdx1 , dBdx2 , dBdx
    REAL*8 :: dBdy1 , dBdy2 , dBdy

    REAL*8 :: max_slope_angle_rad
    REAL*8 :: velocity_ang_release_rad


    R = ( released_volume / pig )**(1.0/3.0)

    ! indexes of the lower-left vertex of the cell in the computational grid
    idx1 = CEILING( ( x_release - x0 ) / dx )
    idy1 = CEILING( ( y_release - y0 ) / dy )

    ! relative coordinates of (x_release,y_release) within the cell
    coeff_x = ( x_release - x_stag(idx1) ) / dx
    coeff_y = ( y_release - y_stag(idy1) ) / dy


    ! x-derivatives at bottom and top edges of the cell
    dBdx1 = ( B_ver(idx1+1,idy1) - B_ver(idx1,idy1) ) / dx
    dBdx2 = ( B_ver(idx1+1,idy1+1) - B_ver(idx1,idy1+1) ) / dx

    ! x-derivative of terrain at (x_release,y_release)
    dBdx = coeff_y * dBdx2 + ( 1.d0 - coeff_y ) * dBdx1
    
    ! y-derivatives at left and right edges of the cell
    dBdy1 = ( B_ver(idx1,idy1+1) - B_ver(idx1,idy1) ) / dy
    dBdy2 = ( B_ver(idx1+1,idy1+1) - B_ver(idx1+1,idy1) ) / dx

    ! y-derivative of terrain at (x_release,y_release)
    dBdy = coeff_x * dBdy2 + ( 1.d0 - coeff_x ) * dBdy1


    
    IF ( DSQRT( (x-x_release)**2 + (y-y_release)**2 ) .LE. R ) THEN

       velocity_ang_release_rad = velocity_ang_release * pig / 180.d0
       
       IF ( velocity_ang_release .GE. 0.D0 ) THEN 
          
          ! angle departing from positive x-axis
          velocity_v_function = velocity_mod_release                            &
               * SIN( velocity_ang_release * ( 2.D0 * pig / 360.d0 ) )
          
       ELSE
          
          ! angle departing from maximum slope direction
          velocity_v_function = velocity_mod_release                            &
               * SIN( max_slope_angle_rad + velocity_ang_release_rad )
          
       END IF

    ELSE

      velocity_v_function = 0.d0

    ENDIF

  END FUNCTION velocity_v_function
  

  !******************************************************************************
  !> Temperature function
  !
  !> This subroutine defines the temperature in the pile and outside as a function 
  !> of the input (x,y) grid point
  !> \date OCTOBER 2016
  !> \param    x           original grid                (\b input)
  !> \param    y           original grid                (\b input)
  !******************************************************************************

  REAL*8 FUNCTION temperature_function(x,y)
  
    USE parameters_2d, ONLY : released_volume , x_release , y_release
    USE parameters_2d, ONLY : T_init , T_ambient
  
    IMPLICIT NONE
    
    REAL*8, INTENT(IN) :: x,y
   
    REAL*8, PARAMETER :: pig = 4.0*ATAN(1.0)
    REAL*8 :: R
   
    R = ( released_volume / pig )**(1.0/3.0)
    
    IF ( DSQRT( (x-x_release)**2 + (y-y_release)**2 ) .LE. R ) THEN

       temperature_function = T_init

    ELSE

       temperature_function = T_ambient

    ENDIF

  END FUNCTION temperature_function
  
END MODULE init_2d
