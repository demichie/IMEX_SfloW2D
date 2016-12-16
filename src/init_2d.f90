!*********************************************************************
!> \brief Initial solution
!
!> This module contains the variables and the subroutine for the
!> initialization of the solution for a Riemann problem.
!*********************************************************************

MODULE init_2d

  IMPLICIT none

  REAL*8, ALLOCATABLE :: q_init(:,:,:)
  
  !> Riemann problem interface relative position. It is a value
  !> between 0 and 1
  REAL*8 :: riemann_interface  


  REAL*8 :: hB_L         !< Left height
  REAL*8 :: u_L          !< Left velocity x
  REAL*8 :: v_L          !< Left velocity y
  
  REAL*8 :: hB_R         !< Right height
  REAL*8 :: u_R          !< Right velocity x
  REAL*8 :: v_R          !< Right velocity y



CONTAINS


  !*********************************************************************
  !> \brief Riemann problem initialization
  !
  !> This subroutine initialize the solution for a Riemann problem. The 
  !> values for the left and right states and the interface location 
  !> are read from the input file.\
  !> \date 26/08/2011
  !*********************************************************************

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

    riemann_int_search:DO j = 1,comp_cells_x

       IF ( x_comp(j) .LT. riemann_interface ) THEN

          i1 = j

       ELSE

          EXIT riemann_int_search

       END IF

    END DO riemann_int_search

    eps = 1.D-10

    ! Left initial state
    !qp(1) = hB_L
    !qp(2) = u_L
    !qp(3) = v_L

    qp(1,1:i1,:) = hB_L

    qp(2,1:i1,:) = u_L

    qp(3,1:i1,:) = v_L

    DO j = 1,i1

       DO k = 1,comp_cells_y

         ! evaluate the vector of conservative variables
         CALL qp_to_qc( qp(:,j,k) , B_cent(j,k) , qj )

         q(1:n_vars,j,k) = qj

         IF ( verbose_level .GE. 1 ) WRITE(*,*) j,k,B_cent(j,k),qp(:,j,k)

       ENDDO

    END DO

    ! Right initial state
    !qp(1) = hB_R
    !qp(2) = u_R
    !qp(3) = v_R

    qp(1,i1+1:comp_cells_x,:) = hB_R

    qp(2,i1+1:comp_cells_x,:) = u_R

    qp(3,i1+1:comp_cells_x,:) = v_R

    DO j = i1+1,comp_cells_x

       DO k = 1,comp_cells_y

         ! evaluate the vector of conservative variables
         CALL qp_to_qc( qp(:,j,k) , B_cent(j,k) , qj )

         q(1:n_vars,j,k) = qj

         IF ( verbose_level .GE. 1 ) WRITE(*,*) j,k,B_cent(j,k),qp(:,j,k)
    
      END DO

    ENDDO

    IF ( verbose_level .GE. 1 ) READ(*,*)

    RETURN

  END SUBROUTINE riemann_problem

  !*********************************************************************
  !> \brief Problem initialization
  !
  !> This subroutine initialize the solution for a a generic problem. The 
  !> values are read from thickness and velocity functions
  !> \date OCTOBER 2016
  !*********************************************************************

  SUBROUTINE initial_conditions

    USE constitutive_2d, ONLY : qp_to_qc

    USE geometry_2d, ONLY : x_comp , comp_cells_x, y_comp , comp_cells_y , B_cent

    ! USE geometry_2d, ONLY : x0 , xN , y0 , yN

    USE parameters_2d, ONLY : n_vars , verbose_level

    USE solver_2d, ONLY : q

    IMPLICIT none

    REAL*8 :: qp(n_vars,comp_cells_x,comp_cells_y) , qj(n_vars)

    INTEGER :: j,k          !< loop counter

    DO j=1,comp_cells_x

      DO k=1,comp_cells_y

         qp(1,j,k) = thickness_function(x_comp(j),y_comp(k),B_cent(j,k))

         qp(2,j,k) = velocity_u_function(x_comp(j),y_comp(k),B_cent(j,k))

         qp(3,j,k) = velocity_v_function(x_comp(j),y_comp(k),B_cent(j,k))

      ENDDO

    ENDDO

    DO j = 1,comp_cells_x

       DO k = 1,comp_cells_y

         ! evaluate the vector of conservative variables
         CALL qp_to_qc( qp(:,j,k) , B_cent(j,k) , qj )

         q(1:n_vars,j,k) = qj

         IF ( verbose_level .GE. 1 ) WRITE(*,*) j,k,B_cent(j,k),qp(:,j,k)
    
      END DO

    ENDDO

    IF ( verbose_level .GE. 1 ) READ(*,*)

    RETURN

  END SUBROUTINE initial_conditions


!---------------------------------------------------------------------------
!> Thickness function
!
!> This subroutine defines thickness height hB=h+B
!> in the input (x,y) grid point
!> \date OCTOBER 2016
!> \param    x           original grid                (\b input)
!> \param    y           original grid                (\b input)
!> \param    Bj          original grid                (\b input)
!---------------------------------------------------------------------------
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

!---------------------------------------------------------------------------
!> Velocity u function
!
!> This subroutine defines x component of the velocity
!> in the input (x,y) grid point
!> \date OCTOBER 2016
!> \param    x           original grid                (\b input)
!> \param    y           original grid                (\b input)
!---------------------------------------------------------------------------
  REAL*8 FUNCTION velocity_u_function(x,y,Bj)
    IMPLICIT NONE
    
    REAL*8, INTENT(IN) :: x,y
    REAL*8, INTENT(IN) :: Bj
    
    ! example 2D from Kurganov and Petrova 2007    
    !IF(ABS(y).LE.0.5)THEN
    !   velocity_u_function=0.5
    !ELSE
    !   velocity_u_function=0.0
    !ENDIF

    velocity_u_function=0.0

  END FUNCTION velocity_u_function

!---------------------------------------------------------------------------
!> Velocity v function
!
!> This subroutine defines y component of the velocity
!> in the input (x,y) grid point
!> \date OCTOBER 2016
!> \param    x           original grid                (\b input)
!> \param    y           original grid                (\b input)
!---------------------------------------------------------------------------
  REAL*8 FUNCTION velocity_v_function(x,y,Bj)
    IMPLICIT NONE
    
    REAL*8, INTENT(IN) :: x,y
    REAL*8, INTENT(IN) :: Bj
   
    ! example 2D from Kurganov and Petrova 2007    
    !velocity_v_function=0.0

    velocity_v_function=0.0

  END FUNCTION velocity_v_function
  

END MODULE init_2d
