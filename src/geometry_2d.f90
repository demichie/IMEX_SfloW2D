!*********************************************************************
!> \brief Grid module
!
!> This module contains the variables and the subroutines related to 
!> the computational grid
!*********************************************************************
MODULE geometry_2d

  USE parameters_2d, ONLY : verbose_level

  IMPLICIT NONE

  !> Location of the centers (x) of the control volume of the domain
  REAL*8, ALLOCATABLE :: x_comp(:)

  !> Location of the boundaries (x) of the control volumes of the domain
  REAL*8, ALLOCATABLE :: x_stag(:)

  !> Location of the centers (y) of the control volume of the domain
  REAL*8, ALLOCATABLE :: y_comp(:)

  !> Location of the boundaries (x) of the control volumes of the domain
  REAL*8, ALLOCATABLE :: y_stag(:)

  !> Topography at the boundaries (x) of the control volumes
  REAL*8, ALLOCATABLE :: B_stag_x(:,:)

  !> Topography at the boundaries (y) of the control volumes
  REAL*8, ALLOCATABLE :: B_stag_y(:,:)

  !> Topography at the vertices of the control volumes
  REAL*8, ALLOCATABLE :: B_ver(:,:)

  !> Topography at the centers of the control volumes 
  REAL*8, ALLOCATABLE :: B_cent(:,:)

  !> Topography slope (x‪ direction) at the centers of the control volumes 
  REAL*8, ALLOCATABLE :: B_prime_x(:,:)

  !> Topography slope (y direction) at the centers of the control volumes 
  REAL*8, ALLOCATABLE :: B_prime_y(:,:)

  !> Solution in ascii grid format (ESRI)
  REAL*8, ALLOCATABLE :: grid_output(:,:)

  !> gravity vector wrt surface coordinates for each cell
  REAL*8, ALLOCATABLE :: grav_surf(:,:)

  !> curvature wrt mixed directions for each cell
  REAL*8, ALLOCATABLE :: curv_xy(:,:)

  REAL*8, ALLOCATABLE :: topography_profile(:,:,:)

  INTEGER :: n_topography_profile_x, n_topography_profile_y

  REAL*8 :: dx                 !< Control volumes size
  REAL*8 :: x0                 !< Left of the physical domain
  REAL*8 :: xN                 !< Right of the physical domain
  REAL*8 :: dy                 !< Control volumes size
  REAL*8 :: y0                 !< Bottom of the physical domain
  REAL*8 :: yN                 !< Top of the physical domain
  REAL*8 :: dx2                !< Half x Control volumes size
  REAL*8 :: dy2                !< Half y Control volumes size
  INTEGER :: comp_cells_x      !< Number of control volumes x in the comp. domain
  INTEGER :: comp_interfaces_x !< Number of interfaces (comp_cells_x+1)
  INTEGER :: comp_cells_y      !< Number of control volumes y in the comp. domain
  INTEGER :: comp_interfaces_y !< Number of interfaces (comp_cells_y+1)
  REAL*8 :: cell_size

CONTAINS

  !*********************************************************************
  !> \brief Finite volume grid initialization
  !
  !> This subroutine initialize the grids for the finite volume solver.
  !> \date 16/08/2011
  !*********************************************************************

  SUBROUTINE init_grid

    USE parameters_2d, ONLY: eps_sing, topography_demfile

    IMPLICIT none

    INTEGER j,k      !> loop counter

    comp_interfaces_x = comp_cells_x+1
    comp_interfaces_y = comp_cells_y+1

    ALLOCATE( x_comp(comp_cells_x) )
    ALLOCATE( x_stag(comp_interfaces_x) )
    ALLOCATE( y_comp(comp_cells_y) )
    ALLOCATE( y_stag(comp_interfaces_y) )

    ALLOCATE( B_cent(comp_cells_x,comp_cells_y) )
    ALLOCATE( B_prime_x(comp_cells_x,comp_cells_y) )
    ALLOCATE( B_prime_y(comp_cells_x,comp_cells_y) )
    ALLOCATE( B_ver(comp_interfaces_x,comp_interfaces_y) )
    ALLOCATE( B_stag_x(comp_interfaces_x,comp_cells_y) )
    ALLOCATE( B_stag_y(comp_cells_x,comp_interfaces_y) )

    ALLOCATE( grid_output(comp_cells_x,comp_cells_y) )

    ALLOCATE( grav_surf(comp_cells_x,comp_cells_y) )

    IF ( comp_cells_x .GT. 1 ) THEN

       dx = cell_size

    ELSE
       
       dx = 1.D0

    END IF

    
    IF ( comp_cells_y .GT. 1 ) THEN

       dy = cell_size
    
    ELSE

       dy = 1.D0

    END IF

    xN = x0 + comp_cells_x * dx
    yN = y0 + comp_cells_y * dy

    dx2 = dx/2.d0
    dy2 = dy/2.d0

    !eps_sing = MIN( dx ** 4.D0,dy ** 4.D0 )
    eps_sing=MIN(MIN( dx ** 4.D0,dy ** 4.D0 ),1.d-16)

    IF ( verbose_level .GE. 1 ) WRITE(*,*) 'eps_sing = ',eps_sing

    x_comp(1) = x0 + 0.5D0 * dx
    x_stag(1) = x0
    y_comp(1) = y0 + 0.5D0 * dy
    y_stag(1) = y0

    ! if topography is defined in file .inp we do a rescaling
    IF ( .NOT.topography_demfile ) THEN

       topography_profile(1,:,:) = x0 + ( xN - x0 ) * topography_profile(1,:,:)

       topography_profile(2,:,:) = y0 + ( yN - y0 ) * topography_profile(2,:,:)

       B_ver(1,1) = topography_profile(3,1,1)

       ! bottom row
       DO j=1,comp_cells_x

          x_stag(j+1) = x_stag(j) + dx

          IF(k.EQ.comp_cells_y)THEN

             ! right-bottom vertex
             B_ver(j+1,1)=topography_profile(3,n_topography_profile_x,1)

          ELSE

             CALL interp_1d_scalar( topography_profile(1,:,1) ,                 &
                  topography_profile(3,:,1) , x_stag(j+1) , B_ver(j+1,1) ) 

          ENDIF

          B_stag_y(j,1)=0.5*(B_ver(j+1,1)+B_ver(j,1))

       ENDDO

       ! left column
       DO k=1,comp_cells_y

          y_stag(k+1) = y_stag(k) + dy

          IF(k.EQ.comp_cells_y)THEN

             ! left-top vertex
             B_ver(1,k+1)=topography_profile(3,1,n_topography_profile_y)

          ELSE

             CALL interp_1d_scalar( topography_profile(2,1,:) ,                 &
                  topography_profile(3,1,:) , y_stag(k+1) , B_ver(1,k+1) ) 

          ENDIF

          B_stag_x(1,k)=0.5*(B_ver(1,k+1)+B_ver(1,k))

       ENDDO

       ! all the other cells
       DO j=1,comp_cells_x

          !x_stag(j+1) = x_stag(j) + dx

          x_comp(j) = 0.5 * ( x_stag(j) + x_stag(j+1) )

          DO k=1,comp_cells_y

             !y_stag(k+1) = y_stag(k) + dy

             y_comp(k) = 0.5 * ( y_stag(k) + y_stag(k+1) )

             ! right column
             IF ( j.EQ.comp_cells_x .AND. k.NE.comp_cells_y ) THEN

                CALL interp_1d_scalar(                                          &
                     topography_profile(2,n_topography_profile_x,:) ,           &
                     topography_profile(3,n_topography_profile_x,:) ,           &
                     y_stag(k+1) , B_ver(j+1,k+1) )

                ! top row
             ELSEIF ( j.NE.comp_cells_x .AND. k.EQ.comp_cells_y ) THEN

                CALL interp_1d_scalar(                                          &
                     topography_profile(1,:,n_topography_profile_y) ,           &
                     topography_profile(3,:,n_topography_profile_y) ,           &
                     x_stag(j+1) , B_ver(j+1,k+1) ) 

                ! right-top vertex
             ELSEIF ( j.EQ.comp_cells_x .AND. k.EQ.comp_cells_y ) THEN

                B_ver(j+1,k+1) = topography_profile( 3, n_topography_profile_x ,&
                     n_topography_profile_y)

                ! internal cells
             ELSE

                CALL interp_2d_scalar( topography_profile(1,:,:) ,              &
                     & topography_profile(2,:,:), topography_profile(3,:,:) ,   &
                     & x_stag(j+1), y_stag(k+1) , B_ver(j+1,k+1) )

             ENDIF

             B_stag_x(j+1,k)=0.5*(B_ver(j+1,k+1)+B_ver(j+1,k))
             B_stag_y(j,k+1)=0.5*(B_ver(j+1,k+1)+B_ver(j,k+1))

             B_cent(j,k) = 0.25 * ( B_stag_x(j,k) + B_stag_x(j+1,k)             &
                  + B_stag_y(j,k) + B_stag_y(j,k+1) )

             B_prime_x(j,k) = ( B_stag_x(j+1,k) - B_stag_x(j,k) ) /             &
                  (  x_stag(j+1) - x_stag(j) )

             B_prime_y(j,k) = ( B_stag_y(j,k+1) - B_stag_y(j,k) ) /             &
                  (  y_stag(k+1) - y_stag(k) )

             IF ( verbose_level .GE. 2 ) THEN

                WRITE(*,*) topography_profile(1,:,:) 
                WRITE(*,*) topography_profile(2,:,:)
                WRITE(*,*) topography_profile(3,:,:)
                WRITE(*,*) x_stag(j+1) , y_stag(k+1) , B_stag_x(j+1,k) ,        &
                     B_stag_y(j,k+1) , x_comp(j) , y_comp(k) , B_cent(j,k) ,    &
                     B_ver(j+1,k+1) ,  B_prime_x(j,k) , B_prime_y(j,k) 
                READ(*,*)

             END IF

          END DO

       ENDDO

    ! topography from larger dem
    ELSE

       DO j=0,comp_cells_x

          IF ( j.GT.0 ) THEN

             x_stag(j+1) = x_stag(j) + dx

             x_comp(j) = 0.5 * ( x_stag(j) + x_stag(j+1) )

          ENDIF

          DO k=0,comp_cells_y

             IF(k.GT.0)THEN

                y_stag(k+1) = y_stag(k) + dy

                y_comp(k) = 0.5D0 * ( y_stag(k) + y_stag(k+1) )

             ENDIF

             CALL interp_2d_scalar( topography_profile(1,:,:) ,                 &
                  topography_profile(2,:,:), topography_profile(3,:,:) ,        &
                  x_stag(j+1), y_stag(k+1) , B_ver(j+1,k+1) )

             IF(k.GT.0)THEN

                B_stag_x(j+1,k)=0.5D0 * (B_ver(j+1,k+1)+B_ver(j+1,k))

             ENDIF

             IF(j.GT.0)THEN

                B_stag_y(j,k+1)=0.5D0 * (B_ver(j+1,k+1)+B_ver(j,k+1))

             ENDIF

             IF( ( j.GT.0 ).AND.( k.GT.0 ) ) THEN

                B_cent(j,k) = 0.25D0 * ( B_stag_x(j,k) + B_stag_x(j+1,k) +      &
                     B_stag_y(j,k) + B_stag_y(j,k+1) )

                B_prime_x(j,k) = ( B_stag_x(j+1,k) - B_stag_x(j,k) ) /          &
                     (  x_stag(j+1) - x_stag(j) )
                B_prime_y(j,k) = ( B_stag_y(j,k+1) - B_stag_y(j,k) ) /          &
                     (  y_stag(k+1) - y_stag(k) )

             ENDIF

          END DO

       ENDDO

    ENDIF

    ! we calculate g components wrt surface reference system Tx, Ty and N
    DO j = 1,comp_cells_x

       DO k=1,comp_cells_y

          grav_surf(j,k) = - ( 1.d0/ DSQRT( 1.d0 + B_prime_x(j,k)**2            & 
               + B_prime_y(j,k)**2) )

       ENDDO

    ENDDO

  END SUBROUTINE init_grid

  !---------------------------------------------------------------------------
  !> Scalar interpolation
  !
  !> This subroutine interpolate the values of the  array f1, defined on the 
  !> grid points x1, at the point x2. The value are saved in f2
  !> \date 13/02/2009
  !> \param    x1           original grid                (\b input)
  !> \param    f1           original values              (\b input)
  !> \param    x2           new point                    (\b output)
  !> \param    f2           interpolated value           (\b output)
  !---------------------------------------------------------------------------

  SUBROUTINE interp_1d_scalar(x1, f1, x2, f2)
    IMPLICIT NONE

    REAL*8, INTENT(IN), DIMENSION(:) :: x1, f1
    REAL*8, INTENT(IN) :: x2
    REAL*8, INTENT(OUT) :: f2
    INTEGER :: n, n1x, t
    REAL*8 :: grad , rel_pos

    n1x = SIZE(x1)

    !
    ! ... locate the grid points near the topographic points
    ! ... and interpolate linearly the profile
    !
    t = 1

    search:DO n = 1, n1x-1

       rel_pos = ( x2 - x1(n) ) / ( x1(n+1) - x1(n) )

       IF ( ( rel_pos .GE. 0.D0 ) .AND. ( rel_pos .LE. 1.D0 ) ) THEN

          grad = ( f1(n+1)-f1(n) ) / ( x1(n+1)-x1(n) )
          f2 = f1(n) + ( x2-x1(n) ) * grad

          EXIT search

       ELSEIF  ( rel_pos .LT. 0.D0 ) THEN

          f2 = f1(n)

       ELSE

          f2 = f1(n+1)

       END IF

    END DO search

    RETURN

  END SUBROUTINE interp_1d_scalar

  !---------------------------------------------------------------------------
  !> Scalar interpolation (2D)
  !
  !> This subroutine interpolate the values of the  array f1, defined on the 
  !> grid points (x1,y1), at the point (x2,y2). The value are saved in f2
  !> \date OCTOBER 2016
  !> \param    x1           original grid                (\b input)
  !> \param    y1           original grid                (\b input)
  !> \param    f1           original values              (\b input)
  !> \param    x2           new point                    (\b output)
  !> \param    y2           new point                    (\b output)
  !> \param    f2           interpolated value           (\b output)
  !---------------------------------------------------------------------------

  SUBROUTINE interp_2d_scalar(x1, y1, f1, x2, y2, f2)
    IMPLICIT NONE

    REAL*8, INTENT(IN), DIMENSION(:,:) :: x1, y1, f1
    REAL*8, INTENT(IN) :: x2, y2
    REAL*8, INTENT(OUT) :: f2

    INTEGER :: n, m, n1x, n2x
    REAL*8 :: p,q

    n1x = SIZE(x1,1)
    n2x = SIZE(x1,2)

    IF ( n1x .EQ. 1 ) THEN 

       n = 1

    ELSE

       n=2
       
       DO WHILE (x2.GT.x1(n,1))
          
          n=n+1
          
       ENDDO
       
    END IF

    IF ( n2x .EQ. 1 ) THEN

       m = 1

    ELSE

       m=2
       
       DO WHILE (y2.GT.y1(1,m))
          
          m=m+1
          
       ENDDO

    END IF
       
    IF ( n1x .EQ. 1 ) THEN

       q = ( y2-y1(n,m-1) ) / ( y1(n,m) - y1(n,m-1) )
       
       f2 = ( 1.d0 - q ) * f1(1,m-1) + q * f1(1,m)
       
    ELSEIF ( n2x .EQ. 1 ) THEN

       p = ( x2-x1(n-1,m) ) / ( x1(n,m) - x1(n-1,m) )
       
       f2 = ( 1.d0 - p ) * f1(n-1,1) + p * f1(n,1)
       
    ELSE
      
       p = ( x2-x1(n-1,m) ) / ( x1(n,m) - x1(n-1,m) )
       
       q = ( y2-y1(n,m-1) ) / ( y1(n,m) - y1(n,m-1) )
       
       f2 = ( 1.d0 - q ) * ( ( 1.d0 - p ) * f1(n-1,m-1) + p * f1(n,m-1) ) +        &
            q * ( ( 1.d0 - p ) * f1(n-1,m) + p * f1(n,m) )
       
    END IF

  END SUBROUTINE interp_2d_scalar


  !------------------------------------------------------------------------------
  !> Topography function
  !
  !> This subroutine generates a point of the topography
  !> from the input (x,y) grid point
  !> \date OCTOBER 2016
  !> \param    x           original grid                (\b input)
  !> \param    y           original grid                (\b input)
  !------------------------------------------------------------------------------
  REAL*8 FUNCTION topography_function(x,y)
    IMPLICIT NONE

    REAL*8, INTENT(IN) :: x,y

    REAL*8, PARAMETER :: pig = 4.0*ATAN(1.0)
    REAL*8, PARAMETER :: eps_dis = 1.d-8
    REAL*8 :: a

    ! example 1D from Kurganov and Petrova 2007    
    !IF(y.LT.0.0)THEN
    !
    !  topography_function = 1.d0
    !
    !ELSEIF(y.GE.0.0.AND.y.LE.0.4)THEN
    !
    !  topography_function = COS(pig*y)**2
    !
    !ELSEIF(y.GT.0.4.AND.y.LE.0.5)THEN 
    !
    !  topography_function = COS(pig*y)**2+0.25*(COS(10.0*pig*(y-0.5))+1)
    !
    !ELSEIF(y.GT.0.5.AND.y.LE.0.6)THEN
    !
    !  topography_function = 0.5*COS(pig*y)**4+0.25*(COS(10.0*pig*(y-0.5))+1)
    !
    !ELSEIF(y.GT.0.6.AND.y.LT.1.0-eps_dis)THEN
    !
    !  topography_function = 0.5*COS(pig*y)**4
    !
    !ELSEIF(y.GE.1.0-eps_dis.AND.y.LE.1.0+eps_dis)THEN
    !
    !  topography_function = 0.25
    !
    !ELSEIF(y.GT.1.0+eps_dis.AND.y.LE.1.5)THEN
    !
    !  topography_function = 0.25*SIN(2*pig*(y-1))
    !
    !ELSE
    !
    !  topography_function = 0.d0
    !
    !ENDIF


    ! example 2D from Kurganov and Petrova 2007    
    IF(ABS(y).LE.0.5.AND.x.LE.(y-1.0)/2.0)THEN

       a=y**2

    ELSEIF(ABS(y).GT.0.5.AND.x.LE.(y-1.0)/2.0)THEN

       a=y**2+0.1*sin(pig*x)

    ELSE

       a=MAX(0.125,y**2+0.1*sin(pig*x))

    ENDIF

    topography_function=7.0/32.0*exp(-8.0*(x-0.3)**2-60.0*(y-0.1)**2)- &
         & 1.0/8.0*exp(-30.0*(x+0.1)**2-90.0*(y+0.2)**2) + a

  END FUNCTION topography_function


END MODULE geometry_2d
