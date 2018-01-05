!*******************************************************************************
!> \brief Numerical solver
!
!> This module contains the variables and the subroutines for the 
!> numerical solution of the equations.  
!
!> \date 07/10/2016
!> @author 
!> Mattia de' Michieli Vitturi
!
!********************************************************************************
MODULE solver_2d

  ! external variables

  USE constitutive_2d, ONLY : implicit_flag

  USE geometry_2d, ONLY : comp_cells_x,comp_cells_y
  USE geometry_2d, ONLY : comp_interfaces_x,comp_interfaces_y

  USE geometry_2d, ONLY : B_cent , B_prime_x , B_prime_y , B_stag_x , B_stag_y
  USE geometry_2d, ONLY : grav_surf

  USE parameters_2d, ONLY : n_eqns , n_vars , n_nh
  USE parameters_2d, ONLY : n_RK
  USE parameters_2d, ONLY : verbose_level

  USE parameters_2d, ONLY : bcW , bcE , bcS , bcN

  IMPLICIT none

  !> Conservative variables
  REAL*8, ALLOCATABLE :: q(:,:,:)        
  !> Conservative variables
  REAL*8, ALLOCATABLE :: q0(:,:,:)        
  !> Solution of the finite-volume semidiscrete cheme
  REAL*8, ALLOCATABLE :: q_fv(:,:,:)     
  !> Reconstructed value at the left of the x-interface
  REAL*8, ALLOCATABLE :: q_interfaceL(:,:,:)        
  !> Reconstructed value at the right of the x-interface
  REAL*8, ALLOCATABLE :: q_interfaceR(:,:,:)
  !> Reconstructed value at the bottom of the y-interface
  REAL*8, ALLOCATABLE :: q_interfaceB(:,:,:)        
  !> Reconstructed value at the top of the y-interface
  REAL*8, ALLOCATABLE :: q_interfaceT(:,:,:)
  !> Local speeds at the left of the x-interface
  REAL*8, ALLOCATABLE :: a_interface_xNeg(:,:,:)
  !> Local speeds at the right of the x-interface
  REAL*8, ALLOCATABLE :: a_interface_xPos(:,:,:)
  !> Local speeds at the bottom of the y-interface
  REAL*8, ALLOCATABLE :: a_interface_yNeg(:,:,:)
  !> Local speeds at the top of the y-interface
  REAL*8, ALLOCATABLE :: a_interface_yPos(:,:,:)
  !> Semidiscrete numerical interface fluxes 
  REAL*8, ALLOCATABLE :: H_interface_x(:,:,:)
  !> Semidiscrete numerical interface fluxes 
  REAL*8, ALLOCATABLE :: H_interface_y(:,:,:)
  !> Physical variables (\f$\alpha_1, p_1, p_2, \rho u, w, T\f$)
  REAL*8, ALLOCATABLE :: qp(:,:,:)

  !> Array defining fraction of cells affected by source term
  REAL*8, ALLOCATABLE :: source_xy(:,:)

  LOGICAL, ALLOCATABLE :: solve_mask(:,:) , solve_mask0(:,:)

  !> Time step
  REAL*8 :: dt

  LOGICAL, ALLOCATABLE :: mask22(:,:) , mask21(:,:) , mask11(:,:) , mask12(:,:)

  !> Butcher Tableau for the explicit part of the Runge-Kutta scheme
  REAL*8, ALLOCATABLE :: a_tilde_ij(:,:)
  !> Butcher Tableau for the implicit part of the Runge-Kutta scheme
  REAL*8, ALLOCATABLE :: a_dirk_ij(:,:)

  !> Coefficients for the explicit part of the Runge-Kutta scheme
  REAL*8, ALLOCATABLE :: omega_tilde(:)
  !> Coefficients for the implicit part of the Runge-Kutta scheme
  REAL*8, ALLOCATABLE :: omega(:)

  !> Explicit coeff. for the hyperbolic part for a single step of the R-K scheme
  REAL*8, ALLOCATABLE :: a_tilde(:)
  !> Explicit coeff. for the non-hyp. part for a single step of the R-K scheme
  REAL*8, ALLOCATABLE :: a_dirk(:)
  !> Implicit coeff. for the non-hyp. part for a single step of the R-K scheme
  REAL*8 :: a_diag

  !> Intermediate solutions of the Runge-Kutta scheme
  REAL*8, ALLOCATABLE :: q_rk(:,:,:,:)
  !> Intermediate hyperbolic terms of the Runge-Kutta scheme
  REAL*8, ALLOCATABLE :: divFlux(:,:,:,:)
  !> Intermediate non-hyperbolic terms of the Runge-Kutta scheme
  REAL*8, ALLOCATABLE :: NH(:,:,:,:)

  !> Intermediate explicit terms of the Runge-Kutta scheme
  REAL*8, ALLOCATABLE :: expl_terms(:,:,:,:)

  !> Local Intermediate hyperbolic terms of the Runge-Kutta scheme
  REAL*8, ALLOCATABLE :: divFluxj(:,:)
  !> Local Intermediate non-hyperbolic terms of the Runge-Kutta scheme
  REAL*8, ALLOCATABLE :: NHj(:,:)
  !> Local Intermediate explicit terms of the Runge-Kutta scheme
  REAL*8, ALLOCATABLE :: expl_terms_j(:,:)

  !> Flag for the normalization of the array q in the implicit solution scheme
  LOGICAL :: normalize_q

  !> Flag for the normalization of the array f in the implicit solution scheme
  LOGICAL :: normalize_f

  !> Flag for the search of optimal step size in the implicit solution scheme
  LOGICAL :: opt_search_NL

  !> Sum of all the terms of the equations except the transient term
  REAL*8, ALLOCATABLE :: residual_term(:,:,:)


CONTAINS

  !******************************************************************************
  !> \brief Memory allocation
  !
  !> This subroutine allocate the memory for the variables of the 
  !> solver module.
  !
  !> \date 07/10/2016
  !> @author 
  !> Mattia de' Michieli Vitturi
  !
  !******************************************************************************

  SUBROUTINE allocate_solver_variables

    IMPLICIT NONE

    REAL*8 :: gamma , delta

    INTEGER :: i,j

    ALLOCATE( q( n_vars , comp_cells_x , comp_cells_y ) , q0( n_vars ,          &
         comp_cells_x , comp_cells_y ) )

    ALLOCATE( q_fv( n_vars , comp_cells_x , comp_cells_y ) )

    ALLOCATE( q_interfaceL( n_vars , comp_interfaces_x, comp_cells_y ) )
    ALLOCATE( q_interfaceR( n_vars , comp_interfaces_x, comp_cells_y ) )
    ALLOCATE( a_interface_xNeg( n_eqns , comp_interfaces_x, comp_cells_y ) )
    ALLOCATE( a_interface_xPos( n_eqns , comp_interfaces_x, comp_cells_y ) )
    ALLOCATE( H_interface_x( n_eqns , comp_interfaces_x, comp_cells_y ) )


    ALLOCATE( q_interfaceB( n_vars , comp_cells_x, comp_interfaces_y ) )
    ALLOCATE( q_interfaceT( n_vars , comp_cells_x, comp_interfaces_y ) )
    ALLOCATE( a_interface_yNeg( n_eqns , comp_cells_x, comp_interfaces_y ) )
    ALLOCATE( a_interface_yPos( n_eqns , comp_cells_x, comp_interfaces_y ) )


    ALLOCATE( H_interface_y( n_eqns , comp_cells_x, comp_interfaces_y ) )

    ALLOCATE( solve_mask( comp_cells_x , comp_cells_y ) )
    ALLOCATE( solve_mask0( comp_cells_x , comp_cells_y ) )

    ALLOCATE( qp( n_vars , comp_cells_x , comp_cells_y ) )

    ALLOCATE( source_xy( comp_cells_x , comp_cells_y ) )

    
    ALLOCATE( a_tilde_ij(n_RK,n_RK) )
    ALLOCATE( a_dirk_ij(n_RK,n_RK) )
    ALLOCATE( omega_tilde(n_RK) )
    ALLOCATE( omega(n_RK) )


    ! Allocate the logical arrays defining the implicit parts of the system
    ALLOCATE( mask22(n_eqns,n_eqns) )
    ALLOCATE( mask21(n_eqns,n_eqns) )
    ALLOCATE( mask11(n_eqns,n_eqns) )
    ALLOCATE( mask12(n_eqns,n_eqns) )

    ! Initialize the logical arrays with all false (everything is implicit)
    mask11(1:n_eqns,1:n_eqns) = .FALSE.
    mask12(1:n_eqns,1:n_eqns) = .FALSE.
    mask22(1:n_eqns,1:n_eqns) = .FALSE.
    mask21(1:n_eqns,1:n_eqns) = .FALSE.

    ! Set to .TRUE. the elements not corresponding to equations and variables to 
    ! be solved implicitly
    DO i = 1,n_eqns

       DO j = 1,n_eqns

          IF ( .NOT.implicit_flag(i) .AND. .NOT.implicit_flag(j) )              &
               mask11(j,i) = .TRUE.
          IF ( implicit_flag(i) .AND. .NOT.implicit_flag(j) )                   &
               mask12(j,i) = .TRUE.
          IF ( implicit_flag(i) .AND. implicit_flag(j) )                        &
               mask22(j,i) = .TRUE.
          IF ( .NOT.implicit_flag(i) .AND. implicit_flag(j) )                   &
               mask21(j,i) = .TRUE.

       END DO

    END DO

    ! Initialize the coefficients for the IMEX Runge-Kutta scheme
    ! Please note that with respect to the schemes described in Pareschi & Russo 
    ! (2000) we do not have the coefficient vectors c_tilde and c, because the 
    ! explicit and implicit terms do not depend explicitly on time.

    ! Explicit part coefficients (a_tilde_ij=0 for j>=i)
    a_tilde_ij = 0.D0

    ! Weight coefficients of the explicit part in the final assemblage
    omega_tilde = 0.D0

    ! Implicit part coefficients (a_dirk_ij=0 for j>i)
    a_dirk_ij = 0.D0

    ! Weight coefficients of the explicit part in the final assemblage
    omega = 0.D0

    gamma = 1.D0 - 1.D0 / SQRT(2.D0)
    delta = 1.D0 - 1.D0 / ( 2.D0 * gamma )

    IF ( n_RK .EQ. 1 ) THEN

       a_tilde_ij(1,1) = 1.D0

       omega_tilde(1) = 1.D0

       a_dirk_ij(1,1) = 0.D0

       omega(1) = 0.D0

    ELSEIF ( n_RK .EQ. 2 ) THEN

       a_tilde_ij(2,1) = 1.0D0

       omega_tilde(1) = 1.0D0
       omega_tilde(2) = 0.0D0

       a_dirk_ij(2,2) = 1.0D0

       omega(1) = 0.D0
       omega(2) = 1.D0

      
    ELSEIF ( n_RK .EQ. 3 ) THEN

       ! Tableau for the IMEX-SSP(3,3,2) Stiffly Accurate Scheme
       ! from Pareschi & Russo (2005), Table IV

       a_tilde_ij(2,1) = 0.5D0
       a_tilde_ij(3,1) = 0.5D0
       a_tilde_ij(3,2) = 0.5D0

       omega_tilde(1) =  1.0D0 / 3.0D0
       omega_tilde(2) =  1.0D0 / 3.0D0
       omega_tilde(3) =  1.0D0 / 3.0D0

       a_dirk_ij(1,1) = 0.25D0
       a_dirk_ij(2,2) = 0.25D0
       a_dirk_ij(3,1) = 1.0D0 / 3.0D0
       a_dirk_ij(3,2) = 1.0D0 / 3.0D0
       a_dirk_ij(3,3) = 1.0D0 / 3.0D0

       omega(1) =  1.0D0 / 3.0D0
       omega(2) =  1.0D0 / 3.0D0
       omega(3) =  1.0D0 / 3.0D0

    ELSEIF ( n_RK .EQ. 4 ) THEN

       ! LRR(3,2,2) from Table 3 in Pareschi & Russo (2000)

       a_tilde_ij(2,1) = 0.5D0
       a_tilde_ij(3,1) = 1.D0 / 3.D0
       a_tilde_ij(4,2) = 1.0D0

       omega_tilde(1) = 0.D0
       omega_tilde(2) = 1.0D0
       omega_tilde(3) = 0.0D0
       omega_tilde(4) = 0.D0

       a_dirk_ij(2,2) = 0.5D0
       a_dirk_ij(3,3) = 1.0D0 / 3.0D0
       a_dirk_ij(4,3) = 0.75D0
       a_dirk_ij(4,4) = 0.25D0

       omega(1) = 0.D0
       omega(2) = 0.D0
       omega(3) = 0.75D0
       omega(4) = 0.25D0

    END IF

    ALLOCATE( a_tilde(n_RK) )
    ALLOCATE( a_dirk(n_RK) )

    ALLOCATE( q_rk( n_vars , comp_cells_x , comp_cells_y , n_RK ) )
    ALLOCATE( divFlux( n_eqns , comp_cells_x , comp_cells_y , n_RK ) )
    ALLOCATE( NH( n_eqns , comp_cells_x , comp_cells_y , n_RK ) )

    ALLOCATE( expl_terms( n_eqns , comp_cells_x , comp_cells_y , n_RK ) )

    ALLOCATE( divFluxj(n_eqns,n_RK) )
    ALLOCATE( NHj(n_eqns,n_RK) )
    ALLOCATE( expl_terms_j(n_eqns,n_RK) )

    ALLOCATE( residual_term( n_vars , comp_cells_x , comp_cells_y ) )

  END SUBROUTINE allocate_solver_variables

  !******************************************************************************
  !> \brief Memory deallocation
  !
  !> This subroutine de-allocate the memory for the variables of the 
  !> solver module.
  !
  !> \date 07/10/2016
  !> @author 
  !> Mattia de' Michieli Vitturi
  !
  !******************************************************************************

  SUBROUTINE deallocate_solver_variables

    DEALLOCATE( q , q0 )

    DEALLOCATE( q_fv )

    DEALLOCATE( q_interfaceL )
    DEALLOCATE( q_interfaceR )
    DEALLOCATE( q_interfaceB )
    DEALLOCATE( q_interfaceT )

    DEALLOCATE( a_interface_xNeg )
    DEALLOCATE( a_interface_xPos )
    DEALLOCATE( a_interface_yNeg )
    DEALLOCATE( a_interface_yPos )

    DEALLOCATE( H_interface_x )
    DEALLOCATE( H_interface_y )

    DEALLOCATE( solve_mask , solve_mask0 )

    Deallocate( qp )

    DEALLOCATE( source_xy )
    
    DEALLOCATE( a_tilde_ij )
    DEALLOCATE( a_dirk_ij )
    DEALLOCATE( omega_tilde )
    DEALLOCATE( omega )

    DEALLOCATE( implicit_flag )

    DEALLOCATE( a_tilde )
    DEALLOCATE( a_dirk )

    DEALLOCATE( q_rk )
    DEALLOCATE( divFlux )
    DEALLOCATE( NH )
    DEALLOCATE( expl_terms )

    DEALLOCATE( divFluxj )
    DEALLOCATE( NHj )
    DEALLOCATE( expl_terms_j )

    DEALLOCATE( mask22 , mask21 , mask11 , mask12 )

    DEALLOCATE( residual_term )

  END SUBROUTINE deallocate_solver_variables


  !******************************************************************************
  !> \brief Masking of cells to solve
  !
  !> This subroutine compute a 2D array of logicals defining the cells where the
  !> systems of equations have to be solved. It is defined according to the 
  !> positive thickness in the cell and in the neighbour cells
  !
  !> \date 20/04/2017
  !> @author 
  !> Mattia de' Michieli Vitturi
  !
  !******************************************************************************

  SUBROUTINE check_solve

    IMPLICIT NONE

    INTEGER :: i

    solve_mask0(1:comp_cells_x,1:comp_cells_y) = .FALSE.

    solve_mask = solve_mask0

    WHERE ( q(1,:,:) - B_cent(:,:) .GT. 0.D0 ) solve_mask = .TRUE.

    DO i = 1,n_RK

       solve_mask(1+i:comp_cells_x,:) =  solve_mask(1+i:comp_cells_x,:) .OR.    &
            solve_mask(1:comp_cells_x-i,:) 

       solve_mask(1:comp_cells_x-i,:) =  solve_mask(1:comp_cells_x-i,:) .OR.    &
            solve_mask(1+i:comp_cells_x,:) 

       solve_mask(:,1+i:comp_cells_y) =  solve_mask(:,1+i:comp_cells_y) .OR.    &
            solve_mask(:,1:comp_cells_y-i) 

       solve_mask(:,1:comp_cells_y-i) =  solve_mask(:,1:comp_cells_y-i) .OR.    &
            solve_mask(:,1+i:comp_cells_y) 

    END DO


  END SUBROUTINE check_solve

  !******************************************************************************
  !> \brief Time-step computation
  !
  !> This subroutine evaluate the maximum time step according to the CFL
  !> condition. The local speed are evaluated with the characteristic
  !> polynomial of the Jacobian of the fluxes.
  !
  !> \date 07/10/2016
  !> @author 
  !> Mattia de' Michieli Vitturi
  !
  !******************************************************************************

  SUBROUTINE timestep

    ! External variables
    USE geometry_2d, ONLY : dx,dy
    USE parameters_2d, ONLY : max_dt , cfl

    ! External procedures
    USE constitutive_2d, ONLY : eval_local_speeds_x, eval_local_speeds_y

    IMPLICIT none

    REAL*8 :: vel_max(n_vars)
    REAL*8 :: vel_min(n_vars)
    REAL*8 :: vel_j         !< maximum speed in the j-th cell
    REAL*8 :: dt_cfl        !< local time step
    REAL*8 :: qj(n_vars)    !< conservative variables

    INTEGER :: j,k          !< loop counter

    dt = max_dt

    IF ( cfl .NE. -1.d0 ) THEN

       DO j = 1,comp_cells_x

          DO k = 1,comp_cells_y

             qj = q( 1:n_vars , j , k )

             IF ( comp_cells_x .GT. 1 ) THEN

                ! x direction
                CALL eval_local_speeds_x( qj , B_cent(j,k) , vel_min , vel_max )
                
                vel_j = MAX( MAXVAL(ABS(vel_min)) , MAXVAL(ABS(vel_max)) )
                
                dt_cfl = cfl * dx / vel_j
                
                dt = MIN( dt , dt_cfl )

             END IF

             IF ( comp_cells_y .GT. 1 ) THEN
                
                ! y direction
                CALL eval_local_speeds_y( qj , B_cent(j,k) , vel_min , vel_max )
                
                vel_j = MAX( MAXVAL(ABS(vel_min)) , MAXVAL(ABS(vel_max)) )
                
                dt_cfl = cfl * dy / vel_j
                                
                dt = MIN( dt , dt_cfl )

             END IF
                
          ENDDO

       END DO

    END IF

  END SUBROUTINE timestep


  !*****************************************************************************
  !> \brief Time-step computation
  !
  !> This subroutine evaluate the maximum time step according to the CFL
  !> condition. The local speed are evaluated with the characteristic
  !> polynomial of the Jacobian of the fluxes.
  !
  !> \date 07/10/2016
  !> @author 
  !> Mattia de' Michieli Vitturi
  !
  !*****************************************************************************

  SUBROUTINE timestep2

    ! External variables
    USE geometry_2d, ONLY : dx,dy
    USE parameters_2d, ONLY : max_dt , cfl

    ! External procedures
    USE constitutive_2d, ONLY : eval_local_speeds2_x, eval_local_speeds2_y

    IMPLICIT none

    REAL*8 :: dt_cfl        !< local time step

    REAL*8 :: a_interface_x_max(n_eqns,comp_interfaces_x,comp_cells_y)
    REAL*8 :: a_interface_y_max(n_eqns,comp_cells_x,comp_interfaces_y)
    REAL*8 :: dt_interface_x, dt_interface_y

    INTEGER :: i,j,k          !< loop counter

    dt = max_dt

    IF ( cfl .NE. -1.d0 ) THEN

       CALL reconstruction

       CALL eval_speeds

       DO i=1,n_vars

          a_interface_x_max(i,:,:) =                                            &
               MAX( a_interface_xPos(i,:,:) , -a_interface_xNeg(i,:,:) )

          a_interface_y_max(i,:,:) =                                            &
               MAX( a_interface_yPos(i,:,:) , -a_interface_yNeg(i,:,:) )

       END DO

       DO j = 1,comp_cells_x

          DO k = 1,comp_cells_y

             dt_interface_x = cfl * dx / MAX( MAXVAL(a_interface_x_max(:,j,k))  &
                  , MAXVAL(a_interface_x_max(:,j+1,k)) )

             dt_interface_y = cfl * dy / MAX( MAXVAL(a_interface_y_max(:,j,k))  &
                  , MAXVAL(a_interface_y_max(:,j,k+1)) )

             dt_cfl = MIN( dt_interface_x , dt_interface_y )

             dt = MIN(dt,dt_cfl)

          ENDDO

       END DO

    END IF

  END SUBROUTINE timestep2

  !******************************************************************************
  !> \brief Runge-Kutta integration
  !
  !> This subroutine integrate the hyperbolic conservation law with
  !> non-hyperbolic terms using an implicit-explicit runge-kutta scheme.
  !> The fluxes are integrated explicitely while the non-hyperbolic terms
  !> are integrated implicitely.
  !
  !> \date 07/10/2016
  !> @author 
  !> Mattia de' Michieli Vitturi
  !
  !******************************************************************************

  SUBROUTINE imex_RK_solver

    USE constitutive_2d, ONLY : eval_nonhyperbolic_terms

    USE constitutive_2d, ONLY : qc_to_qp

    USE geometry_2d, ONLY : dx

    IMPLICIT NONE

    REAL*8 :: q_guess(n_vars) !< initial guess for the solution of the RK step
    INTEGER :: i_RK           !< loop counter for the RK iteration
    INTEGER :: j,k            !< loop counter over the grid volumes

    REAL*8 :: h_new

    ! Initialization of the solution guess
    q0( 1:n_vars , 1:comp_cells_x , 1:comp_cells_y ) =                          &
         q( 1:n_vars , 1:comp_cells_x , 1:comp_cells_y )

    IF ( verbose_level .GE. 2 ) WRITE(*,*) 'solver, imex_RK_solver: beginning'

    ! Initialization of the variables for the Runge-Kutta scheme
    q_rk(1:n_vars,1:comp_cells_x,1:comp_cells_y,1:n_RK) = 0.d0

    divFlux(1:n_eqns,1:comp_cells_x,1:comp_cells_y,1:n_RK) = 0.d0

    NH(1:n_eqns,1:comp_cells_x,1:comp_cells_y,1:n_RK) = 0.d0

    expl_terms(1:n_eqns,1:comp_cells_x,1:comp_cells_y,1:n_RK) = 0.d0

    runge_kutta:DO i_RK = 1,n_RK

       IF ( verbose_level .GE. 2 ) WRITE(*,*) 'solver, imex_RK_solver: i_RK',i_RK

       ! define the explicits coefficients for the i-th step of the Runge-Kutta
       a_tilde = 0.d0
       a_dirk = 0.d0

       ! in the first step of the RK scheme all the coefficients remain to 0
       a_tilde(1:i_RK-1) = a_tilde_ij(i_RK,1:i_RK-1)
       a_dirk(1:i_RK-1) = a_dirk_ij(i_RK,1:i_RK-1)

       ! define the implicit coefficient for the i-th step of the Runge-Kutta
       a_diag = a_dirk_ij(i_RK,i_RK)

       loop_over_xcells:DO j = 1,comp_cells_x

          loop_over_ycells:DO k = 1,comp_cells_y

             IF ( verbose_level .GE. 2 ) THEN

                WRITE(*,*) 'solver, imex_RK_solver: j',j,k

             END IF

             ! initialize the RK step
             IF ( i_RK .EQ. 1 ) THEN

                ! solution from the previous time step
                q_guess(1:n_vars) = q0( 1:n_vars , j , k) 

             ELSE

                ! solution from the previous RK step
                q_guess(1:n_vars) = q_rk( 1:n_vars , j , k , MAX(1,i_RK-1) )

             END IF

             ! RK explicit hyperbolic terms for the volume (i,j)
             divFluxj(1:n_eqns,1:n_RK) = divFlux( 1:n_eqns , j , k , 1:n_RK )

             ! RK implicit terms for the volume (i,j)
             NHj(1:n_eqns,1:n_RK) = NH( 1:n_eqns , j , k , 1:n_RK )

             ! RK additional explicit terms for the volume (i,j)
             Expl_terms_j(1:n_eqns,1:n_RK) = expl_terms( 1:n_eqns,j,k,1:n_RK )

             ! New solution at the i_RK step without the implicit term
             q_fv( 1:n_vars , j , k ) = q0( 1:n_vars , j , k )                  &
                  - dt * (MATMUL( divFluxj(1:n_eqns,1:i_RK)                     &
                  + Expl_terms_j(1:n_eqns,1:i_RK) , a_tilde(1:i_RK) )           &
                  - MATMUL( NHj(1:n_eqns,1:i_RK) , a_dirk(1:i_RK) ) )

            
             IF ( verbose_level .GE. 2 ) THEN

                WRITE(*,*) 'q_guess',q_guess
                CALL qc_to_qp( q_guess , B_cent(j,k) , qp(1:n_vars,j,k) )
                WRITE(*,*) 'q_guess: qp',qp(1:n_vars,j,k)

             END IF

             IF ( a_diag .NE. 0.D0 ) THEN

                IF ( q_fv(1,j,k) .GT.  B_cent(j,k) )  THEN

                   ! Solve the implicit system to find the solution at the 
                   ! i_RK step of the IMEX RK procedure
                   CALL solve_rk_step( B_cent(j,k) , B_prime_x(j,k) ,           &
                        B_prime_y(j,k), grav_surf(j,k) ,                        &
                        q_guess , q0(1:n_vars,j,k ) , a_tilde ,                 &
                        a_dirk , a_diag )

                ELSE

                   q_guess(1:n_vars) = q_fv( 1:n_vars , j , k ) 

                END IF

             END IF

             ! Eval and store the implicit term at the i_RK step
             CALL eval_nonhyperbolic_terms( B_cent(j,k) , B_prime_x(j,k) ,      &
                  B_prime_y(j,k) , grav_surf(j,k) , r_qj = q_guess ,            &
                  r_nh_term_impl = NH(1:n_eqns,j,k,i_RK) ) 

             IF ( q_fv(2,j,k)**2 + q_fv(3,j,k)**2 .EQ. 0.D0 ) THEN

                q_guess(2:3) = 0.D0 

             ! Check if the impl. friction term changed the sign of the velocity
             ELSEIF ( ( q_guess(2)*q_fv(2,j,k) .LT. 0.D0 ) .AND.                &
                  ( q_guess(3)*q_fv(3,j,k) .LT. 0.D0 ) ) THEN

                q_guess(2:3) = 0.D0 
          
             ! Align the velocity vector with previous one
             ELSE

                q_guess(2:3) = DSQRT( q_guess(2)**2 + q_guess(3)**2 ) *         &
                     q_fv(2:3,j,k) / DSQRT( q_fv(2,j,k)**2 + q_fv(3,j,k)**2 ) 

             END IF
          
             IF ( a_diag .NE. 0.D0 ) THEN
                
                ! Update the viscous term with the correction on the new velocity
                NH(2:3,j,k,i_RK) = - ( q_fv(2:3,j,k) - q_guess(2:3) ) /         &
                     ( dt*a_diag )

             END IF

             ! Store the solution at the end of the i_RK step
             q_rk( 1:n_vars , j , k , i_RK ) = q_guess

             ! Check the sign of the flow thickness
             h_new = q_guess(1) - B_cent(j,k)

             IF ( h_new .LT. 0.D0 ) THEN

                IF ( DABS(h_new) .LT. 1.D-8 ) THEN

                   q_rk(1,j,k,i_RK) = B_cent(j,k)
                   q_rk(2:n_vars,j,k,i_RK) = 0.D0
                   
                ELSE

                   IF ( i_RK .EQ. 1 ) THEN
                      
                      WRITE(*,*) 'q0', q0( 1 , j , k) 
                      
                   ELSE
                      
                      WRITE(*,*) 'q0', q_rk( 1 , j , k , MAX(1,i_RK-1) )
                      
                   END IF
                   
                   CALL qc_to_qp(  q_rk( 1:n_vars , j , k , MAX(1,i_RK-1) ) ,   &
                        B_cent(j,k) , qp(1:n_vars,j,k) )

                   WRITE(*,*) 'q1',q_guess(1)
                   WRITE(*,*) 
                   WRITE(*,*) 'j,k,h',j,k,h_new,qp(1,j,k),B_cent(j,k)
                   WRITE(*,*) 
                   WRITE(*,*) 'i_RK,dt',i_RK,dt
                   WRITE(*,*) 
                   WRITE(*,*) 'divFluxj(1,1:n_RK)',divFluxj(1,1:n_RK)
                   WRITE(*,*) 'NHj(1,1:n_RK)',NHj(1,1:n_RK)
                   WRITE(*,*) 'Expl_terms_j(1,1:n_RK)',Expl_terms_j(1,1:n_RK)
                   WRITE(*,*) 
                   
                   
                   WRITE(*,*) 'h^-_j-1/2',q_interfaceL(1,j,k) - B_stag_x(j,k)
                   WRITE(*,*) 'h^+_j-1/2',q_interfaceR(1,j,k) - B_stag_x(j,k)
                   WRITE(*,*) 'h^-_j+1/2',q_interfaceL(1,j+1,k) - B_stag_x(j+1,k)
                   WRITE(*,*) 'h^+_j+1/2',q_interfaceR(1,j+1,k) - B_stag_x(j+1,k)
                   
                   WRITE(*,*)
                   
                   WRITE(*,*) 'hu^-_j-1/2',q_interfaceL(2,j,k)
                   WRITE(*,*) 'hu^+_j-1/2',q_interfaceR(2,j,k)
                   WRITE(*,*) 'hu^-_j+1/2',q_interfaceL(2,j+1,k)
                   WRITE(*,*) 'hu^+_j+1/2',q_interfaceR(2,j+1,k)
                   
                   WRITE(*,*)
                   WRITE(*,*) 'a^-_j-1/2',a_interface_xNeg(1,j,k)
                   WRITE(*,*) 'a^+_j-1/2',a_interface_xPos(1,j,k)
                   WRITE(*,*) 'a^-_j+1/2',a_interface_xNeg(1,j+1,k)
                   WRITE(*,*) 'a^+_j+1/2',a_interface_xPos(1,j+1,k)

                   WRITE(*,*)
                   WRITE(*,*) 'lambda*a^-_j-1/2' , dt / dx                      &
                        * ABS(a_interface_xNeg(1,j,k))

                   WRITE(*,*) 'lambda*a^+_j-1/2' , dt / dx                      &
                        * ABS(a_interface_xPos(1,j,k))

                   WRITE(*,*) 'lambda*a^-_j+1/2' , dt / dx                      &
                        * ABS(a_interface_xNeg(1,j+1,k))
                   WRITE(*,*) 'lambda*a^-_j+1/2' , dt / dx                      &
                        * ABS(a_interface_xPos(1,j+1,k))

                   WRITE(*,*)
                   
                   WRITE(*,*) 'u^-_j-1/2', q_interfaceL(2,j,k) /                & 
                        ( q_interfaceL(1,j,k) - B_stag_x(j,k) )
 
                   WRITE(*,*) 'u^+_j-1/2', q_interfaceR(2,j,k) /                &
                        ( q_interfaceR(1,j,k) - B_stag_x(j,k) )
                
                   WRITE(*,*) 'u^-_j+1/2', q_interfaceL(2,j+1,k) /              &
                        ( q_interfaceL(1,j+1,k) - B_stag_x(j+1,k) )
                
                   WRITE(*,*) 'u^+_j+1/2', q_interfaceR(2,j+1,k) /              &
                        ( q_interfaceR(1,j+1,k) - B_stag_x(j+1,k) )
                   
                   WRITE(*,*) 
                   
                   WRITE(*,*) 'lambda',dt/dx
                   
                   WRITE(*,*) 
                   
                   WRITE(*,*) ( a_interface_xPos(1,j,k) -  q_interfaceR(1,j,k)  &
                        / ( q_interfaceR(2,j,k) - B_stag_x(j,k) ) ) /           &
                        (  a_interface_xPos(1,j,k) - a_interface_xNeg(1,j,k) )
                   
                   WRITE(*,*) (  q_interfaceL(1,j+1,k) / ( q_interfaceL(2,j+1,k)&
                        -  B_stag_x(j+1,k) ) - a_interface_xNeg(1,j+1,k) ) /    &
                        ( a_interface_xPos(1,j+1,k) - a_interface_xNeg(1,j+1,k) )
                   
                   WRITE(*,*) ( a_interface_xPos(1,j+1,k)                       &
                        - q_interfaceR(1,j+1,k) /                               &
                        ( q_interfaceR(2,j+1,k) - B_stag_x(j+1,k) ) ) /         &
                        (  a_interface_xPos(1,j+1,k) - a_interface_xNeg(1,j+1,k))
                   
                   WRITE(*,*) (  q_interfaceL(1,j,k) / ( q_interfaceL(2,j,k) -  &
                        B_stag_x(j,k) ) - a_interface_xNeg(1,j,k) ) /           &
                        ( a_interface_xPos(1,j,k) - a_interface_xNeg(1,j,k) )
                                      
                   WRITE(*,*) 
                                      
                   WRITE(*,*) dt/dx * q_interfaceR(1,j,k)
                   
                   READ(*,*) 
                   
                END IF
                   
             END IF


             IF ( verbose_level .GE. 2 ) THEN

                WRITE(*,*) 'imex_RK_solver: qc',q_guess
                CALL qc_to_qp( q_guess, B_cent(j,k) , qp(1:n_vars,j,k) )
                WRITE(*,*) 'imex_RK_solver: qp',qp(1:n_vars,j,k)
                READ(*,*)

             END IF

          END DO loop_over_ycells

       ENDDO loop_over_xcells

       ! Eval and store the explicit hyperbolic (fluxes) terms
       CALL eval_hyperbolic_terms( q_rk(1:n_vars,1:comp_cells_x,1:comp_cells_y, &
            i_RK) , divFlux(1:n_eqns,1:comp_cells_x,1:comp_cells_y,i_RK) )

       ! Eval and store the other explicit terms (e.g. gravity or viscous forces)
       CALL eval_explicit_terms( q_rk(1:n_vars,1:comp_cells_x,1:comp_cells_y,   &
            i_RK) , expl_terms(1:n_eqns,1:comp_cells_x,1:comp_cells_y,i_RK) )

       IF ( verbose_level .GE. 1 ) THEN

          WRITE(*,*) 'div_flux(2),div_flux(3),expl_terms(2),expl_terms(3)'

          DO j = 1,comp_cells_x

             DO k = 1,comp_cells_y

                WRITE(*,*) divFlux(2,j,k,i_RK) , divFlux(3,j,k,i_RK) ,          &
                     expl_terms(2,j,k,i_RK) , expl_terms(3,j,k,i_RK)

             ENDDO

          END DO

          READ(*,*)

       END IF

    END DO runge_kutta

    DO j = 1,comp_cells_x

       DO k = 1,comp_cells_y

          residual_term(1:n_vars,j,k) = MATMUL( divFlux(1:n_eqns,j,k,1:n_RK)    &
               + expl_terms(1:n_eqns,j,k,1:n_RK) , omega_tilde )                &
               - MATMUL( NH(1:n_eqns,j,k,1:n_RK) , omega )

       ENDDO

    END DO

    DO j = 1,comp_cells_x

       DO k = 1,comp_cells_y

          IF ( verbose_level .GE. 1 ) THEN

             WRITE(*,*) 'cell jk =',j,k
             WRITE(*,*) 'before imex_RK_solver: qc',q0(1:n_vars,j,k)
             CALL qc_to_qp(q0(1:n_vars,j,k) , B_cent(j,k) , qp(1:n_vars,j,k))
             WRITE(*,*) 'before imex_RK_solver: qp',qp(1:n_vars,j,k)

          END IF

          q(1:n_vars,j,k) = q0(1:n_vars,j,k) - dt * residual_term(1:n_vars,j,k)

          IF ( q(1,j,k) .LT. B_cent(j,k) ) THEN

             IF ( DABS( q(1,j,k)-B_cent(j,k) ) .LT. 1.D-10 ) THEN
                
                q(1,j,k) = B_cent(j,k)
                q(2:n_vars,j,k) = 0.D0
                
             ELSE
                
                WRITE(*,*) 'j,k,n_RK',j,k,n_RK
                WRITE(*,*) 'dt',dt
                
                WRITE(*,*) 'before imex_RK_solver: qc',q0(1:n_vars,j,k)
                CALL qc_to_qp(q0(1:n_vars,j,k) , B_cent(j,k) , qp(1:n_vars,j,k))
                WRITE(*,*) 'before imex_RK_solver: qp',qp(1:n_vars,j,k)
                WRITE(*,*) 'h old',q0(1,j,k) - B_cent(j,k)


                WRITE(*,*) 'h^-_j-1/2',q_interfaceL(1,j,k) - B_stag_x(j,k)
                WRITE(*,*) 'h^+_j-1/2',q_interfaceR(1,j,k) - B_stag_x(j,k)
                WRITE(*,*) 'h^-_j+1/2',q_interfaceL(1,j+1,k) - B_stag_x(j+1,k)
                WRITE(*,*) 'h^+_j+1/2',q_interfaceR(1,j+1,k) - B_stag_x(j+1,k)

                WRITE(*,*)
                
                WRITE(*,*) 'hu^-_j-1/2',q_interfaceL(2,j,k)
                WRITE(*,*) 'hu^+_j-1/2',q_interfaceR(2,j,k)
                WRITE(*,*) 'hu^-_j+1/2',q_interfaceL(2,j+1,k)
                WRITE(*,*) 'hu^+_j+1/2',q_interfaceR(2,j+1,k)

                WRITE(*,*)

                WRITE(*,*) 'a^-_j-1/2',a_interface_xNeg(1,j,k)
                WRITE(*,*) 'a^+_j-1/2',a_interface_xPos(1,j,k)
                WRITE(*,*) 'a^-_j+1/2',a_interface_xNeg(1,j+1,k)
                WRITE(*,*) 'a^+_j+1/2',a_interface_xPos(1,j+1,k)

                WRITE(*,*)

                WRITE(*,*) 'u^-_j-1/2', q_interfaceL(1,j,k) / & 
                     ( q_interfaceL(2,j,k) - B_stag_x(j,k) )
                WRITE(*,*) 'u^+_j-1/2', q_interfaceR(1,j,k) / &
                     ( q_interfaceR(2,j,k) - B_stag_x(j,k) )
                WRITE(*,*) 'u^-_j+1/2', q_interfaceL(1,j+1,k) / &
                     ( q_interfaceL(2,j+1,k) - B_stag_x(j+1,k) )
                WRITE(*,*) 'u^+_j+1/2', q_interfaceR(1,j+1,k) / &
                     ( q_interfaceR(2,j+1,k) - B_stag_x(j+1,k) )

                WRITE(*,*) 
              
                WRITE(*,*) 'lambda',dt/dx

                WRITE(*,*) 

                WRITE(*,*) ( a_interface_xPos(1,j,k) -  q_interfaceR(1,j,k) /   &
                     ( q_interfaceR(2,j,k) - B_stag_x(j,k) ) ) /                &
                     (  a_interface_xPos(1,j,k) - a_interface_xNeg(1,j,k) )

                WRITE(*,*) (  q_interfaceL(1,j+1,k) / ( q_interfaceL(2,j+1,k) - &
                     B_stag_x(j+1,k) ) - a_interface_xNeg(1,j+1,k) ) /          &
                     ( a_interface_xPos(1,j+1,k) - a_interface_xNeg(1,j+1,k) )

                WRITE(*,*) ( a_interface_xPos(1,j+1,k) -  q_interfaceR(1,j+1,k) &
                     / ( q_interfaceR(2,j+1,k) - B_stag_x(j+1,k) ) ) /          &
                     (  a_interface_xPos(1,j+1,k) - a_interface_xNeg(1,j+1,k) )

                WRITE(*,*) (  q_interfaceL(1,j,k) / ( q_interfaceL(2,j,k) -     &
                     B_stag_x(j,k) ) - a_interface_xNeg(1,j,k) ) /              &
                     ( a_interface_xPos(1,j,k) - a_interface_xNeg(1,j,k) )


                

                WRITE(*,*) 
  
                WRITE(*,*) 'divFlux',divFlux(1,j,k,1:n_RK)    
                WRITE(*,*) 'expl_terms',expl_terms(1,j,k,1:n_RK)
                WRITE(*,*) 'NH',NH(1,j,k,1:n_RK)
               
                READ(*,*)
                
             END IF

          END IF

          IF ( verbose_level .GE. 1 ) THEN
          
             CALL qc_to_qp(q(1:n_vars,j,k) , B_cent(j,k) , qp(1:n_vars,j,k))
             
             WRITE(*,*) 'after imex_RK_solver: qc',q(1:n_vars,j,k)
             WRITE(*,*) 'after imex_RK_solver: qp',qp(1:n_vars,j,k)
             WRITE(*,*) 'h new',q(1,j,k) - B_cent(j,k)
             READ(*,*)

          END IF


       ENDDO

    END DO

  END SUBROUTINE imex_RK_solver

  !******************************************************************************
  !> \brief Runge-Kutta single step integration
  !
  !> This subroutine find the solution of the non-linear system 
  !> given the a step of the implicit-explicit Runge-Kutta scheme for a
  !> cell:\n
  !> \f$ Q^{(i)} = Q^n - dt \sum_{j=1}^{i-1}\tilde{a}_{j}\partial_x 
  !> F(Q^{(j)}) +  dt \sum_{j=1}^{i-1} a_j  NH(Q^{(j)}) 
  !> + dt a_{diag} NH(Q^{(i)}) \f$\n
  !
  !> \param[in]     Bj        topography at the cell center
  !> \param[in]     Bprimej   topography slope at the cell center
  !> \param[in,out] qj        conservative variables 
  !> \param[in]     qj_old    conservative variables at the old time step
  !> \param[in]     a_tilde   explicit coefficents for the fluxes
  !> \param[in]     a_dirk    explicit coefficient for the non-hyperbolic terms
  !> \param[in]     a_diag    implicit coefficient for the non-hyperbolic terms 
  !
  !> \date 07/10/2016
  !> @author 
  !> Mattia de' Michieli Vitturi
  !
  !******************************************************************************

  SUBROUTINE solve_rk_step( Bj, Bprimej_x, Bprimej_y, grav3_surf,               &
       qj, qj_old, a_tilde , a_dirk , a_diag )

    USE parameters_2d, ONLY : max_nl_iter , tol_rel , tol_abs

    USE constitutive_2d, ONLY : qc_to_qp

    IMPLICIT NONE

    REAL*8, INTENT(IN) :: Bj
    REAL*8, INTENT(IN) :: Bprimej_x
    REAL*8, INTENT(IN) :: Bprimej_y
    REAL*8, INTENT(IN) :: grav3_surf
    REAL*8, INTENT(INOUT) :: qj(n_vars)
    REAL*8, INTENT(IN) :: qj_old(n_vars)
    REAL*8, INTENT(IN) :: a_tilde(n_RK)
    REAL*8, INTENT(IN) :: a_dirk(n_RK)
    REAL*8, INTENT(IN) :: a_diag

    REAL*8 :: qj_init(n_vars)

    REAL*8 :: qj_org(n_vars) , qj_rel(n_vars)

    REAL*8 :: left_matrix(n_eqns,n_vars)
    REAL*8 :: right_term(n_eqns)

    REAL*8 :: scal_f

    REAL*8 :: coeff_f(n_eqns)

    REAL*8 :: qj_rel_NR_old(n_vars)
    REAL*8 :: scal_f_old
    REAL*8 :: desc_dir(n_vars)
    REAL*8 :: grad_f(n_vars)

    INTEGER :: pivot(n_vars)

    REAL*8 :: left_matrix_small22(n_nh,n_nh)
    REAL*8 :: left_matrix_small21(n_eqns-n_nh,n_nh)
    REAL*8 :: left_matrix_small11(n_eqns-n_nh,n_vars-n_nh)
    REAL*8 :: left_matrix_small12(n_nh,n_vars-n_nh)

    REAL*8 :: desc_dir_small2(n_nh)
    INTEGER :: pivot_small2(n_nh)

    REAL*8 :: desc_dir_small1(n_vars-n_nh)

    INTEGER :: ok

    INTEGER :: i 
    INTEGER :: nl_iter

    REAL*8, PARAMETER :: STPMX=100.D0
    REAL*8 :: stpmax
    LOGICAL :: check

    REAL*8, PARAMETER :: TOLF=1.D-10 , TOLMIN=1.D-6
    REAL*8 :: TOLX

    REAL*8 :: qpj(n_vars)

    REAL*8 :: desc_dir2(n_vars)

    REAL*8 :: desc_dir_temp(n_vars)

    normalize_q = .TRUE.
    normalize_f = .FALSE.
    opt_search_NL = .TRUE.

    coeff_f(1:n_eqns) = 1.D0

    qj_init = qj

    ! normalize the functions of the nonlinear system
    IF ( normalize_f ) THEN

       qj = qj_old - dt * ( MATMUL(divFluxj+ Expl_terms_j,a_tilde)              &
            - MATMUL(NHj,a_dirk) )

       CALL eval_f( Bj , Bprimej_x , Bprimej_y , grav3_surf ,                   &
            qj , qj_old , a_tilde , a_dirk , a_diag , coeff_f , right_term ,    &
            scal_f )

       IF ( verbose_level .GE. 3 ) THEN

          WRITE(*,*) 'solve_rk_step: non-normalized right_term'
          WRITE(*,*) right_term
          WRITE(*,*) 'scal_f',scal_f

       END IF

       DO i=1,n_eqns

          IF ( ABS(right_term(i)) .GE. 1.D0 ) coeff_f(i) = 1.D0 / right_term(i)

       END DO

       right_term = coeff_f * right_term

       scal_f = 0.5D0 * DOT_PRODUCT( right_term , right_term )

       IF ( verbose_level .GE. 3 ) THEN                    
          WRITE(*,*) 'solve_rk_step: after normalization',scal_f
       END IF

    END IF

    !---- normalize the conservative variables ------

    IF ( normalize_q ) THEN

       qj_org = qj

       qj_org = MAX( ABS(qj_org) , 1.D-3 )

    ELSE 

       qj_org(1:n_vars) = 1.D0

    END IF

    qj_rel = qj / qj_org

    ! -----------------------------------------------

    newton_raphson_loop:DO nl_iter=1,max_nl_iter

       TOLX = epsilon(qj_rel)

       IF ( verbose_level .GE. 2 ) WRITE(*,*) 'solve_rk_step: nl_iter',nl_iter

       CALL eval_f( Bj , Bprimej_x , Bprimej_y , grav3_surf ,                   &
            qj , qj_old , a_tilde , a_dirk , a_diag , coeff_f , right_term ,    &
            scal_f )

       IF ( verbose_level .GE. 2 ) THEN

          WRITE(*,*) 'solve_rk_step: right_term',right_term

       END IF

       IF ( verbose_level .GE. 2 ) THEN

          WRITE(*,*) 'before_lnsrch: scal_f',scal_f

       END IF

       ! check the residual of the system

       IF ( MAXVAL( ABS( right_term(:) ) ) < TOLF ) THEN

          IF ( verbose_level .GE. 3 ) WRITE(*,*) '1: check',check
          RETURN

       END IF

       IF ( ( normalize_f ) .AND. ( scal_f < 1.D-6 ) ) THEN

          IF ( verbose_level .GE. 3 ) WRITE(*,*) 'check scal_f',check
          RETURN

       END IF

       ! ---- evaluate the descent direction ------------------------------------

       IF ( COUNT( implicit_flag ) .EQ. n_eqns ) THEN

          CALL eval_jacobian(Bj,Bprimej_x,Bprimej_y,grav3_surf,                 &
               qj_rel,qj_org,coeff_f,left_matrix)

          desc_dir_temp = - right_term

          CALL DGESV(n_eqns,1, left_matrix , n_eqns, pivot, desc_dir_temp ,     &
               n_eqns, ok)

          desc_dir = desc_dir_temp

       ELSE

          CALL eval_jacobian(Bj,Bprimej_x,Bprimej_y,grav3_surf,                 &
               qj_rel,qj_org,coeff_f,left_matrix)

          left_matrix_small11 = reshape(pack(left_matrix, mask11),              &
               [n_eqns-n_nh,n_eqns-n_nh]) 

          left_matrix_small12 = reshape(pack(left_matrix, mask12),              &
               [n_nh,n_eqns-n_nh]) 

          left_matrix_small22 = reshape(pack(left_matrix, mask22),              &
               [n_nh,n_nh]) 

          left_matrix_small21 = reshape(pack(left_matrix, mask21),              &
               [n_eqns-n_nh,n_nh]) 


          desc_dir_small1 = pack( right_term, .NOT.implicit_flag )
          desc_dir_small2 = pack( right_term , implicit_flag )

          DO i=1,n_vars-n_nh

             desc_dir_small1(i) = desc_dir_small1(i) / left_matrix_small11(i,i)

          END DO

          desc_dir_small2 = desc_dir_small2 -                                   &
               MATMUL( desc_dir_small1 , left_matrix_small21 )

          CALL DGESV(n_nh,1, left_matrix_small22 , n_nh , pivot_small2 ,        &
               desc_dir_small2 , n_nh, ok)

          desc_dir = unpack( - desc_dir_small2 , implicit_flag , 0.0D0 )        &
               + unpack( - desc_dir_small1 , .NOT.implicit_flag , 0.0D0 )

       END IF

       IF ( verbose_level .GE. 3 ) WRITE(*,*) 'desc_dir',desc_dir

       qj_rel_NR_old = qj_rel
       scal_f_old = scal_f

       IF ( ( opt_search_NL ) .AND. ( nl_iter .GT. 1 ) ) THEN
          ! Search for the step lambda giving a suffic. decrease in the solution 

          stpmax = STPMX * MAX( SQRT( DOT_PRODUCT(qj_rel,qj_rel) ) ,            &
               DBLE( SIZE(qj_rel) ) )

          grad_f = MATMUL( right_term , left_matrix )

          desc_dir2 = desc_dir

          CALL lnsrch( Bj , Bprimej_x , Bprimej_y , grav3_surf ,                &
               qj_rel_NR_old , qj_org , qj_old , scal_f_old , grad_f ,          &
               desc_dir , coeff_f , qj_rel , scal_f , right_term , stpmax ,     &
               check )

       ELSE

          qj_rel = qj_rel_NR_old + desc_dir

          qj = qj_rel * qj_org

          CALL eval_f( Bj , Bprimej_x , Bprimej_y , grav3_surf ,                &
               qj , qj_old , a_tilde , a_dirk , a_diag , coeff_f ,              &
               right_term , scal_f )

       END IF

       IF ( verbose_level .GE. 2 ) WRITE(*,*) 'after_lnsrch: scal_f',scal_f

       qj = qj_rel * qj_org

       IF ( verbose_level .GE. 3 ) THEN

          WRITE(*,*) 'qj',qj
          CALL qc_to_qp( qj , Bj , qpj)
          WRITE(*,*) 'qp',qpj

       END IF

       
       IF ( MAXVAL( ABS( right_term(:) ) ) < TOLF ) THEN

          IF ( verbose_level .GE. 3 ) WRITE(*,*) '1: check',check
          check= .FALSE.
          RETURN

       END IF

       IF (check) THEN

          check = ( MAXVAL( ABS(grad_f(:)) * MAX( ABS( qj_rel(:) ),1.D0 ) /     &
               MAX( scal_f , 0.5D0 * SIZE(qj_rel) ) )  < TOLMIN )

          IF ( verbose_level .GE. 3 ) WRITE(*,*) '2: check',check
          !          RETURN

       END IF

       IF ( MAXVAL( ABS( qj_rel(:) - qj_rel_NR_old(:) ) / MAX( ABS( qj_rel(:)) ,&
            1.D0 ) ) < TOLX ) THEN

          IF ( verbose_level .GE. 3 ) WRITE(*,*) 'check',check
          RETURN

       END IF

    END DO newton_raphson_loop

  END SUBROUTINE solve_rk_step

  !******************************************************************************
  !> \brief Search the descent stepsize
  !
  !> This subroutine search for the lenght of the descent step in order to have
  !> a decrease in the nonlinear function.
  !> \param[in]     Bj               topography at the cell center
  !> \param[in]     Bprimej_x        topography x-slope at the cell center
  !> \param[in]     Bprimej_y        topography y-slope at the cell center
  !> \param[in]     qj_rel_NR_old  
  !> \param[in]     qj_org
  !> \param[in]     qj_old
  !> \param[in]     scal_f_old
  !> \param[in]     grad_f
  !> \param[in,out] desc_dir
  !> \param[in]     coeff_f
  !> \param[out]    qj_rel
  !> \param[out]    scal_f
  !> \param[out]    right_term
  !> \param[in]     stpmax
  !> \param[out]    check
  !> \date 07/10/2016
  !> @author 
  !> Mattia de' Michieli Vitturi
  !******************************************************************************

  SUBROUTINE lnsrch( Bj , Bprimej_x , Bprimej_y , grav3_surf ,                  &
       qj_rel_NR_old , qj_org , qj_old , scal_f_old , grad_f ,                  &
       desc_dir , coeff_f , qj_rel , scal_f , right_term , stpmax , check )

    IMPLICIT NONE

    REAL*8, INTENT(IN) :: Bj

    REAL*8, INTENT(IN) :: Bprimej_x

    REAL*8, INTENT(IN) :: Bprimej_y

    REAL*8, INTENT(IN) :: grav3_surf

    !> Initial point
    REAL*8, DIMENSION(:), INTENT(IN) :: qj_rel_NR_old

    !> Initial point
    REAL*8, DIMENSION(:), INTENT(IN) :: qj_org

    !> Initial point
    REAL*8, DIMENSION(:), INTENT(IN) :: qj_old

    !> Gradient at xold
    REAL*8, DIMENSION(:), INTENT(IN) :: grad_f

    !> Value of the function at xold
    REAL*8, INTENT(IN) :: scal_f_old

    !> Descent direction (usually Newton direction)
    REAL*8, DIMENSION(:), INTENT(INOUT) :: desc_dir

    REAL*8, INTENT(IN) :: stpmax

    !> Coefficients to rescale the nonlinear function
    REAL*8, DIMENSION(:), INTENT(IN) :: coeff_f

    !> Updated solution
    REAL*8, DIMENSION(:), INTENT(OUT) :: qj_rel

    !> Value of the scalar function at x
    REAL*8, INTENT(OUT) :: scal_f

    !> Value of the scalar function at x
    REAL*8, INTENT(OUT) :: right_term(n_eqns)

    !> Output quantity check is false on a normal exit 
    LOGICAL, INTENT(OUT) :: check

    REAL*8, PARAMETER :: TOLX=epsilon(qj_rel)

    INTEGER, DIMENSION(1) :: ndum
    REAL*8 :: ALF , a,alam,alam2,alamin,b,disc
    REAL*8 :: scal_f2
    REAL*8 :: desc_dir_abs
    REAL*8 :: rhs1 , rhs2 , slope, tmplam

    REAL*8 :: scal_f_min , alam_min

    REAL*8 :: qj(n_vars)

    ALF = 1.0d-4

    IF ( size(grad_f) == size(desc_dir) .AND. size(grad_f) == size(qj_rel)      &
         .AND. size(qj_rel) == size(qj_rel_NR_old) ) THEN

       ndum = size(grad_f)

    ELSE

       WRITE(*,*) 'nrerror: an assert_eq failed with this tag:', 'lnsrch'
       STOP 'program terminated by assert_eq4'

    END IF

    check = .FALSE.

    desc_dir_abs = SQRT( DOT_PRODUCT(desc_dir,desc_dir) )

    IF ( desc_dir_abs > stpmax ) desc_dir(:) = desc_dir(:) * stpmax/desc_dir_abs  

    slope = DOT_PRODUCT(grad_f,desc_dir)

    alamin = TOLX / MAXVAL( ABS( desc_dir(:))/MAX( ABS(qj_rel_NR_old(:)),1.D0 ) )

    IF ( alamin .EQ. 0.d0) THEN

       qj_rel(:) = qj_rel_NR_old(:)

       RETURN

    END IF

    alam = 1.0D0

    scal_f_min = scal_f_old

    optimal_step_search: DO

       IF ( verbose_level .GE. 4 ) THEN

          WRITE(*,*) 'alam',alam

       END IF

       qj_rel = qj_rel_NR_old + alam * desc_dir

       qj = qj_rel * qj_org

       CALL eval_f( Bj , Bprimej_x , Bprimej_y , grav3_surf , qj , qj_old ,     &
            a_tilde , a_dirk , a_diag , coeff_f , right_term , scal_f )

       IF ( verbose_level .GE. 4 ) THEN

          WRITE(*,*) 'lnsrch: effe_old,effe',scal_f_old,scal_f
          READ(*,*)

       END IF

       IF ( scal_f .LT. scal_f_min ) THEN

          scal_f_min = scal_f
          alam_min = alam

       END IF

       IF ( scal_f .LE. 0.9 * scal_f_old ) THEN   
          ! sufficient function decrease

          IF ( verbose_level .GE. 4 ) THEN

             WRITE(*,*) 'sufficient function decrease'

          END IF

          EXIT optimal_step_search   

       ELSE IF ( alam < alamin ) THEN   
          ! convergence on Delta_x

          IF ( verbose_level .GE. 4 ) THEN

             WRITE(*,*) ' convergence on Delta_x',alam,alamin

          END IF

          qj_rel(:) = qj_rel_NR_old(:)
          scal_f = scal_f_old
          check = .TRUE.

          EXIT optimal_step_search

          !       ELSE IF ( scal_f .LE. scal_f_old + ALF * alam * slope ) THEN   
       ELSE  

          IF ( alam .EQ. 1.D0 ) THEN

             tmplam = - slope / ( 2.0D0 * ( scal_f - scal_f_old - slope ) )

          ELSE

             rhs1 = scal_f - scal_f_old - alam*slope
             rhs2 = scal_f2 - scal_f_old - alam2*slope

             a = ( rhs1/alam**2.D0 - rhs2/alam2**2.D0 ) / ( alam - alam2 )
             b = ( -alam2*rhs1/alam**2 + alam*rhs2/alam2**2 ) / ( alam - alam2 )

             IF ( a .EQ. 0.D0 ) THEN

                tmplam = - slope / ( 2.0D0 * b )

             ELSE

                disc = b*b - 3.0D0*a*slope

                IF ( disc .LT. 0.D0 ) THEN

                   tmplam = 0.5D0 * alam

                ELSE IF ( b .LE. 0.D0 ) THEN

                   tmplam = ( - b + SQRT(disc) ) / ( 3.D0 * a )

                ELSE

                   tmplam = - slope / ( b + SQRT(disc) )

                ENDIF

             END IF

             IF ( tmplam .GT. 0.5D0*alam ) tmplam = 0.5D0 * alam

          END IF

       END IF

       alam2 = alam
       scal_f2 = scal_f
       alam = MAX( tmplam , 0.5D0*alam )

    END DO optimal_step_search

  END SUBROUTINE lnsrch

  !******************************************************************************
  !> \brief Evaluate the nonlinear system
  !
  !> This subroutine evaluate the value of the nonlinear system in the state 
  !> defined by the variables qj.
  !> \param[in]    Bj          topography at the cell center
  !> \param[in]    Bprimej     topography slope at the cell center
  !> \param[in]    qj          conservative variables 
  !> \param[in]    qj_old      conservative variables at the old time step
  !> \param[in]    a_tilde     explicit coefficients for the hyperbolic terms 
  !> \param[in]    a_dirk      explicit coefficients for the non-hyperbolic terms 
  !> \param[in]    a_diag      implicit coefficient for the non-hyperbolic term
  !> \param[in]    coeff_f     coefficient to rescale the nonlinear functions
  !> \param[out]   f_nl        values of the nonlinear functions
  !> \param[out]   scal_f      value of the scalar function f=0.5*<F,F>
  !> \date 07/10/2016
  !> @author 
  !> Mattia de' Michieli Vitturi
  !******************************************************************************

  SUBROUTINE eval_f( Bj , Bprimej_x , Bprimej_y , grav3_surf , qj , qj_old ,    &
       a_tilde , a_dirk , a_diag , coeff_f , f_nl , scal_f )

    USE constitutive_2d, ONLY : eval_nonhyperbolic_terms

    IMPLICIT NONE

    REAL*8, INTENT(IN) :: Bj
    REAL*8, INTENT(IN) :: Bprimej_x
    REAL*8, INTENT(IN) :: Bprimej_y
    REAL*8, INTENT(IN) :: grav3_surf
    REAL*8, INTENT(IN) :: qj(n_vars)
    REAL*8, INTENT(IN) :: qj_old(n_vars)
    REAL*8, INTENT(IN) :: a_tilde(n_RK)
    REAL*8, INTENT(IN) :: a_dirk(n_RK)
    REAL*8, INTENT(IN) :: a_diag
    REAL*8, INTENT(IN) :: coeff_f(n_eqns)
    REAL*8, INTENT(OUT) :: f_nl(n_eqns)
    REAL*8, INTENT(OUT) :: scal_f

    REAL*8 :: nh_term_impl(n_eqns)
    REAL*8 :: Rj(n_eqns)

    CALL eval_nonhyperbolic_terms( Bj , Bprimej_x , Bprimej_y , grav3_surf ,    &
         r_qj = qj , r_nh_term_impl = nh_term_impl ) 

    Rj = ( MATMUL(divFluxj+ Expl_terms_j,a_tilde) - MATMUL(NHj,a_dirk) )        &
         - a_diag * nh_term_impl

    f_nl = qj - qj_old + dt * Rj

    f_nl = coeff_f * f_nl

    scal_f = 0.5D0 * DOT_PRODUCT( f_nl , f_nl )

  END SUBROUTINE eval_f

  !******************************************************************************
  !> \brief Evaluate the jacobian 
  !
  !> This subroutine evaluate the jacobian of the non-linear system
  !> with respect to the conservative variables.
  !
  !> \param[in]    Bj            topography at the cell center
  !> \param[in]    Bprimej_x     topography x-slope at the cell center
  !> \param[in]    Bprimej_y     topography y-slope at the cell center
  !> \param[in]    grav3_surf
  !> \param[in]    qj_rel        relative variation (qj=qj_rel*qj_org)
  !> \param[in]    qj_org        conservative variables at the old time step
  !> \param[in]    coeff_f       coefficient to rescale the nonlinear functions
  !> \param[out]   left_matrix   matrix from the linearization of the system
  !
  !> \date 07/10/2016
  !> @author 
  !> Mattia de' Michieli Vitturi
  !******************************************************************************

  SUBROUTINE eval_jacobian(Bj , Bprimej_x , Bprimej_y , grav3_surf ,            &
       qj_rel , qj_org , coeff_f, left_matrix)

    USE constitutive_2d, ONLY : eval_nonhyperbolic_terms

    IMPLICIT NONE

    REAL*8, INTENT(IN) :: Bj
    REAL*8, INTENT(IN) :: Bprimej_x
    REAL*8, INTENT(IN) :: Bprimej_y
    REAL*8, INTENT(IN) :: grav3_surf
    REAL*8, INTENT(IN) :: qj_rel(n_vars)
    REAL*8, INTENT(IN) :: qj_org(n_vars)
    REAL*8, INTENT(IN) :: coeff_f(n_eqns)
    REAL*8, INTENT(OUT) :: left_matrix(n_eqns,n_vars)

    REAL*8 :: Jacob_relax(n_eqns,n_vars)
    COMPLEX*16 :: nh_terms_cmplx_impl(n_eqns)
    COMPLEX*16 :: qj_cmplx(n_vars) , qj_rel_cmplx(n_vars)

    REAL*8 :: h

    INTEGER :: i

    h = n_vars * epsilon(1.d0)

    ! initialize the matrix of the linearized system and the Jacobian

    left_matrix(1:n_eqns,1:n_vars) = 0.D0
    Jacob_relax(1:n_eqns,1:n_vars) = 0.D0

    ! evaluate the jacobian of the non-hyperbolic terms

    DO i=1,n_vars

       left_matrix(i,i) = coeff_f(i) * qj_org(i)

       IF ( implicit_flag(i) ) THEN 

          qj_rel_cmplx(1:n_vars) = qj_rel(1:n_vars)
          qj_rel_cmplx(i) = DCMPLX(qj_rel(i), h)

          qj_cmplx = qj_rel_cmplx * qj_org

          CALL eval_nonhyperbolic_terms( Bj , Bprimej_x , Bprimej_y ,           &
               grav3_surf, c_qj = qj_cmplx ,                                    &
               c_nh_term_impl = nh_terms_cmplx_impl ) 

          Jacob_relax(1:n_eqns,i) = coeff_f(i) *                                &
               AIMAG(nh_terms_cmplx_impl) / h

          left_matrix(1:n_eqns,i) = left_matrix(1:n_eqns,i) - dt * a_diag       &
               * Jacob_relax(1:n_eqns,i)

       END IF

    END DO

  END SUBROUTINE eval_jacobian

  !******************************************************************************
  !> \brief Evaluate the explicit terms 
  !
  !> This subroutine evaluate the explicit terms (non-fluxes) of the non-linear 
  !> system with respect to the conservative variables.
  !
  !> \param[in]    q_expl          conservative variables 
  !> \param[out]   expl_terms      explicit terms
  !
  !> \date 07/10/2016
  !> @author 
  !> Mattia de' Michieli Vitturi
  !******************************************************************************

  SUBROUTINE eval_explicit_terms( q_expl , expl_terms )

    USE constitutive_2d, ONLY : eval_expl_terms

    IMPLICIT NONE

    REAL*8, INTENT(IN) :: q_expl(n_vars,comp_cells_x,comp_cells_y)
    REAL*8, INTENT(OUT) :: expl_terms(n_eqns,comp_cells_x,comp_cells_y)

    REAL*8 :: qc(n_vars)      !< conservative variables 
    REAL*8 :: expl_forces_term(n_eqns)      !< conservative variables 

    INTEGER :: j,k

    DO j = 1,comp_cells_x

       DO k = 1,comp_cells_y

          qc = q_expl(1:n_vars,j,k)

          CALL eval_expl_terms( B_cent(j,k), B_prime_x(j,k),                    &
               B_prime_y(j,k), source_xy(j,k) , qc, expl_forces_term )

          expl_terms(1:n_eqns,j,k) =  expl_forces_term

       ENDDO

    END DO

  END SUBROUTINE eval_explicit_terms

  !******************************************************************************
  !> \brief Semidiscrete finite volume central scheme
  !
  !> This subroutine compute the divergence part of the system of the eqns,
  !> with a modified version of the finite volume scheme from Kurganov et al.  
  !> 2001, where the reconstruction at the cells interfaces is applied to a
  !> set of physical variables derived from the conservative vriables.
  !
  !> \param[in]     q_expl        conservative variables
  !> \param[out]    divFlux           divergence term
  !
  !> \date 07/10/2016
  !> @author 
  !> Mattia de' Michieli Vitturi
  !******************************************************************************

  SUBROUTINE eval_hyperbolic_terms( q_expl , divFlux )

    ! External variables
    USE geometry_2d, ONLY : dx,dy
    USE parameters_2d, ONLY : solver_scheme

    IMPLICIT NONE

    REAL*8, INTENT(IN) :: q_expl(n_vars,comp_cells_x,comp_cells_y)
    REAL*8, INTENT(OUT) :: divFlux(n_eqns,comp_cells_x,comp_cells_y)

    REAL*8 :: q_old(n_vars,comp_cells_x,comp_cells_y)

    REAL*8 :: h_new

    INTEGER :: i, j, k      !< loop counters

    q_old = q

    q = q_expl

    ! Linear reconstruction of the physical variables at the interfaces
    CALL reconstruction

    ! Evaluation of the maximum local speeds at the interfaces
    CALL eval_speeds

    ! Evaluation of the numerical fluxes
    SELECT CASE ( solver_scheme )

    CASE ("LxF")

       CALL eval_flux_LxF

    CASE ("GFORCE")

       CALL eval_flux_GFORCE

    CASE ("KT")

       CALL eval_flux_KT

    END SELECT

    ! Advance in time the solution
    DO j = 1,comp_cells_x

       DO k = 1,comp_cells_y
          
          DO i=1,n_eqns
             
             divFlux(i,j,k) = 0.D0
             
             IF ( comp_cells_x .GT. 1 ) THEN
                
                divFlux(i,j,k) = divFlux(i,j,k) +                               &
                     ( H_interface_x(i,j+1,k) - H_interface_x(i,j,k) ) / dx
                
             END IF

             IF ( comp_cells_y .GT. 1 ) THEN
                
                divFlux(i,j,k) = divFlux(i,j,k) +                               &
                     ( H_interface_y(i,j,k+1) - H_interface_y(i,j,k) ) / dy
                
             END IF

             
          END DO

          h_new = q_expl(1,j,k) - dt * divFlux(1,j,k) - B_cent(j,k)

          ! IF ( h_new .NE. 0.D0 ) THEN
          IF ( j .EQ. 0 ) THEN

             WRITE(*,*) 'j,k,h,divF',j,k,h_new, divFlux(1,j,k)
             WRITE(*,*) 'dt',dt

             WRITE(*,*) 'h_interface(j,k) ',q_interfaceL(1,j,k)-B_stag_x(j,k) , &
                  q_interfaceR(1,j,k) - B_stag_x(j,k)

             WRITE(*,*) 'hu_interface(j,k)' , q_interfaceL(2,j,k) ,             &
                  q_interfaceR(2,j,k)

             WRITE(*,*) 'h_interface(j+1,k) ' , q_interfaceL(1,j+1,k)           &
                  - B_stag_y(j+1,k) , q_interfaceR(1,j+1,k) - B_stag_y(j+1,k)

             WRITE(*,*) 'hu_interface(j+1,k)' , q_interfaceL(2,j+1,k) ,         &
                  q_interfaceR(2,j+1,k)

             WRITE(*,*) 'h flux x',H_interface_x(1,j,k) , H_interface_x(1,j+1,k)

             WRITE(*,*) 'h flux y',H_interface_y(i,j,k) , H_interface_y(1,j,k+1)

             WRITE(*,*) 'hu flux x',H_interface_x(2,j,k) , H_interface_x(2,j+1,k)

             WRITE(*,*) 'hu flux y',H_interface_y(2,j,k) , H_interface_y(2,j,k+1)

             READ(*,*) 

          END IF

       ENDDO

    END DO

    q = q_old

  END SUBROUTINE eval_hyperbolic_terms


  !******************************************************************************
  !> \brief Semidiscrete numerical fluxes
  !
  !> This subroutine evaluates the numerical fluxes H at the 
  !> cells interfaces according to Kurganov et al. 2001. 
  !> @author 
  !> Mattia de' Michieli Vitturi
  !> \date 16/08/2011
  !******************************************************************************

  SUBROUTINE eval_flux_KT

    ! External procedures
    USE constitutive_2d, ONLY : eval_fluxes

    IMPLICIT NONE

    REAL*8 :: fluxL(n_eqns)           !< Numerical fluxes from the eqns 
    REAL*8 :: fluxR(n_eqns)           !< Numerical fluxes from the eqns
    REAL*8 :: fluxB(n_eqns)           !< Numerical fluxes from the eqns 
    REAL*8 :: fluxT(n_eqns)           !< Numerical fluxes from the eqns

    REAL*8 :: flux_avg_x(n_eqns)   
    REAL*8 :: flux_avg_y(n_eqns)   

    INTEGER :: j,k                      !< Loop counter
    INTEGER :: i                      !< Loop counter

    IF ( comp_cells_x .GT. 1 ) THEN

       DO j = 1,comp_interfaces_x

          DO k = 1,comp_cells_y

             CALL eval_fluxes( B_stag_x(j,k) ,                                  &
                  r_qj = q_interfaceL(1:n_vars,j,k) , r_flux=fluxL , dir=1 )

             CALL eval_fluxes( B_stag_x(j,k) ,                                  &
                  r_qj = q_interfaceR(1:n_vars,j,k) , r_flux=fluxR , dir=1 )

             CALL average_KT( a_interface_xNeg(:,j,k), a_interface_xPos(:,j,k) ,&
                  fluxL , fluxR , flux_avg_x )

             eqns_loop:DO i=1,n_eqns

                IF ( a_interface_xNeg(i,j,k) .EQ. a_interface_xPos(i,j,k) ) THEN

                   H_interface_x(i,j,k) = 0.D0

                ELSE

                   H_interface_x(i,j,k) = flux_avg_x(i)                         &
                        + ( a_interface_xPos(i,j,k) * a_interface_xNeg(i,j,k) ) &
                        / ( a_interface_xPos(i,j,k) - a_interface_xNeg(i,j,k) ) &
                        * ( q_interfaceR(i,j,k) - q_interfaceL(i,j,k) )             

                END IF

             ENDDO eqns_loop

          END DO

       END DO

    END IF

    IF ( comp_cells_y .GT. 1 ) THEN

       DO j = 1,comp_cells_x

          DO k = 1,comp_interfaces_y

             CALL eval_fluxes( B_stag_y(j,k) ,                                  &
                  r_qj = q_interfaceB(1:n_vars,j,k) , r_flux=fluxB , dir=2 )

             CALL eval_fluxes( B_stag_y(j,k) ,                                  &
                  r_qj = q_interfaceT(1:n_vars,j,k) , r_flux=fluxT , dir=2 )


             CALL average_KT( a_interface_yNeg(:,j,k) ,                         &
                  a_interface_yPos(:,j,k) , fluxB , fluxT , flux_avg_y )

             DO i=1,n_eqns

                IF ( a_interface_yNeg(i,j,k) .EQ. a_interface_yPos(i,j,k) ) THEN

                   H_interface_y(i,j,k) = 0.D0

                ELSE

                   H_interface_y(i,j,k) = flux_avg_y(i)                         &
                        + ( a_interface_yPos(i,j,k) * a_interface_yNeg(i,j,k) ) &
                        / ( a_interface_yPos(i,j,k) - a_interface_yNeg(i,j,k) ) &
                        * ( q_interfaceT(i,j,k) - q_interfaceB(i,j,k) )             

                END IF

             END DO

             IF ( ( j .EQ. 0 ) .AND. ( k .EQ. 3 ) ) THEN

                WRITE(*,*) 'eval_flux_KT',j,k
                WRITE(*,*) 'q_interfaceB',q_interfaceB(1:n_vars,j,k),B_stag_y(j,k)
                WRITE(*,*) 'fluxB',fluxB
                WRITE(*,*) 'q_interfaceT',q_interfaceT(1:n_vars,j,k),B_stag_y(j,k)
                WRITE(*,*) 'fluxT',fluxT
                WRITE(*,*) 'H_interface_y',H_interface_y(:,j,k)
                WRITE(*,*) 'flux_avg_y',flux_avg_y(:) 

             END IF

             
          ENDDO

       END DO

    END IF

  END SUBROUTINE eval_flux_KT

  !******************************************************************************
  !> \brief averaged KT flux
  !
  !> This subroutine compute n averaged flux from the fluxes at the two sides of
  !> a cell interface and the max an min speed at the two sides.
  !> \param[in]     aL            speed at one side of the interface
  !> \param[in]     aR            speed at the other side of the interface
  !> \param[in]     wL            fluxes at one side of the interface
  !> \param[in]     wR            fluxes at the other side of the interface
  !> \param[out]    w_avg         array of averaged fluxes
  !> \date 07/10/2016
  !> @author 
  !> Mattia de' Michieli Vitturi
  !******************************************************************************


  SUBROUTINE average_KT( a1 , a2 , w1 , w2 , w_avg )

    IMPLICIT NONE

    REAL*8, INTENT(IN) :: a1(:) , a2(:)
    REAL*8, INTENT(IN) :: w1(:) , w2(:)
    REAL*8, INTENT(OUT) :: w_avg(:)

    INTEGER :: n
    INTEGER :: i 

    n = SIZE( a1 )

    DO i=1,n

       IF ( a1(i) .EQ. a2(i) ) THEN

          w_avg(i) = 0.5D0 * ( w1(i) + w2(i) )
          w_avg(i) = 0.D0

       ELSE

          w_avg(i) = ( a2(i) * w1(i) - a1(i) * w2(i) ) / ( a2(i) - a1(i) )  

       END IF

    END DO

  END SUBROUTINE average_KT

  !******************************************************************************
  !> \brief Numerical fluxes GFORCE
  !> \date 07/10/2016
  !> @author 
  !> Mattia de' Michieli Vitturi
  !******************************************************************************

  SUBROUTINE eval_flux_GFORCE

    ! to be implemented
    WRITE(*,*) 'method not yet implemented in 2-d case'

  END SUBROUTINE eval_flux_GFORCE

  !******************************************************************************
  !> \brief Numerical fluxes Lax-Friedrichs
  !> \date 07/10/2016
  !> @author 
  !> Mattia de' Michieli Vitturi
  !******************************************************************************

  SUBROUTINE eval_flux_LxF

    ! to be implemented
    WRITE(*,*) 'method not yet implemented in 2-d case'

  END SUBROUTINE eval_flux_LxF


  !******************************************************************************
  !> \brief Linear reconstruction
  !
  !> In this subroutine a linear reconstruction with slope limiters is
  !> applied to a set of variables describing the state of the system, according
  !> to the input parameter reconstr_variables.
  !> @author 
  !> Mattia de' Michieli Vitturi
  !> \date 15/08/2011
  !******************************************************************************

  SUBROUTINE reconstruction

    ! External procedures
    USE constitutive_2d, ONLY : qc_to_qp , qp_to_qc
    USE constitutive_2d, ONLY : qrec_to_qc
    USE parameters_2d, ONLY : limiter

    ! External variables
    USE geometry_2d, ONLY : x_comp , x_stag , y_comp , y_stag , dx , dx2 , dy , &
         dy2

    USE parameters_2d, ONLY : reconstr_variables , reconstr_coeff

    IMPLICIT NONE

    REAL*8 :: qc(n_vars)      !< conservative variables
    REAL*8 :: qrec(n_vars,comp_cells_x,comp_cells_y) 
    REAL*8 :: qrecW(n_vars)     !< recons var at the west edge of the cells
    REAL*8 :: qrecE(n_vars)     !< recons var at the east edge of the cells
    REAL*8 :: qrecS(n_vars)     !< recons var at the south edge of the cells
    REAL*8 :: qrecN(n_vars)     !< recons var at the north edge of the cells
    REAL*8 :: qrec_bdry(n_vars) !< recons variables outside the domain

    REAL*8 :: qrec_stencil(3)   !< recons variables stencil for the limiter
    REAL*8 :: x_stencil(3)    !< grid stencil for the limiter
    REAL*8 :: y_stencil(3)    !< grid stencil for the limiter
    REAL*8 :: qrec_prime_x      !< recons variables slope
    REAL*8 :: qrec_prime_y      !< recons variables slope

    INTEGER :: j,k            !< loop counter (cells)
    INTEGER :: i              !< loop counter (variables)


    ! Compute the variable to reconstruct (phys or cons)
    DO j = 1,comp_cells_x

       DO k = 1,comp_cells_y

          qc = q(1:n_vars,j,k)

          IF ( reconstr_variables .EQ. 'phys' ) THEN
          
             CALL qc_to_qp( qc , B_cent(j,k) , qrec(1:n_vars,j,k) )

          ELSEIF ( reconstr_variables .EQ. 'cons' ) THEN

             qrec(1:n_vars,j,k) = qc
             
          END IF
             
       END DO

    ENDDO

    ! Linear reconstruction

    x_loop:DO j = 1,comp_cells_x

       y_loop:DO k = 1,comp_cells_y

          vars_loop:DO i=1,n_vars

             ! x direction
             IF ( comp_cells_x .GT. 1 ) THEN

                ! west boundary
                IF (j.EQ.1) THEN
                   
                   IF ( bcW(i)%flag .EQ. 0 ) THEN
                      
                      x_stencil(1) = x_stag(1)
                      x_stencil(2:3) = x_comp(1:2)
                      
                      qrec_stencil(1) = bcW(i)%value
                      qrec_stencil(2:3) = qrec(i,1:2,k)
                      
                      CALL limit( qrec_stencil , x_stencil , limiter(i) ,       &
                           qrec_prime_x ) 
                      
                   ELSEIF ( bcW(i)%flag .EQ. 1 ) THEN
                      
                      qrec_prime_x = bcW(i)%value
                      
                   ELSEIF ( bcW(i)%flag .EQ. 2 ) THEN
                      
                      qrec_prime_x = ( qrec(i,2,k) - qrec(i,1,k) ) / dx
                      
                   END IF

                   !east boundary
                ELSEIF (j.EQ.comp_cells_x) THEN
                   
                   IF ( bcE(i)%flag .EQ. 0 ) THEN
                      
                      qrec_stencil(3) = bcE(i)%value
                      qrec_stencil(1:2) = qrec(i,comp_cells_x-1:comp_cells_x,k)
                      
                      x_stencil(3) = x_stag(comp_interfaces_x)
                      x_stencil(1:2) = x_comp(comp_cells_x-1:comp_cells_x)
                      
                      CALL limit( qrec_stencil , x_stencil , limiter(i) ,       &
                           qrec_prime_x ) 
                      
                   ELSEIF ( bcE(i)%flag .EQ. 1 ) THEN
                      
                      qrec_prime_x = bcE(i)%value
                      
                   ELSEIF ( bcE(i)%flag .EQ. 2 ) THEN
                      
                      qrec_prime_x = ( qrec(i,comp_cells_x,k) -                 &
                           qrec(i,comp_cells_x-1,k) ) / dx
                      
                   END IF
                   
                   ! internal x cells
                ELSE
                   
                   x_stencil(1:3) = x_comp(j-1:j+1)
                   qrec_stencil = qrec(i,j-1:j+1,k)
                   
                   CALL limit( qrec_stencil , x_stencil , limiter(i) ,          &
                        qrec_prime_x )
                   
                ENDIF

                qrecW(i) = qrec(i,j,k) - reconstr_coeff * dx2 * qrec_prime_x
                qrecE(i) = qrec(i,j,k) + reconstr_coeff * dx2 * qrec_prime_x

                ! positivity preserving reconstruction for h (i=1, 1st variable)
                IF ( i .EQ. 1 ) THEN
                               
                   IF ( qrecE(i) .LT. B_stag_x(j+1,k) ) THEN
                      
                      qrec_prime_x = ( B_stag_x(j+1,k) - qrec(i,j,k) ) / dx2
                      
                      qrecE(i) = qrec(i,j,k) + dx2 * qrec_prime_x

                      qrecW(i) = 2.D0 * qrec(i,j,k) - qrecE(i) 
                                            
                   ENDIF
                   
                   IF ( qrecW(i) .LT. B_stag_x(j,k) ) THEN
                      
                      qrec_prime_x = ( qrec(i,j,k) - B_stag_x(j,k) ) / dx2
                      
                      qrecW(i) = qrec(i,j,k) - dx2 * qrec_prime_x
                      
                      qrecE(i) = 2.D0 * qrec(i,j,k) - qrecW(i) 
                      
                   ENDIF
                
                END IF
          
             END IF
                
             ! y-direction
             check_comp_cells_y:IF ( comp_cells_y .GT. 1 ) THEN

                ! South boundary
                check_y_boundary:IF (k.EQ.1) THEN
                   
                   IF ( bcS(i)%flag .EQ. 0 ) THEN
                      
                      qrec_stencil(1) = bcS(i)%value
                      qrec_stencil(2:3) = qrec(i,j,1:2)
                      
                      y_stencil(1) = y_stag(1)
                      y_stencil(2:3) = y_comp(1:2)
                      
                      CALL limit( qrec_stencil , y_stencil , limiter(i) ,       &
                           qrec_prime_y ) 
                      
                   ELSEIF ( bcS(i)%flag .EQ. 1 ) THEN
                      
                      qrec_prime_y = bcS(i)%value
                      
                   ELSEIF ( bcS(i)%flag .EQ. 2 ) THEN
                      
                      qrec_prime_y = ( qrec(i,j,2) - qrec(i,j,1) ) / dy 
                      
                   END IF
                   
                   ! North boundary
                ELSEIF ( k .EQ. comp_cells_y ) THEN
                   
                   IF ( bcN(i)%flag .EQ. 0 ) THEN
                      
                      qrec_stencil(3) = bcN(i)%value
                      qrec_stencil(1:2) = qrec(i,j,comp_cells_y-1:comp_cells_y)
                      
                      y_stencil(3) = y_stag(comp_interfaces_y)
                      y_stencil(1:2) = y_comp(comp_cells_y-1:comp_cells_y)
                      
                      CALL limit( qrec_stencil , y_stencil , limiter(i) ,       &
                           qrec_prime_y ) 
                      
                   ELSEIF ( bcN(i)%flag .EQ. 1 ) THEN
                      
                      qrec_prime_y = bcN(i)%value
                      
                   ELSEIF ( bcN(i)%flag .EQ. 2 ) THEN
                      
                      qrec_prime_y = ( qrec(i,j,comp_cells_y) -                 &
                           qrec(i,j,comp_cells_y-1) ) / dy 
                      
                   END IF
                   
                   ! Internal y cells
                ELSE
                   
                   y_stencil(1:3) = y_comp(k-1:k+1)
                   qrec_stencil = qrec(i,j,k-1:k+1)
                   
                   CALL limit( qrec_stencil , y_stencil , limiter(i) ,          &
                        qrec_prime_y )
                   
                ENDIF check_y_boundary
                
                qrecS(i) = qrec(i,j,k) - reconstr_coeff * dy2 * qrec_prime_y
                qrecN(i) = qrec(i,j,k) + reconstr_coeff * dy2 * qrec_prime_y

                ! positivity preserving reconstruction for h
                IF ( i .EQ. 1 ) THEN
                   
                   IF ( qrecN(i) .LT. B_stag_y(j,k+1) ) THEN
                      
                      qrec_prime_y = ( B_stag_y(j,k+1) - qrec(i,j,k) ) / dy2
                      
                      qrecN(i) = qrec(i,j,k) + dy2 * qrec_prime_y
                      
                      qrecS(i) = 2.D0 * qrec(i,j,k) - qrecN(i) 
                      
                   ENDIF
                   
                   IF ( qrecS(i) .LT. B_stag_y(j,k) ) THEN
                      
                      qrec_prime_y = ( qrec(i,j,k) - B_stag_y(j,k) ) / dy2
                      
                      qrecS(i) = qrec(i,j,k) - dy2 * qrec_prime_y
                      
                      qrecN(i) = 2.D0 * qrec(i,j,k) - qrecS(i)
                      
                   ENDIF

                END IF
                
             ENDIF check_comp_cells_y
             
          ENDDO vars_loop


          IF ( reconstr_variables .EQ. 'phys' ) THEN
          
             ! Convert back from physical to conservative variables
             IF ( comp_cells_x .GT. 1 ) THEN

                CALL qp_to_qc( qrecW , B_stag_x(j,k) , q_interfaceR(:,j,k) )
                CALL qp_to_qc( qrecE , B_stag_x(j+1,k) , q_interfaceL(:,j+1,k) )

             END IF

             IF ( comp_cells_y .GT. 1 ) THEN
                
                CALL qp_to_qc( qrecS , B_stag_y(j,k) , q_interfaceT(:,j,k) )
                CALL qp_to_qc( qrecN , B_stag_y(j,k+1) , q_interfaceB(:,j,k+1) )

             END IF
                
          ELSE IF ( reconstr_variables .EQ. 'cons' ) THEN
             
             IF ( comp_cells_x .GT. 1 ) THEN

                CALL qrec_to_qc( qrecW, B_stag_x(j,k) , q_interfaceR(:,j,k) )
                CALL qrec_to_qc( qrecE, B_stag_x(j+1,k) , q_interfaceL(:,j+1,k) )

             END IF

             IF ( comp_cells_y .GT. 1 ) THEN
                
                CALL qrec_to_qc( qrecS, B_stag_y(j,k) , q_interfaceT(:,j,k) )
                CALL qrec_to_qc( qrecN, B_stag_y(j,k+1) , q_interfaceB(:,j,k+1) )

             END IF
                
          END IF
          
          IF ( j.EQ.0 ) THEN

             WRITE(*,*) 'qrecW',qrecW
             WRITE(*,*) 'qrecE',qrecE

          END IF

          ! boundary conditions

          IF ( comp_cells_y .GT. 1 ) THEN
          
             ! South q_interfaceB(:,j,1)
             IF ( k .EQ. 1 ) THEN
                
                DO i=1,n_vars
                   
                   IF ( bcS(i)%flag .EQ. 0 ) THEN
                      
                      qrec_bdry(i) = bcS(i)%value 
                      
                   ELSE
                      
                      qrec_bdry(i) = qrecS(i)
                      
                   END IF
                   
                ENDDO
                
                IF ( reconstr_variables .EQ. 'phys' ) THEN
                   
                   CALL qp_to_qc( qrec_bdry, B_stag_y(j,1), q_interfaceB(:,j,1) )
                   
                ELSEIF ( reconstr_variables .EQ. 'cons' ) THEN
                   
                   CALL qrec_to_qc( qrec_bdry , B_stag_y(j,1) ,                 &
                        q_interfaceB(:,j,1) )
                   
                END IF
                
             ENDIF
             
             ! North qT(i,j,comp_interfaces_y)
             IF(k.EQ.comp_cells_y)THEN
                
                DO i=1,n_vars
                   
                   IF ( bcN(i)%flag .EQ. 0 ) THEN
                      
                      qrec_bdry(i) = bcN(i)%value 
                      
                   ELSE
                      
                      qrec_bdry(i) = qrecN(i)
                      
                   END IF
                   
                ENDDO
                
                IF ( reconstr_variables .EQ. 'phys' ) THEN
                   
                   CALL qp_to_qc( qrec_bdry , B_stag_y(j,comp_interfaces_y) ,   &
                        q_interfaceT(:,j,comp_interfaces_y) )
                   
                ELSEIF ( reconstr_variables .EQ. 'cons' ) THEN
                   
                   CALL qrec_to_qc( qrec_bdry , B_stag_y(j,comp_interfaces_y) , &
                        q_interfaceT(:,j,comp_interfaces_y) )
                   
                END IF
                
             ENDIF

          END IF

          IF ( comp_cells_x .GT. 1 ) THEN
          
             ! West q_interfaceL(:,1,k)
             IF ( j.EQ.1 ) THEN
                
                DO i=1,n_vars
                   
                   IF ( bcW(i)%flag .EQ. 0 ) THEN
                      
                      qrec_bdry(i) = bcW(i)%value 
                      
                   ELSE
                      
                      qrec_bdry(i) = qrecW(i)
                      
                   END IF
                   
                ENDDO
                
                IF ( reconstr_variables .EQ. 'phys' ) THEN
                   
                   CALL qp_to_qc( qrec_bdry, B_stag_x(1,k), q_interfaceL(:,1,k) )
                   
                ELSEIF ( reconstr_variables .EQ. 'cons' ) THEN
                   
                   CALL qrec_to_qc( qrec_bdry , B_stag_x(1,k) ,                 &
                        q_interfaceL(:,1,k) )
                   
                END IF
                
                q_interfaceR(:,1,k) = q_interfaceL(:,1,k)
                
             ENDIF
             
             ! East q_interfaceR(:,comp_interfaces_x,k)
             IF ( j.EQ.comp_cells_x ) THEN
                
                DO i=1,n_vars
                   
                   IF ( bcE(i)%flag .EQ. 0 ) THEN
                      
                      qrec_bdry(i) = bcE(i)%value 
                      
                   ELSE
                      
                      qrec_bdry(i) = qrecE(i)
                      
                   END IF
                   
                ENDDO
                
                IF ( reconstr_variables .EQ. 'phys' ) THEN
                   
                   CALL qp_to_qc( qrec_bdry , B_stag_x(comp_interfaces_x,k) ,   &
                        q_interfaceR(:,comp_interfaces_x,k) )
                   
                ELSEIF ( reconstr_variables .EQ. 'cons' ) THEN
                   
                   CALL qrec_to_qc( qrec_bdry , B_stag_x(comp_interfaces_x,k) , &
                        q_interfaceR(:,comp_interfaces_x,k) )
                   
                END IF
                
                q_interfaceL(:,comp_interfaces_x,k) =                           &
                     q_interfaceR(:,comp_interfaces_x,k)
                
             ENDIF

          END IF
          
       END DO y_loop
       
    END DO x_loop
    
  END SUBROUTINE reconstruction


  !******************************************************************************
  !> \brief Characteristic speeds
  !
  !> This subroutine evaluates the largest characteristic speed at the
  !> cells interfaces from the reconstructed states.
  !> @author 
  !> Mattia de' Michieli Vitturi
  !> \date 16/08/2011
  !******************************************************************************

  SUBROUTINE eval_speeds

    ! External procedures
    USE constitutive_2d, ONLY : eval_local_speeds_x, eval_local_speeds_y 
    USE constitutive_2d, ONLY : eval_local_speeds2_x, eval_local_speeds2_y 

    IMPLICIT NONE

    REAL*8 :: abslambdaL_min(n_vars) , abslambdaL_max(n_vars)
    REAL*8 :: abslambdaR_min(n_vars) , abslambdaR_max(n_vars)
    REAL*8 :: abslambdaB_min(n_vars) , abslambdaB_max(n_vars)
    REAL*8 :: abslambdaT_min(n_vars) , abslambdaT_max(n_vars)
    REAL*8 :: min_r(n_vars) , max_r(n_vars)

    INTEGER :: j,k

    IF ( comp_cells_x .GT. 1 ) THEN
    
       DO j = 1,comp_interfaces_x
          
          DO k = 1, comp_cells_y
             
             CALL eval_local_speeds2_x( q_interfaceL(:,j,k) , B_stag_x(j,k) ,   &
                  abslambdaL_min , abslambdaL_max )
             
             CALL eval_local_speeds2_x( q_interfaceR(:,j,k) , B_stag_x(j,k) ,   &
                  abslambdaR_min , abslambdaR_max )
             
             min_r = MIN(abslambdaL_min , abslambdaR_min , 0.0D0)
             max_r = MAX(abslambdaL_max , abslambdaR_max , 0.0D0)
             
             a_interface_xNeg(:,j,k) = min_r
             a_interface_xPos(:,j,k) = max_r
             
          ENDDO
          
       END DO

    END IF

    IF ( comp_cells_y .GT. 1 ) THEN
    
       DO j = 1,comp_cells_x
          
          DO k = 1,comp_interfaces_y
             
             CALL eval_local_speeds2_y( q_interfaceB(:,j,k) , B_stag_y(j,k) ,   &
                  abslambdaB_min , abslambdaB_max )
             
             CALL eval_local_speeds2_y( q_interfaceT(:,j,k) , B_stag_y(j,k) ,   &
                  abslambdaT_min , abslambdaT_max )
             
             min_r = MIN(abslambdaB_min , abslambdaT_min , 0.0D0)
             max_r = MAX(abslambdaB_max , abslambdaT_max , 0.0D0)
             
             a_interface_yNeg(:,j,k) = min_r
             a_interface_yPos(:,j,k) = max_r
             
          ENDDO
          
       END DO

    END IF
       
  END SUBROUTINE eval_speeds


  !******************************************************************************
  !> \brief Slope limiter
  !
  !> This subroutine limits the slope of the linear reconstruction of 
  !> the physical variables, accordingly to the parameter "solve_limiter":\n
  !> - 'none'     => no limiter (constant value);
  !> - 'minmod'   => minmod slope;
  !> - 'superbee' => superbee limiter (Roe, 1985);
  !> - 'van_leer' => monotonized central-difference limiter (van Leer, 1977)
  !> .
  !> \param[in]     v             3-point stencil value array 
  !> \param[in]     z             3-point stencil location array 
  !> \param[in]     limiter       integer defining the limiter choice
  !> \param[out]    slope_lim     limited slope         
  !> \date 07/10/2016
  !> @author 
  !> Mattia de' Michieli Vitturi
  !******************************************************************************

  SUBROUTINE limit( v , z , limiter , slope_lim )

    USE parameters_2d, ONLY : theta

    IMPLICIT none

    REAL*8, INTENT(IN) :: v(3)
    REAL*8, INTENT(IN) :: z(3)
    INTEGER, INTENT(IN) :: limiter

    REAL*8, INTENT(OUT) :: slope_lim

    REAL*8 :: a , b , c

    REAL*8 :: sigma1 , sigma2

    a = ( v(3) - v(2) ) / ( z(3) - z(2) )
    b = ( v(2) - v(1) ) / ( z(2) - z(1) )
    c = ( v(3) - v(1) ) / ( z(3) - z(1) )

    SELECT CASE (limiter)

    CASE ( 0 )

       slope_lim = 0.D0

    CASE ( 1 )

       ! minmod

       slope_lim = minmod(a,b)

    CASE ( 2 )

       ! superbee

       sigma1 = minmod( a , 2.D0*b )
       sigma2 = minmod( 2.D0*a , b )
       slope_lim = maxmod( sigma1 , sigma2 )

    CASE ( 3 )

       ! generalized minmod

       slope_lim = minmod( c , theta * minmod( a , b ) )

    CASE ( 4 )

       ! monotonized central-difference (MC, LeVeque p.112)

       slope_lim = minmod( c , 2.0 * minmod( a , b ) )

    END SELECT

  END SUBROUTINE limit


  REAL*8 FUNCTION minmod(a,b)

    IMPLICIT none

    REAL*8 :: a , b , sa , sb 

    IF ( a*b .EQ. 0.D0 ) THEN

       minmod = 0.d0

    ELSE

       sa = a / ABS(a)
       sb = b / ABS(b)

       minmod = 0.5D0 * ( sa+sb ) * MIN( ABS(a) , ABS(b) )

    END IF

  END FUNCTION minmod

  REAL*8 function maxmod(a,b)

    IMPLICIT none

    REAL*8 :: a , b , sa , sb 

    IF ( a*b .EQ. 0.d0 ) THEN

       maxmod = 0.d0

    ELSE

       sa = a / ABS(a)
       sb = b / ABS(b)

       maxmod = 0.5D0 * ( sa+sb ) * MAX( ABS(a) , ABS(b) )

    END IF

  END function maxmod

END MODULE solver_2d
