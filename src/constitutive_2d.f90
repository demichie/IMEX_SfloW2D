!********************************************************************************
!> \brief Constitutive equations
!********************************************************************************
MODULE constitutive_2d

  USE parameters_2d, ONLY : n_eqns , n_vars
  USE parameters_2d, ONLY : fischer_flag , rheology_flag

  IMPLICIT none

  LOGICAL, ALLOCATABLE :: implicit_flag(:)

  CHARACTER(LEN=20) :: phase1_name
  CHARACTER(LEN=20) :: phase2_name

  !--------- Constants for the equations of state -----------------------------


  COMPLEX*16 :: h      !< height
  COMPLEX*16 :: u      !< velocity (x direction)
  COMPLEX*16 :: v      !< velocity (y direction)

  !> gravitational acceleration
  REAL*8 :: grav

  !> drag coefficients
  REAL*8 :: mu
  REAL*8 :: xi



CONTAINS

  !******************************************************************************
  !> \brief Initialization of relaxation flags
  !
  !> This subroutine set the number and the flags of the non-hyperbolic
  !> terms.
  !> \date 07/09/2012
  !******************************************************************************

  SUBROUTINE init_problem_param

    USE parameters_2d, ONLY : n_nh

    ALLOCATE( implicit_flag(n_eqns) )

    implicit_flag(1:n_eqns) = .FALSE.
    implicit_flag(2) = .TRUE.
    implicit_flag(3) = .TRUE.

    n_nh = COUNT( implicit_flag )

  END SUBROUTINE init_problem_param

  !******************************************************************************
  !> \brief Physical variables
  !
  !> This subroutine evaluates from the conservative local variables qj
  !> the local physical variables  (\f$\alpha_i, u_i, \rho_i, x_i, \f$)
  !> for the two phases and the mixture denisity \f$ \rho \f$.
  !> \param[in]    r_qj     real conservative variables 
  !> \param[in]    c_qj     complex conservative variables 
  !> \date 15/08/2011
  !******************************************************************************

  SUBROUTINE phys_var(Bj,r_qj,c_qj)
    
    USE COMPLEXIFY
    USE parameters_2d, ONLY : eps_sing
    IMPLICIT none
    
    REAL*8, INTENT(IN) :: Bj
    REAL*8, INTENT(IN), OPTIONAL :: r_qj(n_vars)
    COMPLEX*16, INTENT(IN), OPTIONAL :: c_qj(n_vars)

    COMPLEX*16 :: qj(n_vars)

    IF ( present(c_qj) ) THEN

       qj = c_qj

    ELSE

       qj = DCMPLX(r_qj)

    END IF

    h = qj(1) - DCMPLX( Bj , 0.D0 )

    IF ( REAL( h ) .GT. eps_sing ** 0.25D0 ) THEN

       u = qj(2) / h 

       v = qj(3) / h

    ELSE

       u = DSQRT(2.D0) * h * qj(2) / CDSQRT( h**4 + eps_sing )

       v = DSQRT(2.D0) * h * qj(3) / CDSQRT( h**4 + eps_sing )

    END IF


  END SUBROUTINE phys_var

  !******************************************************************************
  !> \brief Local Characteristic speeds x direction
  !
  !> This subroutine evaluates an the largest positive
  !> and negative characteristic speed for the state qj. 
  !> \date 10/04/2012
  !******************************************************************************

  !SUBROUTINE eval_local_speeds_x(qj,Bj,vel_min,vel_max)
  SUBROUTINE eval_local_speeds_x(qj,Bj,grav3_surf,vel_min,vel_max)
    
    IMPLICIT none
    
    REAL*8, INTENT(IN)  :: qj(n_vars)
    REAL*8, INTENT(IN)  :: Bj
    REAL*8, INTENT(IN)  :: grav3_surf
    REAL*8, INTENT(OUT) :: vel_min(n_vars) , vel_max(n_vars)

    CALL phys_var(Bj,r_qj = qj)
    
    IF (fischer_flag) THEN
    
       vel_min(1:n_eqns) = REAL(u) - DSQRT( - grav * grav3_surf * REAL(h) )
       vel_max(1:n_eqns) = REAL(u) + DSQRT( - grav * grav3_surf * REAL(h) )
       
    ELSE
    
       vel_min(1:n_eqns) = REAL(u) - DSQRT( grav * REAL(h) )
       vel_max(1:n_eqns) = REAL(u) + DSQRT( grav * REAL(h) )
    
    ENDIF
    
  END SUBROUTINE eval_local_speeds_x

  !******************************************************************************
  !> \brief Local Characteristic speeds y direction
  !
  !> This subroutine evaluates an the largest positive
  !> and negative characteristic speed for the state qj. 
  !> \date 10/04/2012
  !******************************************************************************

  SUBROUTINE eval_local_speeds_y(qj,Bj,grav3_surf,vel_min,vel_max)
    
    IMPLICIT none
    
    REAL*8, INTENT(IN)  :: qj(n_vars)
    REAL*8, INTENT(IN)  :: Bj
    REAL*8, INTENT(IN)  :: grav3_surf
    REAL*8, INTENT(OUT) :: vel_min(n_vars) , vel_max(n_vars)
    
    CALL phys_var(Bj,r_qj = qj)
    
    IF (fischer_flag) THEN
    
       vel_min(1:n_eqns) = REAL(v) - DSQRT( - grav * grav3_surf * REAL(h) ) 
       vel_max(1:n_eqns) = REAL(v) + DSQRT( - grav * grav3_surf * REAL(h) )

    ELSE
    
       vel_min(1:n_eqns) = REAL(v) - DSQRT( grav * REAL(h) )
       vel_max(1:n_eqns) = REAL(v) + DSQRT( grav * REAL(h) )
       
    ENDIF

  END SUBROUTINE eval_local_speeds_y

  !******************************************************************************
  !> \brief Local Characteristic speeds
  !
  !> This subroutine evaluates an the largest pos and neg characteristic speeds
  !> from the conservative variables qj, without any change on u and h.
  !> \date 10/04/2012
  !******************************************************************************

  SUBROUTINE eval_local_speeds2_x(qj,Bj,grav3_surf,vel_min,vel_max)
    
    IMPLICIT none
    
    REAL*8, INTENT(IN)  :: qj(n_vars)
    REAL*8, INTENT(IN)  :: Bj
    REAL*8, INTENT(IN)  :: grav3_surf
    REAL*8, INTENT(OUT) :: vel_min(n_vars) , vel_max(n_vars)

    REAL*8 :: h_temp , u_temp

    h_temp = qj(1) - Bj

    IF ( h_temp .NE. 0.D0 ) THEN

       u_temp = qj(2) / h_temp

    ELSE

       u_temp = 0.D0

    END IF

   IF ( fischer_flag ) THEN

        vel_min(1:n_eqns) = REAL(u_temp) - DSQRT( - grav * grav3_surf * REAL(h_temp) ) 
        vel_max(1:n_eqns) = REAL(u_temp) + DSQRT( - grav * grav3_surf * REAL(h_temp) )

    ELSE

        vel_min(1:n_eqns) = REAL(u_temp) - DSQRT( grav * REAL(h_temp) )
        vel_max(1:n_eqns) = REAL(u_temp) + DSQRT( grav * REAL(h_temp) )

    ENDIF

  END SUBROUTINE eval_local_speeds2_x

  !******************************************************************************
  !> \brief Local Characteristic speeds
  !
  !> This subroutine evaluates an the largest pos and neg characteristic speeds
  !> from the conservative variables qj, without any change on v and h.
  !> \date 10/04/2012
  !******************************************************************************

  SUBROUTINE eval_local_speeds2_y(qj,Bj,grav3_surf,vel_min,vel_max)
    
    IMPLICIT none
    
    REAL*8, INTENT(IN)  :: qj(n_vars)
    REAL*8, INTENT(IN)  :: Bj
    REAL*8, INTENT(IN)  :: grav3_surf
    REAL*8, INTENT(OUT) :: vel_min(n_vars) , vel_max(n_vars)

    REAL*8 :: h_temp , v_temp

    h_temp = qj(1) - Bj

    IF ( h_temp .NE. 0.D0 ) THEN

       v_temp = qj(3) / h_temp

    ELSE

       v_temp = 0.D0

    END IF

    IF ( fischer_flag ) THEN

        vel_min(1:n_eqns) = REAL(v_temp) - DSQRT( - grav * grav3_surf * REAL(h_temp) ) 
        vel_max(1:n_eqns) = REAL(v_temp) + DSQRT( - grav * grav3_surf * REAL(h_temp) )

    ELSE

        vel_min(1:n_eqns) = REAL(v_temp) - DSQRT( grav * REAL(h_temp) )
        vel_max(1:n_eqns) = REAL(v_temp) + DSQRT( grav * REAL(h_temp) )

    ENDIF

  END SUBROUTINE eval_local_speeds2_y


  !******************************************************************************
  !> \brief Conservative to physical variables
  !
  !> This subroutine evaluates from the conservative variables qc the 
  !> array of physical variables qp:\n
  !> - qp(1) = \f$ h+B \f$
  !> - qp(2) = \f$ u \f$
  !> .
  !> The physical variables are those used for the linear reconstruction at the
  !> cell interfaces. Limiters are applied to the reconstructed slopes.
  !> \param[in]     qc      conservative variables 
  !> \param[out]    qp      physical variables  
  !> \date 15/08/2011
  !******************************************************************************
  
  SUBROUTINE qc_to_qp(qc,B,qp)
    
    IMPLICIT none
    
    REAL*8, INTENT(IN) :: qc(n_vars)
    REAL*8, INTENT(IN) :: B
    REAL*8, INTENT(OUT) :: qp(n_vars)
    
    CALL phys_var(B,r_qj = qc)
        
    qp(1) = REAL(h+B)
    qp(2) = REAL(u)
    qp(3) = REAL(v)

    
  END SUBROUTINE qc_to_qp

  !******************************************************************************
  !> \brief Physical to conservative variables
  !
  !> This subroutine evaluates the conservative variables qc from the 
  !> array of physical variables qp:\n
  !> - qp(1) = \f$ h + B \f$
  !> - qp(2) = \f$ hu \f$
  !> .
  !> \param[in]    qp      physical variables  
  !> \param[out]   qc      conservative variables 
  !> \date 15/08/2011
  !******************************************************************************

  SUBROUTINE qp_to_qc(qp,B,qc)
    
    USE COMPLEXIFY 
    IMPLICIT none
    
    REAL*8, INTENT(IN) :: qp(n_vars)
    REAL*8, INTENT(IN) :: B
    REAL*8, INTENT(OUT) :: qc(n_vars)
    
    REAL*8 :: r_hB      !> topography + height 
    REAL*8 :: r_u       !> velocity
    REAL*8 :: r_v       !> velocity
    
    r_hB = qp(1)
    r_u = qp(2)
    r_v = qp(3)
    
    qc(1) = r_hB
    qc(2) = ( r_hB - B ) * r_u
    qc(3) = ( r_hB - B ) * r_v

  END SUBROUTINE qp_to_qc

  !******************************************************************************
  !> \brief Conservative to 2nd set of physical variables
  !
  !> This subroutine evaluates from the conservative variables qc the 
  !> array of physical variables qp:\n
  !> - qp(1) = \f$ h \f$
  !> - qp(2) = \f$ u \f$
  !> - qp(3) = \f$ v \f$
  !> .
  !> The physical variables are those used for the linear reconstruction at the
  !> cell interfaces. Limiters are applied to the reconstructed slopes.
  !> \param[in]     qc      conservative variables 
  !> \param[out]    qp      physical variables  
  !> \date 15/08/2011
  !******************************************************************************
  
  SUBROUTINE qc_to_qp2(qc,B,qp)
    
    IMPLICIT none
    
    REAL*8, INTENT(IN) :: qc(n_vars)
    REAL*8, INTENT(IN) :: B
    REAL*8, INTENT(OUT) :: qp(n_vars)
    
    CALL phys_var(B,r_qj = qc)
          
    qp(1) = REAL(h)
    qp(2) = REAL(u)
    qp(3) = REAL(v)

    
  END SUBROUTINE qc_to_qp2

  !******************************************************************************
  !> \brief From 2nd set of physical to conservative variables
  !
  !> This subroutine evaluates the conservative variables qc from the 
  !> array of physical variables qp:\n
  !> - qp(1) = \f$ h \f$
  !> - qp(2) = \f$ u \f$
  !> - qp(2) = \f$ v \f$
  !> .
  !> \param[in]    qp      physical variables  
  !> \param[out]   qc      conservative variables 
  !> \date 15/08/2011
  !******************************************************************************

  SUBROUTINE qp2_to_qc(qp,B,qc)
    
    USE COMPLEXIFY 
    IMPLICIT none
    
    REAL*8, INTENT(IN) :: qp(n_vars)
    REAL*8, INTENT(IN) :: B
    REAL*8, INTENT(OUT) :: qc(n_vars)
    
    REAL*8 :: r_h      !> topography + height 
    REAL*8 :: r_u      !> velocity x direction
    REAL*8 :: r_v      !> velocity y direction
    
    r_h = qp(1)
    r_u = qp(2)
    r_v = qp(3)
    
    qc(1) = r_h + B
    qc(2) = r_h * r_u
    qc(3) = r_h * r_v

  END SUBROUTINE qp2_to_qc


  !******************************************************************************
  !> \brief Hyperbolic Fluxes
  !
  !> This subroutine evaluates the numerical fluxes given the conservative
  !> variables qj, accordingly to the equations for the single temperature
  !> model introduced in Romenki et al. 2010.
  !> \date 01/06/2012
  !> \param[in]     c_qj     complex conservative variables 
  !> \param[in]     r_qj     real conservative variables 
  !> \param[out]    c_flux   complex analytical fluxes    
  !> \param[out]    r_flux   real analytical fluxes    
  !******************************************************************************
  
  !SUBROUTINE eval_fluxes(Bj,c_qj,r_qj,c_flux,r_flux,dir)
  SUBROUTINE eval_fluxes(Bj,grav3_surf,c_qj,r_qj,c_flux,r_flux,dir)
    
    USE COMPLEXIFY
    IMPLICIT none

    REAL*8, INTENT(IN) :: Bj
    REAL*8, INTENT(IN) :: grav3_surf
    COMPLEX*16, INTENT(IN), OPTIONAL :: c_qj(n_vars)
    COMPLEX*16, INTENT(OUT), OPTIONAL :: c_flux(n_eqns)
    REAL*8, INTENT(IN), OPTIONAL :: r_qj(n_vars)
    REAL*8, INTENT(OUT), OPTIONAL :: r_flux(n_eqns)
    INTEGER, INTENT(IN) :: dir

    COMPLEX*16 :: qj(n_vars)
    COMPLEX*16 :: flux(n_eqns)
    COMPLEX*16 :: h_temp , u_temp , v_temp

    INTEGER :: i 

    IF ( present(c_qj) .AND. present(c_flux) ) THEN

       qj = c_qj

    ELSEIF ( present(r_qj) .AND. present(r_flux) ) THEN

       DO i = 1,n_vars

          qj(i) = DCMPLX( r_qj(i) )

       END DO

    ELSE

       WRITE(*,*) 'Constitutive, eval_fluxes: problem with arguments'
       STOP

    END IF


    IF ( dir .EQ. 1 ) THEN

    ! flux F (derivated wrt x in the equations)

       flux(1)=qj(2)

       h_temp = qj(1) - Bj

       IF ( REAL(h_temp) .NE. 0.D0 ) THEN

          u_temp = qj(2) / h_temp
          
          IF ( fischer_flag ) THEN

             flux(2) = h_temp * u_temp**2 - 0.5D0 * grav * grav3_surf * h_temp**2  

          ELSE
          
             flux(2) = h_temp * u_temp**2 + 0.5D0 * grav * h_temp**2  
                     
          ENDIF

          flux(3) = u_temp * qj(3)

       ELSE

          flux(2) = 0.D0

          flux(3) = 0.D0

       ENDIF

    ELSEIF ( dir .EQ. 2 ) THEN

       ! flux G (derivated wrt y in the equations)

       flux(1) = qj(3)

       h_temp = qj(1) - Bj

       IF(REAL(h_temp).NE.0.d0)THEN

          v_temp = qj(3) / h_temp

          flux(2) = v_temp * qj(2)

          IF ( fischer_flag ) THEN

             flux(3) = h_temp * v_temp**2 - 0.5D0 * grav * grav3_surf * h_temp**2
             
          ELSE
          
             flux(3) = h_temp * v_temp**2 + 0.5D0 * grav * h_temp**2

          ENDIF

       ELSE

          flux(2) = 0.0

          flux(3) = 0.0

       ENDIF

    ELSE

       WRITE(*,*) 'Constitutive, eval_fluxes: problem with arguments'

       STOP

    ENDIF

    IF ( present(c_qj) .AND. present(c_flux) ) THEN

       c_flux = flux

    ELSEIF ( present(r_qj) .AND. present(r_flux) ) THEN

       r_flux = REAL( flux )

    END IF

  END SUBROUTINE eval_fluxes

  !******************************************************************************
  !> \brief Non-Hyperbolic terms
  !
  !> This subroutine evaluates the non-hyperbolic terms (relaxation terms
  !> and forces) of the system of equations, both for real or complex 
  !> inputs. These terms are treated implicitely in the DIRK numerical
  !> scheme.
  !> \date 01/06/2012
  !> \param[in]     c_qj            complex conservative variables 
  !> \param[in]     r_qj            real conservative variables 
  !> \param[out]    c_nh_term_impl  complex non-hyperbolic terms     
  !> \param[out]    r_nh_term_impl  real non-hyperbolic terms
  !******************************************************************************

  SUBROUTINE eval_nonhyperbolic_terms( Bj , Bprimej_x , Bprimej_y , grav3_surf ,&
       curvj_x , curvj_y , c_qj , c_nh_term_impl , r_qj , r_nh_term_impl )

    USE COMPLEXIFY 
    IMPLICIT NONE

    REAL*8, INTENT(IN) :: Bj
    REAL*8, INTENT(IN) :: Bprimej_x, Bprimej_y
    REAL*8, INTENT(IN) :: grav3_surf
    REAL*8, INTENT(IN) :: curvj_x
    REAL*8, INTENT(IN) :: curvj_y

    COMPLEX*16, INTENT(IN), OPTIONAL :: c_qj(n_vars)
    COMPLEX*16, INTENT(OUT), OPTIONAL :: c_nh_term_impl(n_eqns)
    REAL*8, INTENT(IN), OPTIONAL :: r_qj(n_vars)
    REAL*8, INTENT(OUT), OPTIONAL :: r_nh_term_impl(n_eqns)

    COMPLEX*16 :: qj(n_vars)

    COMPLEX*16 :: nh_term(n_eqns)

    COMPLEX*16 :: relaxation_term(n_eqns)

    COMPLEX*16 :: forces_term(n_eqns)

    INTEGER :: i

    COMPLEX*16 :: mod_vel
    
    IF ( present(c_qj) .AND. present(c_nh_term_impl) ) THEN

       qj = c_qj

    ELSEIF ( present(r_qj) .AND. present(r_nh_term_impl) ) THEN

       DO i = 1,n_vars

          qj(i) = DCMPLX( r_qj(i) )

       END DO

    ELSE

       WRITE(*,*) 'Constitutive, eval_fluxes: problem with arguments'
       STOP

    END IF

    ! initialize and evaluate the relaxation terms
    relaxation_term(1:n_eqns) = DCMPLX(0.D0) 

    ! initialize and evaluate the forces terms
    forces_term(1:n_eqns) = DCMPLX(0.D0)

    IF (rheology_flag) THEN

       CALL phys_var(Bj,c_qj = qj)
    
       mod_vel = CDSQRT( u**2 + v**2 )
    
       IF ( REAL(mod_vel) .NE. 0.D0 ) THEN 
     
          forces_term(2) = forces_term(2) - (u/mod_vel) *                          &
             ( mu*h * ( - grav * grav3_surf + curvj_x * mod_vel**2 )               &
             + ( grav / xi ) * mod_vel ** 2 )
     
          forces_term(3) = forces_term(3) - (v/mod_vel) *                          &
             ( mu*h * ( - grav * grav3_surf + curvj_y * mod_vel**2 )               &
             + ( grav / xi ) * mod_vel ** 2 )
     
       ENDIF
       
    ENDIF

    nh_term = relaxation_term + forces_term

    IF ( present(c_qj) .AND. present(c_nh_term_impl) ) THEN

       c_nh_term_impl = nh_term

    ELSEIF ( present(r_qj) .AND. present(r_nh_term_impl) ) THEN

       r_nh_term_impl = REAL( nh_term )

    END IF 

  END SUBROUTINE eval_nonhyperbolic_terms

  !******************************************************************************
  !> \brief Explicit Forces term
  !
  !> This subroutine evaluates the forces to be treated explicitely
  !> in the DIRK numerical scheme (e.g. gravity)
  !> \date 01/06/2012
  !> \param[in]     qj                 conservative variables 
  !> \param[out]    expl_forces_term   forces term
  !******************************************************************************

  SUBROUTINE eval_explicit_forces( Bj , Bprimej_x , Bprimej_y , gravj_surf ,    &
       curvj_x , curvj_y , qj , expl_forces_term )
    
    IMPLICIT NONE

    REAL*8, INTENT(IN) :: Bj
    REAL*8, INTENT(IN) :: Bprimej_x
    REAL*8, INTENT(IN) :: Bprimej_y
    REAL*8, INTENT(IN) :: gravj_surf(3)
    REAL*8, INTENT(IN) :: curvj_x
    REAL*8, INTENT(IN) :: curvj_y
    !REAL*8, INTENT(IN) :: curvj_xy

    REAL*8, INTENT(IN) :: qj(n_eqns)                 !< conservative variables 
    REAL*8, INTENT(OUT) :: expl_forces_term(n_eqns)  !< explicit forces 

    expl_forces_term(1:n_eqns) = 0.D0

    CALL phys_var(Bj,r_qj = qj)
    
    IF ( fischer_flag ) THEN

       expl_forces_term(2) = grav * h * gravj_surf(1)
   
       expl_forces_term(3) = grav * h * gravj_surf(2)
    
    ELSE
    
       expl_forces_term(2) = grav * h * Bprimej_x
   
       expl_forces_term(3) = grav * h * Bprimej_y
       
    ENDIF
    
  END SUBROUTINE eval_explicit_forces

END MODULE constitutive_2d

    
