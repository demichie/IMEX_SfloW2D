!********************************************************************************
!> \brief Constitutive equations
!********************************************************************************
MODULE constitutive_2d

  USE parameters_2d, ONLY : n_eqns , n_vars
  USE parameters_2d, ONLY : rheology_flag , rheology_model
  USE parameters_2d, ONLY : temperature_flag

  IMPLICIT none

  LOGICAL, ALLOCATABLE :: implicit_flag(:)

  CHARACTER(LEN=20) :: phase1_name
  CHARACTER(LEN=20) :: phase2_name


  COMPLEX*16 :: h      !< height [m]
  COMPLEX*16 :: u      !< velocity (x direction)
  COMPLEX*16 :: v      !< velocity (y direction)
  COMPLEX*16 :: T      !< temperature

  !> gravitational acceleration
  REAL*8 :: grav

  !> drag coefficients (Voellmy-Salm model)
  REAL*8 :: mu
  REAL*8 :: xi
  
  !> drag coefficients (plastic model)
  REAL*8 :: tau

  !> evironment temperature [K]
  REAL*8 :: T_env

  !> radiative coefficient
  REAL*8 :: rad_coeff

  !> friction coefficient
  REAL*8 :: frict_coeff

  !> fluid density [kg/m3]
  REAL*8 :: rho

  !> reference temperature [K]
  REAL*8 :: T_ref

  !> reference kinematic viscosity [m2/s]
  REAL*8 :: nu_ref

  !> viscosity parameter [K-1] (b in Table 1 Costa & Macedonio, 2005)
  REAL*8 :: visc_par

  !> velocity boundary layer fraction of total thickness
  REAL*8 :: emme

  !> specific heat [J kg-1 K-1]
  REAL*8 :: c_p

  !> atmospheric heat trasnfer coefficient [W m-2 K-1] (lambda in C&M, 2005)
  REAL*8 :: atm_heat_transf_coeff

  !> fractional area of the exposed inner core (f in C&M, 2005)
  REAL*8 :: exp_area_fract

  !> Stephan-Boltzmann constant [W m-2 K-4]
  REAL*8, PARAMETER :: SBconst = 5.67D-8

  !> emissivity (eps in Costa & Macedonio, 2005)
  REAL*8 :: emissivity

  !> thermal boundary layer fraction of total thickness
  REAL*8 :: enne

  !> temperature of lava-ground interface
  REAL*8 :: T_ground

  !> thermal conductivity [W m-1 K-1] (k in Costa & Macedonio, 2005)
  REAL*8 :: thermal_conductivity

  !--- Lahars rheology model parameters

  !> 1st parameter for yield strenght empirical relationship (O'Brian et al, 1993)
  REAL*8 :: alpha2

  !> 2nd parameter for yield strenght empirical relationship (O'Brian et al, 1993)
  REAL*8 :: beta2

  !> 1st parameter for fluid viscosity empirical relationship (O'Brian et al, 1993)
  REAL*8 :: alpha1

  !> 2nd parameter for fluid viscosity empirical relationship (O'Brian et al, 1993)
  REAL*8 :: beta1

  !> Empirical resistance parameter
  REAL*8 :: Kappa

  !> Mannings roughness coefficient ( units: T L^(-1/3) )
  REAL*8 :: n_td
 
  !> Specific weight of water
  REAL*8 :: gamma_w

  !> Specific weight of sediments
  REAL*8 :: gamma_s

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

    IF ( temperature_flag ) implicit_flag(4) = .TRUE.


    n_nh = COUNT( implicit_flag )

  END SUBROUTINE init_problem_param

  !******************************************************************************
  !> \brief Physical variables
  !
  !> This subroutine evaluates from the conservative local variables qj
  !> the local physical variables  (\f$h+B, u, v, T \f$).
  !> \param[in]    r_qj     real conservative variables 
  !> \param[in]    c_qj     complex conservative variables 
  !> @author 
  !> Mattia de' Michieli Vitturi
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

       IF ( temperature_flag ) T = qj(4) / h

    ELSE

       u = DSQRT(2.D0) * h * qj(2) / CDSQRT( h**4 + eps_sing )

       v = DSQRT(2.D0) * h * qj(3) / CDSQRT( h**4 + eps_sing )

       IF ( temperature_flag ) THEN

          T =  DSQRT(2.D0) * h * qj(4) / CDSQRT( h**4 + eps_sing )

       END IF

    END IF
    
  END SUBROUTINE phys_var

  !******************************************************************************
  !> \brief Local Characteristic speeds x direction
  !
  !> This subroutine desingularize the velocities and then evaluates the largest 
  !> positive and negative characteristic speed in the x-direction. 
  !> \param[in]     qj            array of conservative variables
  !> \param[in]     Bj            topography at the cell center
  !> \param[out]    vel_min       minimum x-velocity
  !> \param[out]    vel_max       maximum x-velocity
  !> @author 
  !> Mattia de' Michieli Vitturi
  !> \date 05/12/2017
  !******************************************************************************

  SUBROUTINE eval_local_speeds_x(qj,Bj,vel_min,vel_max)
  
    IMPLICIT none
    
    REAL*8, INTENT(IN)  :: qj(n_vars)
    REAL*8, INTENT(IN)  :: Bj
    REAL*8, INTENT(OUT) :: vel_min(n_vars) , vel_max(n_vars)

    CALL phys_var(Bj,r_qj = qj)
    
    vel_min(1:n_eqns) = REAL(u) - DSQRT( grav * REAL(h) )
    vel_max(1:n_eqns) = REAL(u) + DSQRT( grav * REAL(h) )
        
  END SUBROUTINE eval_local_speeds_x

  !******************************************************************************
  !> \brief Local Characteristic speeds y direction
  !
  !> This subroutine desingularize the velocities and then evaluates the largest 
  !> positive and negative characteristic speed in the y-direction. 
  !> \param[in]     qj            array of conservative variables
  !> \param[in]     Bj            topography at the cell center
  !> \param[out]    vel_min       minimum y-velocity
  !> \param[out]    vel_max       maximum y-velocity
  !> @author 
  !> Mattia de' Michieli Vitturi
  !> \date 05/12/2017
  !******************************************************************************

  SUBROUTINE eval_local_speeds_y(qj,Bj,vel_min,vel_max)
    
    IMPLICIT none
    
    REAL*8, INTENT(IN)  :: qj(n_vars)
    REAL*8, INTENT(IN)  :: Bj
    REAL*8, INTENT(OUT) :: vel_min(n_vars) , vel_max(n_vars)
    
    CALL phys_var(Bj,r_qj = qj)
    
    vel_min(1:n_eqns) = REAL(v) - DSQRT( grav * REAL(h) )
    vel_max(1:n_eqns) = REAL(v) + DSQRT( grav * REAL(h) )
       
  END SUBROUTINE eval_local_speeds_y

  !******************************************************************************
  !> \brief Local Characteristic speeds
  !
  !> This subroutine evaluates an the largest pos and neg characteristic speeds
  !> from the conservative variables qj, without any change on u and h.
  !> \param[in]     qj            array of conservative variables
  !> \param[in]     Bj            topography at the cell center
  !> \param[out]    vel_min       minimum x-velocity
  !> \param[out]    vel_max       maximum x-velocity
  !> @author 
  !> Mattia de' Michieli Vitturi
  !> \date 05/12/2017
  !******************************************************************************

  SUBROUTINE eval_local_speeds2_x(qj,Bj,vel_min,vel_max)
    
    IMPLICIT none
    
    REAL*8, INTENT(IN)  :: qj(n_vars)
    REAL*8, INTENT(IN)  :: Bj
    REAL*8, INTENT(OUT) :: vel_min(n_vars) , vel_max(n_vars)

    REAL*8 :: h_temp , u_temp

    h_temp = qj(1) - Bj

    IF ( h_temp .NE. 0.D0 ) THEN

       u_temp = qj(2) / h_temp

    ELSE

       u_temp = 0.D0

    END IF

    vel_min(1:n_eqns) = u_temp - DSQRT( grav * h_temp )
    vel_max(1:n_eqns) = u_temp + DSQRT( grav * h_temp )

  END SUBROUTINE eval_local_speeds2_x

  !******************************************************************************
  !> \brief Local Characteristic speeds
  !
  !> This subroutine evaluates an the largest pos and neg characteristic speeds
  !> from the conservative variables qj, without any change on v and h.
  !> \param[in]     qj            array of conservative variables
  !> \param[in]     Bj            topography at the cell center
  !> \param[out]    vel_min       minimum y-velocity
  !> \param[out]    vel_max       maximum y-velocity
  !> @author 
  !> Mattia de' Michieli Vitturi
  !> \date 05/12/2017
  !******************************************************************************

  SUBROUTINE eval_local_speeds2_y(qj,Bj,vel_min,vel_max)
    
    IMPLICIT none
    
    REAL*8, INTENT(IN)  :: qj(n_vars)
    REAL*8, INTENT(IN)  :: Bj
    REAL*8, INTENT(OUT) :: vel_min(n_vars) , vel_max(n_vars)

    REAL*8 :: h_temp , v_temp

    h_temp = qj(1) - Bj

    IF ( h_temp .NE. 0.D0 ) THEN

       v_temp = qj(3) / h_temp

    ELSE

       v_temp = 0.D0

    END IF

    vel_min(1:n_eqns) = v_temp - DSQRT( grav * h_temp )
    vel_max(1:n_eqns) = v_temp + DSQRT( grav * h_temp )
    
  END SUBROUTINE eval_local_speeds2_y
  
  !******************************************************************************
  !> \brief Conservative to physical variables
  !
  !> This subroutine evaluates from the conservative variables qc the 
  !> array of physical variables qp:\n
  !> - qp(1) = \f$ h+B \f$
  !> - qp(2) = \f$ u \f$
  !> - qp(3) = \f$ v \f$
  !> - qp(4) = \f$ T \f$
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

    IF ( temperature_flag ) qp(4) = REAL(T)

    
  END SUBROUTINE qc_to_qp

  !******************************************************************************
  !> \brief Physical to conservative variables
  !
  !> This subroutine evaluates the conservative variables qc from the 
  !> array of physical variables qp:\n
  !> - qp(1) = \f$ h + B \f$
  !> - qp(2) = \f$ u \f$
  !> - qp(3) = \f$ v \f$
  !> - qp(4) = \f$ T \f$
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
    REAL*8 :: r_T       !> temperature

    r_hB = qp(1)
    r_u  = qp(2)
    r_v  = qp(3)
    r_T  = qp(4)

    qc(1) = r_hB
    qc(2) = ( r_hB - B ) * r_u
    qc(3) = ( r_hB - B ) * r_v

    IF ( temperature_flag ) qc(4) = ( r_hB - B ) * r_T 

  END SUBROUTINE qp_to_qc

  !******************************************************************************
  !> \brief Reconstructed to conservative variables
  !
  !> This subroutine evaluates the conservative variables qc from the 
  !> array of physical variables qrec:\n
  !> - qrec(1) = \f$ h + B \f$
  !> - qrec(2) = \f$ hu \f$
  !> - qrec(3) = \f$ hv \f$
  !> - qrec(4) = \f$ h \cdot xs \f$
  !> - qrec(5) = \f$ T \f$
  !> .
  !> \param[in]    qrec      physical variables  
  !> \param[out]   qc      conservative variables 
  !> \date 15/08/2011
  !******************************************************************************

  SUBROUTINE qrec_to_qc(qrec,B,qc)
    
    IMPLICIT none
    
    REAL*8, INTENT(IN) :: qrec(n_vars)
    REAL*8, INTENT(IN) :: B
    REAL*8, INTENT(OUT) :: qc(n_vars)

    REAL*8 :: r_hB      !> topography + height 
    REAL*8 :: r_u       !> velocity
    REAL*8 :: r_v       !> velocity
    REAL*8 :: r_T       !> temperature

    ! Desingularization
    CALL phys_var(B,r_qj = qrec)
          
    r_hB = REAL(h) + B
    r_u = REAL(u)
    r_v = REAL(v)

    IF ( temperature_flag ) r_T = REAL(T)

    qc(1) = r_hB
    qc(2) = REAL(h) * r_u
    qc(3) = REAL(h) * r_v
    IF ( temperature_flag ) qc(4) = REAL(h) * r_T 
    
  END SUBROUTINE qrec_to_qc


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
  
  SUBROUTINE eval_fluxes(Bj,c_qj,r_qj,c_flux,r_flux,dir)
    
    USE COMPLEXIFY
    IMPLICIT none

    REAL*8, INTENT(IN) :: Bj
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

       flux(1) = qj(2)

       h_temp = qj(1) - Bj

       IF ( REAL(h_temp) .NE. 0.D0 ) THEN

          u_temp = qj(2) / h_temp
          
          flux(2) = h_temp * u_temp**2 + 0.5D0 * grav * h_temp**2  

          flux(3) = u_temp * qj(3)

          ! Temperature flux in x-direction: U*T*h
          IF ( temperature_flag ) flux(4) = u_temp * qj(4)

       ELSE

          flux(2) = 0.D0

          flux(3) = 0.D0

          IF ( temperature_flag ) flux(4) = 0.D0

       ENDIF

    ELSEIF ( dir .EQ. 2 ) THEN

       ! flux G (derivated wrt y in the equations)

       flux(1) = qj(3)

       h_temp = qj(1) - Bj

       IF(REAL(h_temp).NE.0.d0)THEN

          v_temp = qj(3) / h_temp

          flux(2) = v_temp * qj(2)
          
          flux(3) = h_temp * v_temp**2 + 0.5D0 * grav * h_temp**2

          ! Temperature flux in x-direction: V*T*h
          IF ( temperature_flag ) flux(4) = v_temp * qj(4)

       ELSE

          flux(2) = 0.D0

          flux(3) = 0.D0

          IF ( temperature_flag ) flux(4) = 0.D0

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
       c_qj , c_nh_term_impl , r_qj , r_nh_term_impl )

    USE COMPLEXIFY 
    USE parameters_2d, ONLY : sed_vol_fract
    IMPLICIT NONE

    REAL*8, INTENT(IN) :: Bj
    REAL*8, INTENT(IN) :: Bprimej_x, Bprimej_y
    REAL*8, INTENT(IN) :: grav3_surf

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
    
    COMPLEX*16 :: gamma

    REAL*8 :: radiative_coeff

    COMPLEX*16 :: radiative_term

    REAL*8 :: convective_coeff

    COMPLEX*16 :: convective_term

    COMPLEX*16 :: conductive_coeff , conductive_term

    REAL*8 :: thermal_diffusivity

    REAL*8 :: h_threshold

    !--- Lahars rheology model variables
    
    !> Yield strenght
    COMPLEX*8 :: tau_y

    !> Fluid viscosity
    COMPLEX*8 :: fluid_visc

    !> Sediment volume fraction
    COMPLEX*8 :: sed_vol_fract_cmplx

    !> Specific weight of sediment mixture
    COMPLEX*8 :: gamma_m

    !> Total friction
    COMPLEX*8 :: s_f

    !> Yield slope component of total friction
    COMPLEX*8 :: s_y

    !> Viscous slope component of total Friction
    COMPLEX*8 :: s_v

    !> Turbulent dispersive slope component of total friction
    COMPLEX*8 :: s_td


    IF ( temperature_flag ) THEN

       IF ( ( thermal_conductivity .GT. 0.D0 ) .OR. ( emme .GT. 0.D0 ) ) THEN

          h_threshold = 1.D-10

       ELSE

          h_threshold = 0.D0

       END IF
       
    END IF
       

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
    relaxation_term(1:n_eqns) = DCMPLX(0.D0,0.D0) 

    ! initialize and evaluate the forces terms
    forces_term(1:n_eqns) = DCMPLX(0.D0,0.D0)

    IF (rheology_flag) THEN

       CALL phys_var(Bj,c_qj = qj)
    
       mod_vel = CDSQRT( u**2 + v**2 )
       
       ! Voellmy Salm rheology
       IF ( rheology_model .EQ. 1 ) THEN
       
          IF ( REAL(mod_vel) .NE. 0.D0 ) THEN 
          
             forces_term(2) = forces_term(2) -  ( u / mod_vel ) *               &
                  ( mu * h * ( - grav * grav3_surf )                            &
                  + ( grav / xi ) * mod_vel ** 2 )
          
             forces_term(3) = forces_term(3) -  ( v / mod_vel ) *               &
                  ( mu * h * ( - grav * grav3_surf )                            &
                  + ( grav / xi ) * mod_vel ** 2 )
          
          ENDIF
        
       ! Plastic rheology
       ELSEIF ( rheology_model .EQ. 2 ) THEN
       
          IF ( REAL(mod_vel) .NE. 0.D0 ) THEN 
          
             forces_term(2) = forces_term(2) - tau * (u/mod_vel)
          
             forces_term(3) = forces_term(3) - tau * (v/mod_vel)

          ENDIF

       ! Temperature dependent rheology
       ELSEIF ( rheology_model .EQ. 3 ) THEN

          IF ( REAL(h) .GT. h_threshold ) THEN
    
             ! Equation 6 from Costa & Macedonio, 2005
             gamma = 3.D0 * nu_ref / h * CDEXP( - visc_par * ( T - T_ref ) )

          ELSE

             ! Equation 6 from Costa & Macedonio, 2005
             gamma = 3.D0 * nu_ref / h_threshold * CDEXP( - visc_par            &
                  * ( T - T_ref ) )
             
          END IF
          
          IF ( REAL(mod_vel) .NE. 0.D0 ) THEN 
          
             ! Last R.H.S. term in equation 2 from Costa & Macedonio, 2005
             forces_term(2) = forces_term(2) - gamma * u
          
             ! Last R.H.S. term in equation 3 from Costa & Macedonio, 2005
             forces_term(3) = forces_term(3) - gamma * v

          ENDIF
          
       ! Lahars rheology (O'Brien 1993, FLO2D)
       ELSEIF ( rheology_model .EQ. 4 ) THEN

          h_threshold = 1.D-20

          sed_vol_fract_cmplx = DCMPLX(sed_vol_fract,0.D0)

          ! Convert from mass fraction to volume fraction
          ! sed_vol_fract = xs * gamma_w / ( xs * gamma_w + (  DCMPLX(1.D0,0.D0) - xs ) * gamma_s )

          !IF ( xs .NE. 0.D0 ) THEN
             
             !WRITE(*,*) 'xs',xs
             !WRITE(*,*) 'sed_vol_fract',sed_vol_fract
             !READ(*,*) 
             
          !END IF


          ! Mixture density
          gamma_m = ( DCMPLX(1.D0,0.D0) - sed_vol_fract_cmplx ) * gamma_w       &
               + sed_vol_fract_cmplx * gamma_s 

          ! Yield strength
          tau_y = alpha2 * CDEXP( beta2 * sed_vol_fract_cmplx )

          ! Fluid viscosity
          fluid_visc = alpha1 * CDEXP( beta1 * sed_vol_fract_cmplx )


          IF ( h .GT. h_threshold ) THEN
             
             ! Yield slope component
             s_y = tau_y / ( gamma_m * h )
             
             ! Viscous slope component
             s_v = Kappa * fluid_visc * mod_vel / ( 8.D0 * gamma_m * h**2 )
             
             
             ! Turbulent dispersive component
             s_td = n_td**2 * mod_vel**2 / ( h**(4.D0/3.D0) )
             
          ELSE
             
             ! Yield slope component
             s_y = tau_y / ( gamma_m * h_threshold )
             
             ! Viscous slope component
             s_v = Kappa * fluid_visc * mod_vel / ( 8.D0 * gamma_m *            &
                  h_threshold**2 )
             
             ! Turbulent dispersive components
             s_td = n_td**2 * mod_vel**2 / ( h_threshold**(4.D0/3.D0) )
             
          END IF
          
          ! Total friction slope
          s_f = s_y + s_v + s_td
          
          IF ( mod_vel .GT. 0.D0 ) THEN

             forces_term(2) = forces_term(2) - grav * h * ( u / mod_vel ) * s_f
             forces_term(3) = forces_term(3) - grav * h * ( v / mod_vel ) * s_f
  
          END IF

       ELSEIF ( rheology_model .EQ. 5 ) THEN

          tau = 1.D-3 / ( 1.D0 + 10.D0 * h ) * mod_vel
          
          IF ( REAL(mod_vel) .NE. 0.D0 ) THEN
             
             forces_term(2) = forces_term(2) - tau * ( u / mod_vel )
             forces_term(3) = forces_term(3) - tau * ( v / mod_vel )

          END IF
          
       ENDIF
              
    ENDIF

    IF ( temperature_flag ) THEN

       CALL phys_var(Bj,c_qj = qj)

       IF ( REAL(h) .GT. 0.d0 ) THEN

          ! Equation 8 from Costa & Macedonio, 2005
          radiative_coeff = emissivity * SBconst * exp_area_fract / ( rho * c_p )
       
       ELSE

          radiative_coeff = 0.D0

       END IF

       IF ( REAL(T) .GT. T_env ) THEN

          ! First R.H.S. term in equation 4 from Costa & Macedonio, 2005
          radiative_term = - radiative_coeff * ( T**4 - T_env**4 )

       ELSE

          radiative_term = DCMPLX(0.D0,0.D0)

       END IF

       IF ( REAL(h) .GT. 0.d0 ) THEN
       
          ! Equation 9 from Costa & Macedonio, 2005
          convective_coeff = atm_heat_transf_coeff * exp_area_fract             &
               / ( rho * c_p )

       ELSE

          convective_coeff = 0.D0

       END IF

       IF ( REAL(T) .GT. T_env ) THEN

          ! Second R.H.S. term in equation 4 from Costa & Macedonio, 2005
          convective_term = - convective_coeff * ( T - T_env )

       ELSE

          convective_term =  DCMPLX(0.D0,0.D0)

       END IF

       IF ( REAL(h) .GT. h_threshold ) THEN
    
          thermal_diffusivity = thermal_conductivity / ( rho * c_p ) 

          ! Equation 7 from Costa & Macedonio, 2005
          conductive_coeff = enne * thermal_diffusivity / h

       ELSE

          conductive_coeff =  DCMPLX(0.D0,0.D0)
          conductive_coeff = enne * thermal_diffusivity / DCMPLX(h_threshold,0.D0)

       END IF

       ! Third R.H.S. term in equation 4 from Costa & Macedonio, 2005
       IF ( REAL(T) .GT. T_ground ) THEN

          conductive_term = - conductive_coeff * ( T - T_ground )

       ELSE

           conductive_term = DCMPLX(0.D0,0.D0)

        END IF

       relaxation_term(4) = radiative_term + convective_term + conductive_term

    END IF

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
  !> This subroutine evaluates the non-hyperbolic terms to be treated explicitely
  !> in the DIRK numerical scheme (e.g. gravity,source of mass). The sign of the
  !> terms is taken with the terms on the left-hand side of the equations.
  !> \date 01/06/2012
  !> \param[in]     qj                 conservative variables 
  !> \param[out]    expl_term          explicit term
  !******************************************************************************

  SUBROUTINE eval_expl_terms( Bj, Bprimej_x , Bprimej_y , source_xy , qj ,      &
       expl_term )

    USE parameters_2d, ONLY : source_flag , vel_source , T_source
    
    IMPLICIT NONE

    REAL*8, INTENT(IN) :: Bj
    REAL*8, INTENT(IN) :: Bprimej_x
    REAL*8, INTENT(IN) :: Bprimej_y
    REAL*8, INTENT(IN) :: source_xy
    
    REAL*8, INTENT(IN) :: qj(n_eqns)                 !< conservative variables 
    REAL*8, INTENT(OUT) :: expl_term(n_eqns)         !< explicit forces 

    REAL*8 :: visc_heat_coeff

    expl_term(1:n_eqns) = 0.D0

    CALL phys_var(Bj,r_qj = qj)

    IF ( source_flag ) expl_term(1) = - source_xy * vel_source
    
    expl_term(2) = grav * REAL(h) * Bprimej_x
   
    expl_term(3) = grav * REAL(h) * Bprimej_y

    IF ( temperature_flag .AND. source_flag ) THEN
    
       expl_term(4) = - source_xy * vel_source * T_source

       IF ( rheology_model .EQ. 3 ) THEN
              
          IF ( REAL(h) .GT. 0.D0 ) THEN
             
             ! Equation 10 from Costa & Macedonio, 2005
             visc_heat_coeff = emme * nu_ref / ( c_p * REAL(h) ) 
             
          ELSE
             
             visc_heat_coeff = 0.D0
             
          END IF
          
          ! Viscous heating
          ! Last R.H.S. term in equation 4 from Costa & Macedonio, 2005
          expl_term(4) = expl_term(4) - visc_heat_coeff * ( REAL(u)**2          &
               + REAL(v)**2 ) * DEXP( - visc_par * ( REAL(T) - T_ref ) ) 

       END IF
          
    END IF
           
  END SUBROUTINE eval_expl_terms

END MODULE constitutive_2d

    
