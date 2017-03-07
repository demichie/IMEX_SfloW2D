!********************************************************************************
!> \mainpage   IMEX_SfloW2D - Shallow Water Finite volume solver
!> IMEX_SfloW2D is a FORTRAN90 code designed to solve an hyperbolic 
!> system of partial differential equations with relaxation and source
!> terms.\n 
!> The model is discretized in time with an explicit-implicit Runge-Kutta
!> method where the hyperbolic part is solved explicetely and the other
!> terms (relaxation and surce) are treated implicitely.\n
!> The finite volume solver for the hyperbolic part of the system is based
!> on a semidiscrete central scheme and it is not tied on the specific 
!> eigenstructure of the model.\n
!> The implicit part is solved with a Newton-Raphson method where the 
!> elements of the Jacobian of the nonlinear system are evaluated 
!> numerically with a complex step derivative technique.\n
!> Version 1.0:\n
!> \n 
!> Github project page: http://demichie.github.io/IMEX_SfloW2D/
!> \n
!> \authors Mattia de' Michieli Vitturi (*), Giacomo Lari (**)\n
!> (*) Istituto Nazionale di Geofisica e vulcanologia, sezione di Pisa\n
!>     Via della Faggiola, 36\n
!>     I-56126 Pisa, Italy \n
!>     E-mail: mattia.demichielivitturi@ingv.it \n
!> (**) Dipartimento di Matematica, Università di Pisa\n
!>       E-mail: lari@student.dm.unipi.it
!********************************************************************************

!> \brief Main Program 
PROGRAM IMEX_SfloW_2d

  USE constitutive_2d, ONLY : init_problem_param

  USE geometry_2d, ONLY : init_grid

  USE geometry_2d, ONLY : dx,dy,B_cent

  USE init_2d, ONLY : riemann_problem, initial_conditions

  USE inpout_2d, ONLY : init_param
  USE inpout_2d, ONLY : read_param
  USE inpout_2d, ONLY : output_solution
  USE inpout_2d, ONLY : read_solution
  USE inpout_2d, ONLY : close_units
  USE inpout_2d, ONLY : output_dakota

  USE solver_2d, ONLY : allocate_solver_variables
  USE solver_2d, ONLY : deallocate_solver_variables
  USE solver_2d, ONLY : imex_RK_solver
  USE solver_2d, ONLY : timestep, timestep2

  USE inpout_2d, ONLY : restart

  USE parameters_2d, ONLY : t_start
  USE parameters_2d, ONLY : t_end
  USE parameters_2d, ONLY : t_output
  USE parameters_2d, ONLY : riemann_flag
  USE parameters_2d, ONLY : verbose_level

  USE solver_2d, ONLY : q , dt

  IMPLICIT NONE

  REAL*8 :: t
  REAL*8 :: t1 , t2

  CALL cpu_time(t1)

  CALL init_param

  CALL read_param

  CALL init_grid
  
  CALL init_problem_param
  
  CALL allocate_solver_variables
     
  IF ( restart ) THEN

     CALL read_solution

  ELSE

     IF( riemann_flag .EQV. .TRUE. )THEN

        ! riemann problem defined in file.inp
        CALL riemann_problem

     ELSE

        ! generic problem defined by initial conditions function (in init_2d.f90)
        CALL initial_conditions

     ENDIF

  END IF

  t = t_start
  
  WRITE(*,*) 
  WRITE(*,*) '******** START COMPUTATION *********'
  WRITE(*,*)

  IF ( verbose_level .GE. 1 ) THEN
     
     WRITE(*,*) 'Min q(1,:,:)=',MINVAL(q(1,:,:))
     WRITE(*,*) 'Max q(1,:,:)=',MAXVAL(q(1,:,:))

     WRITE(*,*) 'Min B(:,:)=',MINVAL(B_cent(:,:))
     WRITE(*,*) 'Max B(:,:)=',MAXVAL(B_cent(:,:))

     WRITE(*,*) 'size B_cent',size(B_cent,1),size(B_cent,2)

     WRITE(*,*) 'SUM(q(1,:,:)=',SUM(q(1,:,:))
     WRITE(*,*) 'SUM(B_cent(:,:)=',SUM(B_cent(:,:))
     
  END IF
  

  WRITE(*,*) 't =',t,' dt =',dt, ' h_tot =',dx*dy*(SUM(q(1,:,:)-B_cent(:,:)))

  CALL output_solution(t)

  DO WHILE ( t .LT. t_end )

     CALL timestep
     !CALL timestep2

     IF ( t+dt .GT. t_end ) dt = t_end - t
     IF ( t+dt .GT. t_output ) dt = t_output - t

     CALL imex_RK_solver

     t = t+dt

     WRITE(*,*) 't =',t,' dt =',dt, ' h_tot =',dx*dy*(SUM(q(1,:,:)-B_cent(:,:)))

     IF ( ( t .GE. t_output ) .OR. ( t .GE. t_end ) ) CALL output_solution(t)

  END DO

  ! CALL output_dakota

  CALL deallocate_solver_variables

  CALL close_units

  CALL cpu_time(t2)

  WRITE(*,*) 'Time taken by the code was',t2-t1,'seconds'

END PROGRAM IMEX_Sflow_2d

