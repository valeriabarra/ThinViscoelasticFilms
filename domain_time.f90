MODULE domain_time
  USE nrtype
  IMPLICIT NONE
  SAVE

  !!$  Time domain: [t0, tend] with initial time step, dt_init
  !!$  Have output profile every t_out time unit

  REAL(DP) :: t ! variable t needs to be global so to be used in the check of first iterations in Newton and get_matrix_x_2

  REAL(DP), PARAMETER :: t0=0.0d0 ! initial time
  REAL(DP), PARAMETER :: t_end=1.0d1 ! final time
  REAL(DP), PARAMETER :: dt_init=1.d-1 ! initial dt
  REAL(DP), PARAMETER :: t_out=5.0d0 ! time at which output solution
  REAL(DP), PARAMETER :: t_crit= 9.5d1 ! for debugging purposes, the idea of a critical time was to have finer outputs only after a specific time where solution was varying much
  REAL(DP), PARAMETER :: t_out_small=1.0d-2 ! at the critical time, the output frequency could have been different

END MODULE domain_time
