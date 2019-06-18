! BSD 2-Clause (MIT) License
!
! Copyright (c) [2019] [ThinViscoelasticFilms]
!
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:
!
! The above copyright notice and this permission notice shall be included in all
! copies or substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
! SOFTWARE.

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
