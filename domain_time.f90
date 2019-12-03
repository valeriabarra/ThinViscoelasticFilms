! BSD 2-Clause License
!
! Copyright (c) [2019] [Valeria Barra]
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions are met:
!
! 1. Redistributions of source code must retain the above copyright notice, this
!    list of conditions and the following disclaimer.
!
! 2. Redistributions in binary form must reproduce the above copyright notice,
!    this list of conditions and the following disclaimer in the documentation
!    and/or other materials provided with the distribution.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
! ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
! WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
! DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
! ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
! (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
! LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
! ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
! (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
! SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


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
