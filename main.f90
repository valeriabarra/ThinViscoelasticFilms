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


PROGRAM main
  USE nrtype
  USE domain_space, ONLY: n, nmmax
  USE domain_time
  USE paras, ONLY: fac_inc, dx
  IMPLICIT NONE

  !!$
  !!$ Thin viscoelastic film equation in 2D
  !!$

  INTERFACE
     SUBROUTINE init(res,Qi,Ri,tau_21)
        USE nrtype
        USE domain_space
        USE paras, ONLY: dx, dx2, dx3, dx4, dxh, bb, M, beta, hstar, switch_dom,init_switch, folder
        IMPLICIT NONE
       REAL(DP), DIMENSION(nmmax), INTENT(OUT) :: res,Qi,Ri
	   REAL(DP), DIMENSION(n), INTENT(OUT) ::tau_21
     END SUBROUTINE init

     SUBROUTINE output(res,t,tau_21)
         USE nrtype
         USE domain_time, ONLY: t_end
         USE domain_space, ONLY: x, n, nmmax
         USE paras, ONLY: dx, bb, folder
         IMPLICIT NONE
       REAL(DP), DIMENSION(nmmax), INTENT(IN) :: res
	   REAL(DP), DIMENSION(n), INTENT(IN) ::tau_21
       REAL(DP), INTENT(INOUT) :: t
     END SUBROUTINE output

     SUBROUTINE new_time(res,resold,dt,Qi,Ri,newton_iter,stab_iter,massold,massinit,mass_flag,tau_21)
       USE nrtype
       USE domain_space, ONLY: n, nmmax
       USE paras, ONLY: fac_dec, dt_min, percent, folder
       REAL(DP), DIMENSION(nmmax), INTENT(INOUT) :: res, resold
	   REAL(DP), DIMENSION(n), INTENT(INOUT) ::tau_21
       REAL(DP), DIMENSION(nmmax), INTENT(IN) :: Qi, Ri
       REAL(DP), INTENT(INOUT) :: dt
       REAL(DP), INTENT(INOUT) :: massold
       REAL(DP), INTENT(IN) :: massinit
       INTEGER(I4B), INTENT(OUT) :: newton_iter
       INTEGER(I4B), INTENT(INOUT) :: stab_iter
       INTEGER(I4B), INTENT(INOUT) :: mass_flag
     END SUBROUTINE new_time
  END INTERFACE

  REAL(DP) :: res(nmmax), resold(nmmax),Qi(nmmax),Ri(nmmax), tau_21(n)
  REAL(DP) :: dt, dt_prev, time_output
  REAL, DIMENSION(2) :: tarray_start, tarray_present
  REAL :: run_time_start, run_time_present

  INTEGER(I4B) :: n_time, newton_iter, stab_iter
  REAL(DP) :: massold, massinit
  INTEGER(I4B) :: mass_flag


  CALL ETIME(tarray_start, run_time_start)
  
!  get various parameters
  
  CALL par
  
!  get initial conditions
  
  call init(res,Qi,Ri,tau_21)
  resold(1:n)=res(1:n)
  massold=SUM(res(1:n))*dx
  massinit=massold

  
  WRITE(6,*)'----------'
  WRITE(6,*)'start evolution'
  WRITE(6,*)'----------'
  WRITE(6,99)'n_time','t','dt_did','n_iter','time_1','time_2'

99 FORMAT(A9,2(A11),A9,2(A11))
  
  t=t0
  dt=dt_init
  dt_prev=dt

  CALL output(res,t,tau_21)
  time_output=t+t_out

!!$  We want to output the status every 100 steps
!!$  n_time is the counting on steps

  n_time=0
  stab_iter=0

  DO WHILE (t.LT.t_end)

    CALL new_time(res,resold,dt,Qi,Ri,newton_iter,stab_iter,massold,massinit,mass_flag,tau_21)

    IF (mass_flag.EQ.1) THEN
        CYCLE
    ENDIF

    t=t+dt
    n_time=n_time+1
     
      IF(t.GE.(time_output-1.0d-10))THEN
        CALL output(res,t,tau_21)

        !!   Control if the interface has become steady here !!

     !   IF(MAXVAL(abs(res-resold))<=tolsteady)THEN
     !        WRITE(6,*)'The solution is steady!'
     !        CALL flush(6)
     !       STOP
     !   END IF

        IF(t.GE.(t_crit-1.0d-10))THEN
            IF(dt.GE.t_out_small)THEN
                time_output=time_output+dt
             ELSE
                time_output=time_output+t_out_small ! this is to have more outputs after the critical time
             END IF
        ELSE
           time_output=time_output+t_out
        END IF

     END IF
     
     IF(MOD(n_time,100).EQ.0)THEN
        CALL ETIME(tarray_present, run_time_present)
        WRITE(6,100)n_time, t, dt, newton_iter, run_time_present-run_time_start, &
             tarray_present(2)-tarray_start(2)
        CALL FLUSH(6)
     END IF

100  FORMAT(I9,2(ES11.2),I9,2(ES11.2))

     IF(stab_iter.GE.10)THEN
        dt=dt*fac_inc
        stab_iter=0
     END IF
     IF(dt.GT.DABS(t-t_end)) dt=t_end-t

  END DO
  
! DONE
  
  STOP
END PROGRAM main
