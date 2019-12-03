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


SUBROUTINE get_der_cap_x(res,ffx0,ffxm)
  USE nrtype
  USE domain_space, ONLY: n, nmmax
  USE mobility, ONLY: dmob
  USE paras, ONLY: switch_x0, switch_xL
  IMPLICIT NONE

  !!$  form Jacobian of the nonlinear part for capilary term
  !!$
  !!$  ffx0(k): derivative of (ffx at (k-1/2)) w.r.t. res(k)
  !!$  ffxm(k): derivative of (ffx at (k-1/2)) w.r.t. res(k-1)

  REAL(DP), DIMENSION(nmmax), INTENT(IN) :: res
  REAL(DP), DIMENSION(nmmax), INTENT(OUT) :: ffx0,ffxm

  INTEGER(I4B) np,nm,k,km

  np=n+1
  nm=n-1

  IF(switch_x0.EQ.1)THEN
     ffx0(1)=0.0d0
     ffxm(1)=0.0d0
     
     ffx0(2)=dmob(res(2))/2.0d0
     ffxm(2)=0.0d0
     
  ELSEIF(switch_x0.EQ.2)THEN
     ffx0(1)=dmob(res(1))
     ffxm(1)=0.0d0
     
     ffx0(2)=dmob(res(2))/2.0d0
     ffxm(2)=dmob(res(1))/2.0d0
     

  END IF
  
  DO k=3,nm
     km=k-1
     ffx0(k)=dmob(res(k))/2.0d0
     ffxm(k)=dmob(res(km))/2.0d0
  END DO
  
  IF(switch_xL.EQ.1)THEN
     ffx0(n)=0.0d0
     ffxm(n)=dmob(res(nm))/2.0d0
     
     ffx0(np)=0.0d0
     ffxm(np)=0.0d0
    
  ELSEIF(switch_xL.EQ.2)THEN
     ffx0(n)=dmob(res(n))/2.0d0
     ffxm(n)=dmob(res(nm))/2.0d0
     
     ffx0(np)=0.0d0
     ffxm(np)=dmob(res(n))

     
  END IF

  RETURN
END SUBROUTINE get_der_cap_x
