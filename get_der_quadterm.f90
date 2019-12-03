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


SUBROUTINE get_der_quadterm(res,fderquadx0,fderquadxm)
  USE nrtype
  USE domain_space, ONLY: n, nmmax
  USE mobility, ONLY: dquadterm, dmob
  USE paras, ONLY: switch_x0, switch_xL
  IMPLICIT NONE

  !!$  form Jacobian of the nonlinear part for quadratic term times the third order derivative
  !!$
  !!$  fderquadx0(k): derivative of (fquad at (k-1/2)) w.r.t. res(k)
  !!$  fderquadxm(k): derivative of (fquad at (k-1/2)) w.r.t. res(k-1)

  REAL(DP), DIMENSION(nmmax), INTENT(IN) :: res
  REAL(DP), DIMENSION(nmmax), INTENT(OUT) :: fderquadx0,fderquadxm

  INTEGER(I4B) np,nm,k,km

  np=n+1
  nm=n-1

  IF(switch_x0.EQ.1)THEN
     fderquadx0(1)=dquadterm(res(1))
     fderquadxm(1)=0.0d0
     
     fderquadx0(2)=dquadterm(res(2))/2.0d0
     fderquadxm(2)=dquadterm(res(1))/2.0d0
     
  ELSEIF(switch_x0.EQ.2)THEN
     fderquadx0(1)=dquadterm(res(1))
     fderquadxm(1)=0.0d0
     
     fderquadx0(2)=dquadterm(res(2))/2.0d0
     fderquadxm(2)=dquadterm(res(1))/2.0d0
     
  END IF
  
  DO k=3,nm
     km=k-1
     fderquadx0(k)=dquadterm(res(k))/2.0d0
     fderquadxm(k)=dquadterm(res(km))/2.0d0
  END DO
  
  IF(switch_xL.EQ.1)THEN
     fderquadx0(n)=dquadterm(res(n))/2.0d0
     fderquadxm(n)=dquadterm(res(nm))/2.0d0
     
     fderquadx0(np)=0.0d0
     fderquadxm(np)=dquadterm(res(n))
    
  ELSEIF(switch_xL.EQ.2)THEN
     fderquadx0(n)=dquadterm(res(n))/2.0d0
     fderquadxm(n)=dquadterm(res(nm))/2.0d0
     
     fderquadx0(np)=0.0d0
     fderquadxm(np)=dquadterm(res(n))
     
  END IF

  RETURN
END SUBROUTINE get_der_quadterm
