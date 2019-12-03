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


SUBROUTINE get_der_hhh_x(hxxxp,hxxx0,hxxxm,hxxxmm,hx0,hxm,hxavem,hxave0,hxavep)
  USE nrtype
  USE domain_space, ONLY: n, nmmax
  USE paras, ONLY: switch_x0, switch_xL
  IMPLICIT NONE

  !!$  form Jacobian of the derivative parts
  !!$
  !!$  hxxxp(k): derivative of (hxxx at (k-1/2)) w.r.t. res(k+1)
  !!$  hxxx0(k): derivative of (hxxx at (k-1/2)) w.r.t. res(k)
  !!$  hxxxm(k): derivative of (hxxx at (k-1/2)) w.r.t. res(k-1)
  !!$  hxxxmm(k): derivative of (hxxx at (k-1/2)) w.r.t. res(k-2)
  !!$  hx0(k): derivative of (hx at (k-1/2)) w.r.t. res(k)
  !!$  hxm(k): derivative of (hx at (k-1/2)) w.r.t. res(k-1)
  !!$  hxavem(k): derivative of (hxave at (k-1/2)) w.r.t. res(k-1)
  !!$  hxave0(k): derivative of (hxave at (k-1/2)) w.r.t. res(k)
  !!$  hxavep(k): derivative of (hxave at (k-1/2)) w.r.t. res(k+1)

  REAL(DP), DIMENSION(nmmax), INTENT(OUT) ::  hxxxp,hxxx0,hxxxm,hxxxmm,hx0,hxm,hxavem,hxave0,hxavep

  INTEGER(I4B) :: npp,np,nm,nmm,k
  npp=n+2
  np=n+1
  nm=n-1
  nmm=n-2

  hxxxp(1)=0.0d0
  hxxx0(1)=0.0d0
  hxxxm(1)=0.0d0
  hxxxmm(1)=0.0d0

  IF(switch_x0.EQ.1)THEN
     hxxxp(2)=1.0d0
     hxxx0(2)=-3.0d0
     hxxxm(2)=0.0d0
     hxxxmm(2)=0.0d0

     hxxxp(3)=1.0d0
     hxxx0(3)=-3.0d0
     hxxxm(3)=3.0d0
     hxxxmm(3)=0.0d0
  ELSEIF(switch_x0.EQ.2)THEN
     hxxxp(2)=1.0d0
     hxxx0(2)=-3.0d0
     hxxxm(2)=2.0d0
     hxxxmm(2)=0.0d0
     
     hxxxp(3)=1.0d0
     hxxx0(3)=-3.0d0
     hxxxm(3)=3.0d0
     hxxxmm(3)=-1.0d0
  END IF

  DO k=4,nmm
     hxxxmm(k)=-1.0d0
     hxxxm(k)=3.0d0
     hxxx0(k)=-3.0d0
     hxxxp(k)=1.0d0
  END DO

  IF(switch_xL.EQ.1)THEN
     hxxxp(nm)=0.0d0
     hxxx0(nm)=-3.0d0
     hxxxm(nm)=3.0d0
     hxxxmm(nm)=-1.0d0

     hxxxp(n)=0.0d0
     hxxx0(n)=0.0d0
     hxxxm(n)=3.0d0
     hxxxmm(n)=-1.0d0
  ELSEIF(switch_xL.EQ.2)THEN
     hxxxp(nm)=1.0d0
     hxxx0(nm)=-3.0d0
     hxxxm(nm)=3.0d0
     hxxxmm(nm)=-1.0d0

     hxxxp(n)=0.0d0
     hxxx0(n)=-2.0d0
     hxxxm(n)=3.0d0
     hxxxmm(n)=-1.0d0
  END IF

  hxxxp(np)=0.0d0
  hxxx0(np)=0.0d0
  hxxxm(np)=0.0d0
  hxxxmm(np)=0.0d0

!!$

  hx0(1)=0.0d0
  hxm(1)=0.0d0

  IF(switch_x0.EQ.1)THEN
     hx0(2)=1.0d0
     hxm(2)=0.0d0
  ELSEIF(switch_x0.EQ.2)THEN
     hx0(2)=1.0d0
     hxm(2)=-1.0d0
  END IF

  DO k=3,nm
     hx0(k)=1.0d0
     hxm(k)=-1.0d0
  END DO

  IF(switch_xL.EQ.1)THEN
     hx0(n)=0.0d0
     hxm(n)=-1.0d0
  ELSEIF(switch_xL.EQ.2)THEN
     hx0(n)=1.0d0
     hxm(n)=-1.0d0
  END IF

  hx0(np)=0.0d0
  hxm(np)=0.0d0

  !! This part here is needed to get the derivative of the first derivative averaged hxave.
  !! Needed to find the Jacobian L of term 2

  hxavem(1)=0.0d0
  hxave0(1)=0.0d0
  hxavep(1)=0.0d0

  IF(switch_x0.EQ.1)THEN
    hxavem(2)=0.0d0
    hxave0(2)=0.0d0
    hxavep(2)=0.5d0
  ELSEIF(switch_x0.EQ.2)THEN
    hxavem(2)=-0.5d0
    hxave0(2)=0.0d0
    hxavep(2)=0.5d0
  END IF

  DO k=3,nmm
     hxavem(k)=-0.5d0
     hxave0(k)=0.0d0
     hxavep(k)=0.5d0
  END DO

  IF(switch_xL.EQ.1)THEN
    hxavem(nm)=-0.5d0
    hxave0(nm)=0.0d0
    hxavep(nm)=0.0d0
  ELSEIF(switch_xL.EQ.2)THEN
    hxavem(nm)=-0.5d0
    hxave0(nm)=0.0d0
    hxavep(nm)=0.5d0
  END IF
  
  hxavem(n)=0.0d0
  hxave0(n)=0.0d0
  hxavep(n)=0.0d0

  RETURN
END SUBROUTINE get_der_hhh_x
