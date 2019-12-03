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


SUBROUTINE hhh_cap(res,hxxx,hx,hxxxave,hxave)

  USE nrtype
  USE domain_space, ONLY: n, nmmax
  IMPLICIT NONE

  !!$  get derivatives in x direction of the capillary term
  !!$
  !!$  hxxx(k): hxxx at (k-1/2)
  !!$  hx(k): hx at (k-1/2)
  !!$

  REAL(DP), DIMENSION(nmmax), INTENT(IN) :: res
  REAL(DP), DIMENSION(nmmax), INTENT(OUT) :: hxxx,hx,hxave,hxxxave

  INTEGER(I4B) :: np,nm,nmm,nmmm
  INTEGER(I4B) :: k,kp,km,kmm

  np=n+1
  nm=n-1
  nmm=n-2
  nmmm=n-3

  !!$  hxxx(k): hxxx at (k-1/2)

  hxxx(1)=0.0d0
  hxxx(2)=res(3)-3.0d0*res(2)+2.0d0*res(1)
  DO k=3,nm
     kp=k+1
     km=k-1
     kmm=k-2
     hxxx(k)=res(kp)-3.0d0*res(k)+3.0d0*res(km)-res(kmm)
  END DO
  hxxx(n)=-res(nmm)+3.0d0*res(nm)-2.0d0*res(n)
  hxxx(np)=0.0d0

  !!$  hx(k): hx at (k-1/2)

  hx(1)=0.0d0
  DO k=2,n
     km=k-1
     hx(k)=res(k)-res(km)
  END DO
  hx(np)=0.0d0

  !!$ hxxxave(k) (centered at k) average of hxxx(k)<-- (k - 1/2 actually) and hxxx(k+1) (<-- k + 1/2 actually)
  DO k=1,n
     kp=k+1
     hxxxave(k)=(hxxx(k)+hxxx(kp))/2.0d0
  END DO
  !!$ hxave(k) (centered at k) average of hx(k)<-- (k - 1/2 actually) and hx(k+1) (<-- k + 1/2 actually)

  DO k=1,n
     kp=k+1
     hxave(k)=(hx(k)+hx(kp))/2.0d0
  END DO



  RETURN
END SUBROUTINE hhh_cap
