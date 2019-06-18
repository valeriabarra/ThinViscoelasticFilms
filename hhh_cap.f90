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
