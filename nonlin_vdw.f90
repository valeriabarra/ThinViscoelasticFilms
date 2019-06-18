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

SUBROUTINE nonlin_vdw(res,fvdwx,fvdwxsimple)

  USE nrtype
  USE domain_space, ONLY: n, nmmax
  USE vanderwaals, ONLY: vdw, vdwsimple
  USE mobility, ONLY: mob
  IMPLICIT NONE

  !!$c    get nonlinear part for van der waals terms
  !!$
  !!$  fvdwx(k): fvdwx at (k-1/2) (the one originally in the code, already multiplied by h^3)
  !!$  fvdwxsimple(k): fvdwxsimple at (k-1/2) (needed in Q and R)

  REAL(DP), DIMENSION(nmmax), INTENT(IN) :: res
  REAL(DP), DIMENSION(nmmax), INTENT(OUT) :: fvdwx,fvdwxsimple

  INTEGER(I4B) np,k,km

  np=n+1

  fvdwx(1)=vdw(res(1))
  fvdwxsimple(1)=vdwsimple(res(1)) ! needed in Q and R
  
  DO k=2,n
    km=k-1
    fvdwx(k)=(vdw(res(k))+vdw(res(km)))/2.0d0
    fvdwxsimple(k)=(vdwsimple(res(k)) + vdwsimple(res(km)))/2.0d0 ! needed in Q and R
  END DO
  
  fvdwx(np)=vdw(res(n))
  fvdwxsimple(np)=vdwsimple(res(n)) ! needed in Q and R

  RETURN
END SUBROUTINE nonlin_vdw
