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

SUBROUTINE nonlin_quadterm(res,fquad)
  USE nrtype
  USE domain_space, ONLY: n, nmmax
  USE mobility, ONLY: quadterm
  IMPLICIT NONE

  !!$c	get quadratic part for capilary term and vdW
  !!$
  !!$  fquad(k,l): fquad at (k-1/2, l)
  !!$

  REAL(DP), DIMENSION(nmmax), INTENT(IN) :: res
  REAL(DP), DIMENSION(nmmax), INTENT(OUT) :: fquad

  INTEGER(I4B) np,k,km

  np=n+1

  !!$  fquad(k,l): fquad at (k-1/2, l)

  fquad(1)=quadterm(res(1))
  
  DO k=2,n
    km=k-1
    fquad(k)=(quadterm(res(k))+quadterm(res(km)))/2.0d0
  END DO
  
  fquad(np)=quadterm(res(n))


  RETURN
END SUBROUTINE nonlin_quadterm
