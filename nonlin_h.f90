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

SUBROUTINE nonlin_h(res,ffh)

  USE nrtype
  USE domain_space, ONLY: n, nmmax
  IMPLICIT NONE

  REAL(DP), DIMENSION(nmmax), INTENT(IN) :: res
  REAL(DP), DIMENSION(nmmax), INTENT(OUT) :: ffh

  INTEGER(I4B) np,k,km,kp

  np=n+1

  !!$  ffh(k): h at (k-1/2) (as an average , needed in term_one)

  ffh(1)=(res(1)) ! used the BC h_x=0 which implies h_1=h_0 --> ffh(1)= 2 res(1) / 2 = res(1)
  
  DO k=2,n
    km=k-1
    ffh(k)=(res(k)+res(km))/2.0d0
  END DO
  
  ffh(np)=(res(n)) ! used the BC h_x=0 which implies h_(n+1)=h_n --> ffh(n+1)= 2 res(n) / 2 = res(n)

  RETURN
END SUBROUTINE nonlin_h
