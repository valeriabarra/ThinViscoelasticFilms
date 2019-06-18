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

SUBROUTINE term_c(ffx,hxxx,func_c)
  USE nrtype
  USE domain_space, ONLY: n,nmmax
  USE paras, ONLY: coef_c,dx4
  IMPLICIT NONE
  
  !!$ this subroutine computes the fourth order spatial derivative of the capillary term
  
  REAL(DP), DIMENSION(nmmax), INTENT(IN) :: ffx,hxxx
  REAL(DP), DIMENSION(nmmax), INTENT(OUT) :: func_c

  REAL(DP) :: ccx4
  INTEGER(I4B) :: k,kp
  
  ccx4=coef_c/dx4

  DO k=1,n
     kp=k+1
     func_c(k)=ccx4*(ffx(kp)*hxxx(kp)-ffx(k)*hxxx(k))
  END DO

  RETURN
END SUBROUTINE term_c
