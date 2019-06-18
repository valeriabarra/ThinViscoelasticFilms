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

SUBROUTINE get_matrix_x_vdw(fvdwx,fvdwx0,fvdwxm,hx,hx0,hxm,gp_vdw,g0_vdw,gm_vdw)
  USE nrtype
  USE domain_space, ONLY: n, nmmax
  USE paras, ONLY: coef_vdw,dx2
  IMPLICIT NONE
  
  !!$ this subroutine assembles the entries asscoiated to the van der Waals term in the 
  !!$ pentadiagonal matrix (only three diagonals of it)

  REAL(DP), DIMENSION(nmmax), INTENT(IN) :: fvdwx,fvdwx0,fvdwxm,hx,hx0,hxm
  REAL(DP), DIMENSION(nmmax), INTENT(OUT) :: gp_vdw,g0_vdw,gm_vdw

  INTEGER(I4B) :: k, kp, np
  REAL(DP) :: ccx

  np=n+1

  ccx=coef_vdw/dx2 

  DO k=1,n
     kp=k+1
     gp_vdw(k)=ccx*(fvdwx0(kp)*hx(kp)+fvdwx(kp)*hx0(kp))
     g0_vdw(k)=ccx*(fvdwxm(kp)*hx(kp)+fvdwx(kp)*hxm(kp)-fvdwx0(k)*hx(k)-fvdwx(k)*hx0(k))
     gm_vdw(k)=ccx*(-fvdwxm(k)*hx(k)-fvdwx(k)*hxm(k))
  END DO

  RETURN
END SUBROUTINE get_matrix_x_vdw
