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

SUBROUTINE get_matrix_x_c(ffx,ffx0,ffxm,hxxx,hxxxp,hxxx0,hxxxm,hxxxmm,gpp_c,gp_c,g0_c,gm_c,gmm_c)
  USE nrtype
  USE domain_space, ONLY: n,nmmax
  USE paras, ONLY: coef_c,dx4
  IMPLICIT NONE
  
  !!$ this subroutine assembles the entries asscoiated to the capillary term d/dx(h^3 * hxxx) 
  !!$ in the pentadiagonal matrix

  REAL(DP), DIMENSION(nmmax), INTENT(IN) :: ffx,ffx0,ffxm,hxxx,hxxxp,hxxx0,hxxxm,hxxxmm
  REAL(DP), DIMENSION(nmmax), INTENT(OUT) :: gpp_c,gp_c,g0_c,gm_c,gmm_c

  INTEGER(I4B) :: k, kp, np
  REAL(DP) :: ccx

  np=n+1

  ccx=coef_c/dx4

  DO k=1,n
     kp=k+1
     gpp_c(k)=ccx*(ffx(kp)*hxxxp(kp))
     gp_c(k)=ccx*(ffx0(kp)*hxxx(kp)+ffx(kp)*hxxx0(kp)-ffx(k)*hxxxp(k))
     g0_c(k)=ccx*(ffxm(kp)*hxxx(kp)+ffx(kp)*hxxxm(kp)&
                 -ffx0(k)*hxxx(k)-ffx(k)*hxxx0(k))
     gm_c(k)=ccx*(ffx(kp)*hxxxmm(kp)-ffxm(k)*hxxx(k)-ffx(k)*hxxxm(k))
     gmm_c(k)=ccx*(-ffx(k)*hxxxmm(k))
  END DO

  RETURN
END SUBROUTINE get_matrix_x_c
