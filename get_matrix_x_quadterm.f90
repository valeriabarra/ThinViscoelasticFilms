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

SUBROUTINE get_matrix_x_quadterm(fquad,fderquadx0,fderquadxm,hxxx,hxxxp,hxxx0,hxxxm&
                                ,hxxxmm,gpp_quadterm,gp_quadterm,g0_quadterm,gm_quadterm,gmm_quadterm)
  USE nrtype
  USE domain_space, ONLY: n, nmmax
  USE paras, ONLY: coef_slip,dx4
  IMPLICIT NONE
  
  !!$ this subroutine assembles the entries asscoiated to the d/dx(h^2 * h_xxx) term in the pentadiagonal matrix  

  REAL(DP), DIMENSION(nmmax), INTENT(IN) :: fquad,fderquadx0,fderquadxm,hxxx,hxxxp,hxxx0,hxxxm&
  ,hxxxmm
  REAL(DP), DIMENSION(nmmax), INTENT(OUT) :: gpp_quadterm,gp_quadterm,g0_quadterm&
  ,gm_quadterm,gmm_quadterm

  INTEGER(I4B) :: k, kp, np
  REAL(DP) :: cslipx4

  np=n+1

  cslipx4=1.0d0/dx4

  DO k=1,n
     kp=k+1
     gpp_quadterm(k)=cslipx4*(fquad(kp)*hxxxp(kp))
     gp_quadterm(k)=cslipx4*(fderquadx0(kp)*hxxx(kp)+fquad(kp)*hxxx0(kp)-fquad(k)*hxxxp(k))
     g0_quadterm(k)=cslipx4*(fderquadxm(kp)*hxxx(kp)+fquad(kp)*hxxxm(kp)&
                 -fderquadx0(k)*hxxx(k)-fquad(k)*hxxx0(k))
     gm_quadterm(k)=cslipx4*(fquad(kp)*hxxxmm(kp)-fderquadxm(k)*hxxx(k)-fquad(k)*hxxxm(k))
     gmm_quadterm(k)=cslipx4*(-fquad(k)*hxxxmm(k))
  END DO

  RETURN
END SUBROUTINE get_matrix_x_quadterm
