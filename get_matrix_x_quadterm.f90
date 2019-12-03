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
