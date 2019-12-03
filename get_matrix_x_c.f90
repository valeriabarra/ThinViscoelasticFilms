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
