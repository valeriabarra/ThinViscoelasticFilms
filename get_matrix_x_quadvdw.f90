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


SUBROUTINE get_matrix_x_quadvdw(fquadvdwx,fquadvdwx0,fquadvdwxm,hx,hx0,hxm,gp_quadvdw,g0_quadvdw,gm_quadvdw)
  USE nrtype
  USE domain_space, ONLY: n, nmmax
  USE paras, ONLY: coef_slip,coef_vdw,dx2
  IMPLICIT NONE
  
  !!$ this subroutine assembles the entries asscoiated to the quadratic term times the van der Waals
  !!$ in the pentadiagonal matrix (only three diagonals of it)

  REAL(DP), DIMENSION(nmmax), INTENT(IN) :: fquadvdwx,fquadvdwx0,fquadvdwxm,hx,hx0,hxm
  REAL(DP), DIMENSION(nmmax), INTENT(OUT) :: gp_quadvdw,g0_quadvdw,gm_quadvdw

  INTEGER(I4B) :: k, kp, np
  REAL(DP) :: ccx

  np=n+1

  ccx=(coef_vdw*3.0d0)/dx2  ! I have to multiply by three because in the term (h^2 vdw) there is no 1/3 that instead is in the (h^3/3)*vdw, so it was included in coeff. vdw

  DO k=1,n
     kp=k+1
     gp_quadvdw(k)=ccx*(fquadvdwx0(kp)*hx(kp)+fquadvdwx(kp)*hx0(kp))
     g0_quadvdw(k)=ccx*(fquadvdwxm(kp)*hx(kp)+fquadvdwx(kp)*hxm(kp)-fquadvdwx0(k)*hx(k)-fquadvdwx(k)*hx0(k))
     gm_quadvdw(k)=ccx*(-fquadvdwxm(k)*hx(k)-fquadvdwx(k)*hxm(k))
  END DO

  RETURN
END SUBROUTINE get_matrix_x_quadvdw
