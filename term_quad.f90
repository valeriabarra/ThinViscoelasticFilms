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


SUBROUTINE term_quad(fquad,hxxx,func_quad)
  USE nrtype
  USE domain_space, ONLY: n,nmmax
  USE paras, ONLY: coef_slip,dx4
  IMPLICIT NONE
  
  !!$ this subroutine computes the fourth order spatial derivative of the quadratic term
  
  REAL(DP), DIMENSION(nmmax), INTENT(IN) :: fquad,hxxx
  REAL(DP), DIMENSION(nmmax), INTENT(OUT) :: func_quad

  REAL(DP) :: cslipx4
  INTEGER(I4B) :: k,kp
  
  cslipx4=1.0d0/dx4

  DO k=1,n
    kp=k+1
    func_quad(k)=cslipx4*(fquad(kp)*hxxx(kp)-fquad(k)*hxxx(k))
  END DO

  RETURN
END SUBROUTINE term_quad
