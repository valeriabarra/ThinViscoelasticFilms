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


SUBROUTINE term_coeff_three(coeff_three,QinewExplicit,RinewExplicit,fquad,ffh)
  USE nrtype
  USE domain_space, ONLY: n, nmmax
  USE paras, ONLY:l1, l2
  IMPLICIT NONE

  !!$ this subroutine computes the "mhat" term in the reference paper
  !!$ "Interfacial dynamics of thin viscoelastic films and drops", by Barra, Afkhami, Kondic, JNNFM (2016)
  !!$ (called circled 3 in personal notes)

  REAL(DP), DIMENSION(nmmax), INTENT(IN) :: QinewExplicit,RinewExplicit,fquad,ffh
  REAL(DP), DIMENSION(nmmax), INTENT(OUT) :: coeff_three

  INTEGER(I4B) :: k,kp

  ! this computes the average of the argument in parenthesis of coeff_one, b/c in coeff_three there is no outside derivative d/dx, so I need the grid point 'k'

  DO k=1,n
    kp=k+1
    coeff_three(k)= (l2 - l1)*( (DBLE(1.0d0/2.0d0))*fquad(kp)*QinewExplicit(kp)  - ffh(kp)*RinewExplicit(kp) &
    +(DBLE(1.0d0/2.0d0))*fquad(k)*QinewExplicit(k) - ffh(k)*RinewExplicit(k)  )/2.0d0
  END DO

END SUBROUTINE term_coeff_three
