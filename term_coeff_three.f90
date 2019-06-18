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
