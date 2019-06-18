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
