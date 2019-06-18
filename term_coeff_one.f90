SUBROUTINE term_coeff_one(coeff_one,QinewExplicit,RinewExplicit,fquad,ffh)
  USE nrtype
  USE domain_space, ONLY: n, nmmax
  USE paras, ONLY: dx,l1,l2
  IMPLICIT NONE
   
  !!$ this subroutine computes the first order spatial derivative of the "ghat" term in the reference paper 
  !!$ "Interfacial dynamics of thin viscoelastic films and drops", by Barra, Afkhami, Kondic, JNNFM (2016)
  !!$ (called circled 1 in personal notes)

  REAL(DP), DIMENSION(nmmax), INTENT(IN) :: fquad,ffh,QinewExplicit,RinewExplicit
  REAL(DP), DIMENSION(nmmax), INTENT(OUT) :: coeff_one

  INTEGER(I4B) :: k,kp
  REAL(DP) :: ccx

  ccx=1.0d0/dx

  DO k=1,n
    kp=k+1
    coeff_one(k)= 1.0d0 +(l2 - l1)*ccx*( (DBLE(1.0d0/2.0d0))*fquad(kp)*QinewExplicit(kp)  - ffh(kp)*RinewExplicit(kp) &
    -(DBLE(1.0d0/2.0d0))*fquad(k)*QinewExplicit(k) + ffh(k)*RinewExplicit(k)  )
  END DO

END SUBROUTINE term_coeff_one
