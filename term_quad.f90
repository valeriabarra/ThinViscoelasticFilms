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
