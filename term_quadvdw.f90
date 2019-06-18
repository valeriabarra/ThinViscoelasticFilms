SUBROUTINE term_quadvdw(fquadvdwx,hx,func_quadvdw)
  USE nrtype
  USE domain_space, ONLY: n,nmmax
  USE paras, ONLY: coef_slip,coef_vdw,dx2
  IMPLICIT NONE
  
  !!$ this subroutine computes the second order spatial derivative of the quadratic term times van der Waals term
  
  REAL(DP), DIMENSION(nmmax), INTENT(IN) :: fquadvdwx,hx
  REAL(DP), DIMENSION(nmmax), INTENT(OUT) :: func_quadvdw

  REAL(DP) :: cslipx2
  INTEGER(I4B) :: k,kp
  
  cslipx2=(coef_vdw*3.0d0)/dx2 ! I have to multiply by three because in the term h^2 vdw there is no 1/3 that instead is in the (h^3/3)*vdw, so it was included in coeff. vdw

  DO k=1,n
    kp=k+1
    func_quadvdw(k)=cslipx2*(fquadvdwx(kp)*hx(kp)-fquadvdwx(k)*hx(k))
  END DO

  RETURN
END SUBROUTINE term_quadvdw
