SUBROUTINE term_vdw(fvdwx,hx,func_vdw)
  USE nrtype
  USE domain_space, ONLY: n,nmmax
  USE paras, ONLY: coef_vdw,dx2
  IMPLICIT NONE
  
  !!$ this subroutine computes the second order spatial derivative of the van der Waals term
  
  REAL(DP), DIMENSION(nmmax), INTENT(IN) :: fvdwx,hx
  REAL(DP), DIMENSION(nmmax), INTENT(OUT) :: func_vdw

  REAL(DP) :: ccx2
  INTEGER(I4B) :: k,kp
  
  ccx2=coef_vdw/dx2

  DO k=1,n
    kp=k+1
    func_vdw(k)=ccx2*(fvdwx(kp)*hx(kp)-fvdwx(k)*hx(k))
  END DO

  RETURN
END SUBROUTINE term_vdw
