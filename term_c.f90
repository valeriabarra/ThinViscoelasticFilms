SUBROUTINE term_c(ffx,hxxx,func_c)
  USE nrtype
  USE domain_space, ONLY: n,nmmax
  USE paras, ONLY: coef_c,dx4
  IMPLICIT NONE
  
  !!$ this subroutine computes the fourth order spatial derivative of the capillary term
  
  REAL(DP), DIMENSION(nmmax), INTENT(IN) :: ffx,hxxx
  REAL(DP), DIMENSION(nmmax), INTENT(OUT) :: func_c

  REAL(DP) :: ccx4
  INTEGER(I4B) :: k,kp
  
  ccx4=coef_c/dx4

  DO k=1,n
     kp=k+1
     func_c(k)=ccx4*(ffx(kp)*hxxx(kp)-ffx(k)*hxxx(k))
  END DO

  RETURN
END SUBROUTINE term_c
