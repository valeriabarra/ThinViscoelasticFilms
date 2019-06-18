SUBROUTINE get_matrix_x_quadterm(fquad,fderquadx0,fderquadxm,hxxx,hxxxp,hxxx0,hxxxm&
,hxxxmm,gpp_quadterm,gp_quadterm,g0_quadterm,gm_quadterm,gmm_quadterm)
  USE nrtype
  USE domain_space, ONLY: n, nmmax
  USE paras, ONLY: coef_slip,dx4
  IMPLICIT NONE
  
  !!$ this subroutine assembles the entries asscoiated to the d/dx(h^2 * h_xxx) term in the pentadiagonal matrix  

  REAL(DP), DIMENSION(nmmax), INTENT(IN) :: fquad,fderquadx0,fderquadxm,hxxx,hxxxp,hxxx0,hxxxm&
  ,hxxxmm
  REAL(DP), DIMENSION(nmmax), INTENT(OUT) :: gpp_quadterm,gp_quadterm,g0_quadterm&
  ,gm_quadterm,gmm_quadterm

  INTEGER(I4B) :: k, kp, np
  REAL(DP) :: cslipx4

  np=n+1

  cslipx4=1.0d0/dx4

  DO k=1,n
     kp=k+1
     gpp_quadterm(k)=cslipx4*(fquad(kp)*hxxxp(kp))
     gp_quadterm(k)=cslipx4*(fderquadx0(kp)*hxxx(kp)+fquad(kp)*hxxx0(kp)-fquad(k)*hxxxp(k))
     g0_quadterm(k)=cslipx4*(fderquadxm(kp)*hxxx(kp)+fquad(kp)*hxxxm(kp)&
                 -fderquadx0(k)*hxxx(k)-fquad(k)*hxxx0(k))
     gm_quadterm(k)=cslipx4*(fquad(kp)*hxxxmm(kp)-fderquadxm(k)*hxxx(k)-fquad(k)*hxxxm(k))
     gmm_quadterm(k)=cslipx4*(-fquad(k)*hxxxmm(k))
  END DO

  RETURN
END SUBROUTINE get_matrix_x_quadterm
