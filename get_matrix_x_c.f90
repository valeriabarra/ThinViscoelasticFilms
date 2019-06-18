SUBROUTINE get_matrix_x_c(ffx,ffx0,ffxm,hxxx,hxxxp,hxxx0,hxxxm,hxxxmm,gpp_c,gp_c,g0_c,gm_c,gmm_c)
  USE nrtype
  USE domain_space, ONLY: n,nmmax
  USE paras, ONLY: coef_c,dx4
  IMPLICIT NONE
  
  !!$ this subroutine assembles the entries asscoiated to the capillary term d/dx(h^3 * hxxx) 
  !!$ in the pentadiagonal matrix

  REAL(DP), DIMENSION(nmmax), INTENT(IN) :: ffx,ffx0,ffxm,hxxx,hxxxp,hxxx0,hxxxm,hxxxmm
  REAL(DP), DIMENSION(nmmax), INTENT(OUT) :: gpp_c,gp_c,g0_c,gm_c,gmm_c

  INTEGER(I4B) :: k, kp, np
  REAL(DP) :: ccx

  np=n+1

  ccx=coef_c/dx4

  DO k=1,n
     kp=k+1
     gpp_c(k)=ccx*(ffx(kp)*hxxxp(kp))
     gp_c(k)=ccx*(ffx0(kp)*hxxx(kp)+ffx(kp)*hxxx0(kp)-ffx(k)*hxxxp(k))
     g0_c(k)=ccx*(ffxm(kp)*hxxx(kp)+ffx(kp)*hxxxm(kp)&
                 -ffx0(k)*hxxx(k)-ffx(k)*hxxx0(k))
     gm_c(k)=ccx*(ffx(kp)*hxxxmm(kp)-ffxm(k)*hxxx(k)-ffx(k)*hxxxm(k))
     gmm_c(k)=ccx*(-ffx(k)*hxxxmm(k))
  END DO

  RETURN
END SUBROUTINE get_matrix_x_c
