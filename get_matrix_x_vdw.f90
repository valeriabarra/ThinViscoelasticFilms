SUBROUTINE get_matrix_x_vdw(fvdwx,fvdwx0,fvdwxm,hx,hx0,hxm,gp_vdw,g0_vdw,gm_vdw)
  USE nrtype
  USE domain_space, ONLY: n, nmmax
  USE paras, ONLY: coef_vdw,dx2
  IMPLICIT NONE
  
  !!$ this subroutine assembles the entries asscoiated to the van der Waals term in the 
  !!$ pentadiagonal matrix (only three diagonals of it)

  REAL(DP), DIMENSION(nmmax), INTENT(IN) :: fvdwx,fvdwx0,fvdwxm,hx,hx0,hxm
  REAL(DP), DIMENSION(nmmax), INTENT(OUT) :: gp_vdw,g0_vdw,gm_vdw

  INTEGER(I4B) :: k, kp, np
  REAL(DP) :: ccx

  np=n+1

  ccx=coef_vdw/dx2 

  DO k=1,n
     kp=k+1
     gp_vdw(k)=ccx*(fvdwx0(kp)*hx(kp)+fvdwx(kp)*hx0(kp))
     g0_vdw(k)=ccx*(fvdwxm(kp)*hx(kp)+fvdwx(kp)*hxm(kp)-fvdwx0(k)*hx(k)-fvdwx(k)*hx0(k))
     gm_vdw(k)=ccx*(-fvdwxm(k)*hx(k)-fvdwx(k)*hxm(k))
  END DO

  RETURN
END SUBROUTINE get_matrix_x_vdw
