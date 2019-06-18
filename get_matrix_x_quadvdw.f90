SUBROUTINE get_matrix_x_quadvdw(fquadvdwx,fquadvdwx0,fquadvdwxm,hx,hx0,hxm,gp_quadvdw,g0_quadvdw,gm_quadvdw)
  USE nrtype
  USE domain_space, ONLY: n, nmmax
  USE paras, ONLY: coef_slip,coef_vdw,dx2
  IMPLICIT NONE
  
  !!$ this subroutine assembles the entries asscoiated to the quadratic term times the van der Waals
  !!$ in the pentadiagonal matrix (only three diagonals of it)

  REAL(DP), DIMENSION(nmmax), INTENT(IN) :: fquadvdwx,fquadvdwx0,fquadvdwxm,hx,hx0,hxm
  REAL(DP), DIMENSION(nmmax), INTENT(OUT) :: gp_quadvdw,g0_quadvdw,gm_quadvdw

  INTEGER(I4B) :: k, kp, np
  REAL(DP) :: ccx

  np=n+1

  ccx=(coef_vdw*3.0d0)/dx2  ! I have to multiply by three because in the term (h^2 vdw) there is no 1/3 that instead is in the (h^3/3)*vdw, so it was included in coeff. vdw

  DO k=1,n
     kp=k+1
     gp_quadvdw(k)=ccx*(fquadvdwx0(kp)*hx(kp)+fquadvdwx(kp)*hx0(kp))
     g0_quadvdw(k)=ccx*(fquadvdwxm(kp)*hx(kp)+fquadvdwx(kp)*hxm(kp)-fquadvdwx0(k)*hx(k)-fquadvdwx(k)*hx0(k))
     gm_quadvdw(k)=ccx*(-fquadvdwxm(k)*hx(k)-fquadvdwx(k)*hxm(k))
  END DO

  RETURN
END SUBROUTINE get_matrix_x_quadvdw
