SUBROUTINE get_matrix_x_hxave(hxavem,hxave0,hxavep,gm_hxave,g0_hxave,gp_hxave)
  USE nrtype
  USE domain_space, ONLY: n, nmmax
  USE paras, ONLY: dx
  IMPLICIT NONE

  !!$ this subroutine assembles the entries asscoiated to the hxave term in the 
  !!$ pentadiagonal matrix (only three diagonals of it)

  REAL(DP), DIMENSION(nmmax), INTENT(IN) :: hxavem,hxave0,hxavep
  REAL(DP), DIMENSION(nmmax), INTENT(OUT) :: gm_hxave,g0_hxave,gp_hxave

  INTEGER(I4B) :: k, kp, np
  REAL(DP) :: ccx

  ccx=1.0d0/dx

  DO k=1,n
     gm_hxave(k)=ccx*(hxavem(k))
     g0_hxave(k)=ccx*(hxave0(k))
     gp_hxave(k)=ccx*(hxavep(k))
  END DO

  RETURN
END SUBROUTINE get_matrix_x_hxave
