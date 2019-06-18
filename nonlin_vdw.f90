SUBROUTINE nonlin_vdw(res,fvdwx,fvdwxsimple)

  USE nrtype
  USE domain_space, ONLY: n, nmmax
  USE vanderwaals, ONLY: vdw, vdwsimple
  USE mobility, ONLY: mob
  IMPLICIT NONE

  !!$c    get nonlinear part for van der waals terms
  !!$
  !!$  fvdwx(k): fvdwx at (k-1/2) (the one originally in the code, already multiplied by h^3)
  !!$  fvdwxsimple(k): fvdwxsimple at (k-1/2) (needed in Q and R)

  REAL(DP), DIMENSION(nmmax), INTENT(IN) :: res
  REAL(DP), DIMENSION(nmmax), INTENT(OUT) :: fvdwx,fvdwxsimple

  INTEGER(I4B) np,k,km

  np=n+1

  fvdwx(1)=vdw(res(1))
  fvdwxsimple(1)=vdwsimple(res(1)) ! needed in Q and R
  
  DO k=2,n
    km=k-1
    fvdwx(k)=(vdw(res(k))+vdw(res(km)))/2.0d0
    fvdwxsimple(k)=(vdwsimple(res(k)) + vdwsimple(res(km)))/2.0d0 ! needed in Q and R
  END DO
  
  fvdwx(np)=vdw(res(n))
  fvdwxsimple(np)=vdwsimple(res(n)) ! needed in Q and R

  RETURN
END SUBROUTINE nonlin_vdw
