SUBROUTINE nonlin_cap(res,ffx)
  USE nrtype
  USE domain_space, ONLY: n, nmmax
  USE mobility, ONLY: mob
  IMPLICIT NONE

  !!$c	get nonlinear part for capilary term
  !!$
  !!$  ffx(k,l): ffx at (k-1/2, l)
  !!$

  REAL(DP), DIMENSION(nmmax), INTENT(IN) :: res
  REAL(DP), DIMENSION(nmmax), INTENT(OUT) :: ffx

  INTEGER(I4B) np,k,km,kp

  np=n+1

  !!$  ffx(k,l): ffx at (k-1/2, l)
  !!$  ffnemx(k,l): ffnemx at (k-1/2, l)

  ffx(1)=mob(res(1))
  DO k=2,n
     km=k-1
     ffx(k)=(mob(res(k))+mob(res(km)))/2.0d0
  END DO
  ffx(np)=mob(res(n))

  RETURN
END SUBROUTINE nonlin_cap
