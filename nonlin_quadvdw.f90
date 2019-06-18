SUBROUTINE nonlin_quadvdw(res,fquadvdwx)

  USE nrtype
  USE domain_space, ONLY: n, nmmax
  USE vanderwaals, ONLY: quadvdw
  IMPLICIT NONE

  !!$c	get nonlinear part for van der waals term multiplied by quadratic term
  !!$
  !!$  fquadvdwx(k,l): fquadvdwx at (k-1/2, l)
  !!$

  REAL(DP), DIMENSION(nmmax), INTENT(IN) :: res
  REAL(DP), DIMENSION(nmmax), INTENT(OUT) :: fquadvdwx

  INTEGER(I4B) np,k,km

  np=n+1

  !!$  fquadvdwx(k,l): fquadvdwx at (k-1/2, l)

  fquadvdwx(1)=quadvdw(res(1))
  
  DO k=2,n
    km=k-1
    fquadvdwx(k)=(quadvdw(res(k))+quadvdw(res(km)))/2.0d0
  END DO
  
  fquadvdwx(np)=quadvdw(res(n))


  RETURN
END SUBROUTINE nonlin_quadvdw

