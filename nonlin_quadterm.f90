SUBROUTINE nonlin_quadterm(res,fquad)
  USE nrtype
  USE domain_space, ONLY: n, nmmax
  USE mobility, ONLY: quadterm
  IMPLICIT NONE

  !!$c	get quadratic part for capilary term and vdW
  !!$
  !!$  fquad(k,l): fquad at (k-1/2, l)
  !!$

  REAL(DP), DIMENSION(nmmax), INTENT(IN) :: res
  REAL(DP), DIMENSION(nmmax), INTENT(OUT) :: fquad

  INTEGER(I4B) np,k,km

  np=n+1

  !!$  fquad(k,l): fquad at (k-1/2, l)

  fquad(1)=quadterm(res(1))
  
  DO k=2,n
    km=k-1
    fquad(k)=(quadterm(res(k))+quadterm(res(km)))/2.0d0
  END DO
  
  fquad(np)=quadterm(res(n))


  RETURN
END SUBROUTINE nonlin_quadterm
