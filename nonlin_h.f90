SUBROUTINE nonlin_h(res,ffh)

  USE nrtype
  USE domain_space, ONLY: n, nmmax
  IMPLICIT NONE

  REAL(DP), DIMENSION(nmmax), INTENT(IN) :: res
  REAL(DP), DIMENSION(nmmax), INTENT(OUT) :: ffh

  INTEGER(I4B) np,k,km,kp

  np=n+1

  !!$  ffh(k): h at (k-1/2) (as an average , needed in term_one)

  ffh(1)=(res(1)) ! used the BC h_x=0 which implies h_1=h_0 --> ffh(1)= 2 res(1) / 2 = res(1)
  
  DO k=2,n
    km=k-1
    ffh(k)=(res(k)+res(km))/2.0d0
  END DO
  
  ffh(np)=(res(n)) ! used the BC h_x=0 which implies h_(n+1)=h_n --> ffh(n+1)= 2 res(n) / 2 = res(n)

  RETURN
END SUBROUTINE nonlin_h
