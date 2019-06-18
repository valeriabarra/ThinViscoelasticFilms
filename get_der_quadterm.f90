SUBROUTINE get_der_quadterm(res,fderquadx0,fderquadxm)
  USE nrtype
  USE domain_space, ONLY: n, nmmax
  USE mobility, ONLY: dquadterm, dmob
  USE paras, ONLY: switch_x0, switch_xL
  IMPLICIT NONE

  !!$  form Jacobian of the nonlinear part for quadratic term times the third order derivative
  !!$
  !!$  fderquadx0(k): derivative of (fquad at (k-1/2)) w.r.t. res(k)
  !!$  fderquadxm(k): derivative of (fquad at (k-1/2)) w.r.t. res(k-1)

  REAL(DP), DIMENSION(nmmax), INTENT(IN) :: res
  REAL(DP), DIMENSION(nmmax), INTENT(OUT) :: fderquadx0,fderquadxm

  INTEGER(I4B) np,nm,k,km

  np=n+1
  nm=n-1

  IF(switch_x0.EQ.1)THEN
     fderquadx0(1)=dquadterm(res(1))
     fderquadxm(1)=0.0d0
     
     fderquadx0(2)=dquadterm(res(2))/2.0d0
     fderquadxm(2)=dquadterm(res(1))/2.0d0
     
  ELSEIF(switch_x0.EQ.2)THEN
     fderquadx0(1)=dquadterm(res(1))
     fderquadxm(1)=0.0d0
     
     fderquadx0(2)=dquadterm(res(2))/2.0d0
     fderquadxm(2)=dquadterm(res(1))/2.0d0
     
  END IF
  
  DO k=3,nm
     km=k-1
     fderquadx0(k)=dquadterm(res(k))/2.0d0
     fderquadxm(k)=dquadterm(res(km))/2.0d0
  END DO
  
  IF(switch_xL.EQ.1)THEN
     fderquadx0(n)=dquadterm(res(n))/2.0d0
     fderquadxm(n)=dquadterm(res(nm))/2.0d0
     
     fderquadx0(np)=0.0d0
     fderquadxm(np)=dquadterm(res(n))
    
  ELSEIF(switch_xL.EQ.2)THEN
     fderquadx0(n)=dquadterm(res(n))/2.0d0
     fderquadxm(n)=dquadterm(res(nm))/2.0d0
     
     fderquadx0(np)=0.0d0
     fderquadxm(np)=dquadterm(res(n))
     
  END IF

  RETURN
END SUBROUTINE get_der_quadterm
