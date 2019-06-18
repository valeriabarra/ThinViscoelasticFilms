SUBROUTINE ExplicitQ(Qi,dt,hx,hxxx,QinewExplicit,fvdwxsimple)
  USE nrtype
  USE domain_space, ONLY: n,nmmax
  USE paras, ONLY: coef_vdw,dx,dx3,l2
  IMPLICIT NONE
  REAL(DP), DIMENSION(nmmax), INTENT(IN) :: Qi,hx,hxxx,fvdwxsimple
  REAL(DP), DIMENSION(nmmax), INTENT(INOUT) :: QinewExplicit
  REAL(DP), INTENT(IN) :: dt
  INTEGER(I4B) :: k

  !!$ this routine, solves the algebraic equation for Q (in case l2=0), 
  !!$ or implements the explicit Euler scheme (in case l2 is not zero)

  IF(l2.EQ.0.0d0)THEN
    DO k=1,nmmax
      QinewExplicit(k)= -(1.0d0/dx3)*hxxx(k) - ((coef_vdw*3.0d0)/dx)*fvdwxsimple(k)*hx(k)
    END DO
  
  ELSE

    DO k=1,nmmax
      QinewExplicit(k)= l2*Qi(k) -dt*Qi(k) - dt*((1.0d0/dx3)*hxxx(k)&
      + ((coef_vdw*3.0d0)/dx)*fvdwxsimple(k)*hx(k))
    END DO
    
    QinewExplicit(1:nmmax) = l2*QinewExplicit(1:nmmax)

  END IF

END SUBROUTINE ExplicitQ


SUBROUTINE ExplicitR(Ri,dt,ffh,hx,hxxx,RinewExplicit,fvdwxsimple)
  USE nrtype
  USE domain_space, ONLY: n,nmmax
  USE paras, ONLY: coef_vdw,dx,dx3,l2
  IMPLICIT NONE
  REAL(DP), DIMENSION(nmmax), INTENT(IN) :: Ri,ffh,hx,hxxx,fvdwxsimple
  REAL(DP), DIMENSION(nmmax), INTENT(INOUT) :: RinewExplicit
  REAL(DP), INTENT(IN) :: dt
  INTEGER(I4B) :: k

  !!$ this routine, solves the algebraic equation for R (in case l2=0), 
  !!$ or implements the explicit Euler scheme (in case l2 is not zero)

  IF(l2.EQ.0.0d0)THEN
    DO k=1,nmmax
      RinewExplicit(k)= -(ffh(k))*((1.0d0/dx3)*hxxx(k) + ((coef_vdw*3.0d0)/dx)*fvdwxsimple(k)*hx(k))
    END DO
  
  ELSE

    DO k=1,nmmax
      RinewExplicit(k)= l2*Ri(k) -dt*Ri(k) - dt*(ffh(k))*((1.0d0/dx3)*hxxx(k)&
      + ((coef_vdw*3.0d0)/dx)*fvdwxsimple(k)*hx(k))
    END DO

    RinewExplicit(1:nmmax) = l2*RinewExplicit(1:nmmax)
  END IF

END SUBROUTINE ExplicitR
