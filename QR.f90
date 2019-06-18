! BSD 2-Clause (MIT) License
!
! Copyright (c) [2019] [ThinViscoelasticFilms]
!
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:
!
! The above copyright notice and this permission notice shall be included in all
! copies or substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
! SOFTWARE.

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
