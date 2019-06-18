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

SUBROUTINE get_der_cap_x(res,ffx0,ffxm)
  USE nrtype
  USE domain_space, ONLY: n, nmmax
  USE mobility, ONLY: dmob
  USE paras, ONLY: switch_x0, switch_xL
  IMPLICIT NONE

  !!$  form Jacobian of the nonlinear part for capilary term
  !!$
  !!$  ffx0(k): derivative of (ffx at (k-1/2)) w.r.t. res(k)
  !!$  ffxm(k): derivative of (ffx at (k-1/2)) w.r.t. res(k-1)

  REAL(DP), DIMENSION(nmmax), INTENT(IN) :: res
  REAL(DP), DIMENSION(nmmax), INTENT(OUT) :: ffx0,ffxm

  INTEGER(I4B) np,nm,k,km

  np=n+1
  nm=n-1

  IF(switch_x0.EQ.1)THEN
     ffx0(1)=0.0d0
     ffxm(1)=0.0d0
     
     ffx0(2)=dmob(res(2))/2.0d0
     ffxm(2)=0.0d0
     
  ELSEIF(switch_x0.EQ.2)THEN
     ffx0(1)=dmob(res(1))
     ffxm(1)=0.0d0
     
     ffx0(2)=dmob(res(2))/2.0d0
     ffxm(2)=dmob(res(1))/2.0d0
     

  END IF
  
  DO k=3,nm
     km=k-1
     ffx0(k)=dmob(res(k))/2.0d0
     ffxm(k)=dmob(res(km))/2.0d0
  END DO
  
  IF(switch_xL.EQ.1)THEN
     ffx0(n)=0.0d0
     ffxm(n)=dmob(res(nm))/2.0d0
     
     ffx0(np)=0.0d0
     ffxm(np)=0.0d0
    
  ELSEIF(switch_xL.EQ.2)THEN
     ffx0(n)=dmob(res(n))/2.0d0
     ffxm(n)=dmob(res(nm))/2.0d0
     
     ffx0(np)=0.0d0
     ffxm(np)=dmob(res(n))

     
  END IF

  RETURN
END SUBROUTINE get_der_cap_x
