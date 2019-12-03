! BSD 2-Clause License
!
! Copyright (c) [2019] [Valeria Barra]
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions are met:
!
! 1. Redistributions of source code must retain the above copyright notice, this
!    list of conditions and the following disclaimer.
!
! 2. Redistributions in binary form must reproduce the above copyright notice,
!    this list of conditions and the following disclaimer in the documentation
!    and/or other materials provided with the distribution.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
! ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
! WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
! DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
! ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
! (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
! LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
! ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
! (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
! SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


SUBROUTINE init(res,Qi,Ri,tau_21)
  USE nrtype
  USE domain_space
  USE paras
  IMPLICIT NONE
  
  !!$ this subroutine sets the Initial Condition of h, Q, and R

  REAL(DP), DIMENSION(nmmax), INTENT(OUT) :: res,Qi,Ri
  REAL(DP), DIMENSION(n), INTENT(OUT) ::tau_21

  REAL(DP) :: stepx, h0, Kmax, ls, y0
  INTEGER(I4B) :: k, l

  REAL(DP) :: dom(nmmax)
  
  ! specify grid, init configuration, etc

  ! the grid construction : shifted grid
  stepx=(xmax-x0)/DBLE(n)
  DO k=1,n
     x(k)=x0+stepx*(DBLE(k)-0.5d0)
  END DO

  ! initial reference flat thickness
  h0=1.0d0
 
  ! initial condition depends if it is a new compuation (default init_switch=0) or if it continues an old one (init_switch =1)
  select case (init_switch)
 
  case (0) ! this is the default cosine IC
  IF(switch_dom.EQ.0) THEN
     Kmax= DBLE(2.0d0*PI/xmax)
     WRITE(6,98)'K_init=',Kmax
  ELSEIF(switch_dom.EQ.1)THEN
     Kmax= DBLE(PI/xmax)
     WRITE(6,98)'K_init=',Kmax
  END IF
  WRITE(6,98)'h0=',h0
  WRITE(6,98)'epsilon=',E

  DO k=1,n
        res(k)=h0+E*(COS(Kmax*x(k))) ! for both whole and half domain now since I took care of K already
  END DO

  !  this case is the case of the continuation
  case (1)
  OPEN(11,file=folder//'resLastOut.dat',STATUS='OLD',ACTION='READWRITE')
  DO k=1,n      
      READ(11,101)dom(k), res(k) !  we read the data from file 11 and put it into res
  END DO
  CLOSE(11)

  case(2) !  this is the case of a rivlet/droplet
  
  IF(switch_dom.EQ.0) THEN
     Kmax= DBLE(2.0d0*PI/xmax)
     WRITE(6,98)'K_init=',Kmax
  ELSEIF(switch_dom.EQ.1)THEN
     Kmax= DBLE(PI/xmax)
     WRITE(6,98)'K_init=',Kmax
  END IF
  WRITE(6,98)'h0=',h0
  WRITE(6,98)'epsilon=',E
  WRITE(6,98)'radius=',r_phys
  WRITE(6,98)'init_theta=',theta_init
  
   !  careful about making drop dimensionless:
    y0=-r_phys*COS(theta_init)
    ls=SQRT(r_phys**2-(char_height*hstar-y0)**2)
    DO k=1,n
      IF (x(k).LE.ls/char_length) THEN
  !  careful about making drop dimensionless:
        res(k)=(y0+SQRT(r_phys**2 - (char_length*x(k))**2))/char_height
      ELSE 
        res(k)=hstar
      END IF
    END DO



98 FORMAT(A15,F12.5) 

  end select

!  Initial Conditions for Q and R =0
Qi(1:nmmax)=0.0d0
Ri(1:nmmax)=0.0d0

!  Initial Conditions for tau_21 = 0
tau_21(1:n)=0.0d0


CALL flush(10)
CALL flush(6)

101 FORMAT(2(ES24.16)) 

!  Dirichlet boundary conditions
!  Note : it's fixed for all t for switch=1

  bdx0=res(1)
  bdxL=res(n)

  dx=stepx
  dx2=dx**2
  dx3=dx**3
  dx4=dx**4
  dxh=dx/2.0d0

  RETURN
END SUBROUTINE init
