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


SUBROUTINE newton(res,Qi, Ri,Qinew,Rinew,resold,dt,newton_iter,iflag_newton,tau_21)
  USE nrtype
  USE nr_ban
  USE domain_time
  USE domain_space, ONLY: n, nmmax
  USE paras
  USE vanderwaals
  USE mobility
  IMPLICIT NONE

  !!$ This subroutines assembles all the terms that are needed in the linearized system (Aq = b) to solve for the correction term

  INTERFACE

    SUBROUTINE hhh_cap(res,hxxx,hx,hxxxave,hxave)
        USE nrtype
        USE domain_space, ONLY: n, nmmax
        IMPLICIT NONE
        REAL(DP), DIMENSION(nmmax), INTENT(IN) :: res
        REAL(DP), DIMENSION(nmmax), INTENT(OUT) :: hxxx,hx,hxxxave,hxave
    END SUBROUTINE hhh_cap

    SUBROUTINE nonlin_h(res,ffh)
        USE nrtype
        USE domain_space, ONLY: n, nmmax
        IMPLICIT NONE
        REAL(DP), DIMENSION(nmmax), INTENT(IN) :: res
        REAL(DP), DIMENSION(nmmax), INTENT(OUT) :: ffh
    END SUBROUTINE nonlin_h


    SUBROUTINE nonlin_func2(res,func,func_four,func_six,ffx,hxxx,hx,hxave,hxxxave,fvdwx,fquadvdwx,fquad,fvdwxsimple)
        USE nrtype
        USE domain_space, ONLY: n, nmmax
        USE paras, ONLY: bb
        IMPLICIT NONE
        REAL(DP), DIMENSION(nmmax), INTENT(IN) :: res
        REAL(DP), DIMENSION(nmmax), INTENT(OUT) :: func, ffx, hxxx, hx,hxave,hxxxave,fvdwx, fquadvdwx, fquad,func_four&
       ,func_six,fvdwxsimple
    END SUBROUTINE nonlin_func2

    SUBROUTINE get_der_cap_x(res,ffx0,ffxm)
        USE nrtype
        USE domain_space, ONLY: n, nmmax
        USE mobility, ONLY: dmob
        USE paras, ONLY: switch_x0, switch_xL
        IMPLICIT NONE
        REAL(DP), DIMENSION(nmmax), INTENT(IN) :: res
        REAL(DP), DIMENSION(nmmax), INTENT(OUT) :: ffx0,ffxm
    END SUBROUTINE get_der_cap_x

    SUBROUTINE get_der_quadterm(res,fderquadx0,fderquadxm)
        USE nrtype
        USE domain_space, ONLY: n, nmmax
        USE mobility, ONLY: dquadterm, dmob
        USE paras, ONLY: switch_x0, switch_xL
        IMPLICIT NONE
        REAL(DP), DIMENSION(nmmax), INTENT(IN) :: res
        REAL(DP), DIMENSION(nmmax), INTENT(OUT) :: fderquadx0,fderquadxm
    END SUBROUTINE get_der_quadterm

    SUBROUTINE get_der_vdw_x(res,fvdwx0,fvdwxm)
        USE nrtype
        USE domain_space, ONLY: n, nmmax
        USE vanderwaals, ONLY:  dvdw
        USE paras, ONLY: switch_x0, switch_xL
        REAL(DP), DIMENSION(nmmax), INTENT(IN) :: res
        REAL(DP), DIMENSION(nmmax), INTENT(OUT) :: fvdwx0,fvdwxm
    END SUBROUTINE get_der_vdw_x

    SUBROUTINE get_der_quadvdw_x(res,fquadvdwx0,fquadvdwxm)
        USE nrtype
        USE domain_space, ONLY: n, nmmax
        USE vanderwaals, ONLY: dquadvdw
        USE paras, ONLY: switch_x0, switch_xL
        REAL(DP), DIMENSION(nmmax), INTENT(IN) :: res
        REAL(DP), DIMENSION(nmmax), INTENT(OUT) :: fquadvdwx0,fquadvdwxm
    END SUBROUTINE get_der_quadvdw_x

    SUBROUTINE get_der_hhh_x(hxxxp,hxxx0,hxxxm,hxxxmm,hx0,hxm,hxavem,hxave0,hxavep)
        USE nrtype
        USE domain_space, ONLY: n, nmmax
        USE paras, ONLY: switch_x0, switch_xL
        REAL(DP), DIMENSION(nmmax), INTENT(OUT) ::  hxxxp,hxxx0,hxxxm,hxxxmm,hx0,hxm,hxavem,hxave0,hxavep
    END SUBROUTINE get_der_hhh_x

    SUBROUTINE get_matrix_x2(dt,ma_x,ffx,hxxx,hx,ffx0,ffxm,hxxxp&
              ,hxxx0,hxxxm,hxxxmm,hx0,hxm,hxavem,hxave0,hxavep,fvdwx,fvdwx0&
              ,fvdwxm, fquadvdwx,fquadvdwx0,fquadvdwxm,fquad,fderquadx0,fderquadxm&
              ,coeff_one,coeff_three)
        USE nrtype
        USE domain_time
        USE domain_space, ONLY: n, nmmax
        USE paras, ONLY: theta,l1,bb,l2
        IMPLICIT NONE
        REAL(DP), INTENT(IN) :: dt
        REAL(DP), DIMENSION(nmmax,1:5), INTENT(OUT) :: ma_x
        REAL(DP), DIMENSION(nmmax), INTENT(IN) :: ffx,ffx0,ffxm,hxxxp,hxxx0,hxxxm,hxxxmm,hx0,hxm&
       ,fvdwx,fvdwx0,fvdwxm,fquadvdwx,fquadvdwx0,fquadvdwxm,fquad,fderquadx0,fderquadxm&
       ,hxxx,hx,coeff_one,coeff_three,hxavem,hxave0,hxavep
    END SUBROUTINE get_matrix_x2

    SUBROUTINE get_matrix_x_hxave(hxavem,hxave0,hxavep,gm_hxave,g0_hxave,gp_hxave)
        USE nrtype
        USE domain_space, ONLY: n, nmmax
        USE paras, ONLY: dx
        REAL(DP), DIMENSION(nmmax), INTENT(IN) :: hxavem,hxave0,hxavep
        REAL(DP), DIMENSION(nmmax), INTENT(OUT) :: gm_hxave,g0_hxave,gp_hxave
    END SUBROUTINE get_matrix_x_hxave

    SUBROUTINE ExplicitQ(Qi,dt,hx,hxxx,Qinew,fvdwxsimple)
        USE nrtype
        USE domain_space, ONLY: n,nmmax
        USE paras, ONLY: coef_vdw,dx,dx3,l2
        IMPLICIT NONE
        REAL(DP), DIMENSION(nmmax), INTENT(IN) :: Qi,hx,hxxx,fvdwxsimple
        REAL(DP), DIMENSION(nmmax), INTENT(INOUT) :: Qinew
        REAL(DP), INTENT(IN) :: dt
    END SUBROUTINE ExplicitQ

    SUBROUTINE ExplicitR(Ri,dt,ffh,hx,hxxx,Rinew,fvdwxsimple)
        USE nrtype
        USE domain_space, ONLY: n,nmmax
        USE paras, ONLY: coef_vdw,dx,dx3,l2
        IMPLICIT NONE
        REAL(DP), DIMENSION(nmmax), INTENT(IN) :: Ri,ffh,hx,hxxx,fvdwxsimple
        REAL(DP), DIMENSION(nmmax), INTENT(INOUT) :: Rinew
        REAL(DP), INTENT(IN) :: dt
    END SUBROUTINE ExplicitR

    SUBROUTINE  term_coeff_one(coeff_one,Qinew,Rinew,fquad,ffh)
        USE nrtype
        USE domain_space, ONLY: n, nmmax
        USE paras, ONLY: dx,l1,l2
        IMPLICIT NONE
        REAL(DP), DIMENSION(nmmax), INTENT(IN) ::fquad,ffh,Qinew,Rinew
        REAL(DP), DIMENSION(nmmax), INTENT(OUT) :: coeff_one
    END SUBROUTINE term_coeff_one

    SUBROUTINE  term_coeff_three(coeff_three,Qinew,Rinew,fquad,ffh)
        USE nrtype
        USE domain_space, ONLY: n, nmmax
        USE paras, ONLY:l1, l2
        IMPLICIT NONE
        REAL(DP), DIMENSION(nmmax), INTENT(IN) :: Qinew,Rinew,fquad,ffh
        REAL(DP), DIMENSION(nmmax), INTENT(OUT) :: coeff_three
    END SUBROUTINE term_coeff_three

  END INTERFACE

  REAL(DP), DIMENSION(nmmax), INTENT(INOUT) :: Qinew,Rinew
  REAL(DP), DIMENSION(nmmax), INTENT(INOUT) :: res
  REAL(DP), DIMENSION(n), INTENT(INOUT) :: tau_21

  REAL(DP), DIMENSION(nmmax), INTENT(IN) :: Qi, Ri,resold
  REAL(DP), INTENT(IN) :: dt
  INTEGER(I4B), INTENT(OUT) :: newton_iter, iflag_newton

  REAL(DP) :: ffx(nmmax),ffh(nmmax)
  REAL(DP) :: hxxx(nmmax),hx(nmmax),hxave(nmmax),hxxxave(nmmax)
  REAL(DP) :: res_guess(nmmax),func(nmmax),func_guess(nmmax),b_fix(nmmax),b(nmmax)
  REAL(DP) :: func_four(nmmax),func_six(nmmax),func_four_guess(nmmax),func_six_guess(nmmax),hxave_guess(nmmax)

  REAL(DP) :: err_newton,err_posit,err_resid,err_local,err_trunc,small
  REAL(DP) :: ffx0(nmmax),ffxm(nmmax)

  REAL(DP) :: fquad(nmmax)
  REAL(DP) :: fderquadx0(nmmax),fderquadxm(nmmax)
  REAL(DP) :: fvdwx0(nmmax),fvdwxm(nmmax)
  REAL(DP) :: fvdwx(nmmax),fvdwxsimple(nmmax)
  REAL(DP) :: fquadvdwx(nmmax)
  REAL(DP) :: fquadvdwx0(nmmax), fquadvdwxm(nmmax)
  
  REAL(DP) :: hxxxp(nmmax),hxxx0(nmmax),hxxxm(nmmax),hxxxmm(nmmax),hx0(nmmax),hxm(nmmax)
  REAL(DP) :: hxavem(nmmax),hxave0(nmmax),hxavep(nmmax) ! needed as Jacobian matrix L (hx_ave, centered at the grid point k)
  REAL(DP) :: coeff_one(nmmax),coeff_three(nmmax)!needed as output of term_coeff_one, term_coeff_three

  REAL(DP) :: ma_x(nmmax,1:5)

  REAL(DP) :: flux(n)

  INTEGER(I4B) :: np, k, iflag_penta

  np=n+1

  small=1.0d-20

  newton_iter=0
  iflag_newton=1
  err_newton=1.0d0

  CALL nonlin_h(res,ffh)

  CALL nonlin_func2(res,func,func_four,func_six,ffx,hxxx,hx,hxave,hxxxave,fvdwx,fquadvdwx,fquad,fvdwxsimple)

  CALL ExplicitQ(Qi,dt,hx,hxxx,Qinew,fvdwxsimple)
  CALL ExplicitR(Ri,dt,ffh,hx,hxxx,Rinew,fvdwxsimple)

  CALL term_coeff_one(coeff_one,Qinew,Rinew,fquad,ffh)

  CALL term_coeff_three(coeff_three,Qinew,Rinew,fquad,ffh)

  res_guess(1:n)=res(1:n)
  func_guess(1:n)=func(1:n)
  func_four_guess(1:n)=func_four(1:n)
  func_six_guess(1:n)=func_six(1:n)
  hxave_guess(1:n)=hxave(1:n)

  IF(l2.EQ.0.0d0)THEN ! this IF puts back the code in the original form when l2=0, first order in time equation
    DO k=1,n
        b_fix(k)=coeff_one(k)*res(k)-theta*dt*func(k) +coeff_three(k)*hxave(k) + l1*func_four(k)
    END DO

  ELSE ! this is the version of the code for l2 NOT zero, second order in time equation

    IF(t.EQ.t0)THEN ! This creates a slightly different RHS bfix only for the first time step, n=0, accounting for the initial h^(-1)=h^1 (so there is no resold)

      DO k=1,n
        b_fix(k)=(2.0d0*l2 + dt*coeff_one(k))*res(k)-theta*(dt**2)*func(k) +dt*coeff_three(k)*hxave(k) + l1*dt*func_four(k)&
        +l2*dt*func_six(k)
      END DO

    ELSE
      DO k=1,n
        b_fix(k)=(2.0d0*l2 + dt*coeff_one(k))*res(k)-theta*(dt**2)*func(k) +dt*coeff_three(k)*hxave(k) + l1*dt*func_four(k)&
        +l2*dt*func_six(k) -l2*resold(k) !the only difference in b_fix after the first step is that there is the -l2*resold, before not
      END DO
    END IF

  END IF

  !!$  Start Newton's iteration
  DO WHILE (err_newton.GT.eps_newt)

    IF(l2.EQ.0.0d0)THEN ! this IF puts back the code in the original form when l2=0, first order in DeltaT
      DO k=1,n
        b(k)=-l1*func_four_guess(k)-coeff_one(k)*res_guess(k)-coeff_three(k)*hxave_guess(k)-(1.0d0-theta)*dt*func_guess(k)+b_fix(k)
      END DO

    ELSE ! this is the version of the code for l2 NOT zero, second order in DeltaT

      IF(t.EQ.t0)THEN  ! This creates a slightly different RHS b only for the first time step, n=0, accounting for the initial h^(-1)=h^1 (so there is no resold)
        DO k=1,n
            b(k)=-l1*dt*func_four_guess(k) -l2*dt*func_six_guess(k) - (2.0d0*l2 + dt*coeff_one(k))*res_guess(k)&
            -dt*coeff_three(k)*hxave_guess(k)-(1.0d0-theta)*(dt**2)*func_guess(k) + b_fix(k)
            ! the difference in b for the first time step is that there is - (2.0d0*l2 + dt*coeff_one(k))*res_guess(k), while otherwise
            ! it would be just - (l2 + dt*coeff_one(k))*res_guess(k), without 2 multiplying l2
        END DO
      ELSE
        DO k=1,n
            b(k)=-l1*dt*func_four_guess(k) -l2*dt*func_six_guess(k) - (l2 + dt*coeff_one(k))*res_guess(k)&
            -dt*coeff_three(k)*hxave_guess(k)-(1.0d0-theta)*(dt**2)*func_guess(k) + b_fix(k)
        END DO
      END IF
    END IF


  !!$  Solve for the correction
  CALL get_der_cap_x(res_guess,ffx0,ffxm)

  CALL get_der_quadterm(res_guess,fderquadx0,fderquadxm)

  CALL get_der_vdw_x(res_guess,fvdwx0,fvdwxm)

  CALL get_der_quadvdw_x(res_guess,fquadvdwx0,fquadvdwxm)

  CALL get_der_hhh_x(hxxxp,hxxx0,hxxxm,hxxxmm,hx0,hxm,hxavem,hxave0,hxavep)

  CALL get_matrix_x2(dt,ma_x,ffx,hxxx,hx,ffx0,ffxm,hxxxp&
                    ,hxxx0,hxxxm,hxxxmm,hx0,hxm,hxavem,hxave0,hxavep,fvdwx,fvdwx0&
                    ,fvdwxm, fquadvdwx,fquadvdwx0,fquadvdwxm,fquad,fderquadx0,fderquadxm&
                    ,coeff_one,coeff_three)

  IF(switch_x0.eq.1)THEN
    ma_x(1,3)=1.0d0
    ma_x(1,4)=0.0d0
    ma_x(1,5)=0.0d0
    b(1)=0.0d0
  END IF
     
  IF(switch_xL.eq.1)THEN
    ma_x(n,1)=0.0d0
    ma_x(n,2)=0.0d0
    ma_x(n,3)=1.0d0
    b(n)=0.0d0
  END IF

  ! CALL of the pentadiagonal system solver. Note: the result is stored in b!!!
  CALL penta(ma_x,b,n,iflag_penta)

  err_newton=0.0d0
  err_posit=1.0d2

  DO k=1,n
  !!$  Check if Newton's method converged
    IF(res_guess(k).LT.(1.0d-15))THEN
        err_resid=DABS(b(k)/(1.0d-15))
    ELSE
        err_resid=DABS(b(k)/res_guess(k))
    END IF
    IF(err_resid.GT.err_newton)THEN
        err_newton=err_resid
    END IF

  !!$  Update the guess and check for positivity
    res_guess(k)=res_guess(k)+b(k) ! Since the result of penta was stored in b, we update the sol adding b
    
    IF(res_guess(k).LT.err_posit)THEN
        err_posit=res_guess(k)
    END IF
        
  END DO

  IF(err_posit.LT.small)THEN
    WRITE(6,*)'Negative solution, dt=', dt
    WRITE(6,*)'            err_posit=', err_posit
    CALL FLUSH(6)

    OPEN(9,FILE=folder//'resNeg.dat')
    
    DO k=1,n
        WRITE(9,101)res_guess(k) ! This prints out the last solution at the last time regardless if Newton converged
    END DO
    CALL FLUSH(9)
    STOP
    RETURN
  END IF

  !!$ Update the RHS functions for the truncation error check
  CALL nonlin_func2(res_guess,func_guess,func_four_guess,func_six_guess,ffx,hxxx,hx,hxave,hxxxave,fvdwx,fquadvdwx&
                   ,fquad,fvdwxsimple)

    newton_iter=newton_iter+1

    IF(newton_iter.GT.n_iter_max)THEN
      EXIT
    END IF

  END DO !!! END OF NEWTON's ITERATIVE LOOP

  !!$  Evaluate truncation error
  err_trunc=0.0d0
  DO k=1,n
    err_local=DABS(res_guess(k)-res(k))+dt*(DABS(func_guess(k)-func(k)))! + DABS(g_guess(k)-g(k)) + DABS(hxave_guess(k)-hxave(k)))
    IF(err_local.GT.err_trunc)THEN
       err_trunc=err_local
    END IF
  END DO

  IF(err_newton.LE.eps_newt)THEN
    IF(err_trunc.LE.tol)THEN
       iflag_newton=0
       res(1:n)=res_guess(1:n)
    ELSE
       WRITE(6,98)'Too large local truncation error, dt=', dt, 'err_trunc=', err_trunc
    END IF

  ELSE
    WRITE(6,99)'Newton does not converge, dt=', dt, 'err_newton=', err_newton
  END IF

  IF(t.GT.t0)THEN
    DO k=1,n
      flux(k) = coef_c*mob(res(k))*(1.0d0/dx3)*hxxxave(k) - coef_vdw*vdw(res(k))*(hxave(k)/dx) - bb*quadterm(res(k))*(1.0d0/dx3)&
      *hxxxave(k)-bb*(coef_vdw*3.0d0)*quadvdw(res(k))*(hxave(k)/dx) - (  l1*(quadterm(res(k))/2.0d0)*(1.0d0/dx3)*hxxxave(k)&
      +l1*(coef_vdw*3.0d0/2.0d0)*quadvdw(res(k))*(hxave(k)/dx)  )*((res(k)-resold(k))/dt)
    END DO
  
    !!$ this computes the shear stress tensor component 12 for sanity check 
    DO k=1,n
	  tau_21(k) = (1.0d0/2.0d0)*(res(k))*( (1.0d0/dx3)*hxxxave(k) + (coef_vdw/dx)*vdwsimple(res(k))*(hxave(k)) )
    END DO
  
  END IF

98  FORMAT(A37,ES11.2,A15,ES11.2)
99  FORMAT(A28,ES11.2,A15,ES11.2)
101 FORMAT(ES25.15)


  CALL FLUSH(6)

  RETURN
END SUBROUTINE newton
