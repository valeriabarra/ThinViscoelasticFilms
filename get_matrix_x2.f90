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

SUBROUTINE get_matrix_x2(dt,ma_x,ffx,hxxx,hx,ffx0,ffxm,hxxxp&
                        ,hxxx0,hxxxm,hxxxmm,hx0,hxm,hxavem,hxave0,hxavep,fvdwx,fvdwx0&
                        ,fvdwxm, fquadvdwx,fquadvdwx0,fquadvdwxm,fquad,fderquadx0,fderquadxm&
                        ,coeff_one,coeff_three)
  USE nrtype
  USE domain_time
  USE domain_space, ONLY: n, nmmax
  USE paras, ONLY: theta,l1,bb,l2
  IMPLICIT NONE
  
  !!$ this subroutine assembles all the different contributions in the pentadiagonal matrix

  INTERFACE
     SUBROUTINE get_matrix_x_c(ffx,ffx0,ffxm,hxxx,hxxxp,hxxx0,hxxxm,hxxxmm,gpp_c,gp_c,g0_c,gm_c,gmm_c)
        USE nrtype
        USE domain_space, ONLY: n,nmmax
        USE paras, ONLY: coef_c,dx4
        IMPLICIT NONE
          REAL(DP), DIMENSION(nmmax), INTENT(IN) :: ffx,ffx0,ffxm,hxxx,hxxxp,hxxx0,hxxxm,hxxxmm
          REAL(DP), DIMENSION(nmmax), INTENT(OUT) :: gpp_c,gp_c,g0_c,gm_c,gmm_c
     END SUBROUTINE get_matrix_x_c

     SUBROUTINE get_matrix_x_quadterm(fquad,fderquadx0,fderquadxm,hxxx,hxxxp&
     ,hxxx0,hxxxm,hxxxmm,gpp_quadterm,gp_quadterm,g0_quadterm&
     ,gm_quadterm,gmm_quadterm)
        USE nrtype
        USE domain_space, ONLY: n, nmmax
        USE paras, ONLY: coef_slip,dx4
        IMPLICIT NONE
          REAL(DP), DIMENSION(nmmax), INTENT(IN) :: fquad,fderquadx0,fderquadxm,hxxx,hxxxp,hxxx0&
          ,hxxxm,hxxxmm
          REAL(DP), DIMENSION(nmmax), INTENT(OUT) :: gpp_quadterm,gp_quadterm&
          ,g0_quadterm,gm_quadterm,gmm_quadterm
     END SUBROUTINE get_matrix_x_quadterm

     SUBROUTINE get_matrix_x_vdw(fvdwx,fvdwx0,fvdwxm,hx,hx0,hxm,gp_vdw,g0_vdw,gm_vdw)
        USE nrtype
        USE domain_space, ONLY: n, nmmax
        USE paras, ONLY: coef_vdw,dx2
        IMPLICIT NONE
          REAL(DP), DIMENSION(nmmax), INTENT(IN) :: fvdwx,fvdwx0,fvdwxm,hx,hx0,hxm
          REAL(DP), DIMENSION(nmmax), INTENT(OUT) :: gp_vdw,g0_vdw,gm_vdw
     END SUBROUTINE get_matrix_x_vdw

     SUBROUTINE get_matrix_x_quadvdw(fquadvdwx,fquadvdwx0,fquadvdwxm,hx,hx0,hxm,gp_quadvdw&
     ,g0_quadvdw,gm_quadvdw)
        USE nrtype
        USE domain_space, ONLY: n, nmmax
        USE paras, ONLY: coef_slip,coef_vdw,dx2
        IMPLICIT NONE
          REAL(DP), DIMENSION(nmmax), INTENT(IN) :: fquadvdwx,fquadvdwx0,fquadvdwxm,hx,hx0,hxm
          REAL(DP), DIMENSION(nmmax), INTENT(OUT) :: gp_quadvdw,g0_quadvdw,gm_quadvdw
     END SUBROUTINE get_matrix_x_quadvdw

     SUBROUTINE get_matrix_x_hxave(hxavem,hxave0,hxavep,gm_hxave,g0_hxave,gp_hxave)
       USE nrtype
       USE domain_space, ONLY: n, nmmax
       USE paras, ONLY: dx
       IMPLICIT NONE
         REAL(DP), DIMENSION(nmmax), INTENT(IN) :: hxavem,hxave0,hxavep
         REAL(DP), DIMENSION(nmmax), INTENT(OUT) :: gm_hxave,g0_hxave,gp_hxave
     END SUBROUTINE get_matrix_x_hxave

  END INTERFACE

  REAL(DP), INTENT(IN) :: dt
  REAL(DP), DIMENSION(nmmax,1:5), INTENT(OUT) :: ma_x
  REAL(DP), DIMENSION(nmmax), INTENT(IN) :: ffx,ffx0,ffxm,hxxxp,hxxx0,hxxxm,hxxxmm,hx0,hxm&
  ,fvdwx,fvdwx0,fvdwxm,fquadvdwx,fquadvdwx0,fquadvdwxm,fquad,fderquadx0,fderquadxm&
  ,hxxx,hx,coeff_one,coeff_three,hxavem,hxave0,hxavep

  REAL(DP) :: gpp_c(nmmax),gp_c(nmmax),g0_c(nmmax),gm_c(nmmax),gmm_c(nmmax)
  REAL(DP) :: gpp_quadterm(nmmax), gp_quadterm(nmmax), g0_quadterm(nmmax), gm_quadterm(nmmax)&
   ,gmm_quadterm(nmmax)
  REAL(DP) :: gp_vdw(nmmax),g0_vdw(nmmax),gm_vdw(nmmax)
  REAL(DP) :: gp_quadvdw(nmmax),g0_quadvdw(nmmax),gm_quadvdw(nmmax)
  REAL(DP) :: gm_hxave(nmmax),g0_hxave(nmmax),gp_hxave(nmmax)
  INTEGER(I4B) :: k

  CALL get_matrix_x_c(ffx,ffx0,ffxm,hxxx,hxxxp,hxxx0,hxxxm,hxxxmm,gpp_c,gp_c,g0_c,gm_c,gmm_c)

  CALL get_matrix_x_quadterm(fquad,fderquadx0,fderquadxm,hxxx,hxxxp,hxxx0,hxxxm,hxxxmm&
  ,gpp_quadterm,gp_quadterm,g0_quadterm,gm_quadterm,gmm_quadterm)

   CALL get_matrix_x_vdw(fvdwx,fvdwx0,fvdwxm,hx,hx0,hxm,gp_vdw,g0_vdw,gm_vdw)

  CALL get_matrix_x_quadvdw(fquadvdwx,fquadvdwx0,fquadvdwxm,hx,hx0,hxm&
  ,gp_quadvdw,g0_quadvdw,gm_quadvdw)

  CALL get_matrix_x_hxave(hxavem,hxave0,hxavep,gm_hxave,g0_hxave,gp_hxave)

IF(l2.EQ.0.0d0)THEN ! this IF puts back the code in the original form when l2=0, first order in time equation
    DO k=1,n
     ma_x(k,1)=(1.0d0-theta)*dt*(gmm_c(k) + bb*gmm_quadterm(k)) + l1*(gmm_c(k))

     ma_x(k,2)=(1.0d0-theta)*dt*(gm_c(k)+bb*gm_quadterm(k)+gm_vdw(k)+bb*gm_quadvdw(k))&
     +l1*(gm_c(k) + gm_vdw(k)) + coeff_three(k)*gm_hxave(k)

     ma_x(k,3)=coeff_one(k)+(1.0d0-theta)*dt*(g0_c(k)+bb*g0_quadterm(k)+g0_vdw(k)&
     +bb*g0_quadvdw(k)) + coeff_three(k)*g0_hxave(k) + l1*(g0_c(k) + g0_vdw(k))

     ma_x(k,4)=(1.0d0-theta)*dt*(gp_c(k)+bb*gp_quadterm(k)+gp_vdw(k)+bb*gp_quadvdw(k))&
     +l1*(gp_c(k) + gp_vdw(k)) + coeff_three(k)*gp_hxave(k)

     ma_x(k,5)=(1.0d0-theta)*dt*(gpp_c(k) + bb*gpp_quadterm(k))  + l1*(gpp_c(k))
    END DO

ELSE ! this is the latest version of the code, for l2 NOT zero, second order in time equation

      IF(t.EQ.t0)THEN ! This creates a slightly different matrix only for the first time step, n=0, accounting for the initial h^(-1)=h^1
           DO k=1,n
             ma_x(k,1)=(1.0d0-theta)*(dt**2)*(gmm_c(k) + bb*gmm_quadterm(k))&
             +l1*dt*(gmm_c(k)) + l2*dt*bb*(gmm_quadterm(k))

             ma_x(k,2)=(1.0d0-theta)*(dt**2)*(gm_c(k)+bb*gm_quadterm(k)+gm_vdw(k)+bb*gm_quadvdw(k))&
             +l1*dt*(gm_c(k) + gm_vdw(k))+l2*dt*bb*(gm_quadvdw(k) + gm_quadterm(k))&
             + dt*coeff_three(k)*gm_hxave(k)

             ma_x(k,3)=2.0d0*l2 + dt*coeff_one(k) + (1.0d0-theta)*(dt**2)*(g0_c(k)+bb*g0_quadterm(k)+g0_vdw(k)&
             +bb*g0_quadvdw(k)) + dt*coeff_three(k)*g0_hxave(k) + l1*dt*(g0_c(k) + g0_vdw(k))&
             +l2*dt*bb*(g0_quadvdw(k) + g0_quadterm(k))

             ma_x(k,4)=(1.0d0-theta)*(dt**2)*(gp_c(k)+bb*gp_quadterm(k)+gp_vdw(k)+bb*gp_quadvdw(k))&
             +l1*dt*(gp_c(k) + gp_vdw(k))+l2*dt*bb*(gp_quadvdw(k) + gp_quadterm(k))&
             + dt*coeff_three(k)*gp_hxave(k)

             ma_x(k,5)=(1.0d0-theta)*(dt**2)*(gpp_c(k) + bb*gpp_quadterm(k))&
             +l1*dt*(gpp_c(k)) + l2*dt*bb*(gpp_quadterm(k))
           END DO
      ELSE ! this is the general case, for any time step n>=1, the difference is that there is only l2 in the main diagonal, instead of 2*l2
           DO k=1,n
             ma_x(k,1)=(1.0d0-theta)*(dt**2)*(gmm_c(k) + bb*gmm_quadterm(k))&
             +l1*dt*(gmm_c(k)) + l2*dt*bb*(gmm_quadterm(k))

             ma_x(k,2)=(1.0d0-theta)*(dt**2)*(gm_c(k)+bb*gm_quadterm(k)+gm_vdw(k)+bb*gm_quadvdw(k))&
             +l1*dt*(gm_c(k) + gm_vdw(k))+l2*dt*bb*(gm_quadvdw(k) + gm_quadterm(k))&
             + dt*coeff_three(k)*gm_hxave(k)

             ma_x(k,3)= l2 + dt*coeff_one(k) + (1.0d0-theta)*(dt**2)*(g0_c(k)+bb*g0_quadterm(k)+g0_vdw(k)&
             +bb*g0_quadvdw(k)) + dt*coeff_three(k)*g0_hxave(k) + l1*dt*(g0_c(k) + g0_vdw(k))&
             +l2*dt*bb*(g0_quadvdw(k) + g0_quadterm(k))

             ma_x(k,4)=(1.0d0-theta)*(dt**2)*(gp_c(k)+bb*gp_quadterm(k)+gp_vdw(k)+bb*gp_quadvdw(k))&
             +l1*dt*(gp_c(k) + gp_vdw(k))+l2*dt*bb*(gp_quadvdw(k) + gp_quadterm(k))&
             + dt*coeff_three(k)*gp_hxave(k)

             ma_x(k,5)=(1.0d0-theta)*(dt**2)*(gpp_c(k) + bb*gpp_quadterm(k))&
             +l1*dt*(gpp_c(k)) + l2*dt*bb*(gpp_quadterm(k))
          END DO
     END IF

  END IF

  RETURN
END SUBROUTINE get_matrix_x2
