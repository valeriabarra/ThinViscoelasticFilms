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


SUBROUTINE nonlin_func2(res,func,func_four,func_six,ffx,hxxx,hx,hxave,hxxxave,fvdwx,fquadvdwx,fquad,fvdwxsimple)
  USE nrtype
  USE domain_space, ONLY: n, nmmax
  USE paras, ONLY: bb
  IMPLICIT NONE

  INTERFACE
    SUBROUTINE nonlin_cap(res,ffx)
      USE nrtype
      USE domain_space, ONLY: n, nmmax
      USE mobility, ONLY: mob
      IMPLICIT NONE
      REAL(DP), DIMENSION(nmmax), INTENT(IN) :: res
      REAL(DP), DIMENSION(nmmax), INTENT(OUT) :: ffx
    END SUBROUTINE nonlin_cap

    SUBROUTINE nonlin_quadterm(res,fquad)
      USE nrtype
      USE domain_space, ONLY: n, nmmax
      USE mobility, ONLY: quadterm
      IMPLICIT NONE
      REAL(DP), DIMENSION(nmmax), INTENT(IN) :: res
      REAL(DP), DIMENSION(nmmax), INTENT(OUT) :: fquad
    END SUBROUTINE nonlin_quadterm

    SUBROUTINE nonlin_quadvdw(res,fquadvdwx)
      USE nrtype
      USE domain_space, ONLY: n, nmmax
      USE vanderwaals, ONLY: quadvdw
      IMPLICIT NONE
      REAL(DP), DIMENSION(nmmax), INTENT(IN) :: res
      REAL(DP), DIMENSION(nmmax), INTENT(OUT) :: fquadvdwx
    END SUBROUTINE nonlin_quadvdw

    SUBROUTINE nonlin_vdw(res,fvdwx,fvdwxsimple)
      USE nrtype
      USE domain_space, ONLY: n, nmmax
      USE vanderwaals, ONLY: vdw, vdwsimple
      USE mobility, ONLY: mob
      REAL(DP), DIMENSION(nmmax), INTENT(IN) :: res
      REAL(DP), DIMENSION(nmmax), INTENT(OUT) :: fvdwx,fvdwxsimple
    END SUBROUTINE nonlin_vdw
   
    SUBROUTINE hhh_cap(res,hxxx,hx,hxxxave,hxave)
      USE nrtype
      USE domain_space, ONLY: n, nmmax
      IMPLICIT NONE
      REAL(DP), DIMENSION(nmmax), INTENT(IN) :: res
      REAL(DP), DIMENSION(nmmax), INTENT(OUT) :: hxxx,hx,hxave,hxxxave
    END SUBROUTINE hhh_cap
    
    SUBROUTINE term_c(ffx,hxxx,func_c)
      USE nrtype
      USE domain_space, ONLY: n,nmmax
      USE paras, ONLY: coef_c,dx4
      IMPLICIT NONE
      REAL(DP), DIMENSION(nmmax), INTENT(IN) :: ffx,hxxx
      REAL(DP), DIMENSION(nmmax), INTENT(OUT) :: func_c
    END SUBROUTINE term_c

    SUBROUTINE term_quad(fquad,hxxx,func_quad)
      USE nrtype
      USE domain_space, ONLY: n,nmmax
      USE paras, ONLY: coef_slip,dx4
      IMPLICIT NONE
      REAL(DP), DIMENSION(nmmax), INTENT(IN) :: fquad,hxxx
      REAL(DP), DIMENSION(nmmax), INTENT(OUT) :: func_quad
    END SUBROUTINE term_quad

    SUBROUTINE term_vdw(fvdwx,hx,func_vdw)
      USE nrtype
      USE domain_space, ONLY: n,nmmax
      USE paras, ONLY: coef_vdw,dx2
      IMPLICIT NONE
      REAL(DP), DIMENSION(nmmax), INTENT(IN) :: fvdwx,hx
      REAL(DP), DIMENSION(nmmax), INTENT(OUT) :: func_vdw
    END SUBROUTINE term_vdw

    SUBROUTINE term_quadvdw(fquadvdwx,hx,func_quadvdw)
      USE nrtype
      USE domain_space, ONLY: n,nmmax
      USE paras, ONLY: coef_slip,coef_vdw,dx2
      IMPLICIT NONE
      REAL(DP), DIMENSION(nmmax), INTENT(IN) :: fquadvdwx,hx
      REAL(DP), DIMENSION(nmmax), INTENT(OUT) :: func_quadvdw
    END SUBROUTINE term_quadvdw
  END INTERFACE

  REAL(DP), DIMENSION(nmmax), INTENT(IN) :: res
  REAL(DP), DIMENSION(nmmax), INTENT(OUT) :: func, ffx, hxxx, hx, hxave,hxxxave, fvdwx, fquadvdwx, fquad,func_four,func_six&
  ,fvdwxsimple

  REAL(DP) :: func_c(nmmax), func_vdw(nmmax), func_quad(nmmax),func_quadvdw(nmmax)

  CALL nonlin_cap(res,ffx)
  CALL nonlin_quadterm(res,fquad)
  CALL nonlin_vdw(res,fvdwx,fvdwxsimple)
  CALL nonlin_quadvdw(res,fquadvdwx)
  CALL hhh_cap(res,hxxx,hx,hxxxave,hxave)
  CALL term_c(ffx,hxxx,func_c)
  CALL term_quad(fquad,hxxx,func_quad)! fixed first argument was ffx, wrong!
  CALL term_vdw(fvdwx,hx,func_vdw)
  CALL term_quadvdw(fquadvdwx,hx,func_quadvdw)

  func(1:n)=func_c(1:n)+bb*func_quad(1:n)+func_vdw(1:n)+bb*func_quadvdw(1:n) ! all four pieces needed in f = term 5
  func_four(1:n)=func_c(1:n)+func_vdw(1:n) ! only capillary terms needed in g = term 4
  func_six(1:n)=bb*func_quad(1:n)+bb*func_quadvdw(1:n) ! only quadratic slip terms needed in m = term 6

  RETURN
END SUBROUTINE nonlin_func2
