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
