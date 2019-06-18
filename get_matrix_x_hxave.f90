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

SUBROUTINE get_matrix_x_hxave(hxavem,hxave0,hxavep,gm_hxave,g0_hxave,gp_hxave)
  USE nrtype
  USE domain_space, ONLY: n, nmmax
  USE paras, ONLY: dx
  IMPLICIT NONE

  !!$ this subroutine assembles the entries asscoiated to the hxave term in the 
  !!$ pentadiagonal matrix (only three diagonals of it)

  REAL(DP), DIMENSION(nmmax), INTENT(IN) :: hxavem,hxave0,hxavep
  REAL(DP), DIMENSION(nmmax), INTENT(OUT) :: gm_hxave,g0_hxave,gp_hxave

  INTEGER(I4B) :: k, kp, np
  REAL(DP) :: ccx

  ccx=1.0d0/dx

  DO k=1,n
     gm_hxave(k)=ccx*(hxavem(k))
     g0_hxave(k)=ccx*(hxave0(k))
     gp_hxave(k)=ccx*(hxavep(k))
  END DO

  RETURN
END SUBROUTINE get_matrix_x_hxave
