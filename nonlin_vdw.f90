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


SUBROUTINE nonlin_vdw(res,fvdwx,fvdwxsimple)

  USE nrtype
  USE domain_space, ONLY: n, nmmax
  USE vanderwaals, ONLY: vdw, vdwsimple
  USE mobility, ONLY: mob
  IMPLICIT NONE

  !!$c    get nonlinear part for van der waals terms
  !!$
  !!$  fvdwx(k): fvdwx at (k-1/2) (the one originally in the code, already multiplied by h^3)
  !!$  fvdwxsimple(k): fvdwxsimple at (k-1/2) (needed in Q and R)

  REAL(DP), DIMENSION(nmmax), INTENT(IN) :: res
  REAL(DP), DIMENSION(nmmax), INTENT(OUT) :: fvdwx,fvdwxsimple

  INTEGER(I4B) np,k,km

  np=n+1

  fvdwx(1)=vdw(res(1))
  fvdwxsimple(1)=vdwsimple(res(1)) ! needed in Q and R
  
  DO k=2,n
    km=k-1
    fvdwx(k)=(vdw(res(k))+vdw(res(km)))/2.0d0
    fvdwxsimple(k)=(vdwsimple(res(k)) + vdwsimple(res(km)))/2.0d0 ! needed in Q and R
  END DO
  
  fvdwx(np)=vdw(res(n))
  fvdwxsimple(np)=vdwsimple(res(n)) ! needed in Q and R

  RETURN
END SUBROUTINE nonlin_vdw
