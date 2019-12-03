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


SUBROUTINE nonlin_quadterm(res,fquad)
  USE nrtype
  USE domain_space, ONLY: n, nmmax
  USE mobility, ONLY: quadterm
  IMPLICIT NONE

  !!$c	get quadratic part for capilary term and vdW
  !!$
  !!$  fquad(k,l): fquad at (k-1/2, l)
  !!$

  REAL(DP), DIMENSION(nmmax), INTENT(IN) :: res
  REAL(DP), DIMENSION(nmmax), INTENT(OUT) :: fquad

  INTEGER(I4B) np,k,km

  np=n+1

  !!$  fquad(k,l): fquad at (k-1/2, l)

  fquad(1)=quadterm(res(1))
  
  DO k=2,n
    km=k-1
    fquad(k)=(quadterm(res(k))+quadterm(res(km)))/2.0d0
  END DO
  
  fquad(np)=quadterm(res(n))


  RETURN
END SUBROUTINE nonlin_quadterm
