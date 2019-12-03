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


MODULE mobility
  USE nrtype
  IMPLICIT NONE

  !!$ this small module computes different terms involving powers of h

  CONTAINS

  REAL(DP) FUNCTION mob(h)
    IMPLICIT NONE
    REAL(DP), INTENT(IN) :: h
    
    mob=h**3
    
    RETURN
  END FUNCTION mob

  REAL(DP) FUNCTION dmob(h)
    IMPLICIT NONE
    REAL(DP), INTENT(IN) :: h
    
    dmob=3.0d0*(h**2)
    
    RETURN
  END FUNCTION dmob

  REAL(DP) FUNCTION quadterm(h)
    IMPLICIT NONE
    REAL(DP), INTENT(IN) :: h
    
    quadterm=h**2
    
    RETURN
  END FUNCTION quadterm

  REAL(DP) FUNCTION dquadterm(h)
    IMPLICIT NONE
    REAL(DP), INTENT(IN) :: h
    
    dquadterm=2.0d0*h
    
    RETURN
  END FUNCTION dquadterm

END MODULE mobility
