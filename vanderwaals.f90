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


MODULE vanderwaals
  USE nrtype
  USE paras, ONLY: hstar, nvdw, mvdw
  IMPLICIT NONE
  
  !!$ this small module computes different terms involving the van der Waals term

  CONTAINS

  !!$ This term is the van der Waals term not multiplied by anything, it's need in Q and R. 
  !!$ This is in the notation Pi' (already differentiated).
  REAL(DP) FUNCTION vdwsimple(h)
    IMPLICIT NONE
    REAL(DP), INTENT(IN) :: h 
    
    vdwsimple=-dble(nvdw)*(hstar**nvdw)*(h**(-nvdw-1)) + dble(mvdw)*(hstar**mvdw)*(h**(-mvdw-1))

    RETURN
  END FUNCTION vdwsimple

  !!$ This term is the van der Waals term pre-multiplied by h^3
  REAL(DP) FUNCTION vdw(h)
    IMPLICIT NONE
    REAL(DP), INTENT(IN) :: h
 
    vdw=-dble(nvdw)*(hstar**nvdw)*(h**(-nvdw+2)) + dble(mvdw)*(hstar**mvdw)*(h**(-mvdw+2))

    RETURN
  END FUNCTION vdw

  !!$ This term is the van der Waals term pre-multiplied by the quadratic term h^2
  REAL(DP) FUNCTION quadvdw(h)
    IMPLICIT NONE
    REAL(DP), INTENT(IN) :: h
    
    quadvdw=-dble(nvdw)*(hstar**nvdw)*(h**(-nvdw+1)) + dble(mvdw)*(hstar**mvdw)*(h**(-mvdw+1))
    
    RETURN
  END FUNCTION quadvdw

  !!$ This term is the derivative of the van der Waals term pre-multiplied by h^3
  REAL(DP) FUNCTION dvdw(h)
    IMPLICIT NONE
    REAL(DP), INTENT(IN) :: h
    
    dvdw=dble(nvdw)*(dble(nvdw)-2.0d0)*(hstar**nvdw)*(h**(-nvdw+1)) - dble(mvdw)*(dble(mvdw)-2.0d0)*(hstar**mvdw)*(h**(-mvdw+1))
    
    RETURN
  END FUNCTION dvdw

  !!$ This term is the derivative of the van der Waals term pre-multiplied by the quadratic term h^2
  REAL(DP) FUNCTION dquadvdw(h)
    IMPLICIT NONE
    REAL(DP), INTENT(IN) :: h
    
    dquadvdw=dble(nvdw)*(dble(nvdw)-1.0d0)*(hstar**nvdw)*(h**(-nvdw)) - dble(mvdw)*(dble(mvdw)-1.0d0)*(hstar**mvdw)*(h**(-mvdw))
    
    RETURN
  END FUNCTION dquadvdw

END MODULE vanderwaals

