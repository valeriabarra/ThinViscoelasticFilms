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

