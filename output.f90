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


SUBROUTINE output(res,t,tau_21)
  USE nrtype
  USE domain_time, ONLY: t_end
  USE domain_space, ONLY: x, n, nmmax
  USE paras, ONLY: dx, bb, folder
  IMPLICIT NONE

  REAL(DP), DIMENSION(nmmax), INTENT(IN) :: res
  REAL(DP), DIMENSION(n), INTENT(IN) ::tau_21
  REAL(DP), INTENT(INOUT) :: t

  REAL(DP) :: mass
  INTEGER(I4B) :: k
  INTEGER(I4B), DIMENSION(1) :: k_max, k_min

  !!$  output film height profile

  DO k=1,n
    WRITE(8,99)x(k),res(k)
  END DO

  !!$ output shear stress component tau_21
    DO k=1,n
     WRITE(15,99)x(k),tau_21(k)
  END DO
  
  !!$ write the current output solution in a separate file with more digits
  OPEN(11,file=folder//'resLastOut.dat',STATUS='REPLACE',ACTION='WRITE')
  
  DO k=1,n
    WRITE(11,101)x(k),res(k)
  END DO
  
  CLOSE(11)
101 FORMAT(2(ES24.16))
99 FORMAT(2(ES15.5))

  !!$  output max and min height

  k_max=MAXLOC(res(1:n))
  k_min=MINLOC(res(1:n))
  WRITE(3,98)t,x(k_max),res(k_max),x(k_min),res(k_min)

98 format(5(ES17.9))

  !!$  output total mass

  mass=SUM(res(1:n))*dx
  WRITE(7,97)t,mass

 97 FORMAT(2(ES20.9))

  call flush(3)
  call flush(7)
  call flush(8)
  call flush(15)
  
  RETURN
END SUBROUTINE output
