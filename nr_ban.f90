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


MODULE nr_ban
  INTERFACE
     SUBROUTINE penta(a,b,n,flag)
       USE nrtype
       REAL(DP), DIMENSION(:,:), INTENT(IN) :: a
       REAL(DP), DIMENSION(:), INTENT(INOUT) :: b
       INTEGER(I4B), INTENT(IN) :: n
       INTEGER(I4B), INTENT(OUT) :: flag
     END SUBROUTINE penta
  END INTERFACE
  INTERFACE
     SUBROUTINE bandec(a,n,m1,m2,al,indx,d,flag)
       USE nrtype
       INTEGER(I4B), INTENT(IN) :: n,m1,m2
       INTEGER(I4B), DIMENSION(:), INTENT(OUT) :: indx
       INTEGER(I4B), INTENT(OUT) :: flag
       REAL(DP), INTENT(OUT) :: d
       REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: a
       REAL(DP), DIMENSION(:,:), INTENT(OUT) :: al
     END SUBROUTINE bandec
  END INTERFACE
  INTERFACE
     SUBROUTINE banbks(a,n,m1,m2,al,indx,b)
       USE nrtype
       INTEGER(I4B), INTENT(IN) :: n,m1,m2
       INTEGER(I4B), DIMENSION(:), INTENT(IN) :: indx
       REAL(DP), DIMENSION(:,:), INTENT(IN) :: a,al
       REAL(DP), DIMENSION(:), INTENT(INOUT) :: b
     END SUBROUTINE banbks
  END INTERFACE
END MODULE nr_ban
