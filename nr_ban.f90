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
