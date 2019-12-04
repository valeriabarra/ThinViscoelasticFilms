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


subroutine penta(a22,b,n,iflag_solve)
  USE nrtype; USE nr_ban, ignore_me => penta
  IMPLICIT NONE
  
  !!$ this subroutine solves the banded pentadiagonal matrix with a direct solver

  REAL(DP), DIMENSION(:,:), INTENT(IN) :: a22
  REAL(DP), DIMENSION(:), INTENT(INOUT) :: b
  INTEGER(I4B), INTENT(IN) :: n
  INTEGER(I4B), INTENT(OUT) :: iflag_solve

  REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: a,al
  REAL(DP) :: d
  INTEGER(I4B), ALLOCATABLE, DIMENSION(:) :: indx
  INTEGER(I4B) :: statu

  allocate(a(n,5), stat=statu)
  allocate(al(n,2), stat=statu)
  allocate(indx(n), stat=statu)

  a(1:n,1:5)=a22(1:n,1:5)

  call bandec(a(1:n,1:5),n,2,2,al(1:n,1:5),indx(1:n),d,iflag_solve)

  IF(iflag_solve.EQ.1)THEN
     WRITE(*,*)'penta flag=1, NO SOLUTION'
     STOP
  END IF
  call banbks(a(1:n,1:5),n,2,2,al(1:n,1:5),indx(1:n),b(1:n))

  return
end subroutine penta

!!$
!!$ this is an auxiliary subroutine used by "penta"
!!$

SUBROUTINE bandec(a,n,m1,m2,al,indx,d,flag)
  USE nrtype
  IMPLICIT NONE
  REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: a
  INTEGER(I4B), INTENT(IN) :: m1,m2
  REAL(DP), DIMENSION(:,:), INTENT(OUT) :: al
  INTEGER(I4B), DIMENSION(:), INTENT(INOUT) :: indx
  REAL(DP), INTENT(OUT) :: d
  REAL(DP), PARAMETER :: TINY=1.0d-20
  INTEGER(I4B), INTENT(OUT) :: flag
  INTEGER(I4B) :: i,j,k,l,mm,n
  INTEGER(I4B), DIMENSION(1) :: im
  REAL(DP) :: dum

  mm=m1+m2+1
  flag=0

  l=m1
  do i=1,m1
     a(i,(m1+2-i-l):m1+i)=a(i,(m1+2-i):mm) ! it was up to :mm
     l=l-1
     a(i,(mm-l):mm)=0.0d0
  enddo

  d=1.0d0
  
  do k=1,n
     l=min(m1+k,n)
     im=maxloc(abs(a(k:l,1)))
     i=im(1)+k-1
     dum=a(i,1)
     if (abs(dum) <= TINY) then
        flag=1
        RETURN
     end if
     indx(k)=i
     if (i /= k) then
        d=-d
        do j=1,mm
           dum=a(k,j)
           a(k,j)=a(i,j)
           a(i,j)=dum
        enddo
     end if
     do i=k+1,l
        dum=a(i,1)/a(k,1)
        al(k,i-k)=dum
        a(i,1:mm-1)=a(i,2:mm)-dum*a(k,2:mm)
        a(i,mm)=0.0d0
     end do
  end do
END SUBROUTINE bandec

!!$
!!$ this is another auxiliary subroutine used by "penta"
!!$

SUBROUTINE banbks(a,n,m1,m2,al,indx,b)
  USE nrtype
  IMPLICIT NONE
  REAL(DP), DIMENSION(:,:), INTENT(IN) :: a,al
  INTEGER(I4B), INTENT(IN) :: m1,m2
  INTEGER(I4B), DIMENSION(:), INTENT(IN) :: indx
  REAL(DP), DIMENSION(:), INTENT(INOUT) :: b
  INTEGER(I4B) :: i,k,l,mm,n
  REAL(DP) :: dum

  mm=m1+m2+1
  do k=1,n
     l=min(n,m1+k)
     i=indx(k)
     if (i /= k) then 
        dum=b(k)
        b(k)=b(i)
        b(i)=dum
     end if
     b(k+1:l)=b(k+1:l)-al(k,1:l-k)*b(k)
  end do
  do i=n,1,-1
     l=min(mm,n-i+1)
     b(i)=(b(i)-dot_product(a(i,2:l),b(1+i:i+l-1)))/a(i,1)
  end do
END SUBROUTINE banbks
