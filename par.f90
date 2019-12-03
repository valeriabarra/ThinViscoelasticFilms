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


SUBROUTINE par
  USE domain_space
  USE paras
  IMPLICIT NONE

  if((nmmax-n).lt.1)then
    WRITE(6,*)'NOT ENOUGH GRID POINTS IN X'
    STOP
  end if

  !!$  this subroutine takes care of opening the output files and writing into the status one
 
  OPEN(3,file=folder//'range.dat')
  OPEN(6,file=folder//'status.dat')
  OPEN(8,file=folder//'thickness.dat')
  OPEN(7,file=folder//'mass.dat')
  
  WRITE(3,*)'VARIABLES = t,xmax,"max_value",xmin,"min_value"'

  WRITE(6,97)'switch_x0=',switch_x0
  WRITE(6,97)'switch_xL=',switch_xL
  WRITE(6,97)'switch_dom=',switch_dom
  WRITE(6,97)'init_switch=',init_switch
  WRITE(6,98)'coef_c=',coef_c
  WRITE(6,98)'coef_m=',coef_m
  WRITE(6,98)'coef_vdw=',coef_vdw
  WRITE(6,98)'surf_tens=',surf_tens
  WRITE(6,98)'kappa=',kappa
  WRITE(6,98),'contact=',contact_angle
  WRITE(6,98)'vdwM=',vdwM
  WRITE(6,98)'hstar=',hstar
  WRITE(6,98)'M=',M
  WRITE(6,98)'b=',bb
  WRITE(6,98)'l1=',l1
  WRITE(6,98)'l2=',l2
  WRITE(6,200)'tol_mass=',tolmass
  WRITE(6,98)'percent=',percent
  WRITE(6,98)'x0=',x0
  WRITE(6,98)'xmax=',xmax
  WRITE(6,97)'n=',n
  
  WRITE(7,*)'VARIABLES = time,"mass"'
  
200  FORMAT(A15,1ES11.4)
98   FORMAT(A15,F12.5)
97   FORMAT(A15,I12)
95   FORMAT(A15,F12.2,F12.2,I12)

  CALL flush(1)
  CALL flush(3)
  CALL flush(6)
  CALL flush(7)
  CALL flush(8)
  CALL flush(15)

RETURN
END SUBROUTINE par
