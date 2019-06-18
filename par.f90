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
