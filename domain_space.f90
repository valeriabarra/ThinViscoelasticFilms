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

MODULE domain_space
  USE nrtype
  IMPLICIT NONE
  SAVE

  !!$  Spatial domain in x: [x0, xmax] with n grid points
  !!$  nmmax should be larger than n

  REAL(DP), PARAMETER :: x0=0.0d0
  REAL(DP), PARAMETER :: xmax=10.0d0
  INTEGER(I4B), PARAMETER :: n=1000, nmmax=n+1 !! Tot. no. of points (+ 2 ghost pts)
  REAL(DP) :: x(nmmax)
  REAL(DP) :: bdx0, bdxL

END MODULE domain_space
