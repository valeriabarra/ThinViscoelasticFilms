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

MODULE paras
  USE nrtype
  IMPLICIT NONE
  SAVE

  !!$ this module contains several of the physical parameters that needs to be defined
  
  !!$ output folder
  CHARACTER(len=*), PARAMETER :: folder='/home/user/Output/000/' ! needs to be changed with your global path

  !!$ initialization of Tolerance for Steady State !
  REAL(DP), PARAMETER ::  tolsteady=1.0d-10
  REAL(DP), PARAMETER :: bb=0.0d0 ! bb is the slip coefficient (b=0 means no slip)
  REAL(DP), PARAMETER :: coef_c=1.0d0/3.0d0
  REAL(DP), PARAMETER :: coef_slip=bb
  REAL(DP), PARAMETER :: hstar=5.0d-3
  INTEGER(I4B), PARAMETER :: switch_dom=0  ! 1 indicates half domain, 0 the whole
  INTEGER(I4B), PARAMETER :: init_switch=0 ! 0 indicates default initial condition from zero,
                                           ! 1 indicates that the solution starts as a continuation of old solution
										   ! 2 indicates initial rivulet/droplet
  REAL(DP), PARAMETER :: tolmass=1.0d-1

  !!$ Epsilon (or delta) (not the small parameter, but the small amplitude of initial perturbation, but they can be considered of the same order)
  REAL(DP), PARAMETER :: E=1.0d-2

  !!$ Jeffreys model viscoelastic parameters
  !!$ l1 is the relaxation time
  REAL(DP), PARAMETER :: l1=0.0d0
  !!$ l2 is the retardation time
  REAL(DP), PARAMETER :: l2=0.0d0

  !!$  parameters for Newton's iteration
  REAL(DP), PARAMETER :: dt_max=1.0d1
  REAL(DP), PARAMETER :: dt_min=1.0d-15
  REAL(DP), PARAMETER :: fac_inc=1.01d0 ! time step incremental factor when Newton converges
  REAL(DP), PARAMETER :: fac_dec=2.0d0  ! time step decremental factor (cutback) for when Newton does not converge
  REAL(DP), PARAMETER :: percent=1.0d-5
  INTEGER(I4B) :: n_time=0
  REAL(DP), PARAMETER :: theta=0.5d0 ! this is theta-scheme = 1/2, ie Crank-Nicolson
  REAL(DP), PARAMETER :: eps_newt=1.0d-9
  REAL(DP), PARAMETER :: tol=1.0d-4
  INTEGER(I4B), PARAMETER :: n_iter_max=10 ! max no. iterations for each Newton's itertion it was 10
  INTEGER(I4B) :: iflag_fine=0
  INTEGER(I4B) :: n_increase=10
  INTEGER(I4B) :: nout_save=-1
  INTEGER(I4B) :: ncount_four=300

  !!$  Exponents of step sizes
  REAL(DP) :: dx, dx2, dx3, dx4, dxh

  !!$ droplet for initial condition
  REAL(DP), PARAMETER :: theta_init=45.0d0
  REAL(DP), PARAMETER :: char_height=1.0d0
  REAL(DP), PARAMETER :: char_length=1.0d0
  REAL(DP), PARAMETER :: r_phys=3.0d0

  !!$  constants in magnetic term
  REAL(DP), PARAMETER :: mu0=0.0d0
  REAL(DP), PARAMETER :: chiSI=0.0d0
  REAL(DP), PARAMETER :: surf_tens=1.0d0
  REAL(DP), PARAMETER :: psib=0.0d0
  REAL(DP), PARAMETER :: lambda=0.0d0
  REAL(DP), PARAMETER :: lambdam=0.0d0
  REAL(DP), PARAMETER :: coef_m=0.0d0
  REAL(DP), PARAMETER :: beta=0.0d0
  REAL(DP), PARAMETER :: M=0.0d0

  !!$  constants in van der waals term  
  INTEGER(I4B), PARAMETER :: nvdw=3
  INTEGER(I4B), PARAMETER :: mvdw=2
  REAL(DP), PARAMETER :: vdwM=DBLE(nvdw-mvdw)/DBLE((mvdw-1)*(nvdw-1))
  REAL(DP), PARAMETER :: contact_angle=15.0d0
  REAL(DP), PARAMETER :: kappa=surf_tens*((1.0d0-COS(contact_angle*PI/180.0d0))/(vdwM*hstar*E**2))
  REAL(DP), PARAMETER :: coef_vdw=(1.0d0/3.0d0)*kappa ! this way vdw is already divided by three considering its version multiplied by h^3
                                                      ! but when I need it simple or multiplied by h^2 I multiply it again by 3.0d0


  !!$  boundary conditions
  !!$  switch_x0 = 1 for h=h(t=0), h'=0
  !!$  switch_x0 = 2 for h'=h'''=0

  INTEGER(I4B), PARAMETER :: switch_x0=2
  INTEGER(I4B), PARAMETER :: switch_xL=2




END MODULE paras
