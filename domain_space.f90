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
