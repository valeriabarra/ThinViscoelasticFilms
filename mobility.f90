MODULE mobility
  USE nrtype
  IMPLICIT NONE

  !!$ this small module computes different terms involving powers of h

  CONTAINS

  REAL(DP) FUNCTION mob(h)
    IMPLICIT NONE
    REAL(DP), INTENT(IN) :: h
    
    mob=h**3
    
    RETURN
  END FUNCTION mob

  REAL(DP) FUNCTION dmob(h)
    IMPLICIT NONE
    REAL(DP), INTENT(IN) :: h
    
    dmob=3.0d0*(h**2)
    
    RETURN
  END FUNCTION dmob

  REAL(DP) FUNCTION quadterm(h)
    IMPLICIT NONE
    REAL(DP), INTENT(IN) :: h
    
    quadterm=h**2
    
    RETURN
  END FUNCTION quadterm

  REAL(DP) FUNCTION dquadterm(h)
    IMPLICIT NONE
    REAL(DP), INTENT(IN) :: h
    
    dquadterm=2.0d0*h
    
    RETURN
  END FUNCTION dquadterm

END MODULE mobility
