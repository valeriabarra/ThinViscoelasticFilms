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
