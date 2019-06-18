SUBROUTINE new_time(res,resold,dt,Qi,Ri,newton_iter,stab_iter,massold,massinit,mass_flag,tau_21)
  USE nrtype
  USE domain_space, ONLY: n, nmmax
  USE domain_time
  USE paras, ONLY: fac_dec, dt_min, percent, folder, tolmass, dx
  IMPLICIT NONE

  !!$  This subroutine allows us to advance in time.
  !!$  It solves for the next time step, using Newton's iteration.
  !!$  
  !!$  Input:
  !!$          res: solution at t=t_n
  !!$           dt: trial time step
  !!$  
  !!$  Output:
  !!$          res: solution at t=t_{n+1}
  !!$           dt: time step from t_n to t_{n+1}
  !!$  newton_iter: number of newton's iteration

  INTERFACE
     SUBROUTINE newton(res_guess,Qi, Ri,Qinew,Rinew,resold,dt,newton_iter,iflag_newton,tau_21)
       USE nrtype
       USE nr_ban
       USE domain_time
       USE domain_space, ONLY: n, nmmax
       USE paras
       USE vanderwaals
       USE mobility
       IMPLICIT NONE
       REAL(DP), DIMENSION(nmmax), INTENT(INOUT) :: res_guess
	   REAL(DP), DIMENSION(n), INTENT(INOUT) :: tau_21
       REAL(DP), DIMENSION(nmmax), INTENT(IN) :: Qi, Ri,resold
       REAL(DP), DIMENSION(nmmax), INTENT(INOUT) :: Qinew, Rinew
       REAL(DP), INTENT(IN) :: dt
       INTEGER(I4B), INTENT(OUT) :: newton_iter,iflag_newton
     END SUBROUTINE newton

  END INTERFACE



  REAL(DP), DIMENSION(nmmax), INTENT(INOUT) :: res, resold
  REAL(DP), DIMENSION(nmmax), INTENT(INOUT) :: Qi, Ri
  REAL(DP), DIMENSION(n), INTENT(INOUT) :: tau_21
  REAL(DP), INTENT(INOUT) :: dt
  INTEGER(I4B), INTENT(OUT) :: newton_iter
  INTEGER(I4B), INTENT(INOUT) :: stab_iter
  REAL(DP), INTENT(INOUT) :: massold
  REAL(DP), INTENT(IN) :: massinit
  INTEGER(I4B), INTENT(INOUT) :: mass_flag

  REAL(DP) :: res_guess(nmmax),Qinew(nmmax),Rinew(nmmax)
  INTEGER(I4B) :: iflag_newton,k
  REAL(DP) :: diff_min
  REAL(DP) :: massnew, massdiff, massdiffold, massdiffnew

!!$  Guess for t_{n+1}
  res_guess(1:n)=res(1:n)

  iflag_newton=1


  DO WHILE (iflag_newton.EQ.1)

     iflag_newton=0
     newton_iter=1
	 
    CALL newton(res_guess,Qi, Ri,Qinew,Rinew,resold,dt,newton_iter,iflag_newton,tau_21)
	 
    diff_min=DABS(MINVAL(res_guess(1:n))-MINVAL(resold(1:n)))

     IF(iflag_newton.EQ.1)THEN
        res_guess(1:n)=res(1:n)
        dt=dt/fac_dec
        stab_iter=0
        IF(dt.LE.dt_min)THEN
           WRITE(6,*)'Min. time step reached, dt=', dt
             OPEN(9,FILE=folder//'resMindt.dat')
             DO k=1,n
               WRITE(9,101)res(k)
             END DO  
           CALL FLUSH(9)
           STOP
        END IF
		
!         IF((diff_min.GE.percent))THEN !case of insuccess
!         res_guess(1:n)=res(1:n)
!         dt=dt/fac_dec
!         WRITE(6,99)'H diminished too much, diffmin=', diff_min
!         stab_iter=0
!            IF(dt.LE.dt_min)THEN
!              WRITE(6,*)'Min. time step reached, dt=', dt
!              OPEN(9,FILE=folder//'resMindt.dat')
!              DO k=1,n
!                WRITE(9,101)res(k)
!              END DO  
!            CALL FLUSH(9)
!            STOP
!            END IF
!        END IF
     END IF
	 	  
 END DO

 ! Check for mass loss !
 massnew=SUM(res_guess(1:n))*dx
 massdiffold=(massinit-massold)/massinit
 massdiffnew=(massinit-massnew)/massinit

massdiff=DABS(massnew-massold)/dt
! massdiff=DABS(massdiffold-massdiffnew)/dt

 IF(massdiff.GE.tolmass)THEN
    dt=dt/fac_dec
    WRITE(6,99)'Mass loss too big=',massdiff, 'dt=', dt
    mass_flag = 1
    IF(dt.LE.dt_min)THEN
        WRITE(6,*)'Min. time step reached, dt=', dt
          OPEN(9,FILE=folder//'resMindt.dat')
          DO k=1,n
            WRITE(9,101)res(k)
          END DO
         CALL FLUSH(9)
        STOP
    END IF

 ELSE ! this is the part that all went well
  stab_iter=stab_iter+1


  !!$  Newton converged !!$ So we update all the variables
  resold(1:n)= res(1:n)
  res(1:n)=res_guess(1:n)

  Qi(1:n)=Qinew(1:n)
  Ri(1:n)=Rinew(1:n)

  ! update the mass now
  massold=massnew
  !!$
  mass_flag = 0

  END IF



101 FORMAT(ES25.15)
99  FORMAT(A28,ES19.10,A15,ES11.2)

  RETURN
END SUBROUTINE new_time
