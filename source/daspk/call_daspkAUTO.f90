! 9-Jul-2009 mrharper:
!	Added Lindemann reactions
! 6/26/08 gmagoon:
! filename: call_daspkAUTO.f90;
! changed name from CALL_DASPK to CALL_DASPKLAUTO
! to indicate that this version requires changes to RMG;
! changed from .f to .f90 and from "C" comments to "!"
! comments to address compilation issues; also changed "$" continued line
! to "& /n &" continued line 
PROGRAM CALL_DASPKAUTO

      IMPLICIT NONE

! IN THIS CODE WE ASSUME THAT THE 
! MAXIMUM # SPECIES = 1500 
! MAXIMUM # REACTIONS = 100,000
! MAXIMUM # THIRDBODYREACTIONS = 100
! MAXIMUM # TROEREACTIONS = 100
! MAXIMUM # LINDEREACTIONS = 100

      INCLUDE 'reaction.fh'


      INTEGER INFO(30), LRW, LIW, I, J, IDID, NSTEPS, IMPSPECIES, numiter, & 
     &     AUTOFLAG, ESPECIES, EREACTIONSIZE, SENSFLAG
      DOUBLE PRECISION Y(NEQMAX), YPRIME(NEQMAX), T, TOUT, RTOL, ATOL, &
     &     THERMO(SPMAX), TARGETCONC(50), THRESH, CORETHRESH
     ! 6/26/08 gmagoon:make auto arrays allocatable and double
     ! precision for KVEC
     INTEGER, DIMENSION(:), ALLOCATABLE :: NEREAC,NEPROD
     INTEGER, DIMENSION(:,:), ALLOCATABLE :: IDEREAC, IDEPROD
     DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: KVEC

      IDID = 0 !gmagoon 1/25/10: initialize IDID to zero (no meaning in terms of DASPK idid outputs) to avoid a situation where dassl is never called (and hence IDID is never assigned) due to edge exceeding flux at t=0
      OPEN (UNIT=12, FILE = 'SolverInput.dat', STATUS = 'OLD')

 101  Format(E24.15)

! READ THE NUMBER OF SPECIES; 6/26/08 gmagoon: added autoFlag,
! which will equal 1 for automatic (solve ODE until flux
! exceeds threshhold) or -1 for conventional operation
! (i.e. user-specification of time/conversion steps) 
! gmagoon 080509: added sensflag so that DASPK can use this to decide whether
! to do sensitivity or not; previously, it was using numiter to decide (which is
! OK in most cases, but not if we use auto or don't specify intermediate
! conversions)
      READ(12,*) NSTATE, NEQ, IMPSPECIES, NUMITER, AUTOFLAG, SENSFLAG
      NPARAM = INT(NEQ/NSTATE) - 1
      NSTATE = NSTATE+1
      NEQ = NSTATE*(NPARAM + 1)


      READ(12,*) (TARGETCONC(I), I=1,NUMITER) 


!     READ THE TIME AND TOUT
      READ(12,*) T, TOUT


!  READ THE CONCENTRATIONS
      if (t .eq. 0.0) then
! it is a new job and just have to read the state variables
! read the concentrations
         READ(12,*) (Y(I), i=1,nstate-1)
!     READ THE RATE OF CHANGE OF CONCENTRATIONS
         READ(12,*) (Yprime(I), i=1,nstate-1)
         do i=nstate,neq
            y(i) = 0
            yprime(i) = 0
         end do
!     READ THE INFO ARRAY
         READ(12,*) (INFO(I), I=1,30)
!     READ RTOL AND ATOL
         READ(12,*) RTOL, ATOL
         
!     READ THE THERMO OF SPECIES
         READ(12,*)(THERMO(I), I=1,NSTATE-1)
         
!     READ TEMPERATURE AND PRESSURE
         READ(12,*) TEMPERATURE, PRESSURE

         
         y(nstate) =  8.314*temperature/pressure/1e-6
         do i=1,nstate-1
            y(i) = y(i)*y(nstate)
         end do
            

         
!     READ INFORMATION ABOUT REACTIONS
         READ(12,*) REACTIONSIZE
         READ(12,*) (REACTIONARRAY(I),I=1,9*REACTIONSIZE)
         READ(12,*) (REACTIONRATEARRAY(I),I=1,5*REACTIONSIZE)

!     READ INFORMATION ABOUT THIRDBODYREACTIONS
         READ(12,*) THIRDBODYREACTIONSIZE
         READ(12,*)(THIRDBODYREACTIONARRAY(I),I=1, &
         &   20*THIRDBODYREACTIONSIZE)
         READ(12,*)(THIRDBODYREACTIONRATEARRAY(I),I=1, &
        &    16*THIRDBODYREACTIONSIZE)
         
!     READ INFORMATION ABOUT TROEREACTIONS
         READ(12,*) TROEREACTIONSIZE
         READ(12,*)(TROEREACTIONARRAY(I),I=1,21*TROEREACTIONSIZE)
         READ(12,*)(TROEREACTIONRATEARRAY(I),I=1,21*TROEREACTIONSIZE)

!     READ INFORMATION ABOUT LINDEREACTIONS
         READ(12,*) LINDEREACTIONSIZE
		 READ(12,*)(LINDEREACTIONARRAY(I),I=1,20*LINDEREACTIONSIZE)
		 READ(12,*)(LINDEREACTIONRATEARRAY(I),I=1,17*LINDEREACTIONSIZE)

! 6/26/08 gmagoon: if autoFlag = 1, read in additional information 
! specific to automatic time stepping
	IF (AUTOFLAG .EQ. 1) THEN
		! read the threshhold, corresponding to the value of
		! termination tolerance specified in condition.txt input file
	        !CORETHRESH is the threshhold for moving something from edge to core,
	        !whereas THRESH is the threshhold for terminating the run before the end
	        !THRESH >= CORETHRESH
		READ(12,*) THRESH, CORETHRESH
		! read the number of edge species and edge reactions
		READ(12,*) ESPECIES, EREACTIONSIZE
		! allocate memory for arrays
		ALLOCATE(NEREAC(EREACTIONSIZE), NEPROD(EREACTIONSIZE),&
	&	    IDEREAC(EREACTIONSIZE,3), IDEPROD(EREACTIONSIZE,3),&
	&            KVEC(EREACTIONSIZE))
		! read in the reaction parameters for each reaction;
		! a maximum of 3 products and 3 reactants is assumed
		! for each reaction; parameters read for each reaction
		! are: number of reactants, number of products,
		! three reactant ID's (integers from 1 to NSTATE-1),
		! three product ID's (integers from 1 to ESPECIES),
		! and the rate coefficient k, such that dCi/dt=k*Ca*Cb;
		! note that cases where abs. value of stoic. coeff. 
		! does not equal one are handled by using repeated
		! reactant/product IDs
		! note: use of KVEC rather than Arrhenius parameters
		! requires assumption that system is isothermal
		! (and isobaric for pressure dependence)
		DO I=1, EREACTIONSIZE
			READ(12,*) NEREAC(I),NEPROD(I),IDEREAC(I,1), & 
         &         	  IDEREAC(I,2),IDEREAC(I,3),IDEPROD(I,1), &
         &       	  IDEPROD(I,2), IDEPROD(I,3), KVEC(I)
		! alternative for reading Arrhenius parameters instead
		! of k values; this allows easier extention to
		! non-isothermal systems in the future
		! form: k=A*T^n*e^(Ea/(RT)) *?
		! units: A[=], n[=]dimensionless, Ea[=]
	  !		READ(12,*) NEREAC(I),NEPROD(I),IDEREAC(I,1), &
        ! &               IDEREAC(I,2),IDEREAC(I,3),IDEPROD(I,1), &
        ! &               IDEPROD(I,2),IDEPROD(I,3), AVEC(I), NVEC(I), &
        ! &               EAVEC(I)
	!      for isothermal, isobaric systems, (or just isothermal 
	!      systems if pressure dependence is not considered)
      !      the following may be calculated here once
	!      (versus calculating at every timestep)
	!      ...units may need adjustment:
	!    		KVEC(I)=AVEC(I)*TEMPERATURE**NVEC(I)*EXP(EAVEC(I) &
	!	&       /(8.314*TEMPERATURE))
		END DO
	END IF
        ! read constantConcentration data (if flag = 1 then the concentration of that species will not be integrated)
        ! there is one integer for each species (up to nstate-1), then the last one is for the VOLUME
        READ(12,*) (ConstantConcentration(I), i=1,nstate)
! 6/26/08 gmagoon: if t.ne.0, we are presumably not using AUTO method, since
! with AUTO method we would return to zero each time; alternative would be to
! also read and write AUTO parameters from/to variables.dat if AUTOFLAG = 1;
! not including AUTO parameters (as done here) could give error,
! but an error would indicate a problem with the AUTO method where it is not
! returning to zero after ODE solver finishes
      ELSE
         OPEN(UNIT=13, FILE = 'variables.dat', form='unformatted')
! the first line is neq read that
         read(13) idid
         read(13) t
         read(13) nstate
         read(13) nparam
         read(13) neq
         DO I=1,NEQ
            READ(13) Y(I)
            READ(13) YPRIME(I)
         END DO
         read(13) rtol, atol
         do i=1,30
            read(13) info(i)
         end do

!     c write the reaction info
         read(13) reactionsize
         do i=1,9*reactionsize
            read(13) reactionarray(i)
         end do
         do i=1,5*reactionsize
            read(13) reactionratearray(i)
         end do
         
         read(13) thirdbodyreactionsize
         do i=1,20*thirdbodyreactionsize
            read(13) thirdbodyreactionarray(i)
         end do

         do i=1,16*thirdbodyreactionsize
            read(13) thirdbodyreactionratearray(i)
         end do
         
         read(13) troereactionsize
         do i=1,21*troereactionsize
            read(13) troereactionarray(i)
         end do
         do i=1,21*troereactionsize
            read(13) troereactionratearray(i)
         end do

		 READ(13) LINDEREACTIONSIZE
		 DO I=1,20*LINDEREACTIONSIZE
			READ(13) LINDEREACTIONARRAY(I)
		 END DO
		 DO I=1,17*LINDEREACTIONSIZE
			READ(13) LINDEREACTIONRATEARRAY(I)
		 END DO

         do i=1,nstate-1
            read(13) thermo(i)
         end do
         read(13) temperature, pressure
         
         !gmagoon: even though concentration flags will be in SolverInput file, it will be simplest to write and read them from variables.dat
         DO i=1,nstate
            read(13) ConstantConcentration(i)
         END DO

         close(13)
      END IF
      
      CLOSE(12) !gmagoon 6/11/09: moved location for closing 12
      !(the input file) so that it is not closed before AUTO
      !information is read in
      
    !  if (NUMITER .eq. 1) then
    !     CALL SOLVEODE(Y, YPRIME, T, TOUT, INFO, RTOL, ATOL, IDID, &
    ! &        THERMO, IMPSPECIES, TARGETCONC(1),AUTOFLAG, &
    ! &     THRESH,ESPECIES,EREACTIONSIZE,NEREAC,NEPROD,IDEREAC, &
    ! &     IDEPROD,KVEC)
    !  else
    !     CALL SOLVESEN(Y,YPRIME,T, TOUT, INFO, RTOL, ATOL, IDID,&
    ! &        THERMO, IMPSPECIES, NUMITER, TARGETCONC)
    !  end if
     !gmagoon 080509: switched to use sensitivity flag to decide whether to 
     ! perform sensitivity analysis or not
      if (sensflag .eq. 1) then
         CALL SOLVESEN(Y,YPRIME,T, TOUT, INFO, RTOL, ATOL, IDID,&
     &        THERMO, IMPSPECIES, NUMITER, TARGETCONC)
      else
         CALL SOLVEODE(Y, YPRIME, T, TOUT, INFO, RTOL, ATOL, IDID, &
     &        THERMO, IMPSPECIES, TARGETCONC(1),AUTOFLAG, &
     &     THRESH,CORETHRESH,ESPECIES,EREACTIONSIZE,NEREAC,NEPROD,IDEREAC, &
     &     IDEPROD,KVEC)
      end if

! 6/26/08 gmagoon: deallocate memory from allocatable arrays
! (may or may not be useful)
      IF (AUTOFLAG .EQ. 1) THEN
      	DEALLOCATE (NEREAC, NEPROD, KVEC, &
 	 &	IDEREAC, IDEPROD)
      END IF
         
      END PROGRAM CALL_DASPKAUTO



!*******************************************************************

      SUBROUTINE PSOL()
      END SUBROUTINE PSOL

!*******************************************************************
      SUBROUTINE G_RES()
      END SUBROUTINE G_RES

!*******************************************************************


      SUBROUTINE SOLVEODE(Y, YPRIME, Time, TOUT, INFO, RTOL, ATOL,&
     &     IDID, THERMO, IMPSPECIES, TARGETCONC, AUTOFLAG, &
     &     THRESH,CORETHRESH,ESPECIES,EREACTIONSIZE,NEREAC,NEPROD,IDEREAC, &
     &     IDEPROD,KVEC)

      IMPLICIT NONE
      
!     INITIALIZE VARIABLES IN COMMON BLOC
      INCLUDE 'reaction.fh'


      INTEGER  INFO(30),LIW,LRW, IWORK(41+2*NSTATE), SENPAR(1), IPAR(1),&
     &     IDID, ires, I, J, iter, IMPSPECIES, AUTOFLAG, ESPECIES, &
     &      EREACTIONSIZE, EDGEFLAG, ITER_OUTPT, &
     &      IWORK_OUTPT(41 + 2*NSTATE), EDGEFLAGCORE, CORESTORED
      INTEGER NEREAC(EREACTIONSIZE),NEPROD(EREACTIONSIZE)
      INTEGER IDEREAC(EREACTIONSIZE,3), IDEPROD(EREACTIONSIZE,3)
      DOUBLE PRECISION KVEC(EREACTIONSIZE),MAXRATIO(ESPECIES),HIGHESTRATIO

      DOUBLE PRECISION Y(NEQ), YPRIME(NEQ), Time, TOUT, RTOL, ATOL, & 
    &     RWORK(51+9*NEQ+NSTATE**2), THERMO(SPMAX), TFINAL, & 
    &     TARGETCONC, THRESH, CORETHRESH, Y_OUTPT(NSTATE), &
    &     YPRIME_OUTPT(NSTATE), TIME_OUTPT,  &
    &     RWORK_OUTPT(51 + 9*NEQ + NSTATE**2)

      INTEGER PRUNEVEC(ESPECIES)

      DOUBLE PRECISION RPAR(REACTIONSIZE+THIRDBODYREACTIONSIZE+ &
     &     TROEREACTIONSIZE+LINDEREACTIONSIZE+NSTATE-1), del(neq), CJ, &
     &     TOTALREACTIONFLUX(REACTIONSIZE+THIRDBODYREACTIONSIZE+ &
     &     TROEREACTIONSIZE+LINDEREACTIONSIZE), &
     &     CURRENTREACTIONFLUX(REACTIONSIZE+THIRDBODYREACTIONSIZE+ &
     &     TROEREACTIONSIZE+LINDEREACTIONSIZE), &
     &     PREVREACTIONFLUX(REACTIONSIZE+THIRDBODYREACTIONSIZE+&
     &     TROEREACTIONSIZE+LINDEREACTIONSIZE), PREVTIME, &
     &     TOTALREACTIONFLUX_OUTPT(REACTIONSIZE+THIRDBODYREACTIONSIZE+ &
     &     TROEREACTIONSIZE+LINDEREACTIONSIZE)


      EXTERNAL RES, JAC, PSOL, G_RES


      LIW = 41 + NSTATE + NSTATE
      LRW = 51 + 9*NEQ + NSTATE**2

 101  Format(E24.15)
!     READ RWORK AND IWORK
      IF (time .NE. 0) THEN
         OPEN(UNIT=13, FILE = 'RWORK.DAT', form = 'unformatted')
         read(13) lrw
         do i=1, lrw
            READ(13) RWORK(I)
         end do
         close(13)

         OPEN(UNIT=14, FILE = 'IWORK.DAT', form = 'unformatted')
         read(14) liw
         do i=1,liw
            READ(14) IWORK(I)
         end do
         close(14)
      Else
         do i=1,lrw
            rwork(i) = 0
         end do
         do i=1, liw
            iwork(i) = 0
         end do
      end if



!     ******* INITIALIZE THE REAL ARRAY
      DO I=0,REACTIONSIZE-1
         RPAR(I+1) = REACTIONRATEARRAY(5*I+1)
      END DO

      DO I=0,THIRDBODYREACTIONSIZE-1
         RPAR(REACTIONSIZE+I+1) = THIRDBODYREACTIONRATEARRAY(16*I+1)
      END DO

      DO I=0,TROEREACTIONSIZE-1
         RPAR(REACTIONSIZE+THIRDBODYREACTIONSIZE+I+1) = &
     &        TROEREACTIONRATEARRAY(21*I+1)
      END DO

	  DO I=0,LINDEREACTIONSIZE-1
	     RPAR(REACTIONSIZE+THIRDBODYREACTIONSIZE+TROEREACTIONSIZE+I+1) = &
	 &		  LINDEREACTIONRATEARRAY(17*I+1)
	  END DO

      DO I=0,NSTATE-2
         RPAR(REACTIONSIZE+THIRDBODYREACTIONSIZE+TROEREACTIONSIZE+ &
     &        LINDEREACTIONSIZE+I+1) = THERMO(I+1)
      END DO

!     **********************************

!     initialize the yprimes for state variables and sensitivities
!     of the second kind
      
      do i=1,neq
         del(i) = 0
      end do

      if (time .eq. 0) then
         if (info(20) .ne. 0) then
            ires = 1
            cj = 0.0
            CALL RES(time, Y, del, cj, yprime, ires, RPAR, ipar, senpar)
         else
            call getflux(y, yprime, rpar)
         end if
      end if

! 6/26/08 gmagoon: call EdgeFlux if AUTOFLAG is 1 to
! determine EDGEFLAG and EDGEFLAGCORE (-1 if flux threshhold has not yet been met,
! positive integer otherwise)
	EDGEFLAG = -1
	EDGEFLAGCORE = -1
	!initialize flag CORESTORED to -1 indicating that the output values (when core inclusion threshold is exceeded) have not been stored yet; this will be set to 1 once the values are stored
	CORESTORED = -1
	!initialize MAXRATIO to a vector of all zeroes;
	!this variable will track the maxiumum
	!of rate(i)/Rchar for each edge species i, with respect to time
	MAXRATIO = 0
	!initialize HIGHESTRATIO to 0; this variable (a scalar)
	!tracks the maximum of maxratio; when a new maximum is
	!attained, the present time step creating the new maximum
	!will be stored for output from the ODE solver to the Java code
	HIGHESTRATIO = 0
	IF (AUTOFLAG.eq.1) THEN
		CALL EDGEFLUX(EDGEFLAG,EDGEFLAGCORE,Y,YPRIME,THRESH,CORETHRESH,ESPECIES, &
     &     EREACTIONSIZE,NEREAC,NEPROD,IDEREAC,IDEPROD,KVEC, &
     &     MAXRATIO,PRUNEVEC,NSTATE)
		!update HIGHESTRATIO
		!note that we cannot exit the do loop once the IF statement is
		!caught because the HIGHESTRATIO will not necessarily have the highest value
		DO I=1, ESPECIES
		    IF (MAXRATIO(I) .GT. HIGHESTRATIO) THEN
			HIGHESTRATIO = MAXRATIO(I)
		    END IF
		END DO
		!when core flux threshold is exceeded for the first time, store the values for output
		!we want to use the first point of exceeding the core inclusion threshold
		!(aka error tolerance) to determine which species to add to the core
		!note: TOTALREACTIONFLUX apparently tracks integrated flux, not
		!instantaneous flux, so setting value to be zero here should be OK...
		!in any case, the result doesn't seem to be used by the Java code
		IF(EDGEFLAGCORE .GT. 0 .AND. CORESTORED .EQ. -1) THEN
		    ITER_OUTPT=0
		    TIME_OUTPT=TIME
		    Y_OUTPT = Y
		    YPRIME_OUTPT = YPRIME
		    !IWORK_OUTPT = IWORK
		    !RWORK_OUTPT = RWORK
		    TOTALREACTIONFLUX_OUTPT = 0
		    !PRUNEVEC_OUTPT = PRUNEVEC !note that we do not copy prunevec for the output value; this will be based on the final time integrated to (in case core threshold, but not termination threshold, is exceeded at t=0...it is ok to prune species that have zero flux at t=0
		    CORESTORED = 1
		END IF

	ENDIF

      iter =0
      PREVTIME = TIME
      DO I=1,REACTIONSIZE+THIRDBODYREACTIONSIZE+TROEREACTIONSIZE+ &
	 &		LINDEREACTIONSIZE
         PREVREACTIONFLUX(I) = 0
         TOTALREACTIONFLUX(I) = 0
         CURRENTREACTIONFLUX(I) = 0
      END DO
      CALL REACTION_FLUX(Y, PREVREACTIONFLUX, RPAR)

! 6/26/08 gmagoon: added criteria that edgeflag = -1 for loop to
! continue; calls to DASSL will stop once EDGEFLAG takes on a
! different value
      IF (IMPSPECIES .EQ. -1) THEN

 1       IF (Time .LE. TOUT.AND. EDGEFLAG .EQ. -1) THEN !6/11/09 gmagoon: condensing onto one line to avoid warning
            CALL DDASPK(RES, NEQ, Time, Y, YPRIME, TOUT, INFO, RTOL, &
     &           Atol,IDID,RWORK, LRW, IWORK, LIW, RPAR, IPAR, JAC,PSOL, &
     &           SENPAR, G_RES)
            
            IF (IDID .EQ. -1) THEN
               WRITE(*,*) "Warning: 500 steps were already taken, taking &
     &              another 500 steps"
               INFO(1) = 1
               GO TO 1
            else if (idid .eq. 1) then
               iter = iter+1

               CALL REACTION_FLUX(Y, CURRENTREACTIONFLUX, RPAR)
              DO I=1,REACTIONSIZE+THIRDBODYREACTIONSIZE+TROEREACTIONSIZE+LINDEREACTIONSIZE
                  TOTALREACTIONFLUX(I) = TOTALREACTIONFLUX(I) + &
     &                 (PREVREACTIONFLUX(I) + CURRENTREACTIONFLUX(I))* &
     &                 (TIME - PREVTIME)/2
                  PREVREACTIONFLUX(I) = CURRENTREACTIONFLUX(I)
               END DO
               PREVTIME = TIME
               ! 6/26/08 gmagoon: call EdgeFlux if AUTOFLAG is 1 to
               ! determine EDGEFLAG 
		IF (AUTOFLAG.eq.1) THEN
			CALL EDGEFLUX(EDGEFLAG,EDGEFLAGCORE,Y,YPRIME,THRESH,CORETHRESH,ESPECIES, &
    		 &     EREACTIONSIZE,NEREAC,NEPROD,IDEREAC,IDEPROD,KVEC, &
		 &     MAXRATIO,PRUNEVEC,NSTATE)
			!update HIGHESTRATIO
			DO I=1, ESPECIES
			    IF (MAXRATIO(I) .GT. HIGHESTRATIO) THEN
				HIGHESTRATIO = MAXRATIO(I)
			    END IF
			END DO
			!when core flux threshold is exceeded for the first time, store the values for output
			IF(EDGEFLAGCORE .GT. 0 .AND. CORESTORED .EQ. -1) THEN
			    ITER_OUTPT=ITER
			    TIME_OUTPT=TIME
			    Y_OUTPT = Y
			    YPRIME_OUTPT = YPRIME
			    !IWORK_OUTPT = IWORK
			    !RWORK_OUTPT = RWORK
			    TOTALREACTIONFLUX_OUTPT = TOTALREACTIONFLUX
			    CORESTORED=1
			END IF
		END IF
               
               go to 1
            END IF
         END IF
      ELSE
 2       IF (Y(IMPSPECIES) .GE. TARGETCONC*Y(NSTATE).AND. EDGEFLAG .EQ. -1) &
     &        THEN

            CALL DDASPK(RES, NEQ, Time, Y, YPRIME, TOUT, INFO, RTOL, &
     &           Atol,IDID,RWORK, LRW, IWORK, LIW, RPAR, IPAR, JAC,PSOL,&
     &           SENPAR, G_RES)
            
            IF (IDID .EQ. -1) THEN
               WRITE(*,*) "Warning: 500 steps were already taken, taking &
     &              another 500 steps"
               INFO(1) = 1
               GO TO 2 !gmagoon 6/11/09: fixed to go to 2 (was originally 1,
  ! but I think this was incorrect and it was causing compiler warning)
            else if (idid .eq. 1) then
               iter = iter+1

               CALL REACTION_FLUX(Y, CURRENTREACTIONFLUX, RPAR)
              DO I=1,REACTIONSIZE+THIRDBODYREACTIONSIZE+TROEREACTIONSIZE+LINDEREACTIONSIZE
                  TOTALREACTIONFLUX(I) = TOTALREACTIONFLUX(I) + &
     &                 (PREVREACTIONFLUX(I) + CURRENTREACTIONFLUX(I))*&
     &                 (TIME - PREVTIME)/2
                  PREVREACTIONFLUX(I) = CURRENTREACTIONFLUX(I)
               END DO
!               write(*,*) time
               PREVTIME = TIME
               if (time .le. 1e6) then
                  tout = 200*time
               end if
               ! 6/26/08 gmagoon: call EdgeFlux if AUTOFLAG is 1 to
               ! determine EDGEFLAG 
		IF (AUTOFLAG.eq.1) THEN
			CALL EDGEFLUX(EDGEFLAG,EDGEFLAGCORE,Y,YPRIME,THRESH,CORETHRESH,ESPECIES, &
    		 &     EREACTIONSIZE,NEREAC,NEPROD,IDEREAC,IDEPROD,KVEC, &
		 &     MAXRATIO,PRUNEVEC,NSTATE)
		 	!update HIGHESTRATIO
			DO I=1, ESPECIES
			    IF (MAXRATIO(I) .GT. HIGHESTRATIO) THEN
				HIGHESTRATIO = MAXRATIO(I)
			    END IF
			END DO
			!when core flux threshold is exceeded for the first time, store the values for output
			IF(EDGEFLAGCORE .GT. 0 .AND. CORESTORED .EQ. -1) THEN
			    ITER_OUTPT=ITER
			    TIME_OUTPT=TIME
			    Y_OUTPT = Y
			    YPRIME_OUTPT = YPRIME
			    !IWORK_OUTPT = IWORK
			    !RWORK_OUTPT = RWORK
			    TOTALREACTIONFLUX_OUTPT = TOTALREACTIONFLUX
			    CORESTORED=1
			END IF
		END IF
               go to 2
            END IF
         END IF
      END IF

      !if we are not using AUTOFLAG, we just want to output the final times;
      !another case where we want to use the final time is if the highest ratio
      !never reaches CORETHRESH...if this is the case, we are "converged",
      !as we have reached the target time/conversion without exceeding the
      !core inclusion threshhold, and we therefore want to report the
      !final state back to the Java code so it will recognize that the
      !target time/conversion has been reached
      !in either case, CORESTORED should be -1
      IF (CORESTORED .EQ. -1) THEN
	ITER_OUTPT=ITER
	TIME_OUTPT=TIME
	Y_OUTPT = Y
	YPRIME_OUTPT = YPRIME
	!IWORK_OUTPT = IWORK
	!RWORK_OUTPT = RWORK
	TOTALREACTIONFLUX_OUTPT = TOTALREACTIONFLUX
      END IF
      IWORK_OUTPT = IWORK
      RWORK_OUTPT = RWORK

      OPEN(UNIT=15, FILE='SolverOutput.dat')

 100  Format(E24.15)
      write(15,*) ITER_OUTPT
      WRITE(15,*) (NSTATE-1)*(nparam+1)      
      WRITE(15,*) TIME_OUTPT
      do i=0, nparam
         DO j=1,NSTATE
            WRITE(15,*) Y_OUTPT(I*nstate+j)/Y_OUTPT(NSTATE)
         END DO
      end do
      
   ! 6/26/08 gmagoon: these values are read by RMG as flux; corrected to
   ! include volume changing effects (dV/dt) using quotient rule
   ! note: it appears that we only need to go to NSTATE-1, but extra "1.0" is taken into
   ! account here (and above) with extra lines in RMG's JDASPK: line = br.readLine();
      do i=0, nparam
         DO j=1,NSTATE
            !WRITE(15,*) Yprime(I*nstate+j)/Y(NSTATE)
            WRITE(15,*) (Y_OUTPT(NSTATE)*YPRIME_OUTPT(I*nstate+j)- &
	     &           Y_OUTPT(I*nstate+j)*YPRIME_OUTPT(NSTATE)) &
	     &           /(Y_OUTPT(NSTATE)**2)
         END DO
      end do

      DO I=1,REACTIONSIZE+THIRDBODYREACTIONSIZE+TROEREACTIONSIZE+LINDEREACTIONSIZE
         WRITE(15,*) TOTALREACTIONFLUX_OUTPT(I)
      END DO
      

      DO I=1,30
         WRITE(15,*) INFO(I) !gmagoon 032510: I am assuming this is not modified by ODE solver and that we don't need to keep track of intermediate values this may take on (like we do with Y, YPRIME, etc.)
      END DO

      write(15,100) Y_OUTPT(NSTATE)

      !for autoflag cases, display the edgeflag (if < 0, the ODE solver did not reach the target time/conversion),
      !, the final time integrated to,
      !along with pruning information
      IF (AUTOFLAG .EQ. 1) THEN
	  WRITE(15,*) EDGEFLAG
	  WRITE(15,*) TIME
	  DO I=1, ESPECIES
	      WRITE(15,*) PRUNEVEC(I)
	      WRITE(15,*) MAXRATIO(I)
	  END DO
      END IF

      CLOSE(15)

!gmagoon 032510: this information should only be used in non-auto cases,
!but I have modified to use OUTPT variables (which should be the same as
!"regular" variables in non-auto cases) for consistency
!     write formatted restart data to the file
      open (unit=16, file='variables.dat', form='unformatted')
      write(16) idid
      write(16) TIME_OUTPT
      write(16) nstate
      write(16) nparam
      write(16) neq
      do i=1, neq
         write(16) Y_OUTPT(i)
         write(16) YPRIME_OUTPT(i)
      end do
      write(16) rtol, atol
      do i=1,30
         write(16) info(i)
      end do

      write(16) reactionsize
      do i=1,9*reactionsize
         write(16) reactionarray(i)
      end do
      do i=1,5*reactionsize
         write(16) reactionratearray(i)
      end do

      write(16) thirdbodyreactionsize
      do i=1,20*thirdbodyreactionsize
         write(16) thirdbodyreactionarray(i)
      end do
      do i=1,16*thirdbodyreactionsize
         write(16) thirdbodyreactionratearray(i)
      end do

      write(16) troereactionsize
      do i=1,21*troereactionsize
         write(16) troereactionarray(i)
      end do
      do i=1,21*troereactionsize
         write(16) troereactionratearray(i)
      end do

      WRITE(16) LINDEREACTIONSIZE
      DO i=1,20*LINDEREACTIONSIZE
         WRITE(16) LINDEREACTIONARRAY(I)
      END DO
      DO i=1,17*LINDEREACTIONSIZE
         WRITE(16) LINDEREACTIONRATEARRAY(I)
      END DO

      do i=1,nstate-1
         write(16) thermo(i)
      end do
      write(16) temperature, pressure
     
      DO i=1,nstate
         WRITE(16) ConstantConcentration(i)
      END DO

      close(16)

      OPEN(UNIT=17, FILE='RWORK.DAT', FORM='unformatted')
      write(17) lrw
      do i=1, lrw
         WRITE(17) RWORK_OUTPT(I)
      end do
      CLOSE(17)

      OPEN(UNIT=18, FILE='IWORK.DAT', FORM='unformatted')
      write(18) liw
      do i=1,liw
         WRITE(18) IWORK_OUTPT(I)
      end do
      CLOSE(18)


      
      if (idid .eq. 1 .or. idid .eq. 2 .or. idid .eq. 3 .or. idid .eq. 0) then
         WRITE(*,*) "******ODESOLVER SUCCESSFUL: IDID=",idid
      else
         WRITE(*,*) "******ODESOLVER FAILED : IDID=", idid
      end if

      END SUBROUTINE SOLVEODE
            
!     *****************************************************************


      SUBROUTINE SOLVESEN(Y, YPRIME, Time, TOUT ,INFO, RTOL, ATOL,&
     &     IDID, THERMO, IMPSPECIES, NUMITER, TARGETCONC)

      IMPLICIT NONE
      
!     INITIALIZE VARIABLES IN COMMON BLOC
      INCLUDE 'reaction.fh'


      INTEGER  INFO(30),LIW,LRW, IWORK(41+2*NSTATE), SENPAR(1), IPAR(1), &
     &     IDID, ires, I, J, iter, NSTEPS, K, IMPSPECIES, NUMITER

      DOUBLE PRECISION Y(NEQ), YPRIME(NEQ), Time, TOUT, RTOL, ATOL &
     &     , RWORK(51+9*NEQ+NSTATE**2), THERMO(SPMAX), TSTEP, &
     &     TARGETCONC(numiter)

      DOUBLE PRECISION RPAR(REACTIONSIZE+THIRDBODYREACTIONSIZE+ &
     &     TROEREACTIONSIZE+LINDEREACTIONSIZE+NSTATE-1), del(neq), CJ, &
     &     TOTALREACTIONFLUX(REACTIONSIZE+THIRDBODYREACTIONSIZE+&
     &     TROEREACTIONSIZE+LINDEREACTIONSIZE), &
     &     CURRENTREACTIONFLUX(REACTIONSIZE+THIRDBODYREACTIONSIZE+ &
     &     TROEREACTIONSIZE+LINDEREACTIONSIZE), &
     &     PREVREACTIONFLUX(REACTIONSIZE+THIRDBODYREACTIONSIZE+ &
     &     TROEREACTIONSIZE+LINDEREACTIONSIZE), PREVTIME, TSTEPS(NUMITER)



      EXTERNAL RES, JAC, PSOL, G_RES


      LIW = 41 + NSTATE + NSTATE
      LRW = 51 + 9*NEQ + NSTATE**2

 101  Format(E24.15)

      do i=1,lrw
         rwork(i) = 0
      end do
      do i=1, liw
         iwork(i) = 0
      end do
      

!     ******* INITIALIZE THE REAL ARRAY
      DO I=0,REACTIONSIZE-1
         RPAR(I+1) = REACTIONRATEARRAY(5*I+1)
      END DO

      DO I=0,THIRDBODYREACTIONSIZE-1
         RPAR(REACTIONSIZE+I+1) = THIRDBODYREACTIONRATEARRAY(16*I+1)
      END DO

      DO I=0,TROEREACTIONSIZE-1
         RPAR(REACTIONSIZE+THIRDBODYREACTIONSIZE+I+1) = &
     &        TROEREACTIONRATEARRAY(21*I+1)
      END DO

	  DO I=0,LINDEREACTIONSIZE-1
	     RPAR(REACTIONSIZE+THIRDBODYREACTIONSIZE+TROEREACTIONSIZE+I+1) = &
	 &		  LINDEREACTIONRATEARRAY(17*I+1)
	  END DO

      DO I=0,NSTATE-2
         RPAR(REACTIONSIZE+THIRDBODYREACTIONSIZE+TROEREACTIONSIZE+LINDEREACTIONSIZE+i+1) = &
     &        THERMO(I+1)
      END DO

!     **********************************

!     initialize the yprimes for state variables and sensitivities
!     of the second kind
      
      do i=1,neq
         del(i) = 0
      end do

      ires = 1;
      cj = 0.0
      CALL RES(time, Y, del, cj, yprime, ires, RPAR, ipar, senpar)
      
      OPEN(UNIT=15, FILE='SolverOutput.dat')

      TSTEP = TOUT

      IF (IMPSPECIES .EQ. -1) THEN
         DO I=1,NUMITER
            TSTEPS(I) = TARGETCONC(I)
            TARGETCONC(I) = 0
         END DO
      ELSE
         DO I=1,NUMITER
            TSTEPS(I) = 1E6
         END DO
      END IF

      do k=1,NUMITER
         PREVTIME = TIME
         DO I=1,REACTIONSIZE+THIRDBODYREACTIONSIZE+TROEREACTIONSIZE+LINDEREACTIONSIZE
            PREVREACTIONFLUX(I) = 0
            TOTALREACTIONFLUX(I) = 0
            CURRENTREACTIONFLUX(I) = 0
         END DO
         call reaction_flux(y, prevreactionflux, rpar)
         
 1       IF (Time .LE. TSTEPS(K) .AND. Y(IMPSPECIES) .GE. TARGETCONC(K) &
     &        *Y(NSTATE)) THEN
            CALL DDASPK(RES, NEQ, Time, Y, YPRIME, Tout,INFO, RTOL, &
     &           Atol, IDID, RWORK, LRW, IWORK, LIW, RPAR, IPAR, JAC, &
     &           PSOL, SENPAR, G_RES)
            
            IF (IDID .EQ. -1) THEN
               WRITE(*,*) "Warning: 500 steps were already taken, taking &
     &              another 500 steps"
               INFO(1) = 1
               GO TO 1
            else if (idid .eq. 1) then
               iter = iter+1
               call reaction_flux(y, currentreactionflux, rpar)
              DO I=1,REACTIONSIZE+THIRDBODYREACTIONSIZE+TROEREACTIONSIZE+LINDEREACTIONSIZE
                  TOTALREACTIONFLUX(I) = TOTALREACTIONFLUX(I) + &
     &                 (PREVREACTIONFLUX(I) + CURRENTREACTIONFLUX(I))*&
     &                 (TIME - PREVTIME)/2
                  PREVREACTIONFLUX(I) = CURRENTREACTIONFLUX(I)
               END DO
               PREVTIME = TIME
               go to 1
            END IF
         END IF
         
         WRITE(*,*) "STEP: ", K, " OF ", NUMITER," DONE"
         WRITE(15,*) (NSTATE-1)*(nparam+1)      
         WRITE(15,*) TIME
         do i=0, nparam
            DO j=1,NSTATE
               !for concentrations, normalize by volume
               IF (i .EQ. 0) THEN
                WRITE(15,*) Y(I*nstate+j)/Y(NSTATE)
               !for sensitivity coefficients, we must use quotient rule:
               !Zcj,ki = dcj/dki=d(nj/V)/dki=(V*dnj/dki-nj*dV/dki)/V^2
               ELSE
                WRITE(15,*) (Y(NSTATE)*Y(I*nstate+j)-Y(j)*Y(I*nstate+NSTATE))/(Y(NSTATE)**2)
               END IF
            END DO
         end do
         
         ! 6/26/08 gmagoon: these values are read by RMG as flux; corrected to
         ! include volume changing effects (dV/dt) using quotient rule
         ! 12/21/09 gmagoon: update: "fluxes" for sensitivity coefficients may not be computed correctly, but I don't believe they are used
         do i=0, nparam
            DO j=1,NSTATE
               !WRITE(15,*) Yprime(I*nstate+j)/Y(NSTATE)
               WRITE(15,*) (Y(NSTATE)*YPRIME(I*nstate+j)-Y(I*nstate+j)*YPRIME(NSTATE))/(Y(NSTATE)**2)
            END DO
         end do

         DO I=1,REACTIONSIZE+THIRDBODYREACTIONSIZE+TROEREACTIONSIZE+LINDEREACTIONSIZE
            WRITE(15,*) TOTALREACTIONFLUX(I)
         END DO         

      end do   
      CLOSE(15)

      if (idid .eq. 1 .or. idid .eq. 2 .or. idid .eq. 3) then
         WRITE(*,*) "******ODESOLVER SUCCESSFUL: IDID=",idid
      else
         WRITE(*,*) "******ODESOLVER FAILED : IDID=", idid
      end if

      END SUBROUTINE SOLVESEN
            
! 6/26/08 gmagoon: adding subroutines EDGEFLUX and RCHAR

! EDGEFLUX determines whether core or termination threshhold flux has been
! reached by any species; if so, EDGEFLAGCORE and EDGEFLAG, respectively, are set to a
! positive value corresponding to the offending species ID;
! otherwise, these values remain at whatever they were when passed to the function
! EDGEFLUX also calls RCHAR
! when edgeflux exceeds core threshold, edgeflagcore is set to the offending species ID
! if edgeflux also exceeds the termination threshold flux, edgeflag is set to the offending species ID

SUBROUTINE EDGEFLUX(EDGEFLAG,EDGEFLAGCORE,Y,YPRIME,THRESH,CORETHRESH,ESPECIES, &
    		 &     EREACTIONSIZE,NEREAC,NEPROD,IDEREAC,IDEPROD, &
		&      KVEC, MAXRATIO,PRUNEVEC,NSTATE)
	IMPLICIT NONE
	INTEGER EDGEFLAG, EDGEFLAGCORE, ESPECIES, EREACTIONSIZE, NSTATE, I, J, K
	DOUBLE PRECISION THRESH,CORETHRESH,FLUXRC, RFLUX, Y(NSTATE), YPRIME(NSTATE)
	DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: RATE
	INTEGER NEREAC(EREACTIONSIZE),NEPROD(EREACTIONSIZE)
        INTEGER IDEREAC(EREACTIONSIZE,3), IDEPROD(EREACTIONSIZE,3)
	INTEGER PRUNEVEC(ESPECIES)
        DOUBLE PRECISION KVEC(EREACTIONSIZE),RATIO, MAXRATIO(ESPECIES)
	
	ALLOCATE(RATE(ESPECIES))
	! initialize array to have all zeroes
	RATE=0
	! it is assumed edgeflag = -1 at this point, since
      ! this should be the only way to access this subroutine

      ! calculate characteristic flux, FLUXRC
	CALL RCHAR(FLUXRC, Y, YPRIME, NSTATE)

	!initialize PRUNEVEC; at first, assume everything is prunable
	PRUNEVEC=1

	! calculate the vector of fluxes for each species
	! by summing contributions from each reaction
	ILOOP: DO I=1, EREACTIONSIZE
	      ! calculate reaction flux by multiplying k
		! by concentration(s)
		RFLUX = KVEC(I)
		KLOOP: DO K=1, NEREAC(I)
			RFLUX = RFLUX*Y(IDEREAC(I,K))/Y(NSTATE)
		END DO KLOOP
		! loop over reaction products, adding RFLUX
		JLOOP: DO J=1, NEPROD(I)
			RATE(IDEPROD(I,J))=RATE(IDEPROD(I,J))+ RFLUX
			!if RFLUX is zero, then presumably one of the
			!reactant concentrations is zero (it is assumed
			!that k != 0), and therefore, we don't want to
			!prune the species yet
			IF(RFLUX .EQ. 0.0) THEN
			    PRUNEVEC(IDEPROD(I,J)) = 0
			END IF
		END DO JLOOP
	END DO ILOOP


	! check if any of the edge species fluxes exceed
	! the threshhold and update the MAXRATIO vector
	FLOOP: DO I=1, ESPECIES
		RATIO = RATE(I)/FLUXRC
		!if we reach the error tolerance (aka core tolerance),
		!set EDGEFLAGCORE to the species ID
		!if we also reach the termination tolerance, set EDGEFLAG to the species ID
		IF (RATIO .GE. CORETHRESH) THEN
		    EDGEFLAGCORE = I
		    IF(RATIO .GE. THRESH) THEN
			EDGEFLAG = I
		    END IF
		END IF
		!update the MAXRATIO vector
		IF(RATIO .GT. MAXRATIO(I)) THEN
		    MAXRATIO(I) = RATIO
		END IF
	END DO FLOOP
	DEALLOCATE(RATE)
END SUBROUTINE EDGEFLUX

! RCHAR calculates the characteristic flux (i.e. the 
! L2 norm of the core species fluxes), which it
! returns through the variable FLUXRC); the subroutine
! is called by EDGEFLUX; Y and YPRIME contain the 
! state variables and state variable derivatives,
! respectively; state vars. include the number of moles of
! each species, with the volume as the final state
! variable; NSTATE indicates the number of state
! variables (equal to the number of core species
! plus one)
SUBROUTINE RCHAR(FLUXRC, Y, YPRIME, NSTATE)
       IMPLICIT NONE 
	INTEGER NSTATE, I
	DOUBLE PRECISION SSF, FLUXRC, Y(NSTATE), YPRIME(NSTATE)
	! compute the sum of squared fluxes
	SSF=0.0

        !6/26/08 gmagoon: note that we are implicitly treating only indexing
        !first NSTATE elements of array (i.e. nparam=0) since we do not care
        !about sensitivity here
	DO I=1, NSTATE-1
	! add (dCi/dt)^2 with dCi/dt computed based on
        ! quotient rule...dCi/dt=d(Ni/V)/dt
	! =(V*dNi/dt-Ni*dV/dt)/V^2
		SSF=SSF+((Y(NSTATE)*YPRIME(I)-Y(I)*YPRIME(NSTATE)) &
       &                        /(Y(NSTATE)**2))**2
	END DO
	! compute the square root, corresponding to
      ! L2 norm
	FLUXRC = SSF**0.5
END SUBROUTINE RCHAR


