C     
C Perfectly-stirred reactor model (either IC engine or Constant T,P)
C Uses DASSL to integrate DAE system
C
C Paul Yelvington, MIT (5 April 04)
C  
C System of Equations: 
C
C	For IC engine model (adiabatic, variable-volume reactor),
C
C     dY/dt = V*WDOT*WT , K=1,NSPC
C     Cv*dT/dt = -P*dV/dt - V*SUM(E*WDOT*WT)
C     RHO*R/WTM*T = P
C     VCYL = VCYL(t),  V = VCYL(t)/M
C
C	For Constant T,P reactor,
C
C     dY/dt = V*WDOT*WT , K=1,NSPC
C
C
      SUBROUTINE REACTOR(NEQ,LIN,LPAR,LDAT,LOUT,RWORK,IWORK,
     >                       CWORK,Q,DQ,KSYM,REAC,LRDSL,LIDSL,XML)
 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
      
      DIMENSION RWORK(*), IWORK(*), Q(*), DQ(*), REAC(*)
      DIMENSION ATOL(1), RTOL(1), INFO(25)
      CHARACTER CWORK(*)*(*), KSYM(*)*(*), XERN*8
      LOGICAL KERR, XML

      EXTERNAL RICEN,RCNTP

      COMMON /IPOINT/ NIPAR,NICK,NIDSL
      COMMON /RPOINT/ NRPAR,NREAC,NWT,NU,NWDOT,NCVMS,NC,NQ,NDQ,NRCK,
     >                NRDSL
      COMMON /ICONS/  NSPC,NRXN,NVAR,NT,NY,NP
      COMMON /RCONS/  R,SPDRAD,THETAZERO,RMASS
      COMMON /INPUTS/ RTOL,ATOL,TEMP,PRES,TZERO,TSTOP,SPD,CADZERO, 
     >                RC,RLA,VC

      DATA  KERR/.FALSE./, INFO /25*0/

C     Get species names and molecular weights
      CALL CKWT(IWORK(NICK), RWORK(NRCK), RWORK(NWT))
      CALL CKSYMS(CWORK, LOUT, KSYM, KERR)
      IF (KERR) THEN
		CALL ERRORXML(LDAT,'Error reading species names.')
		STOP
	ENDIF

C     Read keyword or xml input file (reactor.inp or reactorInput.xml)
	ICASE = 0
	IF (XML) THEN
		CALL READXML (NSPC,KSYM,LIN,LOUT,REAC,KERR,ICASE)
	ELSE
		CALL REDKEY (NSPC,KSYM,LIN,LOUT,REAC,KERR,ICASE)
	ENDIF

      IF (KERR) THEN
		CALL ERRORXML(LDAT,'Error reading input file.')
		STOP
	ENDIF

      CLOSE(LIN)

C     Read keyword input for IC engine parameters (parameter.inp)
	IF(ICASE==1) THEN
		CALL REDKEY (NSPC,KSYM,LPAR,LOUT,REAC,KERR,ICASE)
		IF (KERR) THEN
			CALL ERRORXML(LDAT,'Error reading parameters file.')
			STOP
		ENDIF
		CLOSE(LPAR)
	ENDIF

C     Other parameters
      T    = TZERO                      ! Start time for simulation [s]
C      DELT = (TSTOP-TZERO)/10.0       ! Time-step for printout [s]
                 
      R    = 8.314D0                    ! Universal gas constant [m3*Pa/mol-K]
      PI   = 3.14159                    ! Pi

      SPDRAD    = SPD*2.0*PI/60.0       ! Engine speed [radian/sec]
      THETAZERO = CADZERO*2.0*PI/360.0  ! Initial crank angle [radian]

C     Parameters for DASSL
      INFO(3) = 0             ! "1" sets intermediate output mode

C *** Set up consistent initial conditions

C     Temperature, Mass Fractions, and Pressure
      Q(NT) = TEMP   ! [K]
      Q(NP) = PRES   ! [Pa]
      CALL CKCTY(REAC,IWORK(NICK),RWORK(NRCK),Q(NY)) ! convert [mol/cm3] to Yk

	IF (ICASE==1) THEN ! IC engine
C		Calculate cylinder volume and its derivative
		CALL VOLUME(T,VCYL,DVCYLDT)
     
C		Calculate initial mass in cylinder
		CALL CKMMWY(Q(NY),IWORK(NICK),RWORK(NRCK),WTM)
		WTM = WTM*1.D-3                ! [kg/mol]
		RMASS = PRES*VCYL*WTM/(R*TEMP) ! [kg]

C		Temperature and Mass Fractions Derivatives
		CALL FUNZP(T,Q,DQ,RWORK,IWORK)

C		Pressure Derivative (algebraic equation)
		DQ(NP) = Q(NP)*(DQ(NT)/Q(NT) - DVCYLDT/VCYL)

	ELSEIF (ICASE==2) THEN ! Constant T,P
		CALL RCNTP(T,Q(NY),DQ(NY),DQ(NY),IRES,RWORK,IWORK)
	ELSE
		CALL ERRORXML(LDAT,'Problem keyword not found.')
		STOP
	ENDIF

C     Print page headings
C      WRITE(LOUT,*)
C      WRITE(LOUT,7100) (KSYM(K)(:10), K=1,NSPC)

C *** MAIN LOOP

C 10   CONTINUE

C *** Print out solution
C      CAD = CADZERO + T*SPD*6.0  ! [Degrees]

C      WRITE(LOUT, 7101) T, CAD, Q(NT), Q(NP)*1.0D-5, 
C     >     (Q(NY+K-1), K=1,NSPC)

C *** Increment time and call the integrator

C      TOUT = T+DELT

 15	CONTINUE

C	IC Engine Model
	IF (ICASE==1) THEN
		CALL DDASSL(RICEN, NEQ, T, Q, DQ, TSTOP, INFO,
     >            RTOL, ATOL, IDID, RWORK(NRDSL), LRDSL, IWORK(NIDSL),
     >            LIDSL, RWORK, IWORK, JAC)
C	Constant T and P Model
	ELSEIF (ICASE==2) THEN
		CALL DDASSL(RCNTP, NSPC, T, Q(NY), DQ(NY), TSTOP, INFO,
     >            RTOL, ATOL, IDID, RWORK(NRDSL), LRDSL, IWORK(NIDSL),
     >            LIDSL, RWORK, IWORK, JAC)
	ENDIF

C     Print error message from DASSL      
      IF (IDID.LT.0) THEN
         IF (IDID.EQ.-1) THEN
            INFO(1) = 1
            GOTO 15
         ELSE
            WRITE (XERN, '(I8)') IDID
            CALL ERRORXML(LDAT,'DASSL ERROR: IDID = ' // XERN)
            STOP
         END IF
      END IF

      IF (IDID.EQ.1) GO TO 15

C     Exit condition
C      IF (TOUT.LT.TSTOP*0.999) GO TO 10

C *** Print DASSL statistics
      CALL CPUTIM(TIMER)

      WRITE(LOUT,*)
      WRITE(LOUT,*) 'CPU Time (sec):',TIMER
      WRITE(LOUT,*)
      WRITE(LOUT,*) 'Number of steps: ',IWORK(NIDSL+10)
      WRITE(LOUT,*) 'Number of function evaluations: ',IWORK(NIDSL+11)
      WRITE(LOUT,*) 'Number of Jacobian evaluations: ',IWORK(NIDSL+12)
      WRITE(LOUT,*) 'Number of converg. test failures: ',IWORK(NIDSL+14)
      WRITE(LOUT,*) 'Number of error test failures: ',IWORK(NIDSL+13)
      WRITE(LOUT,*) 'SUCCESSFULLY COMPLETED RUN.'
      WRITE(LOUT,*)

C *** Write reactorOutput.xml or reactor.out file to pass back to RMG
      PCGS = 10.0D0*Q(NP) ! convert Pa to dyne/cm2 
      CALL CKYTCP(PCGS,Q(NT),Q(NY),IWORK(NICK),RWORK(NRCK),RWORK(NC)) ! get [Xk]
      CALL CKWC(Q(NT),RWORK(NC),IWORK(NICK),RWORK(NRCK),RWORK(NWDOT)) ! get wdot

	IF (XML) THEN
		WRITE(LDAT, '(A)') '<?xml version="1.0" standalone="no"?>'
		WRITE(LDAT, '(A)') '<!DOCTYPE reactoroutput SYSTEM ' //
     >         '"./documentTypeDefinitions/reactorOutput.dtd">'
		WRITE(LDAT, '(A)') '<reactoroutput>'
		WRITE(LDAT, '(A)') '<header>'
		WRITE(LDAT, '(A)') '<title>Reactor Output File</title>'
		WRITE(LDAT, '(2A)') '<description>External-solver-generated ', 
     >       'file that returns the mixture state to RMG</description>'
		WRITE(LDAT, '(A)') '</header>'
		WRITE(LDAT, '(A)') '<returnmessage>' //
     >                    'SUCCESSFULLY COMPLETED RUN.</returnmessage>'
		WRITE(LDAT, '(A)') '<outputvalues>'
		WRITE(LDAT, 7103) '<time units="sec">',T,'</time>'
		WRITE(LDAT, '(A)') '<systemstate>'
		WRITE(LDAT, 7103) 
     >                 '<temperature units="K">',Q(NT),'</temperature>'
		WRITE(LDAT, 7103) '<pressure units="Pa">',Q(NP),'</pressure>'
		DO K=1,NSPC
			NCHARS = ILASCH(KSYM(K))
			WRITE(LDAT,7104) '<amount units="molPerCm3" speciesid="',
     >               KSYM(K)(:NCHARS), '">', RWORK(NC+K-1), '</amount>'
		ENDDO
		DO K=1,NSPC
			NCHARS = ILASCH(KSYM(K))
			WRITE(LDAT,7104) '<flux units="molPerCm3-Sec" speciesid="',
     >               KSYM(K)(:NCHARS), '">', RWORK(NWDOT+K-1), '</flux>'
		ENDDO
		WRITE(LDAT, '(A)') '</systemstate>'
		WRITE(LDAT, '(A)') '</outputvalues>'
		WRITE(LDAT, '(A)') '</reactoroutput>'
	ELSE
		WRITE(LDAT, '(A)') '// units are sec, K, Pa, mol/cm3'
		WRITE(LDAT,'(A,14X,E16.8E3)') 'time', T
		WRITE(LDAT,'(A,14X,E16.8E3)') 'temp', Q(NT)
		WRITE(LDAT,'(A,14X,E16.8E3)') 'pres', Q(NP)
		DO K=1,NSPC
			WRITE(LDAT,7102) KSYM(K)(:16), RWORK(NC+K-1)
		ENDDO
	ENDIF

C FORMATS
 7100 FORMAT (2X, 'T(sec)', 8X, 'CAD', 11X, 'TMP(K)', 8X, 
     >     'P(bar)', 4X, 1100(2X,A12))
 7101 FORMAT (1100(1X,E13.5E3))
 7102 FORMAT (A16,2X, E16.8E3)
 7103 FORMAT (A,1X,E16.8E3,1X,A)
 7104 FORMAT (3A,1X,E16.8E3,1X,A)

      RETURN
      END
C ***
C ***
C ***
C ***
	SUBROUTINE ERRORXML(LDAT,ERRMSG)
C	Write an error message in XML format for RMG

      IMPLICIT NONE
	INTEGER LDAT,NCHARS,ILASCH
	CHARACTER ERRMSG*(*)

	NCHARS = ILASCH(ERRMSG)

	WRITE(LDAT, '(A)') '<?xml version="1.0" standalone="no"?>'
		WRITE(LDAT, '(A)') '<!DOCTYPE reactoroutput SYSTEM ' //
     >         '"./documentTypeDefinitions/reactorOutput.dtd">'
	WRITE(LDAT, '(A)') '<reactoroutput>'
	WRITE(LDAT, '(A)') '<header>'
	WRITE(LDAT, '(A)') '<title>Reactor Output File</title>'
	WRITE(LDAT, '(2A)') '<description>External-solver-generated ', 
     >       'file that returns  the mixture state to RMG</description>'
	WRITE(LDAT, '(A)') '</header>'
	WRITE(LDAT, '(A)') '<returnmessage>' // ERRMSG(:NCHARS) // 
     >                   '</returnmessage>'
	WRITE(LDAT, '(A)') '</reactoroutput>'

	RETURN
	END
