C
C     Read input from text file. (PEY 9March04)
C
      SUBROUTINE REDKEY (KK,KSYM,LIN,LOUT,REAC,KERR,ICASE)

      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C
      DIMENSION REAC(KK), RTOL(1), ATOL(1), VALUE(1)
C
      CHARACTER KEY*4, KSYM(*)*(*), LINE*80, CKCHUP*4
      LOGICAL KERR, IERR

      COMMON /INPUTS/ RTOL,ATOL,TEMP,PRES,TZERO,TSTOP,SPD,CADZERO, 
     >                RC,RLA,VC

C     defaults (if there are any)

C     read next input line
 90   CONTINUE

      LINE = ' '
      IERR = .FALSE.
      READ  (LIN,  '(A)',END=100) LINE
      WRITE (LOUT, '(A)') LINE
      CALL CKDTAB (LINE)
      KEY = ' '
      KEY = CKCHUP(LINE,4)
      LINE(1:4) = ' '

C     is this keyword commented out
      IF (KEY(1:1) .EQ. '/' .OR. KEY(1:1) .EQ. '!') GO TO 90
	
C	problem choice (ic engine model)
	IF (KEY .EQ. 'ICEN') THEN
		ICASE = 1
C	problem choice (const t,p model)
	ELSEIF (KEY .EQ. 'CNTP') THEN
		ICASE = 2     
C     relative tolerance
      ELSEIF (KEY .EQ. 'RTOL') THEN
         CALL CKXNUM (LINE, 1, LOUT, NVAL, RTOL, IERR)
C     absolute tolerance
      ELSEIF (KEY .EQ. 'ATOL') THEN
         CALL CKXNUM (LINE, 1, LOUT, NVAL, ATOL, IERR)
C     temperature (K)
      ELSEIF (KEY .EQ. 'TEMP') THEN
         CALL CKXNUM (LINE, 1, LOUT, NVAL, TEMP, IERR)
C     pressure (Pa)
      ELSEIF (KEY .EQ. 'PRES') THEN
         CALL CKXNUM (LINE, 1, LOUT, NVAL, PRES, IERR)
C     species concentrations (mol/cm3)
      ELSEIF (KEY .EQ. 'REAC') THEN
         CALL CKSNUM (LINE, 1, LOUT, KSYM, KK, KSPEC, NVAL,
     >        VALUE, IERR)
         IF (.NOT. IERR) REAC(KSPEC) = VALUE(1)
C     initial time
      ELSEIF (KEY .EQ. 'TIM0') THEN
         CALL CKXNUM (LINE, 1, LOUT, NVAL, TZERO, IERR)
C     final time
      ELSEIF (KEY .EQ. 'TIMF') THEN
         CALL CKXNUM (LINE, 1, LOUT, NVAL, TSTOP, IERR)      
C     engine speed (rpm)                  
      ELSEIF (KEY .EQ. 'SPED') THEN
         CALL CKXNUM (LINE, 1, LOUT, NVAL, SPD, IERR)
C     initial CAD (degrees)
      ELSEIF (KEY .EQ. 'CADI') THEN
         CALL CKXNUM (LINE, 1, LOUT, NVAL, CADZERO, IERR)
C     compression ratio
      ELSEIF (KEY .EQ. 'COMP') THEN
         CALL CKXNUM (LINE, 1, LOUT, NVAL, RC, IERR)
C     ratio of connecting rod to crank radius
      ELSEIF (KEY .EQ. 'RLA') THEN
         CALL CKXNUM (LINE, 1, LOUT, NVAL, RLA, IERR)
C     clearance volume [m3]
      ELSEIF (KEY .EQ. 'VOLC') THEN
         CALL CKXNUM (LINE, 1, LOUT, NVAL, VC, IERR)
C     last line
      ELSEIF (KEY .EQ. 'END ') THEN
         GO TO 100
C     keyword not found
      ELSE
	   KERR = .TRUE.
         WRITE (LOUT, '(/2X,A,A,A)')
     >        ' Warning, keyword ', KEY, ' was not found.'
      ENDIF

C     error message
      IF (IERR) THEN
         KERR =.TRUE.
         WRITE (LOUT, '(2X,A,A,A,/)')
     >        ' Error reading data for keyword ', KEY,'.'
      ENDIF
      
      GO TO 90
C     
 100  CONTINUE
C     END OF INPUT

      IF (KERR) RETURN

C     Check that concentrations are consistent with T and P
      CTOT1 = PRES/(8.314D6*TEMP)

      CTOT2 = 0.0D0
      DO K=1,KK
         CTOT2 = CTOT2 + REAC(K)
      ENDDO

      DIFFC = ABS(CTOT1-CTOT2)/(0.5D0*(CTOT1+CTOT2))
      IF(DIFFC .GT. 0.1) THEN
         WRITE (LOUT, *) ' Warning: concentrations, T, P inconsistent'
         KERR =.TRUE.
      ENDIF

C     end of SUBROUTINE REDKEY
      RETURN
      END
C     
C----------------------------------------------------------------------
C     C
      SUBROUTINE CKDTAB (STRING)
C     
C     START PROLOGUE
C     
C     SUBROUTINE CKDTAB (STRING)
C     Replaces any tab character in a character string with one space.
C     
C*****precision > double
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C*****END precision > double
C*****precision > single
C     IMPLICIT REAL (A-H, O-Z), INTEGER (I-N0
C*****END precision > single
C     
      CHARACTER STRING*(*), TAB*1
C     
      TAB = CHAR(9)
 10   CONTINUE
      IND = INDEX(STRING,TAB)
      IF (IND .GT. 0) THEN
         STRING(IND:IND) = ' '
         GO TO 10
      ENDIF
C     
C     end of SUBROUTINE CKDTAB
      RETURN
      END
C     
C----------------------------------------------------------------------C
C     C
      CHARACTER*(*) FUNCTION CKCHUP(ISTR, ILEN)
C     
C     Change lowercase letters to uprrercase
C     
C     
      CHARACTER*(*) ISTR
      CHARACTER*1 LCASE(26), UCASE(26)
      DATA LCASE /'a','b','c','d','e','f','g','h','i','j','k','l','m',
     1            'n','o','p','q','r','s','t','u','v','w','x','y','z'/,
     2     UCASE /'A','B','C','D','E','F','G','H','I','J','K','L','M',
     3            'N','O','P','Q','R','S','T','U','V','W','X','Y','Z'/
C     
      CKCHUP = ' '
      CKCHUP = ISTR(1:ILEN)
      JJ = MIN (LEN(CKCHUP), LEN(ISTR), ILEN)
      DO 10 J = 1, JJ
         DO 05 N = 1,26
            IF (ISTR(J:J) .EQ. LCASE(N)) CKCHUP(J:J) = UCASE(N)
   05    CONTINUE
 10   CONTINUE
C     
C     end of FUNCTION CKCHUP
      RETURN
      END
C
C----------------------------------------------------------------------C
C
