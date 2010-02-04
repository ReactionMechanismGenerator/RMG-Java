C
C     Driver routine for REACTOR code
C
      PROGRAM DRIVER

      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
      PARAMETER (LR=2000000, LI=500000, LC=5000)
      PARAMETER (LRDSL=500000, LIDSL=200000)
      PARAMETER (LIN=9,LPAR=10,LDAT=11,LINKCK=12,LOUT=6)
      DIMENSION RWORK(LR), IWORK(LI) !, TEMP(10), PRES(10)
      CHARACTER CWORK(LC)*16
	LOGICAL XML

      COMMON /ICONS/  NSPC,NRXN,NVAR,NT,NY,NP
      COMMON /IPOINT/ NIPAR,NICK,NIDSL
      COMMON /RPOINT/ NRPAR,NREAC,NWT,NU,NWDOT,NCVMS,NC,NQ,NDQ,NRCK,
     >                NRDSL

      DATA RWORK/LR*0.0D0/, IWORK/LI*0/, XML/.TRUE./

C     Open input and output files
	IF(XML) THEN
		OPEN (LIN,    STATUS='UNKNOWN',FILE='reactorInput.xml')
		OPEN (LDAT,   STATUS='UNKNOWN',FILE='reactorOutput.xml')
	ELSE
		OPEN (LIN,    STATUS='UNKNOWN',FILE='reactor.inp')
	    OPEN (LDAT,   STATUS='UNKNOWN',FILE='reactor.out')
	ENDIF

      OPEN (LINKCK, FORM='UNFORMATTED',STATUS='UNKNOWN',FILE='chem.bin')
      OPEN (LPAR,   STATUS='UNKNOWN',FILE='parameter.inp')
      OPEN (LOUT,   STATUS='UNKNOWN',FILE='diagnostics.out')

C     Get size of workspaces
      CALL CKLEN(LINKCK, LOUT, LENI, LENR, LENC)

C     Determine mechanism size
      IF (LENI.LT.LI .AND. LENR.LT.LR .AND. LENC.LT.LC) THEN
         CALL CKINIT(LI, LR, LC, LINKCK, LOUT, IWORK, RWORK, CWORK)
         CALL CKINDX(IWORK, RWORK, NELM, NSPC, NRXN, NFIT)
      END IF

C     Number of states and number of equations
      NVAR = NSPC + 2
      NEQ  = NVAR

C     Pointers to Q,DQ arrays
      NT    = 1               ! Pointer to temperature variable
      NY    = NT + 1          ! Pointer to mass fraction variable
      NP    = NY + NSPC       ! Pointer to pressure variable

C     Apportion integer workspace
      NNPAR = 1               ! Number of integer parameters (NRXN)
      NIPAR = NNPAR + 1       ! Integer parameters (point to position of A factors)
      NICK  = NIPAR + NRXN    ! Pointer to integer chemkin workspace
      NIDSL = NICK  + LENI    ! DSL48S integer work space
      LITOT = NIDSL + LIDSL-1 ! Length of IWORK

C     Apportion real workspace
      NRPAR = 1               ! Real parameters (A factors)
      NREAC = NRPAR + NRXN    ! Pointer to REAC array
      NWT   = NREAC + NSPC    ! Pointer to molecular weights
      NU    = NWT   + NSPC    ! Pointer to internal energies
      NWDOT = NU    + NSPC    ! Pointer to molar production rates
      NCVMS = NWDOT + NSPC    ! Pointer to Cv's
      NC    = NCVMS + NSPC    ! Pointer to concentrations
      NQ    = NC    + NSPC    ! Pointer to states and sensitivities
      NDQ   = NQ    + NEQ     ! Pointer to derivatives of states and sens
      NRCK  = NDQ   + NEQ     ! Pointer to real chemkin workspace
      NRDSL = NRCK  + LENR    ! DSL48S real work space
      LRTOT = NRDSL + LRDSL-1 ! Length of RWORK

C     Apportion character workspace
      NCCK  = 1              ! Pointer to character chemkin workspace
      NKSYM = NCCK  + LENC   ! Pointer to species names
      LCTOT = NKSYM + NSPC-1 ! Length of CWORK

C     Check if the alloted workspace is large enough
      IF (LITOT.GT.LI .OR. LRTOT.GT.LR .OR. LCTOT.GT.LC) THEN

         WRITE (LOUT, 7000) LI, LITOT, LR, LRTOT, LC, LCTOT
         WRITE (LOUT, '(A)') 'ERROR: NOT ENOUGH WORKSPACE PROVIDED'

 7000    FORMAT (/,'                WORKING SPACE REQUIREMENTS',
     1           /,'                 PROVIDED        REQUIRED ',
     2           /,' INTEGER  ' , 2I15,
     3           /,' REAL     ' , 2I15,
     4           /,' CHARACTER' , 2I15,/)

	   CALL ERRORXML(LDAT,'ERROR: NOT ENOUGH WORKSPACE PROVIDED')
         STOP
      END IF

C     Initialize CHEMKIN workspace
      DO I=1,LR
         RWORK(I) = 0.0D0
      END DO

      DO I=1,LI
         IWORK(I) = 0
      END DO

      REWIND(LINKCK)
      CALL CKINIT(LI, LR, LC, LINKCK, LOUT, IWORK(NICK), RWORK(NRCK), 
     >           CWORK(NCCK))
      CLOSE(LINKCK)

c	begin test
c	do i=1,10
c		temp(i) = 300.0 + float(i-1)*20.0
c		pres(i) = 1.01325e6*2 + float(i-1)*1.01325e6
c	enddo
c	rwork(nc) = 1.0d0
c	do i=1,10
c		do j=1,10
c			call ckwyp(pres(i), temp(j), rwork(nc), iwork(nick), 
c     >                   rwork(nrck), rwork(nwdot))
c	    enddo
c      enddo
c	stop
c	end test

C     Call reactor subroutine
      CALL REACTOR(NEQ,LIN,LPAR,LDAT,LOUT,RWORK,IWORK,CWORK,
     >             RWORK(NQ),RWORK(NDQ),CWORK(NKSYM),RWORK(NREAC),
     >             LRDSL,LIDSL,XML)

      CLOSE(LOUT)
	CLOSE(LDAT)

      STOP
      END
