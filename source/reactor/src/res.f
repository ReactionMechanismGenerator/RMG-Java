      SUBROUTINE RICEN (TIME, Z, ZP, DEL, IRES, RPAR, IPAR)
C	Residual function for the IC engine reactor model

      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
      DIMENSION Z(*), ZP(*), DEL(*), RPAR(*), IPAR(*)

      COMMON /IPOINT/ NIPAR,NICK,NIDSL
      COMMON /RPOINT/ NRPAR,NREAC,NWT,NU,NWDOT,NCVMS,NC,NQ,NDQ,NRCK,
     >                NRDSL
      COMMON /ICONS/  NSPC,NRXN,NVAR,NT,NY,NP
      COMMON /RCONS/  R,SPDRAD,THETAZERO,RMASS

C *** Calculate time derivatives

      CALL FUNZP (TIME, Z, DEL, RPAR, IPAR)

C *** Compute residuals

      DO I=1,NVAR-1
         DEL(I) = DEL(I) - ZP(I)
      END DO

C *** Enforce ideal gas law
      CALL CKMMWY(Z(NY),IPAR(NICK),RPAR(NRCK),WTM)
      WTM = WTM*1.0D-3                               ! [kg/mol]

      CALL VOLUME(TIME, VCYL, DVCYLDT)
      V   = VCYL/RMASS                               ! [m3/kg]

      DEL(NP) = Z(NP)*V - R/WTM*Z(NT)

      RETURN
      END
C ***
C ***
C ***
C ***
      SUBROUTINE FUNZP(TIME, Q, DQ, RPAR, IPAR)
C	Calculates derivatives for the IC engine model

      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
      DIMENSION Q(*),DQ(*),RPAR(*),IPAR(*)

      COMMON /IPOINT/ NIPAR,NICK,NIDSL
      COMMON /RPOINT/ NRPAR,NREAC,NWT,NU,NWDOT,NCVMS,NC,NQ,NDQ,NRCK,
     >                NRDSL
      COMMON /ICONS/  NSPC,NRXN,NVAR,NT,NY,NP
      COMMON /RCONS/  R,SPDRAD,THETAZERO,RMASS
      COMMON /INPUTS/ RTOL,ATOL,TEMP,PRES,TZERO,TSTOP,SPD,CADZERO, 
     >                RC,RLA,VC

      T    = Q(NT)
      P    = Q(NP)

C     Calculate cylinder volume
      CALL VOLUME(TIME,VCYL, DVCYLDT)

      V    = VCYL/RMASS              ! [m3/kg]
      DVDT = DVCYLDT/RMASS           ! [m3/kg-s]

      RHO  = 1.0D0/V*1.0D-3          ! [g/cm3]

C     Call CHEMKIN subroutines
      CALL CKYTCR (RHO, T, Q(NY), IPAR(NICK), RPAR(NRCK), RPAR(NC))
      CALL CKWC (T, RPAR(NC), IPAR(NICK), RPAR(NRCK), RPAR(NWDOT))
      CALL CKCVBS (T, Q(NY), IPAR(NICK), RPAR(NRCK), CVB)
      CALL CKUMS (T, IPAR(NICK), RPAR(NRCK), RPAR(NU))

      CVB = CVB*1.0D-4         ! convert dyne/gm  to J/kg-K

C     Form governing equations

      HSRC = 0.0
      DO K=1, NSPC
         U    = RPAR(NU    + K-1)*1.0D-4                     ! [J/kg]
         WDOT = RPAR(NWDOT + K-1)*1.0D6                      ! [mol/m3-s]
         WT   = RPAR(NWT   + K-1)*1.0D-3                     ! [kg/mol]

         DQ(NY + K-1) = V*WDOT*WT
         HSRC = HSRC + U*WDOT*WT
      END DO

      DQ(NT) = - (V/CVB)*HSRC - (P/CVB)*DVDT
      
      RETURN
      END
C ***
C ***
C ***
C ***
      SUBROUTINE VOLUME(TIME, VOL, DVOLDT)
C     Returns the extensive volume[m3] and its derivative[m3/s] given time[s] 
C	for IC engine model.

      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)

      COMMON /RCONS/  R,SPDRAD,THETAZERO,RMASS
      COMMON /INPUTS/ RTOL,ATOL,TEMP,PRES,TZERO,TSTOP,SPD,CADZERO, 
     >                RC,RLA,VC

      THETA   = THETAZERO + TIME*SPDRAD                    ! [Radians]

      VOL     = VC*(1.0 + 0.5*(RC-1.0)*                    ! [m3]
     >          (RLA + 1.0 - COS(THETA) - SQRT(RLA**2 - SIN(THETA)**2)))

      DVOLDT  = VC*SPDRAD*0.5*(RC-1.0)*SIN(THETA)*         ! [m3/s]
     >              (1.0 + COS(THETA)/SQRT(RLA**2 - SIN(THETA)**2))

      RETURN
      END
C ***
C ***
C ***
C ***
      SUBROUTINE RCNTP(TIME, Z, ZP, DEL, IRES, RPAR, IPAR)
C	Residual subroutine for a constant T and P reactor

      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER (I-N)
      DIMENSION Z(*), ZP(*), DEL(*), RPAR(*), IPAR(*)

      COMMON /IPOINT/ NIPAR,NICK,NIDSL
      COMMON /RPOINT/ NRPAR,NREAC,NWT,NU,NWDOT,NCVMS,NC,NQ,NDQ,NRCK,
     >                NRDSL
      COMMON /ICONS/  NSPC,NRXN,NVAR,NT,NY,NP
      COMMON /RCONS/  R,SPDRAD,THETAZERO,RMASS
      COMMON /INPUTS/ RTOL,ATOL,TEMP,PRES,TZERO,TSTOP,SPD,CADZERO, 
     >                RC,RLA,VC

C	Constant temperature and pressure
      T = TEMP
      P = PRES
      PCGS = 10.0D0*P ! convert Pa to dyne/cm2

C     Call CHEMKIN subroutines
      CALL CKWYP (PCGS, T, Z, IPAR(NICK), RPAR(NRCK), RPAR(NWDOT))
	CALL CKRHOY (PCGS, T, Z, IPAR(NICK), RPAR(NRCK), RHO)
	VSP = 1./RHO

C     Form species equations
      DO K=1, NSPC
         WDOT = RPAR(NWDOT + K-1)
         WT   = RPAR(NWT   + K-1)
         DEL(K) = VSP*WDOT*WT - ZP(K)
      END DO
  
      RETURN
      END
