C     ************************************************
      SUBROUTINE REACTION_FLUX(TEMPY, REACTIONFLUX, RPAR) 

      IMPLICIT NONE

C     INITIALIZE VARIABLES IN COMMON BLOC
      INTEGER NSTATE, REACTIONSIZE, THIRDBODYREACTIONSIZE,
     $     TROEREACTIONSIZE, LINDEREACTIONSIZE
      COMMON /SIZE/ NSTATE, REACTIONSIZE, THIRDBODYREACTIONSIZE,
     $     TROEREACTIONSIZE, LINDEREACTIONSIZE

      INTEGER REACTIONARRAY(9*100000), THIRDBODYREACTIONARRAY(20*100),
     $     TROEREACTIONARRAY(21*100), LINDEREACTIONARRAY(20*100)

      DOUBLE PRECISION REACTIONRATEARRAY(5*100000),
     $     THIRDBODYREACTIONRATEARRAY(16*100),
     $     TROEREACTIONRATEARRAY(21*100), temperature, pressure,
     $     LINDEREACTIONRATEARRAY(17*100)

      COMMON /REAC/ REACTIONRATEARRAY, ThirdbodyREACTIONRATEARRAY,
     $     TROEREACTIONRATEARRAY,temperature, pressure, 
     $     REACTIONARRAY, THIRDBODYREACTIONARRAY, TROEREACTIONARRAY,
     $     LINDEREACTIONRATEARRAY, LINDEREACTIONARRAY


      
      INTEGER I, J, RNUM, PNUM, NUMCOLLIDER, DIRECTION
      DOUBLE PRECISION TEMPY(*), REACTIONFLUX(*), RPAR(*), FLUX, 
     $     M,FCENT, D, PR, LOGFCENT, LOGPR, N, C, INSIDE, LOGF, F,ALPHA,
     $     TSTAR, T2STAR, T3STAR, LOWRATE, RATE, sumy, sumyprime,
     $     y(nstate), FRATE, RRATE, INERTEFFICIENCY

      

      DO I=1,REACTIONSIZE+THIRDBODYREACTIONSIZE+TROEREACTIONSIZE
     $     +LINDEREACTIONSIZE
         REACTIONFLUX(I) = 0
      END DO

c     tempy is the total number of moles, get concentration y from it
      do i=1,nstate-1
         y(i) = tempy(I)/tempy(nstate)
      end do
      y(nstate) = tempy(nstate)


C     CALCULATE THE FLUX DUE TO REACTIONS
      DO I=0,REACTIONSIZE-1
  
         FRATE = RPAR(I+1)
         IF (REACTIONARRAY(9*I+9) .EQ. 1) THEN
            RRATE = RPAR(I+1)/REACTIONRATEARRAY(5*I+5)
         ELSE
            RRATE = 0
         END IF

         RNUM = REACTIONARRAY(9*I+1)
         PNUM = REACTIONARRAY(9*I+2)
  
         DO J=1,RNUM
            FRATE = FRATE*Y(REACTIONARRAY(9*I+2+J))
         END DO

         DO J=1,PNUM
            RRATE = RRATE*Y(REACTIONARRAY(9*I+5+J))
         END DO
         
         REACTIONFLUX(I+1) = FRATE-RRATE

      END DO
      

C     CALCULATE THE FLUX DUE TO THIRDBODYREACTIONS
      DO I=0,THIRDBODYREACTIONSIZE-1
         FRATE = RPAR(REACTIONSIZE+I+1)


C     *****CALCULATE THE THIRDBODY EFFICIENCY


         inertefficiency = pressure*1e-6/8.314/temperature

         NUMCOLLIDER = THIRDBODYREACTIONARRAY(I*20+10)
         DO J=1,NUMCOLLIDER
            INERTEFFICIENCY = INERTEFFICIENCY + Y(THIRDBODYREACTIONARRAY
     $           (I*20+10+J))
     $           *(THIRDBODYREACTIONRATEARRAY(16*I+6+J)-1)
         END DO
C     *************************************

         FRATE = FRATE * INERTEFFICIENCY

         IF (THIRDBODYREACTIONARRAY(20*I+9) .EQ. 1) THEN
            RRATE = fRATE / THIRDBODYREACTIONRATEARRAY(16*I+5)
         ELSE
            RRATE = 0
         END IF

         RNUM = THIRDBODYREACTIONARRAY(I*20+1)
         PNUM = THIRDBODYREACTIONARRAY(I*20+2)
         DO J=1,RNUM
            FRATE = FRATE*Y(THIRDBODYREACTIONARRAY(20*I+2+J))
         END DO

         DO J=1,PNUM
            RRATE = RRATE*Y(THIRDBODYREACTIONARRAY(20*I+5+J))
         END DO

         REACTIONFLUX(REACTIONSIZE+I+1) = FRATE-RRATE

      END DO


C     ******CALCULATE THE FLUX DUE TO TROE REACTIONS
      DO I=0,TROEREACTIONSIZE -1

C     FIRST CALCULATE THE RATE OF TROEREACTION

         m = pressure*1e-6/8.314/temperature

         NUMCOLLIDER = TROEREACTIONARRAY(I*21 + 10)
         DO J=1,NUMCOLLIDER
            M = M + Y(TROEREACTIONARRAY(I*21+10+J))
     $           *(TROEREACTIONRATEARRAY(21*I+6+J)-1)
         END DO

         ALPHA = TROEREACTIONRATEARRAY(21*I + 17)
         TSTAR = TROEREACTIONRATEARRAY(21*I + 18)
         T2STAR = TROEREACTIONRATEARRAY(21*I + 19)
         T3STAR = TROEREACTIONRATEARRAY(21*I + 20)
         LOWRATE = TROEREACTIONRATEARRAY(21*I + 21)
         RATE = RPAR(REACTIONSIZE+THIRDBODYREACTIONSIZE
     $        +I+1)
         DIRECTION = TROEREACTIONARRAY(21*I+9)
c         write(*,*) temperature
         IF (TROEREACTIONARRAY(21*I + 21) .EQ. 0) THEN
            FCENT = (1-ALPHA)*EXP(-TEMPERATURE/T3STAR) + ALPHA
     $           *EXP(-TEMPERATURE/TSTAR) + EXP(-T2STAR/TEMPERATURE)
         ELSE if  (TROEREACTIONARRAY(21*I + 21) .EQ. 1) THEN
            FCENT = (1-ALPHA)*EXP(-TEMPERATURE/T3STAR) + ALPHA
     $           *EXP(-TEMPERATURE/TSTAR)
         else
            stop
         END IF

         D = 0.14
         PR = LOWRATE * M/RATE
         IF (FCENT .GE. 1E-30) THEN
            LOGFCENT = LOG10(FCENT)
         ELSE
            LOGFCENT = -30
         END IF


         IF (PR .GE. 1E-30) THEN
            LOGPR = LOG10(PR)
         ELSE
            LOGPR = -30
         END IF
      
         N = 0.75 - 1.27*LOGFCENT
         C = -0.4 - 0.67*LOGFCENT

         INSIDE = (LOGPR + C)/(N - D*(LOGPR + C))
         LOGF = LOGFCENT/(1 + INSIDE**2)
         F = 10**LOGF

C      9-Jul-2009: MRH
C      For unimolecular/recombination rxns, effective rate constant is
C           keff = kinf * (Pr/(1+Pr)) * F
C      For chemically-activated bimolecular rxns, effective rate constant is
C           keff = kzero * (1/(1+Pr)) * F
         IF (TROEREACTIONARRAY(I*21+1) .GT. 1 .AND.
     $            TROEREACTIONARRAY(I*21+2) .GT. 1) THEN
            FRATE = LOWRATE * (1/(1+PR)) * F
         ELSE
            FRATE = RATE * (PR/(1+PR)) * F
         END IF

         IF (TROEREACTIONARRAY(21*I+9) .EQ. 1) THEN
            RRATE = fRATE/TROEREACTIONRATEARRAY(21*I+5)
         ELSE
            RRATE = 0
         END IF

C     **********************************

         RNUM = TROEREACTIONARRAY(I*21+1)
         PNUM = TROEREACTIONARRAY(I*21+2)
         DO J=1,RNUM
            FRATE = FRATE*Y(TROEREACTIONARRAY(21*I+2+J))
         END DO

         DO J=1,PNUM
            RRATE = RRATE*Y(TROEREACTIONARRAY(21*I+5+J))
         END DO

         REACTIONFLUX(REACTIONSIZE+THIRDBODYREACTIONSIZE+I+1) = FRATE
     $        - RRATE

      END DO

C     ******CALCULATE THE FLUX DUE TO LINDEMANN REACTIONS
      DO I=0,LINDEREACTIONSIZE -1

C     FIRST CALCULATE THE RATE OF LINDEREACTION

         m = pressure*1e-6/8.314/temperature

         NUMCOLLIDER = LINDEREACTIONARRAY(I*20 + 10)
         DO J=1,NUMCOLLIDER
            M = M + Y(LINDEREACTIONARRAY(I*20+10+J))
     $           *(LINDEREACTIONRATEARRAY(17*I+6+J)-1)
         END DO

         LOWRATE = LINDEREACTIONRATEARRAY(17*I + 17)
         RATE = RPAR(REACTIONSIZE+THIRDBODYREACTIONSIZE
     $        +TROEREACTIONSIZE+I+1)
         DIRECTION = LINDEREACTIONARRAY(20*I+9)
c         write(*,*) temperature

         PR = LOWRATE * M/RATE

C      9-Jul-2009: MRH
C      For unimolecular/recombination rxns, effective rate constant is
C           keff = kinf * (Pr/(1+Pr)) * F
C      For chemically-activated bimolecular rxns, effective rate constant is
C           keff = kzero * (1/(1+Pr)) * F
         IF (LINDEREACTIONARRAY(I*20+1) .GT. 1 .AND.
     $            LINDEREACTIONARRAY(I*20+2) .GT. 1) THEN
            FRATE = LOWRATE * (1/(1+PR))
         ELSE
            FRATE = RATE * (PR/(1+PR))
         END IF

         IF (LINDEREACTIONARRAY(20*I+9) .EQ. 1) THEN
            RRATE = fRATE/LINDEREACTIONRATEARRAY(17*I+5)
         ELSE
            RRATE = 0
         END IF

C     **********************************

         RNUM = LINDEREACTIONARRAY(I*20+1)
         PNUM = LINDEREACTIONARRAY(I*20+2)
         DO J=1,RNUM
            FRATE = FRATE*Y(LINDEREACTIONARRAY(20*I+2+J))
         END DO

         DO J=1,PNUM
            RRATE = RRATE*Y(LINDEREACTIONARRAY(20*I+5+J))
         END DO

         REACTIONFLUX(REACTIONSIZE+THIRDBODYREACTIONSIZE+
     $      TROEREACTIONSIZE+I+1) = FRATE - RRATE

      END DO

      SUMYPRIME = 0
      SUMY = 0

      END SUBROUTINE REACTION_FLUX
