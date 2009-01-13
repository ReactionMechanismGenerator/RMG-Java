C     ************************************************
      SUBROUTINE GETFLUX(TEMPY, DEL, RPAR) 

      IMPLICIT NONE

C     INITIALIZE VARIABLES IN COMMON BLOC
      INTEGER NSTATE, REACTIONSIZE, THIRDBODYREACTIONSIZE,
     $     TROEREACTIONSIZE
      COMMON /SIZE/ NSTATE, REACTIONSIZE, THIRDBODYREACTIONSIZE,
     $     TROEREACTIONSIZE

      INTEGER REACTIONARRAY(9*100000), THIRDBODYREACTIONARRAY(20*100),
     $     TROEREACTIONARRAY(21*100), ConstantConcentration(1501)

      DOUBLE PRECISION REACTIONRATEARRAY(5*100000),
     $     THIRDBODYREACTIONRATEARRAY(16*100),
     $     TROEREACTIONRATEARRAY(21*100), temperature, pressure

      COMMON /REAC/ REACTIONRATEARRAY, ThirdbodyREACTIONRATEARRAY,
     $     TROEREACTIONRATEARRAY,temperature, pressure, 
     $     REACTIONARRAY, THIRDBODYREACTIONARRAY, TROEREACTIONARRAY,
     $     ConstantConcentration


      
      INTEGER I, J, RNUM, PNUM, NUMCOLLIDER, DIRECTION
      DOUBLE PRECISION TEMPY(*), DEL(*), RPAR(*), FLUX, INERTEFFICIENCY,
     $     M,FCENT, D, PR, LOGFCENT, LOGPR, N, C, INSIDE, LOGF, F,ALPHA,
     $     TSTAR, T2STAR, T3STAR, LOWRATE, RATE, sumy, sumyprime,
     $     y(nstate), FRATE, RRATE

      

      DO I=1,NSTATE
         DEL(I) = 0
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
         
         DO J=1,RNUM
            DEL(REACTIONARRAY(9*I+2+J)) = DEL(REACTIONARRAY(9*I+2+J))
     $           - FRATE + RRATE
         END DO

         DO J=1,PNUM
            DEL(REACTIONARRAY(9*I+5+J)) = DEL(REACTIONARRAY(9*I+5+J))
     $           + FRATE - RRATE
         END DO

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

         DO J=1,RNUM
            DEL(THIRDBODYREACTIONARRAY(20*i+2+J)) =
     $           DEL(THIRDBODYREACTIONARRAY(20*I+2+J)) - FRATE + RRATE
         END DO

         DO J=1,PNUM
            DEL(THIRDBODYREACTIONARRAY(20*I+5+J)) =
     $           DEL(THIRDBODYREACTIONARRAY(20*I+5+J)) + FRATE - RRATE
         END DO

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

         FRATE = RATE * (PR/(1+PR)) * F
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

         DO J=1,RNUM
            DEL(TROEREACTIONARRAY(21*i+2+J)) = DEL(TROEREACTIONARRAY
     $           (21*I+2+J)) - FRATE + RRATE
         END DO

         DO J=1,PNUM
            DEL(TROEREACTIONARRAY(21*I+5+J)) = DEL(TROEREACTIONARRAY
     $           (21*I+5+J)) + FRATE - RRATE
         END DO

      END DO

      SUMYPRIME = 0
      SUMY = 0

c del contains rate of change of concentrations, change it to rate of
c change of moles
      DO I=1, NSTATE-1
         DEL(I) = DEL(I) * tempY(NSTATE)
      END DO
      
c   NO CHANGE IN FIRST SPECIES CONCENTRATION!!  
c   (eg. dissolved oxygen in liquid phase)
c      DEL(1)=0

      DO I=1,NSTATE-1
         SUMYPRIME = sumyprime + DEL(I)
      END DO

c this is the rate of change of volume
      DEL(NSTATE) = sumyprime*8.314*temperature/pressure/1e-6
      
c   NO CHANGE IN VOLUME!
c      DEL(NSTATE) = 0
      DO I=1,NSTATE
         IF (constantConcentration(I) .EQ. 1) THEN
            DEL(I) = 0
         END IF
      END DO
      
      
c      del(nstate) = 0
c      write(*,*) sumyprime
      END SUBROUTINE GETFLUX
