      PROGRAM CALL_DASPK

      IMPLICIT NONE

C IN THIS CODE WE ASSUME THAT THE 
C MAXIMUM # SPECIES = 1500 
C MAXIMUM # REACTIONS = 100,000
C MAXIMUM # THIRDBODYREACTIONS = 100
C MAXIMUM # TROEREACTIONS = 100

      INCLUDE 'reaction.fh'


      INTEGER INFO(30), LRW, LIW, I, J, IDID, NSTEPS, IMPSPECIES,numiter
      DOUBLE PRECISION Y(NEQMAX), YPRIME(NEQMAX), T, TOUT, RTOL, ATOL,
     $     THERMO(SPMAX), TARGETCONC(50)

      OPEN (UNIT=12, FILE = 'SolverInput.dat', STATUS = 'OLD')

 101  Format(E24.15)

C READ THE NUMBER OF SPECIES
      READ(12,*) NSTATE, NEQ, IMPSPECIES, NUMITER
      NPARAM = INT(NEQ/NSTATE) - 1
      NSTATE = NSTATE+1
      NEQ = NSTATE*(NPARAM + 1)


      READ(12,*) (TARGETCONC(I), I=1,NUMITER) 


c     READ THE TIME AND TOUT
      READ(12,*) T, TOUT


C  READ THE CONCENTRATIONS
      if (t .eq. 0.0) then
c it is a new job and just have to read the state variables
c read the concentrations
         READ(12,*) (Y(I), i=1,nstate-1)
C     READ THE RATE OF CHANGE OF CONCENTRATIONS
         READ(12,*) (Yprime(I), i=1,nstate-1)
         do i=nstate,neq
            y(i) = 0
            yprime(i) = 0
         end do
c     READ THE INFO ARRAY
         READ(12,*) (INFO(I), I=1,30)
C     READ RTOL AND ATOL
         READ(12,*) RTOL, ATOL
         
C     READ THE THERMO OF SPECIES
         READ(12,*)(THERMO(I), I=1,NSTATE-1)
         
C     READ TEMPERATURE AND PRESSURE
         READ(12,*) TEMPERATURE, PRESSURE

         
         y(nstate) =  8.314*temperature/pressure/1e-6
         do i=1,nstate-1
            y(i) = y(i)*y(nstate)
         end do
            

         
C     READ INFORMATION ABOUT REACTIONS
         READ(12,*) REACTIONSIZE
         READ(12,*) (REACTIONARRAY(I),I=1,9*REACTIONSIZE)
         READ(12,*) (REACTIONRATEARRAY(I),I=1,5*REACTIONSIZE)

C     READ INFORMATION ABOUT THIRDBODYREACTIONS
         READ(12,*) THIRDBODYREACTIONSIZE
         READ(12,*)(THIRDBODYREACTIONARRAY(I),I=1,20*THIRDBODYREACTIONSI
     $        ZE)
         READ(12,*)(THIRDBODYREACTIONRATEARRAY(I),I=1,16*THIRDBODYREACTI
     $        ONSIZE)
         
C     READ INFORMATION ABOUT TROEREACTIONS
         READ(12,*) TROEREACTIONSIZE
         READ(12,*)(TROEREACTIONARRAY(I),I=1,21*TROEREACTIONSIZE)
         READ(12,*)(TROEREACTIONRATEARRAY(I),I=1,21*TROEREACTIONSIZE)
         CLOSE(12)
      ELSE
         OPEN(UNIT=13, FILE = 'variables.dat', form='unformatted')
c the first line is neq read that
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

c     c write the reaction info
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
         do i=1,nstate-1
            read(13) thermo(i)
         end do
         read(13) temperature, pressure
         close(13)
      END IF


      if (NUMITER .eq. 1) then
         CALL SOLVEODE(Y, YPRIME, T, TOUT, INFO, RTOL, ATOL, IDID,
     $        THERMO, IMPSPECIES, TARGETCONC(1))
      else
         CALL SOLVESEN(Y,YPRIME,T, TOUT, INFO, RTOL, ATOL, IDID,
     $        THERMO, IMPSPECIES, NUMITER, TARGETCONC)
      end if
         
      END PROGRAM CALL_DASPK



c*******************************************************************

      SUBROUTINE PSOL()
      END SUBROUTINE PSOL

c*******************************************************************
      SUBROUTINE G_RES()
      END SUBROUTINE G_RES

c*******************************************************************


      SUBROUTINE SOLVEODE(Y, YPRIME, Time, TOUT, INFO, RTOL, ATOL,
     $     IDID, THERMO, IMPSPECIES, TARGETCONC)

      IMPLICIT NONE
      
C     INITIALIZE VARIABLES IN COMMON BLOC
      INCLUDE 'reaction.fh'


      INTEGER  INFO(30),LIW,LRW, IWORK(41+2*NSTATE), SENPAR(1), IPAR(1),
     $     IDID, ires, I, J, iter, IMPSPECIES

      DOUBLE PRECISION Y(NEQ), YPRIME(NEQ), Time, TOUT, RTOL, ATOL
     $     , RWORK(51+9*NEQ+NSTATE**2), THERMO(SPMAX), TFINAL,
     $     TARGETCONC

      DOUBLE PRECISION RPAR(REACTIONSIZE+THIRDBODYREACTIONSIZE+
     $     TROEREACTIONSIZE+NSTATE-1), del(neq), CJ,
     $     TOTALREACTIONFLUX(REACTIONSIZE+THIRDBODYREACTIONSIZE+
     $     TROEREACTIONSIZE),
     $     CURRENTREACTIONFLUX(REACTIONSIZE+THIRDBODYREACTIONSIZE+
     $     TROEREACTIONSIZE),
     $     PREVREACTIONFLUX(REACTIONSIZE+THIRDBODYREACTIONSIZE+
     $     TROEREACTIONSIZE), PREVTIME


      EXTERNAL RES, JAC, PSOL, G_RES


      LIW = 41 + NSTATE + NSTATE
      LRW = 51 + 9*NEQ + NSTATE**2

 101  Format(E24.15)
c     READ RWORK AND IWORK
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



C     ******* INITIALIZE THE REAL ARRAY
      DO I=0,REACTIONSIZE-1
         RPAR(I+1) = REACTIONRATEARRAY(5*I+1)
      END DO

      DO I=0,THIRDBODYREACTIONSIZE-1
         RPAR(REACTIONSIZE+I+1) = THIRDBODYREACTIONRATEARRAY(16*I+1)
      END DO

      DO I=0,TROEREACTIONSIZE-1
         RPAR(REACTIONSIZE+THIRDBODYREACTIONSIZE+I+1) =
     $        TROEREACTIONRATEARRAY(21*I+1)
      END DO

      DO I=0,NSTATE-2
         RPAR(REACTIONSIZE+THIRDBODYREACTIONSIZE+TROEREACTIONSIZE+i+1) =
     $        THERMO(I+1)
      END DO

C     **********************************

c     initialize the yprimes for state variables and sensitivities
C     of the second kind
      
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

      iter =0
      PREVTIME = TIME
      DO I=1,REACTIONSIZE+THIRDBODYREACTIONSIZE+TROEREACTIONSIZE
         PREVREACTIONFLUX(I) = 0
         TOTALREACTIONFLUX(I) = 0
         CURRENTREACTIONFLUX(I) = 0
      END DO
      CALL REACTION_FLUX(Y, PREVREACTIONFLUX, RPAR)

      IF (IMPSPECIES .EQ. -1) THEN

 1       IF (Time .LE. TOUT)
     $        THEN
            CALL DDASPK(RES, NEQ, Time, Y, YPRIME, TOUT, INFO, RTOL, 
     $           Atol,IDID,RWORK, LRW, IWORK, LIW, RPAR, IPAR, JAC,PSOL,
     $           SENPAR, G_RES)
            
            IF (IDID .EQ. -1) THEN
               WRITE(*,*) "Warning: 500 steps were already taken, taking
     $              another 500 steps"
               INFO(1) = 1
               GO TO 1
            else if (idid .eq. 1) then
               iter = iter+1
               
               CALL REACTION_FLUX(Y, CURRENTREACTIONFLUX, RPAR)
              DO I=1,REACTIONSIZE+THIRDBODYREACTIONSIZE+TROEREACTIONSIZE
                  TOTALREACTIONFLUX(I) = TOTALREACTIONFLUX(I) +
     $                 (PREVREACTIONFLUX(I) + CURRENTREACTIONFLUX(I))*
     $                 (TIME - PREVTIME)/2
                  PREVREACTIONFLUX(I) = CURRENTREACTIONFLUX(I)
               END DO
               PREVTIME = TIME
               
               go to 1
            END IF
         END IF
      ELSE
 2       IF (Y(IMPSPECIES) .GE. TARGETCONC*Y(NSTATE))
     $        THEN
            CALL DDASPK(RES, NEQ, Time, Y, YPRIME, TOUT, INFO, RTOL, 
     $           Atol,IDID,RWORK, LRW, IWORK, LIW, RPAR, IPAR, JAC,PSOL,
     $           SENPAR, G_RES)
            
            IF (IDID .EQ. -1) THEN
               WRITE(*,*) "Warning: 500 steps were already taken, taking
     $              another 500 steps"
               INFO(1) = 1
               GO TO 1
            else if (idid .eq. 1) then
               iter = iter+1
               
               CALL REACTION_FLUX(Y, CURRENTREACTIONFLUX, RPAR)
              DO I=1,REACTIONSIZE+THIRDBODYREACTIONSIZE+TROEREACTIONSIZE
                  TOTALREACTIONFLUX(I) = TOTALREACTIONFLUX(I) +
     $                 (PREVREACTIONFLUX(I) + CURRENTREACTIONFLUX(I))*
     $                 (TIME - PREVTIME)/2
                  PREVREACTIONFLUX(I) = CURRENTREACTIONFLUX(I)
               END DO
               PREVTIME = TIME
               tout = 200*time
               go to 2
            END IF
         END IF
      END IF

      OPEN(UNIT=15, FILE='SolverOutput.dat')

 100  Format(E24.15)
      write(15,*) iter
      WRITE(15,*) (NSTATE-1)*(nparam+1)      
      WRITE(15,*) TIME
      do i=0, nparam
         DO j=1,NSTATE
            WRITE(15,*) Y(I*nstate+j)/Y(NSTATE)
         END DO
      end do

      do i=0, nparam
         DO j=1,NSTATE
            WRITE(15,*) Yprime(I*nstate+j)/Y(NSTATE)
         END DO
      end do

      DO I=1,REACTIONSIZE+THIRDBODYREACTIONSIZE+TROEREACTIONSIZE
         WRITE(15,*) TOTALREACTIONFLUX(I)
      END DO
      

      DO I=1,30
         WRITE(15,*) INFO(I)
      END DO

      write(15,100) y(NSTATE)
      CLOSE(16)


c     write formatted restart data to the file
      open (unit=16, file='variables.dat', form='unformatted')
      write(16) idid
      write(16) time
      write(16) nstate
      write(16) nparam
      write(16) neq
      do i=1, neq
         write(16) y(i)
         write(16) yprime(i)
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
      do i=1,nstate-1
         write(16) thermo(i)
      end do
      write(16) temperature, pressure
      close(16)

      OPEN(UNIT=17, FILE='RWORK.DAT', FORM='unformatted')
      write(17) lrw
      do i=1, lrw
         WRITE(17) RWORK(I)
      end do
      CLOSE(17)

      OPEN(UNIT=18, FILE='IWORK.DAT', FORM='unformatted')
      write(18) liw
      do i=1,liw
         WRITE(18) IWORK(I)
      end do
      CLOSE(18)


      
      if (idid .eq. 1 .or. idid .eq. 2 .or. idid .eq. 3) then
         WRITE(*,*) "******ODESOLVER SUCCESSFUL: IDID=",idid
      else
         WRITE(*,*) "******ODESOLVER FAILED : IDID=", idid
      end if

      END SUBROUTINE SOLVEODE
            
c     *****************************************************************


      SUBROUTINE SOLVESEN(Y, YPRIME, Time, TOUT ,INFO, RTOL, ATOL,
     $     IDID, THERMO, IMPSPECIES, NUMITER, TARGETCONC)

      IMPLICIT NONE
      
C     INITIALIZE VARIABLES IN COMMON BLOC
      INCLUDE 'reaction.fh'


      INTEGER  INFO(30),LIW,LRW, IWORK(41+2*NSTATE), SENPAR(1), IPAR(1),
     $     IDID, ires, I, J, iter, NSTEPS, K, IMPSPECIES, NUMITER

      DOUBLE PRECISION Y(NEQ), YPRIME(NEQ), Time, TOUT, RTOL, ATOL
     $     , RWORK(51+9*NEQ+NSTATE**2), THERMO(SPMAX), TSTEP,
     $     TARGETCONC(numiter)

      DOUBLE PRECISION RPAR(REACTIONSIZE+THIRDBODYREACTIONSIZE+
     $     TROEREACTIONSIZE+NSTATE-1), del(neq), CJ,
     $     TOTALREACTIONFLUX(REACTIONSIZE+THIRDBODYREACTIONSIZE+
     $     TROEREACTIONSIZE),
     $     CURRENTREACTIONFLUX(REACTIONSIZE+THIRDBODYREACTIONSIZE+
     $     TROEREACTIONSIZE),
     $     PREVREACTIONFLUX(REACTIONSIZE+THIRDBODYREACTIONSIZE+
     $     TROEREACTIONSIZE), PREVTIME, TSTEPS(NUMITER)



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
      

C     ******* INITIALIZE THE REAL ARRAY
      DO I=0,REACTIONSIZE-1
         RPAR(I+1) = REACTIONRATEARRAY(5*I+1)
      END DO

      DO I=0,THIRDBODYREACTIONSIZE-1
         RPAR(REACTIONSIZE+I+1) = THIRDBODYREACTIONRATEARRAY(16*I+1)
      END DO

      DO I=0,TROEREACTIONSIZE-1
         RPAR(REACTIONSIZE+THIRDBODYREACTIONSIZE+I+1) =
     $        TROEREACTIONRATEARRAY(21*I+1)
      END DO

      DO I=0,NSTATE-2
         RPAR(REACTIONSIZE+THIRDBODYREACTIONSIZE+TROEREACTIONSIZE+i+1) =
     $        THERMO(I+1)
      END DO

C     **********************************

c     initialize the yprimes for state variables and sensitivities
C     of the second kind
      
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
         DO I=1,REACTIONSIZE+THIRDBODYREACTIONSIZE+TROEREACTIONSIZE
            PREVREACTIONFLUX(I) = 0
            TOTALREACTIONFLUX(I) = 0
            CURRENTREACTIONFLUX(I) = 0
         END DO
         call reaction_flux(y, prevreactionflux, rpar)
         
 1       IF (Time .LE. TSTEPS(K) .AND. Y(IMPSPECIES) .GE. TARGETCONC(K)
     $        *Y(NSTATE)) THEN
            CALL DDASPK(RES, NEQ, Time, Y, YPRIME, Tout,INFO, RTOL, 
     $           Atol, IDID, RWORK, LRW, IWORK, LIW, RPAR, IPAR, JAC, 
     $           PSOL, SENPAR, G_RES)
            
            IF (IDID .EQ. -1) THEN
               WRITE(*,*) "Warning: 500 steps were already taken, taking
     $              another 500 steps"
               INFO(1) = 1
               GO TO 1
            else if (idid .eq. 1) then
               iter = iter+1
               call reaction_flux(y, currentreactionflux, rpar)
              DO I=1,REACTIONSIZE+THIRDBODYREACTIONSIZE+TROEREACTIONSIZE
                  TOTALREACTIONFLUX(I) = TOTALREACTIONFLUX(I) +
     $                 (PREVREACTIONFLUX(I) + CURRENTREACTIONFLUX(I))*
     $                 (TIME - PREVTIME)/2
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
               WRITE(15,*) Y(I*nstate+j)/Y(NSTATE)
            END DO
         end do
         
         do i=0, nparam
            DO j=1,NSTATE
               WRITE(15,*) Yprime(I*nstate+j)/Y(NSTATE)
            END DO
         end do

         DO I=1,REACTIONSIZE+THIRDBODYREACTIONSIZE+TROEREACTIONSIZE
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
            


