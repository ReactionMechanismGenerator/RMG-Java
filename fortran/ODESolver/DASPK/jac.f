      SUBROUTINE JAC(T, Y, YPRIME, PD, CJ, RPAR, IPAR, SENPAR, IJAC)

      IMPLICIT NONE
      INCLUDE 'reaction.fh'


      DOUBLE PRECISION T, Y(NEQ), YPRIME(NEQ), PD(Nstate,Nstate),
     $     CJ, RPAR(*), DEL(NEQ), AJAC(nstate*nstate)
      INTEGER IPAR(*), SENPAR(*), IJAC, NJAC,IROW(nstate*nstate), 
     $     JCOL(nstate*nstate), IDPWRK(nstate*nstate), IRES, I, J

      IF (IJAC .EQ. 0) THEN

         call RESAD(T, Y, YPRIME, CJ, DEL, RPAR, IPAR, SENPAR,
     $        AJAC, NJAC, IROW, JCOL, IDPWRK)

         
         DO I=1,NJAC
            PD(IROW(I), JCOL(I))= AJAC(I)
         END DO
         
         DO I=1,NSTATE
            PD(I,I) = PD(I,I) - CJ
         END DO

      ELSE 
         WRITE(*,*) "IJAC NOT EQUAL TO 0", IJAC
      END IF


      END SUBROUTINE JAC
