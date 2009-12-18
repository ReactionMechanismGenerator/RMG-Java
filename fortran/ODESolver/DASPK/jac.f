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
         
         oloop: DO I=1,NSTATE
            PD(I,I) = PD(I,I) - CJ
            !gmagoon 12/18/09: set rows for equations for constant concentration species equal to all zeroes
            IF (ConstantConcentration(I) .EQ. 1) THEN
                iloop: DO J=1,NSTATE
                    PD(I,J) = 0
                END DO iloop
            END IF
         END DO oloop

      ELSE 
         WRITE(*,*) "IJAC NOT EQUAL TO 0", IJAC
      END IF


      END SUBROUTINE JAC
