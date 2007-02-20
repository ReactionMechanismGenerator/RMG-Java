      SUBROUTINE RES(T, Y, YPRIME, CJ, DEL, RPAR, IPAR, SENPAR)

      IMPLICIT NONE
      INCLUDE 'reaction.fh'
      
      INTEGER IPAR(*), IRES, I, ijac, j, k, senpar
      DOUBLE PRECISION T, Y(NEQ), YPRIME(NEQ), CJ, DEL(NEQ),
     $     RPAR(REACTIONSIZE+THIRDBODYREACTIONSIZE+
     $     TROEREACTIONSIZE+nstate-1)


      CALL GETFLUX(Y, DEL, RPAR)
      
      DO I=1,Nstate
         DEL(I) = DEL(I) - YPRIME(I)
      END DO

      END


