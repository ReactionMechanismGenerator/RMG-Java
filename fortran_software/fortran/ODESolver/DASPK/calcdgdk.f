      SUBROUTINE calcdgdk(T, Y, YPRIME, PD, CJ, RPAR, IPAR, SENPAR)

      IMPLICIT NONE
      INCLUDE 'reaction.fh'


      DOUBLE PRECISION T, Y(NEQ), YPRIME(NEQ), PD(Nstate,nparam),
     $     CJ, RPAR(*), DEL(Neq), AJAC(nstate*nparam)

      INTEGER IPAR(*), SENPAR(*), IJAC, NJAC,IROW(nstate*nparam), 
     $     JCOL(nstate*nparam), IDPWRK(nstate+nparam), IRES, I, J


c      NSTATE = IPAR(1)

      call RES_rpar(T, Y, YPRIME, CJ, DEL, RPAR, IPAR,SENPAR,
     $     AJAC, NJAC, IROW, JCOL, IDPWRK)
      
      
      do i=1, nstate
         do j=1,nparam
            pd(i,j) = 0
         end do
      end do
      
      DO I=1,NJAC
         PD(IROW(I), JCOL(I))= AJAC(I)
      END DO
      
      !gmagoon 12/21/09: set rows for equations for constant concentration species equal to all zeroes, similar to what was done in jac; I'm not exactly sure what res_rpar does, so I don't really have a good justification for this, but assuming res_rpar calculates dg/dk, then I think this should be accurate
      oloop: DO I=1,NSTATE
         IF (ConstantConcentration(I) .EQ. 1) THEN
             iloop: DO J=1,NPARAM
                 PD(I,J) = 0
             END DO iloop
         END IF
      END DO oloop
      
      END SUBROUTINE calcdgdk
