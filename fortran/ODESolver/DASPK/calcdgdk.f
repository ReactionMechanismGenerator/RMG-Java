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
      
      
      END SUBROUTINE calcdgdk
