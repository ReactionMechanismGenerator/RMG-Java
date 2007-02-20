      SUBROUTINE RES(T, Y, YPRIME, CJ, DEL, IRES, RPAR, IPAR, SENPAR)

      IMPLICIT NONE
      INCLUDE 'reaction.fh'
      
      INTEGER IPAR(*), IRES, I, ijac, j, k, senpar
      DOUBLE PRECISION T, Y(NEQ), YPRIME(NEQ), CJ, DEL(NEQ),
     $     RPAR(REACTIONSIZE+THIRDBODYREACTIONSIZE+
     $     TROEREACTIONSIZE+nstate-1), SUMY, SUMYPRIME,
     $     Jacobian(NSTATE,NSTATE), tempcj,
     $     dgdk(NSTATE, REACTIONSIZE+THIRDBODYREACTIONSIZE+
     $     TROEREACTIONSIZE+NSTATE-1)


      do i=1,nstate
         do j=1,nstate
            jacobian(i,j) = 0
         end do
      end do

      do i=1,nstate
         do j=1,reactionsize+thirdbodyreactionsize+troereactionsize+
     $        nstate-1
            dgdk(i,j) = 0
         end do
      end do

      if (ires .eq. 0) then
         CALL GETFLUX(Y, DEL, RPAR)

         DO I=1,Nstate
            DEL(I) = DEL(I) - YPRIME(I)
         END DO

      else if (ires .eq. 1) then

c calculate the res for state variables
         CAll GETFlux(y, DEL, RPAR)

c calculate the res for sensitivity variables
         ijac = 0
         tempcj =0
         call jac(t, y, yprime, jacobian, tempcj, rpar, ipar, senpar,
     $        ijac)
         call calcdgdk(t, y, yprime, dgdk, tempcj, rpar, ipar, senpar)
         do j=1,reactionsize+thirdbodyreactionsize+troereactionsize+
     $        nstate-1
c     parameter number
            do i=1,nstate
c     the equation
               do k=1, nstate
c     the state variable
                  del(j*nstate + i) = del(j*nstate +i) + jacobian(i,k)*
     $                 y(j*nstate + k) 
               end do
               del (j*nstate + i) = del(j*nstate+i) + dgdk(i,j)
            end do
         end do



         do i=1, neq
            del(i) = del(i) - yprime(i)
         end do
      end if

      END


