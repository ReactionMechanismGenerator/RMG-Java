



C***********************************************************************
C**************************SUBROUTINE BS ******************************
C***********************************************************************
c
c dmm 20000310
c     Stripped out of QCPE routine (QCPE557?) and hindrot routine,
c     simplified and turned into 'pure' BS: single T array, frequencies
c     only.
c
c     dmm 2000318  Not debugged yet; removed DEGEN and will code up for
c     a set of classical 1-d rotors + frequencies.  This requires us to
c     use Asholt's method as described in Gilbert and Smith (see binder
c     notes, 3/16/00).  Note that this code returns the SUM of states
c     in the TBS vector; the NUMBER of states would be TBS(I) - TBS(I-1)
c     , which would be the number of states in width DE.
c     
c
c     NMAX = max no. of grains to count up to
c     NV = number of vibrational modes
c     VIB = vibrational mode vector (cm-1)
c     nrot = number of classical independent 1-d rotors
c     BR = vector of rotational constants (cm-1)
c     BDG = rotational mode degeneracies
c     DE = grainsize (cm-1)
c**********************************************************************
      SUBROUTINE SBS(NMAX,NV,VIB,nrot,BR,BDG,DE)

      IMPLICIT REAL*8(A-H,O-Z)     
      
      include 'straightBS.fh'

      integer NMAX, NV, nrot
      real*8 VIB(maxmodes)
      real*8 BR(maxmodes)
      real*8 BDG(maxmodes)


C     set EMAX

      EMAX = DBLE(NMAX-1)*DE  

c***********************************************************************
C     Initialize the TBS array.  If there are 1-d rotors, initialize
c     with the sum of states of a classical set of convolved 1-d rotors
c     (see RN II p. 113)
c***********************************************************************

      if(nrot.gt.0) then

c     Check degeneracies of rotors; if they are 0, set them to 1.
c     Calculate the sum-of-states constant

         Bprod=1.0
         do 30 i = 1, nrot
            if(BDG(i).le.0.0D0) BDG(i) = 1.D0
            Bprod = Bprod*BR(i)*BDG(i)
 30      continue

c jmg fix 
         Crot = (2.0/dfloat(nrot))*exp(-gammln( (dfloat(nrot)/2.0) ))*( 
     $        (3.14159)**((dfloat(nrot)/2.0)) )*(1.0/(Bprod**0.5))

c     Initialize the TBS vector, assuming W(0) (the number of rotational
c     states at E = 0) is 1; changed to 0.99 to prevent numerical
c     occurence of negative DOS when small BS grainsize is used ...

         TBS(1) = 0.99D0
         do 35 i = 2, NMAX
cjmg fix
            TBS(i) = Crot*( (DE*dfloat(i-1))**( dfloat(nrot)/2.0 ))
c           TBS(i) = Crot*((DE*dfloat(i-1))**(dfloat(nrot/2)))
 35      continue

      else
c     There are no rotors to convolve ...

         TBS(1) = 1.0D0
         
         DO 10 I= 2,NMAX     
           TBS(I) = 1.D0
c            TBS(I) = 0.0D0
 10      CONTINUE 
      endif

C---- >COUNT VIBRATIONAL STATES
      DO 99 IV = 1,NV
         IF (DMOD(VIB(IV),DE).EQ.0.D00) THEN
            KV = nint(VIB(IV)/DE)
            ISTART = KV + 1
            DO 20 I = ISTART,NMAX
               TBS(I) = TBS(I) + TBS(I-KV)
 20         CONTINUE
         ELSE
                                                             
c     DO 60 I = 1,NMAX
c     TBS(I) = AT(I)
c     60         CONTINUE
            
            RV = VIB(IV)/DE
            KMAX = int(EMAX/VIB(IV))
            KV = int(VIB(IV)/DE)

 
               
c     RKV = DBLE(K)*RV
c     KV = RKV
c     IF (RKV-DBLE(KV).GE.0.5D00) KV = KV + 1
            
            ISTART = KV + 1
            DO 70 I = ISTART,NMAX                                   
               TBS(I) = TBS(I) + TBS(I-KV)
 70         CONTINUE                

         ENDIF   
 99   CONTINUE  
      
      
      RETURN
      END                                                                       

c $Id$
c $Author$
c $Date$
c $Source$
c $Revision$
c $Log$
c Revision 1.1  2007-02-20 23:10:23  sandeeps
c Initial revision
c
c Revision 1.8  2003/04/23 19:17:12  dmmatheu
c test revision checkin
c
c Revision 1.6  2001/10/18 20:34:02  dmmatheu
c before examining effect of changes to W(0) = 1 assumption (changing to
c prevent numerical difficulty with a first rho as negative)
c
c Revision 1.5  2000/03/24 02:03:29  dmmatheu
c Debugged -- erroneous KMAX not needed ... use 'nint' for counting of
c vibrational states ...
c
c Revision 1.4  2000/03/21 22:15:28  dmmatheu
c Tested on SRtest.XS4, sheet3, athena/research directory, with one free
c rotor.  Works great.
c
c Revision 1.3  2000/03/14 15:39:58  dmmatheu
c tested once; moved on to initialization of T vector with 1-d rotors
c ...
c
c Revision 1.2  2000/03/09 16:46:40  dmmatheu
c First try at straightBS ... untested
c
c Revision 1.1  2000/03/09 16:19:40  dmmatheu
c Initial revision
c


















