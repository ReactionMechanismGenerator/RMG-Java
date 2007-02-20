c***INITRHO******************************
      subroutine initRho
c    initializes block common common rhoVar
      
      implicit none
      include 'cdparams.fh'
      include 'cdisprop.fh'
      include 'cdwell0.fh'
      include 'cdrhovar.fh'
      
c    local variables
      real*8 sumdlg,rhoOf1,rhoOf0, pi
      integer iWell,ifreq,index
      real*8 rhoOfN,dlgam
c    end declarations

      pi = 3.14159265d0
      
c    calculate prodGamS; also initialize e's and rho's
c    note we keep a set of every variable for each channel;
c    and 0 for the bottom of well

      do 120 iWell = 1, nWells
         sumdlg = 0.d0 
         do 100 ifreq = 1,nfreqs(iWell)
            sumdlg = sumdlg + dlgam(degen(iWell,ifreq))
100      continue

c       fold in active 1-D rotations; see note in calQ
         if (nRot.gt.0) then
            prodGamS(iWell) = dexp(sumdlg + dlgam(.5d0*dfloat(nRot)))
     2         * pi**(-0.5d0*dfloat(nRot))
         else
            prodGamS(iWell) = dexp(sumdlg)
         endif

c       we introduce rhoOf 0 cause we can no longer assume it
c       is unity due to possible rotations
         rhoOf0 = rhoOfN(0,iWell)	 
         rhoOf1 = rhoOfN(1,iWell)	
	 
         do 120 index = 0, nWells + nProds(iWell)	 
            nLower(iWell,index) = 0
            rhoLower(iWell,index) = rhoOf0
            rhoUpper(iWell,index) = rhoOf1
120   continue

      return
      end

c***OLDCALQ******************************
c*** tal 20000621
c*** from old version of chemdis
c*** provided to allow switching between
c*** old and new versions of density
c*** of states calculations

      function oldcalQ(T,iWell)
      
      implicit none
      include 'cdparams.fh'
      include 'cdisprop.fh'
      include 'cdrhovar.fh'
      
c    local variables
      real*8 oldcalQ
      integer iWell,ifreq
      real*8 T,sumfrac, pi, R
c    end declarations

      pi = 3.14159265d0
      R = 1.987d-3

      sumfrac = 0.d0
   
      do 100 ifreq = 1,nfreqs(iWell)
         sumfrac = sumfrac + degen(iWell,ifreq)*
     2      dlog(1.d0-dexp(-1.439d0*freq(iWell,ifreq)/T))
100   continue

c    fold in active rotations 1-D; note we factor out the 
c    Sum [1/(sigma B^1/2)] term both here and in getRho, as
c    these cancel;  
c    (we don't need B for the calculation so we don't ask)
c    also note we use R instead of kBoltzman cause our E's
c    are in kcal and our rho's are per kcal

      if (nRot.gt.0) then
         oldcalQ = dexp(-sumFrac)*(pi*R*T)**(0.5d0*dfloat(nRot))
      else
         oldcalQ = dexp(-sumFrac)
      endif
      
      return
      end
c*** tal 20000621


c***CALQ******************************
      function calQ(T,iWell, dE)
c dmm 20000511
c
c     Completely rewrote this function to calculate Q by summing up
c     xNOS results * boltzmann factor.  This should deal with the
c     obscure disagreement between xNOS and whatever it is that calRho
c     computes.  
      
      implicit real*8 (a-h, o-z)
      include 'cdparams.fh'
      include 'cdisprop.fh'
      include 'straightBS.fh'
      include 'cdnos.fh'
      
c    local variables
      real*8 calQ
      integer iWell,ifreq
      real*8 T, dE, R, sumQ, energy
c    end declarations


      R = 1.987d-3
      energy = 0.0D0
      width = dE
      Qval = 0.0D0
c     set pfflag to true:  calculating partition function
      pfflag = .true.
      tol = 1.0D-14
      zintcheck = 1.0D0
      do 10 while(zintcheck.gt.tol)
         xstates = xNOS(energy, width, iWell)     
c         Qval = Qval + xstates*exp(-(energy + 0.5*width)/(R*T))
         Qval = Qval + xstates*exp(-(energy + 0.5*width)/(R*T))

         energy = energy + width
c integration check
         zintcheck = (xstates*exp(-(energy + 0.5*width)/(R*T)))/Qval

 10   continue
      pfflag = .false.
      calQ = Qval
      
      return
      end
           
c***RHOOFN******************************
      function rhoOfN(N,iWell)

      implicit none
      include 'cdparams.fh'
      include 'cdisprop.fh'
      include 'cdrhovar.fh'
      
c    local variables
      real*8 rhoOfN
      integer N,iWell
      real*8 Eavail
      real*8 recur1,termVR
c    end declarations

c    NOTE:  we assume FORTRAN 77 no-trip do's
c    we incorporate a series of recursive subroutines recur[1]..recur[mxfreqs-1]
c    because this is fortran, we must define each one of these separately

c    NOTE:  prodGamS is calculated by initRho

c    set rho to zero if N is negative
      if (N.lt.0) then
         rhoOfN = 0.d0

c    (1-freq case, just call term)
      else if (nfreqs(iWell).eq.1) then
            rhoOfN = termVR(N,degen(iWell,1),nRot) / prodGamS(iWell)
	 
c    (n-freq case, call recur routines) 
c    note we add half quanta to E to center sum on proper energy level
      else
         Eavail = 2.859d-3*(dfloat(N)+.5d0)*freq(iWell,1) 
         rhoOfN = recur1(Eavail,iWell) / prodGamS(iWell)
      endif
      
      return
      end
           
c***CALRHO******************************
      function calRho(E,dE,iWell,index)
c    this interpolates values of rho from calculated values stored in
c    rhoVar; also calls rhoOfN to enter new values into rhoVar when
c    necessary
      
      implicit none
      include 'cdparams.fh'
      include 'cdisprop.fh'
      include 'cdrhovar.fh'
      
c    local variables
      integer nLessor,iWell,index
      real*8 calRho,rhoOfN
      real*8 E,dE,value
c    end declarations

c    note with array index we can maintain several sets of rho's simultaneously
c    in memory; this should save computation time

      if (E.lt.0.d0) then
         calRho = 0.d0
	 return
      endif
     
      nLessor = int(E/(2.859d-3*freq(iWell,1)))

c    note default (nLessor=nLower) is to do nothing 
c    (leave var's unchanged)         
      if (nLessor.eq.nLower(iWell,index)+1) then
         nLower(iWell,index) = nLower(iWell,index)+1
         rhoLower(iWell,index) = rhoUpper(iWell,index)
         rhoUpper(iWell,index) = rhoOfN(nLower(iWell,index)+1,iWell)
	 
      else if (nLessor.eq.nLower(iWell,index)-1) then
         nLower(iWell,index) = nLower(iWell,index)-1
         rhoUpper(iWell,index) = rhoLower(iWell,index)
         rhoLower(iWell,index) = rhoOfN(nLower(iWell,index),iWell)

      else if (nLessor.ne.nLower(iWell,index)) then
         nLower(iWell,index) = nLessor
         rhoUpper(iWell,index) = rhoOfN(nLower(iWell,index)+1,iWell)
         rhoLower(iWell,index) = rhoOfN(nLower(iWell,index),iWell)
      endif
      
      value = rhoLower(iWell,index) + 
     2   (rhoUpper(iWell,index) - rhoLower(iWell,index)) * 
     3   (E / (2.859d-3*freq(iWell,1)) - dfloat(nLower(iWell,index)) )

c    renormalize to proper dE (this is the sum to integral conversion)
      calRho = value * dE*
     2   (2.859d-3*freq(iWell,1))**(dfloat(nRot)/2.d0 - 1.d0) 

      return
      end
      
c****TERMV*******************************************************
      function termV(n,rs)
      implicit none
      real*8 rs
      integer n
      real*8 termV, dlgam

c    this computes (n+rs-1)!/(n)!
c    the extra factor of gam(rs) is incorporated in prodGamS 

c    this routine returns gam(rs) = (rs-1)! if n = 0, since dlgam(1) = 0
c    however, we must force to zero for n<0
    
      if (n.ge.0) then 
         termV = dexp(dlgam(dfloat(n) + rs) - dlgam(dfloat(n) + 1.d0))	 
      else
         termV = 0.d0
      endif
      
      return
      end
      
c****TERMVR*******************************************************
      function termVR(n,rs,nR)
      implicit none
      real*8 rs,sum,power
      integer n,nR
      integer i
      real*8 termVR,termV

c    this incorporates internal rotations into the last sum for rho of E

c    Normalizations are partially accounted for in calQ; and the factor
c    gam(nR/2) is incorporated in prodGamS; other factors not included
c    are expected to cancel in both K(E) and p(E)/Q; see note in calQ
    
c dmm test 
c      write(6,*)  rs
      if (nR.eq.0) then 
         termVR = termV(n,rs)

      else
         sum = 0.d0
         power = .5d0*dfloat(nR) - 1.d0

c       sum just to n-1 because Rot rho(0) should equal 1, not E^power
         do 10 i = 1,n-1
            sum = sum + termV(i,rs) * (dfloat(n-i))**power
10       continue

c       do last n term separately, letting rot part be 1
         termVR = sum + termV(n,rs)
      endif
            
      return
      end
      
c****************************************************************
c now we have 2 versions of dlgam - one adapted from joe which
c uses stirling's formula directly and one from numeric recipes
c we also use the relation: n*gamma(n) = gamma(n+1) to evaluate
c for arguments less than 1
c presently version 2 is used (it returns closer to 0 for z = 1 and 2)

C*** DLGAM (VER 2)***********************************************

c****************************************************************
c dmm 20000225
c
c altered to use lookup array
c
c 
c****************************************************************

      function dlgam(z)
c    adapted from numerical recipies; press et al. <ayc 3/93>
      implicit none

      include 'cdgamma.fh'

      real*8 dlgam,cof(6),stp,x,tmp,ser,gamln,z,xx
      integer j
      
c dmm 20000224
c      real*8 zz
c      save zz
c dmm 20000224

      data cof,stp/76.18009173d0,-86.50532033d0,24.01409822d0,
     2   -1.231739516d0,0.120858003d-2,-0.536382d-5,2.50662827465d0/

c dmm 20000224
c      zz = dmax1(z, zz)
c      if (zz.eq.z) then
c         write(*, 10) zz
c      endif
c 10   format(1X, 'Max dlgam arg: ', E14.7)
c dmm 20000224

c dmm 20000225
      if( (z.le.gzmax).and.(z.ge.gzinc)) then
         dlgam = dlgtab( nint(z/gzinc))
         return
      else

         if (z.ge.1.d0) then
            xx = z
         else
            xx = z + 1.d0
         endif
         
         x = xx - 1.d0
         tmp = x + 5.5d0
         tmp =(x + 0.50d0)*dlog(tmp) - tmp
         ser = 1.d0
         do 11 j = 1,6
            x = x + 1.d0
            ser = ser + cof(j)/x
 11      continue
         
         gamln = tmp + dlog(stp*ser)
         
         if (z.gt.1.d0) then
            dlgam = gamln
         else
            dlgam = gamln - dlog(z)
         endif
         
         return
      endif
      end
c***********************************************************************
c dmm 20000224
c
c     Gamma function table generation
c     
c     This function generates a lookup table for the gamma function
c     above, dlgam
c***********************************************************************

      subroutine dlgamgen

      implicit real*8 (a-h, o-z)      
      include 'cdgamma.fh'
      
      lgam = 105
c     If the GAMMA keyword did not appear (user did not specify gamma
c     table parameters), then use default
      if(.not.gammaspec) then
         gzinc = 0.01
         gzmax = 500.0
      endif
      
      open(lgam, file = 'cdlgtable.dat', status = 'new')
      write(lgam, 199) gzinc, gzmax
      do 10 i = 1, (nint(gzmax/gzinc) + 1)
         z = z + gzinc
         value = dlgam(z)
         write(lgam, 200) i, z, value
c         write(*, 200) i, z, value

 10   continue
 199  format(F9.4, 1X, F9.4)
 200  format(I5, 1X, F8.3, 1X, E16.9)
      close(lgam, status = 'keep')
      
      end
      
c***********************************************************************
c dmm 20000224
c
c     loads in gamma function into vector
c***********************************************************************
      subroutine dlgamload


      implicit real*8 (a-h, o-z)
      include 'cdgamma.fh'
      
      lgam = 105

 10   continue

      open(lgam, file = 'cdlgtable.dat', status = 'old')
      rewind lgam
      read(lgam, *) xinc, xmax

c     gzinc and gzmax should have been read in already by getparams if
c     GAMMA keyword was used
      if(gammaspec) then
         if ( (xinc.ne.gzinc).or.(xmax.ne.gzmax) ) then
c            gzinc = xinc
c            gzmax = xmax
            close(lgam, status = 'delete')
            call dlgamgen
            goto 10
         endif
      endif

      if(.not.gammaspec) then
         gzinc = xinc
         gzmax = xmax
      endif

      do 20 i = 1, nint(xmax/xinc)
         read(lgam, *) j, z, value
         dlgtab(j) = value
 20   continue

      end

C*** DLGAM (VER 1)*******************************************
      function dlgamo(z)
      implicit none
      real*8 dlgamo,x,z,squarx,pi,dlgamx,term1,term2,term3,series
      
c    we are experiencing overflow problems possible due to this
c    subroutine - we try reformulation <ayc 3/93>

      pi=3.141592654d0
	 
      if (z.gt.1.d0) then
         x = z
      else
         x = z + 1.d0
      endif
	 
      term1 = (x-0.5d0)*dlog(x)
      term2 = 0.5d0*dlog(2.d0*pi)
      squarx = x*x
      series = - 571.d0/(2488320.d0*squarx*squarx)
     2   - 139.d0/(51840.d0*X*squarx)
     3   + 1.d0/(288.d0*squarx) + 1.d0/(12.d0*X) + 1.d0 
      term3 = dlog(series)
      dlgamx = term1 + term2 + term3 - x
	 
      if (z.gt.1.d0) then
         dlgamo = dlgamx
      else
         dlgamo = dlgamx - dlog(z)
      endif
	 
      return
      end


c***********************************************************************

c $Id$
c $Author$
c $Date$
c $Source$
c $Revision$
c $Log$
c Revision 1.1  2007-02-20 23:10:23  sandeeps
c Initial revision
c
c Revision 1.10  2000/06/23 17:43:34  dmmatheu
c after Tom lada merge
c
c Revision 1.9  2000/06/23 14:37:25  dmmatheu
c After fix to calQ bug ... seems to be working better and more
c consistently.  Passes HP limit test for HONO2 example.
c
c Revision 1.8  2000/06/22 17:13:52  dmmatheu
c Before modification to fix "greater than HP limit bug" for dissoc.
c dE used in Qcal needs to be consistent with that used for fofEnum, and
c since this grainsize is clumsily large (1 kcal) at default, also use
c "half-width" technique ...
c
c Revision 1.7  2000/05/18 02:46:16  dmmatheu
c new calQ version ...
c
c Revision 1.6  2000/05/11 15:55:22  dmmatheu
c before modifying calQ to agree with our xNOS.  Still can't figure out
c source of large differences between our xNOS and calRho ... will have
c to do this eventually if we want to compare DOS methods ...
c
c Revision 1.5  2000/03/21 22:43:30  dmmatheu
c before commenting out genRhos to deal with testing ...
c
c Revision 1.4  2000/03/09 14:33:47  dmmatheu
c before adding BS algorithm and calcRhos routine
c
c Revision 1.3  2000/02/25 16:31:20  dmmatheu
c before lookup table alteration
c
c Revision 1.2  2000/02/24 22:14:47  dmmatheu
c before dlgamgen
c
