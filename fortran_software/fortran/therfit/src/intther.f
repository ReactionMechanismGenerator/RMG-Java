c***********************************************************************
c***********************************************************************
c***********************************************************************
c     intther.f:
c
c     Collection of Modified Routines for THERFIT
c     dmmatheu 20000529
c
c     These routines allow modification of the THERFIT code to allow
c     integer-only frequency degeneracies for fitting the 3-frequency
c     model to heat capacity data.  It uses the mrqmin routine family
c     from Numerical Recipes for the marquardt minimization.
c
c     should be compiled and linked with therfit.f ...
c***********************************************************************
c***********************************************************************
c***********************************************************************


c***********************************************************************
c***********************************************************************
c     SUBROUTINE cpfunc
c
c     This subroutine returns predicted heat capacity y for a
c     temperature x, and the vector of derivatives dy/db given the
c     parameter vector b.  na is the number of these parameters
c
c     x = temperature, K
c     y = Cp, cal/mol*K, *CONST(1) which is a unit conversion factor
c         (usually 1)
c     b = coefficients:  b(1) = 1st mode degeneracy
c                        b(2) = 2nd mode degeneracy
c                        b(3) = 1st mode frequency (cm-1)
c                        b(4) = 2nd mode frequency
c                        b(5) = 3rd mode frequency
c     3rd mode degeneracy given by CONST(2) - b(1) - b(2)
c     
c     na = number of coefficients (total -- obviously must be 5!)
c***********************************************************************
c***********************************************************************

      subroutine cpfunc(x, b, y, dydb, nb)

      implicit real*8 (a-h, o-z)
      integer nb, i, j
      real*8 x, y, b(nb), dydb(nb)

      include 'const.fh'

c      real  const(3)
c      integer nterms

c      common /const/ CONST, NTERMS

c     First, call FOFX to get value of y; first y is dummy arg.

      i = 1
      call fofx(i, x, b, ydummy, y)

c     Next, compute dydb vector (see binder notes 20000529)
c     hkB = hc/kB in K/cm-1
c     Rconst = gas constant in cal/mol K

      Rconst = 1.987*CONST(1)
      sigma3 = CONST(2)-b(1)-b(2)
      
c     dy/db (1) = derivative of y wrt 1st mode degeneracy

      cv1 = 0.0
      cv2 = 0.0
      cv3 = 0.0
      theta = 0.0

      call cvvib(b(3), x, cv1, theta)
      call cvvib(b(4), x, cv2, theta)
      call cvvib(b(5), x, cv3, theta)

      dydb(1) = (cv1 - cv3)*Rconst
      dydb(2) = (cv2 - cv3)*Rconst

c     dy/db(3) = derivative of y wrt 1st mode frequency

      dydb(3) = dCvibdNu(x, b(3))*b(1)*Rconst
      dydb(4) = dCvibdNu(x, b(4))*b(2)*Rconst
      dydb(5) = dCvibdNu(x, b(5))*sigma3*Rconst

      end

c***********************************************************************
c***********************************************************************
c     FUNCTION dCvibdNu
c
c     This function returns the derivative of Cvib (as defined by stat
c     mech for quantum HO's) wrt frequency.  The gas constant R is not
c     included, and so the answer must be multiplied by R in the desired
c     units, as well as the mode degeneracy.
c     
c     x = temperature (K)
c     freq = frequency, cm-1
c
c     Reference:  dmmatheu's Maple file 'cpfitderivative.mws' in
c     research dir.
c***********************************************************************
c***********************************************************************

      function dCvibdNu(x, freq)

      implicit real*8 (a-h, o-z)
      real*8 x, freq

c     hokB = hc/kB in K/cm-1
      hokB = 1.438659929
      hokBT = hokB/x

c     dimensionless temperature, or h*nu/kT
      dimT = hokB*freq/x

      topex = exp(dimT)
      botex = exp(dimT) - 1.0

      term1 = 2.0*freq*(hokBT**2.0)*topex/(botex**2.0)
      term2 = (hokBT**3.0)*(freq**2.0)*topex/(botex**2.0)
      term3 = -2.0*(hokBT**3.0)*(freq**2.0)*(topex**2.0)/(botex**3.0)

      dCvibdNu = term1 + term2 + term3

      return
      end

c***********************************************************************
c***********************************************************************
c     SUBROUTINE IDCPFIT
c
c     This subroutine is intended as a modular replacement for THCPFI in
c     therfit.f.  It uses the same Marquardt method, this time adapted
c     from the general one of Numerical Recipes.  
c
c     This method has the advantage of giving back some fitting criteria
c     automatically, and much more importantly, allows us to 'fix'
c     parameters by turning on or off the integers of the ia vector.
c
c     The code to force the frequency degeneracies to be integers is
c     described below.  When the IDG flag is set to 1, the degeneracies
c     will be forced to be integers, and the total number of modes will
c     be also set to an integer.  Note that this doesn't change the high
c     -T limit Cp.
c     
c     NT = number of temperature data points (except high-T limit)
c     CPL = vector of Cp at various T (cal/mol*K)
c     NATOM = number of atoms in species
c     TEMPR = temperature vector (K)
c     CPINFI = high-T limit Cp
c     icp = integer flag for when an actual temperature has been
c     provided for the high-T limit.
c     thigh = given high T for the high-T limit ... otherwise this will
c     be 5000 K
c***********************************************************************
c***********************************************************************

      subroutine IDCPFIT(NT, CP, NATOM, bparm, TEMPR, CPINF, icp, thigh)

      implicit real*8 (a-h, o-z)

      include 'intdegen.fh'
      include 'const.fh'
      include 'params.fh'

      external cpfunc

      real*8 CP(maxdat), TEMPR(maxdat), CPINF, thigh
      real*8 xdat(maxdat), ydat(maxdat), aparm(maxpar), alpha(maxpar,
     $     maxpar), covar(maxpar, maxpar), sigma(maxdat), bparm(maxpar)
     
      real*8 bmin(maxpar), bmax(maxpar)

      integer ia(maxpar)
      integer maxmqit, minit
      real*8 bbparm(maxpar), rparm(maxpar), csqv(4)
      logical goodcomb
 

c***********************************************************************      
c     Copy data points to xdat and ydat arrays.  Assume std. dev. on
c     each data point is the same (!).  
c***********************************************************************

      stdev = 1.0
      do 10 i = 1, NT
         ydat(i) = CP(i)
         xdat(i) = TEMPR(i)
         sigma(i) = stdev
 10   continue

c***********************************************************************
c     Set up CPINFI as high-T data point; ndata = total number of data
c     points.  
c
c     Note that because they use so few data points (from benson format
c     -- 8 data points for 5 fitted pars), it matters inordinately what
c     Tinf is.
c***********************************************************************

      ndata = NT + 1
      Tinf = 99999.0

      if(icp.eq.1) then
         xdat(ndata) = thigh
      else
         xdat(ndata) = Tinf
      endif
      ydat(ndata) = CPINF
      sigma(ndata) = stdev

c***********************************************************************
c     Default settings -- it's not a really good idea to set these here,
c     but that is the way it's done in the original -- hard to avoid.
c***********************************************************************     

      NTERMS=3

c     Note -- don't set tolerance too tight! can cause
c     other problems ... frequencies migrate too high.  Fit will be good
c     anyway
      
c      reltol1 = 1.0D-2
      reltol1 = 5.0D-3

c      idg = 0
      maxmqit = 200
      nparm = maxpar
      ntest = nparm
 
c***********************************************************************
c     Set ia(i) matrix.  If ia(i) = 1, 'fit' aparm(i); if 0, consider
c     this value 'fixed'.  Default is to fit all except the nonexistent '
c     6th' parameter
c***********************************************************************
     
      do 20 i = 1, maxpar-1
         ia(i) = 1
 20   continue
      ia(maxpar) = 0



c***********************************************************************
c     Initial call to imqsolve;  if not fixing integer degeneracies,
c     this will be sufficient.  If we are, this is an initialization
c     call.
c***********************************************************************

c     in-house guess for aparm ... may need to be revised as per old THERFIT
      aparm(1) = 0.3*CONST(2)
      aparm(2) = 0.3*CONST(2)
      aparm(3) = 600.0
      aparm(4) = 950.0
c      aparm(4) = 1200.0
      aparm(5) = 1500.0

c     set mins and maxes ...

      bmin(1) = 1.0
      bmin(2) = 1.0
      bmin(3) = 250.0
      bmin(4) = 600.0
      bmin(5) = 800.0

      bmax(1) = 0.9*CONST(2)
      bmax(2) = 0.9*CONST(2)
      bmax(3) = 1500.0
      bmax(4) = 2500.0
      bmax(5) = 3600.0

c reset for Jing's cases ...
      bmax(5) = 5000.0

c     vinyl + O2 guess ...

c      aparm(1) = 5.052
c      aparm(2) = 5.087
c     aparm(3) = 452.6
c     aparm(4) = 1049.8
c     aparm(5) = 2915.9


      call imqsolve(xdat, ydat, sigma, ndata, aparm, ia, nparm,
     $     covar, alpha, ntest, chisq, cpfunc, alamda, maxmqit, relerr,
     $     reltol1, bmin, bmax)


c***********************************************************************
c     If fitting integer degeneracies, enter this block
c***********************************************************************

c***********************************************************************   
c     If IDG is not 0, then we are fixing the frequency degeneracies as
c     integers, and CONST(2) needs to be reset to an integer total
c     number of modes (this means the CPINF was still eval'd w/correct
c     limit for hindered rotors)
c***********************************************************************

c     idg = 1
      if(IDG.gt.0) then
         CONST(2) = float(int(CONST(2)))
      endif



      if(idg.gt.0) then
         it = 1
         ia(1) = 0
         ia(2) = 0
         csqmin = 1.0D50

         do 40 it = 1, 4
            goodcomb = .true.
            call getcombo(it, aparm, rparm)
            call chkcombo(rparm, aparm, CONST(2), goodcomb)

            if(goodcomb) then
                              
               call imqsolve(xdat, ydat, sigma, ndata, rparm, ia, nparm,
     $              covar, alpha, ntest, chisq, cpfunc, alamda, maxmqit,
     $              relerr,reltol1, bmin, bmax)

               csqv(it) = chisq

               if(chisq.lt.csqmin) then
                  csqmin = chisq
                  minit = it
                  do 50 j = 1, nparm
                     bbparm(j) = rparm(j)
 50               continue
               endif

            endif

 40      continue

c     Copy results to aparm ...
         do 60 i = 1, nparm
            aparm(i) = bbparm(i)
 60      continue


      endif

c***********************************************************************
c     End integer degeneracy fitting block
c***********************************************************************

c***********************************************************************
c     Copy results into bparm vector
c***********************************************************************

      bparm(1) = aparm(1)
      bparm(2) = aparm(2)
      bparm(3) = CONST(2) - aparm(1) - aparm(2)
      bparm(4) = aparm(3)
      bparm(5) = aparm(4)
      bparm(6) = aparm(5)

c***********************************************************************
c     Print out data points, fit to xmgr-compatible data plot
c***********************************************************************

c (placeholder ... for now in therfit.f)


      return
      end

c***********************************************************************
c***********************************************************************
c     SUBROUTINE imqsolve
c
c     This subroutine constitutes the inner marquardt iteration loop.
c
c     Uses Marquardt routine to find best fit for the unfixed
c     parameters of aparm (those for which ia(i) > 0).  Criteria for
c     term. is relative chisq tolerance between two iterations (not nec.
c     good -- should eventually change to 2-norm of update vector?).
c     
c     alamda.gt.oalamda indicates last step was "wrong direction", and
c     while relerr1 will be 0 this is because that's what mrqmin returns
c     for wrong direction search, so stay in loop (terminate with step
c     in the right chisq direction.
c     
c***********************************************************************     
c***********************************************************************

      subroutine imqsolve(xdat, ydat, sigma, ndata, aparm, ia, nparm,
     $     covar, alpha, ntest, chisq, cpfunc, alamda, maxmqit, relerr,
     $     reltol, bmin, bmax)
      
      implicit real*8 (a-h, o-z)

      include 'params.fh'

      real*8 xdat(maxdat), ydat(maxdat), sigma(maxdat), aparm(maxpar),
     $     alpha(maxpar,maxpar), covar(maxpar,maxpar), chisq, alamda,
     $     relerr, reltol, bmin(maxpar), bmax(maxpar)

      integer ndata, nparm, ntest, maxmqit
      integer ia(maxpar)

      external cpfunc

c***********************************************************************
c     Initialize loop variables, and chisq
c***********************************************************************

      chisq = 1.0D50
      i = 1
      relerr = 1.0D50
      alamda = -1.0
      oalamda = 0.0

c***********************************************************************
c     Fit parameters iteratively
c***********************************************************************

      do 30 while ( (i.lt.maxmqit).and.(((relerr.gt.reltol).or.relerr
     $     .lt.0.0).or.(alamda.gt.oalamda)))

         oalamda = alamda
         chiold = chisq
         call mrqmin(xdat, ydat, sigma, ndata, aparm, ia, nparm, covar,
     $        alpha, ntest, chisq, cpfunc, alamda, bmin, bmax)         
         i = i + 1
         relerr = (chiold - chisq)/chisq
         
 30   continue

c***********************************************************************
c     Final call to mrqmin -- find covariance matrix
c***********************************************************************

      alamda = 0.0
      call mrqmin(xdat, ydat, sigma, ndata, aparm, ia, nparm, covar,
     $     alpha, ntest, chisq, cpfunc, alamda, bmin, bmax)

c***********************************************************************
c     Errors and warnings
c***********************************************************************

      if(i.ge.maxmqit) then 
         write (*, 100), i
 100     format ('WARNING::imqsolve::intther.f:  exceeded max. allowed
     $        iterations; i = ', I5)
      endif


c     placeholder -- beyond max it. error handling; goodness of fit
c     checks, etc.

c     placeholder -- write out covariance at least for fitted pars.

      return
      end

c***********************************************************************
c***********************************************************************
c     SUBROUTINE getcombo
c
c     For the three-frequency model, there are 4 combinations of
c     rounding of the first 2 degeneracies when the total number of
c     vibrational modes has already been rounded down to the nearest int
c     .  This subroutine choses one for each iteration i when it is
c     called.  
c
c***********************************************************************
c***********************************************************************

      subroutine getcombo(it, aparm, rparm)
      
      implicit real*8 (a-h, o-z)
      include 'params.fh'
      
      integer it
      real*8 aparm(maxpar), rparm(maxpar)

c     Down-Down Combo
      if(it.eq.1) then
         rparm(1) = dfloat(int(aparm(1)))
         rparm(2) = dfloat(int(aparm(2)))
     
c     Down-Up Combo
      else if(it.eq.2) then
         rparm(1) = dfloat(int(aparm(1)))
         rparm(2) = dfloat(int(aparm(2)) + 1)

c     Up-Down Combo
      else if(it.eq.3) then
         rparm(1) = dfloat(int(aparm(1)) + 1)
         rparm(2) = dfloat(int(aparm(2)))

c     Up-Up Combo
      else if(it.eq.4) then
         rparm(1) = dfloat(int(aparm(1)) + 1)
         rparm(2) = dfloat(int(aparm(2)) + 1)

c     Invalid Combo
      else
         write(*,100), it
 100     format('ERROR::getcombo::intther.f: getcombo called with bad
     $ iteration number, it = ', I4)
      endif
      
c     Set the rest of rparm
      rparm(3) = aparm(3)
      rparm(4) = aparm(4)
      rparm(5) = aparm(5)

      return
      end

c***********************************************************************
c***********************************************************************
c     SUBROUTINE chkcombo
c
c     This routine checks whether the rounding combo is 'good', e.g.
c     third degeneracy is a round up or down ... this subroutine may not
c     actuall be necessary ...
c     
c
c***********************************************************************
c***********************************************************************

      subroutine chkcombo(rparm, aparm, vtot, goodcomb)

      implicit real*8 (a-h, o-z)
      include 'params.fh'
      
      real*8 rparm(maxpar), aparm(maxpar), vtot
      logical goodcomb

      rs3 = vtot - rparm(1) - rparm(2)
      s3 = vtot - aparm(1) - aparm(2)

      sck = dabs(rs3 - s3)
      if(sck.gt.1.0) then
         goodcomb = .false.
      else
         goodcomb = .true.
      endif

c     20031107  dmm A new problem with O-O bonds.  We sometimes get bad
C     combos which would give the third degeneracy as 0.  This gives a
C     singular  matrix in the solver, so we must throw such combos out.

      vsum = rparm(1) + rparm(2)
      if(dabs(vsum-vtot).lt.0.01) then
         goodcomb = .false.
      endif

      return
      end

c***********************************************************************
c***********************************************************************
c     SUBROUTINE printfit
c
c     This routine prints out the xmgr-compatible data table of heat
c     capacity data points, + the fit in 20 K intervals (this to check
c     against fit weirdness)
c***********************************************************************
c***********************************************************************

      subroutine printfit(NT, CP, bparm, tempr, CPINF, icp, thigh, SPEC
     $     ,irec)

      implicit real*8 (a-h, o-z)

      include 'intdegen.fh'
      include 'params.fh'
      include 'const.fh'

      integer irec, NT, icp
      character*16 SPEC
      real*8 CP(maxdat), bparm(maxpar), tempr(maxdat), CPINF, thigh
      real*8 aparm(maxpar)


      write(21, 100) SPEC(:16), 'DATA'

      irec = irec + 1
      
      do 10 i = 1, NT

         write(21, 200) tempr(i), CP(i)
         irec = irec + 1
 10   continue

      tmax = 2000.0
      twidth = 50.0
      imax = int(tmax/twidth)
      
      aparm(1) = bparm(1)
      aparm(2) = bparm(2)
      aparm(3) = bparm(4)
      aparm(4) = bparm(5)
      aparm(5) = bparm(6)

      xt = 300.0
      write(21, 100) SPEC(:16), 'PRED'
      do 20 i = 1, imax
        
         call fofx(i, xt, aparm, ydummy, yval)
         write(21,200) xt, yval
         xt = xt + twidth
         
 20   continue
 100   format(A16, 1X, A4, 15X)
 200   format(F6.1, 1X, F8.2, 15X)
 300   format(//)
      return
      end

c***********************************************************************
c***********************************************************************
c     SUBROUTINE checkbounds
c
c     This subroutine checks whether the new parameters exceed their
c     bounds ... if bounds are le 0 then these are not checked against .
c     ..
c
c
c***********************************************************************

      subroutine checkbounds(atry, amin, amax, ma, ibound)

      implicit real*8 (a-h, o-z)
      integer ma, ibound
      real*8 atry(ma), amin(ma), amax(ma)

      ibound = 0

c     examine bounds, increment ibound for each over bound

      do 20 j = 1, ma
         if( (atry(j).lt.amin(j)).and.(amin(j).gt.0.01)) then
            ibound = ibound + 1
         else if ((atry(j).gt.amax(j)).and.(amax(j).gt.0.01)) then
            ibound = ibound + 1
         endif
 20   continue

      return
      end





c $Id$
c $Author$
c $Date$
c $Source$
c $Revision$
c $Log$
c Revision 1.1  2007-02-20 23:10:22  sandeeps
c Initial revision
c
c Revision 1.11  2000/06/07 14:43:51  dmmatheu
c Works on henning molecules, but get weird numbers on vinyl + O2
c species ... reset tolerance to 0.01 as in previous code.  This sort of
c fixes problem ... can't quite get to old THERFIT numbers on
c no-integer-constraint option, but close (reasonable and good agreement
c in any case w/Cp data)
c
c Before putting in bounds-checking ...
c
c Revision 1.10  2000/06/02 16:14:47  dmmatheu
c A working integer-fit version!!! confirmed working nicely for ALL of
c Henning's test molecules.
c
c Revision 1.9  2000/06/02 01:49:57  dmmatheu
c intermediate version -- building integer loop.  Modularized inner
c marquardt loop.  Compiles.
c
c Revision 1.8  2000/06/02 00:50:49  dmmatheu
c Tests well against old algorithm for Henning's case.  Fixed all loop
c problems (hopefully).  Before turning mrqmin loop into separate
c function ...
c
c Revision 1.7  2000/05/31 20:46:35  dmmatheu
c loop termination criteria corrected
c
c Revision 1.6  2000/05/31 20:25:40  dmmatheu
c tested OK ... loop error from not setting ia(6) = 0!
c
c Revision 1.5  2000/05/31 18:49:45  dmmatheu
c fixed memory problem but now different answer than before !
c
c Revision 1.4  2000/05/31 16:07:24  dmmatheu
c solved loop problem w/raw integer ...
c
c Revision 1.3  2000/05/31 15:27:45  dmmatheu
c mystery loop problem -- loop commented out and works ok ... with loop
c gives 'singular matrix error' ????
c
c Revision 1.2  2000/05/30 00:01:40  dmmatheu
c checked against FOFX; OK
c
c Revision 1.1  2000/05/29 19:28:27  dmmatheu
c Initial revision
c
