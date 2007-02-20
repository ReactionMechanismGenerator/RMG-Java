 
c***********************************************************************
c***********************************************************************
c dmm 20000330
c
c     genRhos subroutine
c
c     This subroutine calculates and loads values into the xrhoofn
c     tables, replacing the recursion and gamma-table based rhoOfN
c     calculations used previously.  
c     
c     The sums-of-states algorithm used, SBS in 'straightBS.f',
c     automatically convolves with any number of 1d-rotors, assuming
c     them classical after method of Astholz (see Gilbert and Smith).
c     IN this particular implementation, we can set all rotational
c     constants and symmetries to 1.0, since they will cancel (all
c     implementations are ratios).  Unit issues for the energies should
c     also cancel.
c
c     There is a lot of complicated machinery in here to deal with
c     noninteger degeneracies from the frequency fitting.  If we could
c     force integer frequencies, lots of this wouldn't be needed.
c
c     eMax and frequencies are expected in cm-1.  The grainsize should
c     be in cdnos.fh, its default set (10cm-1) or read in by cdgets.f
c
c     MODIFICATIONS
c     
c     dmm 20000627 
c     Adding 'wellscalcd' vector, which is initialized to
c     0.0 and stores the maximum temperature for which the xrhooofn
c     vector for a particular array was last calculated.  If new T is
c     greater than this T, recalc the well, otherwise, skip (should be
c     OK).  This should not affect the operation of the maxBSindex
c     vector, which lets us generate a warning when a call to xNOS
c     exceeds the energy levels to which we have calculated rho.
c
c***********************************************************************
c***********************************************************************

      subroutine genRhos(eMax, dE, temp)

      implicit real*8 (a-h, o-z)
      include 'cdparams.fh'
      include 'straightBS.fh'
      include 'cdnos.fh'
      include 'cdisprop.fh'
      include 'cdwell0.fh'
      include 'cdwell1.fh'
      include 'cdcontrl.fh'

      integer numfreq(mxWells, mxFreqs)
      integer iWell, Isize, imax1, ip1, iw1, nfrot, ndos
      real*8 VIB(maxmodes), BR(maxmodes), BDG(maxmodes)
      real*8 remaind(mxFreqs)
      real*8 Rcm 
      real*8 Egrain1, Emax2, Viblo, Ebarmax
      integer maxnmax = 0


c     Set partition function flag.  See cdnos.fh.  Will be problematic
c     here if NOS's are pre-calc'd somewhere else

      pfflag = .false.
c     

c     Boltz. constant cm-1/mol/K
      Rcm = 0.694995954

c     Set all BR, BDG to 1.0 as per discussion above
      do 40 i=1, maxmodes
         BR(i) = 20.0D0
c         BR(i) = 1.0D0
         BDG(i) = 1.0D0
 40   continue

c----------------------------------------------------------------------c
c *** RRKM 1 k(E) Thread *** START                                     c
c----------------------------------------------------------------------c
c     Loop over wells; go one more if RRKM option is turned on to get
c     the RRKM TS sum of states.  
      
      if(RRKM) then
         ndos = nWells + 1
      else
         ndos = nWells   
      endif
c----------------------------------------------------------------------c
c *** RRKM 1 k(E) Thread *** END                                       c
c----------------------------------------------------------------------c

      do 10 iWell = 1, ndos

c***********************************************************************
c     Begin If-statement block (goes all the way to '10 continue') which
c     is entered if the density of states for this well needs to be
c     calculated.
c***********************************************************************

c     dmm 20000908 this logic is screwed up and inpwonly actually isn't
c     available here anyway.

c     dmm 20010612 commented out below because I don't understand why '
c     inpwonly' logic is there ...
c$$$         if(temp.gt.wellscalcd(iWell) .or. (.not. inpwonly)) then
c$$$            wellscalcd(iWell) = temp
         if(temp.gt.wellscalcd(iWell)) then

c     set nfrot, number of free internal rotors, to 0
            nfrot = 0

         
c     2.  Find sum of degeneracies degensum
c     3.  Find the integer degeneracies for the BS
c     4.  Find the sum of the integer degeneracies ntotal
c     5.  Load the VIB vector for the BS algorithm
            
            
            sumdegen = 0.0
            ntotal = 0
            jcount = 1
            
            do 20 ifreq = 1, nFreqs(iWell)
                              
               sumdegen = sumdegen + degen(iWell, ifreq)
               numfreq(iWell, ifreq) = nint(degen(iWell, ifreq))
               remaind(ifreq) = degen(iWell, ifreq) -
     $              dfloat(numfreq(iWell,ifreq))
               
               ntotal = ntotal + numfreq(iWell, ifreq)

c     Avoid loading the VIB vector if the frequencies have been read
c     from a file
               if(.not.FrRtRd(iWell)) then
                  do 15 j = 1, numfreq(iWell, ifreq)
                     VIB(jcount) = freq(iWell, ifreq)
                     jcount = jcount + 1
 15               continue
               endif
 20         continue


c***********************************************************************
c     Cleanup:  Check whether nint(sumdegen) is the same as summation of
c     frequences ntotal.  If not (if ntotal is less) then add that
c     frequency whose degeneracy is closest to the next integer number
c     .This complicated machinery is only necessary if we plan on using
c     noninteger degeneracies ...
c***********************************************************************

c     It shouldn't matter, but avoid this block if we have read the
c     frequencies from a file ...

            if(.not.FrRtRd(iWell)) then
               if( nint(sumdegen - 0.01).gt.ntotal) then
                  rmax = 0.0D0
c     Look for largest remainder
                  do 30 ifreq = 1, nFreqs(iWell)
                     if (remaind(ifreq).gt.rmax) then
                        rmax = remaind(ifreq)
                        nmax = ifreq
                     endif
 30               continue
c     Add one more frequency corresponding to that mode with the
c     largest remainder (see RN II, 3/30/00)
                  VIB(jcount) = freq(iWell, nmax)
                  ntotal = ntotal + 1
               endif
            endif

c***********************************************************************
c     If frequencies and rotors have been read from file, set this up;
c     because other parts of code use degen vector, this must be set ...
c***********************************************************************
            if(FrRtRd(iWell)) then
               do 32 jcount = 1, nFreqs(iWell)
                  VIB(jcount) = freq(iWell,jcount)
                  degen(iWell,jcount) = 1.0
 32            continue
               sumdegen = dfloat(nFreqs(iWell))
               nfrot = n1drot(iWell)
               ntotal = nFreqs(iWell)

               
            endif

c----------------------------------------------------------------------c
c *** RRKM 1 k(E) Thread *** START                                     c
c----------------------------------------------------------------------c
c RRKM substitution lines -- use same as for well ...
c skip setting of maxBSindex if we are on the 'TS' -- last well and 
c RRKM is set ...
               if(RRKM.and.(iWell.eq.ndos)) then
                  maxBSindex(iWell) = maxBSindex(i1RRKM)
                  goto 31
               endif
c----------------------------------------------------------------------c
c *** RRKM 1 k(E) Thread *** END                                       c
c----------------------------------------------------------------------c


c***********************************************************************
c     Calculation of max energy level for SBS
c     
c     1.  Send along B = 1.0 for the 1d rotors since these technically
c     don't matter (they should cancel)  1.0*dE is added to agree with
c     end of compchem function in cdchem.f
c     
c     2.  For chemical activation:  go to 4.5 kT above entrance channel
c     ,or 1.5kT above eMax when 4.5kT above entrance does not exceed
c     eMax + 1.5kT
c     
c     3.  For dissoc:  go to 4.5 kT above well bottom; or 1.5 kT above
c     highest exit eMax if 4.5kT does not exceed eMax + 1.5 kT
c     
c     Note:  AYC is Modular Man!  No special handling needed for dissoc.
c     Emax will be from bottom of well, and Ewell(inpWell) = 0 for
c     dissoc case; everything works out the same

            xnatoms = (sumdegen + 6.0D0 + dfloat(nrot))/(3.0)
            
c     eDist = 2.0*(1.5*xnatoms*Rcm*temp) + eMax

c     dmm 20020429  After some problems running Preeti cases in XMG,
c     changed this to be 1/2 of eDist above eMax for criteria (Preeti
c     calc ends up exceeding nmax or maxBSindex)

            eDist = 6.0*(1.5*xnatoms*Rcm*temp)
            
            if ( eDist.ge.(eMax + 0.5*eDist)) then
               nmax = nint((eDist - Ewell(iWell)/2.859D-3)/bsgs) + 1
            else
               nmax = nint( (eMax + 0.5*eDist-Ewell(iWell)/2.859D-3)
     $              /bsgs)+ 1
            endif
            maxBSindex(iWell) = nmax     





c***********************************************************************
c     BARKER OUTPUT Setting of nmax ...
c***********************************************************************

c***********************************************************************
c     For printing Barker output, need nmax to be large enough to agree
c     with 'barker.inp' file ...
            if(barker) then
               open(111, file='barker.inp', status='old')
               read(111,*) Egrain1, imax1, Emax2, Isize, Viblo
               close(111)

c     Notes on 'nmax' for Barker code generation:  Barker code
c     multiwell-1.01 computes k(E) with E = 0 defined as TS ZPE.  This
c     means if calculation goes up to Emax2, then max rho necessary will
c     actually be Emax2 + Ebarmax (maximum barrier out of this well).
c     It's dumb since the actual wells are only calculated to Emax2, but
c     the Barker code (foolishly) follows this format with E=0 at TS
c     ZPE for the k(E) values as well.  These super-high k(E)'s can't
c     actually be used by multiwell (because the rho's only go up to
c     Emax, not Emax +Ebarmax!) but MW barfs w/o them. 
c     
c     Decided not to rewrite Barker code to fix this inefficiency, so
c     here lies much wasted CPU time until Barker fixes it himself ...
c     
c     Note small chance that nmax calc'd this way will not be large
c     enough for good CHEMACT or DISSOC calc, but in practice probably
c     OK.

c***********************************************************************
c     Loops to find max barrier in kcal/mol
c***********************************************************************

               Ebarmax = 0.0
c     First check product channels ...
               do 300 ip1 = 1, nProds(iWell)
                  if(EProd(iWell, ip1).gt.Ebarmax) then
                     Ebarmax = EProd(iWell, ip1)
                  endif
 300           continue

c     Now check isomer channels ...
               do 301 iw1 = 1, mxWells
                  if(EIsom(iWell, iw1).gt.Ebarmax) then
                     Ebarmax = Eisom(iWell, iw1)
                  endif
 301           continue

c***********************************************************************
c     End loops to find max barrier ... convert Ebarmax to kcals!
c***********************************************************************
               
               Emax2 = Emax2 + Ebarmax*349.775
               nmax = nint(Emax2/bsgs) + nint(2.0*(Emax2/dfloat(Isize))
     $              /bsgs) + 100

               maxBSindex(iWell) = nmax

            endif
            
c***********************************************************************
c     End BARKER OUTPUT setting of nmax
c***********************************************************************


c***********************************************************************

c     For printing code: find the max nmax
            if (nmax.gt.maxnmax) then
               maxnmax = nmax
            endif
c     Error checking
            if(nmax.gt.maxNOS) then
               write(*,*)
     $   ' ERROR::genRhos.f: Number of grains for sum of states greater'
               write(*,*)
     $   '                   than maxNOS.  Choose larger grainsize or '
               write(*,*)
     $   '                   increase maxNOS in straightBS.fh'

               write(*,*) 'maxNOS = ', maxNOS, ' nmax = ', nmax
               write(*,*) ' ... exiting ...'
               stop
            endif

c----------------------------------------------------------------------c
c *** RRKM 1 k(E) Thread *** START                                     c
c----------------------------------------------------------------------c
c     re-entry point below for case where iWell is the RRKM TS

 31         continue
c----------------------------------------------------------------------c
c *** RRKM 1 k(E) Thread *** END                                       c
c----------------------------------------------------------------------c

c----------------------------------------------------------------------c
c *** RRKM 1 k(E) Thread *** START                                     c
c----------------------------------------------------------------------c
c     Set the rotational constants BR(1) for the k-modes if these have
c     been specified (i.e. krot = 1)
            if(RRKM.and.(rotk(iWell).gt.0.0)) then
               BR(1) = rotk(iWell)
               BDG(1) = 1.0
            else
               BR(1) = 20.0
               BDG(1) = 0.0
            endif
c----------------------------------------------------------------------c
c *** RRKM 1 k(E) Thread *** END                                       c
c----------------------------------------------------------------------c

            call SBS(nmax, ntotal, VIB, (nRot + nfrot), BR,BDG,bsgs)
c     Copy results to xrhoofn
            do 60 i=1, nmax
               xrhoofn(i, iWell) = TBS(i)
 60         continue 
            
         endif
c***********************************************************************
c     End if-statement block for recalc of well, continue on to next
c     well 
c***********************************************************************

 10   continue

c***********************************************************************
c Printing results for testing (comment out for normal usage)
c
c$$$      open(100, FILE='genRho.out')
c$$$      do 110 i=1, maxnmax
c$$$         write(100, 120) dfloat(i-1)*bsgs, (xrhoofn(i,iWell), iWell=1,
c$$$     $        nWells)
c$$$
c$$$
c$$$
c$$$
c$$$ 110  continue
c$$$ 120  format(F9.1, 6(1X,5E16.9,/))
c***********************************************************************


c***********************************************************************
c Print results for Barker MEQ density file
c
      if(barker) then
         call PrintBarkerDens(maxnmax)
      endif
c***********************************************************************



      return
      end


c***********************************************************************
c dmm 20000404
c
c     FUNCTION xNOS
c
c     This function evaluates the number of states in a width dE given
c     that xrhoofn arrays have been previously evaluated.  
c
c     Note that this algorithm EXPECTS bsgs to be considerably lower than 
c     dE.  If it is not, it issues a warning but ploughs ahead, at the
c     user's peril.  The results become a bit unpredictable for bsgs
c     larger than the dE grain, jsut as for master equation cases. 
c
c     E   = energy in kcal/mol above ZPE of molecule
c     dE  = width, in kcal/mol
c     iWell = well number
c***********************************************************************

      function xNOS(energy, width, iWell)
      
      implicit real*8 (a-h, o-z)
      include 'cdparams.fh'
      include 'straightBS.fh'
      include 'cdnos.fh'

      integer nE1, nE2

c     Set xNOS to 0.0 if called with a negative energy.  This is natural
c     when coming from getKofEmatrix (signifies no reaction poss. to
c     that well/product at energy)
     
      if(energy.lt.0.D0) then
         xNOS = 0.0D0
         return
      endif

c     Convert energy and width into cm-1; convert into grainsize index

c     dmm 20011102  Added offset of 1 to nE1 and nE2 ... probably will
c     not affect most results but this is more accurate as xrhoofn
c     begins as xrhoofn(1) at energy level 0.0 cm-1

      xen = energy*349.755
      xwd = width*349.755
      nE1 = nint(xen/bsgs) + 1
      nE2 = nint( (xen + xwd)/bsgs) + 1
c pey 15aug04
c      write(*,*) 'xnos = ', ne1, ne2, maxBSindex(iWell)
c pey
c     Warning check: issue if dE < 2*bsgs and not calculating Q.
      if((xwd.lt.2.0*bsgs).and.(.not.pfflag)) then
         write(*,10)
         write(*,20)
         write(*,30)
      endif

 10   format(/,'WARNING::genRhos.f:xNOS: width for calculation order')
 20   format('           of bsgs. Reset bsgs to lower value or use ')
 30   format('           larger dE.')
c     Warning if nE2 exceeds max. energy of xrhoofn for this well
      if(nE2.gt.maxBSindex(iWell)) then
         write(*,40)
         write(*,50)
         write(*,60)
         write(*,70), iWell, nE2, maxBSindex(iWell)
c         pause
      endif
 40   format(/,'WARNING::genRhos.f:xNOS: xrhoofn index higher than max')
 50   format('           calculated index for this well.  Results will')
 60   format('           be unpredictable.  Use larger nmax.')
 70   format(/,'iWell = ',I4, ' nE2 = ', I10, ' maxBSindex = ', I10)
 
c     Return xNOS
      xNOS = xrhoofn(nE2, iWell) - xrhoofn(nE1, iWell)


      return
      end
c dmm 2000404
c***********************************************************************


c***********************************************************************
c***********************************************************************
c     FUNCTION xNOScm
c
c     Same as xNOS but assumes energy in cm-1
c***********************************************************************
c***********************************************************************

      function xNOScm(energy, width, iWell)
      
      implicit real*8 (a-h, o-z)
      include 'cdparams.fh'
      include 'straightBS.fh'
      include 'cdnos.fh'

      integer nE1, nE2

c     Set xNOS to 0.0 if called with a negative energy.  This is natural
c     when coming from getKofEmatrix (signifies no reaction poss. to
c     that well/product at energy)
     
      if(energy.lt.0.D0) then
         xNOS = 0.0D0
         return
      endif

c     Convert energy and width into cm-1; convert into grainsize index

c     dmm 20011102  See above for explanation of offset on nE1 and nE2
      xen = energy
      xwd = width
      nE1 = nint(xen/bsgs) + 1
      nE2 = nint( (xen + xwd)/bsgs) + 1

c     Warning check: issue if dE < 2*bsgs and not calculating Q.
      if((xwd.lt.2.0*bsgs).and.(.not.pfflag)) then
         write(*,10)
         write(*,20)
         write(*,30)
      endif

 10   format(/,'WARNING::genRhos.f:xNOS: width for calculation order')
 20   format('           of bsgs. Reset bsgs to lower value or use ')
 30   format('           larger dE.')
c     Warning if nE2 exceeds max. energy of xrhoofn for this well
      if(nE2.gt.maxBSindex(iWell)) then
         write(*,40)
         write(*,50)
         write(*,60)
         write(*,70), iWell, nE2, maxBSindex(iWell)
      endif
 40   format(/,'WARNING::genRhos.f:xNOS: xrhoofn index higher than max')
 50   format('           calculated index for this well.  Results will')
 60   format('           be unpredictable.  Use larger nmax.')
 70   format(/,'iWell = ',I4, ' nE2 = ', I10, ' maxBSindex = ', I10)
 
c     Return xNOS
      xNOScm = xrhoofn(nE2, iWell) - xrhoofn(nE1, iWell)


      return
      end
c dmm 2000404
c***********************************************************************



c***********************************************************************
c***********************************************************************
c dmm 20000627
c
c     SUBROUTINE InitGenRhos
c
c     All this subroutine does is zero out the wellscalcd vector (see
c     genRhos, cdnos.fh).  It is probably not necessary, but it helps
c     make sure.  Called from chemdis2.f ... it's clunky but it
c     works for now.  A more robust and elegant solution escapes me. 
c     
c***********************************************************************
c***********************************************************************

      subroutine InitGenRhos

      implicit none
      include 'cdparams.fh'
      include 'straightBS.fh'
      include 'cdnos.fh'

      integer iWell

      do 10 iWell = 1, mxWells
         wellscalcd(iWell) = 0.0D0
 10   continue

      return
      end

c***********************************************************************
c***********************************************************************
c dmm 20011216
c
c     SUBROUTINE WRSumDens
c
c     This function computes the Whitten Rabinovitch approx to the sums
c     and densities of states; note that Albert & orig. authors have
c     their own version.  I don't trust its units consistency with
c     direct count, hence this function.  Plan is for this to be called
c     above the "break" energy of the double-array density of states
c     files produced by PrintBarkerDens for the multiwell code.
c
c     PLACEHOLDER:  possible serious issue at higher T and P by using WR
c     approx for DOS/SOS but 3-freq direct count for QRRK k(E)s.
c     Probably better to use the smoother WR where possible anyway,
c     since 3-freq approx is only meaningful in a "smooth" sense, and
c     spikiness/noise resulting from direct count is artificial.
c     
c     Evr:  energy in cm-1 above ZPE
c     dens:  density of states returned in 1/cm-1 (molecular)
c
c***********************************************************************
c***********************************************************************

      SUBROUTINE WRSumDens(iWell, Evr, Sum, Dens)

      implicit none

      integer iWell, i, j, ndg
      real*8 Evr, Sum, Dens, piv, Eprime, Ezero
      real*8 vmsq, vsqm, vfsqsum, vfsum, betaR, xnrot, xnvib, Qrot, Bk,
     $     Bdegen, wfac, afrac, lprodf, lWRsum, gamarg, dwdEp, lWRdens

      double precision gammln

      include 'cdparams.fh'
      include 'cdisprop.fh'

      piv = 3.14159

c***********************************************************************
c     Begin constants evaluation -- these do not change given same
c     molecule
c***********************************************************************

c     Calculate zero-point energy

      Ezero = 0.0
      vfsum = 0.0
      vfsqsum = 0.0
      lprodf = 0.0
      xnvib = 0.0
      do 10 i = 1, nFreqs(iWell)
         ndg = nint(degen(iWell,i))
         do 15 j = 1, ndg
            xnvib = xnvib + 1.0
            Ezero = Ezero + 0.5*freq(iWell,i)
            vfsum = vfsum + freq(iWell,i)
            vfsqsum = vfsqsum + (freq(iWell,i))**(2.0)
            lprodf = lprodf + LOG(freq(iWell,i))
 15      continue
 10   continue

c     Total (float) number of rotational modes
      xnrot = dfloat(nRot + n1drot(iWell))

c     Calculate vmsq, vsqm
      vsqm = vfsqsum/xnvib
      vmsq = (vfsum/xnvib)**(2.0)

c     betaR for the WR approx
      betaR = ((xnvib-1.0)/(xnvib))*((xnvib+0.5*xnrot)/(xnvib))*(
     $    vsqm/vmsq)

c     Evaluate rotational partition function: k-mode
c     Note: rotational constant default here of 20.0 MUST agree with that
c     set above in genRhos!
c     PLACEHOLDER:  there is currently no handling of rotational mode
c     degeneracies.  These are all assumed 1 by the code (should be
c     changed at some point) since QRRK uses only ratios of DOS anyway,
c     where the degeneracy would cancel out.

      Qrot = 1.0
      if(nRot.gt.0) then
         if(rotk(iWell).le.0.0) then
            Bk = 20.0
            Bdegen = 1.0
         else
            Bk = rotk(iWell)
            Bdegen = 1.0
         endif
         Qrot = Qrot*((piv/(Bk*Bdegen))**(0.5))
      endif

c     Evaluate other rotors.  PLACEHOLDER:  Again, rotational modes are never
c     actually read in for internal rotors, nor stored, since the QRRK
c     k(E) formula causes them to cancel along with degeneracies.  This
c     is probably a deficiency which should be corrected.  Constant
c     here should be same as the default set around line 70 of genRhos.f
     
      if(n1drot(iWell).gt.0) then
         do 20 i = 1, n1drot(iWell)
            Bk = 20.0
            Bdegen = 1.0
            Qrot = Qrot*((piv/(Bk*Bdegen))**(0.5))
 20      continue
      endif

c***********************************************************************
c End constants evaluation
c***********************************************************************

c     Evaluate Eprime
      Eprime = Evr/Ezero
c     Test function limits
      if(Eprime.lt.0.1) then
         write(*,900) 
         stop
      else if(Eprime.gt.8.0) then
         write(*,901)
      endif
c     Calculate wfac
      if((Eprime.lt.1.0)) then
         wfac = (5.0*Eprime + 2.73*(Eprime**0.5) + 3.51)**(-1.0)
         dwdEp = -(5.0 + 1.365*(Eprime**(-0.5)))*(wfac**2.0)
      else
         wfac = EXP(-2.4191*(Eprime**0.25))
         dwdEp = -(0.60478*(Eprime**-0.75))*wfac
      endif
c     Calculate afrac
      afrac = 1.0 - betaR*wfac
c     Calculate gammln function argument
      gamarg = xnvib + 0.5*xnrot + 1.0
c     Calculate log of WR sum of states and dens of states
      lWRsum = LOG(Qrot) - gammln(gamarg) + (gamarg-1.0)*LOG(Evr + afrac
     $     *Ezero) - lprodf
      Sum = EXP(lWRsum)

      gamarg = xnvib+0.5*xnrot
      lWRdens = (LOG(Qrot) - gammln(gamarg) + (gamarg-1.0)*LOG(Evr +
     $     afrac*Ezero) - lprodf) + LOG(1.0 - betaR*dwdEp)

      Dens = EXP(lWRdens)


 900  format(1x,'ERROR:  WR called with energy too low for accuracy',1x,
     $    ' ... stopping')
 901  format(1x
     $     ,'WARNING:  WR called with energy greater than 8 times ZPE')

      return
      end
     
      
c***********************************************************************
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
c Revision 1.16  2001/12/16 20:44:45  dmmatheu
c before adding Whitten-Rabinovitch subroutine, to get smooth DOS/SOS
c profile above high energy for Barker file input ...
c
c Revision 1.15  2001/12/02 17:51:06  dmmatheu
c Before offset fix:  it seems that nE1 and nE2 should be incremented by
c 1 automatically, since xrhoofn(1) is the same as energy level 0.  Not
c sure how this error got past for so long.  This gives "look-ahead" dos
c in xNOScm:  call at energy, width gives number of states from energy
c to energy+width.
c
c Revision 1.14  2001/10/19 16:02:04  dmmatheu
c before changes to generate SOS for the 1 RRKM rate constant
c
c Revision 1.13  2001/06/13 14:12:16  dmmatheu
c After changes to turn Barker file outputs on/off, major debug of
c mxWells and mxProds (zeroing loops of MakeShadow were exceeding array
c bounds).  Also added xNOScm from ASA-MasterEqn.. directory ... merging
c to one chemdis file
c
c Revision 1.12  2000/08/09 18:37:31  dmmatheu
c after tlada bug fix on dissoc
c
c Revision 1.11  2000/06/27 22:12:16  dmmatheu
c Changed to avoid recalculating density of states for a given well,
c when it's not necessary.  Corresponding changes to chemdis2.f,
c cdnos.fh.  Checked on standard vinyl+O2 9P14T example, looks good (no
c ASA attempted).
c
c Revision 1.10  2000/06/23 17:54:34  dmmatheu
c before Tlada merge ...
c
c Revision 1.9  2000/06/20 18:47:28  dmmatheu
c Bug in previous version never discovered ... things seem to be working
c fine now (??).  Very strange.
c
c Revision 1.8  2000/06/20 18:18:47  dmmatheu
c Before working on very strange bug -- commenting out the genRho.out
c file production lines causes different results in fort.65 for dissoc
c case
c
c Revision 1.7  2000/05/09 18:05:30  dmmatheu
c Works for vinyl + O2 example from May 4.
c Before changes with dissoc implementation.
c
c Revision 1.6  2000/04/04 21:13:50  dmmatheu
c Applied 3kT rule to estimate max for density of states ... need to
c build safety AND warning for exceeding this into code ...
c
c Revision 1.5  2000/04/04 20:38:51  dmmatheu
c Code compiles and runs without core-dumps.  eMax/max height for genRho
c is not handled correctly.  Need to estimate eMax.
c
c Revision 1.4  2000/04/04 15:12:04  dmmatheu
c Check in before adding xNOS function
c
c Revision 1.3  2000/03/30 23:43:18  dmmatheu
c Successful implementation in CHEMDIS on vinyl + O2 (generation of
c xrhoofn only)
c
c Revision 1.2  2000/03/30 20:32:23  dmmatheu
c pre-compiled version
c





