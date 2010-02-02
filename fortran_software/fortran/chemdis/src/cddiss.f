c ******COMPDISS******************
      subroutine compDISS(lecho,ldiag,dE) 

      implicit none
      include 'cdparams.fh'
      include 'cdwell0.fh'
      include 'cdwell1.fh'
      include 'cdwell2.fh' 
      include 'cdisprop.fh'
      include 'cdcolls.fh'
      include 'cdrange.fh'
      include 'cdrates.fh'
      include 'cdcontrl.fh'
      include 'cdASATempPres.fh'

c     dmm 20000811
c     for printing out pops and well names we include labels for now
      include 'cdlabels.fh'
c     dmm 20000811

      integer mxSteps
      parameter (mxSteps = 1000)
      
c   local variables
      logical divert(mxWells)
      logical rkovar,inAcces(mxWells),noExit
      integer it,ip,iWell,iProd
      integer nEint,lecho,ldiag,nTdone,nElast(mxTpts),imostd,nvibs,
     2   ierr,ipfill
      real*8 E,fofEnum,Emin,rNoExNm,Edenom,Emax,dE,Estart,ACsave,
     2   Eoff,expEoff,rhoNum
      real*8 Qpart(mxTpts),stab(mxWells,mxTpts),ACpop(mxWells),
     2   ACarray(mxWells,mxWells),rateP(mxWells,0:mxProds),
     3   rateI(mxWells,mxWells),B(mxWells)    
      real*8 rKofEI(mxWells,mxWells),rKofEP(mxWells,mxProds),
     2   ACDepI(mxWells),ACDepP(mxWells)
     
      real*8 rks,calRho,xNOS, calQ, tkinf, tksum, oldcalQ
      integer ipop
      real*8 tpop(mxWells),M
c tal 20000621
      equivalence (rkovar, rkold)
c tal 20000621
     
c******end declarations ***********

c    NOTE - this is just a modified version of compchem - however, our zero's 
c    are different here E = 0 from bottom of Well

c       now energy must be defined with reference to the bottom of each
c       Well - however variable E in loop is referenced with respect to
c       offset energy Estart - hopefully this will enable us to avoid
c       underflows of the exponential terms  - renormalization will
c       occur at end <ayc 3/93>

c    setup resultant matrix for solver; init noDExit flag
      do 5 iWell = 1,nWells
         B(iWell) = 0.0d0
         noDExit(iWell) = .false.
 5    continue
      B(inpWell) = -1.d0
      
c     initialize density of states data block
      tksum = 0.0D0
c     dmm 20000510 commented out
c     call initRho
c     dmm 20000510
      
c     tal 20000621 adding switch
      if (oldrho) call initRho
c     tal 20000621

c     dmm 20000908 attempt to fix ACpop init bug ... zero out acpop
      call initCDcalc(ACpop, RateP)
c     dmm 20000908

c    find max barrier height in system, deepest Well depth
c    find minimum barrier to first exit (Product) channel
      imostd = 1
      Emax = 0.d0
      Emin = Exmax(inpWell)
      do 10 iWell = 1,nWells      
         if (idpWell(iWell).gt.idpWell(imostd)) imostd = iWell         
         Emax = dmax1(Exmax(iWell) + EWell(iWell),Emax)
         do 10 iProd = 1, NProds(iWell)
            Emin = dmin1(EbarLP(iWell) + EProd(iWell,iProd),Emin)
 10   continue
         
c     we start at minimum exit channel
      
c     dmm 20000510 commented out for new xNOS number-of-states
c     calculation.  No longer necessary.
c      Estart = Exmin(inpWell)
c dmm 20000510

c tal 20000621 adding switch
      if (oldrho) then
         Estart = Exmin(inpWell) + .5d0*dE
      else
         Estart = Exmin(inpWell)
      end if
c tal 20000621



c    IMPORTANT: for the low p limit, this subroutine must switch if
c    calculating a stabilization channel and there are no Product
c    channels at that E - this results in a different P-dependence
c    than the high E terms

c    we try adding an offset energy to avoid overflows in our exponentials
      Eoff = .5d0*(Emax - Estart)
 
      write(lecho,15) Estart,Emax,dE,dE/2.859d-3
 15   format(/,1x,'DISSOC calculations from E of ',
     2     f6.2,' to at least ',f6.2,' kcals,',/,
     3     ' in ',f5.2,' kcal (',f6.1,' cm-1) steps.')       

c    initialize variables before E sum
      do 30 it = 1,nTemps
         nElast(it) = mxsteps
         do 30 iWell = 1,nWells
            do 30 iProd = 0,nProds(iWell)
               do 30 ip = -1,nPres+1
                  RK(iWell,iProd,it,ip) = 0.0d0
30    continue 

c***********************************************************************
c dmm 20000510
c
c     Generate xrhoofn, or sum-of-states vector, for every well.  See
c     cdchem.f comments.  Skip if chemact calc has already been done
c     (assumes chemact calc will be done first if both are selected)


c      if(.not.chemact) then
c         call genRhos((eMax/2.859D-3),(dE/2.859D-3), temp(nTemps))
c      endif
c dmm 20000510


c tal 20000621
c      adding switch
      if ((.not.chemact).and.(.not.oldrho)) then
         call genRhos((eMax/2.859D-3),(dE/2.859D-3), temp(nTemps))
      endif
c tal 20000621





c***********************************************************************

c    calculate stabilization rate which is a function of T only
c    parameter rkovar will set up constant FE, original omega calculation
c    we use lowest E barrier for Well as Enot
c    we redirect output for stabilization parameters to ldiag

      do 35 iWell = 1,nWells
         write(ldiag,*) 'Calculating stabilization rate parameters for',
     2        ' Well #',iWell
         if (.not.rkovar) call rksset(ldiag,iWell,Exmin(iWell),
     2        nvibs,Edenom)
         do 35 it = 1,nTemps
            stab(iWell,it) = rks(ldiag,rkovar,temp(it),nvibs,
     2           Exmin(iWell),Edenom,iWell)
35    continue

c    compute partition function vs T
      do 37 it = 1,nTemps
c tal 20000621 for old vs. new chemdis
         if (oldrho) then
            Qpart(it) = oldcalQ(temp(it),inpWell)
         else
            Qpart(it) = calQ(temp(it),inpWell,dE)
         end if
c tal 20000621
 37   continue
        
      nEint = - 1
      nTdone = 0
      
c****start E loop************
c    loop is terminated when nTdone = nTemps
c    E is kcal units and is referenced from bottom of inpWell
c    note each step represents energy between E-dE/2 and E+dE/2  
c    nEint is now the total # of iterations     


c***********************************************************************
c     20000809 tal for printing out population vectors
c      open(23, file='ac.dat',status='unknown')
c     20000809 tal
c
c     20000809 dmm for printing out detailed balance terms
c      open(24, file='detail.txt', status = 'unknown')
c      write(24,27) temp(1), (ISname(ipop), ipop=1,nwells)
c 27   format (F8.2, 1x, 8(A10, 3x))
c***********************************************************************

c     tal 20000815 - for new rmin for dissoc
      rdiss = 0.d0
c     tal 20000815

 40   continue 
      nEint = nEint + 1
      E = Estart + nEint*dE

c     get density of states for numerator
c     rhoNum = calRho(E,dE,inpWell,0)  dmm 20000510
c     rhoNum =  xNOS(E,dE,inpWell)
c     tal 20000621
      if (oldrho) then
         rhonum = calRho(E,dE,inpWell,0)
      else
         rhoNum =  xNOS(E,dE,inpWell)
      endif
c tal 20000621


c       compute all k(E)'s at this E
      call getKofEMatrix(E,dE,rKofEI,rKofEP)
        
c       we use EbarLP to flag inaccessible channels at this E
      do 330 iWell = 1,NWells
         if (E.le.EbarLP(iWell)) then
            inAcces(iWell) = .true.
         else
            inAcces(iWell) = .false.
         endif
 330  continue

c       start T loop
      do 500 it = 1,nTemps
         
c     variable nElast starts at mxsteps but is downsized
c     as determined directly below
         
         if (nEint.le.nElast(it)) then
            
c     calculate new f(E,T)  
c     if using old rho method, do not offset by half-width for
c     exponential; if using new method, need offset in exponential
c     term
            if (oldrho) then
               fofEnum = rhoNum / Qpart(it)*
     $              dexp(-(E-Eoff)/(1.987d-3*temp(it)))
            else
c     calculate new f(E,T) with half-width ... 
               fofEnum = rhoNum / Qpart(it)*
     2              dexp(-(E + 0.5*dE -Eoff)/(1.987d-3*temp(it)))
            endif
c***********************************************************************
c     dmm 20000623
c     TEMPORARY CALQ debugging lines ... these lines evaluate the HP
c     limit for dissociation and dump results to fort.91.  tksum will
c     be the integrated estimate of the HP limit for dissociation fed to
c     the code in fort.10.  It should not be higher than the limit!
c     
c     
c     
c      if(.not.oldrho) then
c         tkinf = Aprod(1,1)*xNOS(E-Ewell(1)-Eprod(1,1),dE
c     $        ,inpWell)*dE*exp(-(E+0.5)/(1.987d-3*temp(it)))*
c     $        (temp(it)**rNprod(1,1))/Qpart(it)
c         
c         
c         tksum = tksum + tkinf
c         write(91, 111) temp(it), tkinf, tksum
c 111     format(F8.2, 1X, E16.9, 1X, E16.9)
c      endif
c     dmm 20000623
c***********************************************************************

c     setup T-dependent portion of RateI and RateP
            call fillAr1(it,rKofEI,RKofEP,RateI,RateP,
     2           ACDepI,ACDepP)   
            
c     start P loop   
            do 400 ip = 0,nPres+1
               
c     complete activated complex matrix AC and RateP
c     note M is factored in limits		 
c     also we do an incredible kludge in fillAr so the
c     low P dissoc limit works out for noExit conditions
               call fillAr2(it,ip,ACDepP,ACDepI,stab,ACarray,
     2              RateI,RateP,inAcces,noExit,rNoExNm)     
               
c     this is for debugging matrix solver
c     if ((nEint.eq.0).or.(nEint.eq.1)) then
c     write(lecho,43) ip,it,ACarray(1,1),ACarray(1,2)
c     write(lecho,43) ip,it,ACarray(2,1),ACarray(2,2)
c     43                  format(1x,i2,' ',i2,' ',1pe14.7,' ',1pe14.7)
c     endif

c***********************************************************************
c     20000810 dmm code to print detailed balance terms
c               if((ip.gt.0).and.(ip.lt.nPres+1)) then
c                  do 101 ipop = 1, nwells
c                     tpop(ipop) = RateI(inpWell, ipop)*xNOS(E, dE
c     $                    ,inpWell)
c 101              continue
c                  write(24,102) E, (tpop(ipop), ipop=1, nWells)
c 102              format(F6.2, 1x, 8(E12.5, 1x))
c               endif
c     20000811 dmm
c***********************************************************************

		  
c     now we are almost ready to solve matrix
c     hopefully, we've eliminated all singularities
               call solve(NWells,ACarray,B,idpWell,
     2               idpWell(imostd),ACpop,inAcces,ierr,lecho,E)


c***********************************************************************
c     20000809 tal code to print out population vectors (commented out by pey, 25/3/05)
c
c               if((ip.gt.0).and.(ip.lt.nPres+1)) then
c                  expEoff = dexp(-Eoff/(1.987d-3*temp(it)))
c                  do 111 ipop = 1, nwells   
c
c                     tpop(ipop) = ACpop(ipop)*fofenum*expEoff
c     $                    *stab(inpWell,it)*pres(it,ip)/(82.1d0*temp(it)
c     $                    )
c 111              continue
c                  ipop =1
c                  write(23,112) E, (tpop(ipop),ipop=1,nwells)
c 112              format(F6.2, 1x, 8(E12.5, 1x))
c               endif
c***********************************************************************
c     20000809 tal code end

c     here's our kludge!!! recognizing the solution
c     vector goes as one over our very small M factor, we
c     multiply it out by same factor; also trip noDExit
c     flag and divert flag for latter

               if ((ip.eq.0).and.(noExit)) then                     
                  do 340 iWell = 1,NWells
                     ACpop(iWell) = ACpop(iWell) * rNoExNm
                     if (.not.inAcces(iWell)) then
                        noDExit(iWell) = .true.
                        divert(iWell) = .true.
                     else
                        divert(iWell) = .false.
                     endif
 340              continue
                  
               else
                  do 345 iWell = 1,NWells
                     divert(iWell) = .false.
 345              continue
               endif

c                we save parameter for convergence criterion     
               if (ip.eq.nPres+1) ACsave = ACpop(iwcmax)
               
c     more matrix debugging stuff		  
c     if ((nEint.eq.0).or.(nEint.eq.1))
c     2               write(lecho,*) ACpop(1),ACpop(2)
               
               if (ierr.eq.1) then
                  write(lecho,*) ' ' 
                  write(lecho,*) 'Failure at nEint, E, it, ip: ',
     2                 nEint,E,it,ip
                  stop
               endif
               
c     finally compute contribution to RK and keep
c     running sum over all levels
               do 350 iWell = 1,nWells
                  do 350 iProd = 0,nProds(iWell)
                     ipfill = ip
                     if ((ip.eq.0).and.divert(iWell).and.
     2                    (iProd.eq.0)) ipfill = -1
                     RK(iWell,iProd,it,ipfill) = RK(iWell,iProd,
     2                     it,ipfill) + fofEnum*RateP(iWell,iProd)
     3                    *ACpop(iWell)
 350           continue                                        
 400        continue
c  end P loop

c     tal 20000815
c     setting up new criteria for dissoc ASA cases
            if (ASA) then
               M = pres(icurtemp(1),icurpres(1)) / (82.1 * 
     2              temp(icurtemp(1)))
               rdiss = rdiss + fofenum*stab(inpwell,it)*M*
     2			dexp(-Eoff/(1.987d-3*temp(it)))
            end if
c     tal 20000815

c             check to see if this is our last sum vs T                  
c             don't need to sum any more terms if E is greater
c             and E is greater then Emax and there is convergence on
c             highest Product channel - we believe high P limit is
c             slowest to converge

            if (((E-1.d0*dE.gt.Emax).and.(RK(iwcmax,ipcmax,it,
     2           nPres+1).gt.1.d-175)
     3           .and.(fofEnum*RateP(iwcmax,ipcmax)*ACsave/
     4           RK(iwcmax,ipcmax,it,nPres+1).lt.2.d-3)).or.
     5           (nEint.ge.mxsteps)) then
               nElast(it) = nEint
               nTdone = nTdone + 1
            endif
            
c     when N = nElast(T) there will be no more iterations
c     on this T; when nTdone = nTemps we are done
                   
         endif           
 500  continue
c end T loop    
      if (nTdone.lt.nTemps) goto 40    
c*****E loop done
         close(23)
c    finally we normalize RK by stabilization rate and account for offset
c    we also factor M from the high and low p limits  
      do 700 it = 1,nTemps
         expEoff = dexp(-Eoff/(1.987d-3*temp(it)))                
                 
         do 700 iWell = 1,nWells
            do 700 iProd = 0,nProds(iWell) 
              do 700 ip = -1,nPres+1 	      
                  if ((ip.le.0).or.(ip.eq.nPres+1)) then
                     RK(iWell,iProd,it,ip) =  RK(iWell,iProd,it,ip)
     2                  *expEoff*stab(inpWell,it)
                  else
                     RK(iWell,iProd,it,ip) = 
     2                  RK(iWell,iProd,it,ip)*expEoff*stab(inpWell,it)
     3                  * pres(it,ip) / (82.1d0*temp(it))

c***********************************************************************
c     dmm 20000823 print out ks*M (temporary)
                     
                     if ((iWell.eq.1).and.(iProd.eq.1)) then
                        write(*,701) temp(it), pres(it,ip), stab(inpWell
     $                       ,it)*pres(it,ip)/(82.1d0*temp(it))
 701                    format('T = ', F8.2, 1x,'P = ', F8.2, 1x,
     $                       'ks[M] = ', E13.3)

                     endif
c***********************************************************************

                  endif
                  if (RK(iWell,iProd,it,ip).lt.1.d-200) 
     2               RK(iWell,iProd,it,ip) = 1.d-200
700   continue

c    echo out some integration parameters
      write(lecho,710)
710   format(/,1x,'INTEGRATION STEP RANGE',
     2   /,1x,' T (K) ',2x,'# steps',2x,'Estart (kcal) Eend (kcal)')

      do 760 it = 1,nTemps
         write(lecho,750) temp(it),nElast(it)+1,Estart,
     2      Estart+dE*nElast(it)
750      format(1x,f7.1,4x,i3,4x,2x,f7.2,4x,f7.2)
760   continue



      return
      end
                  

c $Author$
c $Date$
c $Source$
c $Revision$
c $Log$
c Revision 1.1  2007-02-20 23:10:23  sandeeps
c Initial revision
c
c Revision 1.11  2000/08/18 13:11:06  dmmatheu
c after merge w/tlada code ... changes include pop vector labelling and
c printing, detailed balance labeling and printing
c
c Revision 1.10  2000/08/09 18:20:48  dmmatheu
c after tlada changes to print out population vectors
c
c Revision 1.9  2000/06/26 14:00:50  dmmatheu
c successful merge with tlada code
c
c Revision 1.8  2000/06/23 17:08:14  dmmatheu
c after merge with Tom Lada code
c
c Revision 1.7  2000/06/23 17:00:08  dmmatheu
c before merge with Tom Lada code
c
c Revision 1.6  2000/06/22 21:59:43  dmmatheu
c dissoc side fixed -- Q consistent with fofEnum and fofEnum fixed ...
c OK!
c
c Revision 1.5  2000/06/22 17:10:46  dmmatheu
c Note:  Before fixing "greater rate constant at high pressure than HP
c limit rate constant" in dissoc problem.  Basically this comes from the
c fact that fofEnum is evaluated with a dE that's usually too big (1
c kcal) so that either the exponential term must either be evaluated halfway
c through (at E + 0.5*dE) or the partition function must be evaluated
c "erroneously" (e.g. with grainsize dE = 1 kcal), or even both.  It's
c not clear what's best.
c
c Revision 1.4  2000/05/11 15:36:47  dmmatheu
c Doesn't work well ... final results ~2 OOM higher than they should
c be.  Likely due to use of calQ and Qpart which may be inconsistent
c with our xNOS.  There is some serious unit disagreement ugliness here
c -- I can't seem to get to the bottom of it!  Someday may need to iron
c out calRho/xNOS differences once and for all.  For now try hack
c solution (running denom approach as in cdchem.f)
c
c Revision 1.3  2000/05/10 19:20:43  dmmatheu
c before altering calRhos to put in xNOS calls ...
c







