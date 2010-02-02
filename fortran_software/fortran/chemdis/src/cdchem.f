c*****12/12/97  more changes via Jeff Grenda
c*********12/9/97  working more on setting emax properly for case where all exit wells
c                   lower than entrance
C*********6/12/94--INITIALIZED eMax to be 10.0 kcal to try to handle case where all exit
c   channels lower than entrance.  (Otherwise eMax=0 and go into never-never land.)
c ******COMPCHEM******************
      subroutine compCHEM(lecho,ldiag,dE) 
 
      implicit none
      include 'cdparams.fh'
      include 'cdwell0.fh'
      include 'cdwell1.fh'
      include 'cdwell2.fh' 
      include 'cdisprop.fh'
      include 'cdcolls.fh'
      include 'cdrange.fh'
      include 'cdrates.fh'
c     tal 20000621 for density of states switching
      include 'cdcontrl.fh'
c     tal 20000621
      
      integer mxSteps
      real tarray(2)
      real tutime, tstime
      save tutime, tstime
      parameter (mxSteps = 2000)
      
c   local variables
      logical inAcces(mxWells),noExit
      integer it,ip,iWell,iProd
      integer nEint,lecho,ldiag,nTdone,nElast(mxTpts),imostd,nvibs,ierr
      real*8 E,fofEnum,RKin,Edenom,eMax,dE,Estart,ACsave,Eoff,rNoExNm,
     2   rhoNum
      real*8 denom(mxTpts),stab(mxWells,mxTpts),ACpop(mxWells),
     2   ACarray(mxWells,mxWells),rateP(mxWells,0:mxProds),
     3   rateI(mxWells,mxWells),B(mxWells)    
      real*8 rKofEI(mxWells,mxWells),rKofEP(mxWells,mxProds),
     2   ACDepI(mxWells),ACDepP(mxWells)
     
      real*8 rks,calRho, xNOS

     
c******end declarations ***********

c    NOTE - structure of this subroutine is determined by the fact that 
c    we want efficient code but don't want to have to store arrays vs E; 
c    hence E is the outside loop

c    setup resultant matrix for solver; also noDExit which will be
c    set to false but will never be changed for the chemact calc.
      tutime = 0.0
      tstime = 0.0
      do 5 iWell = 1,nWells
         B(iWell) = 0.0d0
         noDExit(iWell) = .false.
5     continue
      B(inpWell) = -1.d0

c    initialize density of states data block

c     20000509 dmm commented out
c      call initRho
c     20000509 dmm
c     20000620 tal added switch
      if (oldrho) call initRho
c     20000620

c     20000908 dmm bug fix -- zero out matrices
      call initCDcalc(ACpop, RateP)
c     20000908 dmm

c    find max barrier height in system, deepest Well depth
      imostd = 1
      eMax = 0.d0
      do 10 iWell = 1,nWells      
         if (idpWell(iWell).gt.idpWell(imostd)) imostd = iWell         
         eMax = dmax1(Exmax(iWell)+EWell(iWell),eMax)
10    continue
c*******checking emax to make sure it is >0.
         if (eMax.le.0.d0) then
             eMax= 10.
             write(lecho,155) eMax
 155         format(/,1x,'Set minimum integration limit above
     * entrance well: eMax = ',f7.2 )
c*********above approach still doesn't keep calc. from going to never never land
c  for time being will simply halt code with a message to add
c   a dummy high temperature channel
c        write(lecho, 154)
c154    format(/,1x,'Need to add a dummy high energy exit channel to get
c     *  calc. to work properly')
c       stop         
         else
         endif
c*********above approach still doesn't keep calc. from going to never never land
c  for time being will simply halt code witha message to add
c   a dummy high temperature channel


c    we try adding an offset energy to avoid overflows in our exponentials
c    note this is automatically cancelled in RK by denom
      Eoff = .5d0*eMax
 
c    we start at dE/2
c     dmm 20000404 commented out for new xNOS number-of-states
c     calculation.  No longer necessary.
c      Estart = .5d0*dE
c      Estart = 0.0D0
c     dmm 20000404

c     tal 20000620
c     changed to give option of xNOS vs. old calcs
      if (oldrho) then
         Estart = .5d0*dE
      else
         Estart = 0.0d0
      end if
c     tal 20000620




      write(lecho,15) Estart,eMax,dE,dE/2.859d-3
15    format(/,1x,'CHEMACT calculations from E of ',
     2   f6.2,' to at least ',f6.2,' kcals,',/,
     3   ' in ',f5.2,' kcal (',f6.1,' cm-1) steps.')       

c    initialize variables before E sum
      do 30 it = 1,nTemps
         denom(it) = 0.d0
         nElast(it) = mxsteps
         do 30 iWell = 1,nWells
            do 30 iProd = 0,nProds(iWell)
               do 30 ip = 0,nPres+1
                  RK(iWell,iProd,it,ip) = 0.0d0
30    continue 


c***********************************************************************
c dmm 20000330
c     Generate xrhoofn, or sum-of-states vector, for every well.  
c     NOTE:  As now formulated max temperature must be last temperature
c     (temps in ascending order)
c
c     timing routine
c     genRho expects units of cm-1
c      call dtime(tarray)


c tal 20000621
c     adding switch to make possible old and new CHEMDIS
c     call genRhos((eMax/2.859D-3),(dE/2.859D-3), temp(nTemps))
      if (.not.(oldrho)) then
      call genRhos((eMax/2.859D-3),(dE/2.859D-3), temp(nTemps))
      end if
c tal 20000621


c      call dtime(tarray)
c      write(*,949), tarray(1)
c      write(*,959), tarray(2)
 949  format(1X, 'BS Code User Time:    ', E16.9)
 959  format(1X, 'BS Code System Time:  ', E16.9)


c dmm 20000330
c***********************************************************************


c    calculate stabilization rate which is a function of T only
c    parameter rkold will set up constant FE, original omega calculation
c    we use lowest E barrier for Well as Enot
c    we redirect output for stabilization parameters to ldiag

      do 35 iWell = 1,nWells
         write(ldiag,*) 'Calculating stabilization rate parameters for',
     2      ' Well #',iWell
         if (.not.rkold) call rksset(ldiag,iWell,Exmin(iWell),
     2      nvibs,Edenom)
     
         do 35 it = 1,nTemps
            stab(iWell,it) = rks(ldiag,rkold,temp(it),nvibs,
     2         Exmin(iWell),Edenom,iWell)
35    continue
        
      nEint = - 1
      nTdone = 0
      
c****start E loop************
c    loop is terminated when nTdone = nTemps
c    E is kcal units and is referenced from input barrier height
c    note each step represents energy between E-dE/2 and E+dE/2  
c    nEint is now the total # of iterations     

40    continue 
         nEint = nEint + 1
         E = Estart + nEint*dE

c       get density of states for numerator; offset cancels due to def of 0

c     dmm 20000404 replace calRho with xNOS
c         rhoNum = calRho(E,dE,inpWell,nWells+inpChan)
c         rhoNum = xNOS(E,dE,inpWell)
c     dmm 20000404	

c     tal 20000621 add switch to facilitate dens. of states method choice
         if (oldrho) then
            rhoNum = calRho(E,dE,inpWell,nWells+inpChan)
         else
            rhoNum = xNOS(E,dE,inpWell)
         end if
c     tal 20000621



c       compute all k(E)'s at this E\
c     Begin timing of getKofEMatrix (dmm 20000219)

c         call dtime(tarray)
         call getKofEMatrix(E,dE,rKofEI,rKofEP)
c     End timing of getKofEMatrix (dmm 20000219)
c         call dtime(tarray)
c         write(*,909), tarray(1)
c         write(*,919), tarray(2)
c         tutime = tutime + tarray(1)
c         tstime = tstime + tarray(2)
c         write(*, 929), tutime
c         write(*, 939), tstime
c         tarray(1) = 0.0
c         tarray(2) = 0.0
 909  format(1X, 'getKofEMatrix User Time:    ', E16.9)
 919  format(1X, 'getKofEMatrix System Time:  ', E16.9)
 929  format(1X, '              Accum. User:  ', E16.9)
 939  format(1X, '              Accum. System:', E16.9)
c       we use EbarLP to flag inaccessible channels at this E
         do 330 iWell = 1,NWells
            if (E.le.EbarLP(iWell)) then
               inAcces(iWell) = .true.
            else
               inAcces(iWell) = .false.
            endif
330      continue

c       start T loop
         do 500 it = 1,nTemps
c*************12/12/97 via Jeff Grenda:
cjmg fix                                           
      if (mod(nEint-1,10).eq.0)           
     *write(6,203)'T = ',temp(it),'  E= ',E  
 203  format(1x,2(a,2x,f12.5))              
  
c          variable nElast starts at mxsteps but is downsized
c          as determined directly below

            if (nEint.le.nElast(it)) then

c     dmm 20000622
c             calculate new f(E,T)  
c     Use old method for old gamma-function DOS method, else evaluate
c     exponential at midpoint of width
      if (oldrho) then
         fofEnum = rhoNum*dexp(-(E-Eoff)/(1.987d-3*temp(it)))
      else
         fofEnum = rhoNum*dexp(-(E + 0.5*dE -Eoff)/(1.987d-3
     $        *temp(it)))
      endif

c             keep running sum on denom's     
               denom(it) = denom(it) + fofEnum

c             setup T-dependent portion of RateI and RateP
               call fillAr1(it,rKofEI,RKofEP,RateI,RateP,
     2            ACDepI,ACDepP)     

c             start P loop   
               do 400 ip = 0,nPres+1
         
c                complete activated complex matrix AC and RateP
c                note M is factored in limits
                     call fillAr2(it,ip,ACDepP,ACDepI,stab,ACarray,
     2                  RateI,RateP,inAcces,noExit,rNoExNm)     
		 
c                now we are almost ready to solve matrix
c                hopefully, we've eliminated all singularities


c                     call dtime(tarray)
                     call solve(NWells,ACarray,B,idpWell,
     2                    idpWell(imostd),ACpop,inAcces,ierr,lecho,0.0)
                     
c                     call dtime(tarray)
c                     write(*,909), tarray(1)
c                     write(*,919), tarray(2)
c                     tutime = tutime + tarray(1)
c                     tstime = tstime + tarray(2)
c                     write(*, 929), tutime
c                     write(*, 939), tstime
c                     tarray(1) = 0.0
c                     tarray(2) = 0.0


c                we save parameter for convergence criterion     
                  if (ip.eq.nPres+1) ACsave = ACpop(iwcmax)

                  if (ierr.eq.1) then
                     write(lecho,*) ' ' 
                     write(lecho,*) 'Failure at nEint, E, it, ip: ',
     2                  nEint,E,it,ip
                     stop
                  endif

c                finally compute contribution to RK and keep
c                running sum over all levels
                  do 350 iWell = 1,nWells
                     do 350 iProd = 0,nProds(iWell)
                        RK(iWell,iProd,it,ip) = RK(iWell,iProd,
     2                     it,ip) + fofEnum*RateP(iWell,iProd)
     3                     *ACpop(iWell)
350               continue                                        
400            continue

c             check to see if this is our last sum vs T                  
c             don't need to sum any more terms if E is greater
c             and E is greater then eMax and there is convergence on
c             highest Product channel - we believe high P limit is
c             slowest to converge

               if (((E-1.d0*dE.gt.eMax).and.
     2            (RK(iwcmax,ipcmax,it,nPres+1).gt.1.d-75)
     3            .and.(fofEnum*RateP(iwcmax,ipcmax)*ACsave/
     4            RK(iwcmax,ipcmax,it,nPres+1).lt.2.d-3)).or.
     5            (nEint.ge.mxsteps)) then
                     nElast(it) = nEint
                     nTdone = nTdone + 1
               endif
                  
c              when N = nElast(T) there will be no more iterations
c              on this T; when nTdone = nTemps we are done
                   
           endif           
500      continue

      if (nTdone.lt.nTemps) goto 40    
c*****E loop done
         
c    finally we normalize RK by input rate and account for denom's      
      do 700 it = 1,nTemps                 
         RKin = Ain*temp(it)**rNin*dexp(-alphaIn*temp(it))*
     2      dexp(-Ein/(1.987d-3*temp(it)))  
                 
         do 700 iWell = 1,nWells
            do 700 iProd = 0,nProds(iWell) 
               do 700 ip = 0,nPres+1     
                 RK(iWell,iProd,it,ip) = 
     2              RKin*RK(iWell,iProd,it,ip)/denom(it)
                 if (RK(iWell,iProd,it,ip).lt.1.d-200) 
     2              RK(iWell,iProd,it,ip) = 1.d-200
700   continue

c    echo out some integration parameters
      write(lecho,710)
cjmg fix
c710   format(/,1x,'INTEGRATION STEP RANGE',
c     2   /,x,' T (K) ',2x,'# steps',2x,'Estart (kcal) Eend (kcal)')
710   format(/,1x,'INTEGRATION STEP RANGE',
     2   /,1x,' T (K) ',2x,'# steps',2x,'Estart (kcal) Eend (kcal)')

      do 760 it = 1,nTemps
         write(lecho,750) temp(it),nElast(it)+1,Estart,
     2      Estart+dE*nElast(it)
750      format(1x,f7.1,4x,i3,4x,2x,f7.2,4x,f7.2)
760   continue

      return
      end
                  
c $Id$
c $Author$
c $Date$
c $Source$
c $Revision$
c $Log$
c Revision 1.1  2007-02-20 23:10:23  sandeeps
c Initial revision
c
c Revision 1.7  2000/06/23 16:54:28  dmmatheu
c After Tom Lada code merge
c
c Revision 1.6  2000/06/23 15:32:50  dmmatheu
c Before merge with Tom Lada code
c
c Revision 1.5  2000/05/09 13:41:57  dmmatheu
c Working version, used for walkthrough/May 4 example.  3kT rule
c adjusted.  Before commenting out initRho, updating cddiss.f chain to
c use xNOS instead of gamma functions.
c
c Revision 1.4  2000/04/04 21:51:19  dmmatheu
c Runs Without Crashing.  Tested OK for vinyl+O2.  Used 3kT rule for
c estimating nmax in sum of states vector calculations.
c
c Revision 1.3  2000/04/04 17:29:46  dmmatheu
c before replacement of calRho with xNOS
c
c Revision 1.2  2000/03/30 23:47:07  dmmatheu
c successful generation
c
c Revision 1.1  2000/02/19 20:28:13  dmmatheu
c Initial revision
c
c
c
