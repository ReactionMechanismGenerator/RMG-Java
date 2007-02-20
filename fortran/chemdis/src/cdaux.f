c *****getKofEMatrix******************
c    fills rate matricies

      subroutine getKofEMatrix(E,dE,rKofEI,rKofEP)
      implicit none
      
c    get header files
      include 'cdparams.fh'
      include 'cdwell0.fh'
      include 'cdwell1.fh'
      include 'cdwell2.fh'
      include 'cdisprop.fh'
      include 'cdcontrl.fh'
      
c   local variables
      integer iWell,ioWell,iProd
      real*8 E,dE,rhoDenom
      real*8 rKofEI(mxWells,mxWells),rKofEP(mxWells,mxProds)
      real*8 calRho, xNOS
c    end declaratons

c    note we kludge it if the ratio goes awry 

      do 100 iWell = 1,nWells
      
c       make sure bottom of well is above E
         if (E.gt.EWell(iWell)) then    

c dmm 20000404 replace calRho with xNOS	 
c            rhoDenom = calRho(E-EWell(iWell),dE,iWell,0)	    
c            rhoDenom = xNOS(E-Ewell(iWell), dE, iWell)

c tal 20000621 switch
             if (oldrho) then
                rhoDenom = calRho(E-EWell(iWell),dE,iWell,0)
             else
                rhoDenom = xNOS(E-Ewell(iWell), dE, iWell)
             end if
c tal 20000621

c          first, do Product channels
            do 60 iProd = 1,nProds(iWell)


c tal 20000621
             if (oldrho) then
               rKofEP(iWell,iProd) = 
     2            dmin1(calRho(E-EWell(iWell)-EProd(iWell,iProd),dE,
     3               iWell,nWells+iProd) / rhoDenom, 1.d0)
             else
               rKofEP(iWell,iProd) = 
     2            dmin1(xNOS(E-EWell(iWell)-EProd(iWell,iProd),dE,
     3               iWell) / rhoDenom, 1.d0)
             end if
c tal 20000621

60           continue 
                               
c          second, do isomerization channels
            do 80 ioWell = 1, nWells 
               if (Aisom(iWell,ioWell).gt.1.d-175) then

c tal 20000621
                 if (oldrho) then
                  rKofEI(iWell,ioWell) = 
     2               dmin1(calRho(E-EWell(iWell)-EIsom(iWell,ioWell),dE,
     3                  iWell,ioWell) / rhoDenom, 1.d0)
                 else
                  rKofEI(iWell,ioWell) = 
     2               dmin1(xNOS(E-EWell(iWell)-EIsom(iWell,ioWell),dE,
     3                  iWell) / rhoDenom, 1.d0)
                 end if
               else
                  rKofEI(iWell,ioWell) = 0.d0
               endif
c temporary k(E) printing line (commented out by pey, 25/3/04
c               write(99,81) iWell, ioWell, E, EWell(iWell), Eisom(iWell
c     $              ,ioWell), rKofEI(iWell,ioWell)
80          continue
 81         format(I2,1x,I2,1x,F7.2,1x, F7.2, 1x, F7.2, 1x, E13.4)
c     dmm 20000404 end replacement of calRho with xNOS

c      zero out for well whose depths are above E
         else
            do 90 iProd = 1,nProds(iwell)
               rkofEP(iWell,iProd) = 0.d0
90          continue

            do 95 ioWell = 1,nWells
               rkofEI(iWell,ioWell) = 0.d0
95          continue
         endif
	 
100   continue

      return
      end
       
c ******FILLAR1******************
c    fills rate vectors

      subroutine fillAr1(it,rKofEI,RKofEP,RateI,RateP,ACDepI,ACDepP)   
      implicit none
      
c    get header files
      include 'cdparams.fh'
      include 'cdwell0.fh'
      include 'cdwell1.fh'
      include 'cdwell2.fh'
c     include 'cdisprop.fh'
      include 'cdrange.fh'
      
c   local variables
      integer iWell,ioWell,iProd,it
      real*8 rateP(mxWells,0:mxProds),rateI(mxWells,mxWells) 
      real*8 rKofEI(mxWells,mxWells),rKofEP(mxWells,mxProds),
     2   ACDepI(mxWells),ACDepP(mxWells)

c    end declaratons

      do 100 iWell = 1,nWells
         ACdepI(iWell) = 0.d0
         ACdepP(iWell) = 0.d0
               
c       first, do Product channels
         do 60 iProd = 1,nProds(iWell)
            RateP(iWell,iProd) = Aprod(iWell,iProd)*
     2         temp(it)**rNProd(iWell,iProd)*
     3         dexp(-alphaProd(iWell,iProd)*temp(it))*
     4         rKofEP(iWell,iProd)
     
            ACdepP(iWell) = ACdepP(iWell) +
     2        RateP(iWell,iProd)
60        continue 
                               
c       second, do isomerization channels
         do 80 ioWell = 1, nWells 
            if (Aisom(iWell,ioWell).gt.1.d-175) then
               RateI(iWell,ioWell) = AIsom(iWell,ioWell)*
     2            temp(it)**rNIsom(iWell,ioWell)*
     3            dexp(-alphaIsom(iWell,ioWell)*temp(it))*
     4            rKofEI(iWell,ioWell)
     
               ACdepI(iWell) = ACdepI(iWell) +
     2            RateI(iWell,ioWell)
            else
               RateI(iWell,ioWell) = 0.d0
            endif
80       continue
100   continue

      return
      end
       
c ******FILLAR2******************
c    fills ACarray

      subroutine fillAr2(it,ip,ACDepP,ACDepI,stab,ACarray,
     2   RateI,RateP,noAcces,noExit,Rleak)
     
      implicit none  
      
c    get header files
      include 'cdparams.fh'
      include 'cdwell0.fh'
      include 'cdwell1.fh'
      include 'cdwell2.fh'
      include 'cdisprop.fh'
      include 'cdrange.fh'
      
c   local variables
      logical noAcces(mxWells),noexit
      integer iWell,ioWell,it,ip
      real*8 ACdepI(mxWells),ACdepP(mxWells),rateP(mxWells,0:mxProds),
     2   rateI(mxWells,mxWells),ACarray(mxWells,mxWells),
     3   stab(mxWells,mxTpts),Rleak
c    end declaratons

c    we provide some leakage current of M for the low P limit
c    hopefully, this will help avoid singularities; compare this
c    to typical M's:  P(atm) /(80*T)
      Rleak = 1.d-70
      noexit = .false.

c    first case (ip=0) is low p limit
      if (ip.eq.0) then                                    
         do 130 iWell = 1,nWells                                                            
c          first write input channels       
            do 120 ioWell = 1,nWells
               ACarray(iWell,ioWell) = RateI(ioWell,iWell)
120         continue
                        
c          then write depletion channel but
c          don't count stabilization rate for low p limit
            RateP(iWell,0) = stab(iWell,it)
            ACarray(iWell,iWell) =  -ACdepP(iWell) -ACdepI(iWell)
     2          -Rleak*stab(iWell,it)
130      continue

c    second case is for general p
      else if (ip.le.npres) then                 
         do 160 iWell = 1,nWells                                       
c          first write input channels            
            do 150 ioWell = 1,nWells
                ACarray(iWell,ioWell) = RateI(ioWell,iWell)
150         continue

c          then write depletion channel
            RateP(iWell,0) = stab(iWell,it) * pres(it,ip) /
     2         (82.1d0*temp(it))
            ACarray(iWell,iWell) =  -ACdepP(iWell) -ACdepI(iWell)
     2         - RateP(iWell,0)
160      continue

c       third case (ip=npres+1) is h p limit, here only count 
c       input channels from Wells of one less depth
         else
            do 190 iWell = 1,nWells                     
               do 180 ioWell = 1,nWells
                   if (idpWell(ioWell).eq.idpWell(iWell)-1) then
                      ACarray(iWell,ioWell) = RateI(ioWell,iWell)
                   else
                      ACarray(iWell,ioWell) = 0.d0
                   endif
180             continue
                  
c           in high P limit deactivation channel domimates
             RateP(iWell,0) = stab(iWell,it)
             ACarray(iWell,iWell) = -RateP(iWell,0)
190       continue
       endif 
       
c     NEW!!! since we can make arbitrary linear combination of 
c     our equations, we change the equation corresponding to inpWell
c     to be the sum of all equations - the right hand side does not
c     change because all other entries are zero; on the matrix side,
c     all isomerization rates must cancel, leaving only the stabilization
c     and Product channels.  we put these in explicitly.
c     this should make things better behaved at low pressures


c IMPORTANT:  we put an amazing kludge to make the low pressure rate
c             come out in the low P case for dissoc - if there are 
c             no exit channels the stabilization rates are the only 
c             things across the inpWell line - in this case we can't
c             put in zero but we put in a small number times stab
c             and flag it;
c             above, we must cancel it 

       if (ip.eq.0) then
          noexit = .true.
          do 200 iWell = 1,nWells
             if ((ACdepP(iWell).gt.1.d-175).and.(.not.noAcces(iWell)))
     2          noexit = .false.
200       continue

          if (noexit) then
             do 210 iWell = 1,nWells
                ACarray(inpWell,iWell) = -Rleak*stab(iWell,it)
210          continue
          else
             do 220 iWell = 1,nWells
                ACarray(inpWell,iWell) =  -ACdepP(iWell)
     2             -Rleak*stab(iWell,it)
220          continue
          endif

       else if (ip.le.npres) then
          do 250 iWell = 1,nWells
             ACarray(inpWell,iWell) =  -ACdepP(iWell) -RateP(iWell,0)
250       continue
      endif
       
c    we leave hp limit alone - this matrix should be okay by itself,
c    plus it doesn't satisfy the same summation condition anyway

      return
      end

c***SOLVE*******************
      subroutine solve(N,A,B,id,Nd,X,elim,ierr,lecho, E)
c    id corresponds to idpWell; we need this for elimination procedure
c    Nd is largest value of id
c    IMPORTANT: parameter mxWells must be identical to that of main

      implicit none
      include 'cdparams.fh'
      
      real rNegThresh
      parameter (rNegThresh = 1.d-10)
      
c    local variables
      integer i,j
      integer N,Nsubmit,index,jndex,id(mxWells),Nd,ierr,lecho,lout
      logical elim(mxWells),neg
      real*8 A(mxWells,mxWells),Asubmit(mxWells,mxWells),
     2   B(mxWells),X(mxWells),Bsubmit(mxWells),Xout(mxWells), E
c    end local variables

      neg = .false.

c temp debug line
c      write(*,*) 'mxWells in solve', mxWells

c    we use input vector elim to get rid of inactive channels
c    note this could still leave some singularities near threshold
c    in the low p limit where there may not be any exit channels
c    to an isomer

c    now eliminate rows and columns corresponding to inactive channels
      index = 0
      do 200 i = 1,N
         if (.not.elim(i)) then
            index = index+1
             Bsubmit(index) = B(i)
	     
             jndex = 0
             do 190 j = 1,N
               if (.not.elim(j)) then
                  jndex = jndex+1
                  Asubmit(index,jndex) = A(i,j)
               endif
190         continue

         endif	 
200   continue
      Nsubmit = index

c    now get solution to Asubmit * Xout = Bsubmit

      if (Nsubmit.gt.0) 
     2   call getB(Nsubmit,Asubmit,Bsubmit,Xout,ierr,lout, E)
      
c    reconstruct solution, note X's have already been zeroed for default
c    NOTE:  now - we try something - if any of the X's come out really
c    negative (more than rNegThresh: this is unphysical and we zero them
c    all out; else we just zero out that channel
    
      index = 0
      do 300 i=1,N
         if (.not.elim(i)) then
            index = index+1
	    if (Xout(index).lt.-rNegThresh) neg = .true.
            X(i) = Xout(index)
            if (X(i).lt.0.d0) X(i) = 0.d0
         endif
300   continue

      if (neg) then
         do 350 i = 1,N
            X(i) = 0.d0
350      continue
      endif

      return
      end

c****GETB1**********************************************
c    here, we attempt to link to gauss elim routines
c    from Numerical recipes, by press et al.
      subroutine getB(N,A,B,X,ierr,lout, E)
      
      implicit none     
      include 'cdparams.fh'
      include 'cdwell2.fh'
      include 'cdwell0.fh'
            
c    local variables
      integer i
      integer N,ierr,lout
      real*8 A(mxWells,mxWells),B(mxWells),X(mxWells), E
      character*8 name
c    end local variables

c***********************************************************************
c     20000809 tlada code for printing out A arrays ... commented out
c     for now ...
c      name = char(int((E - Exmin(inpwell))/10)+48) //
c     2       char(int(mod((E-Exmin(inpwell)),10)) + 48)
c      open (78, file=name, status='unknown')



c      write(78,10) E
c 10   format('E = ', D10)
c      if (N .eq. 1) then
c         write(78,25) A(1,1)
c      else if (N .eq. 2) then
c         write(78,26) A(1,1),A(1,2),A(2,1),A(2,2)
c      else if (N.eq. 3) then
c         write (78,27) A(1,1),A(1,2),A(1,3),A(2,1),A(2,2),A(2,3),
c     2        A(3,1),A(3,2),A(3,3)
c      else 
c         write (78,28) A(1,1),A(1,2),A(1,3),A(1,4),A(2,1),A(2,2),
c     2        A(2,3),A(2,4),A(3,1),A(3,2),A(3,3),A(3,4),A(4,1),A(4,2),
c     3        A(4,3),A(4,4)
c      end if
c 25   format('A = [', D10, ']')
c 26   format('A = [', D10, D10, ' ; ', D10, D10, ']')
c 27   format('A = [',D10,D10,D10,';',D10,D10,D10,';',D10,D10,D10,']')
c 28   format('A = [',D10,d10,d10,d10,';',d10,d10,d10,d10,';',
c     2     d10,d10,d10,d10,';',d10,d10,d10,d10,']')
c***********************************************************************

      ierr = 0

c    get gaussian elimin solver
      call gaussj(A,N,mxWells,B,1,1)

c    return X's
      do 300 i=1,N
         X(i) = B(i)
300   continue

c***********************************************************************
c     20000809 tlada code for printing matrices ...
c      if (N.eq.1) then
c         write(78,305) B(1)
c      else if (N.eq.2) then
c         write(78,306) B(1),B(2)
c         else if (N.eq.3) then
c            write(78,307) B(1),B(2),B(3)
c            else 
c               write(78,308) B(1),B(2),B(3), B(4)
c      end if
c 305  format('B = [',D10,']')
c 306  format('B = [',D10,';',D10,']')
c 307  format('B = [',D10,';',D10,';',D10,']')
c 308  format('B= [',D10,';',D10,';',D10,';',D10,']')
c      close(78)
c***********************************************************************

      return
      end


c****GETB2**********************************************
c    here we attempt link to ldu routines 
c    from NUMERICAL RECIPES by Press et al.
      subroutine getB2(N,A,B,X,ierr,lout)
      
      implicit none          
      include 'cdparams.fh'
      
c    local variables
      integer i,j
      integer N,ierr,lout,Ndim,Idum(mxWells)
      real*8 A(mxWells,mxWells),Adum(mxWells,mxWells),Bdum(mxWells),
     2   B(mxWells),X(mxWells),Vdum(mxWells),D
c    end local variables

c    we use dummy arrays to leave originals uncorrupted
c    we initialize D and Vdum just to surpress warnings

      Ndim = mxWells
      D = 0.d0
      ierr = 0
      
c    initialize dummy arrays for input matrices
      do 260 i = 1,N
	 Vdum(i) = 0.d0
         Bdum(i) = B(i)
         do 260 j = 1,N
            Adum(i,j) = A(i,j)
260   continue

c    get LU decomposed Adum
      call ludcmp(Adum,N,Ndim,Idum,D,Vdum,ierr,lout)

c    proceed if no errors encountered
      if(ierr.eq.0) then     
c       get solution for X
         call lubksb(Adum,N,Ndim,Idum,Bdum)
      
c      reconstruct solution   
         do 300 i = 1,N
            X(i) = Bdum(i)
300      continue
      endif

      return
      end

c***REMOVED GETB3 for successful linking on SGI <ayc 1/94>

c***CLEANUP**********************
      subroutine cleanup(option)      
      implicit none

c    get header files      
      include 'cdparams.fh'
      include 'cdwell0.fh'
      include 'cdwell1.fh'
      include 'cdwell2.fh'
      include 'cdrates.fh'

c    local variables
      character*8 option
      integer iWell,iProd
c    end declarations 

c       get integer powers for lindemann forms 
c       (depends on depth of Well and whether 
c       stabilization or Product channel)
            
      if (option.eq.'Chemact') then
            
         do 500 iWell = 1,nWells
            do 500 iProd = 0,nProds(iWell)
            
               nPr(iWell,iProd) = idpWell(iWell)
                       
               if (iProd.eq.0) then
                  nRKinf(iWell,iProd) = 1 - idpWell(iWell)
               else
                  nRKinf(iWell,iProd) = -idpWell(iWell)
               endif
500      continue
   
      else              
c       all dissoc channels  
   
         do 600 iWell = 1,nWells
            do 600 iProd = 0,nProds(iWell)
            
               nPr(iWell,iProd) = idpWell(iWell)
                       
               if (iProd.eq.0) then
                  nRKinf(iWell,iProd) = 2 - idpWell(iWell)
               else
                  nRKinf(iWell,iProd) = 1 - idpWell(iWell)
               endif
600      continue
         
      endif
       
      return
      end

c***********************************************************************
c***********************************************************************
c dmm 20000908
c     initCDcalc
c
c     After debugging we have found that certain arrays need to be 0'd
c     out at the beginning of each CHEMDIS run.  This function does that
c     .  Called from cddiss.f, cdchem.f
c
c***********************************************************************
c***********************************************************************

      subroutine initCDcalc(ACpop, RateP)

      implicit none
      include 'cdparams.fh'
      include 'cdrates.fh'
      
      real*8 ACpop(mxWells), RateP(mxWells, 0:mxProds)
      integer iWell, iProd, it, ip

      do 10 iWell = 1, mxWells
         ACpop(iWell) = 0.0D0
         do 20 iProd = 0, mxProds
            RateP(iWell, iProd) = 0.0D0
c$$$            do 30 it = 1, mxTPts
c$$$               do 40 ip = -1, mxPPts+1
c$$$                  RK(iWell, iProd, it, ip) = 0.0D0
c$$$ 40            continue
c$$$ 30         continue
 20      continue
 10   continue

      return 
      end

c----------------------------------------------------------------------c
c *** RRKM 1 k(E) Thread *** START                                     c
c----------------------------------------------------------------------c
c***********************************************************************
c***********************************************************************
c     dmm 20011022
c     
c     kofERRKM
c
c     This function uses the sum of states stored in xrhoofn as 'well'
c     iTS and the sum of states for iWell to compute a simple version
c     of the RRKM k(E) for printing to the Barker file.
c
c     Function evalues the density of states for the denominator as the
c     average of DOS over Evr + Ebar to Evr + Ebar + Egrain.  Sum of
c     states evaluated for the arithmetic middle of the grain ... 
c***********************************************************************
c***********************************************************************

      subroutine kofERRKM(Evr, Ebar, Egrain, iWell, iTS, xkofE, xkrev)

      implicit none

      include 'cdparams.fh'
      include 'straightBS.fh'
      include 'cdnos.fh'
      include 'cdwell0.fh'
      include 'cdwell1.fh'

      real*8 xNOScm, xkrev, rhode2, Ebar2
      real*8 Evr, Ebar, Egrain, xkofE, hplanc, pdgn, sumTS, rhodenom
      integer iWell, iTS, nETS

c     Degeneracy is needed here, since this would normally be embedded
c     in the A-factor for QRRK but is never available explicitly

      pfflag = .true.
      pdgn = dfloat(idgRRKM)
      hplanc = (6.6260755D-34/(1000.0*4.18))*349.755*6.02D+23
      rhodenom = xNOScm( (Evr + Ebar),Egrain,i1RRKM)/Egrain
c     dmm 20011102  added +1 offset to nETS
      nETS = nint((Evr + 0.5*Egrain)/bsgs) + 1
      sumTS = xrhoofn(nETS, nWells+1)
      if(sumTS.le.0.0D0) then
         sumTS = 1.0D0
      endif
      xkofE = pdgn*sumTS/(hplanc*rhodenom)
      pfflag = .false.

c***********************************************************************
c detailed balance k(E) calculation
c***********************************************************************

      Ebar2 = Eisom(i2RRKM,i1RRKM)*349.755
      rhode2 = xNOScm( (Evr + Ebar2), Egrain, i2RRKM)/Egrain
      
      xkrev = xkofE*rhodenom/rhode2
c      write(91,100) Evr, xkrev
 100   format(1x, F9.2, 1x, E16.9)
      return
      end
c----------------------------------------------------------------------c
c *** RRKM 1 k(E) Thread *** END                                       c
c----------------------------------------------------------------------c

c $Id$
c $Author$
c $Date$
c $Source$
c $Revision$
c $Log$
c Revision 1.1  2007-02-20 23:10:23  sandeeps
c Initial revision
c
c Revision 1.8  2001/11/21 03:16:45  dmmatheu
c after RRKM changes (successful?); before adding reverse rate constant feature
c
c Revision 1.7  2001/10/22 17:55:44  dmmatheu
c before adding RRKM k(E) function
c
c Revision 1.6  2001/09/20 18:03:19  dmmatheu
c before modifying initCDcalc to see whether we can do without
c time-consuming zeroing of the very large RK matrix
c
c Revision 1.5  2000/08/09 16:57:29  dmmatheu
c tlada changes -- matrix printing code (for matlab-based solver
c debugging) is incomplete and here commented out (dmm)
c
c Revision 1.4  2000/06/23 15:16:51  dmmatheu
c test of VC after merging ...
c
c Revision 1.3  2000/06/23 14:59:39  dmmatheu
c before merging with Tom's code ... working debugged version
c
c Revision 1.2  2000/04/04 17:25:02  dmmatheu
c before replacing calRho with xNOS
c
c Revision 1.1  2000/04/04 17:24:43  dmmatheu
c Initial revision
c

     
