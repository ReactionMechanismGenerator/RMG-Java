c
c
c
c*   now print out desnity of state info to separate file (filedef 31) - fort.31 removed by pey (30/3/04)
c***RKSSET***************************
      subroutine rksset(lout,iwell,Enot,nvibs,Edenom)
c    calculates Edenom following troe JPC 1979; also returns nvibs
c
c    PV modifications on 06/02/1995 for printing out density of states
c    corrected Whitten-Rabinovitch DOS code for integer and non-integer
c    cases  6/5/95
      implicit none
      include 'cdparams.fh'
      include 'cdisprop.fh'

      
c    local variables
      integer ifreq,iwell
      integer nvibs,lout
      integer ivib,ilev,nlev
      real*8 beta,w,Ezero,Enot,Edenom,sumdegen,sumf,sumf2,sumdif
      real*8 roe,rh,Er,Er0,rlfrq,Ebot,rfac
      character*1 tab
c    end declarations
	tab=char(9)
      
      sumdegen= 0.d0
      do 50 ifreq = 1,nfreqs(iwell) 
        sumdegen = sumdegen + degen(iwell,ifreq)
50    continue
      nvibs = int(sumdegen+0.5d0)
      write(lout,126) sumdegen, nvibs
126   format(1x, 'sum of THERM degeneracies = ',f7.2,' nvibs = ', i5)  
      sumf = 0.d0
      sumf2 = 0.d0 
      rlfrq = freq(iwell,1)
      do 100 ifreq = 1,nfreqs(iwell)
         if(rlfrq.gt.freq(iwell,ifreq))rlfrq = freq(iwell,ifreq)
         sumf = sumf + degen(iwell,ifreq)*freq(iwell,ifreq)
         sumf2 = sumf2 + degen(iwell,ifreq)*freq(iwell,ifreq)*
     2      freq(iwell,ifreq)
100   continue
      write(lout,123) rlfrq
 123  format('lowest frequency = ', f7.1, 'cm^-1')
c    removed factor: (sumdegen-1)/sumdegen <ayc 2/22/94>
      beta = (sumdegen-1.d0)*sumf2 / (sumf*sumf)
      Ezero = .5d0*2.859d-3*sumf

      if (Enot.gt.Ezero) then
         w = 10.d0**(-1.0506d0*(Enot/Ezero)**0.25d0)
      else
         w = 1.d0/(5.d0*(Enot/Ezero) + 2.73d0*dsqrt(Enot/Ezero) + 
     2      3.51d0)
      endif

      Edenom = Enot + (1.d0-beta*w)*Ezero

      write(lout,120) Enot,Ezero,nvibs,w,beta,Edenom
120   format(' Enot: ',f6.2,'  Ezero: ',f6.2,
     2   ' nsum: ',i3,'  w: ',f5.3,'  betav: ',f5.3,'  Edenom: ',f6.2)
c
      Er=Enot
      rh=2.859d-03
      roe=(Edenom**(sumdegen-1))/rfac(sumdegen-1)
      do 88 ifreq =1,nfreqs(iwell)
      roe=roe/((rh*freq(iwell,ifreq))**degen(iwell,ifreq))
88    continue
      write(lout,121)Er,tab,roe
121   format(1x,'At Energy Enot = ',f7.1,a1,
     #'Vibrational Density of States (Rho-Vib) = ',1pe10.3)
c	following line comment out by pey, 30/3/04
c      write(31, 1221) tab
cjmg fix
c1221 format (' Energy (kcal)'a1,'Density of States From Well Bottom')     
 1221 format (' Energy (kcal)',a1,'Density of States From Well Bottom')     
c
      Ebot=0.0d0
      nlev=1001
      do 999 ilev=1,nlev
      Er0=Ebot+(dfloat(ilev-1))*rlfrq*rh
      if(Er0.gt.300.0d0)return
      if (Er0.gt.Ezero) then
         w = 10.d0**(-1.0506d0*(Er0/Ezero)**0.25d0)
      else
         w = 1.d0/(5.d0*(Er0/Ezero) + 2.73d0*dsqrt(Er0/Ezero) +
     2      3.51d0)
      endif
      Er = Er0 + (1.d0-beta*w)*Ezero
      roe=(Er**(sumdegen-1))/rfac(sumdegen-1)
      do 881 ifreq =1,nfreqs(iwell)
      roe=roe/((rh*freq(iwell,ifreq))**degen(iwell,ifreq))
881   continue
c	following line comment out by pey, 30/3/04
c      write(31,122)Er0,tab,roe 
999   continue
122   format(1x,f7.1,a1,1pe10.3)
c
      return
      end
c
      double precision function rfac(x)
      implicit double precision (a-h,o-z)
      y=x
      rfac=1.0d0
      b1=-0.577191652
      b2=0.988205891
      b3=-0.897056937
      b4=0.918206857
      b5=-0.756704078
      b6=0.482199394
      b7=-0.193527818
      b8=0.035868343
1     if((y.lt.1).and.(y.ge.0))goto 2
      rfac=rfac*y
      y=y-1.0d0 
      goto 1
2     rfac=rfac*(1+
     #     b1*y+
     #     b2*y*y+
     #     b3*y*y*y+
     #     b4*y*y*y*y+
     #     b5*y*y*y*y*y+
     #     b6*y*y*y*y*y*y+
     #     b7*y*y*y*y*y*y*y+
     #     b8*y*y*y*y*y*y*y*y) 
      return
      end

c*RKS*******************
      function rks(lout,old,T,nvibs,Enot,Edenom,iwell)

      implicit none
      include 'cdparams.fh'
      include 'cdisprop.fh'
      include 'cdcolls.fh'
c     dmm to add dEdcm control ...
      include 'cdcontrl.fh'
c    local variables
      real*8 rks   
      logical old
      integer icoll
      integer iwell,lout,nvibs
      real*8 T,redmass,avgsig,ekmean,pi,alpha,FkT,result,Enot,
     2   Edenom,FE,partial,omega,delta,delt1,diff12,Eabs,
     3   FkT0,FEeff,Teff,betacol
      real*8 calFE,TFE 
      data pi /3.14159265d0/
c    end declarations

c    (orginal method: fixed FE; alpha is a function of T)    
c    parameters FkT,FkT0,delt1&diff12 are used to calculate beta
      if (old) then
         FE = 1.15d0
         FkT = FE*1.987d-3*T
         FkT0 = FkT
         delt1 = 1.d0
         diff12 = 1.d0

c    (or calculate FE; nvibs,Enot and Edenom have been returned
c    previously by rksset; we also use constant alpha, not deltE)
      else         
         FE = calFE(nvibs,Edenom,T)       
c       we freeze out effective FE at 1d6; get Teff
         if (FE.gt.1.d6) then
            Teff = TFE(nvibs,Edenom,T,1.d6)
            FEeff = calFE(nvibs,Edenom,Teff)
         else
            Teff = T
            FEeff = FE
         endif
         
c       get delt1,diff12 using Teff and FEeff
         call caldelt(FEeff,Enot,Edenom,nvibs,Teff,delt1,diff12)
         FkT = FEeff*1.987d-3*Teff
         
c       here, alpha is calculated using room T values  
c       changed from 295 to 300 ayc 3/94
         FkT0 = calFE(nvibs,Edenom,300.d0)*1.987d-3*300.d0
      endif

c    now perform sum over all colliders
      result = 0.d0
      do 100 icoll = 1,ncolls

         avgsig = 0.5d-8*(sig(iwell)+sig2(icoll))
         redmass = rmass*rmass2(icoll) /(rmass+rmass2(icoll))
         ekmean = dsqrt(ek(iwell)*ek2(icoll))
     
c       (option to calculate beta; remember DeltaE is in cal not kcal)
c       see gilbert, luther & troe, ber. bunsenges 87 (1983), p. 169.

         if (.not.consbeta(icoll)) then
            Eabs = dabs(deltaE(icoll))/1.d3

c          added exponent     <ayc 3/94>

c     dmm 20010813 added this if-else structure to allow user to specify
c     deltaEdown directly, in cm-1, instead of needing to specify dEall
c     in cal/mol.  At high temperatures this can be important, since
c     conversion method usually used is good only at low T.  Note that
c     Albert believes deltaEd should be at 300 K; we should check
c     Gilbert's paper on this.

            if(.not.dEdcm) then
               alpha = (dsqrt(.25d0*Eabs*Eabs + Eabs*FkT0) + .5d0*Eabs)*
     2              (T/300.d0)**dEExp(iColl)
            else
               alpha = dabs(deltaE(icoll))/349.755
            endif
    
            delta = (alpha*delt1 + FkT*diff12) /(alpha+FkT)
            if (delta.lt.1.d-75) delta = 1.d-175
            betacol = (alpha*alpha)/((alpha+FkT)*(alpha+FkT)*delta)
            
c     constrain beta to be less than 1
            if (betacol.gt.1.d0) betacol = 1.d0
            
c     (or we use fixed beta from input)         
         else
            betacol = beta(icoll)
         endif         
           
c       (omega term: original method; Forst p. 185) 
         if (old) then
            omega = 2.708d0*(ekmean/T)**0.333d0
            
c       (or use Reid, Prausnitz & Poling, pg. 393)
         else
            omega = 1.16145d0*(T/ekmean)**(-0.14874d0) + 
     2         0.52487d0*dexp(-0.77320d0*T/ekmean) + 
     3         2.16178d0*dexp(-2.43787d0*T/ekmean) 
         endif 

         partial = betacol*6.023d23*pi*avgsig*avgsig*omega*
     2      dsqrt(8.0d0*8.314d7*T/(pi*redmass))
     
         write(lout,80) T,icoll,FE,betacol,omega,xcoll(icoll),partial
80       format(1x,f5.0,' K: #',i2,'  FE = ',1pe8.2,'  beta = ', 0pf5.3,
     2      '  omega = ',f4.2,'  rks = ',f4.2,' x ',1pe8.2)

         if (.not.consbeta(icoll))
     2      write(lout,50) FkT,alpha,delt1,1.d0/delta
50       format('    FEkT: ',1pe9.3,'  alpha: ',1pe9.3,'  delta1: ',
     2      1pe9.3,'  1/delta: ',1pe9.3) 

         result = result + xcoll(icoll)*partial
 100  continue
 
      rks = result

c  amamo changes 20010801 (commented out by pey, 25/3/04)
cjmg fix
c     open(6789, file='feparam')
c     write(6789, *) FE
c      open(68, file='feparam')
c      write(68, *) FE

      return
      end

c***CALFE****************************
      function calFE(nvibs,Edenom,T)
      implicit none
      real*8 calFE

c    local variables
      integer nvibs 
      real*8 Edenom,T,a
      real*8 rint 
c    end declarations

c    this is integral of Whitten-Rabinovitch with constant a(E), evaluated
c    at Enot; comes from troe JPC 1979; consistent with grad.& ryshik

      a = -1.d0/(1.987d-3*T)
      calFE = rint(Edenom,nvibs-1,a)

      return
      end

c***TFE*****************************
      function TFE(nvibs,Edenom,Tinit,FE)
c    returns T corresponding to FE value
c    Tinit = initial guess for T
      implicit none
      real*8 TFE

c    local variables
      integer nvibs 
      real*8 Edenom,Tinit,FE,Tlow,Thigh,Tguess,Ftest,Terr
      real*8 calFE
c    end declarations

c    error tolerance in K
      Terr = 10.

      Tlow = 0.d0
      Thigh = 1.d4
      if (calFE(nvibs,Edenom,Tinit-500.d0).lt.FE) Tlow = Tinit-500.d0
      if (calFE(nvibs,Edenom,Tinit+500.d0).gt.FE) Thigh = Tinit+500.d0
      Tguess = Tinit

10    continue
         Ftest = calFE(nvibs,Edenom,Tguess)
         if (Ftest.lt.FE) then
            Tlow = Tguess
         else
            Thigh = Tguess
         endif
c        write(20,*) 'T: ',Tguess,'  FE: ',Ftest
         Tguess = .5d0*(Thigh+Tlow)
      if ((Thigh-Tlow).gt.Terr) goto 10

c    should be good within Terr/2  
      TFE = Tguess

      return
      end

c***CALDELT**************************
      subroutine caldelt(fE,Enot,Eeff,nvibs,T,delt1,diff12)
c    returns parameters for calculating high T denominator for beta factor
      implicit none

c    local variables
      real*8 FE,Enot,Eeff,aEzero,delup,dello,diff12,difup,diflo,delt1,
     2   a1,T,Eoff
      integer nvibs
      real*8 rint,dfint
c    end declarations 

c    as with calFE we integrate Whitten-Rabinovitch with constant a(E),
c    evaluating at Enot

      aEzero = Eeff - Enot
      a1 = -1.d0/(1.987d-3*T)
      Eoff = Enot + aEzero

c    we try this difference formalism because we encoutered roundoff with
c    more straight-forward approach (unfortunately we still have errors
c    in extreme cases)
c    we have factored x**n exp(a*x) /a from rint and dfint
c    leaving sums 1 + ...
c    note -dello is the normalization

      delup =  rint(Enot+aEzero,nvibs-1,a1)
      dello =  rint(aEzero,nvibs-1,a1)
      difup = dfint(Enot+aEzero,nvibs-1,a1,FE,Eoff)
      diflo = dfint(aEzero,nvibs-1,a1,FE,Eoff)
      
      delt1 = dexp(a1*Enot + dfloat(nvibs-1)*
     2   dlog((Enot+aEzero)/aEzero))*delup/(-dello) + 1.d0
      diff12 = (dexp(a1*Enot + dfloat(nvibs-1)*
     2   dlog((Enot+aEzero)/aEzero))*difup - diflo)/(-dello)

c     write(20,*) ' '
c     write(20,*) 'delup:  ',delup,'  dello:  ',dello 
c     write(20,*) 'difup:  ',difup,'  diflo:  ',diflo 
     
      return
      end
      
c***RINT******************************
      function rint(x,n,a)
c    evaluates integral of x^n exp (ax)
c    from grad & ryshik
      implicit none
      real*8 rint

c    local variables
      real*8 x,a,sum
      integer n,k
      real*8 facrat
c    end declarations
      
      sum = 0.d0    
      do 10 k = 0,n
         sum = sum + facrat(n,n-k) * (-a*x)**(-k)
10    continue

c    we factor out x**n exp(a*x) /a
      rint = sum
      
      return
      end
      
c***DFINT******************************
      function dfint(x,n,a,FE,offset)
c    evaluates differnce of two integrals of x^n exp (ax)
c    (see rint) - this routine is necessary to avoid roundoff errors
c    a2 = a*(1-1/FE)
      implicit none
      real*8 dfint

c    local variables
      real*8 x,a,FE,sum,offset
      integer n,k
      real*8 facrat
c    end declarations

      sum = 0.d0
      do 10 k = 0,n
         sum = sum + facrat(n,n-k) * (-a*x)**(-k) *        
     2      (1.d0 - (1.d0-1.d0/FE)**(-(k+1))*dexp(-a*(x-offset)/FE))
10    continue
     
c    above exponential term should reduce to dexp(0) in upper
c    limit of call and dexp(-Enot/(FEKT)) in lower limit      

c    we leave out common factor x**n exp(a*x) /a
      dfint = sum
      
      return
      end
     
c***FACRAT****************************
      function facrat(n,k)
c    returns n!/k! if n > k; else returns 1.

      implicit none
      real*8 facrat

      integer i,n,k
      real*8 prod

      prod = 1.d0
      if (n.gt.k) then
         do 10 i=k+1,n
            prod = prod*dfloat(i)
10       continue
      endif

      facrat = prod

      return
      end

