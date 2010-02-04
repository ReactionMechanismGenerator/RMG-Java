c***RKLIM**************************
c    fun to generate limiting rate constants from the kFit common block
      function RKlim(t,index)
      
      implicit none
      include 'cdparams.fh'
      include 'cdlimfit.fh'

c    local variables
      integer index,icoef
      real*8 t,basis(mxKfit),RKfln
      real*8 RKlim
c    end declarations

      RKfln = 0.d0
      call fPoly(t,basis,Nkfit)	    
      do 50 icoef = 1,Nkfit
         RKfln = Rkfln + CoefLim(icoef,index)*basis(icoef)
50    continue

      RKlim = dexp(Rkfln)
      
      return
      end
      
c***FPOLY**************************
c    extern basis function for polynomial lsq fit to limiting rates
      subroutine fPoly(x,VecCoefs,nCoefs)
      
      implicit none

      integer nCoefs,i    
      real*8 x,VecCoefs(*),var

      var = dsqrt(1.d3/x)    
      VecCoefs(1) = 1.d0
      
      do 10 i = 2,nCoefs
         VecCoefs(i) = VecCoefs(i-1)*var
10    continue
      
      return
      end      
 
c***KLIMFIT************************
c    this sub does all the 4+ parameter lsq fitting of
c    the high and low p limits
      subroutine kLimFit(option,lecho,lkout,iwell,iprod)
      
      implicit none
      include 'cdparams.fh'
      include 'cdwell0.fh'
      include 'cdlabels.fh'
      include 'cdrange.fh'
      include 'cdrates.fh'
      include 'cdlimfit.fh'

c    local variables
      character option*8,RCname*20
      integer it,index,icoef,istart,ipout
      integer lkout,lecho,iwell,iprod
      real*8 RKlog(mxTpts),RKfit,error,avgerr,errpar,
     2   Coefs(mxKfit)
      real*8 RKlim
      character tab,str*2,fitID*80
c    end declarations

      tab = char(9)
      
c    we either call fPoly once to get Nkfit or we set it here
      Nkfit = 5
c     call fPoly(1000.d0,basis,Nkfit)    

c    set ID string for fit
      fitID = 
     2   'Fit to log K with 5 coefs (x^0,..,x^4) in (1000/T)^1/2'
      
c    set high and low p exponents (0 = lp, 1 = hp; -1 = special
c    add'l dissoc stab term)
      nlim(0) = nRKinf(iwell,iprod) + nPr(iwell,iprod)
      nlim(1) = nRKinf(iwell,iprod)
      
      if (addTerm) then
         nlim(-1) = nlim(0) -1
         istart = -1
      else
         istart = 0
      endif

      RCname = PDname(inpwell,inpchan)
      
      do 300 index = istart,1           
         
         write(lkout,10) option,inpwell,iwell,iprod,
     2      RCname,PDname(iwell,iprod)
10       format(//,1x,a,'(',i2,') -> ',i2,':',i2,'  ',a,' = ',a)

         if (index.eq.-1) then 
            str = ' l'
            ipout = -1
            write(lkout,5) nlim(-1)
5          format(/,' Low pressure limit (add''l term): ',
     2         'multiply rate by M**(',i2,')')
     
         else if (index.eq.0) then 
            str = ' l'
            ipout = 0
            write(lkout,15) nlim(0)
15          format(/,' Low pressure limit: ',
     2         'multiply rate by M**(',i2,')')
     
         else 
            str = ' h'
            ipout = npres+1
            write(lkout,25) nlim(1)
25          format(/,' High pressure limit: ',
     2         'multiply rate by M**(',i2,')')
         endif

         do 100 it = 1,ntemps
            RKlog(it) = dlog(RK(iwell,iprod,it,ipout))
100      continue

         call polylin(lecho,ntemps,temp,RKlog,Nkfit,
     2      Coefs,errpar)
     
         do 110 icoef = 1,nKfit
            CoefLim(icoef,index) = Coefs(icoef)
110      continue

         write (lkout,115) fitID
115      format(1x,a,/)
            
         do 130 icoef = 1,Nkfit
            write(lkout,120) icoef-1,CoefLim(icoef,index)
cjmg fix
c120         format(1x,' Kcoef(',i2'):  ',1pe11.4) 
120         format(1x,' Kcoef(',i2,'):  ',1pe11.4) 
130      continue
    
         write(lkout,140) tab,tab,tab,tab,tab,tab
140      format(/,1x,'  T (K)',a,' 1000/T',a,'  P (atm) ',a,
     2     '     k    ',a,' log k ',a,' log kf',a,' % error')

         avgerr = 0.d0
         do 200 it = 1,ntemps

c          note we can regenerate solution with out knowing
c          what the model really is at this point
            RKfit = RKlim(temp(it),index)
     
            error = 100.d0* (RK(iwell,iprod,it,ipout) -
     2         RKfit) / dmin1(RK(iwell,iprod,it,ipout),RKfit)
     
            avgerr = avgerr + dabs(error)
     
            write(lkout,160) temp(it),tab,1.d3/temp(it),tab,str,
     2         tab,RK(iwell,iprod,it,ipout),tab,
     3         dlog10(RK(iwell,iprod,it,ipout)),tab,
     4         dlog10(RKfit),tab,error
160         format(1x,f7.1,a,f7.4,a,a,'.p.limit',a,1pe10.3,a,
     2         0pf7.3,a,f7.3,a,f8.2)

200      continue

         avgerr = avgerr/ntemps

         write(lkout,220) temp(1),temp(ntemps),avgerr
220      format(/,' fit from ',f6.1,' K to ',f6.1,' K with avg error ',
     2      'of ',f7.2,' %')
     
300   continue

      return
      end
      
c*POLYLIN*******************************
      subroutine PolyLin(lerr,nData,xin,yin,nCoefs,A,chiSq)
     
      implicit none
      include 'cdparams.fh'

c    polynomial fitting routine - we do this just like the new non-
c    arrhenius routine

c    local variables
      integer i
      integer lerr,nData,nCoefs
      real*8 xin(mxTpts),yin(mxTpts),A(nCoefs),chiSq
      integer listA(mxKfit)
      real*8 sig(mxTpts),covar(mxKfit,mxKfit)   
      external fPoly
c    end declarations

c    perform some initializing
      do 20 i = 1,nData
         sig(i) = 1.d0
20    continue

c    listA is a vector telling what coefficients to adjust;
c    we want all of them
      do 50 i = 1,nCoefs
         listA(i) = i
50    continue      

c    note we added add'l argument lerr to this sub
      call lfit (xin,yin,sig,nData,A,nCoefs,listA,nCoefs,covar,
     2   mxKfit,chiSq,fPoly,lerr)
        
      return
      end
      
