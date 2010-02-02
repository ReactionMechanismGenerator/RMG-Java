c***INITX************************************      
c    sub to return initial-guess solution and scaling parameters
      subroutine initX(Npar,X0,Xscale,Fscale,IDstr,option)
      implicit none
      
c    local variables
      integer Npar
      real*8 X0(*),Xscale(*),Fscale
      character IDstr*(*),option*8
c    end declarations

c    IMPORTANT: these options are set in getparams

      if (option.eq.'Troe') then
         call initTroe(Npar,X0,Xscale,Fscale,IDstr)
      
      else if (option.eq.'SRI') then
         call initSRI(Npar,X0,Xscale,Fscale,IDstr)
      
      else
         call initDEF(Npar,X0,Xscale,Fscale,IDstr)
              
      endif
      
      return
      end

c***GETFFIT************************************      
c    routine to return fitted F
      function getFfit(Npar,X,T,P,option)
      implicit none
      
c    local variables
      integer Npar
      real*8 X(*),T,P
      real*8 getFfit,getTroe,getSRI,getDef
      character option*8
c    end declarations

      if (option.eq.'Troe') then
         getFfit = getTroe(Npar,X,T,P)
      
      else if (option.eq.'SRI') then
         getFfit = getSRI(Npar,X,T,P)
      
      else
         getFfit = getDef(Npar,X,T,P)
         
      endif
      
      return
      end
      
c***INITDEF************************************      
c    companion routine to getDef for initializing variables
      subroutine initDef(Npar,X0,Xscale,Fscale,IDstr)
      implicit none
      
c    local variables
      integer iter
      integer Npar,iamp,ioff,iwidth,iasym
      real*8 X0(*),Xscale(*),Fscale
      character IDstr*(*)
c    end declarations

c    make consistent with getDef 
      data iamp,iwidth,iasym,ioff/1,6,9,12/
      Npar = 12
      
c    no better idea so set Xscale to 1 - see IMSL description;
c    also zero out all X0's before setting non-zero values
c    ok - let's try not quite zeroing out (we want init vector
c    all the same order of magnitude)

c    enter FIT function ID string here
      IDstr = 'G (12): 5 amp, 3 width, 3 asym, 1 shift, exp = 1/2'
            
      do 10 iter = 1,Npar
         Xscale(iter) = 1.d0
         X0(iter) = 1.d-1
10    continue
      Fscale = 1.d0      
      
c    Now we reset the amplitude and width
      X0(iamp) = 1.d0
      X0(iwidth) = 2.d0

      return
      end
      
c***GETDEF************************************      
c    function to return defaulf F fitting function
      function getDef(Npar,X,T,P)
      implicit none
      
c    local variables
      integer Npar,iamp,ioff,iwidth,iasym
      real*8 X(*),T,P,FcentLn,FLShape,Plgoff,Width,asym,Plg,
     2   arg,Tnm
      real*8 getDef
c    end declarations
     
      data iamp,iwidth,iasym,ioff/1,6,9,12/

c    IMPORTANT: data statement must be consistent with initDef
c               IDstr set in initDef

c    we invent new scaled T
      Tnm = (1.d-3*T)**0.5
      
c    use base 10 log for P
      Plg = dlog10(P)
      
c    a: we use first 5 parameters for amplitude function
c    and we try fitting to cubic in 1000./T      
      FcentLn = X(iamp) + X(iamp+1)/(Tnm) + X(iamp+2)/(Tnm*Tnm)
     2   + X(iamp+3)/(Tnm*Tnm*Tnm) + X(iamp+4)/(Tnm*Tnm*Tnm*Tnm)

c    the rest are lineshape parameters:
c    b: we use 3 parameters for Half-width
      Width = X(iwidth) + X(iwidth+1)/(Tnm) + X(iwidth+2)/(Tnm*Tnm)
     
c    c: we use 3 parameters for asymmetry parameter     
      asym = X(iasym) + X(iasym+1)/(Tnm) + X(iasym+2)/(Tnm*Tnm)

c    d: we use 1 parameter for offset from 0
      Plgoff =  X(ioff)
      
c    now the lineshape function

c    Lorentzian: we don't use it but we keep it just in case
c     FLShape = 1.d0 /
c    2   ( 1.d0 + ((Plg + Plgoff)/(Width + asym*Plg))**2 )

c    Gaussian:
      FLShape = dexp( -((Plg + Plgoff)/(Width + asym*Plg))**2 )

c    finally return the solution; but constrain argument
c    to exponential to avoid overflows
         
      arg = -FLShape*FcentLn
      if (arg.gt.100.d0) arg = 100.d0
      if (arg.lt.-100.d0) arg = -100.d0

      getDef = dexp(arg)
      
      return
      end
      
c***INITTROE************************************      
c    companion routine to getTroe for initializing variables
      subroutine initTroe(Npar,X0,Xscale,Fscale,IDstr)
      implicit none
      
c    local variables
      integer iter
      integer Npar
      real*8 X0(*),Xscale(*),Fscale
      character IDstr*(*)
c    end declarations

      Npar = 4
     
c    no better idea so set Xscale to 1 - see IMSL description;
c    set all x's to order 1 (we define T's as x's * 1000)

c    enter FIT function ID string here
      IDstr = '4-Parameter Troe Fit in T/1000 (w/atan)'
            
      do 10 iter = 1,Npar
         Xscale(iter) = 1.d0
         X0(iter) = 5.d-1
10    continue
      Fscale = 1.d0      
      
c    Now we reset the a factor
      X0(1) = 4.d-1

      return
      end
      
c***GETTROE************************************      
c    function to return defaulf F fitting function
      function getTroe(Npar,X,T,P)
      implicit none
      
c    local variables
      integer Npar,icoef
      real*8 X(*),T,P,TStar(3),a,Fcent,
     2   FLShape,Rn,Rc,Rd,Plg
      real*8 getTroe
c    end declarations
     
      data Rd /0.14/

c    we try atan function to constrain a between 0 and 1
      a = datan(X(1))/3.14159265359d0 + .5d0

c    we wish to constrain T's to be positive
      do 10 icoef = 1,3
         Tstar(icoef) = 1.d3*dabs(X(icoef+1))
10    continue

c    use base 10 log for P
      Plg = dlog10(P)

c    this is right out of the Chemkin II manual, pg 23      
      Fcent = (1.d0 - a)*dexp(-T/TStar(3)) + a*dexp(-T/TStar(1))
     2   + dexp(-TStar(2)/T)
     
      Rc = -0.4d0 - 0.67d0*dlog10(Fcent)
      Rn = 0.75 - 1.27*dlog10(Fcent)
      
      FLShape = 1.d0 /
     2   ( 1.d0 + ( (Plg + Rc)/ (Rn - Rd*(Plg + Rc)) )**2 )

      getTroe = Fcent**(FLShape)
      
      return
      end
      
c***INITSRI************************************      
c    companion routine to getSRI for initializing variables
      subroutine initSRI(Npar,X0,Xscale,Fscale,IDstr)
      implicit none
      
c    local variables
      integer iter
      integer Npar
      real*8 X0(*),Xscale(*),Fscale
      character IDstr*(*)
c    end declarations

      Npar = 3
     
c    no better idea so set Xscale to 1 - see IMSL description;

c    enter FIT function ID string here
      IDstr = '3-Parameter SRI Fit in T/1000'
            
      do 10 iter = 1,Npar
         Xscale(iter) = 1.d0
         X0(iter) = 1.d0
10    continue
      Fscale = 1.d0      
      
c    Now we reset the a factor
      X0(1) = 4.d-1

      return
      end
      
c***GETSRI************************************      
c    function to return defaulf F fitting function
      function getSRI(Npar,X,T,P)
      implicit none
      
c    local variables
      integer Npar
      real*8 X(Npar),T,P,a,b,c,Fcent,FLShape,Plg
      real*8 getSRI
c    end declarations

c    (tried an equiv statement here, but compiler didn't
c    like it - probably because of passed variable)

      a = dabs(X(1))
      b = 1.d3*dabs(X(2))
      c = 1.d3*dabs(X(3))
     
c    use base 10 log for P
      Plg = dlog10(P)

c    this is right out of the Chemkin II manual, pg 24 -
c    but we have not included the optional d and e factors,
c    because we don't see that kind of limiting behavior
     
      Fcent = a*dexp(-b/T) + dexp(-T/c)
     
      FLShape = 1.d0 / ( 1.d0 + Plg**2 )

      getSRI = Fcent**(FLShape)
      
      return
      end
      
c  we tried these subs but we don't use them anymore
c**BOUNCR************************************ 
c    func to reflect fitting parameter to right of of boundary
      function bouncR(xval,xlim,xzone)
      implicit none

      real*8 xval,xlim,xzone
      real*8 bouncR

c    bouncR is the reflection of xval to the other side
c    of xlim but within the region defined by xzone

c    xzone should be small compared to expected magnitude of
c    xval (because legitimate xval's could be stuck in the 
c    zone); still xzone should be large enough that things
c    put there don't simple vanish into the noise

c    the test condition (xval < xlim) is assumed

      bouncR = xlim + xzone/(1.d0 + (xlim - xval))
      
      return
      end

c**BOUNCL************************************ 
c    func to reflect fitting parameter to left of of boundary
      function bouncL(xval,xlim,xzone)
      implicit none
      
      real*8 xval,xlim,xzone
      real*8 bouncL

c    bouncL is the reflection of xval to the other side
c    of xlim but within the region defined by xzone

c    xzone should be small compared to expected magnitude of
c    xval (because legitimate xval's could be stuck in the 
c    zone); still xzone should be large enough that things
c    put there don't simple vanish into the noise

c    the test condition (xval > xlim) is assumed

      bouncL = xlim - xzone/(1.d0 + (xval - xlim))
      
      return
      end

