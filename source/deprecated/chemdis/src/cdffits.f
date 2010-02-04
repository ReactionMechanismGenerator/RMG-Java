c***FITF************************************      
      subroutine FitF(lecho)
      
      implicit none
      include 'cdparams.fh'
      include 'cdrange.fh'
      include 'cdffunc.fh'
      include 'cdffit.fh'

c    local variables
      integer iter
      integer Iparam(7),lecho,mxIter,mxFEval,mxGEval
      real*8 Rparam(7),Xguess(mxFitP),Xscale(mxFitP),Fvalue,
     2   X(mxFitP),WK(mxFitP*(mxFitP+8)),Fscale
      external FCN
c    end declarations           

c    get initial guess for solution; scaling info              
      call initX(Nfit,Xguess,Xscale,Fscale,funstr,Fopt)

c    intialize minimizer 
c    flag to reset parameters 
      Iparam(1) = 1  
c      call DU4INF(Iparam,Rparam)

c    change some parameters for minimizer      
      call resetP(Iparam,Rparam,mxIter,mxFEval,mxGEval)

c    U2INF is the IMSL minimization fitting routine
c    NOTE U2INF is really for single precision but we
c    always compile with opt CRAY real*8 ->> single precision
c    (changed to du2inf for SGI  <ayc 3/94>)

c      call DU2INF(FCN,Nfit,Xguess,Xscale,Fscale,Iparam,Rparam,
c     2   X,Fvalue,WK) 
     
c    test for convergence
      if (Iparam(3).ge.mxIter) write(lecho,100) mxIter
100   format(/,' ERROR: No convergence in ',i5,' iterations.')

      if (Iparam(4).ge.mxFEval) write(lecho,110) mxFEval
110   format(/,' ERROR: No convergence in ',i5,
     2   ' function evaluations.')
     
      if (Iparam(5).ge.mxGEval) write(lecho,120) mxGEval
120   format(/,' ERROR: No convergence in ',i5,
     2   ' gradient evaluations.')

c    fill out fit parameters into common variables
      do 150 iter = 1,Nfit
         FitVec(iter) = X(iter)
150   continue

      return
      end

c***RESETP************************************           
      subroutine resetP(iparam,Rparam,mxI,mxF,mxG)
c    sets up parameters for minimizer - see IMSL doc.

      parameter (mxIter = 1000, mxFEval = 4000, mxGEval = 4000)
      
      integer Iparam(7),mxI,mxF,mxG
      real*8 Rparam(7)

c    all these things have defaults, but some
c    we reset

c    iparam(1) was control flag which we already invoked

c    iparam(2) is # of machine-dependent good digits 

c    #3 is max # of iterations (def: 100)         
      Iparam(3) = mxIter
      mxI = mxIter
      
c    #4 is max # of function evaluations (def: 400)
      Iparam(4) = mxFEval   
      mxF = mxFEval
      
c    #5 is max # of gradient evaluations (def: 400)
      Iparam(5) = mxGEval   
      mxG = mxGEval
      
c    iparam(6) def = 0 - init Hessian to Identity

c    iparam(7) is not used

c    Rparam defaults look ok - based on machine
c    dependent properties - i don't think we can
c    do better - so we leave alone this time

      return
      end

c***FCN************************************      
c    subroutine to define objective function
c    see IMSL documentation

      subroutine FCN(N, X, F)
      
      implicit none
      include 'cdparams.fh'
      include 'cdrange.fh'
      include 'cdffunc.fh'
      include 'cdffit.fh'

c    local variables      
      integer N
      integer it, ip
      real*8 X(N),error,cumerr,Ffit,F
      real*8 getFfit
c    end declarations      

      cumerr = 0.d0
            
      do 100 it = 1,ntemps
         do 100 ip = 1,npres

            Ffit = getFfit(N,X,temp(it),Pr(it,ip),Fopt)            
            error = (dlog10(Fcalc(it,ip)) - dlog10(Ffit))**2
                    
            cumerr = cumerr + error 
100   continue

      F = cumerr
      
      return
      end

c***FITFOUT************************************      
      subroutine FitFout(option,lfout1,lfout2,iwell,iprod)
c    fill out output files with Fcalc vs T, P
      
      implicit none
      include 'cdparams.fh'
      include 'cdwell0.fh'
      include 'cdlabels.fh'
      include 'cdrange.fh'
      include 'cdffunc.fh'
      include 'cdffit.fh'
                         
c    local variables
      character option*8,RCname*20
      integer it,ip,iwell,iprod,iter
      integer lfout1,lfout2,itEmax,ipEmax,nFlt90,nEgt10
      real*8 Ffit,error,errormx,toterr
      real*8 getFfit
      character tab
c    end declarations

      tab = char(9)
      
      RCname = PDname(inpwell,inpchan)
      write(lfout1,50) option,inpwell,iwell,iprod,
     2   RCname,PDname(iwell,iprod)
50    format(//,1x,a,'(',i2,') -> ',i2,':',i2,'  ',a,' = ',a)

      write(lfout1,150) tab,tab,tab,tab,tab,tab,tab,tab
150   format(/,1x,' T (K) ',a,' 1000/T',a,' log P ',a,' log Pr',a,
     2   ' Fcalc ',a,'-lnFcal',a,' F fit ',a,'-lnFfit',a,' %error')

      nFlt90 = 0
      nEgt10 = 0
      errormx = 0.d0
      toterr = 0.d0
           
      do 350 it = 1,ntemps     
         do 300 ip = 1,npres
            Ffit = getFfit(Nfit,FitVec,temp(it),Pr(it,ip),Fopt)
            if (Ffit.lt.9.d-1) nFlt90 = nFlt90 + 1
            
            error = (Fcalc(it,ip)-Ffit) / dmin1(Fcalc(it,ip),Ffit)
            if (dabs(error).gt.1.d-1) nEgt10 = nEgt10 + 1
            
            toterr = toterr + dabs(error)
            
            if (dabs(error).gt.errormx) then
               errormx = dabs(error)
               itEmax = it
               ipEmax = ip
            endif
            
            write(lfout1,280) temp(it),tab,1000.d0/temp(it),
     2         tab,dlog10(pres(it,ip)),tab,dlog10(Pr(it,ip)),
     3         tab,Fcalc(it,ip),tab,-dlog(Fcalc(it,ip)),
     4         tab,Ffit,tab,-dlog(Ffit),tab,error*100.d0
280         format(' ',f7.2,a,f7.4,a,f7.3,a,f7.3,a,f7.4,a,f7.3,
     2         a,f7.4,a,f7.3,a,f7.2)
300      continue
350   continue

c    now write summary output to lfout2
      write(lfout2,50) option,inpwell,iwell,iprod,
     2   RCname,PDname(iwell,iprod)
      
      write(lfout2,400) funstr
400   format(/,' Fit: ',a,/)

      do 430 iter = 1,Nfit
         write(lfout2,410) iter,FitVec(iter)
410      format(' Param(',i2,'): ',1pe14.7)
430   continue

      write(lfout2,450) npres*ntemps,nFlt90,nEgt10,
     2   100.d0*toterr/(npres*ntemps)
450   format(/,' Total points: ',i3,'  # F''s < 0.9: ',i3,
     2   '  # err > 10%: ',i3,'  avg err%: ',f5.2)

      Ffit = getFfit(Nfit,FitVec,temp(itEmax),
     2   Pr(itEmax,ipEmax),Fopt)
     
      write(lfout2,460) errormx*100.d0,temp(itEmax),
     2   dlog10(pres(itEmax,ipEmax)),Fcalc(itEmax,ipEmax),Ffit
460   format(/,' Worst error%: ',f5.1'  T: ',f6.1,'  logP: ',f5.2,
     2   '  Fcalc: ',f7.4,'  Ffit: ',f7.4)
    
      return
      end
      
      
