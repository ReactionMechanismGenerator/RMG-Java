c***CALCF**********************
      subroutine calcF(option,iwell,iprod)
      
      implicit none
      include 'cdparams.fh'
      include 'cdrange.fh'
      include 'cdrates.fh'
      include 'cdlimfit.fh'
      include 'cdffunc.fh'

c    local variables
      character*8 option
      integer it,ip,iwell,iprod
      real*8 RKlind,concM
c    end declarations 
      
      do 500 it = 1,ntemps
         do 500 ip = 1,npres
         
            concM = pres(it,ip) / (82.1d0*temp(it))

c          get lindemann forms (depends on depth of well and
c          whether stabilization or product channel);  also do
c          special dissoc case

            if (addTerm) then
               Pr(it,ip) = (RK(iwell,iprod,it,0) + 
     2            RK(iwell,iprod,it,-1)/concM)
     3            *concM**nPr(iwell,iprod) / 
     4            RK(iwell,iprod,it,npres+1)        
c          or, do general case
            else              
               Pr(it,ip) = RK(iwell,iprod,it,0)
     2            *concM**nPr(iwell,iprod) / 
     3            RK(iwell,iprod,it,npres+1)
            endif
                 
            RKlind = RK(iwell,iprod,it,npres+1)
     2         *concM**nRKinf(iwell,iprod) 
     3         *Pr(it,ip) / (1.d0 + Pr(it,ip))
                     
            Fcalc(it,ip) = RK(iwell,iprod,it,ip) / RKlind
500   continue
       
      return
      end
       
c***FOUTPUT************************
      subroutine Foutput(option,lfout,iwell,iprod)
c    fill out output files with Fcalc vs T, P

      implicit none
      include 'cdparams.fh'
      include 'cdwell0.fh'
      include 'cdlabels.fh'
      include 'cdrange.fh'
      include 'cdffunc.fh'
                         
c    local variables
      character option*8,RCname*20
      integer it,ip,iwell,iprod
      integer lfout,ipmax
      real*8 Flognm
      character tab
c    end declarations

      tab = char(9)
      
      RCname = PDname(inpwell,inpchan)
      write(lfout,50) option,inpwell,iwell,iprod,
     2   RCname,PDname(iwell,iprod)
50    format(//,1x,a,'(',i2,') -> ',i2,':',i2,'  ',a,' = ',a)

      write(lfout,150) tab,tab,tab,tab,tab
150   format(/,1x,' T (K) ',a,' log P ',a,' log Pr',a,' Fcalc ',
     2   a,' -ln F ',a,'ln F |n')
     
      do 350 it = 1,ntemps     
         write(lfout,270) tab,tab,tab,temp(it),tab,temp(it),tab,temp(it)
270      format(' ',a,a,3(a,f5.0,' K'))

c        get lineshape maxima
         ipmax = 1
         do 200 ip = 1,npres
            if (Fcalc(it,ip).le.Fcalc(it,ipmax)) ipmax = ip
200      continue

         Flognm = -dlog(Fcalc(it,ipmax))

c       make sure there are no divide overflows
         if (Flognm.gt.0) then

            do 300 ip = 1,npres
               write(lfout,280) temp(it),tab,dlog10(pres(it,ip)),tab,
     2            dlog10(Pr(it,ip)),tab,Fcalc(it,ip),tab,
     3            -dlog(Fcalc(it,ip)),tab,-dlog(Fcalc(it,ip))/Flognm
280            format(' ',f7.2,a,f7.3,a,f7.3,a,f7.4,a,f7.4,a,f7.4)
300         continue

         else
            do 320 ip = 1,npres
               write(lfout,280) temp(it),tab,dlog10(pres(it,ip)),tab,
     2            dlog10(Pr(it,ip)),tab,Fcalc(it,ip),tab,
     3            -dlog(Fcalc(it,ip))
320         continue
         endif

350   continue

      return
      end
      
c***TOUTPUT************************
      subroutine Toutput(option,ltout,iwell,iprod)
c    fill out output files with Fcalc vs T, P

      implicit none
      include 'cdparams.fh'
      include 'cdwell0.fh'
      include 'cdlabels.fh'
      include 'cdrange.fh'
      include 'cdffunc.fh'

c    local variables
      character option*8,RCname*20
      integer it,ip,iwell,iprod
      integer ltout,ipmax,ipleft,ipright
      real*8 slope,PRLlog,PRRlog,PRsum,PRdif
      character tab
c    end declarations

      tab = char(9)
      
      RCname = PDname(inpwell,inpchan)
      write(ltout,50) option,inpwell,iwell,iprod,
     2   RCname,PDname(iwell,iprod)
50    format(//,1x,a,'(',i2,') -> ',i2,':',i2,'  ',a,' = ',a)

      write(ltout,160) tab,tab,tab,tab,tab,tab
160   format(/,1x,' T (K) ',a,' log P ',a,' log Pr',a,' Fcalc ',
     2   a,' -ln F ',a,'log Pr+',a,'log Pr-')
     
      do 400 it = 1,ntemps
          
c        get lineshape maxima (F minima)
         ipmax = 1
         do 200 ip = 1,npres
            if (Fcalc(it,ip).le.Fcalc(it,ipmax)) ipmax = ip
200      continue

c       proceed if Fmin is not too close to unity (otherwise problems)
         if (-dlog10(Fcalc(it,ipmax)).gt.1.d-4) then

c          get lineshape left and right HWHM
            ipleft = 1
            ipright = npres

            do 230 ip = 1,npres
               if ((dlog(Fcalc(it,ip))/dlog(Fcalc(it,ipmax)).lt.0.5d0)
     2            .and.(ip.lt.ipmax)) ipleft = ip
               if ((dlog(Fcalc(it,ip))/dlog(Fcalc(it,ipmax)).gt.0.5d0)
     2            .and.(ip.gt.ipmax)) ipright = ip
230         continue

c          we do some interpolation (could be trouble if F is constant)
            slope = (dlog10(PR(it,ipleft+1))-dlog10(PR(it,ipleft)))/
     2         (dlog(Fcalc(it,ipleft+1))-dlog(Fcalc(it,ipleft)))
            PRLlog = dlog10(PR(it,ipleft))+slope*
     2         (0.5d0*dlog(Fcalc(it,ipmax))-dlog(Fcalc(it,ipleft)))
            if (ipright.eq.npres) ipright = npres-1
            slope = (dlog10(PR(it,ipright+1))-dlog10(PR(it,ipright)))/
     2         (dlog(Fcalc(it,ipright+1))-dlog(Fcalc(it,ipright)))
            PRRlog = dlog10(PR(it,ipright))+slope*
     2         (0.5d0*dlog(Fcalc(it,ipmax))-dlog(Fcalc(it,ipright)))

c          note PRLlog is negative     
            PRsum = PRRlog - PRLlog 
            PRdif = PRRlog + PRLlog 

            write(ltout,250) temp(it),tab,dlog10(pres(it,ipmax)),tab,
     2         dlog10(PR(it,ipmax)),tab,Fcalc(it,ipmax),tab,
     3         -dlog(Fcalc(it,ipmax)),tab,PRsum,tab,PRdif
250         format(' ',f7.2,a,f7.3,a,f7.3,a,f7.4,a,f7.4,a,f7.3,a,f7.3)

         else
            write(ltout,250) temp(it),tab,dlog10(pres(it,ipmax)),tab,
     2         dlog10(PR(it,ipmax)),tab,Fcalc(it,ipmax),tab,
     3         -dlog(Fcalc(it,ipmax))
         endif

400   continue

      return
      end

