c***** tal 20000713
c created longtabout for full tree priniting

c****** tal 20000628
c modified to produce both text editor and Excel friendly output
c to load in Excel - select delimited format, with spaces as delimiters

c*******4/14/98  modified kfitout.f so that fit done over smaller T range (if needed)
c       to get better fit
c**********12/13/97  added leading tab to output
c********11/25/97  added tabs to chemkin output
c********12/6/95  modified chemkin output routine to suppress extra pressure printout

C***********************************************************************
c***********************************************************************
c  tabXMG  dmm 20020314
c
c     This routine prints out tables that are more easily read by xmg
c     than the painful long-string outputs.  When XMG logical is true,
c     should be called ...
c	
c
c***********************************************************************
c***********************************************************************

      subroutine tabXMG(option, lout)
      
      implicit none
      include 'cdparams.fh'
      include 'cdwell0.fh'
      include 'cdlabels.fh'
      include 'cdrange.fh'
      include 'cdrates.fh'   
      include 'cdcontrl.fh'

      character option*8
      integer lout
      integer it, ip, iwell, iprod, siowell
      real*8 sum

c***********************************************************************
c     Write output file header ...
c***********************************************************************

      write(lout,10)
 10   format('OUTPUT for PDEP NETWORK via CHEMDIS-XMG')

      if(option.eq.'Chemact') then
         write(lout, 15) PDname(inpWell, inpchan)
 15      format('CHEMACT case, input channel ', a25)
      else
         write(lout, 16) PDName(inpWell, 0)
 16      format('DISSOC case, input well ', a25)
      endif

c-----------------------------------------------------------------------
c***********************************************************************
c     BEGIN TEMPERATURE-PRESSURE  LOOP FOR RESULTS OUTPUT
c***********************************************************************
c-----------------------------------------------------------------------

      do 200 it = 1,ntemps
         do 201 ip = 1, npres
            
c***********************************************************************
c     BEGIN loop on wells and products
c***********************************************************************
            
            do 170 iwell = 1,nwells
               do 171 iprod = 0,nprods(iwell)
                  
c     For the dissoc case, don't print the results for the well
c     dissociating back to itself; otherwise, go ahead
                  if (.not.((option.eq.'Dissoc').and.
     2                 (iwell.eq.inpwell).and.(iprod.eq.0))) then
                     
                     if (ASA) then
c                        call MatchProductToSIsom(iwell, iprod,
c     $                       siowell)
                     else
                        siowell = -1
                     end if
                     
                     if ((siowell .lt. 0) .or. (iProd .eq. 0)) then
c     Complicated if-test removed here; unsure of effect if do this with
c     ASA on as a results.
                        
c     Finally, print the results ...
                        if(iprod.eq.0) then
                           write(lout, 20) temp(it), pres(it,ip),
     $                          PDname(iwell,iprod), RK(iwell,iprod,
     $                          it, ip)
                        else
                           write(lout, 21) temp(it), pres(it,ip),
     $                          PDname(iwell,iprod), RK(iwell,iprod,
     $                          it, ip)
                        endif					
                        
                     endif
c     Ends if((siowell ..)) test
                     
                  endif
c     Ends if(.not dissoc) test ...
                  
 171           continue
 170        continue
c     Closes loop over products and wells
 201     continue
 200  continue
c     Write an end line
      write(lout, 30)
 30   format('END')

 20   format(F6.1, 1x, E9.2, 1x, 'ISOMER  ', A25, 1x, E10.2)
 21   format(F6.1, 1x, E9.2, 1x,'PRODUCT ',A25,1x,E10.2)

      end

c***RMGOUT*************************
c print values of k(T,P) to be read by RMG -- pey 29/3/04
c similar to tabXMG, but better suited for printing values at many T and P's

      subroutine rmgout(option, lout)
      
      implicit none
      include 'cdparams.fh'
      include 'cdwell0.fh'
      include 'cdlabels.fh'
      include 'cdrange.fh'
      include 'cdrates.fh'   
      include 'cdcontrl.fh'
	include 'cdcheb.fh'

	logical isInNetwork
      character option*8, substr*8
      integer lout, nstrings
      integer it, ip, iwell, iprod
	real*8 rkleak(mxTPts, mxPPts), rkmax(mxTPts, mxPPts), rkCurr
	character addspc(mxTPts,mxPPts)*12

c	initialize rkleak, rkmax
	do it=1,mxTPts
		do ip=1,mxPPts
			rkleak(it,ip) = 0.0d0
			rkmax(it,ip) = 0.0d0
		enddo
	enddo

c	write file header
	if (cheby) then
		write(lout,'(2A)') 'OUTPUT for PDEP NETWORK via CHEMDIS-RMG',
     >                          ' (Chebyshev-polynomial format)'
	else
		write(lout,'(2A)') 'OUTPUT for PDEP NETWORK via CHEMDIS-RMG',
     >                          ' (Lookup-table format)'		
	endif

      if(option.eq.'Chemact') then
         write(lout, 15) PDname(inpWell, inpchan)
 15      format('CHEMACT case, input channel ', a25)
      else
         write(lout, 16) PDName(inpWell, 0)
 16      format('DISSOC case, input well ', a25)
      endif

c	temperatures and pressures
	if(cheby) then
		write(lout,'(A)') 'Temperature range'
		write(lout,22) tmin, tmax
		write(lout,'(A)') 'Pressure range'
		write(lout,22) pmin, pmax
	else
		write(lout,'(A)') 'Temperature'
		do it=1,ntemps
			write(lout,22) (temp(it), ip=1,npres)
		enddo
		write(lout,'(A)') 'Pressure'
		do it=1,ntemps
			write(lout,23) (pres(it,ip), ip=1,npres)
		enddo
	endif

	if (inclks) then ! print only k's included in the present pdep network

		do 180 iwell = 1,nwells ! loop over wells and products
			do 181 iprod = 0,nprods(iwell)
				
c	For the dissoc case, don't print the results for the well
c	dissociating back to itself; otherwise, go ahead
				if (.not.((option.eq.'Dissoc').and.
     2				 (iwell.eq.inpwell).and.(iprod.eq.0))) then

c					if a pathway is a PRODUCT channel and leads to only one species,
c					it's not in the pressure dependent network
					isInNetwork = .true.
					call parstr(PDname(iwell,iprod),2,substr,nstrings)
					if (nstrings.eq.1 .and. iprod.ne.0) then
						isInNetwork = .false.
					endif
					
					if(isInNetwork) then ! product is in network

					if(iprod.eq.0) then
						write(lout, 24) PDname(iwell,iprod)
					else
						write(lout, 25) PDname(iwell,iprod)
					endif

					if(cheby) then
						write(lout,'(A)')
     >                            'Cannot use CHEBYSHEV and INCLKS'
					else
						do it=1,ntemps
							write(lout,23) (RK(iwell,iprod,it,ip), 
     >                                        ip=1,npres)
						enddo
					endif ! ends cheby test

					else ! product is NOT in network

c					calculate rleak and record the non-included isomer with the max flux
					do it=1,mxTPts
						do ip=1,mxPPts
							rkCurr = RK(iwell,iprod,it,ip)
							rkleak(it,ip) = rkleak(it,ip) + rkCurr

							if (rkCurr.gt.rkmax(it,ip)) then
								rkmax(it,ip)  = rkCurr
								addspc(it,ip) = PDname(iwell,iprod)
							endif

						enddo
					enddo
					endif ! Ends if(isInNetwork) test

				endif ! Ends if(.not dissoc...) test
				
 181			continue
 180		continue

c		print rkleak, rkmax and addspc
		write(lout,'(A)') 'Kleak'
		do it=1,ntemps
			write(lout,23) (rkleak(it,ip), ip=1,npres)
		enddo
		write(lout,'(A)') 'NextIsomToAdd'
		do it=1,ntemps
			write(lout,26) (addspc(it,ip), ip=1,npres)
		enddo
		write(lout,'(A)') 'MaxNonInclRateConstant'
		do it=1,ntemps
			write(lout,23) (rkmax(it,ip), ip=1,npres)
		enddo
		
	else  ! print all k's

		do 190 iwell = 1,nwells ! loop over wells and products
			do 191 iprod = 0,nprods(iwell)
				
c	For the dissoc case, don't print the results for the well
c	dissociating back to itself; otherwise, go ahead
				if (.not.((option.eq.'Dissoc').and.
     2				 (iwell.eq.inpwell).and.(iprod.eq.0))) then
				   
					if(iprod.eq.0) then
						write(lout, 24) PDname(iwell,iprod)
					else
						write(lout, 25) PDname(iwell,iprod)
					endif

					if(cheby) then
						call chebyshev(iwell,iprod)
						do it=1,nchebT
						write(lout,27) (chebcoeffs(it,ip),ip=1,nchebP)
						enddo
					else
						do it=1,ntemps
							write(lout,23) (RK(iwell,iprod,it,ip), 
     >										ip=1,npres)
						enddo
					endif ! ends if(cheby) test

				endif
c	Ends if(.not dissoc) test ...
				
 191			continue
 190		continue		

 	endif ! end check if inclks = true or false

c     Write an end line
      write(lout, 30)
 30   format('END')

c	Formats
 22	format(500(F8.1))
 23	format(500(E11.3))
 24	format('ISOMER ', A25)
 25	format('PRODUCT', A25)
 26	format(1x,500(A12))
 27	format(10(E18.8))
      end
		
c***MATOUT*************************
c print matlab friendly table of k(T,P) -- pey 25/3/04
      subroutine matout(option, lout)
      
      implicit none
      include 'cdparams.fh'
      include 'cdwell0.fh'
      include 'cdlabels.fh'
      include 'cdrange.fh'
      include 'cdrates.fh'   
      include 'cdcontrl.fh'

      character option*8
      integer lout, nrxns, i
      integer it, ip, iwell, iprod, siowell
      real*8 sum

c	write file header
      write(lout,10)
 10   format('OUTPUT for PDEP NETWORK via CHEMDIS-XMG')

      if(option.eq.'Chemact') then
         write(lout, 15) PDname(inpWell, inpchan)
 15      format('CHEMACT case, input channel ', a25)
      else
         write(lout, 16) PDName(inpWell, 0)
 16      format('DISSOC case, input well ', a25)
      endif

	nrxns = 0
c	print list of products and isomers
	do iwell = 1,nwells
		do iprod = 0,nprods(iwell)

c     For the dissoc case, don't print the results for the well
c     dissociating back to itself; otherwise, go ahead
			if (.not.((option.eq.'Dissoc').and.
     2			(iwell.eq.inpwell).and.(iprod.eq.0))) then

				nrxns = nrxns+1
				if(iprod.eq.0) then
					write(lout,20) nrxns, PDname(iwell,iprod)
				else
					write(lout,21) nrxns, PDname(iwell,iprod)
				endif
			endif

		enddo
	enddo

 20	format('Column:',1x,I3,3x,'ISOMER ',A25)
 21   format('Column:',1x,I3,3x,'PRODUCT',A25)

c	column headings
	write(lout,24)
	write(lout,19) (i, i=1,nrxns)
 19	format(1x,'T(K)',3x,'P(atm)',4x,500(:,'Col.',I4,2x))
		
c	print table
      do 200 ip = 1, npres
      do 201 it = 1,ntemps

		write(lout, 22) temp(it), pres(it,ip)
		            
		do 170 iwell = 1,nwells
		do 171 iprod = 0,nprods(iwell)

c     For the dissoc case, don't print the results for the well
c     dissociating back to itself; otherwise, go ahead
			if (.not.((option.eq.'Dissoc').and.
     2		   (iwell.eq.inpwell).and.(iprod.eq.0))) then 
                        
				write(lout,23) RK(iwell,iprod,it,ip)				
                     
			endif
                  
 171		continue
 170		continue

	write(lout,24) ! carriage return
 201  continue
 200  continue

 22	format(F6.1, 1x, E9.2,\)
 23   format(E10.2,\)
 24	format(' ')

      end


c***KFITOUT************************
c    this sub does all the non-arrhenius fitting
      subroutine Kfitout(option,lecho,lkout,iwell,iprod,ipout)

      implicit none
      include 'cdparams.fh'
      include 'cdwell0.fh'
      include 'cdlabels.fh'
      include 'cdrange.fh'
      include 'cdrates.fh'
      include 'cdkfit.fh'
c*******added following line 4/14/98      
      include 'cdcolls.fh'

c    local variables
      character option*8,RCname*20
      integer it
      integer lkout,lecho,iwell,iprod,ipout
c*******added following line 4/14/98      
      integer icoll
      real*8 RKlog(mxTpts),RKfit,error,avgpres,errpar
c*******added following lines 4/14/98      
      real*8 Rkl(mxTpts), ttemp(mxTpts), tlow, thigh
      integer ntfit
c*******          4/14/98      
      character tab,str*2
c    end declarations

      tab = char(9)
      
      RCname = PDname(inpwell,inpchan)
      write(lkout,10) option,inpwell,iwell,iprod,
     2   RCname,PDname(iwell,iprod)
10    format(//,1x,a,'(',i2,') -> ',i2,':',i2,'  ',a,' = ',a)

c    low p limit
      if (ipout.eq.-1) then
         str = ' l'
         write(lkout,15) nRKinf(iwell,iprod) + nPr(iwell,iprod) - 1
15       format(1x,'Low pressure limit (add''l term): ',
     2      'multiply rate by M**(',i2,')')
    
      else if (ipout.eq.0) then 
         str = ' l'
         write(lkout,30) nRKinf(iwell,iprod) + nPr(iwell,iprod)
30       format(1x,'Low pressure limit: ',
     2      'multiply rate by M**(',i2,')')
 
c    all finite p    
      else if (ipout.le.npres) then
         avgpres = 0.d0
         do 50 it = 1,ntemps
            avgpres = avgpres + pres(it,ipout)/ntemps
50       continue      
         write(lkout,60) ipout,avgpres
60       format(1x,'Pressure sequence # ',i2,
     2      ';  P (nom): ',1pe10.3,' atm') 

c    high p limit
      else
         str = ' h'
         write(lkout,70) nRKinf(iwell,iprod)
70       format(1x,'High pressure limit: ',
     2      'multiply rate by M**(',i2,')')
      endif
         
      do 100 it = 1,ntemps
         RKlog(it) = dlog(RK(iwell,iprod,it,ipout))
100   continue

      call arrnlin(lecho,0.d0,ntemps,temp,RKlog,ntlow(ipout),
     2   nthigh(ipout),afit(ipout),rnfit(ipout),efit(ipout),
     3   errpar)
     
c    note Efit is now in kcal 
 
      write(lkout,140) tab,tab,tab,tab,tab,tab
140   format(/,1x,'  T (K)',a,' 1000/T',a,'  P (atm) ',a,
     2   '     k    ',a,' log k ',a,' log kf',a,' % error')


      avgerr(ipout) = 0.d0
      do 200 it = 1,ntemps
         
         RKfit = Afit(ipout)*temp(it)**rNfit(ipout)*
     2      dexp(-efit(ipout)/(1.987d-3*temp(it)))
     
         error = 100.d0* (RK(iwell,iprod,it,ipout) - RKfit) /
     2      dmin1(RK(iwell,iprod,it,ipout),RKfit)
     
         avgerr(ipout) = avgerr(ipout) + dabs(error)
     
         if ((ipout.le.0).or.(ipout.eq.npres+1)) then
            write(lkout,160) temp(it),tab,1.d3/temp(it),tab,str,tab,
     2         RK(iwell,iprod,it,ipout),tab,
     3         RKlog(it)/dlog(10.d0),tab,dlog10(RKfit),tab,error
160         format(1x,f7.1,a,f7.4,a,a,'.p.limit',a,1pe10.3,a,
     2         0pf7.3,a,f7.3,a,f8.2)
         else
            write(lkout,170) temp(it),tab,1.d3/temp(it),tab,
     2         pres(it,ipout),tab,RK(iwell,iprod,it,ipout),tab,
     3         RKlog(it)/dlog(10.d0),tab,dlog10(RKfit),tab,error
170         format(1x,f7.1,a,f7.4,a,1pe10.3,a,1pe10.3,a,0pf7.3,a,
     2         f7.3,a,f8.2)
         endif

200   continue

      avgerr(ipout) = avgerr(ipout)/ntemps
            write(lkout,210) afit(ipout),rnfit(ipout),efit(ipout)
210   format(/,' Fit:  A = ',1pe10.3,'  n = ',0pf6.2,'  E = ',f7.2,
     2   ' kcals')
      write(lkout,220) temp(ntlow(ipout)),temp(nthigh(ipout)),
     2   avgerr(ipout)
220   format(' over range ',f6.1,' K to ',f6.1,' K with avg error ',
     2   'of ',f7.2,' %')


c**********added 4/14/98********************************************************      
      if ((ipout.le.0).or.(ipout.gt.npres)) go to 208
      thigh=temp(ntemps)
c  checking to see if fit good enough (less than 15%)
      if  (avgerr(ipout).gt.15.0)  then
c      write (55, 222) avgerr(ipout)
c222   format('first check:  avgerr =', f7.2, '%')      
       avgerr(ipout) = 0.d0
c****** removing lowest T point       
       ntfit=ntemps-1
c      write(55,124) ntfit
c124   format('removing lowest T', '# T points=',i3)      
       do 101 it = 2,ntemps
        ttemp(it-1)=temp(it)
         RKl(it-1) = dlog(RK(iwell,iprod,it,ipout))
101   continue
      call arrnlin(lecho,0.d0,ntfit,ttemp,RKl,ntlow(ipout),
     2   nthigh(ipout),afit(ipout),rnfit(ipout),efit(ipout),
     3   errpar)
      do 201 it = 2,ntemps
         RKfit = Afit(ipout)*temp(it)**rNfit(ipout)*
     *      dexp(-efit(ipout)/(1.987d-3*temp(it)))
         error = 100.d0* (RK(iwell,iprod,it,ipout) - RKfit) /
     *      dmin1(RK(iwell,iprod,it,ipout),RKfit)
         avgerr(ipout) = avgerr(ipout) + dabs(error)
201   continue
      tlow=temp(2)    
      avgerr(ipout) = avgerr(ipout)/ntfit
c     write(55,601) tab,RCname,PDname(iwell,iprod),tab, afit(ipout),
c    2   tab, rnfit(ipout),tab,1.d3*efit(ipout),tab,
c    3   pres(ntlow(ipout),ipout), int(tlow+.5d0),     
c    4   int(thigh+.5d0),int(avgerr(ipout)+.5d0),
c    5   (xcoll(icoll),collNM(icoll),icoll=1,ncolls)
      else
c      write(55,134) 
c134   format('using all T')   
      tlow=temp(1)
      write(55,601) tab,RCname,PDname(iwell,iprod),tab, afit(ipout),
     2   tab, rnfit(ipout),tab,efit(ipout),tab,
     3   pres(ntlow(ipout),ipout), int(tlow+.5d0),     
     4   int(thigh+.5d0),int(avgerr(ipout)+.5d0),
     5   (xcoll(icoll),collNM(icoll),icoll=1,ncolls)
      go to 208
      endif 
      if((avgerr(ipout).gt.15.0).and.(ipout.gt.0).and.(ipout.le.npres))
     *  then
c      write (55, 223) avgerr(ipout)
c223   format('second check:  avgerr =', f7.2, '%')      
       avgerr(ipout) = 0.d0
c****** removing next lowest T point       
       ntfit=ntemps-2
c      write(55,154) ntfit
c154   format('removing 2nd lowest T', '# T points=',i3)      
       do 102 it = 3,ntemps
        ttemp(it-2)=temp(it)
         RKl(it-2) = dlog(RK(iwell,iprod,it,ipout))
102   continue
      call arrnlin(lecho,0.d0,ntfit,ttemp,RKl,ntlow(ipout),
     2   nthigh(ipout),afit(ipout),rnfit(ipout),efit(ipout),
     3   errpar)
      do 202 it = 3,ntemps
         RKfit = Afit(ipout)*temp(it)**rNfit(ipout)*
     *      dexp(-efit(ipout)/(1.987d-3*temp(it)))
         error = 100.d0* (RK(iwell,iprod,it,ipout) - RKfit) /
     *      dmin1(RK(iwell,iprod,it,ipout),RKfit)
         avgerr(ipout) = avgerr(ipout) + dabs(error)
202   continue
      tlow=temp(3)   
      avgerr(ipout) = avgerr(ipout)/ntfit
c     write(55,601) tab,RCname,PDname(iwell,iprod),tab, afit(ipout),
c    2   tab, rnfit(ipout),tab,1.d3*efit(ipout),tab,
c    3   pres(ntlow(ipout),ipout), int(tlow+.5d0),     
c    4   int(thigh+.5d0),int(avgerr(ipout)+.5d0),
c    5   (xcoll(icoll),collNM(icoll),icoll=1,ncolls)
      else
c      write(55,135) 
c135   format('fit OK if lowest T dropped')   
       tlow=temp(2)
      write(55,601) tab,RCname,PDname(iwell,iprod),tab, afit(ipout),
     2   tab, rnfit(ipout),tab,efit(ipout),tab,
     3   pres(ntlow(ipout),ipout), int(tlow+.5d0),     
     4   int(thigh+.5d0),int(avgerr(ipout)+.5d0),
     5   (xcoll(icoll),collNM(icoll),icoll=1,ncolls)
      go to 208
      endif
      if((avgerr(ipout).gt.15.0).and.(ipout.gt.0).and.(ipout.le.npres))
     *  then
c      write (55, 233) avgerr(ipout)
c233   format('third check:  avgerr =', f7.2, '%')      
       avgerr(ipout) = 0.d0
c****** removing next lowest T point       
       ntfit=ntemps-3
c      write(55,144) ntfit
c144   format('removing 3rd lowest T', '# T points=',i3)      
       do 103 it = 4,ntemps
        ttemp(it-3)=temp(it)
         RKl(it-3) = dlog(RK(iwell,iprod,it,ipout))
103   continue
      call arrnlin(lecho,0.d0,ntfit,ttemp,RKl,ntlow(ipout),
     2   nthigh(ipout),afit(ipout),rnfit(ipout),efit(ipout),
     3   errpar)
      do 203 it = 4,ntemps
          RKfit = Afit(ipout)*temp(it)**rNfit(ipout)*
     *      dexp(-efit(ipout)/(1.987d-3*temp(it)))
          error = 100.d0* (RK(iwell,iprod,it,ipout) - RKfit) /
     *      dmin1(RK(iwell,iprod,it,ipout),RKfit)
         avgerr(ipout) = avgerr(ipout) + dabs(error)
203   continue
      tlow=temp(4)
      avgerr(ipout) = avgerr(ipout)/ntfit
      write(55,601) tab,RCname,PDname(iwell,iprod),tab, afit(ipout),
     2   tab, rnfit(ipout),tab,efit(ipout),tab,
     3   pres(ntlow(ipout),ipout), int(tlow+.5d0),     
     4   int(thigh+.5d0),int(avgerr(ipout)+.5d0),
     5   (xcoll(icoll),collNM(icoll),icoll=1,ncolls)
      else
c       write(55,145) 
c145   format('fit OK if 2 lowest T dropped')   
      tlow=temp(3)
      write(55,601) tab,RCname,PDname(iwell,iprod),tab, afit(ipout),
     2   tab, rnfit(ipout),tab,efit(ipout),tab,
     3   pres(ntlow(ipout),ipout), int(tlow+.5d0),     
     4   int(thigh+.5d0),int(avgerr(ipout)+.5d0),
     5   (xcoll(icoll),collNM(icoll),icoll=1,ncolls)
601   format(1x,a1,a,' = ',a,a1,1pe9.2,a1,0pf6.2,a1,f7.2,a1,
     2   '!',1pe8.2,' atm, ',0pi4,'-',0pi4,' K,'
     3   0pi4,'% err, ',5(f4.2,' x ',a4,1x))
      endif
208   continue
c**********added 4/14/98********************************************************      
     
      return
      end
      
c*ARRNLIN*******************************
      subroutine arrnlin(lerr,dummy1,nData,xin,yin,nlow,nhigh,afit,
     2      rnfit,efit,chiSq)
     
      implicit none
      include 'cdparams.fh'
      
      integer nCoefs
      parameter (nCoefs = 3)      

c    non-linear arrhenius fitting routine - for now, we keep the same
c    input format as the original routine which we don't understand
c    and are completely replacing - also we use subs from num. recipies
c    (even though we don't understand these either, at least we know
c    where they are documented)    -   <ayc  4/93>

c    local variables
      integer i
      integer lerr,nData,nlow,nhigh
      real*8 xin(mxTpts),yin(mxTpts),afit,rnfit,efit,dummy1,chiSq
      integer listA(nCoefs)
      real*8 sig(mxTpts),A(nCoefs),covar(nCoefs,nCoefs)   
      external fArrn
c    end declarations

c    nlow and nhigh are part of old autoranging stuff - here, we fix and 
c    leave alone
      nlow = 1
      nhigh = nData 
      
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
     2   nCoefs,chiSq,fArrn,lerr)
        
c    set parameters for output
      afit = dexp(A(1))
      rnfit = A(2)
      efit = -A(3)*1.987d0
      
      return
      end
      
c***FARRN**************************
c    extern basis functions for lsq fit of arrhenius coefficients
      subroutine fArrn(x,VecCoefs,nCoefs)
      implicit none

      integer nCoefs      
      real*8 x,VecCoefs(nCoefs)
      
      VecCoefs(1) = 1.d0
      VecCoefs(2) = dlog(x)
      VecCoefs(3) = 1000.d0/x      
      
      return
      end
 
c***CKOUTPUT*******************************      
c**********4/14/98*************
c  commented out write statement--now in cdkouts
      subroutine CKoutput(option,lckout,iwell,iprod,ipout)

      implicit none
      include 'cdparams.fh'
      include 'cdwell0.fh'
      include 'cdlabels.fh'
      include 'cdrange.fh'
      include 'cdkfit.fh'
      include 'cdcolls.fh'
      
c    local variables
      character option*8,RCname*20
      character tab
      integer lckout,iwell,iprod,ipout,icoll

c    end declarations

      tab = char(9)




c    writes rates in chemkin compatible form - note chemkin form not
c    compatible with high and low pressure rate limits
       RCname = PDname(inpwell,inpchan)
c*******12/6/95  modified below:
c      write(lckout,60) RCname,PDname(iwell,iprod),afit(ipout),
c     2   rnfit(ipout),1.d3*efit(ipout),pres(ntlow(ipout),ipout),      
c     3   pres(nthigh(ipout),ipout),int(temp(ntlow(ipout))+.5d0),
c     4   int(temp(nthigh(ipout))+.5d0),int(avgerr(ipout)+.5d0),
c     5   (xcoll(icoll),collNM(icoll),icoll=1,ncolls)
c60    format(1x,a,' <=> ',a,1pe9.2,'  ',0pf6.2,'   ',f7.0,5x,
c     2   '!  ',1pe8.2,' (',1pe8.2,') atm, ',0pi4,'-',0pi4,' K,'
c     3   0pi4,'% err, ',5(f4.2,' x ',a4,1x))
c     write(lckout,60) RCname,PDname(iwell,iprod),afit(ipout),
c     2   rnfit(ipout),1.d3*efit(ipout),pres(ntlow(ipout),ipout),      
c     3   int(temp(ntlow(ipout))+.5d0),
c     4   int(temp(nthigh(ipout))+.5d0),int(avgerr(ipout)+.5d0),
c     5   (xcoll(icoll),collNM(icoll),icoll=1,ncolls)
c60    format(1x,a,' <=> ',a,1pe9.2,'  ',0pf6.2,'   ',f7.0,1x,
c     2   '!',1pe8.2,' atm, ',0pi4,'-',0pi4,' K,'
c     3   0pi4,'% err, ',5(f4.2,' x ',a4,1x))
c********11/25/97  modified to include tabs
c$$$       write(lckout,60) tab,RCname,PDname(iwell,iprod),tab, afit(ipout),
c$$$     2   tab, rnfit(ipout),tab,1.d3*efit(ipout),tab,
c$$$     3   pres(ntlow(ipout),ipout), int(temp(ntlow(ipout))+.5d0),     
c$$$     4   int(temp(nthigh(ipout))+.5d0),int(avgerr(ipout)+.5d0),
c$$$     5   (xcoll(icoll),collNM(icoll),icoll=1,ncolls)
c$$$ 60    format(1x,a1,a,' = ',a,a1,1pe9.2,a1,0pf6.2,a1,f7.0,a1,
c$$$     2   '!',1pe8.2,' atm, ',0pi4,'-',0pi4,' K,'
c$$$     3   0pi4,'% err, ',5(f4.2,' x ',a4,1x))

      return
      end
      


c***********************************************************************
c***********************************************************************
c
c  dmm 20010925
c
c     PrintCKLine
c
c     This function prints a line to fort.98 (for now), which can be
c     cut and pasted to a chemkin mechanism ... expects ratek to be in
c     normal, not log, units.  Different than the Albert Chang version
c     (CKOutput) in that fitting is not required; can be called after
c     ASA; less sophisticated in that p-dep fitted forms are not used.
c     This just puts out the rate constant as the A-factor.
c
c***********************************************************************     
c***********************************************************************     

      SUBROUTINE PrintCKinp(option, iWell, iProd, it, ip, ratek)
      
      implicit none
      include 'cdparams.fh'
      include 'cdlabels.fh'
      include 'cdwell0.fh'
      
      character option*8
      integer iWell, iProd, it, ip
      real*8 ratek

      integer lprod

      if(option.eq.'Chemact') then
         write(98,100) reName, pdName(iWell,iProd), ratek
         if(iProd.ne.0) then
c     search for whether a duplicate channel exists ... it will need to
c     be labeled if it does ...
            do 10 lprod = 1, nprods(iWell)
               if((pdName(iWell,lprod).eq.pdName(iWell,iProd)).and
     $              .(lprod.ne.iProd)) then
                  write(98,101) 
               endif
 10         continue
         endif
      else if(option.eq.'Dissoc') then
         
      endif
      return
 100  format(a,'=',a,T31,E7.2)
 101  format('DUPLICATE')
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
c Revision 1.11  2002/03/08 22:45:27  dmmatheu
c before changes to produce chemdis-xmg version with easier-to-read
c output for reading into xmg ...
c
c Revision 1.10  2001/09/26 22:38:53  dmmatheu
c before changing units of CK printout to kcal/mol ...
c
c Revision 1.9  2001/09/25 22:52:52  dmmatheu
c before changes to add PrintCKLine function
c
c Revision 1.8  2000/08/29 00:36:55  dmmatheu
c after tlada merge ... no changes
c
c Revision 1.7  2000/08/18 13:17:54  dmmatheu
c Major changes to output by tlada ... I don't know them all!
c
c Revision 1.6  2000/08/09 21:21:43  dmmatheu
c after tlada changes, addition of GetRK function ...
c
c Revision 1.5  2000/07/17 13:20:23  dmmatheu
c after tlada changes to output
c
c Revision 1.4  2000/07/17 13:18:02  dmmatheu
c before tlada merge ...
c
c Revision 1.3  2000/07/07 19:52:20  dmmatheu
c product 'isomer' channel flagging
c
c Revision 1.2  2000/06/29 13:29:59  dmmatheu
c tlada version
c
