c********12/6/95  having trouble getting code to recognize N2, CH4 etc as colliders,
c                 but works with AR and HE--suspect problem might be the routine
c                 'upcase' in "COLL' section--will comment out and observe effect
c***GETPARAMS*************************
      subroutine getparams(linput,lecho,dE)

      implicit none
      include 'cdparams.fh'
      include 'cdwell0.fh'
      include 'cdwell1.fh'
      include 'cdwell2.fh'
      include 'cdisprop.fh'
      include 'cdcolls.fh'
      include 'cdlabels.fh'
      include 'cdrange.fh'
      include 'cdffit.fh'
      include 'cdcontrl.fh'
	include 'cdcheb.fh' ! added by pey 4/6/04

c dmm 20000224
      include 'cdgamma.fh'
c dmm 20000224
c dmm 20000330 
      include 'straightBS.fh'
      include 'cdnos.fh'
c dmm 20000330

c    local variables
      real*8 degendef(mxfreqs),freqdef(mxfreqs),sigdef,ekdef,xtot,dE
      integer iWell,jWell,ioWell,ifrq,it,ip,icoll, iTS, isp
      integer linput,lecho,numRC,nfdef
	real*8 deltat, deltap
      character instr*50,cmdstr*20,Wname*15, fname*15
      logical readin,readma,readdf,readdc,tpcorr,fillp,noinfo
      logical readfr(mxWells),readcp(mxWells), trange, prange

	real*8 pi
	parameter (pi = 3.14159265)
c    end declarations

c    these variables are for checking input
      numRC = 0
      readin = .false.
      readma = .false.
      readdf = .false.
      readdc = .false.
      fillp = .false.
   
      
c    controls reassignment      
      do 10 iWell = 1,mxWells     
         readcp(iWell) = .false.
         readfr(iWell) = .false.
10    continue 
      
c    defaults
      fitrange = .false.
      fitglobal = .false.
      showline = .false.
      tpcorr = .false.
      chemact = .false.
      dissoc = .false.
      inpwonly = .false.
      isomers = .false.
      rkold = .false.
      table = .true.
      caseStr = ' '
      Fopt = 'default'
      dE = 1.d0
      nRot = 1
c dmm 20000331
      gammaspec = .false.
      bsgs = 10.0
c dmm 20000331      
c tal 20000620-20000822
      oldrho = .false.
      ASA = .false.
      ASActrl = .false.
      noASAstepout = .false.
      nosumisom = .false.
      absrmin = .false.
c tal 20000822
c dmm 20010612 Barker file output
      barker = .false.
c dmm 20010612
      steady = .false.
c dmm 20010813
      dEdcm = .false.
c dmm 20011019
      RRKM = .false.
c dmm 20020314
      XMG = .false.
c pey 29/3/04
	inclks = .false.
	cheby  = .false.
	trange = .false.
	prange = .false.

c    we only zero out variables as are completely necessary    
      ncolls = 0
      nWells = 0
      xtot = 0.d0
      do 50 iWell = 1,mxWells
         nfreqs(iWell) = 0
         nProds(iWell) = 0
         idpWell(iWell) = 0
         do 45 ioWell = 1,mxWells
            Aisom(iWell,ioWell) = 0.0d0
            Eisom(iWell,ioWell) = 0.0d0
            rNisom(iWell,ioWell) = 0.0d0
            alphaIsom(iWell,ioWell) = 0.0d0
45       continue
50    continue

c1   loop for reading 1st level commands

c    main entry point
100   continue
         read(linput,'(a)',end=450) instr
         call ljust(instr)
         cmdstr = instr
         call upcase(cmdstr)
          
c    entry point for inner loops terminated on unknown command  
110   continue

c       comment or blank: do nothing
         if ((cmdstr(1:1).eq.'#').or.(cmdstr(1:1).eq.' ')) then
            goto 100
        
c       get casestr (tag)
         else if (cmdstr(1:3).eq.'TAG') then
         read(linput,'(a)') instr
         call ljust(instr)
         caseStr = instr
        
c     tal 20000620
c       use old density of states instead of fast
         else if (cmdstr(1:3).eq.'OLD') then
            oldrho = .true.
c     tal 20000620

c     tal 20000627
c       activates ASA algorithm
         else if (cmdstr(1:3).eq.'ASA') then
            ASA = .true.
            ASActrl = .true.
c     tal 20000705
c       allows specification of rmin by using ASA RMIN keywords
            call rjust(cmdstr)
            rminmult = -1.d0
c               write(*,*) cmdstr(len(cmdstr)-6:len(cmdstr))
            if (cmdstr(len(cmdstr)-3:len(cmdstr)) .eq. 'RMIN') then
               read(linput,*,err=999) rminmult
c     tal 20000822
c        allows specification of absolute RMIN
c        (units differ between chemact and dissoc!!!)
            else if (cmdstr(len(cmdstr)-6:len(cmdstr)).eq.'RMINABS')then
               absrmin = .true.
               read(linput,*,err=999) rminmult
c     tal 20000822
            else
               rminmult = 0.01
            end if
c     tal 20000705
c     tal 20000627

c     tal 20000705
c        allows turning off of step output files for ASA.
         else if (cmdstr(1:6) .eq. 'NOSTEP') then
            if (.not. ASA) goto 999
            noASAstepout = .true.
c     tal 20000705

c     dmm 20010612
c     if 'BARKER' keyword used, turn on printing of *.rke, *.dens and
c     TSkey.dat files for Barker MEQ code  input.  Note this requires
c     barker.inp file.

         else if (cmdstr(1:6) .eq. 'BARKER') then
            barker = .true.
c     dmm 20010612

c----------------------------------------------------------------------c
c *** RRKM 1 k(E) Thread *** START                                     c
c----------------------------------------------------------------------c
c     dmm 20011019 RRKM option:  allows substitution of one isom rate
c     constant as an RRKM rate constant (hack job, for H + cycloalkenes
c     test).  

         else if (cmdstr(1:4).eq. 'RRKM') then
            RRKM = .true.
            read(linput,*,err=999) i1RRKM, i2RRKM, idgRRKM
c----------------------------------------------------------------------c
c *** RRKM 1 k(E) Thread *** END                                       c
c----------------------------------------------------------------------c
            

c     dmm 20010702
         else if (cmdstr(1:6).eq.'STEADY') then
            steady = .true.
c     dmm 20010702

c     dmm 20020314 
         else if (cmdstr(1:3).eq.'XMG') then
            XMG = .true.

c     dmm 20010813
c     These lines cause CHEMDIS to read the collisional energy transfer
c     energy as deltaEdown, in cm-1, rather than deltaEall in cal/mol

         else if (cmdstr(1:6).eq.'DEDOWN') then
            dedcm = .true.

c     tal 20000814
c         allows choice of RMIN criteria for CHEMACT
c         default: sum of isomer fluxes < RMIN
c         nosumisom: all unique fluxes < RMIN
         else if (cmdstr(1:9) .eq. 'NOSUMISOM') then
            if (.not. ASA)  goto 999
            nosumisom = .true.
c     tal 20000814

c       chemact
         else if (cmdstr(1:4).eq.'CHEM') then
            chemact = .true.
        
c       dissoc
         else if (cmdstr(1:4).eq.'DISS') then
            dissoc = .true.

c     dmm 20000510
c     dissociation on just one input well (avoid cycling through wells)
c     inpWspec holds well number for dissoc calculation
c     tal 20000705 inpWspec now holds well identifier (species name)
         else if (cmdstr(1:8).eq.'INPWONLY') then
            inpwonly = .true.
            read(linput,'(a)',err=999) inpWspec
            call ljust(inpWspec)
            call upcase(inpWspec)
c     dmm 20000510


c       fitglobal
         else if (cmdstr(1:4).eq.'FITG') then
            fitglobal = .true.
        
c       surpress rate table output
         else if (cmdstr(1:5).eq.'NOTAB') then
            table = .false.

c	  write only k's included in the pdep network to chemdis-xmg.out (pey 29/3/04)
c	  only has an effect when printing k's at multiple values of T,P
	  else if (cmdstr(1:6).eq.'INCLKS') then
		  inclks = .true.

c	  fit k(T,P) to a Chebyshev polynomial (pey 4/6/04)
	  else if (cmdstr(1:9).eq.'CHEBYSHEV') then
		  read(linput,*,err=999) nchebT, nchebP
		  if (nchebT.gt.mxChebT) then
		     write(lecho,'(A)') 'Max chebyshev order in T exceeded!'
			 stop
	      elseif (nchebP.gt.mxChebP) then
	         write(lecho,'(A)') 'Max chebyshev order in P exceeded!'
		     stop
		  endif
	      cheby = .true.

c       change default fit option
         else if (cmdstr(1:4).eq.'TROE') then
            Fopt = 'Troe'
        
c       change default fit option
         else if (cmdstr(1:3).eq.'SRI') then
            Fopt = 'SRI'
        
c       show lineshape parameters (don't invoke unless full p range)
         else if (cmdstr(1:4).eq.'SHOW') then
            showline = .true.

c       use original calculation method for stabilization
         else if (cmdstr(1:3).eq.'RKO') then
            rkold = .true.

c       fitrange
         else if (cmdstr(1:4).eq.'FITR') then
            fitrange = .true.
        
c       norot - supress 1 active rotation
         else if (cmdstr(1:3).eq.'NOR') then
            nRot = 0
        
c       change default integration interval (kcal)
         else if (cmdstr(1:3).eq.'INT') then
            read(linput,*,err=999) dE

         else if (cmdstr(1:5).eq.'BSGS') then
            read(linput,*,err=999) bsgs
c       input rate constant
         else if (cmdstr(1:5).eq.'INPUT') then
            read(linput,*,err=999) Ain,rNin,alphaIn,Ein
            readin = .true.
          
c       mass for all isomers  
         else if (cmdstr(1:4).eq.'MASS') then
            read(linput,*,err=999) rmass
            readma = .true.
        
c       temperature values 
         else if (cmdstr(1:4).eq.'TEMP') then
            read(linput,*,err=999) ntemps,(temp(it), 
     2         it = 1,min(ntemps,mxtpts))
     
            if (ntemps.gt.mxtpts) then
               write(lecho,120) mxtpts
120            format(/,' ERROR: max # of ',i2,' temperatures ',
     2            'exceeded.')
               stop
            endif

c	   accept temperature input in form: Tmin, Tmax, nT  (PEY 29/3/04)
		else if(cmdstr(1:6).eq.'TRANGE') then
			read(linput,*,err=999) tmin, tmax, ntemps

			if (ntemps.gt.mxtpts) then
				write(lecho,120) mxtpts
				stop
			endif
		    trange = .true. ! fill temp array below

c       pressure values 
         else if(cmdstr(1:4).eq.'PRES') then
            read(linput,*,err=999) npres,(pres(1,ip),
     2         ip = 1,min(npres,mxPpts))
     
            if (npres.gt.mxPpts) then
               write(lecho,125) mxPpts
125            format(/,' ERROR: max # of ',i2,' pressures ',
     2            'exceeded.')
               stop
            endif
            fillp = .true.
c          later we must fill out entire p array

c	   accept pressure input in the form: Pmin, Pmax, nP (PEY 29/3/04)
		else if(cmdstr(1:6).eq.'PRANGE') then
			read(linput,*,err=999) pmin, pmax, npres

			if (npres.gt.mxPpts) then
				write(lecho,125) mxPpts
				stop
			endif
		    prange = .true. ! fill pres array below
			fillp = .true.

c        correlated t,p option
         else if (cmdstr(1:4).eq.'CORR') then
            read(linput,*,err=999) ntemps,npres
            if (ntemps.gt.mxtpts) then
               write(lecho,120) mxtpts
               stop
            endif
            if (npres.gt.mxPpts) then
               write(lecho,125) mxPpts
               stop
            endif
            do 130 it = 1,ntemps
               read(linput,*,err=999) temp(it),(pres(it,ip),ip=1,npres)
130         continue
        
c       default freq parameters
         else if (cmdstr(1:4).eq.'FREQ') then
            read(linput,*,err=999) nfdef, (freqdef(ifrq),degendef(ifrq),
     2         ifrq = 1,min(nfdef,mxfreqs))
     
            if (nfdef.gt.mxfreqs) then
               write(lecho,145) mxfreqs
145            format(/,' ERROR: max # of ',i2,' ind. frequencies ',
     2            'exceeded.')
               stop
            endif
            readdf=.true.
        
c       default collision parameters
         else if (cmdstr(1:5).eq.'PARAM') then
            read(linput,*,err=999) sigdef,ekdef
            readdc = .true.
 
c       option for colliders
         else if (cmdstr(1:4).eq.'COLL') then
            ncolls = ncolls+1
            if (ncolls.le.mxcolls) then
               read(linput,'(a)') instr
               call ljust(instr)
               collNM(ncolls) = instr
c********12/6/95   comment out next line
c               call upcase(collNM(ncolls))
               if (collNM(ncolls)(1:1).eq.'!') then
                  read(linput,*,err=999) xcoll(ncolls),rmass2(ncolls),
     2               sig2(ncolls),ek2(ncolls),deltaE(ncolls)
               else
                  read(linput,*,err=999) xcoll(ncolls)
                  call lookup(collNM(ncolls),rmass2(ncolls),
     2               sig2(ncolls),ek2(ncolls),deltaE(ncolls),noinfo) 
                  if (noinfo) then
                     write(lecho,580) collNM(ncolls)
580                  format(/,' ERROR: No info on COLLIDER: ',a)
                     stop
                  endif 
               endif

               xtot = xtot + xcoll(ncolls)
            else
               write(lecho,140) mxcolls
140            format(/,' ERROR: max # of ',i2,' colliders ',
     2            'exceeded.')
               stop
            endif

c          set default options           
            consbeta(ncolls) = .false.
            dEExp(ncolls) = 0.d0
	    
c          checks if beta or exponent specified;
c          if not, return to main loop with command string loaded

            read(linput,'(a)',end=450) instr
            call ljust(instr)
            cmdstr = instr
            call upcase(cmdstr)

c          beta option	       
            if (cmdstr(1:4).eq.'BETA') then
               read(linput,*,err=999) beta(ncolls)
               consbeta(ncolls) = .true.

c          exponent option             
            elseif (cmdstr(1:3).eq.'EXP') then
               read(linput,*,err=999) dEExp(ncolls)

            else 
               goto 110
            endif
        
c       get Well name, parameters       
         else if (cmdstr(1:4).eq.'WELL') then
            read(linput,'(a)') instr
            call ljust(instr)
            Wname = instr
            call upcase(Wname)
            
c          check to see if iWell aready asigned for this isom; if not
c          increase nWells
            iWell = 0
            if (nWells.gt.0) then 
c for debugging ...
c               write(*,*) Wname
               do 150 jWell=1,nWells
                  if (Wname.eq.ISname(jWell)) then
                     iWell = jWell
c                     write(*,*) ISname(jWell), iWell, jWell
                  endif
150             continue
            endif
            if (iWell.eq.0) then
               nWells = nWells+1
               iWell = nWells
               if (nWells.gt.mxWells) then
c                  write(*,*) "xxx"
                  write(lecho,160) mxWells
160               format(/,' ERROR:  max # of Wells of ',i3,
     $               ' exceeded.')
                  pause
               endif
               ISname(iWell) = Wname
            endif
            
c          we now have iWell                    
            PDname(iWell,0) = ISname(iWell)
          
c          entry point for Well loop
200         continue
               read(linput,'(a)',end=450) instr
               call ljust(instr)
               cmdstr = instr
               call upcase(cmdstr)

c             comment or blank: do nothing
               if ((cmdstr(1:1).eq.'#').or.(cmdstr(1:1).eq.' ')) then
                  goto 200

c             get reactant channel name, rate               
               else if (cmdstr(1:4).eq.'REAC') then
                  numRC = numRC+1
                  inpWell = iWell
                  nProds(iWell)=nProds(iWell)+1
                  if (nProds(iWell).gt.mxProds) then
                     write(lecho,260) mxProds,ISname(iWell)
cjmg fix
c260                  format(/,' ERROR:  max # of Products for one',
c     2                  ' Well of ',i2,' exceeded on ',a'.')
260                  format(/,' ERROR:  max # of Products for one',
     2                  ' Well of ',i2,' exceeded on ',a,'.')
                     stop
                  endif
                  inpchan = nProds(iWell)
                  read(linput,'(a)') instr
                  call ljust(instr)
                  PDname(iWell,nProds(iWell)) = instr
                  call upcase(PDname(iWell,nProds(iWell)))
                  read(linput,*,err=999) AProd(iWell,nProds(iWell)),
     2               rNProd(iWell,nProds(iWell)),
     3               alphaProd(iWell,nProds(iWell)),
     4               EProd(iWell,nProds(iWell))
          
c             get Product channel name, rate
               else if (cmdstr(1:4).eq.'PROD') then
                  nProds(iWell)=nProds(iWell)+1
                  if (nProds(iWell).gt.mxProds) then
                     write(lecho,260) mxProds,ISname(iWell)
                     stop
                  endif
                  read(linput,'(a)') instr
                  call ljust(instr)
                  PDname(iWell,nProds(iWell)) = instr
                  call upcase(PDname(iWell,nProds(iWell)))
                  read(linput,*,err=999) AProd(iWell,nProds(iWell)),
     2               rNProd(iWell,nProds(iWell)),
     3               alphaProd(iWell,nProds(iWell)),
     4               EProd(iWell,nProds(iWell))
      
c             get isomerization channel name, rate    
               else if (cmdstr(1:4).eq.'ISOM') then
                  isomers = .true.
                  read(linput,'(a)') instr
                  call ljust(instr)
                  Wname = instr
                  call upcase(Wname)
                  
c                again, check if ioWell already assigned
                  ioWell = 0
                  do 300 jWell=1,nWells
                     if (Wname.eq.ISname(jWell)) ioWell = jWell
300               continue
                  if (ioWell.eq.0) then
                     nWells = nWells+1
                     if (nWells.gt.mxWells) then
                        write(lecho,160) mxWells
                        stop
                     endif
                     ioWell = nWells
                     ISname(ioWell) = Wname
                  endif
                  
c                we now have proper ioWell
                  if (ioWell.eq.iWell) then
                     write(lecho,310) ISname(iWell)
                     write(lecho,*) iWell, ioWell
310                  format(/,' ERROR:  ',a,' isomerizes to itself.')
                     stop
                  endif
                 
c                everything ok so continue input
                  read(linput,*,err=999) Aisom(iWell,ioWell),
     2               rNisom(iWell,ioWell),alphaIsom(iWell,ioWell),
     3               Eisom(iWell,ioWell)
                         
c             override default freq parameters - freq
               else if (cmdstr(1:4).eq.'FREQ') then
                  read(linput,*,err=999) nfreqs(iWell),
     2               (freq(iWell,ifrq),degen(iWell,ifrq),
     3               ifrq = 1,min(nfreqs(iWell),mxfreqs))
     
                  if (nfreqs(iWell).gt.mxfreqs) then
                     write(lecho,145) mxfreqs
                     stop
                  endif
                  readfr(iWell) =.true.
                  FrRtRd(iWell) = .false.

c     20010808 dmm
c     read in frequency file and no. of free rotors 
               else if (cmdstr(1:7).eq.'FR-READ') then
                  call VibRotRead(iWell)
                  readfr(iWell) = .true.
                  FrRtRd(iWell) = .true.

c             override default collision parameters
               else if (cmdstr(1:5).eq.'PARAM') then
                  read(linput,*,err=999) sig(iWell),ek(iWell)
                  readcp(iWell) =.true.

c             terminate Well loop; no pass        
               else if (cmdstr(1:3).eq.'END') then
                  goto 100 

c             terminate Well loop; pass command to upper level
               else
                  goto 110

               endif

            goto 200
c          end Well loop
   
c          terminate input loop and allow remainder of file for comments
            else if (cmdstr(1:4).eq.'COMM') then
               goto 450 

            else
               write(lecho,360) instr
360               format(/,' ERROR: command not recognized: ',a,'.')
               stop

            endif 

         goto 100
450      continue       
c    end command loop

c----------------------------------------------------------------------c
c *** RRKM 1 k(E) Thread *** START                                     c
c----------------------------------------------------------------------c
c     dmm 20011018 entry point for separately reading TS data for 1 k(E)
c     RRKM substitution ...
         if(RRKM) then
            iTS = nWells + 1
            fname=ISname(i1RRKM)
            isp = INDEX(fname, ' ')
            if((isp.ne.0).and.(isp.le.(LEN(fname)-2))) then
               fname(isp:(isp+1))='TS'
            else
               fname((LEN(fname)-1):LEN(fname)) = 'TS'
            endif
            ISname(iTS) = fname
            call VibRotRead(iTS)
            FrRtRd(iTS) = .true.
         endif
c----------------------------------------------------------------------c
c *** RRKM 1 k(E) Thread *** END                                       c
c----------------------------------------------------------------------c

c    all done reading
c    fill out default Well parameters if not specified

      do 550 iWell = 1,nWells
         if (.not.readfr(iWell)) then
            if (readdf) then
               nfreqs(iWell) = nfdef
               do 510 ifrq = 1,nfreqs(iWell)
                  freq(iWell,ifrq) = freqdef(ifrq)
                  degen(iWell,ifrq) = degendef(ifrq)
510            continue
            else
               write(lecho,530) ISname(iWell)
530            format(/,' ERROR: FREQUENCIES not specified for ',a)
               stop
            endif
         endif
            
         if (.not.readcp(iWell)) then
            if (readdc) then
               sig(iWell) = sigdef
               ek(iWell) = ekdef
            else 
               write(lecho,540) ISname(iWell)
540            format(/,' ERROR: COLLISION parameters not ',
     2            'specified for ',a)
               stop
            endif
         endif
550   continue

c	begin pey (4/6/04)
c	fill in temp and pres arrays if TRANGE and PRANGE keywords used
c	(also fill in the tcheb and pcheb arrays if using Chebyshev fitting)

	if(cheby) then !space points at Gauss-Chebyshev points in 1/T and log10(P)

		if (trange.and.prange) then

			if (ntemps.lt.nChebT .or. npres.lt.nChebP) then
			   write(lecho,'(2A)') 'Error: Number T and P points ',
     >           'specified in TRANGE,PRANGE must be >= CHEBYSHEV vals.'
	           stop
			endif

      		do it=1,ntemps
				tcheb(it) = cos(pi*(float(it) - 0.5)/float(ntemps))
				temp(it) = 1.0/((tcheb(it)*(1.0/tmin - 1.0/tmax)
     >                           + 1.0/tmin + 1.0/tmax)/2.0)
			enddo

			do ip=1,npres
				pcheb(ip) = cos(pi*(float(ip) - 0.5)/float(npres))
				pres(1,ip)=10.0**((pcheb(ip)*(log10(pmax)-log10(pmin))
     >                           + log10(pmin) + log10(pmax))/2.0)
			enddo

		else
			write(lecho,'(2A)') 'TRANGE and PRANGE keywords must be',
     >		 ' used with CHEBYSHEV keyword.' 
              stop
		endif
      else ! not cheby
		if(trange) then ! space temperature points linearly
			deltat = (tmax-tmin)/float(ntemps-1)
			do it=1,ntemps
				temp(it) = tmin + deltat*float(it-1)
			enddo
		endif
		if(prange) then ! space pressure points linearly
			deltap = (pmax-pmin)/float(npres-1)
			do ip=1,npres
				pres(1,ip) = pmin + deltap*float(ip-1)
			enddo
		endif
	endif

c	end pey

c     fill out remainder of p matrix
      if (fillp) then
         do 570 it = 1,ntemps
            do 570 ip = 1,npres
               pres(it,ip) = pres(1,ip)
570      continue
      endif
      
c    finish collider specification and recompute total mole fractions
      if (ncolls.gt.0) then
c       remormalize x's
         do 620 icoll = 1,ncolls
            xcoll(icoll) = xcoll(icoll)/xtot
620      continue

c    otherwise if no input do default option
      else
         collNM(1) = 'N2'
         xcoll(1) = 1.
         ncolls = 1
         call lookup(collNM(1),rmass2(1),sig2(1),ek2(1),deltaE(1),
     2      noinfo) 
         if (noinfo) then
             write(lecho,580) collNM(1)
             stop
         endif
         collNM(1) = 'DEFAULT'
         beta(1) = 0.2d0
         consbeta(1) = .true.
      endif
         
c    checking for other errors
      if (chemact.and.(.not.readin)) then
         write(lecho,660) 
660      format(/,1x,'ERROR: INPUT rate required for CHEMACT ',
     2      'calculation.')
         stop
      endif
      if (chemact.and.(numRC.ne.1)) then
         write(lecho,670) numRC
cjmg fix
c670      format(/,1x,'ERROR: ',i2,' REACTANT channels specified '
c     2      'for CHEMACT calculation.')
670      format(/,1x,'ERROR: ',i2,' REACTANT channels specified '
     2      ,'for CHEMACT calculation.')
         stop
      endif
      if (nWells.lt.1) then
         write(lecho,680)
680      format(/,1x,'ERROR: no Well info specified.')
         stop
      endif
      
      return

c    flag read errors      
999   continue
      write(lecho,1000) instr
cjmg fix
c1000  format(/,' ERROR:  reading input file after ',a'.')
1000  format(/,' ERROR:  reading input file after ',a,'.')
      stop
      
      end
      
c *************** [add'l routines removed to cdsets.f]*********
      

c ECHOPARAMS******************
      subroutine echoparams(lecho)
      
      implicit none
      include 'cdparams.fh'
      include 'cdwell0.fh'
      include 'cdwell1.fh'
      include 'cdwell2.fh'
      include 'cdisprop.fh'
      include 'cdcolls.fh'
      include 'cdlabels.fh'
      include 'cdcontrl.fh'

c    local variables
      integer iWell,ioWell,iProd,ifrq,icoll,n,iter
      integer lecho,nloops,nstop,Niter
c    end declarations

      write(lecho,10) caseStr
10    format(/,1x,'TAG ID: ',a)

      write(lecho,40) rmass
40    format(/,1x,'ISOMER MASS: ',f5.2,' amu')

      write(lecho,75)
75    format(/,1x,'FREQUENCIES (cm-1)')

      do 100 iWell = 1,nWells
      
c       we output 3 on a line
         if (mod(nfreqs(iWell),3).eq.0) then
            Niter = nfreqs(iWell) / 3
         else
            Niter = nfreqs(iWell) / 3 + 1
         endif
         
         write(lecho,80) iWell,ISname(iWell),
     2      (degen(iWell,ifrq),freq(iWell,ifrq),
     3      ifrq = 1,min(3,nfreqs(iWell)))
80       format(1x,i2,2x,a20,2x,3(2x,f5.2,' x',f7.1)) 

c       no trip do's
         do 90 iter = 2,Niter       
           write(lecho,85) (degen(iWell,ifrq),freq(iWell,ifrq),
     2         ifrq = 3*iter-2,min(3*(iter),nfreqs(iWell)))
85          format(17x,3(2x,f5.2,' x',f7.1)) 
90       continue   
100   continue

      if (nRot.ne.0) write (lecho,120)
120   format(/,1x,'One active rotation assumed.')

      if (chemact) then
         nloops = (nWells-1)/5 +1
         write(lecho,240) 
240      format(/,1x,'ISOMERIZATION CHANNEL A''S')
         do 290 n = 1,nloops
            nstop = min(5*n,nWells)
            write(lecho,250) (isname(ioWell), ioWell = 5*n-4,nstop)
250         format(1x,10x,5(2x,a10))
            do 280 iWell = 1,nWells
               write(lecho,270) ISname(iWell),(Aisom(iWell,ioWell),
     2            ioWell = 5*n-4,nstop)
270            format(1x,a20,5(2x,1pe10.4))
280         continue
290      continue
      
         write(lecho,300) 
300      format(/,1x,'ISOMERIZATION CHANNEL N''S')
         do 330 n = 1,nloops
            nstop = min(5*n,nWells)
            write(lecho,250) (isname(ioWell), ioWell = 5*n-4,nstop)
            do 320 iWell = 1,nWells
               write(lecho,310) ISname(iWell),(rNisom(iWell,ioWell),
     2            ioWell = 5*n-4,nstop)
310            format(1x,a20,5(2x,1pe10.4))
320         continue
330      continue

         write(lecho,335) 
335      format(/,1x,'ISOMERIZATION CHANNEL alpha''S')
         do 360 n = 1,nloops
            nstop = min(5*n,nWells)
            write(lecho,250) (isname(ioWell), ioWell = 5*n-4,nstop)
            do 350 iWell = 1,nWells
               write(lecho,340) ISname(iWell),(alphaIsom(iWell,ioWell),
     2            ioWell = 5*n-4,nstop)
340            format(1x,a10,5(2x,1pe10.4))
350         continue
360      continue

         write(lecho,380) 
380      format(/,1x,'ISOMERIZATION CHANNEL E''S (kcal)')
         do 400 n = 1,nloops
            nstop = min(5*n,nWells)
            write(lecho,250) (isname(ioWell), ioWell = 5*n-4,nstop)
            do 390 iWell = 1,nWells
               write(lecho,370) ISname(iWell),(Eisom(iWell,ioWell),
     2            ioWell = 5*n-4,nstop)
370            format(1x,a10,5(2x,f7.3,3x))
390         continue
400      continue

      endif

      write(lecho,420)
420   format(/,1x,'PRODUCT CHANNELS ',24x,
     2      '   A   ',4x,'   n',4x,'  alpha    ',1x,' E (kcal)')     
      do 450 iWell = 1,nWells
         if (nProds(iWell).gt.0) then
            do 440 iProd = 1,nProds(iWell)
               write(lecho,430) iWell,iProd,ISname(iWell),
     2            PDname(iWell,iProd),AProd(iWell,iProd),
     3            rnProd(iWell,iProd),alphaProd(iWell,iProd),
     4            EProd(iWell,iProd)
cjmg fix
c430            format(1x,i2':',i2,2x,a,' = ',a,1pe10.4,2x,
c     2            0pf6.3,2x,1pe9.3,3x,0pf8.3)
430            format(1x,i2,':',i2,2x,a,' = ',a,1pe10.4,2x,
     2            0pf6.3,2x,1pe9.3,3x,0pf8.3)
440         continue
         endif
450   continue
 
      write(lecho,580)
cjmg fix
c580   format(/,1x,'COLLIDERS',5x,'  x  ',4x,'  m  ',4x,' sigma',3x,
c     2   '  e/k ',3x,'dEavg(300)',2x,'exp(dEd)',3x' beta')
580   format(/,1x,'COLLIDERS',5x,'  x  ',4x,'  m  ',4x,' sigma',3x,
     2   '  e/k ',3x,'dEavg(300)',2x,'exp(dEd)',3x,' beta')
      do 700 icoll = 1,ncolls      
         if (consbeta(icoll)) then
            write(lecho,620) collNM(icoll),xcoll(icoll),rmass2(icoll),
     2         sig2(icoll),ek2(icoll),beta(icoll)
620         format(1x,a,3x,f5.3,4x,f5.1,4x,f6.2,4x,f6.2,5x,21x,f5.3)
         else
            write(lecho,650) collNM(icoll),xcoll(icoll),rmass2(icoll),
     2         sig2(icoll),ek2(icoll),deltaE(icoll),dEExp(icoll)
650         format(1x,a,3x,f5.3,4x,f5.1,4x,f6.2,4x,f6.2,5x,f6.1,
     2         5x,f6.3)
         endif
700   continue

      if (rkold) then
         write(lecho,710) 
710      format(/,' Stabilization rates computed with RKOLD option.')
      endif
      
      return
      end

c***********************************************************************
c***********************************************************************
c     Subroutine VibRotRead
c     
c     dmm 20010808
c
c     This subroutine is called from getparams when FFILE keyword
c     appears where the "FREQ" keyword normally would in fort.10.  It
c     reads a set of frequencies from an input file of the name '
c     ISname(iWell)'.dat.  The format of the file should be
c
c     line 1  no. of freq   no. of free rotors   no. of hindered rotors
c     lines 2-xx frequency table as per sumathy's tstuni.f
c     input, typically a 3-column table of frequencies (free format)
c
c     number of hindered rotors is read in but ignored for now;
c     eventually, should read in fourier forms for the surfaces,
c     geometry/external modes, and figure out cross-correlations?, then
c     solve Sch. eq'n for energy levels, somehow send these to density
c     -of-states routines
c
c     ikrot = flag to indicate whether k-rotor MOMI will be read or not
c***********************************************************************
c***********************************************************************

      SUBROUTINE VibRotRead(iWell)
      
      implicit none
      include 'cdparams.fh'
      include 'cdisprop.fh'
      include 'cdlabels.fh'
      include 'cdcontrl.fh'

      integer iWell, ilen, isp, ifin
      integer ihrot, ifrot, i, ikrot
      logical isthere
      character fname*20
      
c     Construct the correct input file name and open

      fname = ISname(iWell)
      isp = INDEX(fname, ' ')
      if(isp.ne.0) then
         ifin = isp - 1
      else
         ifin = LEN(fname)
      endif

c     see if this file exists before proceeding
      inquire(file = ''//fname(1:ifin)//'.dat', EXIST = isthere)
      if(.not.isthere) then
         write(*, 110), fname(1:ifin)
         stop
      endif

c     open and read data
      open(100, file = ''//fname(1:ifin)//'.dat')
      read(100,*) nFreqs(iWell), n1drot(iWell), ihrot, ikrot
      if(nFreqs(iWell).gt.mxFreqs) then
         write(*, 111), ISname(iWell)
         stop
      endif

      read(100,*) (freq(iWell,i), i=1,nFreqs(iWell))

c----------------------------------------------------------------------c
c *** RRKM 1 k(E) Thread *** START                                     c
c----------------------------------------------------------------------c
c     if ikrot flag is set, then read the moment of inertia for this
c     rotor (normally should only be needed for the RRKM hack option).
c     This will be in cm-1

      if(ikrot.eq.1) then
         read(100,*) rotk(iWell)
      else
         rotk(iWell) = -1.0
      endif
c----------------------------------------------------------------------c
c *** RRKM 1 k(E) Thread *** END                                       c
c----------------------------------------------------------------------c

c     placeholder:  may need to set degen vector to 1.0 values
c     everywhere ...

 110  format(1X,'<<ERROR>>: failed to find freq-rotor input file '
     $     ,A20, 1x, ' ..., stopping')

 111  format(1X, '<<ERROR>>: number of freqs read for ', A20,' more ',
     $   'than mxFreqs')

      end
 
c***********************************************************************
c***********************************************************************


c $Id$
c $Author$
c $Date$
c $Revision$
c $Log$
c Revision 1.1  2007-02-20 23:10:23  sandeeps
c Initial revision
c
c Revision 1.17  2001/10/18 20:58:35  dmmatheu
c before changes to allow 1 RRKM k(E) table to be calculated a la BARKER
c option ...
c
c Revision 1.16  2001/10/18 19:52:16  dmmatheu
c before editing to generate a k(E) from RRKM estimate (hack)
c
c Revision 1.15  2001/09/20 00:00:48  dmmatheu
c before changing Wname length to 15 from 20 (originally 10???)  trying
c to get H+cyclohexene to work.  QRRK.input file generated with 15
c character strings for the WELL lines but 20 character strings for the
c ISOMER lines (!!!!) ARRRGH! hard to fix.
c
c Revision 1.14  2001/09/19 00:03:51  dmmatheu
c during changes to try to get longer wellnames for debugging in
c H+cyclohexene
c
c Revision 1.13  2001/08/08 19:32:20  dmmatheu
c Version before adding VibRotRead, to read a list of frequencies and
c the number of free rotors (hindered rotors not yet supported) for the
c well from a file.
c
c Revision 1.12  2001/06/13 14:07:02  dmmatheu
c After changes to turn Barker file outputs on/off, major debug of
c mxWells and mxProds (zeroing loops of MakeShadow were exceeding array
c bounds).
c
c Added barker logical variable to this file, and lines to check for
c BARKER keyword in input file, which when used causes barker files to
c be printed
c
c Revision 1.11  2000/08/29 00:35:19  dmmatheu
c after tlada changes -- added absolute Rmin criteria, new printing
c line.
c
c Revision 1.10  2000/08/18 13:13:34  dmmatheu
c tlada added NOSUMISOM keyword to specify that individual pathway
c fluxes and not the sum of fluxes to all nonincluded isomers should be
c used to decide whether to add more isomers.
c
c Revision 1.9  2000/07/07 01:01:10  dmmatheu
c Checking in Tom's changes (lots) ... rmin criteria flag, input well
c specification fix (specify by well name for dissoc, not number)
c
c Revision 1.7  2000/06/27 18:20:07  dmmatheu
c before ASA additions ...
c
c Revision 1.6  2000/06/23 17:39:35  dmmatheu
c after Tom Lada merge
c
c Revision 1.5  2000/06/23 17:37:05  dmmatheu
c before merge with Tom Lada code
c
c Revision 1.4  2000/03/30 18:05:00  dmmatheu
c set Beyer-Sweinhart grainsize (bsgs) default to 10 cm-1
c
c Revision 1.3  2000/03/30 17:59:56  dmmatheu
c added bsgs variable and cdnos.fh
c
c Revision 1.2  2000/02/25 01:32:58  dmmatheu
c before adding gamma table check ...
c
