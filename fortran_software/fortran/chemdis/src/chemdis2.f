      program chemdis2


cc
c    4/29/95  modified cdstab to print out density of states AMD
c    version 5.15  <ayc 4/94>
c    compiled with -mp -mips2 -limsl -limslblas on SGI;
c    fitting seems to be working but not checked extensively <ayc 4/94>

c    NOTE: there routines are somewhat sloppy with regards to case;
c          compile case insensitive  (e.g. FOLD TO LOWER CASE)

c*********program chronology:
 
c <v1>    chemact for arbitrary arrangement of Wells and Products,   <ayc 10/92>
c         plus dissoc model built-in; computes high and low pressure 
c         limits; also lindemann-like F forms

c (v1.5>  new stabilization rate model following gunther et al       <ayc 11/92>

c <v2>    corrected methodology to be consistent with paper          <ayc  1/93>
c         new gamma and matrix routines                              <ayc  3/93>

c <v3>    implemented n-freq model with recursive subroutines        <ayc  4/93>
c         (removed 3f routines); slight correction in calrate  

c         changed chemact convergence criteria to 0.2% of rate of    <ayc  4/93>
c         least E-accessible Product channel in hp limit -
c         tweaked dissoc to same consistency

c <v3.5>  fitting for F routines added - gaussian with polynomomial  <ayc  4/93>
c         parameters or Troe and SRI forms if requested; increased
c         default P points to 61 (1,2,5 from 1e-10 to 1e10);
c         added routines for fitting high and low p limits plus
c         routines for reconstructing solution from fit coefficients;
c         replaced old lsq fitting procedure

c <v4>    New Dissoc model - dissoc's allowed to all channels in     <ayc  5/93>
c         system!  Note: we observe that for diss calcs at low E's where
c         there are no Product channels (only isomerization) the stabili-
c         zation Products have a different p-dependence in the low p limit.
c         we thus must factor the low p cases into two separate cases;
c         with ip = 0 the normal case, and ip = -1 the weird case.  we note
c         this weirdness by setting flag noDExit(iWell) = .true., with the 
c         understanding that this only affects the iProd = 0 (stab. channel).
c         For these cases there are two low P limit equations to be fit;
c         the 'normal' one (when there are exit channels), and an additional
c         term which must be carried along and propagated into the reduced 
c         pressure, etc. 
c         [WE STILL HAVE FITTING PROBLEMS ON F FUNCTIONS IN THESE CASES]

c         surprisingly, gaussian elim. works best for the matrix solver
c         (ldu becomes ill behaved under extreme conditions, and sdv
c         isn't happy at all!!)   [use GetB1 in cdaux]               <ayc  5/93>          

c         NOTE: to avoid singularities, we allow in the matrix equation
c         small leakage out, corresponding to an [M] = Rleak of 1e-70
    
c <v5>    switched to INCLUDE file headers for all common blocks     <ayc  8/93>
c         ('standard', non-standard fortran); 
c         also switched to implicit NONE (from char or logical) for 
c         better var checking; also relatively standard;
c         abandon 1st-6 significant character rule for identifiers

c         overhauled cdchem, cddiss, cdaux, cdrecur & cdsubs to      <ayc  8/93>

c         a: interpolate pho based on n-freq model; integration interval
c            now inputted as separate parameter (keyword INT (def 1 kcal)
c         b: calculate all rates based on A, n and E; note this changes 
c            input file format
c         c: corrected n-freq model by increasing sums to E+nhv/2, with
c            int truncation on innermost E term
c         d: added 1-D active rotation in rho(E);
c            use NOROT to suppress

c <v5.1>  compiled on SGI machine: had to rearrange some common blocks
c         to avoid misalignment warnings (reals/ints/log/char)       <ayc  1/94>

c <v5.15> added exponent to dEDown                                   <ayc  3/94>

c***************** more info
c    NOTE: for SGI version - remove commented call statements
c    to fitF, fitFout and reConst in main below - Ffits not implemented
c    on Mac yet (need non-linear fitter or func. minimizer)
      
c    this file: chemdis2.f
c         this includes the driver program as Well as key sub-
c         routines the user may wish to modify: those controlling
c         the file handling, the default T and P ranges and the
c         library of collision parter info
 
c    link chemdis2.f with:
c         cdgets.f:   reads input file and echos parameters
c         cdsets.f:   sets up variables for chemact and diss calc
c         cdchem.f:   does chemact calculations
c         cddiss.f:   does dissoc calculations
c         cdaux.f:    auxillary routines for above
c         cdsubs.f:   support routines for rates, f(E)'s, etc.
c         cdrecur.f:  recursive routines for rate functions
c         cdstab.f:   support routines for stabilization calcs
c         cdkouts.f:  output routines for rates
c         cdbits.f:   auxillary routines for fitting, matrix solving
c         cdstr.f:    miscellaneous string handling routines

c     and F function routines; F fitting uses IMSL
c         cdkfits.f:  add'l fitting routines for rates
c         cdfouts.f:  output routines for lineshape functions
c         cdffits.f:  fitting routines for f surface - CRAY ver. only
c         cdffuns.f:  fitting functions for f - CRAY ver. only
c         cdrecon.f:  routine to reconstruct fit - CRAY ver. only

c     include header files used for block common declarations
c         cdparams.fh: parameters used to dimension arrays in commons 
c         cdwell0.fh:  globals nWells, nProds, inpWell, and inpChan
c         cdwell1.fh:  inputted A's, E's and I's
c         cdwell2.fh:  misc Well stats needed for calculation
c         cdisprop.fh: isomer properties including freqs, collision parameters
c         cdcolls.fh:  collision partner parameters
c         cdlabels.fh: names of all things
c         cdrange.fh:  T and P variables
c         cdrates.fh:  calculated integrated rates for chemact or dissoc loop
c         cdffunc.fh:  calculated F function and Pr array for single channel 
c         cdffit.fh:   fit parameters on F
c         cdkfit.fh:   arrhenius fit parameters; all channels
c         cdlimfit.fh: fit to high and low p limit for single channel
c         cdcontrl.fh: misc control and option specification variables
c         cdrhovar.fh: storage for density function calculations

c*****start here

c    non-standard implcit, but almost universal; 
c    if problems switch to character or logical
      implicit none
      
c    grab parameter statments for commons
      include 'cdparams.fh'

c     grab only those commons needed for main      
	include 'cdwell0.fh'
      include 'cdrange.fh'
      include 'cdrates.fh'
      include 'cdlimfit.fh'
      include 'cdcontrl.fh'
      include 'cdlabels.fh'

c	pey 4/6/04
	include 'cdcheb.fh'
c	end pey

c     dmm 20000224
c     gamma table common block
      include 'cdgamma.fh'
c     dmm 20000224

c     dmm 20000219
      real tarray(2)
c     dmm 20000219

c    delcare local variables
      character option*8,fstamp*120
      integer iWell,iProd,ip,n
      integer Nluout
      real*8 dE
      logical proceed, isthere
      
c     dmm 20020314 xmg output file indicator
      integer lxmgout, lmtlbout


c    assign lu's; Nluout is in data statement    
      integer linput,lecho,lkout,ldiag,lfout,ltout,lckout,ltabout,
     2   ltablong,lffit,lfsum,lrecon,luarr(11)
      equivalence (luarr(1),lecho),(luarr(2),lkout),(luarr(3),ldiag),
     2   (luarr(4),lfout),(luarr(5),ltout),(luarr(6),lckout),
     3   (luarr(7),ltabout),(luarr(8),ltablong),(luarr(9),lffit),
     4   (luarr(10),lfsum),(luarr(11),lrecon)

      data linput,lecho,lkout,ldiag,lfout,ltout,lckout,ltabout,
     2   ltablong,lffit,lfsum,lrecon, Nluout
c	begin pey 30/3/04 -- changed Nluout so lots of empty output files aren't created
c     3   /10, 20,25,30,35,45,55,65,73,75,85,95,11/ 
     3   /10, 20,25,30,35,45,55,65,73,75,85,95,3/ 
c	end pey
     
c******end declarations ***********

c    default T and P range assigned through block data getrange 
c    these defaults can be overwritten if option specified

c     Begin timing point:  whole code case,  dmm 20000219
      call dtime(tarray)
c     dmm 20000219


c    open and read file header from input file
      open(linput, file = 'fort.10', status = 'old')
      read(linput,'(a)') fstamp
      
c    open all output files and stamp header      
      call initfs(Nluout,luarr,fstamp)
      
c    read all other input parameters from input file and close
      call getparams(linput,lecho,dE)
      close(linput, status = 'keep')
 
c    echo generic parameters to file 
      call echoparams(lecho)

c dmm 20000627 zero out variables for genRhos.f
      call InitGenRhos
c dmm 20000627

c    loop through remainder of program doing chemact calculation first
c    then do dissoc calculations if desired

      do 200 n = 0,nWells
      
c       set switches     
         if ((n.eq.0).and.(chemact)) then    
            option = 'Chemact'
            proceed = .true.

         else if ((n.ge.1).and.(dissoc).and.(.not.inpwonly)) then
            option = 'Dissoc'
            inpWell = n
            inpchan = 0
            proceed = .true.

c dmm 20000510
c     Added block to allow specification of particular input well as
c     source for dissoc calculations (instead of going through all of
c     them).  This is prep for conversion to routine called by XMG.

         else if ((n.gt.0).and.(isname(n).eq.inpWspec(1:len(isname(n))))
     $           .and.(dissoc).and.(inpwonly)) then
            option = 'Dissoc'
            inpWell = n
            inpchan = 0
            proceed = .true.
c dmm 20000510
         else
            proceed = .false.
         endif

         if (proceed) then
         
c          get array idpWell giving depth of complex on isomer chart;
c          also compute Well depths in kcal; and Emax
            write(*,*) 'setcalc'
            call setCalc(lecho,option)  

c tal 20000705 changing for a softer, gentler ASA
c dmm 20000627 temporary flow control
c            if(.not.ASA) then
c               write(*,*) 'echo2'
c               call echo2(lecho,option)
c            else
c               call ASAScreen(option,lecho, ldiag, dE)
c            end if
c dmm 20000627

            call echo2(lecho, option)
            if (ASA) then
c eliminated ASA portion for PC build 20031105
c               call ASAScreen(option, lecho, ldiag, dE)
            else

c          now do actual calculation
               if (option.eq.'Chemact') then
                  write(*,*) 'compChem'

c     dmm 20000411 timing of CHEMACT calculation
c               call dtime(tarray) 
                  call compCHEM(lecho,ldiag,dE)
c               call dtime(tarray)
c     dmm 20000411 end timing
               else
                  write(*,*) 'compDiss'
                  call compDISS(lecho,ldiag,dE)
               endif
            end if
c tal 20000705

c          finish calculation
            write(*,*) 'cleanup'
            call cleanup(option)         
            
            write(*,*) 'output'
c          now we output results; first we print table
            if (table) then
c               call tabout(option,ltabout)
c              if (ASA) call longtabout(option, ltablong)
               
c dmm 20020314 XMG-formatted table output (pey 30/3/04 - replaced tabXMG with rmgout)
               if (XMG.and.(.not.ASA)) then
                  lxmgout = 66
				lmtlbout = 67
                  open(lxmgout, file = 'chemdis-rmg.out')
                  call rmgout(option, lxmgout)
c				open(lmtlbout, file = 'chemdis-mtlb.out')
c				call matout(option, lmtlbout)
               endif
               
            end if

c          now we loop over all reactions and output and or fit
c          remember iProd = 0 is stabilization channel for chemact

c          note common blocks kfit and Ffit only store one reaction's
c          worth of parameters at a time

            do 100 iWell = 1,nWells
               do 80 iProd = 0,NProds(iWell)

c                the boolean addTerm indicates whether we must include
c                secondary low p limit term
                  if ((option.eq.'Dissoc').and.(noDExit(iWell))
     2               .and.(iProd.eq.0)) then
                     addTerm = .true.
                  else
                     addTerm = .false.
                  endif
	       
c             don't output dissoc channels where dissociating isomer
c             is Product - this is a usless rate and we integrated for 
c             an E > 0 so we didn't calculate this correctly anyway
               if (.not.((option.eq.'Dissoc').and.(iWell.eq.inpWell)
     2            .and.(iProd.eq.0))) then

c                   fit limiting arrhenius factors

c                     if ((fitglobal.or.fitrange).and.(nTemps.gt.2)) then
c                        if (addTerm) 
c     2                     call Kfitout(option,lecho,lkout,iWell,
c     3                     iProd,-1)
c                        call Kfitout(option,lecho,lkout,iWell,iProd,0)
c                        call Kfitout(option,lecho,lkout,iWell,iProd,
c     2                     npres+1)
c                     endif

c                   fit arrhenius factors over p range
                     if (fitrange.and.(nTemps.gt.2)) then
                        do 70 ip = 1,npres
                           call Kfitout(option,lecho,lkout,iWell,
     2                        iProd,ip)
                           call CKoutput(option,lckout,iWell,iProd,ip)
70                      continue
                     endif

c                   fit lineshape factors; make sure enough parameters
                     if ((fitglobal).and.(nTemps.gt.3).and.
     2                  (nPres.gt.5)) then
                        call kLimFit(option,lecho,lfsum,iWell,iProd)
                        call calcF(option,iWell,iProd)
   
                        if (showline) then 
                           call Foutput(option,lfout,iWell,iProd)
                           call Toutput(option,ltout,iWell,iProd)
                        endif
     
c      IMPORTANT!!  these subroutines presently only exist in SGI vers.
C                        call FitF(lecho)
C                       call FitFout(option,lffit,lfsum,iWell,iProd)
C                        call reConst(option,lrecon,iWell,iProd)
                     endif

                  endif
80             continue
100         continue
         endif
200   continue

c    close all output files
999   continue
      call closfs(Nluout,luarr)
c     End whole code timing; dmm 20000219
      call dtime(tarray)
c     Setup call for printing out Barker kofE files
      if(barker) then
         if (chemact.or.(dissoc.and.inpWonly)) then
            call PrintBarkerkofE
         endif
      endif
c     Setup call for printing out steady-state calculator inputs
      if(steady) then
         if((chemact.and.dissoc).or.(dissoc.and.(.not.inpwonly))) then
            write(*,*) '<<WARNING>>: STEADY on with both chemact and',
     $           ' dissoc, or multi-dissoc.  Using last calcd'
            option='dissoc'
         else if(chemact) then
            option='chemact'
         else if(dissoc) then
            option='dissoc'
         endif
         call PrintSteadyInput(option)
      endif

      write(*,909), tarray(1)
      write(*,919), tarray(2)
 909  format(1X, 'Total Code User Time:    ', E16.9)
 919  format(1X, 'Total Code System Time:  ', E16.9)
      stop
      end
        

c***INITFS****************************
      subroutine initfs(Nlu,lu,string)
c    opens and clears out all output files and stamps header
      implicit none
 
c    local variables
      integer Nlu,lu(Nlu)
      integer i
      character string*(*),num*2,fname*8
c    end declarations

      do 20 i = 1,Nlu
         write(num,'(i2)') lu(i)
         call concat('fort.',num,fname)
         open(lu(i), file = fname, status = 'unknown')
         close(lu(i), status = 'delete')
         open(lu(i), file = fname, status = 'new')
         if ((lu(i) .eq. 65) .or. (lu(i) .eq. 73)) then
           write(lu(i),15) string
         else
         write(lu(i),'(1x,a)') string
         end if
15       format('"',a,'"')
20    continue

      return
      end
      
c***CLOSFS****************************
      subroutine closfs(Nlu,lu)
c    closes output files
      implicit none
     
c    local variables
      integer Nlu,lu(Nlu)
      integer i
c    end declarations

      do 20 i = 1,Nlu
         close(lu(i), status = 'keep')
20    continue

      return
      end
      
c GETRANGE******************
      block data getrange
      
      implicit none
      include 'cdparams.fh'
      include 'cdrange.fh'

c    local variables
c      integer i,j
c    end declarations

c$$$      data ntemps,npres/13,61/
c$$$      
c$$$      data (temp(i),i=1,13) /300.,400.,500.,600.,800.,1000.,1200.,1500.,
c$$$     2   2000.,2500.,3200.,4000.,4800./
c$$$     
c$$$      data ((pres(i,j),i=1,13),j=1,21) /13*1.d-10,13*2.d-10,13*5.d-10,
c$$$     2   13*1.d-9,13*2.d-9,13*5.d-9,13*1.d-8,13*2.d-8,13*5.d-8,13*1.d-7,
c$$$     3   13*2.d-7,13*5.d-7,13*1.d-6,13*2.d-6,13*5.d-6,13*1.d-5,13*2.d-5,
c$$$     4   13*5.d-5,13*1.d-4,13*2.d-4,13*5.d-4/
c$$$     
c$$$      data ((pres(i,j),i=1,13),j=22,42) /13*1.d-3,13*2.d-3,13*5.d-3,
c$$$     2   13*1.d-2,13*2.d-2,13*5.d-2,13*1.d-1,13*2.d-1,13*5.d-1,13*1.d0,
c$$$     3   13*2.d0,13*5.d0,13*1.d1,13*2.d1,13*5.d1,13*1.d2,13*2.d2,
c$$$     4   13*5.d2,13*1.d3,13*2.d3,13*5.d3/
c$$$     
c$$$      data ((pres(i,j),i=1,13),j=43,61) /13*1.d4,13*2.d4,13*5.d4,
c$$$     2   13*1.d5,13*2.d5,13*5.d5,13*1.d6,13*2.d6,13*5.d6,13*1.d7,
c$$$     3   13*2.d7,13*5.d7,13*1.d8,13*2.d8,13*5.d8,13*1.d9,13*2.d9,
c$$$     4   13*5.d9,13*1.d10/
c$$$      
      end
      
c***LOOKUP*******************************
      subroutine lookup(sname,rmass,sig,ek,deltaE,noinfo)
      implicit none
      
      integer mxdata
      parameter (mxdata = 10) 
      
c    local variables
      logical noinfo
      character sname*10
      integer n
      integer ndata
      real*8 rmass,sig,ek,deltaE
      character species(mxdata)*10
      real*8 smass(mxdata),ssig(mxdata),sek(mxdata),sdeltaE(mxdata)
      
c    do assignments here
      data ndata/6/     
      data species(1) /'N2'/
      data smass(1),ssig(1),sek(1),sdeltaE(1)/28.0, 3.621, 97.5, 830./     
      data species(2) /'AR'/
      data smass(2),ssig(2),sek(2),sdeltaE(2)/40.0, 3.33, 136.5, 630./
      data species(3) /'HE'/
      data smass(3),ssig(3),sek(3),sdeltaE(3)/4.0, 2.6, 10.2, 431./
      data species(4) /'CH4'/
      data smass(4),ssig(4),sek(4),sdeltaE(4)/16.0, 3.746, 141.4, 2100./
      data species(5) /'C3H8'/
      data smass(5),ssig(5),sek(5),sdeltaE(5)/44.0, 4.98, 266.8, 4200./
      data species(6) /'SF6'/
      data smass(6),ssig(6),sek(6),sdeltaE(6)/146.0, 5.13, 222.1, 3400./
      
      noinfo = .true.

      do 10 n = 1,ndata
         if (sname.eq.species(n)) then
            rmass = smass(n)
            sig = ssig(n)
            ek = sek(n)
            deltaE = sdeltaE(n)
            noinfo = .false.
         endif
10    continue

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
c Revision 1.16  2003/02/14 19:55:16  dmmatheu
c before changes to Kfitout calls to eliminate (for XMG fitting
c purposes) the high and low pressure limit calls
c
c Revision 1.15  2001/06/13 14:08:34  dmmatheu
c After changes to turn Barker file outputs on/off, major debug of
c mxWells and mxProds (zeroing loops of MakeShadow were exceeding array
c bounds).
c
c Added PrintBarkerkofE call to this file
c
c Revision 1.14  2000/10/09 23:04:27  dmmatheu
c no change of note
c
c Revision 1.13  2000/07/07 01:13:42  dmmatheu
c tlada changes checked in ...
c 1.  debugging of isomerization loops
c 2.  addition of inpWspec as character (not integer) so as to specify
c single well for dissoc by name
c 3.  debugging of rmin criteria
c 4.  user specification of rmin option
c 5.  user option to supress staged output files
c
c Revision 1.11  2000/06/29 13:34:01  dmmatheu
c tlada changes checked in -- cdkouts, chemdis2 and cdasa.f  for writing
c output files at every step and fixing formatting to work better with
c text editors.
c
c Revision 1.10  2000/06/27 22:13:43  dmmatheu
c added InitGenRhos to allow more elegant control of when a well's
c density of states must be recalc'd.  See genRhos.f.  Passes vinyl+O2
c example test.
c
c Revision 1.9  2000/06/27 18:15:34  dmmatheu
c Before changes to add ASA
c
c Revision 1.8  2000/06/23 17:54:03  dmmatheu
c after Tlada merge
c
c Revision 1.7  2000/06/23 17:52:52  dmmatheu
c before Tom lada merge ... after debugging to fix dissoc problem
c
c Revision 1.6  2000/04/04 21:41:52  dmmatheu
c speedup of chemical activation compiled and tested for vinyl+O2 case.
c calRho replaced by xNOS.  Seems to work fine.
c
c Revision 1.5  2000/03/30 20:47:11  dmmatheu
c Before adding genRho ...
c
c Revision 1.4  2000/02/24 22:28:05  dmmatheu
c before modification to run gamma lookup table
c
c Revision 1.3  2000/02/19 20:26:29  dmmatheu
c whole code timing
c
c Revision 1.2  2000/02/19 20:17:09  dmmatheu
c Before adding timing routines (Perf. An. suggests k(E) calcs as
c bottleneck)
c

