c***********************************************************************
c***********************************************************************
c     Subroutine PrintBarkerDens
c
c     This subroutine prints out "__________.dens" type files for use
c     with Barker's Master Equation Code.  It reads the file "barker
c     .inp" for the parameters which make up the third header line in
c     the barker dens files
c
c***********************************************************************
c***********************************************************************

      SUBROUTINE PrintBarkerDens(maxn)

      implicit none
      include 'cdparams.fh'
      include 'straightBS.fh'
      include 'cdnos.fh'
      include 'cdisprop.fh'
      include 'cdwell0.fh'
      include 'cdlabels.fh'

      integer maxn, ifreq, iWell, imax1, Isize, nbslength, ispace, ic
      integer i, nE
      real*8 Egrain1, Emax2, Viblo, energy, denstates, sumstates
      real*8 xNOScm, Egrain2, Emax1, dum1, dum2
      character bkname(mxWells)*10


c     Open Barker data file

      open(101, file='barker.inp', status = 'old')
      read(101,*) Egrain1, imax1, Emax2, Isize, Viblo
      close(101)

c     Pad extra spaces with underscores

      do 10 iWell=1, nWells
         bkname(iWell) = isname(iWell)
         call PadBlanks(bkname(iWell),'_')
 10   continue

      do 20 iWell = 1, nWells
c     Determine the lowest frequency for this well

         Viblo = 1.0E6
         do 25 ifreq = 1, nFreqs(iWell)
            if (freq(iWell,ifreq).lt.Viblo) then
               Viblo = freq(iWell,ifreq)
            endif
 25      continue

c     Open the well "dens" file and write the headers
         open(102, file=''//bkname(iWell)//'.dens')
         write(102,*) bkname(iWell)
         write(102,*) 
c         write(102,*)
         write(102,901) Egrain1, imax1, Emax2, Isize, Viblo
         write(102,900)
 900     format(5x,'No.',5x,'(cm-1)', 5x,'Density', 5x, 'Sum')
 901     format(1x,f10.1,1x,I10,1x,f10.1,1x,I10,1x,f10.1)

c     Begin first part of double array: 1 to imax1

c     Set Emax1 for WR code
c         Emax1 = dfloat(imax1)*Egrain1
         Emax1 = 4000.0
         energy = 0.0
         do 30 i = 1, imax1

            if(energy.ge.Emax1) then
               call WRSumDens(iWell,energy, sumstates, denstates)
            else
               denstates = xNOScm(energy, Egrain1, iWell)/Egrain1

c     dmm20011102  offset problem fix ... xhroofn(1) is energy = 0!
               nE = nint( (energy)/bsgs) + 1
c     nE = nint( (energy + Egrain1)/bsgs)

               sumstates = xrhoofn(nE,iWell)
            endif
            write(102,903) i, energy, denstates, sumstates

            energy = energy + Egrain1
 903        format(1x, I5, 3x, F9.2, 3x, E10.4, 3x, E10.4)
 30      continue

c     Set Emax1 for WR code
c     Emax1 = energy - Egrain1


c     Begin second part of double array:  imax1+1 to Isize
c     Check whether Emax2 is too large ...

         if (nint(Emax2/bsgs).gt.maxBSindex(iWell)) then
            write(*,*) 
     $ ' ERROR::  PrintBarkerD: Emax2 greater than maximum energy in', 
     $ ' xrhoofn array'
            stop
         endif

         Egrain2 = Emax2/(dfloat(Isize-(imax1+1)))
         energy = 0.0
         do 40 i = imax1+1, Isize


c     WR option -- uses Whitten-Rabinovtich code for energies above the
c     double-array 'split' of the barker files ...

            if(energy.ge.Emax1) then
               call WRSumDens(iWell,energy, sumstates, denstates)
            else
               denstates = xNOScm(energy, Egrain2, iWell)/Egrain2
               nE = nint(energy/bsgs) + 1
c     dmm20011102 offset problem fix (see above)
c            nE = nint((energy + Egrain2)/bsgs)
               sumstates = xrhoofn(nE, iWell)
            endif
            write(102,903) i, energy, denstates, sumstates
            energy = energy + Egrain2
            
 40      continue
c     close the well dens file
         close(102)
 20   continue

      return
      end

c***********************************************************************


c***********************************************************************
c***********************************************************************
c     SUBROUTINE PrintBarkerkofE
c
c     This subroutine prints the Barker-style '.rke' files for each
c     pathway.
c
c***********************************************************************
c***********************************************************************

      SUBROUTINE PrintBarkerkofE
      
      implicit none

      include 'cdparams.fh'
      include 'straightBS.fh'
      include 'cdnos.fh'
      include 'cdwell0.fh'
      include 'cdwell1.fh'
      include 'cdwell2.fh'
      include 'cdlabels.fh'
      include 'cdisprop.fh'
      include 'cdrange.fh'
c----------------------------------------------------------------------c
c *** RRKM 1 k(E) Thread *** START                                     c
c----------------------------------------------------------------------c
      include 'cdcontrl.fh'
c----------------------------------------------------------------------c
c *** RRKM 1 k(E) Thread *** END                                       c
c----------------------------------------------------------------------c



      integer iWell, iProd, imax1, Isize, nchan2, ndum, JJ, Jr
      integer ifrom, ito, i, ioWell, nSum,count
      real*8 Egrain1, Emax2, Viblo, RR, Afactor, Ebarr, energy
      real*8 ebarcm, rhdenom, rkofE, Egrain2, ebarcm2
      real*8 xNOScm, rkrev
      character name1*10, ctot, tsnames(mxWells,mxWells), tstr1*2
      character tstr2*2, tmpname*10


c     Set "parameters for transistion states and reactions" for barker
c     -compatible lines in the TSkey.DAT file

      RR = 1.0
      ndum = 1
      JJ = 2
      Jr = 0

c     Open Barker data file and read input data

      open(101, file='barker.inp', status = 'old')
      read(101,*) Egrain1, imax1, Emax2, Isize, Viblo
      close(101)

c     Open file for TS Key and write out labels for each well and product

c      ctot = 'A'
      open(102, file='TSkey.DAT')

      nchan2 = 1
      count=0

c     Egrain1, isize1, Imax, Emax2, IDUM
      write(102, 971) Egrain1, imax1, Isize, Emax2
 971  format(F3.0,I9,I9,F10.0,'        0')
c     no IDUM in barker.inp
      nSum=0
      do 972 iWell=1,nWells
         nSum=nSum+nProds(iWell)
 972  continue
      write(102,*),nWells,nSum

      do 10 iWell = 1, nWells

c     Note -- below is for old Barker format ... comment out to use

         tmpname = isName(iWell)
         call PadBlanks(tmpname, '_')
         
c***** OLD BARKER FORMAT TSKEY *****************************************
         write(102,901) iWell, tmpname, Ewell(iWell), 1.0, 1, 1, 1
     $        , 1, 57, 1, 0.0, 0.0, 0.0
 901     format(I2,T5,'''',A10,'''',TR5,F9.2,1x,F3.1,1x,I2,1x,I2,1x,I2
     $        ,1x,I2,1x,I3
     $        ,1x,I2, 1x, F3.1, 1x, F3.1, 1x, F3.1, 1x)
c***********************************************************************


c         write(102,901) iWell, isName(iWell)
c         write(102,911) Ewell(iWell), 1.0, 1, 1, 1
c     $        , 1, 57, 1, 0.0, 0.0, 0.0
c
c 901     format(I2,','A10)
c 911     format(F9.2,1x,F3.1,1x,I2,1x,I2,1x,I2,1x,I2,1x,I3
c     $        ,1x,I2, 1x, F3.1, 1x, F3.1, 1x, F3.1, 1x)

         nchan2 = nchan2 + 1
 10   continue

      do 15 iWell = 1, nWells
         do 17 iProd = 1, nProds(iWell)
            tmpname = pdName(iWell,iProd)
            call PadBlanks(tmpname, '_')
            write(102, 901) nchan2, tmpname
            nchan2 = nchan2 + 1
 17      continue
 15   continue

c***********************************************************************
c     BEGIN LOOP OVER PRODUCT CHANNELS
c***********************************************************************
      
      nchan2 = nWells + 1
      do 20 iWell = 1, nWells
         do 30 iProd = 1, nProds(iWell)

            ifrom = iWell
            ito = nchan2
         
c     Build "TS" filename

            name1 = 'TS'
            write(tstr1, '(I2)') iWell
            write(tstr2, '(I2)') ito
            name1(3:4) = tstr1
            name1(5:6) = 'to'
            name1(7:8) = tstr2
            
            call PadBlanks(name1, '_')
            
c            tsnames(iWell,iProd) = name1

c     Calculate Afactor with temperature factors ...

            Afactor = Aprod(iWell,iProd)*temp(1)**rNProd(iWell,iProd)
     $           *dexp(-alphaProd(iWell,iProd)*temp(1))
            Ebarr = EProd(iWell,iProd)

c     From Barker input file desc.:
c     RR = 2D external moment of inertia (should not be needed)
c     j, k, l (ndum, ndum, ndum) = ext sym, elec sym, optical isomer
c     numbers
c     JJ = 2 for read kofe files
c     Jr = 0 (neglect reverse; necessary where JJ = 2)
c     AA = afactor, EE = critical energy (barrier in CHEMDIS case)

c     Write line to TSkey.DAT file

            
c***** OLD BARKER FORMAT ***********************************************
            write(102, 902) ifrom, ito, name1, RR, ndum, ndum, ndum, JJ
     $           ,Jr, Afactor, Ebarr
 902        format(I2,TR6,I2,TR6,'''',A10,'''',TR6,F4.1,TR6,I2,TR6,I2
     $           ,TR6,I2,TR6
     $           ,I2,TR6,I2,TR6,E8.2,1X,F7.2)
            count=count+1
c***********************************************************************

c            write(102, 902) ifrom, ito, name1
c            write(102, 912) RR, ndum, ndum, ndum, JJ
c     $           ,Jr, Afactor, Ebarr


c 902        format(I2,','I2','A10)


 912        format(F4.1,TR6,I2,TR6,I2,TR6,I2,TR6
     $           ,I2,TR6,I2,TR6,E8.2,1X,F7.2)
            
c           ctot = CHAR(ICHAR(ctot) + 1)
            nchan2 = nchan2 + 1

c     open the .rke file, print header lines

            open(103,file = ''//name1//'.rke')
            write(103,903) Egrain1, imax1, Emax2, Isize
 903        format(//,1x,f10.1,1x,I10,1x,f10.1,1x,I10,1x,f10.1)
            write(103,*)

c     Begin first part of double array:  1 to imax1.  Energy is grained
c     in cm-1
            energy = 0.0d0
            ebarcm = Ebarr*349.755
            do 40 i = 1, imax1

             
c     check whether calling energy will be greater than that available
c     for the well

               if( nint((energy+ebarcm)/bsgs).gt.maxBSindex(iWell)) then
c     may want to change this to max possbile value, e.g. set denom to
c     largest avail as barker does (!)
                  rkofE = 0.0
                  
c                  write (*,*) 'Exceeded well DOS'
c                  write (*,*) iWell, iProd, energy, nint((energy+Ebarr)
c     $                 /bsgs)
                  goto 41
               endif
               
               rhdenom = xNOScm(energy+ebarcm, Egrain1,iWell)
               rkofE = Afactor*xNOScm(energy, Egrain1, iWell)/rhdenom

 41            continue

               write(103, 904) i, energy, RR, rkofE
 904           format(1x, I5, 3x, F9.2, 3x, F9.2, 3x, G10.4)
               energy = energy + Egrain1

 40         continue

c     Second part of double array
            
            Egrain2 = Emax2/(dfloat(Isize-(imax1+1)))
            energy = 0.0
            do 50 i = imax1+1, Isize
               if( nint((energy+ebarcm)/bsgs).gt.maxBSindex(iWell)) then
                  rkofE = 0.0
                  goto 51
               endif

              rhdenom = xNOScm(energy+ebarcm,Egrain2,iWell)
              rkofE = Afactor*xNOScm(energy, Egrain2, iWell)/rhdenom
 51           continue
              write(103,904) i, energy, RR, rkofE
              energy = energy + Egrain2
 50        continue

c     Close the current .rke file
           close(103)
 30     continue     !close product loop
 20   continue       !close isomer loop
      
c***********************************************************************
c     END LOOP OVER PRODUCT CHANNELS; BEGIN LOOP OVER ISOMs
c***********************************************************************

      do 60 iWell = 1, nWells
         do 70 ioWell = 1, nWells
            ifrom = iWell
            ito = ioWell
            
c     Determine whether a channel exists; if so, enter block
            if(Aisom(iWell,ioWell).gt.1.D-175) then
               
c     Build "TS" filename
               name1 = 'TS'
               write(tstr1, '(I2)') iWell
               write(tstr2, '(I2)') ito
               name1(3:4) = tstr1
               name1(5:6) = 'to'
               name1(7:8) = tstr2
            
               call PadBlanks(name1, '_')
c            tsnames(iWell,iProd) = name1

               Afactor = Aisom(iWell,ioWell)*temp(1)**rNIsom(iWell
     $              ,ioWell)*dexp(-alphaIsom(iWell,ioWell)*temp(1))
               Ebarr = EIsom(iWell,ioWell)

c     Write line to TSkey.DAT file

c***** OLD BARKER FORMAT ***********************************************
            write(102, 902) ifrom, ito, name1, RR, ndum, ndum, ndum, JJ
     $           ,Jr, Afactor, Ebarr
            count=count+1
c***********************************************************************

c            write(102, 902) ifrom, ito, name1
c            write(102, 912) RR, ndum, ndum, ndum, JJ
c     $           ,Jr, Afactor, Ebarr



c     Open the .rke file, print header lines
               open(103,file = ''//name1//'.rke')
               write(103,903) Egrain1, imax1, Emax2, Isize
               write(103,*)

c     Begin first part of double array:  1 to imax1
               energy = 0.0
               ebarcm = Ebarr*349.755
               do 80 i = 1, imax1
c     check whether calling energy will be greater than that available
c     for the well
                  
                  if( nint((energy+ebarcm)/bsgs).gt.maxBSindex(iWell))
     $                 then
c     may want to change this to max possbile value, e.g. set denom to
c     largest avail as barker does (!)
                     rkofE = 0.0
                     
c     write (*,*) 'Exceeded well DOS'
c     write (*,*) iWell, iProd, energy, nint((energy+Ebarr)
c     $                 /bsgs)
                     goto 81
                  endif

c----------------------------------------------------------------------c
c *** RRKM 1 k(E) Thread *** START                                     c
c----------------------------------------------------------------------c
                  if((RRKM).and.(i1RRKM.eq.iWell).and.(i2RRKM.eq.ioWell)
     $                 ) then
                     call kofERRKM(energy, ebarcm, Egrain1, iWell,
     $                    ioWell, rkofE, rkrev)
c     this section forces reverse k(E) to be calc'd by detailed balance
c     ... requires rotor moment and freqs to be specified for other
c     isomer (i2RRKM) as well as the first.  Comment out if desired.
c     Problem here is that density of states ratio used for conversion
c     probably is not right -- one is 3-freq, other may be explicit, not
c     clear that we can form this ratio  Possible units issue
c     
c     
c                  else if((RRKM).and.(i1RRKM.eq.ioWell).and.(i2RRKM.eq
c     $                    .iWell)) then
c                     ebarcm2 = 349.755*Eisom(i1RRKM, i2RRKM)
c                     call kofERRKM(energy, ebarcm2, Egrain1, i1RRKM,
c     $                    i2RRKM, rkofE, rkrev)
c                     rkofE = rkrev
c----------------------------------------------------------------------c
c *** RRKM 1 k(E) Thread *** END                                       c
c----------------------------------------------------------------------c
                  else
                     rhdenom = xNOScm(energy+ebarcm, Egrain1, iWell)
                     rkofE = Afactor*xNOScm(energy, Egrain1, iWell)
     $                    /rhdenom
                  endif
                  
 81               continue
                  write(103, 904) i, energy, RR, rkofE
                  energy = energy + Egrain1
 80            continue
               
c     Second part of double array
               
               Egrain2 = Emax2/(dfloat(Isize-(imax1+1)))
               energy = 0.0
               do 90 i = imax1+1, Isize
                  if( nint((energy+ebarcm)/bsgs).gt.maxBSindex(iWell))
     $                 then
                     rkofE = 0.0
                     goto 91
                  endif

c----------------------------------------------------------------------c
c *** RRKM 1 k(E) Thread *** START                                     c
c----------------------------------------------------------------------c
                  if((RRKM).and.(i1RRKM.eq.iWell).and.(i2RRKM.eq.ioWell)
     $                 ) then
                     call kofERRKM(energy, ebarcm, Egrain2, iWell,
     $                    ioWell, rkofE, rkrev)
c     this section forces reverse k(E) to be calc'd by detailed balance
c     ... requires rotor moment and freqs to be specified for other
c     isomer (i2RRKM) as well as the first
c$$$                  else if((RRKM).and.(i1RRKM.eq.ioWell).and.(i2RRKM.eq
c$$$     $                    .iWell)) then
c$$$                     ebarcm2 = 349.755*Eisom(i1RRKM, i2RRKM)
c$$$                     call kofERRKM(energy, ebarcm2, Egrain1, i1RRKM,
c$$$     $                    i2RRKM, rkofE, rkrev)
c$$$                     rkofE = rkrev
c
c----------------------------------------------------------------------c
c *** RRKM 1 k(E) Thread *** END                                       c
c----------------------------------------------------------------------c
                  else
                     rhdenom = xNOScm(energy+ebarcm, Egrain2, iWell)
                     rkofE = Afactor*xNOScm(energy, Egrain2, iWell)
     $                    /rhdenom
                  endif

 91               continue
                  write(103,904) i, energy, RR, rkofE
                  energy = energy + Egrain2
 90            continue

c     Close file
               close(103)
            endif
            
 70      continue
 60   continue
      write(102,*) 'Count'
      write(102,*) count
      write(102,*) 'END'
      write(102,*)

      return
      end

c***********************************************************************
c***********************************************************************
c     SUBROUTINE PadBlanks
c
c     This subroutine pads all the trailing string1 'blanks' with pchar,
c     assumes length of 10 (stupid FORTRAN!)
c     
c***********************************************************************
c***********************************************************************
      
      SUBROUTINE PadBlanks(string1, pchar)

      character pchar*1
      character string1*10
      integer ispace, nbslength, ic
      
      nbslength = 10
      ispace = INDEX(string1, ' ')
      do 10 ic = 1, nbslength
         if(string1(ic:ic).eq.' ') then            
            string1(ic:ic) = pchar
         endif
 10   continue
      
      return
      end

c$Id$
c$Author$
c$Date$
c$Revision$
c$Log$
cRevision 1.1  2007-02-20 23:10:23  sandeeps
cInitial revision
c
cRevision 1.4  2001/12/02 18:17:35  dmmatheu
cBefore changes to correspond to offset problem.  Printing of sum of
cstates should be offset by 1 and not by Egrain1.  Should just give
cnumber of states at this energy!
c
cRevision 1.3  2001/12/01 01:30:26  dmmatheu
cbefore changes to correct apparent density of states printout error:
chave been (!!!) printing out NUMBER of states in width, not density
c... does MW even use this or just Sum of States ???
c
cRevision 1.2  2001/07/07 18:12:31  dmmatheu
cbefore adding amamo changes ...
c
cRevision 1.1  2001/07/02 21:00:44  dmmatheu
cInitial revision
c


