c***********************************************************************
c***********************************************************************
c     
c     Subroutine PrintSteadyInput
c
c     Prints output for calculation of steady state concentration ratios
c     for purpose of testing agreement with Barker MEQ code
c     
c     dmm 20010201
c***********************************************************************
c***********************************************************************

      SUBROUTINE PrintSteadyInput(option)

cjmg fix
c     character option*8
cjmg fix
      implicit none

      include 'cdparams.fh'
      include 'cdwell0.fh'
      include 'cdlabels.fh'
      include 'cdrange.fh'
      include 'cdrates.fh'   
      include 'cdcontrl.fh'

cjmg fix
c      implicit none
      integer ip, it
      integer numchan
      integer iWell, iProd
      character tmpname*10
cjmg fix
      character option*8

      ip = 1
      it = 1
c     Warnings for multiple T and P
      if((nPres.gt.1).or.(nTemps.gt.1)) then
         write(*,*)
     $        '<<INFO>>:  Multiple T and/or P for STEADY:  only 1st T', 
     $         'and P will be used'

      endif

      if(option.eq.'chemact') then
         open(100, file='saCHEMACT.inp', status = 'unknown')      
         write(100,510) 
         write(100,520) PDname(inpwell,inpchan)
      else if(option.eq.'dissoc') then
         tmpname = PDname(inpwell,inpchan)
         call PadBlanks(tmpname, '_')
         open(100, file='ss'//tmpname//'.inp')
         write(100,540)
         write(100,520) PDname(inpwell,inpchan)
      else
         write(*,*) '<<ERROR>> PrintSteady called without option'
         pause
      endif

      call ChanCount(numchan)
      write(100,530) numchan
      write(100,530) nWells
c     First write out rate constants to other wells
      do 10 iWell = 1, nWells
         write(100,550) PDname(iWell,0), dlog10(RK(iWell,0, it, ip))
 10   continue

c     Next write out rate constants to product channels

      do 20 iWell = 1, nWells
         do 30 iProd = 1, nProds(iWell)
            
            write(100,550) PDname(iWell,iProd), dlog10(RK(iWell,iProd,it
     $           ,ip))
 30      continue
 20   continue



 510  format(1x, 'CHEMACT ENTRANT')
 520  format(1x,'''', A20, '''')
 530  format(1x,I4)
 540  format('DISSOC INPWELL')
 550  format('''',A20,'''', 5x, F10.5)
      end

c***********************************************************************
c***********************************************************************


c***********************************************************************
c***********************************************************************
c  ChanCount
c
c     Counts total number of channels for which rate constants have been
c     predicted
c     
c     dmm 20010702
c
c***********************************************************************
c***********************************************************************

      SUBROUTINE ChanCount(numchan)

      include 'cdparams.fh'
      include 'cdwell0.fh'

      integer i,j
      integer numchan = 0


      do 10 i = 1, nWells
         do 20 j = 0, nProds(i)
            numchan = numchan + 1
 20      continue
 10   continue


      return
      end
c***********************************************************************
c***********************************************************************


c$Id$
c$Author$
c$Date$
c$Revision$
c$Log$
cRevision 1.1  2007-02-20 23:10:23  sandeeps
cInitial revision
c
cRevision 1.2  2003/04/23 18:59:55  dmmatheu
cbefore JMG revision merge
c
cRevision 1.1  2001/07/02 21:00:52  dmmatheu
cInitial revision
c
