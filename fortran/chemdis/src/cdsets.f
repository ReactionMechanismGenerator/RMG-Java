c******6/10/98  corrected output formatting error in statement #30   amd
c********12/12/97  modivied as suggested by Jeff Grenda to address problem
c   of code blowing up when exit channel lower than entrance
c    (for only 1 well and 1 exit channel)
c SETCALC*******************
      subroutine setCalc(lecho,option)
c    finds how deep the wells are on the chart; needed for high p limit 
c    also determines actual well depths in kcal

      implicit none
      include 'cdparams.fh'
      include 'cdwell0.fh'
      include 'cdwell1.fh'
      include 'cdwell2.fh'
      include 'cdlabels.fh'
      include 'cdcontrl.fh'

c    local variables
      integer iwell,iowell,iprod,icount
      integer idepth,lecho,ioccur
      real*8 Epmax
      character option*8
cjmg fix - already declared
c     logical isomers
c   end declarations

c    check isomerization matrix
      do 50 iwell = 1,nwells
         do 50 iowell = iwell,nwells
            if (.not.(((Aisom(iwell,iowell).gt.1.d-175).and.
     2         (Aisom(iowell,iwell).gt.1.d-175)).or.
     3         ((Aisom(iwell,iowell).lt.1.d-175).and.(Aisom(iowell,
     4         iwell).lt.1.d-175)))) then
                  write(lecho,40) ISname(iwell),ISname(iowell)
40                format(/,' ERROR:  no reciprocity between isomers ',
     2               a,' and ',a)
                  stop
            endif            
50    continue
 
c    get minimum and maximum E exit channel for each well
      do 100 iwell = 1,nwells
         Exmin(iwell) =  1.d175
         Exmax(iwell) = -1.d175
         if (nprods(iwell).gt.0) then
            do 70 iprod = 1, nprods(iwell)
               Exmax(iwell) = dmax1(Exmax(iwell),Eprod(iwell,iprod))
               Exmin(iwell) = dmin1(Exmin(iwell),Eprod(iwell,iprod))
70          continue
         endif
         if (isomers) then
            do 80 iowell = 1,nwells
               if (Aisom(iwell,iowell).gt.1.d-175) then
                  Exmax(iwell) = dmax1(Exmax(iwell),
     2               Eisom(iwell,iowell))
                  Exmin(iwell) = dmin1(Exmin(iwell),
     2               Eisom(iwell,iowell))
               endif
80          continue
         endif
100   continue

c    input well is entrance and has depth 1; find those connected to depth 1
c    and assign them depth 2; repeat etc.; first initialize

      do 170 iwell = 1,nwells
         idpwell(iwell) = nwells+1
170   continue 
        
      idpwell(inpwell) = 1
      
      do 200 idepth = 2,nwells
         do 190 iwell=1,nwells
            do 180 iowell = 1,nwells
               if ((Aisom(iowell,iwell).gt.1.d-175).and.(idpwell(iowell)
     2            .eq.idepth-1).and.(idepth.lt.idpwell(iwell)))
     3            idpwell(iwell) = idepth
180         continue
190      continue
200   continue

c    check for unconnected isomers
      do 250 iwell = 1,nwells
         if (idpwell(iwell).gt.nwells) then
            write(lecho,240) ISname(iwell)
240         format(/,1x,'ERROR: isomer ',a,' not connected to input ',
     2         'channel.')
            stop
         endif
250   continue            

c dmm 20000627
c     If ASA has been specified, and this is our 'first time' through
c     setCalc, then break out here for ASAScreen
      if(ASActrl) then
         write(*, 252)
cjmg fix
c252     format(1X, 'INFO::cdsets.f: setCalc check complete on whole',1X
c    $   ,'network.  ASA specified, '/, 16X' leaving setCalc.')
 252     format(1X, 'INFO::cdsets.f: setCalc check complete on whole',1X
     $   ,'network.  ASA specified, '/, 16X,' leaving setCalc.')
         return
      endif
c dmm 20000627

c    now calculate well depths in kcal and EbarHP - the total barrier
c    to isomerization in the hp limit (one that takes the shortest
c    path)

      if (option.eq.'Chemact') then
         Ewell(inpwell) = -Eprod(inpwell,inpchan)
      else
         Ewell(inpwell) = 0.d0
      endif

      EbarHP(inpwell) = Ewell(inpwell)
      
      do 300 idepth = 2,nwells
         do 290 iwell = 1,nwells
	 
            if (idpwell(iwell).eq.idepth) then
               do 280 iowell = 1,nwells
	       
                  ioccur = 0
                  if ((idpwell(iowell).eq.idepth-1).and.
     2               (Aisom(iwell,iowell).gt.1.d-175)) then
                     ioccur = ioccur+1	     
     
                     if (ioccur.eq.1) then
                        Ewell(iwell) = Ewell(iowell) + 
     2                     Eisom(iowell,iwell) - Eisom(iwell,iowell)
                        EbarHP(iwell) = dmax1(EbarHP(iowell), 
     2                     Ewell(iowell) + Eisom(iowell,iwell))
     
                     else if (ioccur.gt.1) then
	                EbarHP(iwell) = dmin1(EbarHP(iwell), 
     2                     dmax1(EbarHP(iowell), 
     3                     Ewell(iowell) + Eisom(iowell,iwell)))                   
                     endif

                  endif
280            continue
            endif
290      continue
300   continue

c    now actually find most isolated output channel in hp limit
      Epmax = 0.d0
c*********added 12/12/97  as suggested by Jeff Grenda
cjmg                            
      Epmax = -1.d10            
      do 400 iwell = 1,nwells
         if (EbarHP(iwell).gt.Epmax) then
            Epmax = EbarHP(iwell)
            iwcmax = iwell
            ipcmax = 0
         endif
         
         do 400 iprod = 1,Nprods(iwell)
            if ((EbarHP(iwell) + Eprod(iwell,iprod)).gt.Epmax) then
               Epmax = EbarHP(iwell) + Ewell(iwell) + Eprod(iwell,iprod)
               iwcmax = iwell
               ipcmax = iprod
            endif              
400   continue         

c    now find EbarLP - the barrier to isomerization in the low pressure
c    limit - first initialize
      do 500 iwell = 1,nwells
         EbarLP(iwell) = EbarHP(iwell)
500   continue

c    we have to repeat this nwell times in order that we have everything
c    updated properly - this is the same algorithm as EbarHP, except
c    we don't care what the depth of the connecting well is

      do 600 icount = 1,nwells
         do 600 iwell = 1,nwells
            do 600 iowell = 1,nwells
               if (Aisom(iwell,iowell).gt.1.d-175)
     2            EbarLP(iwell) = dmin1(EbarLP(iwell), 
     3            dmax1(EbarHP(iowell), 
     4            Ewell(iowell) + Eisom(iowell,iwell))) 
600    continue

c    finally check consistency of E's
      do 660 iwell = 1,nwells
         do 660 iowell = 1,nwells
            if ((Aisom(iwell,iowell).gt.1.d-175).and.
     2         abs( (Eisom(iwell,iowell)-Eisom(iowell,iwell)) -
     3         (Ewell(iowell)-Ewell(iwell)) ).ge.1.d-2) then
               write(lecho,655) ISname(iwell),ISname(iowell)
655            format(/,1x,'ERROR: E''s on ',a,' and ',a,' are ',
     2            'not consistent.')
               stop
            endif
660   continue

      return
      end
            
c ECHO2******************
      subroutine echo2(lecho,option)

      implicit none
      include 'cdparams.fh'
      include 'cdwell0.fh'
      include 'cdwell1.fh'
      include 'cdwell2.fh'
      include 'cdisprop.fh'
      include 'cdlabels.fh'
                             
c    local variables
      integer iwell
      integer lecho
      character option*8
c    end declarations

      if (option.eq.'Chemact') then
         write(lecho,10)
10       format(//,1x,'**** CHEMACT PARAMETERS ****')
         write(lecho,20)
20       format(/,1x,'REACTANT CHANNEL ',26x,
     2      '   A   ',4x,'   n',4x,'  alpha    ',1x,' E (kcal)')     
         write(lecho,30) inpwell,inpchan,PDname(inpwell,inpchan),
     2      ISname(inpwell),Ain,rNin,alphaIn,Ein     
c******6/10/98  corrected output formatting error in statement below:
c 30       format(1x,i2':-',i2,2x,a,' => ',a,1pe10.4,2x,
c     2      0pf6.3,2x,1pe9.3,3x,f8.3)
cjmg fix
c30       format(1x,i2':-',i2,2x,a,' => ',a,1pe10.4,2x,
c     2      0pf6.3,2x,1pe9.3,3x,0pf8.3)
30       format(1x,i2,':-',i2,2x,a,' => ',a,1pe10.4,2x,
     2      0pf6.3,2x,1pe9.3,3x,0pf8.3)
     
      else
         write(lecho,35) inpwell
35       format(//,1x,'**** DISSOC (',i2,') PARAMETERS ****')
         write(lecho,40)
40       format(/,1x,'REACTANT CHANNEL ')
         write(lecho,45) inpwell,inpchan,ISname(inpwell),ISname(inpwell)
cjmg fix
c45       format(1x,i2': ',i2,2x,a,' => *',a)
45       format(1x,i2,': ',i2,2x,a,' => *',a)
      endif

      write(lecho,50)
50    format(/,1x,'WELL INFO ',6x,'  kcal ',2x,
     2   '# depth',2x,' sigma ',1x,' e/k  ',3x,'Ebar HP',
     3   3x,'Ebar LP')
     
      do 70 iwell = 1,nwells
         write(lecho,60) iwell,ISname(iwell),Ewell(iwell),
     2         idpwell(iwell),sig(iwell),ek(iwell),EbarHP(iwell),
     3         EbarLP(iwell)
60       format(1x,i2,2x,a,2x,f7.2,5x,i2,4x,f5.2,3x,f6.2,4x,f6.2,
     2      4x,f6.2)
70    continue

      write(lecho,110)
110   format(/,1x,'Well depth (kcal) is determined from input',
     2   ' barrier heights.')  
      write(lecho,120)
120   format(1x,'Chart depth is determined by isomerization',
     2   ' linkage.')     

      return
      end

c $Id$
c $Author$
c $Date$
c $Log$
c Revision 1.1  2007-02-20 23:10:23  sandeeps
c Initial revision
c
c Revision 1.1  2003/04/23 19:11:12  dmmatheu
c Initial revision
c
