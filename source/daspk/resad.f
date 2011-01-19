C     - DAEPACK Version 1.1 - Copyright (C) Numerica Technology, LLC
C     - DERIVATIVE COMPUTATION -
C     -   GENERATED: Wed Jan 31 02:02:16 2007
      subroutine resad(t,y,yprime,cj,del,rpar,ipar,senpar,zzzderiv,
     $   zzzne,zzzirn,zzzjcn,zzziw)
!!independent { y(1:nstate) }
!!dependent   { del(1:nstate) }
!DEC$ IF DEFINED(DLLEXP)
!DEC$ ATTRIBUTES DLLEXPORT :: resad
!DEC$ END IF
      implicit none
      integer neqmax
      parameter(neqmax=500*25000)
      integer spmax
      parameter(spmax=500)
      integer rmax
      parameter(rmax=25000)
      integer tbrmax
      parameter(tbrmax=50)
      integer troemax
      parameter(troemax=50)
      integer lindemax
      parameter(lindemax=50)
      integer i
      integer j
      integer k
      double precision t
      double precision temperature
      integer nstate
      double precision pressure
      integer troereactionsize
      integer lindereactionsize
      integer ires
      double precision cj
      integer ijac
      integer senpar
      integer thirdbodyreactionsize
      integer nparam
      integer neq
      integer reactionsize
      double precision thirdbodyreactionratearray
      integer thirdbodyreactionarray
      double precision y(neq)
      double precision reactionratearray
      integer ipar(*)
      double precision del(neq)
      double precision yprime(neq)
      double precision rpar(reactionsize+thirdbodyreactionsize+
     $     troereactionsize+lindereactionsize+nstate-1)
      integer reactionarray
      integer troereactionarray
      integer lindereactionarray
      double precision troereactionratearray
      double precision lindereactionratearray
      common /size/ nstate,neq,nparam,reactionsize,troereactionsize,
     $   thirdbodyreactionsize,lindereactionsize
      common /reac/ reactionratearray(5*rmax),
     $   thirdbodyreactionratearray(16*tbrmax),troereactionratearray(
     $   21*troemax),temperature,pressure,reactionarray(9*rmax),
     $   thirdbodyreactionarray(20*tbrmax),troereactionarray(21*
     $   troemax),lindereactionratearray(17*lindemax),
     $   lindereactionarray(20*lindemax)
C        zzzderiv - partial derivative array stored in sparse matrix
C                   format, pattern contained in zzzirn and zzzjcn.
C        zzzne    - number of entries in zzzderiv.
C        zzzirn   - row numbers for entries in zzzderiv.
C        zzzjcn   - column numbers for entries in zzzderiv.
      double precision zzzderiv(*)
      integer zzzne,zzzirn(zzzne),zzzjcn(zzzne)
C        zzziw - integer workspace array.
      integer zzziw(*)
C        zzzimodel - unique key for this model.
      integer zzzimodel
C     
C     Define elementary variables and adjoints
      double precision zzzv1,zzzvbar1
      double precision zzzv2,zzzvbar2
      double precision zzzv3,zzzvbar3
C     
C     Active variable offsets and counters
      integer deloft
      integer yoft
      integer reactionratearrayoft
      integer thirdbodyreactionratearrayoft
      integer troereactionratearrayoft
      integer lindereactionratearrayoft
      integer temperatureoft
      integer pressureoft
      integer zzzindx
C     Declare number of independent and dependent variables
      integer zzzn
      integer zzzm
C     Loop control variables
      integer izzz1
C     
C     Variable offsets for independent, dependent, and 
C     non-common block local active variables.
      deloft=0
      yoft=deloft+neq
C     Variable offsets for COMMON BLOCK size
C     Variable offsets for COMMON BLOCK reac
      reactionratearrayoft=0
      thirdbodyreactionratearrayoft=reactionratearrayoft+(5*rmax)
      troereactionratearrayoft=thirdbodyreactionratearrayoft+(16*tbrmax)
      lindereactionratearrayoft=troereactionratearrayoft+(21*troemax)
      temperatureoft=lindereactionratearrayoft+(17*lindemax)
      pressureoft=temperatureoft+1
C     Number of independent and dependent variables.
      zzzn=nstate
      zzzm=nstate
C     Create unique identifier for this model.
      call GETMDLID(zzzimodel,'resad')
C      call SETTABLESIZE()
C      call SETSEGMENTSIZE()
C      call SETNUMSEGMENTS()
C     Create SVM for SUBROUTINE/FUNCTION res
      call CREATESV(zzzimodel,1,zzzn)
C     Create SVM for SUBROUTINE/FUNCTION getflux
      call CREATESV(zzzimodel,2,zzzn)
C     Create SVM for COMMON BLOCK reac
      call CREATESV(zzzimodel,3,zzzn)
C     Create SVM for COMMON BLOCK size
      call CREATESV(zzzimodel,4,zzzn)
C     Construct mapping between independent variable offsets and indices
      zzzindx=0
      do izzz1=1,nstate
         zzzindx=zzzindx+1
         zzziw(zzzindx)=yoft+izzz1
      end do
C     Initialize independent variable gradients to Cartesian basis vectors.
      call DSETIVCV(1,zzzn,zzziw)
C     
C     Save SVMs of arguments
      zzziw(1)=1
      zzziw(2)=1
C     Save current offsets for use within embedded subroutine
      zzziw(3)=yoft
      zzziw(4)=deloft
C     Call modifed embedded subroutine
      call getfluxad(y,del,rpar,zzziw)
      do i=1,nstate
         zzzv3=del(i)-yprime(i)
         zzzvbar1=1.0d0
         del(i)=zzzv3
         call DSVM1(deloft+i,1,zzzvbar1,deloft+i,1)
C     
      end do
C     Construct derivative matrix and sparsity pattern before returning.
C     Construct mapping between dependent variable offsets and indices
      zzzindx=0
      do izzz1=1,nstate
         zzzindx=zzzindx+1
         zzziw(zzzindx)=deloft+izzz1
      end do
      call DCPRVS(1,zzzm,zzziw,zzzderiv,zzzne,zzzirn,zzzjcn)
C     UNCOMMENT BELOW TO SAVE STATISTICS TO FILE
C           call writesvmstats(zzzimodel,'svm.stats')
C     Delete SVM for SUBROUTINES/FUNCTIONS
      call DELETESV(zzzimodel,1)
      call DELETESV(zzzimodel,2)
C     Delete SVM for COMMON BLOCK reac
      call DELETESV(zzzimodel,3)
C     Delete SVM for COMMON BLOCK size
      call DELETESV(zzzimodel,4)
C     Done
      end
      subroutine getfluxad(tempy,del,rpar,zzziw)
      implicit none
      integer neqmax
      parameter(neqmax=500*25000)
      integer spmax
      parameter(spmax=500)
      integer rmax
      parameter(rmax=25000)
      integer tbrmax
      parameter(tbrmax=50)
      integer troemax
      parameter(troemax=50)
      integer lindemax
      parameter(lindemax=50)
      integer numcollider
      double precision logpr
      double precision lowrate
      double precision frate
      double precision c
      double precision d
      double precision f
      integer i
      integer j
      double precision m
      double precision n
      double precision temperature
      double precision flux
      integer nstate
      double precision pressure
      double precision pr
      integer troereactionsize
      integer lindereactionsize
      double precision dg
      integer thirdbodyreactionsize
      integer nparam
      double precision alpha
      double precision logf
      double precision sumyprime
      integer direction
      integer neq
      double precision tstar
      double precision rate
      double precision rrate
      double precision keq
      integer rnum
      double precision fcent
      integer pnum
      double precision t3star
      double precision t2star
      integer reactionsize
      double precision sumy
      double precision inertefficiency
      double precision logfcent
      double precision inside
      double precision thirdbodyreactionratearray
      integer thirdbodyreactionarray
      double precision y(neq)
      double precision reactionratearray
      double precision tempy(*)
      double precision thermo(spmax)
      double precision del(*)
      double precision rpar(*)
      integer reactionarray
      integer troereactionarray
      integer lindereactionarray
      double precision troereactionratearray
      double precision lindereactionratearray
      common /size/ nstate,neq,nparam,reactionsize,troereactionsize,
     $   thirdbodyreactionsize,lindereactionsize
      common /reac/ reactionratearray(5*rmax),
     $   thirdbodyreactionratearray(16*tbrmax),troereactionratearray(
     $   21*troemax),temperature,pressure,reactionarray(9*rmax),
     $   thirdbodyreactionarray(20*tbrmax),troereactionarray(21*
     $   troemax),lindereactionratearray(17*lindemax),
     $   lindereactionarray(20*lindemax)
C     Additional arguments for partial derivative construction
C        zzziw - integer workspace
      integer zzziw(*)
C     Define elementary variables and adjoints
      double precision zzzv1,zzzvbar1
      double precision zzzv2,zzzvbar2
      double precision zzzv3,zzzvbar3
      double precision zzzv4,zzzvbar4
      double precision zzzv5,zzzvbar5
      double precision zzzv6,zzzvbar6
      double precision zzzv7,zzzvbar7
      double precision zzzv8,zzzvbar8
      double precision zzzv9,zzzvbar9
      double precision zzzv10,zzzvbar10
      double precision zzzv11,zzzvbar11
      double precision zzzv12,zzzvbar12
      double precision zzzv13,zzzvbar13
      double precision zzzv14,zzzvbar14
      double precision zzzv15,zzzvbar15
      double precision zzzv16,zzzvbar16
      double precision zzzv17,zzzvbar17
      double precision zzzv18,zzzvbar18
      double precision zzzv19,zzzvbar19
      double precision zzzv20,zzzvbar20
      double precision zzzv21,zzzvbar21
      double precision zzzv22,zzzvbar22
      double precision zzzv23,zzzvbar23
C     
C     Active variable offsets and counters
      integer tempyoft
      integer deloft
      integer reactionratearrayoft
      integer thirdbodyreactionratearrayoft
      integer troereactionratearrayoft
      integer lindereactionratearrayoft
      integer temperatureoft
      integer pressureoft
      integer yoft
      integer frateoft
      integer keqoft
      integer rrateoft
      integer inertefficiencyoft
      integer moft
      integer alphaoft
      integer tstaroft
      integer t2staroft
      integer t3staroft
      integer lowrateoft
      integer fcentoft
      integer proft
      integer logfcentoft
      integer logproft
      integer noft
      integer coft
      integer insideoft
      integer logfoft
      integer foft
      integer sumyprimeoft
      integer zzzindx
      integer zzzsvm1
      integer zzzsvm2
      zzzsvm1=zzziw(1)
      zzzsvm2=zzziw(2)
C     
C     Compute offsets for arguments and local non-common 
C     block active variables
      tempyoft=zzziw(3)
      deloft=zzziw(4)
C     
      yoft=0
      frateoft=yoft+neq
      keqoft=frateoft+1
      rrateoft=keqoft+1
      inertefficiencyoft=rrateoft+1
      moft=inertefficiencyoft+1
      alphaoft=moft+1
      tstaroft=alphaoft+1
      t2staroft=tstaroft+1
      t3staroft=t2staroft+1
      lowrateoft=t3staroft+1
      fcentoft=lowrateoft+1
      proft=fcentoft+1
      logfcentoft=proft+1
      logproft=logfcentoft+1
      noft=logproft+1
      coft=noft+1
      insideoft=coft+1
      logfoft=insideoft+1
      foft=logfoft+1
      sumyprimeoft=foft+1
C     
C     Variable offsets for COMMON BLOCK size
C     Variable offsets for COMMON BLOCK reac
      reactionratearrayoft=0
      thirdbodyreactionratearrayoft=reactionratearrayoft+(5*rmax)
      troereactionratearrayoft=thirdbodyreactionratearrayoft+(16*tbrmax)
      lindereactionratearrayoft=troereactionratearrayoft+(21*troemax)
      temperatureoft=lindereactionratearrayoft+(17*lindemax)
      pressureoft=temperatureoft+1
C     
      do i=1,neq
         del(i)=0
         call SVZ(deloft+i,zzzsvm2)
C     
      end do
      do i=1,nstate-1
         zzzv3=tempy(i)/tempy(nstate)
         zzzvbar2=-zzzv3/tempy(nstate)
         zzzvbar1=1.0d0/tempy(nstate)
         y(i)=zzzv3
         call DSVM2(yoft+i,2,zzzvbar2,tempyoft+nstate,zzzsvm1,
     $      zzzvbar1,tempyoft+i,zzzsvm1)
C     
      end do
      zzzvbar1=1.0d0
      y(nstate)=tempy(nstate)
      call DSVM1(yoft+nstate,2,1.0d0,tempyoft+nstate,zzzsvm1)
C     
      do i=1,nstate-1
         thermo(i)=rpar(reactionsize+thirdbodyreactionsize+
     $      troereactionsize+lindereactionsize+i)
      end do
      do i=0,reactionsize-1
         frate=rpar(i+1)
         call SVZ(frateoft+1,2)
C     
         rnum=reactionarray(9*i+1)
         pnum=reactionarray(9*i+2)
         dg=0
         do j=1,rnum
            dg=dg-thermo(reactionarray(9*i+2+j))
         end do
         do j=1,pnum
            dg=dg+thermo(reactionarray(9*i+5+j))
         end do
         zzzv2=-dg
         zzzv4=zzzv2*4184
         zzzv6=zzzv4/8.314
         zzzv8=zzzv6/temperature
         zzzv9=exp(zzzv8)
         zzzv11=82.053*temperature
         zzzv14=rnum-pnum
         zzzv15=zzzv11**zzzv14
         zzzv16=zzzv9*zzzv15
         zzzvbar15=zzzv9
         zzzvbar11=zzzv14*zzzv11**(zzzv14-1.0d0)*zzzvbar15
         zzzvbar9=zzzv15
         zzzvbar8=zzzv9*zzzvbar9
         zzzvbar7=-zzzv8/temperature*zzzvbar8+82.053*zzzvbar11
         keq=zzzv16
         call DSVM1(keqoft+1,2,zzzvbar7,temperatureoft+1,3)
C     
         if(reactionarray(9*i+9)==1) then
            zzzv3=rpar(i+1)/keq
            zzzvbar2=-zzzv3/keq
            rrate=zzzv3
            call DSVM1(rrateoft+1,2,zzzvbar2,keqoft+1,2)
C     
         else
            rrate=0
            call SVZ(rrateoft+1,2)
C     
         end if
         do j=1,rnum
            zzzv3=frate*y(reactionarray(9*i+2+j))
            zzzvbar2=frate
            zzzvbar1=y(reactionarray(9*i+2+j))
            frate=zzzv3
            call DSVM2(frateoft+1,2,zzzvbar2,yoft+reactionarray(9*i+2+
     $         j),2,zzzvbar1,frateoft+1,2)
C     
         end do
         do j=1,pnum
            zzzv3=rrate*y(reactionarray(9*i+5+j))
            zzzvbar2=rrate
            zzzvbar1=y(reactionarray(9*i+5+j))
            rrate=zzzv3
            call DSVM2(rrateoft+1,2,zzzvbar2,yoft+reactionarray(9*i+5+
     $         j),2,zzzvbar1,rrateoft+1,2)
C     
         end do
         do j=1,rnum
            zzzv3=del(reactionarray(9*i+2+j))-frate
            zzzv5=zzzv3+rrate
            zzzvbar4=1.0d0
            zzzvbar2=-1.0d0
            zzzvbar1=1.0d0
            del(reactionarray(9*i+2+j))=zzzv5
            call DSVM3(deloft+reactionarray(9*i+2+j),zzzsvm2,zzzvbar4,
     $         rrateoft+1,2,zzzvbar2,frateoft+1,2,zzzvbar1,deloft+
     $         reactionarray(9*i+2+j),zzzsvm2)
C     
         end do
         do j=1,pnum
            zzzv3=del(reactionarray(9*i+5+j))+frate
            zzzv5=zzzv3-rrate
            zzzvbar4=-1.0d0
            zzzvbar2=1.0d0
            zzzvbar1=1.0d0
            del(reactionarray(9*i+5+j))=zzzv5
            call DSVM3(deloft+reactionarray(9*i+5+j),zzzsvm2,zzzvbar4,
     $         rrateoft+1,2,zzzvbar2,frateoft+1,2,zzzvbar1,deloft+
     $         reactionarray(9*i+5+j),zzzsvm2)
C     
         end do
      end do
      do i=0,thirdbodyreactionsize-1
         frate=rpar(reactionsize+i+1)
         call SVZ(frateoft+1,2)
C     
         rnum=thirdbodyreactionarray(i*20+1)
         pnum=thirdbodyreactionarray(i*20+2)
         dg=0
         do j=1,rnum
            dg=dg-thermo(thirdbodyreactionarray(20*i+2+j))
         end do
         do j=1,pnum
            dg=dg+thermo(thirdbodyreactionarray(20*i+5+j))
         end do
         zzzv2=-dg
         zzzv4=zzzv2*4184
         zzzv6=zzzv4/8.314
         zzzv8=zzzv6/temperature
         zzzv9=exp(zzzv8)
         zzzv11=82.053*temperature
         zzzv14=rnum-pnum
         zzzv15=zzzv11**zzzv14
         zzzv16=zzzv9*zzzv15
         zzzvbar15=zzzv9
         zzzvbar11=zzzv14*zzzv11**(zzzv14-1.0d0)*zzzvbar15
         zzzvbar9=zzzv15
         zzzvbar8=zzzv9*zzzvbar9
         zzzvbar7=-zzzv8/temperature*zzzvbar8+82.053*zzzvbar11
         keq=zzzv16
         call DSVM1(keqoft+1,2,zzzvbar7,temperatureoft+1,3)
C     
         zzzv3=pressure*1e-6
         zzzv5=zzzv3/8.314
         zzzv7=zzzv5/temperature
         zzzvbar6=-zzzv7/temperature
         zzzvbar5=1.0d0/temperature
         zzzvbar3=1.0d0/8.314*zzzvbar5
         zzzvbar1=1e-6*zzzvbar3
         inertefficiency=zzzv7
         call DSVM2(inertefficiencyoft+1,2,zzzvbar6,temperatureoft+1,3
     $      ,zzzvbar1,pressureoft+1,3)
C     
         numcollider=thirdbodyreactionarray(i*20+10)
         do j=1,numcollider
            zzzv5=thirdbodyreactionratearray(16*i+6+j)-1
            zzzv6=y(thirdbodyreactionarray(i*20+10+j))*zzzv5
            zzzv7=inertefficiency+zzzv6
            zzzvbar5=y(thirdbodyreactionarray(i*20+10+j))
            zzzvbar3=zzzvbar5
            zzzvbar2=zzzv5
            zzzvbar1=1.0d0
            inertefficiency=zzzv7
            call DSVM3(inertefficiencyoft+1,2,zzzvbar3,
     $         thirdbodyreactionratearrayoft+16*i+6+j,3,zzzvbar2,yoft+
     $         thirdbodyreactionarray(i*20+10+j),2,zzzvbar1,
     $         inertefficiencyoft+1,2)
C     
         end do
         zzzv3=frate*inertefficiency
         zzzvbar2=frate
         zzzvbar1=inertefficiency
         frate=zzzv3
         call DSVM2(frateoft+1,2,zzzvbar2,inertefficiencyoft+1,2,
     $      zzzvbar1,frateoft+1,2)
C     
         if(thirdbodyreactionarray(20*i+9)==1) then
            zzzv3=frate/keq
            zzzvbar2=-zzzv3/keq
            zzzvbar1=1.0d0/keq
            rrate=zzzv3
            call DSVM2(rrateoft+1,2,zzzvbar2,keqoft+1,2,zzzvbar1,
     $         frateoft+1,2)
C     
         else
            rrate=0
            call SVZ(rrateoft+1,2)
C     
         end if
         do j=1,rnum
            zzzv3=frate*y(thirdbodyreactionarray(20*i+2+j))
            zzzvbar2=frate
            zzzvbar1=y(thirdbodyreactionarray(20*i+2+j))
            frate=zzzv3
            call DSVM2(frateoft+1,2,zzzvbar2,yoft+
     $         thirdbodyreactionarray(20*i+2+j),2,zzzvbar1,frateoft+1,
     $         2)
C     
         end do
         do j=1,pnum
            zzzv3=rrate*y(thirdbodyreactionarray(20*i+5+j))
            zzzvbar2=rrate
            zzzvbar1=y(thirdbodyreactionarray(20*i+5+j))
            rrate=zzzv3
            call DSVM2(rrateoft+1,2,zzzvbar2,yoft+
     $         thirdbodyreactionarray(20*i+5+j),2,zzzvbar1,rrateoft+1,
     $         2)
C     
         end do
         do j=1,rnum
            zzzv3=del(thirdbodyreactionarray(20*i+2+j))-frate
            zzzv5=zzzv3+rrate
            zzzvbar4=1.0d0
            zzzvbar2=-1.0d0
            zzzvbar1=1.0d0
            del(thirdbodyreactionarray(20*i+2+j))=zzzv5
            call DSVM3(deloft+thirdbodyreactionarray(20*i+2+j),zzzsvm2
     $         ,zzzvbar4,rrateoft+1,2,zzzvbar2,frateoft+1,2,zzzvbar1,
     $         deloft+thirdbodyreactionarray(20*i+2+j),zzzsvm2)
C     
         end do
         do j=1,pnum
            zzzv3=del(thirdbodyreactionarray(20*i+5+j))+frate
            zzzv5=zzzv3-rrate
            zzzvbar4=-1.0d0
            zzzvbar2=1.0d0
            zzzvbar1=1.0d0
            del(thirdbodyreactionarray(20*i+5+j))=zzzv5
            call DSVM3(deloft+thirdbodyreactionarray(20*i+5+j),zzzsvm2
     $         ,zzzvbar4,rrateoft+1,2,zzzvbar2,frateoft+1,2,zzzvbar1,
     $         deloft+thirdbodyreactionarray(20*i+5+j),zzzsvm2)
C     
         end do
      end do
      do i=0,troereactionsize-1
         rnum=troereactionarray(i*21+1)
         pnum=troereactionarray(i*21+2)
         dg=0
         do j=1,rnum
            dg=dg-thermo(troereactionarray(21*i+2+j))
         end do
         do j=1,pnum
            dg=dg+thermo(troereactionarray(21*i+5+j))
         end do
         zzzv2=-dg
         zzzv4=zzzv2*4184
         zzzv6=zzzv4/8.314
         zzzv8=zzzv6/temperature
         zzzv9=exp(zzzv8)
         zzzv11=82.053*temperature
         zzzv14=rnum-pnum
         zzzv15=zzzv11**zzzv14
         zzzv16=zzzv9*zzzv15
         zzzvbar15=zzzv9
         zzzvbar11=zzzv14*zzzv11**(zzzv14-1.0d0)*zzzvbar15
         zzzvbar9=zzzv15
         zzzvbar8=zzzv9*zzzvbar9
         zzzvbar7=-zzzv8/temperature*zzzvbar8+82.053*zzzvbar11
         keq=zzzv16
         call DSVM1(keqoft+1,2,zzzvbar7,temperatureoft+1,3)
C     
         zzzv3=pressure*1e-6
         zzzv5=zzzv3/8.314
         zzzv7=zzzv5/temperature
         zzzvbar6=-zzzv7/temperature
         zzzvbar5=1.0d0/temperature
         zzzvbar3=1.0d0/8.314*zzzvbar5
         zzzvbar1=1e-6*zzzvbar3
         m=zzzv7
         call DSVM2(moft+1,2,zzzvbar6,temperatureoft+1,3,zzzvbar1,
     $      pressureoft+1,3)
C     
         numcollider=troereactionarray(i*21+10)
         do j=1,numcollider
            zzzv5=troereactionratearray(21*i+6+j)-1
            zzzv6=y(troereactionarray(i*21+10+j))*zzzv5
            zzzv7=m+zzzv6
            zzzvbar5=y(troereactionarray(i*21+10+j))
            zzzvbar3=zzzvbar5
            zzzvbar2=zzzv5
            zzzvbar1=1.0d0
            m=zzzv7
            call DSVM3(moft+1,2,zzzvbar3,troereactionratearrayoft+21*i
     $         +6+j,3,zzzvbar2,yoft+troereactionarray(i*21+10+j),2,
     $         zzzvbar1,moft+1,2)
C     
         end do
         zzzvbar1=1.0d0
         alpha=troereactionratearray(21*i+17)
         call DSVM1(alphaoft+1,2,1.0d0,troereactionratearrayoft+21*i+
     $      17,3)
C     
         zzzvbar1=1.0d0
         tstar=troereactionratearray(21*i+18)
         call DSVM1(tstaroft+1,2,1.0d0,troereactionratearrayoft+21*i+
     $      18,3)
C     
         zzzvbar1=1.0d0
         t2star=troereactionratearray(21*i+19)
         call DSVM1(t2staroft+1,2,1.0d0,troereactionratearrayoft+21*i+
     $      19,3)
C     
         zzzvbar1=1.0d0
         t3star=troereactionratearray(21*i+20)
         call DSVM1(t3staroft+1,2,1.0d0,troereactionratearrayoft+21*i+
     $      20,3)
C     
         zzzvbar1=1.0d0
         lowrate=troereactionratearray(21*i+21)
         call DSVM1(lowrateoft+1,2,1.0d0,troereactionratearrayoft+21*i
     $      +21,3)
C     
         rate=rpar(reactionsize+thirdbodyreactionsize+i+1)
         direction=troereactionarray(21*i+9)
         if(troereactionarray(21*i+21)==0) then
            zzzv3=1-alpha
            zzzv5=-temperature
            zzzv7=zzzv5/t3star
            zzzv8=exp(zzzv7)
            zzzv9=zzzv3*zzzv8
            zzzv11=zzzv5/tstar
            zzzv12=exp(zzzv11)
            zzzv13=alpha*zzzv12
            zzzv14=zzzv9+zzzv13
            zzzv16=-t2star
            zzzv17=zzzv16/temperature
            zzzv18=exp(zzzv17)
            zzzv19=zzzv14+zzzv18
            zzzvbar17=zzzv18
            zzzvbar16=1.0d0/temperature*zzzvbar17
            zzzvbar15=-zzzvbar16
            zzzvbar12=alpha
            zzzvbar11=zzzv12*zzzvbar12
            zzzvbar10=-zzzv11/tstar*zzzvbar11
            zzzvbar8=zzzv3
            zzzvbar7=zzzv8*zzzvbar8
            zzzvbar6=-zzzv7/t3star*zzzvbar7
            zzzvbar5=1.0d0/t3star*zzzvbar7+1.0d0/tstar*zzzvbar11
            zzzvbar4=-zzzvbar5-zzzv17/temperature*zzzvbar17
            zzzvbar3=zzzv8
            zzzvbar2=-zzzvbar3+zzzv12
            fcent=zzzv19
            call DSVM5(fcentoft+1,2,zzzvbar15,t2staroft+1,2,zzzvbar10,
     $         tstaroft+1,2,zzzvbar6,t3staroft+1,2,zzzvbar4,
     $         temperatureoft+1,3,zzzvbar2,alphaoft+1,2)
C     
         else
            if(troereactionarray(21*i+21)==1) then
               zzzv3=1-alpha
               zzzv5=-temperature
               zzzv7=zzzv5/t3star
               zzzv8=exp(zzzv7)
               zzzv9=zzzv3*zzzv8
               zzzv11=zzzv5/tstar
               zzzv12=exp(zzzv11)
               zzzv13=alpha*zzzv12
               zzzv14=zzzv9+zzzv13
               zzzvbar12=alpha
               zzzvbar11=zzzv12*zzzvbar12
               zzzvbar10=-zzzv11/tstar*zzzvbar11
               zzzvbar8=zzzv3
               zzzvbar7=zzzv8*zzzvbar8
               zzzvbar6=-zzzv7/t3star*zzzvbar7
               zzzvbar5=1.0d0/t3star*zzzvbar7+1.0d0/tstar*zzzvbar11
               zzzvbar4=-zzzvbar5
               zzzvbar3=zzzv8
               zzzvbar2=-zzzvbar3+zzzv12
               fcent=zzzv14
               call DSVM4(fcentoft+1,2,zzzvbar10,tstaroft+1,2,
     $            zzzvbar6,t3staroft+1,2,zzzvbar4,temperatureoft+1,3,
     $            zzzvbar2,alphaoft+1,2)
C     
            else
               stop
            end if
         end if
         d=0.14
         zzzv3=lowrate*m
         zzzv5=zzzv3/rate
         zzzvbar3=1.0d0/rate
         zzzvbar2=lowrate*zzzvbar3
         zzzvbar1=m*zzzvbar3
         pr=zzzv5
         call DSVM2(proft+1,2,zzzvbar2,moft+1,2,zzzvbar1,lowrateoft+1,
     $      2)
C     
         if(fcent>=1e-30) then
            zzzv2=log10(fcent)
            zzzvbar1=1.0d0/(fcent*log(10.0))
            logfcent=log10(fcent)
            call DSVM1(logfcentoft+1,2,zzzvbar1,fcentoft+1,2)
C     
         else
            logfcent=-30
            call SVZ(logfcentoft+1,2)
C     
         end if
         if(pr>=1e-30) then
            zzzv2=log10(pr)
            zzzvbar1=1.0d0/(pr*log(10.0))
            logpr=log10(pr)
            call DSVM1(logproft+1,2,zzzvbar1,proft+1,2)
C     
         else
            logpr=-30
            call SVZ(logproft+1,2)
C     
         end if
         zzzv4=1.27*logfcent
         zzzv5=0.75-zzzv4
         zzzvbar4=-1.0d0
         zzzvbar3=1.27*zzzvbar4
         n=zzzv5
         call DSVM1(noft+1,2,zzzvbar3,logfcentoft+1,2)
C     
         zzzv4=0.67*logfcent
         zzzv5=-0.4-zzzv4
         zzzvbar4=-1.0d0
         zzzvbar3=0.67*zzzvbar4
         c=zzzv5
         call DSVM1(coft+1,2,zzzvbar3,logfcentoft+1,2)
C     
         zzzv3=logpr+c
         zzzv6=d*zzzv3
         zzzv7=n-zzzv6
         zzzv8=zzzv3/zzzv7
         zzzvbar7=-zzzv8/zzzv7
         zzzvbar6=-zzzvbar7
         zzzvbar4=zzzvbar7
         zzzvbar3=1.0d0/zzzv7+d*zzzvbar6
         zzzvbar2=zzzvbar3
         zzzvbar1=zzzvbar3
         inside=zzzv8
         call DSVM3(insideoft+1,2,zzzvbar4,noft+1,2,zzzvbar2,coft+1,2,
     $      zzzvbar1,logproft+1,2)
C     
         zzzv5=inside**2
         zzzv6=1+zzzv5
         zzzv7=logfcent/zzzv6
         zzzvbar6=-zzzv7/zzzv6
         zzzvbar5=zzzvbar6
         zzzvbar3=2*inside*zzzvbar5
         zzzvbar1=1.0d0/zzzv6
         logf=zzzv7
         call DSVM2(logfoft+1,2,zzzvbar3,insideoft+1,2,zzzvbar1,
     $      logfcentoft+1,2)
C     
         zzzv3=10**logf
         zzzvbar2=zzzv3*log(10.0)
         f=zzzv3
         call DSVM1(foft+1,2,zzzvbar2,logfoft+1,2)
C     
         if (rnum .gt. 1 .and. pnum .gt. 1) then
            zzzv4=1+pr
            zzzv5=1.0d0/zzzv4
            zzzv6=lowrate*zzzv5
            zzzv8=zzzv6*f
            zzzvbar7=zzzv6
            zzzvbar6=f
            zzzvbar5=lowrate*zzzvbar6
            zzzvbar4=-zzzv5/zzzv4*zzzvbar5
            zzzvbar2=zzzvbar4
            frate=zzzv8
            call DSVM2(frateoft+1,2,zzzvbar7,foft+1,2,
     $         zzzvbar2,proft+1,2)
         else
            zzzv4=1+pr
            zzzv5=pr/zzzv4
            zzzv6=rate*zzzv5
            zzzv8=zzzv6*f
            zzzvbar7=zzzv6
            zzzvbar6=f
            zzzvbar5=rate*zzzvbar6
            zzzvbar4=-zzzv5/zzzv4*zzzvbar5
            zzzvbar2=1.0d0/zzzv4*zzzvbar5+zzzvbar4
            frate=zzzv8
            call DSVM2(frateoft+1,2,zzzvbar7,foft+1,2,
     $         zzzvbar2,proft+1,2)
         end if
C     
         if(troereactionarray(21*i+9)==1) then
            zzzv3=frate/keq
            zzzvbar2=-zzzv3/keq
            zzzvbar1=1.0d0/keq
            rrate=zzzv3
            call DSVM2(rrateoft+1,2,zzzvbar2,keqoft+1,2,zzzvbar1,
     $         frateoft+1,2)
C     
         else
            rrate=0
            call SVZ(rrateoft+1,2)
C     
         end if
         do j=1,rnum
            zzzv3=frate*y(troereactionarray(21*i+2+j))
            zzzvbar2=frate
            zzzvbar1=y(troereactionarray(21*i+2+j))
            frate=zzzv3
            call DSVM2(frateoft+1,2,zzzvbar2,yoft+troereactionarray(21
     $         *i+2+j),2,zzzvbar1,frateoft+1,2)
C     
         end do
         do j=1,pnum
            zzzv3=rrate*y(troereactionarray(21*i+5+j))
            zzzvbar2=rrate
            zzzvbar1=y(troereactionarray(21*i+5+j))
            rrate=zzzv3
            call DSVM2(rrateoft+1,2,zzzvbar2,yoft+troereactionarray(21
     $         *i+5+j),2,zzzvbar1,rrateoft+1,2)
C     
         end do
         do j=1,rnum
            zzzv3=del(troereactionarray(21*i+2+j))-frate
            zzzv5=zzzv3+rrate
            zzzvbar4=1.0d0
            zzzvbar2=-1.0d0
            zzzvbar1=1.0d0
            del(troereactionarray(21*i+2+j))=zzzv5
            call DSVM3(deloft+troereactionarray(21*i+2+j),zzzsvm2,
     $         zzzvbar4,rrateoft+1,2,zzzvbar2,frateoft+1,2,zzzvbar1,
     $         deloft+troereactionarray(21*i+2+j),zzzsvm2)
C     
         end do
         do j=1,pnum
            zzzv3=del(troereactionarray(21*i+5+j))+frate
            zzzv5=zzzv3-rrate
            zzzvbar4=-1.0d0
            zzzvbar2=1.0d0
            zzzvbar1=1.0d0
            del(troereactionarray(21*i+5+j))=zzzv5
            call DSVM3(deloft+troereactionarray(21*i+5+j),zzzsvm2,
     $         zzzvbar4,rrateoft+1,2,zzzvbar2,frateoft+1,2,zzzvbar1,
     $         deloft+troereactionarray(21*i+5+j),zzzsvm2)
C     
         end do
      end do

      do i=0,lindereactionsize-1
         rnum=lindereactionarray(i*20+1)
         pnum=lindereactionarray(i*20+2)
         dg=0
         do j=1,rnum
            dg=dg-thermo(lindereactionarray(20*i+2+j))
         end do
         do j=1,pnum
            dg=dg+thermo(lindereactionarray(20*i+5+j))
         end do
         zzzv2=-dg
         zzzv4=zzzv2*4184
         zzzv6=zzzv4/8.314
         zzzv8=zzzv6/temperature
         zzzv9=exp(zzzv8)
         zzzv11=82.053*temperature
         zzzv14=rnum-pnum
         zzzv15=zzzv11**zzzv14
         zzzv16=zzzv9*zzzv15
         zzzvbar15=zzzv9
         zzzvbar11=zzzv14*zzzv11**(zzzv14-1.0d0)*zzzvbar15
         zzzvbar9=zzzv15
         zzzvbar8=zzzv9*zzzvbar9
         zzzvbar7=-zzzv8/temperature*zzzvbar8+82.053*zzzvbar11
         keq=zzzv16
         call DSVM1(keqoft+1,2,zzzvbar7,temperatureoft+1,3)
C     
         zzzv3=pressure*1e-6
         zzzv5=zzzv3/8.314
         zzzv7=zzzv5/temperature
         zzzvbar6=-zzzv7/temperature
         zzzvbar5=1.0d0/temperature
         zzzvbar3=1.0d0/8.314*zzzvbar5
         zzzvbar1=1e-6*zzzvbar3
         m=zzzv7
         call DSVM2(moft+1,2,zzzvbar6,temperatureoft+1,3,zzzvbar1,
     $      pressureoft+1,3)
C     
         numcollider=lindereactionarray(i*20+10)
         do j=1,numcollider
            zzzv5=lindereactionratearray(17*i+6+j)-1
            zzzv6=y(lindereactionarray(i*20+10+j))*zzzv5
            zzzv7=m+zzzv6
            zzzvbar5=y(lindereactionarray(i*20+10+j))
            zzzvbar3=zzzvbar5
            zzzvbar2=zzzv5
            zzzvbar1=1.0d0
            m=zzzv7
            call DSVM3(moft+1,2,zzzvbar3,lindereactionratearrayoft+17*i
     $         +6+j,3,zzzvbar2,yoft+lindereactionarray(i*20+10+j),2,
     $         zzzvbar1,moft+1,2)
C     
         end do

         zzzvbar1=1.0d0
         lowrate=lindereactionratearray(20*i+20)
         call DSVM1(lowrateoft+1,2,1.0d0,lindereactionratearrayoft+17*i
     $      +17,3)
C     
         rate=rpar(reactionsize+thirdbodyreactionsize+troereactionsize
     $      +i+1)
         direction=lindereactionarray(20*i+9)
C     
         if (rnum .gt. 1 .and. pnum .gt. 1) then
            zzzv4=1+pr
            zzzv5=1.0d0/zzzv4
            zzzvbar5=lowrate
            zzzvbar4=-zzzv5/zzzv4*zzzvbar5
            zzzvbar2=zzzvbar4
            frate=zzzv8
            call DSVM1(frateoft+1,2,zzzvbar2,proft+1,2)
         else
            zzzv4=1+pr
            zzzv5=pr/zzzv4
            zzzvbar5=rate
            zzzvbar4=-zzzv5/zzzv4*zzzvbar5
            zzzvbar2=1.0d0/zzzv4*zzzvbar5+zzzvbar4
            frate=zzzv8
            call DSVM1(frateoft+1,2,zzzvbar2,proft+1,2)
         end if
C     
         if(lindereactionarray(20*i+9)==1) then
            zzzv3=frate/keq
            zzzvbar2=-zzzv3/keq
            zzzvbar1=1.0d0/keq
            rrate=zzzv3
            call DSVM2(rrateoft+1,2,zzzvbar2,keqoft+1,2,zzzvbar1,
     $         frateoft+1,2)
C     
         else
            rrate=0
            call SVZ(rrateoft+1,2)
C     
         end if
         do j=1,rnum
            zzzv3=frate*y(lindereactionarray(20*i+2+j))
            zzzvbar2=frate
            zzzvbar1=y(lindereactionarray(20*i+2+j))
            frate=zzzv3
            call DSVM2(frateoft+1,2,zzzvbar2,yoft+lindereactionarray(20
     $         *i+2+j),2,zzzvbar1,frateoft+1,2)
C     
         end do
         do j=1,pnum
            zzzv3=rrate*y(lindereactionarray(20*i+5+j))
            zzzvbar2=rrate
            zzzvbar1=y(lindereactionarray(20*i+5+j))
            rrate=zzzv3
            call DSVM2(rrateoft+1,2,zzzvbar2,yoft+lindereactionarray(20
     $         *i+5+j),2,zzzvbar1,rrateoft+1,2)
C     
         end do
         do j=1,rnum
            zzzv3=del(lindereactionarray(20*i+2+j))-frate
            zzzv5=zzzv3+rrate
            zzzvbar4=1.0d0
            zzzvbar2=-1.0d0
            zzzvbar1=1.0d0
            del(lindereactionarray(20*i+2+j))=zzzv5
            call DSVM3(deloft+lindereactionarray(20*i+2+j),zzzsvm2,
     $         zzzvbar4,rrateoft+1,2,zzzvbar2,frateoft+1,2,zzzvbar1,
     $         deloft+lindereactionarray(20*i+2+j),zzzsvm2)
C     
         end do
         do j=1,pnum
            zzzv3=del(lindereactionarray(20*i+5+j))+frate
            zzzv5=zzzv3-rrate
            zzzvbar4=-1.0d0
            zzzvbar2=1.0d0
            zzzvbar1=1.0d0
            del(lindereactionarray(20*i+5+j))=zzzv5
            call DSVM3(deloft+lindereactionarray(20*i+5+j),zzzsvm2,
     $         zzzvbar4,rrateoft+1,2,zzzvbar2,frateoft+1,2,zzzvbar1,
     $         deloft+lindereactionarray(20*i+5+j),zzzsvm2)
C     
         end do
      end do

      sumyprime=0
      call SVZ(sumyprimeoft+1,2)
C     
      sumy=0
      do i=1,nstate-1
         zzzv3=del(i)*tempy(nstate)
         zzzvbar2=del(i)
         zzzvbar1=tempy(nstate)
         del(i)=zzzv3
         call DSVM2(deloft+i,zzzsvm2,zzzvbar2,tempyoft+nstate,zzzsvm1,
     $      zzzvbar1,deloft+i,zzzsvm2)
C     
      end do
      do i=1,nstate-1
         zzzv3=sumyprime+del(i)
         zzzvbar2=1.0d0
         zzzvbar1=1.0d0
         sumyprime=zzzv3
         call DSVM2(sumyprimeoft+1,2,zzzvbar2,deloft+i,zzzsvm2,
     $      zzzvbar1,sumyprimeoft+1,2)
C     
      end do
      zzzv3=sumyprime*8.314
      zzzv5=zzzv3*temperature
      zzzv7=zzzv5/pressure
      zzzv9=zzzv7/1e-6
      zzzvbar7=1.0d0/1e-6
      zzzvbar6=-zzzv7/pressure*zzzvbar7
      zzzvbar5=1.0d0/pressure*zzzvbar7
      zzzvbar4=zzzv3*zzzvbar5
      zzzvbar3=temperature*zzzvbar5
      zzzvbar1=8.314*zzzvbar3
      del(nstate)=zzzv9
      call DSVM3(deloft+nstate,zzzsvm2,zzzvbar6,pressureoft+1,3,
     $   zzzvbar4,temperatureoft+1,3,zzzvbar1,sumyprimeoft+1,2)
C     
C     Done
      end
