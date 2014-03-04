C     - DAEPACK Version 1.1 - Copyright (C) Numerica Technology, LLC
C     - DERIVATIVE COMPUTATION -
C     -   GENERATED: Wed Jan 31 02:09:44 2007
      subroutine res_rpar(t,y,yprime,cj,del,rpar,ipar,senpar,
     $   zzzderiv,zzzne,zzzirn,zzzjcn,zzziw)
!!independent { rpar }
!!dependent   { del(1:nstate) }
!DEC$ IF DEFINED(DLLEXP)
!DEC$ ATTRIBUTES DLLEXPORT :: res_rpar
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
     $   21*troemax),temperature,pressure,reactionarray(10*rmax),
     $   thirdbodyreactionarray(21*tbrmax),troereactionarray(22*
     $   troemax),lindereactionarray(21*lindemax),
     $   lindereactionratearray(17*lindemax)
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
      integer rparoft
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
      rparoft=deloft+neq
C     Variable offsets for COMMON BLOCK size
C     Variable offsets for COMMON BLOCK reac
      reactionratearrayoft=0
      thirdbodyreactionratearrayoft=reactionratearrayoft+(5*rmax)
      troereactionratearrayoft=thirdbodyreactionratearrayoft+(16*tbrmax)
      lindereactionratearrayoft=troereactionratearrayoft+(21*troemax)
      temperatureoft=lindereactionratearrayoft+(17*lindemax)
      pressureoft=temperatureoft+1
C     Number of independent and dependent variables.
      zzzn=reactionsize+thirdbodyreactionsize+troereactionsize+
     $   lindereactionsize+nstate-1
      zzzm=nstate
C     Create unique identifier for this model.
      call GETMDLID(zzzimodel,'res_rpar')
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
      do izzz1=1,reactionsize+thirdbodyreactionsize+troereactionsize+
     $   lindereactionsize+nstate-1
         zzzindx=zzzindx+1
         zzziw(zzzindx)=rparoft+izzz1
      end do
C     Initialize independent variable gradients to Cartesian basis vectors.
      call DSETIVCV(1,zzzn,zzziw)
C     
C     Save SVMs of arguments
      zzziw(1)=1
      zzziw(2)=1
C     Save current offsets for use within embedded subroutine
      zzziw(3)=deloft
      zzziw(4)=rparoft
C     Call modifed embedded subroutine
      call getflux_rpar(y,del,rpar,zzziw)
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
      subroutine getflux_rpar(tempy,del,rpar,zzziw)
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
     $   21*troemax),temperature,pressure,reactionarray(10*rmax),
     $   thirdbodyreactionarray(21*tbrmax),troereactionarray(22*
     $   troemax),lindereactionratearray(17*lindemax),
     $   lindereactionarray(21*lindemax)
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
      integer deloft
      integer rparoft
      integer reactionratearrayoft
      integer thirdbodyreactionratearrayoft
      integer troereactionratearrayoft
      integer lindereactionratearrayoft
      integer temperatureoft
      integer pressureoft
      integer thermooft
      integer frateoft
      integer dgoft
      integer keqoft
      integer rrateoft
      integer inertefficiencyoft
      integer moft
      integer alphaoft
      integer tstaroft
      integer t2staroft
      integer t3staroft
      integer lowrateoft
      integer rateoft
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
      deloft=zzziw(3)
      rparoft=zzziw(4)
C     
      thermooft=0
      frateoft=thermooft+spmax
      dgoft=frateoft+1
      keqoft=dgoft+1
      rrateoft=keqoft+1
      inertefficiencyoft=rrateoft+1
      moft=inertefficiencyoft+1
      alphaoft=moft+1
      tstaroft=alphaoft+1
      t2staroft=tstaroft+1
      t3staroft=t2staroft+1
      lowrateoft=t3staroft+1
      rateoft=lowrateoft+1
      fcentoft=rateoft+1
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
         call SVZ(deloft+i,zzzsvm1)
C     
      end do
      do i=1,nstate-1
         y(i)=tempy(i)/tempy(nstate)
      end do
      y(nstate)=tempy(nstate)
      do i=1,nstate-1
         zzzvbar1=1.0d0
         thermo(i)=rpar(reactionsize+thirdbodyreactionsize+
     $      troereactionsize+lindereactionsize+i)
         call DSVM1(thermooft+i,2,1.0d0,rparoft+reactionsize+
     $      thirdbodyreactionsize+troereactionsize+lindereactionsize
     $      +i,zzzsvm2)
C     
      end do
      do i=0,reactionsize-1
         zzzvbar1=1.0d0
         frate=rpar(i+1)
         call DSVM1(frateoft+1,2,1.0d0,rparoft+i+1,zzzsvm2)
C     
         rnum=reactionarray(10*i+1)
         pnum=reactionarray(10*i+2)
         dg=0
         call SVZ(dgoft+1,2)
C     
         do j=1,rnum
            zzzv3=dg-thermo(reactionarray(10*i+2+j))
            zzzvbar2=-1.0d0
            zzzvbar1=1.0d0
            dg=zzzv3
            call DSVM2(dgoft+1,2,zzzvbar2,thermooft+reactionarray(10*i+
     $         2+j),2,zzzvbar1,dgoft+1,2)
C     
         end do
         do j=1,pnum
            zzzv3=dg+thermo(reactionarray(10*i+5+j))
            zzzvbar2=1.0d0
            zzzvbar1=1.0d0
            dg=zzzv3
            call DSVM2(dgoft+1,2,zzzvbar2,thermooft+reactionarray(10*i+
     $         5+j),2,zzzvbar1,dgoft+1,2)
C     
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
         zzzvbar6=1.0d0/temperature*zzzvbar8
         zzzvbar4=1.0d0/8.314*zzzvbar6
         zzzvbar2=4184*zzzvbar4
         zzzvbar1=-zzzvbar2
         keq=zzzv16
         call DSVM2(keqoft+1,2,zzzvbar7,temperatureoft+1,3,zzzvbar1,
     $      dgoft+1,2)
C     
         if(reactionarray(10*i+10)==1) then
            zzzv3=rpar(i+1)/keq
            zzzvbar2=-zzzv3/keq
            zzzvbar1=1.0d0/keq
            rrate=zzzv3
            call DSVM2(rrateoft+1,2,zzzvbar2,keqoft+1,2,zzzvbar1,
     $         rparoft+i+1,zzzsvm2)
C     
         else
            rrate=0
            call SVZ(rrateoft+1,2)
C     
         end if
         do j=1,rnum
            zzzv3=frate*y(reactionarray(10*i+2+j))
            zzzvbar1=y(reactionarray(10*i+2+j))
            frate=zzzv3
            call DSVM1(frateoft+1,2,zzzvbar1,frateoft+1,2)
C     
         end do
         do j=1,pnum
            zzzv3=rrate*y(reactionarray(10*i+5+j))
            zzzvbar1=y(reactionarray(10*i+5+j))
            rrate=zzzv3
            call DSVM1(rrateoft+1,2,zzzvbar1,rrateoft+1,2)
C     
         end do
         do j=1,rnum
            zzzv3=del(reactionarray(10*i+2+j))-frate
            zzzv5=zzzv3+rrate
            zzzvbar4=1.0d0
            zzzvbar2=-1.0d0
            zzzvbar1=1.0d0
            del(reactionarray(10*i+2+j))=zzzv5
            call DSVM3(deloft+reactionarray(10*i+2+j),zzzsvm1,zzzvbar4,
     $         rrateoft+1,2,zzzvbar2,frateoft+1,2,zzzvbar1,deloft+
     $         reactionarray(10*i+2+j),zzzsvm1)
C     
         end do
         do j=1,pnum
            zzzv3=del(reactionarray(10*i+5+j))+frate
            zzzv5=zzzv3-rrate
            zzzvbar4=-1.0d0
            zzzvbar2=1.0d0
            zzzvbar1=1.0d0
            del(reactionarray(10*i+5+j))=zzzv5
            call DSVM3(deloft+reactionarray(10*i+5+j),zzzsvm1,zzzvbar4,
     $         rrateoft+1,2,zzzvbar2,frateoft+1,2,zzzvbar1,deloft+
     $         reactionarray(10*i+5+j),zzzsvm1)
C     
         end do
      end do
      do i=0,thirdbodyreactionsize-1
         zzzvbar1=1.0d0
         frate=rpar(reactionsize+i+1)
         call DSVM1(frateoft+1,2,1.0d0,rparoft+reactionsize+i+1,
     $      zzzsvm2)
C     
         rnum=thirdbodyreactionarray(i*21+1)
         pnum=thirdbodyreactionarray(i*21+2)
         dg=0
         call SVZ(dgoft+1,2)
C     
         do j=1,rnum
            zzzv3=dg-thermo(thirdbodyreactionarray(21*i+2+j))
            zzzvbar2=-1.0d0
            zzzvbar1=1.0d0
            dg=zzzv3
            call DSVM2(dgoft+1,2,zzzvbar2,thermooft+
     $         thirdbodyreactionarray(21*i+2+j),2,zzzvbar1,dgoft+1,2)
C     
         end do
         do j=1,pnum
            zzzv3=dg+thermo(thirdbodyreactionarray(21*i+5+j))
            zzzvbar2=1.0d0
            zzzvbar1=1.0d0
            dg=zzzv3
            call DSVM2(dgoft+1,2,zzzvbar2,thermooft+
     $         thirdbodyreactionarray(21*i+5+j),2,zzzvbar1,dgoft+1,2)
C     
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
         zzzvbar6=1.0d0/temperature*zzzvbar8
         zzzvbar4=1.0d0/8.314*zzzvbar6
         zzzvbar2=4184*zzzvbar4
         zzzvbar1=-zzzvbar2
         keq=zzzv16
         call DSVM2(keqoft+1,2,zzzvbar7,temperatureoft+1,3,zzzvbar1,
     $      dgoft+1,2)
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
         numcollider=thirdbodyreactionarray(i*21+11)
         do j=1,numcollider
            zzzv5=thirdbodyreactionratearray(16*i+6+j)-1
            zzzv6=y(thirdbodyreactionarray(i*21+11+j))*zzzv5
            zzzv7=inertefficiency+zzzv6
            zzzvbar5=y(thirdbodyreactionarray(i*21+11+j))
            zzzvbar3=zzzvbar5
            zzzvbar1=1.0d0
            inertefficiency=zzzv7
            call DSVM2(inertefficiencyoft+1,2,zzzvbar3,
     $         thirdbodyreactionratearrayoft+16*i+6+j,3,zzzvbar1,
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
         if(thirdbodyreactionarray(21*i+10)==1) then
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
            zzzv3=frate*y(thirdbodyreactionarray(21*i+2+j))
            zzzvbar1=y(thirdbodyreactionarray(21*i+2+j))
            frate=zzzv3
            call DSVM1(frateoft+1,2,zzzvbar1,frateoft+1,2)
C     
         end do
         do j=1,pnum
            zzzv3=rrate*y(thirdbodyreactionarray(21*i+5+j))
            zzzvbar1=y(thirdbodyreactionarray(21*i+5+j))
            rrate=zzzv3
            call DSVM1(rrateoft+1,2,zzzvbar1,rrateoft+1,2)
C     
         end do
         do j=1,rnum
            zzzv3=del(thirdbodyreactionarray(21*i+2+j))-frate
            zzzv5=zzzv3+rrate
            zzzvbar4=1.0d0
            zzzvbar2=-1.0d0
            zzzvbar1=1.0d0
            del(thirdbodyreactionarray(21*i+2+j))=zzzv5
            call DSVM3(deloft+thirdbodyreactionarray(21*i+2+j),zzzsvm1
     $         ,zzzvbar4,rrateoft+1,2,zzzvbar2,frateoft+1,2,zzzvbar1,
     $         deloft+thirdbodyreactionarray(21*i+2+j),zzzsvm1)
C     
         end do
         do j=1,pnum
            zzzv3=del(thirdbodyreactionarray(21*i+5+j))+frate
            zzzv5=zzzv3-rrate
            zzzvbar4=-1.0d0
            zzzvbar2=1.0d0
            zzzvbar1=1.0d0
            del(thirdbodyreactionarray(21*i+5+j))=zzzv5
            call DSVM3(deloft+thirdbodyreactionarray(21*i+5+j),zzzsvm1
     $         ,zzzvbar4,rrateoft+1,2,zzzvbar2,frateoft+1,2,zzzvbar1,
     $         deloft+thirdbodyreactionarray(21*i+5+j),zzzsvm1)
C     
         end do
      end do
      do i=0,troereactionsize-1
         rnum=troereactionarray(i*22+1)
         pnum=troereactionarray(i*22+2)
         dg=0
         call SVZ(dgoft+1,2)
C     
         do j=1,rnum
            zzzv3=dg-thermo(troereactionarray(22*i+2+j))
            zzzvbar2=-1.0d0
            zzzvbar1=1.0d0
            dg=zzzv3
            call DSVM2(dgoft+1,2,zzzvbar2,thermooft+troereactionarray(
     $         21*i+2+j),2,zzzvbar1,dgoft+1,2)
C     
         end do
         do j=1,pnum
            zzzv3=dg+thermo(troereactionarray(22*i+5+j))
            zzzvbar2=1.0d0
            zzzvbar1=1.0d0
            dg=zzzv3
            call DSVM2(dgoft+1,2,zzzvbar2,thermooft+troereactionarray(
     $         21*i+5+j),2,zzzvbar1,dgoft+1,2)
C     
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
         zzzvbar6=1.0d0/temperature*zzzvbar8
         zzzvbar4=1.0d0/8.314*zzzvbar6
         zzzvbar2=4184*zzzvbar4
         zzzvbar1=-zzzvbar2
         keq=zzzv16
         call DSVM2(keqoft+1,2,zzzvbar7,temperatureoft+1,3,zzzvbar1,
     $      dgoft+1,2)
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
         numcollider=troereactionarray(i*22+11)
         do j=1,numcollider
            zzzv5=troereactionratearray(21*i+6+j)-1
            zzzv6=y(troereactionarray(i*22+11+j))*zzzv5
            zzzv7=m+zzzv6
            zzzvbar5=y(troereactionarray(i*22+11+j))
            zzzvbar3=zzzvbar5
            zzzvbar1=1.0d0
            m=zzzv7
            call DSVM2(moft+1,2,zzzvbar3,troereactionratearrayoft+21*i
     $         +6+j,3,zzzvbar1,moft+1,2)
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
         zzzvbar1=1.0d0
         rate=rpar(reactionsize+thirdbodyreactionsize+i+1)
         call DSVM1(rateoft+1,2,1.0d0,rparoft+reactionsize+
     $      thirdbodyreactionsize+i+1,zzzsvm2)
C     
         direction=troereactionarray(22*i+10)
         if(troereactionarray(22*i+22)==0) then
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
            if(troereactionarray(22*i+22)==1) then
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
         zzzvbar4=-zzzv5/rate
         zzzvbar3=1.0d0/rate
         zzzvbar2=lowrate*zzzvbar3
         zzzvbar1=m*zzzvbar3
         pr=zzzv5
         call DSVM3(proft+1,2,zzzvbar4,rateoft+1,2,zzzvbar2,moft+1,2,
     $      zzzvbar1,lowrateoft+1,2)
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
         IF (RNUM .GT. 1 .AND. PNUM .GT. 1) THEN
            zzzv4=1+pr
            zzzv5=1.0d0/zzzv4
            zzzv6=lowrate*zzzv5
            zzzv8=zzzv6*f
            zzzvbar7=zzzv6
            zzzvbar6=f
            zzzvbar5=lowrate*zzzvbar6
            zzzvbar4=-zzzv5/zzzv4*zzzvbar5
            zzzvbar2=zzzvbar4
            zzzvbar1=zzzv5*zzzvbar6
            frate=zzzv8
            call DSVM3(frateoft+1,2,zzzvbar7,foft+1,2,zzzvbar2,
     $         proft+1,2,zzzvbar1,lowrateoft+1,2)
         ELSE
            zzzv4=1+pr
            zzzv5=pr/zzzv4
            zzzv6=rate*zzzv5
            zzzv8=zzzv6*f
            zzzvbar7=zzzv6
            zzzvbar6=f
            zzzvbar5=rate*zzzvbar6
            zzzvbar4=-zzzv5/zzzv4*zzzvbar5
            zzzvbar2=1.0d0/zzzv4*zzzvbar5+zzzvbar4
            zzzvbar1=zzzv5*zzzvbar6
            frate=zzzv8
            call DSVM3(frateoft+1,2,zzzvbar7,foft+1,2,zzzvbar2,
     $      proft+1,2,zzzvbar1,rateoft+1,2)
         END IF

         zzzv4=1+pr
         zzzv5=pr/zzzv4
         zzzv6=rate*zzzv5
         zzzv8=zzzv6*f
         zzzvbar7=zzzv6
         zzzvbar6=f
         zzzvbar5=rate*zzzvbar6
         zzzvbar4=-zzzv5/zzzv4*zzzvbar5
         zzzvbar2=1.0d0/zzzv4*zzzvbar5+zzzvbar4
         zzzvbar1=zzzv5*zzzvbar6
         frate=zzzv8
         call DSVM3(frateoft+1,2,zzzvbar7,foft+1,2,zzzvbar2,proft+1,2,
     $      zzzvbar1,rateoft+1,2)
C     
         if(troereactionarray(22*i+10)==1) then
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
            zzzv3=frate*y(troereactionarray(22*i+2+j))
            zzzvbar1=y(troereactionarray(22*i+2+j))
            frate=zzzv3
            call DSVM1(frateoft+1,2,zzzvbar1,frateoft+1,2)
C     
         end do
         do j=1,pnum
            zzzv3=rrate*y(troereactionarray(22*i+5+j))
            zzzvbar1=y(troereactionarray(22*i+5+j))
            rrate=zzzv3
            call DSVM1(rrateoft+1,2,zzzvbar1,rrateoft+1,2)
C     
         end do
         do j=1,rnum
            zzzv3=del(troereactionarray(22*i+2+j))-frate
            zzzv5=zzzv3+rrate
            zzzvbar4=1.0d0
            zzzvbar2=-1.0d0
            zzzvbar1=1.0d0
            del(troereactionarray(22*i+2+j))=zzzv5
            call DSVM3(deloft+troereactionarray(22*i+2+j),zzzsvm1,
     $         zzzvbar4,rrateoft+1,2,zzzvbar2,frateoft+1,2,zzzvbar1,
     $         deloft+troereactionarray(22*i+2+j),zzzsvm1)
C     
         end do
         do j=1,pnum
            zzzv3=del(troereactionarray(22*i+5+j))+frate
            zzzv5=zzzv3-rrate
            zzzvbar4=-1.0d0
            zzzvbar2=1.0d0
            zzzvbar1=1.0d0
            del(troereactionarray(22*i+5+j))=zzzv5
            call DSVM3(deloft+troereactionarray(22*i+5+j),zzzsvm1,
     $         zzzvbar4,rrateoft+1,2,zzzvbar2,frateoft+1,2,zzzvbar1,
     $         deloft+troereactionarray(22*i+5+j),zzzsvm1)
C     
         end do
      end do

      do i=0,lindereactionsize-1
         rnum=lindereactionarray(i*21+1)
         pnum=lindereactionarray(i*21+2)
         dg=0
         call SVZ(dgoft+1,2)
C     
         do j=1,rnum
            zzzv3=dg-thermo(lindereactionarray(21*i+2+j))
            zzzvbar2=-1.0d0
            zzzvbar1=1.0d0
            dg=zzzv3
            call DSVM2(dgoft+1,2,zzzvbar2,thermooft+lindereactionarray(
     $         20*i+2+j),2,zzzvbar1,dgoft+1,2)
C     
         end do
         do j=1,pnum
            zzzv3=dg+thermo(lindereactionarray(21*i+5+j))
            zzzvbar2=1.0d0
            zzzvbar1=1.0d0
            dg=zzzv3
            call DSVM2(dgoft+1,2,zzzvbar2,thermooft+lindereactionarray(
     $         20*i+5+j),2,zzzvbar1,dgoft+1,2)
C     
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
         zzzvbar6=1.0d0/temperature*zzzvbar8
         zzzvbar4=1.0d0/8.314*zzzvbar6
         zzzvbar2=4184*zzzvbar4
         zzzvbar1=-zzzvbar2
         keq=zzzv16
         call DSVM2(keqoft+1,2,zzzvbar7,temperatureoft+1,3,zzzvbar1,
     $      dgoft+1,2)
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
         numcollider=lindereactionarray(i*21+11)
         do j=1,numcollider
            zzzv5=lindereactionratearray(17*i+6+j)-1
            zzzv6=y(lindereactionarray(i*21+11+j))*zzzv5
            zzzv7=m+zzzv6
            zzzvbar5=y(lindereactionarray(i*21+11+j))
            zzzvbar3=zzzvbar5
            zzzvbar1=1.0d0
            m=zzzv7
            call DSVM2(moft+1,2,zzzvbar3,lindereactionratearrayoft+17*i
     $         +6+j,3,zzzvbar1,moft+1,2)
C     
         end do

         zzzvbar1=1.0d0
         lowrate=lindereactionratearray(17*i+17)
         call DSVM1(lowrateoft+1,2,1.0d0,lindereactionratearrayoft+17*i
     $      +17,3)
C     
         zzzvbar1=1.0d0
         rate=rpar(reactionsize+thirdbodyreactionsize+troereactionsize
     $      +i+1)
         call DSVM1(rateoft+1,2,1.0d0,rparoft+reactionsize+
     $      thirdbodyreactionsize+troereactionsize+i+1,zzzsvm2)
C     
         direction=lindereactionarray(21*i+10)

         zzzv3=lowrate*m
         zzzv5=zzzv3/rate
         zzzvbar4=-zzzv5/rate
         zzzvbar3=1.0d0/rate
         zzzvbar2=lowrate*zzzvbar3
         zzzvbar1=m*zzzvbar3
         pr=zzzv5
         call DSVM3(proft+1,2,zzzvbar4,rateoft+1,2,zzzvbar2,moft+1,2,
     $      zzzvbar1,lowrateoft+1,2)
C     
         IF (RNUM .GT. 1 .AND. PNUM .GT. 1) THEN
            zzzv4=1+pr
            zzzv5=1.0d0/zzzv4
            zzzvbar5=lowrate
            zzzvbar4=-zzzv5/zzzv4*zzzvbar5
            zzzvbar2=zzzvbar4
            zzzvbar1=zzzv5
            frate=zzzv5*zzzvbar5
            call DSVM2(frateoft+1,2,zzzvbar2,proft+1,2,
     $         zzzvbar1,lowrateoft+1,2)
         ELSE
            zzzv4=1+pr
            zzzv5=pr/zzzv4
            zzzvbar5=rate
            zzzvbar4=-zzzv5/zzzv4*zzzvbar5
            zzzvbar2=1.0d0/zzzv4*zzzvbar5+zzzvbar4
            zzzvbar1=zzzv5
            frate=zzzv5*zzzvbar5
            call DSVM2(frateoft+1,2,zzzvbar2,proft+1,2,
     $         zzzvbar1,rateoft+1,2)
         END IF
C     
         if(lindereactionarray(21*i+10)==1) then
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
            zzzv3=frate*y(lindereactionarray(21*i+2+j))
            zzzvbar1=y(lindereactionarray(21*i+2+j))
            frate=zzzv3
            call DSVM1(frateoft+1,2,zzzvbar1,frateoft+1,2)
C     
         end do
         do j=1,pnum
            zzzv3=rrate*y(lindereactionarray(21*i+5+j))
            zzzvbar1=y(lindereactionarray(21*i+5+j))
            rrate=zzzv3
            call DSVM1(rrateoft+1,2,zzzvbar1,rrateoft+1,2)
C     
         end do
         do j=1,rnum
            zzzv3=del(lindereactionarray(21*i+2+j))-frate
            zzzv5=zzzv3+rrate
            zzzvbar4=1.0d0
            zzzvbar2=-1.0d0
            zzzvbar1=1.0d0
            del(lindereactionarray(21*i+2+j))=zzzv5
            call DSVM3(deloft+lindereactionarray(21*i+2+j),zzzsvm1,
     $         zzzvbar4,rrateoft+1,2,zzzvbar2,frateoft+1,2,zzzvbar1,
     $         deloft+lindereactionarray(21*i+2+j),zzzsvm1)
C     
         end do
         do j=1,pnum
            zzzv3=del(lindereactionarray(21*i+5+j))+frate
            zzzv5=zzzv3-rrate
            zzzvbar4=-1.0d0
            zzzvbar2=1.0d0
            zzzvbar1=1.0d0
            del(lindereactionarray(21*i+5+j))=zzzv5
            call DSVM3(deloft+lindereactionarray(21*i+5+j),zzzsvm1,
     $         zzzvbar4,rrateoft+1,2,zzzvbar2,frateoft+1,2,zzzvbar1,
     $         deloft+lindereactionarray(21*i+5+j),zzzsvm1)
C     
         end do
      end do

      sumyprime=0
      call SVZ(sumyprimeoft+1,2)
C     
      sumy=0
      do i=1,nstate-1
         zzzv3=del(i)*tempy(nstate)
         zzzvbar1=tempy(nstate)
         del(i)=zzzv3
         call DSVM1(deloft+i,zzzsvm1,zzzvbar1,deloft+i,zzzsvm1)
C     
      end do
      do i=1,nstate-1
         zzzv3=sumyprime+del(i)
         zzzvbar2=1.0d0
         zzzvbar1=1.0d0
         sumyprime=zzzv3
         call DSVM2(sumyprimeoft+1,2,zzzvbar2,deloft+i,zzzsvm1,
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
      call DSVM3(deloft+nstate,zzzsvm1,zzzvbar6,pressureoft+1,3,
     $   zzzvbar4,temperatureoft+1,3,zzzvbar1,sumyprimeoft+1,2)
C     
C     Done
      end
