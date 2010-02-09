c*******9/14/98  added new output file that simply gives the 3 freq and degen
c                for each species in chemdis format
c*******6/8/95  modifeid output format for hf298 stuff
c*******5/19/95  added option to specify number of rotors and # oscillators (to better
c                fit Cpinf)
c********5/17/95  increased dimensions so that could use alternative input and fit larger
c                 number of heat capacity points             amd
c***********************************************************************
c                 format for "nonlst" input:
c    line 1:  cflag, spec(:16)   [a1,a16]   (flag and species name)
c    line 2:  (ie(i),ne(i), i=1,4)   [(4(a2,i3,1x)]  element and #)
c    line 3:  ref(1), ref(2), ref3  [(2a4,a6)]  (id)
c    line 4:  phase,nrot  [a1,i2]   (PHASE, NUMBER OF ROTORS, NEGATIVE IF LINEAR)
c    line 5:  hf298,s298,nt    (HF298, S298, NUMBER OF T,CP DATA--29 max)
c    line 6:  t(1),cpl(1)     (T, Cp)
c    line 6 + nt -1:    t(nt,cp(nt)
c    same sequence for species 2, etc. 
c    END       [start in column 1]    (end of data set)
c**************************************************************************
c**********5/9/95  combined some of the output files and added more fit info to them  amd
c *************also added flag to specify high T Cp input instead of computing kinf
c              now all output except nasa polynomials sent to lu 25 (log file)
c              polynoials sent to lu 54
      PROGRAM THERFIT 
c CRAY VERSION (ADDED THSB TO BOTTOM) 
c     works on SGI - fixed some errors in common blocks <ayc 1/94> 
c 
c     <2/8/92 ayc>  prints hoe fit parameters to lu 25
c                   prints exp fit parameters to lu 30
c                   log output file is lu 20
c                   nasa coefficients are written to lu 54
c                   input file is lu 1
c
c                   cflag = '@' supresses output of nasa coefficients
c                      put this before 1st character of species name
c                   cflag = '*' still skips all fits but prints out
c                         nasa coefficients
c
c                   For polyatomics (# atoms > 2):
c                      use NROT<0 (last 2 cols in input file - after 'G')
c                         to specify a linear molecule
c                      otherwise assumed non-linear with NROT int rots
c
c                   note no fits are performed on monotonics and only
c                   exponential fit option is available for diatomics
c
C********************************************************************** THE00040
C     THIS PROGRAM CREATES A NASA FORMAT THERMO DATABASE FILE           THE00050
C                                                                       THE00060
C       TWO SETS OF COEFFICIENTS ARE GENERATED AFTER FITTING CP DATA    THE00070
C       USING A MODIFIED VERSION OF THE CPFIT PROGRAM                   THE00080
C                                                                       THE00090
C          TMP DATA FILE (UNIT 45) IS USED TO PASS FILENAMES DIRECTLY   THE00100
C          TO THIS PROCEDURE FROM THERM.  IF THE FILE "THERMFIT.TMP"    THE00110
C          DOES NOT EXIST, THEN THE USER IS PROMPTED FOR INPUT/OUTPUT   THE00120
C          FILENAMES.                                                   THE00130
C                                                                       THE00140
C ********************************************************************* THE00150
c      2/08/93 ayc  fixed linear molecules/ added fit bypass option @;
c                   also added exp fit output and more error reporting
C      6/29/92 ayc  prints out hoe parameters and redirects screen ouput
c
C *****   5/25/90 MODIFIED TO RUN ON FLIC.  ALL OPEN AND CLOSES   ***** THE00160
C *****   WERE REPLACED BY FILEDEFS.  ROUTINES TO CAPITOLIZE      ***** THE00170
C *****   STRINGS WERE REMOVED, AS WELL AS ROUTINES TO            ***** THE00180
C *****   MANIPULATE FILENAMES.  THE SUBROUTINE NASAIN WAS        ***** THE00190
C *****   REMOVED AND IN ITS PLACE, THE OUTPUT FILE IS OPENED     ***** THE00200
C *****   IN APPEND MODE.  ALSO THE OUTPUT FILE HAS BEEN CHANGED  ***** THE00210
C *****   FROM A RANDOM ACCESS FILE TO A SEQUENTIAL FILE.  THE    ***** THE00220
C *****   EXEC TO RUN THIS PROGRAM IS CALLED  THERFIT EXEC.       ***** THE00230
C *****   THE INPUT FILE IS #1, THE OUTPUT FILE IS #54.    (JCB)  ***** THE00240
C ********************************************************************* THE00250
C                                                                       THE00260
C                                                                       THE00270
C      9/19/89 : modified to include Wilhoit polynomial                 THE00280
C                and Ritter exponential polynomial fitting methods      THE00290
C                as options.  In addition, if this program is run directTHE00300
C                from DOS, the user is given a menu to choose the fittinTHE00310
C                method.                                                THE00320
C                If HOE (default) fit fails (CPFIT failure) then the    THE00330
C                program defaults to the Wilhoit polynomial method      THE00340
C                to extrapolate Cp data and an appropriate warning messaTHE00350
C                is generated.                                          THE00360
C                A XCHECK procedure is included to allow entries in the THE00370
C                output file to be overwritten if they appear in the inpTHE00380
C                file.                                                  THE00390
C                Also, when run directly from DOS, an option is includedTHE00400
C                for non-LST format input.  In these cases the required THE00410
C                input format are as follows:                           THE00420
C********************************************************************** THE00730
      IMPLICIT REAL*8 (A-H,O-Z)                                         THE00740
        DIMENSION X0(2),WX(2),S(2)                                      THE00750
      DIMENSION CPL(30),TRA(3),A(5),B(6),C(5),F(14),FS(3,14),           THE00760
     +TH(3),TL(3),TM(3),SSRS(3),PERMX(3),A1(5)                          THE00770
      DIMENSION TFIT(30),CPFIT(30),PWORK(120),D(14),CONST(3)            THE00780
      DIMENSION B1(5),T(30),TEST(10),cphoe(30), cppoly(30),CPEST(30)                     

c     dmm 20000602

      include 'intdegen.fh'

      INTEGER Iccount
c     Iccount is to make sure READIN gets called only once ... this is
c     the XMG version ..
        REAL*8 TMID                                                     THE00810
        REAL*8 DELPCT(30)                                               THE00820
      CHARACTER*6 TITLE(5),HEADIN(21),REF3                              THE00830
      CHARACTER*1 PHASE,CFLAG,ANS                                       THE00840
      CHARACTER*2 IE(4)                                                 THE00850
      INTEGER NE(4),NE1(4,300),NX(300),ICNT,HOEFAI                      THE00860
      CHARACTER*4 INAME,REF(2),IFLAG,IE1(4,300)                         THE00870
        CHARACTER*16 FUNITS,NUL,SS(300),INPUT,SPEC                      THE00880
        CHARACTER*70 FILIN,FILOUT,COMMAN,DUMMY,S1,COMS                  THE00890
        LOGICAL EXST,OPNED,linear
        EXTERNAL FX                                                     THE00910
        COMMON/DATA/CPEX1,T1,IT,NT,NS                                   THE00920
        COMMON/CONST/CONST,NTERMS                                       THE00930
        COMMON/TOL/TOL,TSTEP                                            THE00940
        COMMON/TEMFIT/TFIT,F                                            THE00950
        COMMON/FLAG/IREPEA,IFLG, DELPCT                                 THE00960
        COMMON/NROTOR/NROT                                              THE00970
        COMMON/EXPAR/BEXPSV                                             THE00980
      DATA NPTS/7/                                                      THE00990
      DATA IFLAG/'END'/, REF3/'      '/                                 THE01000
      DATA IC1,IC2,IC3,IC4/1,2,3,4/                                     THE01010
        TEST(1)=300.                                                    THE01020
        TEST(2)=400.                                                    THE01030
        TEST(3)=500.                                                    THE01040
        TEST(4)=600.                                                    THE01050
        TEST(5)=700.                                                    THE01060
        TEST(6)=800.                                                    THE01070
        TEST(7)=900.                                                    THE01080
        TEST(8)=1000.                                                   THE01090
        TEST(9)=1100.                                                   THE01100
        TEST(10)=1200.                                                  THE01110
        NEWFIL=1                                                        THE01120
        T0=298.15                                                       THE01130
        TLO=300.0                                                       THE01140
        THI=5000.0                                                      THE01150
        TMID=1500.0                                                     THE01160
        TOL=1.0E-4                                                      THE01170
        TSTEP=0.0                                                       THE01180
3       FORMAT(A70)                                                     THE01190
14      FORMAT(A16)                                                     THE01200
        
        Iccount = 0
2       CONTINUE                                                        THE01210
        COMMAN=' '                                                      THE01220
        INPUT='LST'                                                     THE01230
         FILIN=' '                                                      THE01240
         FILOUT=' '                                                     THE01250
  311    CONTINUE                                                       THE01480
c        WRITE(*,*)'        CP FITTING / EXTRAPOLATION OPTIONS'          THE01490
c        WRITE(*,*)' '                                                   THE01500
c        WRITE(*,*)'                 1  HOE FIT'                         THE01510
c        WRITE(*,*)'                 2  EXP FIT'                         THE01520
c        WRITE(*,*)'                 3  WILHOI FIT'                      THE01530
c        WRITE(*,*)' '                                                   THE01540
c        WRITE(*,*)'              ENTER OPTION {1} (Q = QUIT)'           THE01550
c                READ(*,'(A)')ANS                                        THE01560
         ANS='1'
                IF(ANS.EQ.'Q'.OR.ANS.EQ.'Q')THEN                        THE01570
                        STOP 'THRFIT TERMINATED BY USER'                THE01580
                ELSEIF(ANS.EQ.'2')THEN                                  THE01590
                        COMMAN='EXP'                                    THE01600
                ELSEIF(ANS.EQ.'3')THEN                                  THE01610
                        COMMAN='WILHOI'                                 THE01620
                ELSE                                                    THE01630
                        COMMAN='HOE'                                    THE01640
                ENDIF                                                   THE01650
                COMS=COMMAN                                             THE01660

c***********************************************************************
c dmm 20000602 Integer fitting user option + open plot file
c***********************************************************************

c       write(*,*) '        INTEGER FREQUENCY DEGENERACY OPTIONS'
c       write(*,*) '                                            '
c       write(*,*) ' Enter {1} to force mode degeneracies to be integers'
c       write(*,*) ' Enter any other number for normal case -- noninteger
c     $ degeneracies'

c       read(*, '(A)')ANS
       ANS='1'
       if (ANS.eq.'1') then
          idg = 1
       else
          idg = 0
       endif   
       irec = 1
c       open(unit=21, file = 'cpplot.out', form = 'formatted', status
c     $      ='new', access =  'direct', RECL = 30)



c***********************************************************************
c dmm 20000602 end integer fitting user option
c***********************************************************************

                
C              CALL CLRSCR                                              THE01670
C                                                                       THE01680
C***** input format options menu                                        THE01690
C                                                                       THE01700
c        WRITE(*,*)'                     INPUT OPTIONS'                  THE01710
c        WRITE(*,*)' '                                                   THE01720
c        WRITE(*,*)'                 1 THERM "LST" FILE FORMAT'          THE01730
c        WRITE(*,*)'                 2 "nonlst" format  '                THE01750
c        WRITE(*,*)' '                                                   THE01760
c        WRITE(*,*)'              ENTER OPTION {1} (Q = QUIT)'           THE01770
c                READ(*,'(A)')ANS                                        THE01780
       ANS='1'
                IF(ANS.EQ.'Q')THEN                        
                        GOTO 311                                        THE01800
                 ELSEIF(ANS.EQ.'2')THEN                                  
                        INPUT='NONLST'                                  THE01840
                ELSE                                                    THE01850
                        INPUT='LST'                                     THE01860
                ENDIF                                                   THE01870
c                                                                     THE01940
c   write(*,*) ' Enter 0 to use default Cp(inf.) for all species'
c   write(*,*) ' Enter 1 to use specify  Cp(inf.) for some species'
c   read(*,*) iskip
c begin pey -- changed value of NEWFIL to make reading NASA fits from RMG easier
c      NEWFIL=1
    NEWFIL=0
c end pey
        iskip=0
4       CONTINUE                                                        THE02920
        IF(NEWFIL.EQ.1)THEN                                             THE02930
              WRITE(54,99)                                              THE02940
              WRITE(54,98)                                              THE02950
        ENDIF                                                           THE02960
   99 FORMAT('THERMO ',73(' '))                                         THE02970
   98 FORMAT('   300.000  1500.000  5000.000 ',49(' '))                 THE02980
    1 ICOUNT=1                                                          THE02990
        IF(INPUT.EQ.'LST')THEN                                          THE03000
c                READ(1,3)DUMMY                                          THE03010
c                CALL DBLNK(DUMMY,LEN,0)                                 THE03020
                                                                        THE03030
                FUNITS='KCAL'                                           THE03040
c                IF(DUMMY(:6).EQ.'UNITS:')THEN                           THE03050
c                        IF(DUMMY(7:8).EQ.'KJ')THEN                      THE03060
c                                FUNITS='KJ'                             THE03070
c                        ENDIF                                           THE03080
c                        READ(1,3)DUMMY                                  THE03090
c                ENDIF                                                   THE03100
   95           FORMAT(5A6/,21A6)                                       THE03110
c                READ(1,14)NUL                                           THE03120
        ENDIF                                                           THE03130
20      CONTINUE                                                        THE03140
        COMMAN=COMS                                                     THE03150
        DO 203 IK1=1,20                                                 THE03160
203     DELPCT(IK1)=0.0                                                 THE03170
        YERRS=1.0E+7                                                    THE03180
        IREPEA=0                                                        THE03190
        IFLG=0                                                          THE03200
        IPOLY=0                                                         THE03210
        IPOLYC=0                                                        THE03220
        CFLAG=' '                                                       THE03230
        BEXPSV=0.0                                                      THE03240
        TLO=300.0                                                       THE03250
        THI=5000.0                                                      THE03260
        TMID=1500.0                                                     THE03270

        if(Iccount.eq.0) then
c   pey added COMMAN to the following argument list
           CALL READIN(SPEC,INPUT,HF298,S298,CPL,T,REF,REF3,IE,NE,PHASE,
     +          NROT,NT,FUNITS,CFLAG,ICOUNT,COMMAN)

           Iccount = 1
        else
           INPUT = 'DONE'
        endif

C**IMPORTANT** we allow NROT to be negative to flag linear molecules
C              in list file input
C              we now retranslate flag and reset NROT   <ayc 2/93>
C              (we do this as soon as possible so that NROT can be set
C              to its correct value of zero before we get into trouble)

        if (nrot.lt.0) then
           linear = .true.
           nrot = 0
        else
           linear = .false.
        endif

        HOEFAI=0                                                        THE03310
c   write (*,7621) input
c 7621  format(1x,' input =',a14)
        IF(INPUT.EQ.'DONE')GOTO 900                                     THE03400
        IF(FUNITS.EQ.'KJ')THEN                                          THE03410
                CONST(1)=4.184                                          THE03420
        ELSE                                                            THE03430
                CONST(1)=1.0                                            THE03440
        ENDIF                   
c                                        
c     adding ichck initialization; this seems to be never used in this
C     version for XMG (20031105)
        ICHCK = 0
        IF(NEWFIL.EQ.0.AND.ICHCK.NE.0)THEN                              THE03460
                ISAME=0                                                 THE03470
                CALL XCHEX(SPEC,10,SS,NX,ICNT,IADD,ISAME)               THE03480
                IF(ISAME.EQ.1)THEN                                      THE03490
                        JOUT=NX(IADD)                                   THE03500
                ELSE                                                    THE03510
                        JOUT=NX(ICNT)+4                                 THE03520
                        ICNT=ICNT+1                                     THE03530
                        SS(ICNT)=SPEC                                   THE03540
                        NX(ICNT)=JOUT                                   THE03550
                ENDIF                                                   THE03560
        ELSEIF(NEWFIL.EQ.0.AND.ICHCK.EQ.0)THEN                          THE03570
                JOUT=3+(ICNT*4)                                         THE03580
                ICNT=ICNT+1                                             THE03590
        ELSEIF(NEWFIL.EQ.1)THEN                                         THE03600
                IF(ICOUNT.EQ.1)THEN                                     THE03610
                        JOUT=3                                          THE03620
                ELSE                                                    THE03630
                        JOUT=JOUT+4                                     THE03640
                ENDIF                                                   THE03650
        ENDIF                                                           THE03660
        ICOUNT=ICOUNT+1                                                 THE03670
        TRA(1)=300.                                                     THE03680
        NS=0                                                            THE03690
        NATOM=0                                                         THE03700
        ICARBO=0                                                        THE03710
        NCARBO=0                                                        THE03720
        DO 765 IIJ=1,4                                                  THE03730
        IF(IE(IIJ).EQ.'C '.OR.IE(IIJ).EQ.' C')THEN                      THE03740
                ICARBO=1                                                THE03750
                NCARBO=NE(IIJ)                                          THE03760
        ENDIF                                                           THE03770
765     NATOM=NATOM+NE(IIJ)                                             THE03780
986     CONTINUE                                                        THE03790
C       BYPASS FITTING ROUTINE IF ATOM OR Cp <> F(T)                    THE03800
        DIF1=DABS(CPL(1)-CPL(3))                                        THE03810
        DIF2=DABS(CPL(1)-CPL(6))                                        THE03820
        DIFF=MAX(DIF1,DIF2)                                             THE03830
        IF(FUNITS.EQ.'KJ')DIFF=DIFF/4.184                               THE03840
        IF(NATOM.NE.1.AND.DIFF.GT.0.01)THEN                             THE03850
                IF ((linear).OR.(NATOM.EQ.2)) then
                        NS=(3*NATOM)-5                                  THE03880
                        CONST(2)=FLOAT(NS)                              THE03890
                        CONST(3)=1.0                                    THE03900
                        CPINFI=1.987*(FLOAT(3*NATOM)-1.5)*CONST(1)      THE03910
                        CPZERO=(7./2.)*1.987*CONST(1)                   THE03920
                ELSE                                                    THE03930
                        NS=(3*NATOM)-6                                  THE03940
                        CONST(3)=0.0                                    THE03950
C               count hindered rotor as 1/2 degree of freedom           THE03960
                        CONST(2)=FLOAT(NS)-(FLOAT(NROT)/2.0)            THE03970
                CPINFI=1.987*(FLOAT(3*NATOM)-(2+(FLOAT(NROT)/2.0)))*    THE03980
     +CONST(1)                                                          THE03990
                CPZERO=4.*1.987*CONST(1)                                THE04000
                ENDIF                                                   THE04010
                IF(CFLAG.NE.'*')THEN                                    THE04030
                        IF(CPL(7).EQ.0.0.AND.INPUT.EQ.'LST')THEN        THE04040
                                NT=6                                    THE04050
                        ELSEIF(INPUT.EQ.'LST')THEN                      THE04060
                                NT=7                                    THE04070
                        ENDIF                                           THE04080
                ENDIF                                                   THE04090

c   if the molecule is diatomic, the following if-block overrides the 
c   user's choice of COMMAN and sets it to EXP. This was observed to
c   give very bad results in certain cases. Now the user's option
c   (RMG always uses WILHOI) is perserved. <pey 7/04>
c                                                                        THE04100
c        IF(NATOM.EQ.2)THEN                                              THE04110
c                NT1000=0                                                THE04120
c                NT1500=0                                                THE04130
c                NT2000=0                                                THE04140
c                DO 5465 IJ1=1,NT                                        THE04150
c                        IF(T(IJ1).EQ.1000.)NT1000=IJ1                   THE04160
c                        IF(T(IJ1).EQ.1500.)NT1500=IJ1                   THE04170
c                        IF(T(IJ1).EQ.2000.)NT2000=IJ1                   THE04180
c5465            CONTINUE                                                THE04190
c                COMMAN='EXP'                                            THE04200
c               CFLAG=' '                                               THE04210
c              we are going to comment out the next section
c              cause it makes no sense    <ayc 2/93>
c              OKAY i get it, its a kludge for when the data exceeds
c              the expected high T limit  <ayc 1/94>
c              The next section causes therfit to fail for bad Cp data
c              where Cp(1500) < Cp(1000), which can happen in RMG <pey 7/04>
c                IF(NT1000.NE.0.AND.NT1500.NE.0)THEN                     THE04220
c                        TERM=(CPL(NT1500)-CPL(NT1000))/0.37             THE04230
c                        CPINFI=CPL(NT1000)+TERM                         THE04240
c                        IF(NT2000.EQ.0)THEN                             THE04250
c                                T(NT+1)=2000.                           THE04260
c                                TERM=(CPINFI-CPL(NT))*0.37              THE04270
c                                CPL(NT+1)=CPL(NT)+TERM                  THE04280
c                                NT=NT+1                                 THE04290
c                        ENDIF                                           THE04300
c                        CONST(2)=(CPINFI/(1.987*CONST(1)))-3.5          THE04310
c                ENDIF                                                   THE04320
c        ENDIF                                                           THE04330
                                                                        THE04340
                                                                        THE04350
                IF(CFLAG.EQ.'*')THEN                                    THE04360
                        IF(CPL(7).EQ.0.0.AND.INPUT.EQ.'LST')THEN        THE04370
                                NT=6                                    THE04380
                        ELSEIF(INPUT.EQ.'LST')THEN                      THE04390
                                NT=7                                    THE04400
                        ENDIF                                           THE04410
                DO 69696 II=1,NT                                        THE04420
                IF(CPL(II).GT.CPINFI)THEN                               THE04430
                        CPINFI=1.1*CPL(II)                              THE04440
                ENDIF                                                   THE04450
69696           CONTINUE                                                THE04460
                IF(T(NT).EQ.1000.)THEN                                  THE04470
                        NT=NT+1                                         THE04480
                        T(NT)=1500.                                     THE04490
                        CPL(NT)=CPL(NT-1)*(0.37*(CPINFI-CPL(NT-1)))     THE04500
                ELSEIF(T(NT).EQ.1500.)THEN                              THE04510
                        NT=NT+1                                         THE04520
                        T(NT)=2000.                                     THE04530
                        CPL(NT)=CPL(NT-1)*(0.37*(CPINFI-CPL(NT-1)))     THE04540
                ENDIF                                                   THE04550
                DO 329 II=1,NT                                          THE04560
329             CPL(II)=CPL(II)/(1.987*CONST(1))                        THE04570
                CALL LLSQ(T,CPL,C,5,NT,SSR,AVGERR)                      THE04580
C                                                                       THE04590
                DO 6969 II=1,5                                          THE04600
                        F(II)=C(II)                                     THE04610
6969                    F(II+7)=C(II)                                   THE04620
                        F(6)=F6(0,(HF298/1.987),T0,F)                   THE04630
      F(7)=(S298/1.987)-F(1)*DLOG(T0)-F(2)*T0-F(3)/2.*T0**2             THE04640
     1-F(4)/3.*T0**3-F(5)/4.*T0**4                                      THE04650
                        F(13)=F(6)                                      THE04660
                        F(14)=F(7)                                      THE04670
                        TMID=1000.                                      THE04680
                        THI=2000.                                       THE04690
                        GOTO 12987                                      THE04700
                ENDIF                                                   THE04710
    if (iskip.eq.0) then 
c       write(*,*) '  Using all defaults'
         go to 1313
        else
        endif
    icp=0
    write(*,4321) spec(:16)
 4321   format(1x,'next species is: ',a16)
    write(*,*)'enter 0 to compute default Cp(Tinf) or 1 to specify'
    read (*,*) icp
    if (icp.eq.1) then
        write(*,6521) SPEC(:16) 
 6521   format(1x, 'enter T(high), Cp (hIgh) for: ',a16)
        read(*,*) thigh, cpinfi
c********option to specify number of vib. and # rot
        write(*,*) ' enter 0 for default vib or 1 to specify'
        read(*,*) ivib
        if(ivib.eq.1) then
        write(*,*) 'enter # vib. and # int. rot '
        read(*,*) ns,nrot
        const(2)=ns-0.5*nrot
        endif
    else
        write(*,*)  ' Using default Cp inf and default vib'
    endif                                                                        
C-----------------------------------------------------------------------THE04730
C                                                                       THE04740
C          EXTRAPOLATE/FIT INPUT CP DATA                                THE04750
C                                                                       THE04760
C-----------------------------------------------------------------------THE04770
1313            CONTINUE                                                THE04780
989     CONTINUE                                                        THE04790

c             calculate fits and list to file

c             Wilhoit polynomial fit
                IF (COMMAN.EQ.'WILHOI')THEN                             THE04810
                   CALL WILHOI(NT,CPL,T,CPZERO,CPINFI,A1,B2)
                   write(25,*)
                   write(25,2500) spec(:16)
2500               FORMAT('--> ',A16,' ---------',
     2                '--- WILHOIT FIT ------------')
                   WRITE(25,2505) (A1(ifit),ifit=1,4)
2505               FORMAT(' a0 = ',1pe10.3,' a1 = ',
     2                1pe10.3,' a2 = ',1pe10.3,' a3 = ', 1pe10.3)
                   WRITE(25,2510) (B2)
2510               FORMAT(' B = ',1pe10.3)
                   write(25,*)
                   IF (CONST(3).EQ.1)  THEN
                      WRITE(25,2580) NS
                   ELSE 
                      WRITE(25,2590) NS,NROT
                   ENDIF
                   write(25,2596) cpzero,cpinfi
c               EXP fit
                ELSEIF(COMMAN.EQ.'EXP')THEN
                   CALL CPEXFT(NT,CPL,T,CPZERO,CPINFI,A1)
                   write(25,*)
                   write(25,2520) spec(:16)
2520               FORMAT('--> ',A16,' ---------',
     2                '--- EXPONENTIAL FIT ------------')
                   write(25,2525)
2525               format('    Cp = Cp(0) + [Cp(inf)-Cp(0)] ',
     2                '[1 - exp(- K)]')
                   WRITE(25,2530) (a1(ifit),ifit=1,4)
2530               FORMAT('    K = ',1pe10.3,' + ',
     2                1pe10.3,' T + ',1pe10.3,' T^2 + ',
     3                1pe10.3,' T^3')
                   write(25,*)
                   IF (CONST(3).EQ.1)  THEN
                      WRITE(25,2580) NS
                   ELSE 
                      WRITE(25,2590) NS,NROT
                   ENDIF
                   write(25,2596) cpzero,cpinfi

c               HOE fit
                ELSE
c     dmm 20000531 commented out THCPFI, replaced with integer-degen capable IDCPFIT

c                   CALL THCPFI(NT,CPL,NATOM,B,T,CPINFI,icp,thigh)
                   call idcpfit(NT, CPL, NATOM, B, T, CPINFI, icp,thigh)

c                   call printfit(NT, CPL, B, T, CPINFI, icp, thigh, SPEC
c     $                  ,irec)

c     dmm 20000531 end modifications

                   WRITE(25,*)
                   WRITE(25,2540) SPEC(:16)
2540               FORMAT('--> ',A16,' ---------',
     2                ' HARMONIC OSCIL FIT ------------')
    if(ivib.eq.1) then
    write(25,*)'# vib. and # int. rot. specified at run time'
    endif    
                   DO 2560 IVIBE=1,3
                      WRITE(25,2550) IVIBE,B(IVIBE),B(IVIBE+3)
2550                  FORMAT('       VIBRATIONS #',I1,':  MODES = ',
     2                   F6.1,'  FREQUENCY = ',F7.1,' CM-1')
2560               CONTINUE
                   VTOT = B(1)+B(2)+B(3)
                   FREQM = B(4)**(B(1)/VTOT)*B(5)**(B(2)/VTOT)
     2                *B(6)**(B(3)/VTOT)
                   WRITE(25,2570) VTOT,FREQM
2570               FORMAT(' GEO MEAN VIBRATION:  MODES = ',F6.3,
     2                '  FREQUENCY = ',F7.1,' CM-1')
                   write(25,*)
                   IF (CONST(3).EQ.1)  THEN
                      WRITE(25,2580) NS, nrot
2580                  FORMAT(' CALCULATIONS FOR A LINEAR',
     2                   ' MOLECULE WITH ',I3,' INT MODES (', I2,
     3                   ' INT ROTS)')
                   ELSE 
                      WRITE(25,2590) NS,NROT
2590                  FORMAT(' CALCULATIONS FOR A NON-LINEAR',
     2                   ' MOLECULE WITH ',I3,' INT MODES (', I2,
     3                   ' INT ROTS)')
                   ENDIF
c*************************9/14/98***************************************
c  added printout of only frequency info in chemdis format
    write(26,2640)  spec(:16)
 2640 format(1x, a16)
      write(26,2650) b(4),b(1), b(5), b(2), b(6), b(3)
 2650 format(1x,'3', 3(f7.1,2x,f6.3))          
c*************************9/14/98***************************************

        if (icp.eq.1) then
                   write(25,2595) cpzero,thigh, cpinfi
2595               format(' (IN CHOSEN UNITS:    Cp(0): ',f6.3,
     2                '   Cp(',f6.0,'): ',f7.3,')')
        else
                   write(25,2596) cpzero,cpinfi
2596               format(' (IN CHOSEN UNITS:    Cp(0): ',f6.3,
     2                '   Cp(inf): ',f7.3,')')
        endif
                   IF(IFLG.EQ.1)THEN
                      COMMAN='EXP'
                      write(25,*)
                      write(25,*) '*****WARNING: GOODNESS-OF-FIT',
     2                   ' CRITERIA NOT SATISFIED*****'
                      GOTO 989
                   ENDIF
                ENDIF                                                   THE04910
                                                                        THE04920
1994            CONTINUE                                                THE04940
                                                                        THE04950
        ENDIF

c      added new cflag option - surpress polynomial fit
        if (cflag.eq.'@') goto 2220

C******************************************************************     THE04960
73590   CONTINUE                                                        THE04980
        DO 431 JK=1,2                                                   THE04990
        BEXPSV=0.0                                                      THE05000
        IF(NATOM.NE.1.AND.DIFF.GT.0.01)THEN                             THE05010
C------------------------------------------------------------------     THE05020
C      SECOND POLYNOMIAL FIT / SECOND OVERLAP REGION                    THE05030
C------------------------------------------------------------------     THE05040
        IF(IPOLY.EQ.1)THEN                                              THE05050
                IF(JK.EQ.1)THEN                                         THE05060
                TK=253.3333                                             THE05070
                DO 1290 I8=1,15                                         THE05080
                        TK=TK+46.6667                                   THE05090
                        TFIT(I8)=TK                                     THE05100
                        TFIT(I8+15)=1000.+(71.4236*FLOAT(I8-1))         THE05110
1290            CONTINUE                                                THE05120
C       set lower temperature below room temp so that 298 is not        THE05130
C        outside the range of the regression for Cp                     THE05140
                TFIT(1)=273.15                                          THE05150
                ELSE                                                    THE05160
                TK=928.5764                                             THE05170
                DO 1291 I8=1,15                                         THE05180
                        TFIT(I8)=1000.+(71.4236*FLOAT(I8-1))            THE05190
                        TFIT(I8+15)=2000.+(266.667*FLOAT(I8))           THE05200
1291            CONTINUE                                                THE05210
                ENDIF                                                   THE05220
C------------------------------------------------------------------     THE05230
C      FIRST POLYNOMIAL FIT / FIRST OVERLAP REGION                      THE05240
C------------------------------------------------------------------     THE05250
                                                                        THE05260
        ELSEIF(IPOLY.EQ.0)THEN                                          THE05270
                IF(JK.EQ.1)THEN                                         THE05280
                TK=236.6667                                             THE05290
                DO 21290 I8=1,15                                        THE05300
                        TK=TK+63.3333                                   THE05310
                        TFIT(I8)=TK                                     THE05320
                        TFIT(I8+15)=1250.+(35.714*FLOAT(I8-1))          THE05330
21290           CONTINUE                                                THE05340
C       set lower temperature below room temp so that 298 is not        THE05350
C        outside the range of the regression for Cp                     THE05360
                TFIT(1)=273.15                                          THE05370
                ELSE                                                    THE05380
                DO 21291 I8=1,15                                        THE05390
                        TFIT(I8)=1250.+(35.714*FLOAT(I8-1))             THE05400
                        TFIT(I8+15)=1700.+(283.333*FLOAT(I8))           THE05410
21291           CONTINUE                                                THE05420
                ENDIF                                                   THE05430
C------------------------------------------------------------------     THE05440
C       FOURTH POLYNOMIAL FIT/ FOURTH OVERLAP REGION                    THE05450
C------------------------------------------------------------------     THE05460
                                                                        THE05470
        ELSEIF(IPOLY.EQ.3)THEN                                          THE05480
                IF(JK.EQ.1)THEN                                         THE05490
                TK=263.333                                              THE05500
                DO 22290 I8=1,15                                        THE05510
                        TK=TK+36.667                                    THE05520
                        TFIT(I8)=TK                                     THE05530
                        TFIT(I8+15)=850.+(35.714*FLOAT(I8-1))           THE05540
22290           CONTINUE                                                THE05550
C       set lower temperature below room temp so that 298 is not        THE05560
C        outside the range of the regression for Cp                     THE05570
                TFIT(1)=273.15                                          THE05580
                ELSE                                                    THE05590
                DO 22291 I8=1,15                                        THE05600
                        TFIT(I8)=850.+(35.714*FLOAT(I8-1))              THE05610
                        TFIT(I8+15)=1350.+(310.0*FLOAT(I8))             THE05620
22291           CONTINUE                                                THE05630
                ENDIF                                                   THE05640
C------------------------------------------------------------------     THE05650
C       THIRD POLYNOMIAL FIT/ THIRD OVERLAP REGION                      THE05660
C------------------------------------------------------------------     THE05670
                                                                        THE05680
        ELSEIF(IPOLY.EQ.2)THEN                                          THE05690
                IF(JK.EQ.1)THEN                                         THE05700
                TK=263.333                                              THE05710
                DO 23290 I8=1,15                                        THE05720
                        TK=TK+36.667                                    THE05730
                        TFIT(I8)=TK                                     THE05740
                        TFIT(I8+15)=850.+(21.429*FLOAT(I8-1))           THE05750
23290           CONTINUE                                                THE05760
C       set lower temperature below room temp so that 298 is not        THE05770
C        outside the range of the regression for Cp                     THE05780
                TFIT(1)=273.15                                          THE05790
                ELSE                                                    THE05800
                DO 23291 I8=1,15                                        THE05810
                                                                        THE05820
                        TFIT(I8)=850.+(21.429*FLOAT(I8-1))              THE05830
                        TFIT(I8+15)=1150.+(323.333*FLOAT(I8))           THE05840
23291           CONTINUE                                                THE05850
                ENDIF                                                   THE05860
                                                                        THE05870
                                                                        THE05880
        ENDIF                                                           THE05890
C-------------------------------------------------------------------    THE05900
C  OVERLAP REGION IS DETERMINED:  BEGIN REGRESSION FOR PAIR             THE05910
C     OF POLYNOMIALS                                                    THE05920
*-------------------------------------------------------------------    THE05930
                                                                        THE05940
      TRA(1)=300.                                                       THE05950
      TRA(2)=1000.                                                      THE05960
      DO 30 I=1,30                                                      THE05970
        TIN=TFIT(I)                                                     THE05980
        IF(COMMAN.EQ.'WILHOI'.OR.(COMMAN.EQ.'HOE'.AND.                  THE05990
     +HOEFAI.EQ.1))THEN                                                 THE06000
        CPFIT(I)=CPWH(A1,B2,TIN,CPZERO,CPINFI)/(1.987*CONST(1))         THE06010
                                                                        THE06020
        ELSEIF(COMMAN.EQ.'EXP')THEN                                     THE06030
                CPFIT(I)=CPEXP(A1,TIN,CPZERO,CPINFI)/(1.987*CONST(1))   THE06040
C       WRITE(25,'(1X,I2,2X,2(F12.3))')I,TIN,CPFIT(I)*1.987*CONST(1)    THE06050
        ELSE                                                            THE06060
                CPFIT(I)=CPSERI(B,TIN)                                  THE06070
        ENDIF                                                           THE06080
C        WRITE(25,'(1X,I2,2X,2(F12.3))')I,TIN,CPFIT(I)*1.987*CONST(1)    THE06090
   30 CONTINUE                                                          THE06100
      MDEG=4                                                            THE06110
      RSQ=100.                                                          THE06120
      IER=0                                                             THE06130
        NPTS=30                                                         THE06140
              CALL LLSQ(TFIT,CPFIT,C,5,NPTS,SSR,AVGERR)                 THE06150
                                                                        THE06160
      ELSE                                                              THE06170
C               ATOM:  USE DIFFERENT METHOD FOR Cp                      THE06180
C    BYPASS FITTING ROUTINE IF CP <> f(T)                               THE06190
                SUMCP=0.0                                               THE06200
                DATPNT=0.0                                              THE06210
                DO 1209 IJK=1,7                                         THE06220
                IF(CPL(IJK).NE.0.0)THEN                                 THE06230
                        SUMCP=SUMCP+CPL(IJK)                            THE06240
                        DATPNT=DATPNT+1.0                               THE06250
                ENDIF                                                   THE06260
1209            CONTINUE                                                THE06270
              C(1)=SUMCP/(1.98717*DATPNT)                               THE06280
                IF(FUNITS.EQ.'KJ')C(1)=C(1)/4.184                       THE06290
              DO 15 I=2,5                                               THE06300
              C(I)=0.0                                                  THE06310
   15         CONTINUE                                                  THE06320
                TMID=1000.                                              THE06330
      ENDIF                                                             THE06340
        IF(JK.EQ.1)THEN                                                 THE06350
                DO 75 I=1,5                                             THE06360
                F(I)=C(I)                                               THE06370
                                                                        THE06380
   75           CONTINUE                                                THE06390
        ELSEIF(JK.EQ.2)THEN                                             THE06400
                DO 85 J=1,5                                             THE06410
85              F(J+7)=C(J)                                             THE06420
        ENDIF                                                           THE06430
431     CONTINUE                                                        THE06440
C      CALCULATE CONST C6 IN ENTHALPY EXPRESSION                        THE06450
                                                                        THE06460
C     CALC OF Hf & S SO THAT 1st DERIVATIVE AND                         THE06470
C     FUNCTION VALUES MATCH AT T BREAK POINT                            THE06480
C                                                                       THE06490
C    TO DO THIS WE HOLD (C1-C5, C8-C12) CONSTANT AND SOLVE FOR C6,C13   THE06500
C    SUCH THAT Hfl(Tmid)=Hfh(Tmid)                                      THE06510
C      MUST INTEGRATE CP BETWEEN 298 AND Tmid                           THE06520
C                                                                       THE06530
C    WE MUST PERFORM A SIMILAR DETERMINATION FOR S TO GET C7,C14        THE06540
C                                                                       THE06550
C                                                                       THE06560
C                                                                       THE06570
        C6=F6(0,(HF298/1.987),T0,F)                                     THE06580
C                                                                       THE06590
C      CALCULATE CONST C7 FROM ENTROPY EXPRESSION                       THE06600
C                                                                       THE06610
      C7=(S298/1.987)-F(1)*DLOG(T0)-F(2)*T0-F(3)/2.*T0**2               THE06620
     1-F(4)/3.*T0**3-F(5)/4.*T0**4                                      THE06630
                                                                        THE06640
C       USE F(6)  TO DETERMINE HF(Tmid)                                 THE06650
        F(6)=C6                                                         THE06660
        F(7)=C7                                                         THE06670
        H0=HF(0,T0,F)*1.987*T0                                          THE06680
        S0=ST(0,T0,F)*1.987                                             THE06690
                                                                        THE06700
C       DETERMINE BREAK POINT TEMPERATURE                               THE06710
        IF(NATOM.NE.1.AND.CFLAG.NE.'*')THEN                             THE06720
                TINIT=(TFIT(1)+TFIT(15))/2.0                            THE06730
                CALL TBREAK (TMID,TINIT)                                THE06740
                TMID=FLOAT(INT(TMID+0.5))                               THE06750
        ELSE                                                            THE06760
                TMID=1000.                                              THE06770
        ENDIF                                                           THE06780
C*      THESE MUST BE DETERMINED TO SATISFY H(Tmid)l = H(Tmid)h.        THE06790
        HFTMID=HF(0,TMID,F)*1.987*TMID                                  THE06800
        STMID=ST(0,TMID,F)*1.987                                        THE06810
        F(13)=F6(7,(HFTMID/1.987),TMID,F)                               THE06820
                                                                        THE06830
      C14=(STMID/1.987)-F(8)*DLOG(TMID)-F(9)*TMID-F(10)/2.*TMID**2      THE06840
     1-F(11)/3.*TMID**3-F(12)/4.*TMID**4                                THE06850
        F(14)=C14                                                       THE06860
        HCKTM=HF(7,TMID,F)*1.987*TMID                                   THE06870
        SCKTM=ST(7,TMID,F)*1.987                                        THE06880
        DELTAC=100.*((CPT(7,TMID,F)-CPT(0,TMID,F))/CPT(0,TMID,F))       THE06890
        DDELCP=DCPT(7,TMID,F)-DCPT(0,TMID,F)                            THE06900
        DELTAS=100.*((ST(7,TMID,F)-ST(0,TMID,F))/ST(0,TMID,F))          THE06910
        DDELS=DST(7,TMID,F)-DST(0,TMID,F)                               THE06920
        DELTAH=100.*((HF(7,TMID,F)-HF(0,TMID,F))/HF(0,TMID,F))          THE06930
        DDELH=DHF(7,TMID,F)-DHF(0,TMID,F)                               THE06940
        YERR=0.0                                                        THE06950
        IF(ABS(DELTAH).GT.(50.*TOL).OR.ABS(DELTAS).GT.(50.*TOL)         THE06960
     +.OR.ABS(DDELH).GT.(20.*TOL).OR.ABS(DDELS).GT.(20.*TOL).OR.        THE06970
     +ABS(DELTAC).GT.(50.*TOL).OR.ABS(DDELCP).GT.(20.*TOL))THEN         THE06980
        IF(ABS(DELTAC).GT.(50.*TOL))THEN                                THE06990
                YERR=ABS(DELTAC)                                        THE07000
                WRITE(25,1960)ABS(DELTAC)                               THE07010
1960    FORMAT(1X,'TMID MISMATCH: % DELTA CP(TMID) = ',1PE10.3)         THE07020
        ENDIF                                                           THE07030
        IF(ABS(DDELCP).GT.(25.*TOL))THEN                                THE07040
                IF(ABS(DDELCP).GT.YERR)YERR=ABS(DDELCP)                 THE07050
                WRITE(25,1961)ABS(DDELCP)                               THE07060
1961    FORMAT(1X,'TMID MISMATCH:  DELTA DCP(TMID) = ',1PE10.3)         THE07070
        ENDIF                                                           THE07080
        IF(ABS(DELTAH).GT.(50.*TOL))THEN                                THE07090
                IF(ABS(DELTAH).GT.YERR)YERR=ABS(DELTAH)                 THE07100
                WRITE(25,1962)ABS(DELTAH)                               THE07110
1962    FORMAT(1X,'TMID MISMATCH: % DELTA  HF(TMID) = ',1PE10.3)        THE07120
        ENDIF                                                           THE07130
        IF(ABS(DDELH).GT.(25.*TOL))THEN                                 THE07140
                IF(ABS(DDELH).GT.YERR)YERR=ABS(DDELH)                   THE07150
                WRITE(25,1963)ABS(DDELH)                                THE07160
1963    FORMAT(1X,'TMID MISMATCH:  DELTA DHF(TMID) = ',1PE10.3)         THE07170
        ENDIF                                                           THE07180
        IF(ABS(DELTAS).GT.(50.*TOL))THEN                                THE07190
                IF(ABS(DELTAS).GT.YERR)YERR=ABS(DELTAS)                 THE07200
                WRITE(25,1964)ABS(DELTAS)                               THE07210
1964    FORMAT(1X,'TMID MISMATCH: % DELTA S(TMID) = ',1PE10.3)          THE07220
        ENDIF                                                           THE07230
        IF(ABS(DDELS).GT.(25.*TOL))THEN                                 THE07240
                IF(ABS(DDELS).GT.YERR)YERR=ABS(DDELS)                   THE07250
                WRITE(25,1965)ABS(DDELS)                                THE07260
1965    FORMAT(1X,'TMID MISMATCH:  DELTA DS(TMID) = ',1PE10.3)          THE07270
        ENDIF                                                           THE07280
        IF(IPOLY.EQ.0.AND.IPOLYC.EQ.0)THEN                              THE07290
                WRITE(25,*)' REPEATING POLYNOMIAL FITTING ROUTINE'      THE07300
                WRITE(25,*)' WITH A DIFFERENT OVERLAP REGION'           THE07310
                YERRS=YERR                                              THE07320
                YERR1=YERR                                              THE07330
                IPOLYC=IPOLYC+1                                         THE07340
                IPOLY=1                                                 THE07350
                SSRS(1)=SSR                                             THE07360
                TM(1)=TMID                                              THE07370
                DO 1357 IK1=1,14                                        THE07380
1357            FS(1,IK1)=F(IK1)                                        THE07390
                GOTO 73590                                              THE07400
        ELSEIF(IPOLY.EQ.1)THEN                                          THE07410
                WRITE(25,*)' SECOND REGRESSION FAILED TO MATCH'         THE07420
                WRITE(25,*)' PROPERTY AT THE MID-POINT'                 THE07430
        WRITE(25,*)' MAKING A THIRD ATTEMPT TO FIT THE DATA'            THE07440
                DO 13571 IK1=1,14                                       THE07450
13571           FS(2,IK1)=F(IK1)                                        THE07460
                IPOLY=2                                                 THE07470
                YERR2=YERR                                              THE07480
                SSRS(2)=SSR                                             THE07490
                TM(2)=TMID                                              THE07500
                GOTO 73590                                              THE07510
        ELSEIF(IPOLY.EQ.2)THEN                                          THE07520
                WRITE(25,*)' THIRD ATTEMPT TO MATCH PROPERTIES '        THE07530
                WRITE(25,*)' AT THE MIDPOINT HAS FAILED'                THE07540
                WRITE(25,*)' LAST REGRESSION IN PROGRESS'               THE07550
                IPOLY=3                                                 THE07560
                SSRS(3)=SSR                                             THE07570
                YERR3=YERR                                              THE07580
                TM(3)=TMID                                              THE07590
                DO 13572 IK1=1,14                                       THE07600
13572           FS(3,IK1)=F(IK1)                                        THE07610
                GOTO 73590                                              THE07620
        ELSEIF(IPOLY.EQ.3)THEN                                          THE07630
                IF(YERR1.LT.YERR2)THEN                                  THE07640
                        IF(YERR1.LT.YERR3)THEN                          THE07650
                                KEEP=1                                  THE07660
                        ELSE                                            THE07670
                                KEEP=3                                  THE07680
                        ENDIF                                           THE07690
                ELSE                                                    THE07700
                        IF(YERR2.LT.YERR3)THEN                          THE07710
                                KEEP=2                                  THE07720
                        ELSE                                            THE07730
                                KEEP=3                                  THE07740
                        ENDIF                                           THE07750
                ENDIF                                                   THE07760
                WRITE(25,*)' BEST REGRESSION IS # ',KEEP                THE07770
                DO 24 IK2=1,14                                          THE07780
24              F(IK2)=FS(KEEP,IK2)                                     THE07790
                TMID=TM(KEEP)                                           THE07800
        ENDIF                                                           THE07810
        ENDIF                                                           THE07820
C                                                                       THE07830
C       F ( 1-7 )  : LOW TEMPERATURE COEFFICIENTS                       THE07840
C       F ( 8-14 ) : HIGH TEMPERATURE COEFFICIENTS                      THE07850
C                                                                       THE07860
C      REMEMBER:  NASA FORMAT  FIRST 7 COEFFICIENTS (HIGH TEMP)         THE07870
C                              SECOND 7 COEFFICIENTS (LOW TEMP)         THE07880
C                                                                       THE07890
12987   CONTINUE                                                        THE07900
c   begin pey (18may04) change to chemkin format
c      WRITE(54,1010) SPEC(:16),(REF(I),I=1,2),REF3,(IE(I),              THE07910
    WRITE(54,1010) SPEC(:16),'    ',REF3,(IE(I),

c     1NE(I),I=1,4),PHASE,TLO,THI,TMID,NROT,IC1                          THE07920
     1NE(I),I=1,4),PHASE,TLO,THI,TMID,IC1                               THE07920
      WRITE(54,112) F(8),F(9),F(10),F(11),F(12),IC2                     THE07930
      WRITE(54,112)F(13),F(14),F(1),F(2),F(3),IC3                       THE07940
      WRITE(54,113)F(4),F(5),F(6),F(7),IC4                              THE07950
c 1010 FORMAT(A10,2A4,A6,4(A2,I3),A1,2F10.3,F9.3,3X,I2,I1)               THE07960
 1010 FORMAT(A16,A4,A4,4(A2,I3),A1,2F10.3,F9.3,5X,I1)                   THE07960
c   end pey
  112 FORMAT(5(1PE15.8),I5)                                             THE07970
  113 FORMAT(4(1PE15.8),15X,I5)                                         THE07980
c
c *** comparing fits to input values
c
    write(25,2594)
 2594   format(1x,/,'  HF298 input  HF298 fit  S298 input  S298 fit')  
C---   WANT TO COMPUTE H AND S AT 298 WITH THE LOW TEMP PARAMETERS      
      Tf=298.                                                            
      H=1.98717*Tf*((((F(5)*Tf/5.+F(4)/4.)*Tf+F(3)/3.)*Tf
     * +F(2)/2.)*Tf+F(1)+F(6)/Tf)                                                   
      Sfit=1.98717*((((F(5)*Tf/4.+F(4)/3.)*Tf+F(3)/2.)*Tf 
     *  +F(2))*Tf+F(1)*DLOG(Tf)+F(7))                                             
    write(25, 2593) HF298,H,S298,SFIT
 2593   format(1x,f9.0,4x,f9.0,3x,f7.2,4x,f7.2)
 2220   write (25,2597)
 2597 format(1x,/,'  T(K)',5x,'Cp input',3x,'Cp (HOE)',3x,'Delta',
     * 3x,'Cp (poly)',3x,'Delta')
c
c *** hoe fit
c
c *** if @, no polynomial fit, so skip next section
    if(cflag.eq.'@') go to 2498
c *** if * or diatomic, no harmonic oscillator fit, so only do poly fit
    if(cflag.eq.'*'.or.natom.eq.2.or.(.not.(COMMAN.eq.'HOE'))) then
        do 2699 it=1,nt
                IF (t(it).LE.tmid) THEN
          CPpoly(it)=1.98717*((((F(5)*t(it)+F(4))*t(it)+F(3))*t(it)
     *      +F(2))*t(it)+F(1))         
            else if (t(it).gt.thi) THEN
            go to 2498
                ELSE
          CPpoly(it)=1.98717*((((F(12)*t(it)+F(11))*t(it)+F(10))*t(it)
     *      +F(9))*t(it)+F(8))        
                ENDIF
        dif2=cppoly(it) - cpl(it)
        write(25,2698) t(it), cpl(it), cppoly(it), dif2  
 2698       format(1x,f7.0,3x,f7.2,20x,  2(3x,f7.2))
 2699       continue
c***  end of only poly fit, jump over h.o.e. fit
        go to 4432 
        endif
 2498       continue
c***  start h.o.e. fit comparison
c*** if * or diatomic, no harmonic oscillator fit, so skip
    if(cflag.eq.'*'.or.natom.eq.2) go to 4432
        do 2599 it=1,nt
        cphoe(it)=1.987*cpseri(b,t(it))
                dif1=cphoe(it)-cpl(it)
c***  no poly fit if @, so write h.o.e. comparison only
        if(cflag.eq.'@') then
            write(25,2298) t(it), cpl(it), cphoe(it), dif1
 2298           format(1x,f7.0,3(3x,f7.2))
            go to 2599
             elseIF (t(it).LE.tmid) THEN
          CPpoly(it)=1.98717*((((F(5)*t(it)+F(4))*t(it)+F(3))*t(it)
     *      +F(2))*t(it)+F(1))         
        else if (t(it).gt.thi) THEN
        go to 2799
              ELSE
          CPpoly(it)=1.98717*((((F(12)*t(it)+F(11))*t(it)+F(10))*t(it)
     *      +F(9))*t(it)+F(8))        
              ENDIF
    dif2=cppoly(it) - cpl(it)
        write(25,2598) t(it), cpl(it), cphoe(it), dif1,
     *  cppoly(it), dif2  
 2598   format(1x,f7.0, 5(3x,f7.2))
 2599   continue
 2799   continue
c*****getting high temp values
    if(thi.le.t(nt)) go to 4432
    hihoe=1.987*cpseri(b,thi)
        if(cflag.eq.'@') then
            write(25,2198) thi,hihoe
 2198           format(1x,f7.0,13x,f7.2)
        go to 4432
        endif   
        hipoly=1.98717*((((F(12)*thi  +F(11))*thi  +F(10))*thi  
     *       +F(9))*thi  +F(8))
        write(25,2398) thi, hihoe, hipoly 
 2398   format(1x,f7.0,13x,f7.2,13x,f7.2)
 4432   ICOUNT=ICOUNT+1                                                   
      GO TO 20                                                          THE08430
  900 CONTINUE                                                          THE08440
c   begin pey -- don't print END at end of thermo file
c 1012         FORMAT('END',77(' '))                                     THE08450
c      IF(NEWFIL.EQ.1)THEN                                               THE08460
c                JLAST=JOUT+4                                            THE08470
c                WRITE(54,1012)                                          THE08480
c        ELSEIF(NEWFIL.EQ.0.AND.ICHCK.EQ.1)THEN                          THE08490
c                JLAST=NX(ICNT)+4                                        THE08500
c                WRITE(54,1012)                                          THE08510
c        ELSE                                                            THE08520
c                WRITE(54,1012)                                          THE08530
c        ENDIF                                                           THE08540
c   end pey
C        CLOSE(1,STATUS='KEEP')                                         THE08550
C        CLOSE(54,STATUS='KEEP')                                        THE08560
      STOP ' *** Therfit Job Complete ***'
      END                                                              
C ***
C ***
C ***
C ***                                                                   THE08590
      SUBROUTINE THCPFI(NT,CP,NATOM,AI,TEMPIN,CPINF,icp,thigh)                  
                                                                        THE08610
C                                                                       THE08620
      IMPLICIT REAL*8(A-H,O-Z)                                          THE08630
      CHARACTER*8  FMT1                                                 THE08640
      CHARACTER*39 FMT2                                                 THE08650
      CHARACTER*33 FMT3A(5)                                             THE08660
      CHARACTER*46 FMT3B                                                THE08670
      CHARACTER*21 FMT4A(5)                                             THE08680
      CHARACTER*13 FMT4B                                                THE08690
      CHARACTER*1 PLTCHR,DIGIT(10)                                      THE08700
      CHARACTER*72 TITLE(20),TEXT                                       THE08710
      CHARACTER*4 BLANK,TEST                                            THE08720
      CHARACTER*72 FIN, FOUT                                            THE08730
      LOGICAL DONE, NEG ,bigerr                                         THE08740
      REAL*8 NORML, MEAN                                                THE08750
      DIMENSION B(6),BMIN(6),BMAX(6),Y(30),Z(30),PDERIV(6,30),          THE08760
     @          BI(6),ZI(30),CONV(5),X(1,30),NORMAL(7,7),WORK(7),       THE08770
     @          SCALE(6),GRAD(6),DIAG(7),CONST(3),VAR(6),STD(6),        THE08780
     @          IPIVOT(6), COVAR(6,6), NORML(7,7), DIAG2(6),            THE08790
     @          IX(180), XTEMP(180),CP(30),AI(6),DELPCT(30),BS(5),      THE08800
     @          BSMIN(5),BSMAX(5),TEMPIN(*)                             THE08810
      REAL*8 SUMSQ,B,BMIN,BMAX,BI,NORMAL,GRAD,WORK,DIAG,CONST           THE08820
     @,SCALE,BTEMP,MLARGE,PLARGE,AVG,COEF2,COEF3,AI,CP,BS,BSMIN,BSMAX,  THE08830
     @BTMP,TEMPIN                                                       THE08840
      INTEGER PARM,DATA,ERROR,ERROR2,INPUT,OUTPUT,                      THE08850
     %        OPTION,FREQ,MAXPAR,MAXIND,                                THE08860
     %        MAXDAT,NTITLE,MAXIT,NIND,OPT2,INTMSK,                     THE08870
     %        PRTOPT,PLTOPT,MAXCON,NT,BFLAG                             THE08880
      EQUIVALENCE (IPIVOT(1), SCALE(1)), (COVAR(1,1), PDERIV(1,1)),     THE08890
     @            (NORML(1,1), NORMAL(1,1)), (DIAG2(1), DIAG(1)),       THE08900
     @            (IOVER1, OVERF1), (IOVER2, OVERF2),                   THE08910
     @            (IX(1), PDERIV(1,1)), (XTEMP(1), PDERIV(1,1))         THE08920
                                                                        THE08930
                                                                        THE08940
                                                                        THE08950
      COMMON /DINFO/ IPARM1, IDATA, IPARM, INDPT                        THE08960
      COMMON /CONST/ CONST,NTERMS                                       THE08970
      COMMON /PDV/ PDERIV                                               THE08980
      COMMON/FLAG/IREPEA,IFLAG,DELPCT                                   THE08990
        SAVE BSMAX,BSMIN,BS                                             THE09000
      DATA INPUT/5/,OUTPUT/6/,BLANK/'    '/,MAXPAR/30/,MAXIND/10/       THE09010
      DATA MAXDAT/400/,MAXCON/40/                                       THE09020
      DATA FMT1/'(8F10.5)'/,DONE/.TRUE./                                THE09030
      DATA FMT2/'(3X  /(9H       B(,I2,4H)  =,1PE16.7))'/               THE09040
      DATA FMT3A/'(3X  /8H        ,1(1HX,I1,12X)\)',                     THE09050
     &           '(3X  /8H        ,2(1HX,I1,12X)\)',                    THE09060
     &           '(3X  /8H        ,3(1HX,I1,12X)\)',                    THE09070
     &           '(3X  /8H        ,4(1HX,I1,12X)\)',                    THE09080
     &           '(3X  /8H        ,5(1HX,I1,12X)\)'/                    THE09090
      DATA FMT3B/'(2HX,13X,2HY,13X,2HYC,11X,4HYC-Y,10X,3HPCT//) '/      THE09100
      DATA FMT4A/'(5H     ,4(1PE14.5)\)','(5H     ,5(1PE14.5)\)',       THE09110
     &           '(5H     ,6(1PE14.5)\)','(5H     ,7(1PE14.5)\)',       THE09120
     &           '(5H     ,8(1PE14.5)\)'/, FMT4B/'(1H ,0PF10.2)'/       THE09130
      DATA DIGIT/'0','1','2','3','4','5','6','7','8','9'/,PLTCHR/'O'/   THE09140
       MLARGE = -1.6D+75                                                THE09150
       PLARGE =  1.6D+75                                                THE09160
       bigerr = .false.
                                                                        THE09170
       CALL RESET                                                       THE09180
                                                                        THE09190
      IPARM = MAXPAR                                                    THE09200
      IDATA = MAXDAT                                                    THE09210
      IPARM1 = MAXPAR + 1                                               THE09220
      INDPT = MAXIND                                                    THE09230
100   NTITLE = 1                                                        THE09240
110   CONTINUE                                                          THE09250
1     FORMAT(A72)                                                       THE09260
                                                                        THE09270
120     CONTINUE                                                        THE09280
                                                                        THE09290
        FREQ=1000                                                       THE09300
        MAXIT=65                                                        THE09310
        maxrep=20
        CONV(1)=1.0E-2                                                  THE09320
        CONV(2)=1.0E-2                                                  THE09330
        CONV(3)=0.0                                                     THE09340
        CONV(3)=0.0                                                     THE09350
        CONV(4)=0.0                                                     THE09360
        CONV(5)=0.0                                                     THE09370
C       DEFINE  DEFAULTS                                                THE09380
        IREPEA=0                                                        THE09390
        IREPNU=0                                                        THE09400
        OPTION=1                                                        THE09410
        OPT2=1                                                          THE09420
        NPRT=1                                                          THE09430
        PRTOPT=1                                                        THE09440
        PLTOPT=1                                                        THE09450
        ICALL=0                                                         THE09460
        PARM=5                                                          THE09470
        NTERMS=3                                                        THE09480
        NIND=1                                                          THE09490
        ISTART=1                                                        THE09500
                                                                        THE09510
8     FORMAT(3D10.8)                                                    THE09520
        ICONST=3                                                        THE09530
        FVALUE=0.                                                       THE09540
        TVALUE=0.                                                       THE09550
32    FORMAT(7E10.0,8X,I2)                                              THE09560
11    FORMAT(8D10.0)                                                    THE09570
                                                                        THE09580
C       CONST(1): GAS CONSTANT CONVERSION FACTOR;   CONST(1)=1.0 (CAL)  THE09590
C                                                   CONST(1)=4.184 (J)  THE09600
C        CONST(2) = TOTAL NUMBER OF VIBRATIONAL DEGREES OF FREEDOM      THE09610
C                                                                       THE09620
C       CONST(3) = 0., IF NON-LINEAR MOLECULE; = 1., IF LINEAR          THE09630
170   NIND1 = NIND + 1                                                  THE09640
                                                                        THE09650
C                                                                       THE09660
C             FILL DATA POINT ARRAY WITH CP VALUES.                     THE09670
C                                                                       THE09680
180   CONTINUE                                                          THE09690
        DO 190 IK1=1,NT                                                 THE09700
        Y(IK1)=CP(IK1)                                                  THE09710
        X(1,IK1)=TEMPIN(IK1)                                            THE09720
190   CONTINUE                                                          THE09730
        DATA=NT+1                                                       THE09740
195     CONTINUE                                                        THE09750
                                                                        THE09760
       if(icp.eq.1) then
        x(1,data)=thigh
    else     
                       X(1,DATA)=99999.
        endif                                
                        Y(DATA)=CPINF                                   THE09780
        CPRATI=0.5                                                      THE09790
        IREPEA=0                                                        THE09800
        IREPNU=0                                                        THE09810
        IF(IFLAG.EQ.1)THEN                                              THE09820
                DO 7634 IK1=1,5                                         THE09830
                B(IK1)=BS(IK1)                                          THE09840
                BMAX(IK1)=BSMAX(IK1)                                    THE09850
                BMIN(IK1)=BSMIN(IK1)                                    THE09860
7634            CONTINUE                                                THE09870
        GOTO 9865                                                       THE09880
        ENDIF                                                           THE09890
9268    CONTINUE                                                        THE09900
        IF(IREPEA.NE.0)THEN                                             THE09910
                IREPNU=IREPNU+1                                         THE09920
                IREPEA=0                                                THE09930
                IF(IREPNU.GE.maxrep)THEN
                        WRITE(25,*) ' HOE FIT FAILURE '                 THE09950
                        WRITE(25,*)' defaulting to EXP fit'             THE09960
                        IFLAG=1                                         THE09970
                        RETURN                                          THE09980
C                       IF(IFLAG.EQ.0)THEN                              THE09990
C               WRITE(20,*) ' POLYNOMIAL SMOOTHING OF DATA IN PROGRESS' THE10000
C                                                                       THE10010
C                       ELSE                                            THE10020
C       WRITE(20,*)' TEMPERATURE RANGE DEFAULTING TO 300 - 2000 K'      THE10030
                                                                        THE10040
                ENDIF                                                   THE10050
                IF(IREPNU.GT.1)THEN                                     THE10060
                        IBS=0                                           THE10070
                        DO 9655 I=1,5                                   THE10080
                        DELTAB=(ABS(BS(I)-B(I))/BS(I))                  THE10090
                        IF(DELTAB.GT.0.01)IBS=1                         THE10100
9655                    CONTINUE                                        THE10110
                        IF(IBS.EQ.0)THEN                                THE10120
c                               IREPNU=IREPEA                           THE10130
                                IREPEA=2                                THE10140
                                GOTO 600                                THE10150
                        ENDIF                                           THE10160
                ENDIF                                                   THE10170
                IF(DELPCT(1).LE.0.0.AND.IREPNU.EQ.1)THEN                THE10180
                        CPRATI=0.75                                     THE10190
                        DO 9651 I=1,5                                   THE10200
                        BS(I)=B(I)                                      THE10210
                        BSMIN(I)=BMIN(I)                                THE10220
                        BSMAX(I)=BMAX(I)                                THE10230
9651                    CONTINUE                                        THE10240
                        BFLAG=0                                         THE10250
                ELSE IF(DELPCT(1).GT.0.0.AND.IREPNU.EQ.1)THEN           THE10260
                        CPRATI=0.25                                     THE10270
                        DO 9652 I=1,5                                   THE10280
                        BS(I)=B(I)                                      THE10290
                        BSMIN(I)=BMIN(I)                                THE10300
                        BSMAX(I)=BMAX(I)                                THE10310
9652                    CONTINUE                                        THE10320
                        BFLAG=1                                         THE10330
                ELSEIF((DELPCT(1).LT.0.0.AND.DELPCT(2).LT.0.0)          THE10340
     +                  .AND.IREPNU.GT.1)THEN                           THE10350
                        IF(ABS(B(3)-BMIN(3)).LE.0.1)THEN                THE10360
                                BMIN(3)=BMIN(3)-50.                     THE10370
                        ELSE                                            THE10380
                                IF(BFLAG.EQ.1)THEN                      THE10390
                                        DO 9653 I=1,5                   THE10400
                                        BTMP=B(I)                       THE10410
                                        B(I)=(B(I)+BS(I))/2.0           THE10420
                                        BS(I)=BTMP                      THE10430
                                IF(BMAX(I).LT.BSMAX(I))BMAX(I)=BSMAX(I) THE10440
                                IF(BMIN(I).GT.BSMIN(I))BMIN(I)=BSMIN(I) THE10450
9653                            CONTINUE                                THE10460
                                ENDIF                                   THE10470
                        ENDIF                                           THE10480
                        BFLAG=0                                         THE10490
                        GOTO 9865                                       THE10500
                ELSEIF((DELPCT(1).GT.0.0.AND.DELPCT(2).GT.0.0)          THE10510
     +                  .AND.IREPNU.GT.1)THEN                           THE10520
                        IF(ABS(B(3)-BMIN(3)).LE.0.1)THEN                THE10530
                                BMIN(3)=BMIN(3)+50.                     THE10540
                                B(3)=BMIN(3)+50.                        THE10550
                        ELSE                                            THE10560
                                                                        THE10570
                                IF(BFLAG.EQ.0)THEN                      THE10580
                                        DO 9654 I=1,5                   THE10590
                                        BTMP=B(I)                       THE10600
                                        B(I)=(B(I)+BS(I))/2.0           THE10610
                                        BS(I)=BTMP                      THE10620
                                IF(BMAX(I).LT.BSMAX(I))BMAX(I)=BSMAX(I) THE10630
                                IF(BMIN(I).GT.BSMIN(I))BMIN(I)=BSMIN(I) THE10640
9654                            CONTINUE                                THE10650
                                BLFAG=1                                 THE10660
                                ELSE                                    THE10670
                                        BMIN(1)=0.                      THE10680
                                        BMIN(2)=0.                      THE10690
                                        BMIN(3)=100.                    THE10700
                                        BMIN(4)=300.                    THE10710
                                        BMIN(5)=700.                    THE10720
                                        BMAX(1)=CONST(2)                THE10730
                                        BMAX(2)=CONST(2)                THE10740
                                        BMAX(3)=999.                    THE10750
                                        BMAX(4)=1999.                   THE10760
                                        BMAX(5)=3600.                   THE10770
                                ENDIF                                   THE10780
                                                                        THE10790
                        ENDIF                                           THE10800
                        GOTO 9865                                       THE10810
                ENDIF                                                   THE10820
                                                                        THE10830
        ENDIF                                                           THE10840
                                                                        THE10850
                                                                        THE10860
C       set initial parameter values, max & min                         THE10870
C       RANGE DETERMINED BY CPRATI                                      THE10880
C               CPRATI < 0.4  :  VERY HIGH FREQUENCY WEIGHTED           THE10890
C         0.4 < CPRATI < 0.65 :  MID FREQUENCY WEIGHTED                 THE10900
C               CPRATI > 0.65 :  LOW FREQUENCY WEIGHTED                 THE10910
        IF(CPRATI.LT.0.4)THEN                                           THE10920
C               starting point : high frequency weighted                THE10930
                B(1)=0.2*CONST(2)                                       THE10940
                B(2)=0.5*CONST(2)                                       THE10950
                B(3)=600.                                               THE10960
                B(4)=1500.                                              THE10970
                B(5)=3000.                                              THE10980
                BMIN(1)=0.                                              THE10990
                BMIN(2)=0.                                              THE11000
                BMIN(3)=350.                                            THE11010
                BMIN(4)=900.                                            THE11020
                BMIN(5)=1500.                                           THE11030
                BMAX(1)=0.6*CONST(2)                                    THE11040
                BMAX(2)=CONST(2)                                        THE11050
                BMAX(3)=1500.                                           THE11060
                BMAX(4)=2500.                                           THE11070
                BMAX(5)=3600.                                           THE11080
                IF(IREPNU.NE.1)THEN                                     THE11090
                        DO 911 IKJ=1,5                                  THE11100
911                     BS(IKJ)=B(IKJ)                                  THE11110
                ENDIF                                                   THE11120
                                                                        THE11130
        ELSE IF(CPRATI.GE.0.4.AND.CPRATI.LE.0.63)THEN                   THE11140
C               starting point : spread over range                      THE11150
                B(1)=0.3*CONST(2)                                       THE11160
                B(2)=0.3*CONST(2)                                       THE11170
                B(3)=600.                                               THE11180
                B(4)=1200.                                              THE11190
                B(5)=2400.                                              THE11200
                BMIN(1)=0.                                              THE11210
                BMIN(2)=0.                                              THE11220
                BMIN(3)=250.                                            THE11230
                BMIN(4)=800.                                            THE11240
                BMIN(5)=900.                                            THE11250
                BMAX(1)=CONST(2)                                        THE11260
                BMAX(2)=CONST(2)                                        THE11270
                BMAX(3)=1500.                                           THE11280
                BMAX(4)=2500.                                           THE11290
                BMAX(5)=4000.                                           THE11300
                IF(IREPNU.NE.1)THEN                                     THE11310
                        DO 9121 IKJ=1,5                                 THE11320
9121                    BS(IKJ)=B(IKJ)                                  THE11330
                ENDIF                                                   THE11340
        ELSE IF(CPRATI.GT.0.63)THEN                                     THE11350
C               starting point: low frequency weighted                  THE11360
C               parameters are constrained to low frequencies;          THE11370
C               with lower minimum frequencies acceptable.              THE11380
                B(1)=0.5*CONST(2)                                       THE11390
                B(2)=0.4*CONST(2)                                       THE11400
                B(3)=300.                                               THE11410
                B(4)=500.                                               THE11420
                B(5)=700.                                               THE11430
                BMIN(1)=CONST(2)/10.                                    THE11440
                BMIN(2)=CONST(2)/10.                                    THE11450
                BMIN(3)=100.                                            THE11460
                BMIN(4)=401.                                            THE11470
                BMIN(5)=601.                                            THE11480
                BMAX(1)=CONST(2)                                        THE11490
                BMAX(2)=CONST(2)                                        THE11500
                BMAX(3)=450.                                            THE11510
                BMAX(4)=1000.                                           THE11520
                BMAX(5)=2000.                                           THE11530
                IF(IREPNU.NE.1)THEN                                     THE11540
                        DO 913 IKJ=1,5                                  THE11550
913                     BS(IKJ)=B(IKJ)                                  THE11560
                ENDIF                                                   THE11570
        ENDIF                                                           THE11580
                                                                        THE11590
9865    IOUTLI=0                                                        THE11600
9866    CONTINUE                                                        THE11610
        IREDUC=0                                                        THE11620
      DO 200 ITER = ISTART,MAXIT                                        THE11630
         CALL SOLVE (PARM,DATA,ERROR,ERROR2,B,BMIN,BMAX,Y,Z,            THE11640
     &               SUMSQ,CONV,NORMAL,SCALE,                           THE11650
     &               WORK,GRAD,BI,ZI,DIAG,X )                           THE11660
        IF(SUMSQ.LT.CONV(1)) GOTO 600                                   THE11670
C                                                                       THE11680
                                                                        THE11690
C             CHECK FOR TERMINAL ERROR CODES FROM SUBROUTINE            THE11700
C             SOLVE.                                                    THE11710
C                                                                       THE11720
         IF (ERROR .EQ. -1)   GO TO 900                                 THE11730
         IF (ERROR .EQ. -2)   GO TO 900                                 THE11740
         IF (ERROR .EQ. -3)  THEN                                       THE11750
                IF(IOUTLI.EQ.1)GOTO 900                                 THE11760
                IOUTLI=1                                                THE11770
                BMIN(1)=0.                                              THE11780
                BMIN(2)=0.                                              THE11790
                BMIN(3)=75.                                             THE11800
                BMIN(4)=150.                                            THE11810
                BMIN(5)=300.                                            THE11820
                BMAX(1)=CONST(2)                                        THE11830
                BMAX(2)=CONST(2)                                        THE11840
                BMAX(3)=1500.                                           THE11850
                BMAX(4)=3500.                                           THE11860
                BMAX(5)=4500.                                           THE11870
                GOTO 9866                                               THE11880
        ENDIF                                                           THE11890
         IF (ERROR .EQ. 0)    GO TO 600                                 THE11900
                                                                        THE11910
20001   CONTINUE                                                        THE11920
                                                                        THE11930
200      CONTINUE                                                       THE11940
      ITER=MAXIT                                                        THE11970
      GO TO 600                                                         THE11980
C                                                                       THE11990
C             SOLUTION FOUND                                            THE12000
C                                                                       THE12010
600     CONTINUE                                                        THE12020
        AI(1)=B(1)                                                      THE12030
        AI(2)=B(2)                                                      THE12040
        AI(3)=CONST(2)-(B(1)+B(2))                                      THE12050
        AI(4)=B(3)                                                      THE12060
        AI(5)=B(4)                                                      THE12070
        AI(6)=B(5)                                                      THE12080
                                                                        THE12090
      IF (NPRT .LT. NIND)  NPRT = NIND                                  THE12100
      IF (NPRT .GT. 5)       NPRT = 5                                   THE12110
      SUMP = 0.0
      pctmax = 0.
      DO 630 I = 1,DATA                                                 THE12130
         ZI(I) = Z(I) - Y(I)                                            THE12140
         IF (Y(I) .EQ. 0.0)  GO TO 620                                  THE12150
         PCT = 100.0 * ZI(I) / Y(I)                                     THE12160
         GO TO 625                                                      THE12170
620      PCT = 0.0                                                      THE12180
625     CONTINUE
        if (pct.gt.pctmax) pctmax = pct

c      i think irepea = 2 means fit values aren't changing
c      thus if error still too big we quit with warning <ayc 2/93>

        IF(ABS(PCT).GT.4.5.AND.IREPEA.NE.2)THEN                         THE12210
                IREPEA=1                                                THE12220
                DELPCT(I)=PCT                                           THE12230
        ELSEIF(ABS(PCT).GT.4.5.AND.IREPEA.EQ.2)THEN
                bigerr = .true.                                         THE12270
        ENDIF                                                           THE12290
         SUMP = SUMP + ABS(PCT)                                         THE12300
630     CONTINUE                                                        THE12310
        if (bigerr) then
           irepea = 0
           write(25,666) pctmax,irepnu
666        format('  WARNING: HOE FIT WITH MAX ERROR OF ',f4.1,
     2        '% ON ITERATION ',i2)
           write(25,*) ' defaulting to EXP fit'
           iflag = 1
        endif
        IF(IREPEA.NE.0)THEN                                             THE12320
                CALL RESET                                              THE12330
                GOTO 9268                                               THE12340
        ENDIF                                                           THE12350
C if irepeat = 0 then adequate solution has been found.                 THE12360
C     this is a normal return                                           THE12370
        RETURN                                                          THE12380
900     CONTINUE                                                        THE12390
        WRITE(25,*)' HOE FIT FAILED'                                    THE12400
        WRITE(25,*)' DEFAULTING TO EXP FIT'                             THE12410
        IFLAG=1                                                         THE12420
        RETURN                                                          THE12430
      END                                                               THE12500
                                                                        THE12510
      SUBROUTINE SOLVE (NPARM,NDATA,ERROR,ERROR2,B,BMIN,BMAX,           THE12520
     #                  Y,Z,SUMSQ,CONV,NORMAL,SCALE,                    THE12530
     %                  WORK,GRAD,BI,ZI,DIAG,X )                        THE12540
                                                                        THE12550
        IMPLICIT REAL*8(A-H,O-Z)                                        THE12560
      COMMON /DINFO/ IPARM1, IDATA, IPARM, INDPT                        THE12570
      COMMON /LINITL/ INITL                                             THE12580
      COMMON /PDV/ PDERIV                                               THE12590
      INTEGER ERROR,ERROR2                                              THE12600
      REAL*8 NUFAC,LAMBDA                                               THE12610
      DOUBLE PRECISION SUMSQ,B,BMIN,BMAX,BI,NORMAL,SCALEF,GRAD,WORK,    THE12620
     `                 DIFF,SCALE,SUMSQR,DABS,DSQRT,DEL,DIAG, RATIO,    THE12630
     &                 POIN01, POI001, CONV1, CONV2, CONV3              THE12640
      DIMENSION B(6),BMIN(6),BMAX(6),Y(*),Z(*),                         THE12650
     %          PDERIV(6,30),BI(6),ZI(*),CONV(*),                       THE12660
     $          X(1,30),NORMAL(7,7),WORK(*),                            THE12670
     *          SCALE(*),GRAD(*),DIAG(*)                                THE12680
      LOGICAL INITL                                                     THE12690
      IF (INITL)   GO TO 210                                            THE12700
C                                                                       THE12710
C             INITIALIZE AND CHECK INPUT PARAMETERS                     THE12720
C                                                                       THE12730
C             THE NUMBER OF PARAMETERS NPARM MUST BE LESS               THE12740
C             THAN OR EQUAL TO THE NUMBER OF DATA POINTS, NDATA         THE12750
C             IF NOT, SET ERROR TO -1 AND RETURN                        THE12760
C                                                                       THE12770
      IF (NPARM .LE. NDATA)  GO TO 10                                   THE12780
      ERROR = -1                                                        THE12790
      RETURN                                                            THE12800
10    DO 20 J = 1,NPARM                                                 THE12810
         IF (B(J) .GE. BMIN(J) .AND. B(J) .LE. BMAX(J))    GO TO 20     THE12820
         ERROR = -3                                                     THE12830
         RETURN                                                         THE12840
20       CONTINUE                                                       THE12850
      NPARM1 = NPARM + 1                                                THE12860
      INITL = .TRUE.                                                    THE12870
C                                                                       THE12880
C             NUFAC  -  THE NU FACTOR IN THE SOURCE PAPER               THE12890
C             LAMBDA  -  THE LAMBDA FACTOR IN THE SOURCE PAPER          THE12900
C                                                                       THE12910
      IF (CONV(4) .LE. 0.0)  CONV(4) = 10.0                             THE12920
      NUFAC = CONV(4)                                                   THE12930
      IF (CONV(5) .LE. 0.0)  CONV(5) = 0.01                             THE12940
      LAMBDA = CONV(5)                                                  THE12950
      IF (CONV(3) .LE. 0.0)  CONV(3) = 0.01                             THE12960
      POIN01 = CONV(3)                                                  THE12970
      POI001 = POIN01 / 10.0 D 00                                       THE12980
C                                                                       THE12990
C             IF THE CONVERGENCE CRITERIA IS LESS THAN OR               THE13000
C             EQUAL TO ZERO, SET IT TO .00001.                          THE13010
C                                                                       THE13020
      IF (CONV(1) .LE. 0.0)  CONV(1) = .00001                           THE13030
      IF (CONV(2) .LE. 0.0)  CONV(2) = .00001                           THE13040
      CONV1 = CONV(1)                                                   THE13050
      CONV2 = CONV(2)                                                   THE13060
      CONV3 = CONV1 + 1.0 D 00                                          THE13070
C                                                                       THE13080
C             COMPUTE THE SUM OF SQUARES ((Y - Z) ** 2) FOR THE         THE13090
C             FIRST SET OF B VALUES.                                    THE13100
C                                                                       THE13110
      SUMSQ = 0.0 D 00                                                  THE13120



      DO 200 J = 1,NDATA                                                THE13130
         CALL FOFX (J,X(1,J),B,Y(J),Z(J))                               THE13140
 
200      SUMSQ = SUMSQ + (Y(J) - Z(J)) ** 2                             THE13150
      IF (SUMSQ .EQ. 0.0D00)      GO TO 530                             THE13160
C                                                                       THE13170
                                                                        THE13180
                                                                        THE13190
C*  CALL PARTIAL DERIVATIVE SUBROUTINE (E. RITTER  12/30/88)            THE13200
C        DETERMINE DERIV OF FOFX W/R/T PARAMETERS                       THE13210
C       COMPUTE NORMAL EQUATIONS                                        THE13220
                                                                        THE13230
                                                                        THE13240
210     CALL PD(NDATA,X,NPARM,B,Y,ZZ,Z,PDERIV,NORMAL,WORK,POIN01)       THE13250
                                                                        THE13260
                                                                        THE13270
C                                                                       THE13280
C             START THE PROCEDURE TO ESTIMATE THE PARAMETERS.           THE13290
C                                                                       THE13300
      LAMBDA = LAMBDA / NUFAC                                           THE13310
      DO 312 J1 = 1, NPARM                                              THE13320
         J2 = J1 + 1                                                    THE13330
         DO 310 J3 = J2,NPARM1                                          THE13340
                                                                        THE13350
310         NORMAL(J1,J3) = NORMAL(J3,J1)                               THE13360
312      DIAG(J1) = NORMAL(J1,J1)                                       THE13370
      DIAG(NPARM1) = NORMAL(NPARM1,NPARM1)                              THE13380
      GO TO 315                                                         THE13390
300   DO 313 J1 = 1, NPARM                                              THE13400
         J2 = J1 + 1                                                    THE13410
         DO 314 J3 = J2,NPARM1                                          THE13420
314         NORMAL (J3,J1) = NORMAL (J1,J3)                             THE13430
313      NORMAL (J1,J1) = DIAG (J1)                                     THE13440
      NORMAL (NPARM1,NPARM1) = DIAG (NPARM1)                            THE13450
C                                                                       THE13460
C             CALCULATE AND STORE THE MATRIX SCALING FACTORS            THE13470
C             IN THE SCALE VECTOR                                       THE13480
C                                                                       THE13490
315   DO 320 J = 1,NPARM                                                THE13500
320      SCALE(J) = DSQRT(NORMAL(J,J))                                  THE13510
C                                                                       THE13520
C             SCALE THE MATRIX AND SAVE THE GRADIENT VECTOR             THE13530
C             FOR FUTURE IN COMPUTING THE ANGLE BETWEEN THE GRADIENT    THE13540
C             AND THE GAUSS - NEWTON - MARQUARDT  VECTOR.               THE13550
C                                                                       THE13560
      DO 340 J1 = 1,NPARM                                               THE13570
         DO 330 J2 = J1,NPARM                                           THE13580
330         NORMAL(J2,J1) = NORMAL(J2,J1) / (SCALE(J1) * SCALE(J2))     THE13590
         NORMAL(NPARM1,J1) = NORMAL(NPARM1,J1) / SCALE(J1)              THE13600
340      GRAD(J1) = NORMAL(NPARM1,J1)                                   THE13610
C                                                                       THE13620
C             NORMAL EqUATIONS NOW SCALED.  THE GRAD VECTOR             THE13630
C             CONTAINS THE GRADIENT VECTOR.                             THE13640
C                                                                       THE13650
      DO 350 J1 = 1,NPARM                                               THE13660
350      NORMAL(J1,J1) = NORMAL(J1,J1) + LAMBDA                         THE13670
      ERROR2 = 0                                                        THE13680
      DO 356 I = 1,NPARM                                                THE13690
         J1 = I + 1                                                     THE13700
         DO 356 J2 = 1,I                                                THE13710
            IF (NORMAL(J2,J2) .GT. 0.0 D 00)     GO TO 351              THE13720
C                                                                       THE13730
C             A SINGULAR MATRIX WAS ENCOUNTERED                         THE13740
C             GENERATE ERROR CODE -2 AND TERMINATE                      THE13750
C                                                                       THE13760
            ERROR2 = I                                                  THE13770
            ERROR = -2                                                  THE13780
            RETURN                                                      THE13790
351         RATIO = NORMAL (J1,J2) / NORMAL (J2,J2)                     THE13800
            DO 354 J = J1,NPARM1                                        THE13810
354            NORMAL (J,J1) = NORMAL (J,J1) - RATIO * NORMAL (J,J2)    THE13820
356         NORMAL (J1,J2) = RATIO                                      THE13830
      WORK (NPARM) = NORMAL (NPARM1,NPARM)                              THE13840
      IF (NPARM .LE. 1)      GO TO 360                                  THE13850
      J2 = NPARM - 1                                                    THE13860
      I = J2                                                            THE13870
      J = NPARM1                                                        THE13880
      DO 359 J3 = 2,NPARM                                               THE13890
         WORK (J2) = NORMAL (J,I)                                       THE13900
                                                                        THE13910
         L = NPARM                                                      THE13920
         DO 358 J4 = 2,J3                                               THE13930
            J = J - 1                                                   THE13940
            WORK (J2) = WORK (J2) - WORK (L) * NORMAL (J,I)             THE13950
358         L = L - 1                                                   THE13960
         J2 = J2 - 1                                                    THE13970
         J = NPARM1                                                     THE13980
359      I = I - 1                                                      THE13990
C                                                                       THE14000
C             NORMAL EQUATIONS SOLVED.  CORRECTION VECTOR               THE14010
C             IS STORED IN WORK.                                        THE14020
C                                                                       THE14030
C             COMPUTE THE COSINE OF THE ANGLE BETWEEN THE               THE14040
C             GRADIENT VECTOR GRAD AND THE CORRECTION VECTOR ,WORK      THE14050
C                                                                       THE14060
360   SA = 0.0                                                          THE14070
      SB = 0.0                                                          THE14080
      SC = 0.0                                                          THE14090
      DO 370 J = 1,NPARM                                                THE14100
         SA = SA + WORK(J) * GRAD(J)                                    THE14110
         SB = SB + WORK(J) ** 2                                         THE14120
370      SC = SC + GRAD(J) ** 2                                         THE14130
      CANGLE = ABS(SA / (SQRT(SB * SC)))                                THE14140
C                                                                       THE14150
C             INITIALIZE A SCALE FACTOR SCALEF, TO BE USED IN           THE14160
C             THE EVENT THE ANGLE BETWEEN THE VECTORS IS                THE14170
C             GREATER THAN 30 DEGREES.  ( COSINE OF THE ANGLE           THE14180
C             LESS THAN .866 )                                          THE14190
C                                                                       THE14200
      SCALEF = 1.0 D 00                                                 THE14210
380   DO 400 J = 1,NPARM                                                THE14220
         BI(J) = B(J) + SCALEF * WORK(J) / SCALE(J)                     THE14230
         IF (BI(J) .LT. BMIN(J))  GO TO 430                             THE14240
         IF (BI(J) .GT. BMAX(J))  GO TO 430                             THE14250
400      CONTINUE                                                       THE14260
      SUMSQR = 0.0 D 00                                                 THE14270
      DO 405 J = 1,NDATA                                                THE14280
         CALL FOFX (J,X(1,J),BI,Y(J),ZI(J))                             THE14290
405      SUMSQR = SUMSQR + (Y(J) - ZI(J)) ** 2                          THE14300
C                                                                       THE14310
C             NOW SEE IF THE NEW ESTIMATES HAVE REDUCED THE             THE14320
C             RESIDUAL SUM OF SQUARES.  IF NOT, AND THE ANGLE           THE14330
C             IS LESS THAN 30 DEGREES, ADD ONLY 1/2 OF THE              THE14340
C             CORRECTION PREVIOUSLY USED AND RECOMPUTE THE              THE14350
C             RESIDUAL SUM OF SQUARES.  IF THE ANGLE EXCEEDS            THE14360
C             30 DEGREES, MULTIPLY LAMBDA BY NUFAC AND SOLVE            THE14370
C             NORMAL EQUATIONS AGAIN, AND THEN RECOMPUTE THE            THE14380
C             RESIDUAL SUM OF SQUARES.                                  THE14390
C                                                                       THE14400
      IF (SUMSQR / SUMSQ .LE. CONV3)     GO TO 450                      THE14410
      IF (CANGLE .LE.0.866)  GO TO 440                                  THE14420
C                                                                       THE14430
C             ANGLE LESS THAN 30 DEGREES.  GO BACK AND ADD              THE14440
C             HALF THE CORRECTION.                                      THE14450
C                                                                       THE14460
430   SCALEF = SCALEF / 2.0 D 00                                        THE14470
      GO TO 380                                                         THE14480
C                                                                       THE14490
C             ANGLE GREATER THAN 30 DEGREES INCREASE LAMBDA             THE14500
C             AND RESOLVE                                               THE14510
C                                                                       THE14520
440   LAMBDA = LAMBDA * NUFAC                                           THE14530
      GO TO 300                                                         THE14540
C                                                                       THE14550
C             NEW ESTIMATES HAVE BEEN FOUND WHICH REDUCE THE            THE14560
C             RESIDUAL SUM OF SQUARES, AND THESE ESTIMATES              THE14570
C             ARE IN THE BI ARRAY.  THE CORRESPONDING RESIDUALS         THE14580
C             AND THE RESIDUAL SUM OF SQUARES ARE IN THE ZI             THE14590
C             ARRAY AND SUMSQR, RESPECTIVELY.                           THE14600
C                                                                       THE14610
C             PERFORM THE CONVERGENCE TESTS.                            THE14620
C                                                                       THE14630
450   DO 460 J = 1,NDATA                                                THE14640
460      Z(J) = ZI(J)                                                   THE14650
C                                                                       THE14660
C             ERROR2 WILL CONTAIN THE NUMBER OF PARAMETERS              THE14670
C             THAT FAIL THE CONVERGENCE TEST.  ERROR WILL BE            THE14680
C             SET TO  -5, IF ERROR2 IS GREATER THAN 0.                  THE14690
C                                                                       THE14700
      DO 520 J = 1,NPARM                                                THE14710
         DIFF = DABS(BI(J) - B(J))                                      THE14720
         B(J) = BI(J)                                                   THE14730
         IF (DIFF / (1.0D-20 + DABS(B(J))) .GE. CONV2)                  THE14740
     :        ERROR2 = ERROR2 + 1                                       THE14750
520      CONTINUE                                                       THE14760
      DIFF = DABS (SUMSQR - SUMSQ)                                      THE14770
      SUMSQ = SUMSQR                                                    THE14780
      IF (SUMSQ .EQ. 0.0 D 00)    GO TO 530                             THE14790
      ERROR = 0                                                         THE14800
      IF (DIFF / SUMSQ .GE. CONV1)   ERROR = -4                         THE14810
      IF (ERROR2 .GT. 0)     ERROR = -5                                 THE14820
      RETURN                                                            THE14830
530   ERROR = 0                                                         THE14840
      RETURN                                                            THE14850
      END                                                               THE14860
C                                                                       THE14870
C             CALL RESET FOR MULTIPLE CASES                             THE14880
C                                                                       THE14890
      SUBROUTINE RESET                                                  THE14900
                                                                        THE14910
      LOGICAL INITL                                                     THE14920
      COMMON /LINITL/ INITL                                             THE14930
      INITL = .FALSE.                                                   THE14940
      RETURN                                                            THE14950
      END                                                               THE14960
                                                                        THE14970
      SUBROUTINE FACTOR (A, W, IPIVOT, D, N, IDIM, IFLAG)               THE14980
                                                                        THE14990
        IMPLICIT REAL*8(A-H,O-Z)                                        THE15000
        DOUBLE PRECISION D                                              THE15010
      DIMENSION A(7,7), W(7,7), IPIVOT(*), D(*)                         THE15020
      IFLAG = 1                                                         THE15030
C     INITIALIZE W, IPIVOT, D                                           THE15040
      DO 10 I = 1, N                                                    THE15050
         IPIVOT(I) = I                                                  THE15060
         ROWMAX = 0.                                                    THE15070
         DO 9 J = 1, N                                                  THE15080
            W(I,J) = A(I,J)                                             THE15090
            DWIJ = ABS(W(I,J))                                          THE15100
9           ROWMAX = DMAX1 (ROWMAX, DWIJ)                               THE15110
         IF (ROWMAX .EQ. 0.0)  GO TO 999                                THE15120
10       D(I) = ROWMAX                                                  THE15130
C     GAUSS ELIMINATION WITH SCALED PARTIAL PIVOTING.                   THE15140
      NM1 = N - 1                                                       THE15150
      IF (NM1 .EQ. 0)    RETURN                                         THE15160
      DO 20 K = 1, NM1                                                  THE15170
         J = K                                                          THE15180
         KP1 = K + 1                                                    THE15190
         IP = IPIVOT(K)                                                 THE15200
         COLMAX = ABS(W(IP,K)) / D(IP)                                  THE15210
         DO 11 I = KP1, N                                               THE15220
            IP = IPIVOT(I)                                              THE15230
            AWIKOV = ABS (W(IP,K)) / D(IP)                              THE15240
            IF (AWIKOV .LE. COLMAX)     GO TO 11                        THE15250
            COLMAX = AWIKOV                                             THE15260
            J = I                                                       THE15270
11          CONTINUE                                                    THE15280
         IF (COLMAX .EQ. 0.0)    GO TO 999                              THE15290
C                                                                       THE15300
         IPK = IPIVOT(J)                                                THE15310
         IPIVOT(J) = IPIVOT(K)                                          THE15320
         IPIVOT(K) = IPK                                                THE15330
         DO 20 I = KP1, N                                               THE15340
            IP = IPIVOT(I)                                              THE15350
            W(IP,K) = W(IP,K) / W(IPK,K)                                THE15360
            RATIO = -W(IP,K)                                            THE15370
            DO 20 J = KP1, N                                            THE15380
20          W(IP, J) = RATIO * W(IPK, J) + W(IP, J)                     THE15390
      IF (W(IP,N) .NE. 0.0)   RETURN                                    THE15400
999   IFLAG = 2                                                         THE15410
      RETURN                                                            THE15420
      END                                                               THE15430
                                                                        THE15440
      SUBROUTINE SUBST (W, B, IPS, IPIVOT, N, IDIM)                     THE15450
        IMPLICIT REAL*8(A-H,O-Z)                                        THE15460
      DIMENSION W(7,7), B(*), X(180), IPIVOT(*)                         THE15470
      COMMON /PDV/ X                                                    THE15480
      IF (N .GT. 1)  GO TO 10                                           THE15490
         X(IPS) = B(1) / W(1,1)                                         THE15500
         RETURN                                                         THE15510
C                                                                       THE15520
10    IP = IPIVOT (1)                                                   THE15530
      X(IPS) = B(IP)                                                    THE15540
      DO 15 K = 2, N                                                    THE15550
         IP = IPIVOT(K)                                                 THE15560
         KM1 = K - 1                                                    THE15570
         SUM = 0.0                                                      THE15580
         DO 14 J = 1, KM1                                               THE15590
14       SUM = W(IP, J) * X(IPS+J-1) + SUM                              THE15600
15    X(IPS+K-1) = B(IP) - SUM                                          THE15610
C                                                                       THE15620
      X(IPS+N-1) = X(IPS+N-1) / W(IP, N)                                THE15630
      K = N                                                             THE15640
      DO 20 NP1MK = 2, N                                                THE15650
         KP1 = K                                                        THE15660
         K = K - 1                                                      THE15670
         IP = IPIVOT(K)                                                 THE15680
         SUM = 0.0                                                      THE15690
         DO 19 J = KP1, N                                               THE15700
19       SUM = W(IP,J) * X(IPS+J-1) + SUM                               THE15710
20    X(IPS+K-1) = (X(IPS+K-1) - SUM) / W(IP, K)                        THE15720
      RETURN                                                            THE15730
      END                                                               THE15740
                                                                        THE15750
        FUNCTION AVG(A,B)                                               THE15760
        IMPLICIT REAL*8 (A-H,O-Z)                                       THE15770
        DIMENSION A(1),B(1)                                             THE15780
        AVG=(A(1)+B(1))/2.0                                             THE15790
        RETURN                                                          THE15800
        END                                                             THE15810
                                                                        THE15820
                                                                        THE15830
                                                                  
      SUBROUTINE FOFX(I,X,B,Y,Z)                                     
      
c     dmm 20000530  reset X(*) to be simple DP variable X
        IMPLICIT REAL*8(A-H,O-Z)                                       
      DOUBLE PRECISION B,CONST,XIN,THETA,CV,TERM(3),SUM, X              
      DIMENSION B(*),CONST(3)
c     dmm 20000530 end changes

      COMMON /CONST/ CONST,NTERMS                                       THE15900
        XIN=X                                                        
C       write(*,*)x(1),y                                                THE15920
        Z=0.0                                                           THE15930
        DO 100 II=1,NTERMS                                              THE15940
        THETA=0.0                                                       THE15950
        IF(NTERMS.EQ.1)THEN                                             THE15960
                CALL CVVIB(B(2),XIN,CV,THETA)                           THE15970
        ELSE                                                            THE15980
                CALL CVVIB(B(NTERMS+II-1),XIN,CV,THETA)                 THE15990
        ENDIF                                                           THE16000
        IF(II.NE.NTERMS.OR.NTERMS.EQ.1)THEN                             THE16010
                TERM(II)=B(II)*CV                                       THE16020
        ELSE IF(II.EQ.NTERMS.AND.NTERMS.GT.1)THEN                       THE16030
                SUM=0.0                                                 THE16040
                DO 234 IK=1,NTERMS-1                                    THE16050
234             SUM=SUM+B(IK)                                           THE16060
                TERM(II)=(CONST(2)-SUM)*CV                              THE16070
        ENDIF                                                           THE16080
        Z=Z+TERM(II)                                                    THE16090
100     CONTINUE                                                        THE16100
        IF(CONST(3).EQ.0.)THEN                                          THE16110
                                                                        THE16120
C               non-linear molecule                                     THE16130
                Z=(Z+4.0)*1.987*CONST(1)                                THE16140
        ELSE                                                            THE16150
C               linear molecule                                         THE16160
                Z=(Z+(7./2.))*1.987*CONST(1)                            THE16170
        ENDIF                                                           THE16180
C       write(*,*)'xin ',xin,' z= ',z                                   THE16190
C       write(*,*)'theta= ',theta                                       THE16200
      RETURN                                                            THE16210
      END                                                               THE16220
                                                                        THE16230
        SUBROUTINE PD(NDATA,X,NPARM,B,Y,ZZ,Z,PDERIV,NORMAL,WORK,POIN01) THE16240
                                                                        THE16250
        IMPLICIT REAL*8(A-H,O-Z)                                        THE16260
        DOUBLE PRECISION NORMAL(7,7),WORK(*),B(*),DEL,POIN01,           THE16270
     +CONST(3),U,U1,DU,DU1,V,V1,DV,DV1,ALPHA,PRD,TERM1,TERM2,TERM3,     THE16280
     +TERM4,HP,KB,C                                                     THE16290
        DIMENSION X(1,30),Y(*),PDERIV(6,30),Z(*)                        THE16300
        INTEGER NDATA,NPARM                                             THE16310
        COMMON/CONST/CONST,NTERMS                                       THE16320
        HP=6.62517E-34                                                  THE16330
        KB=1.3804E-23                                                   THE16340
        C=2.998E+10                                                     THE16350
                                                                        THE16360
C             INITIALIZE NORMAL TO 0.0,  ONLY THE PART OF               THE16370
C             THE ARRAY THAT WILL BE USED                               THE16380
C                                                                       THE16390
         DO 101 I = 1,NPARM+1                                           THE16400
         DO 100 J = I,NPARM+1                                           THE16410
         NORMAL(J,I) = 0.0                                              THE16420
100     CONTINUE                                                        THE16430
101     CONTINUE                                                        THE16440
C********************************************************************   THE16450
C          REPLACE THIS SECTION WITH EXPLICIT JACOBIAN                  THE16460
C                                                                       THE16470
C                                                                       THE16480
C             CALCULATE THE PARTIAL DERIVATIVES OF THE FUNCTION         THE16490
C             VALUES WITH RESPECT TO THE PARAMETER VALUES FOR           THE16500
C             EACH SET OF DATA POINTS.  THESE DERIVATIVES WILL BE       THE16510
C             STORED IN PDERIV.                                         THE16520
C                                                                       THE16530
         DO 200 J = 1,NDATA                                             THE16540
        ALPHA=(HP*C)/(KB*X(1,J))                                        THE16550
        U2=EXP(ALPHA*B(5))                                              THE16560
        TERM3=(ALPHA*B(5))**2                                           THE16570
        TERM4=(U2/((U2-1.0)**2))                                        THE16580
        PRD=TERM3*TERM4                                                 THE16590
        DO 150 I=1,2                                                    THE16600
        U1=EXP(ALPHA*B(I+2))                                            THE16610
        TERM1=(ALPHA*B(I+2))**2                                         THE16620
                                                                        THE16630
        TERM2=(U1/((U1-1.0)**2))                                        THE16640
        PDERIV(I,J)=(TERM1*TERM2-PRD)*1.987*CONST(1)                    THE16650
150     CONTINUE                                                        THE16660
        DO 1000 I=3,5                                                   THE16670
        U=(ALPHA*B(I))**2                                               THE16680
        DU=2.0*(ALPHA**2)*B(I)                                          THE16690
        U1=EXP(ALPHA*B(I))                                              THE16700
        DU1=ALPHA*U1                                                    THE16710
        V1=(U1-1.0)**2                                                  THE16720
        DV1=(2.0*(U1-1.0)*DU1)                                          THE16730
        V=U1/V1                                                         THE16740
        DV=((1.0/V1)*DU1)-((U1/(V1**2))*DV1)                            THE16750
        PDERIV(I,J)=(U*DV)+(V*DU)                                       THE16760
        IF(I.LT.5)THEN                                                  THE16770
                PDERIV(I,J)=PDERIV(I,J)*B(I-2)*1.987*CONST(1)           THE16780
        ELSE                                                            THE16790
                PDERIV(I,J)=PDERIV(I,J)*(CONST(2)-(B(1)+B(2)))*1.987*   THE16800
     +CONST(1)                                                          THE16810
        ENDIF                                                           THE16820
1000    CONTINUE                                                        THE16830
                                                                        THE16840
C                PDERIV(J,I) = (ZZ - Z(I)) / DEL                        THE16850
200     CONTINUE                                                        THE16860
C                                                                       THE16870
C                                                                       THE16880
C********************************************************************   THE16890
C                                                                       THE16900
C             COMPUTE THE  NORMAL EQUATIONS                             THE16910
C                                                                       THE16920
      DO 240 J1 = 1,NDATA                                               THE16930
         DO 230 J2 = 1,NPARM                                            THE16940
230         WORK(J2) = PDERIV(J2,J1)                                    THE16950
         WORK(NPARM+1) = Y(J1) - Z(J1)                                  THE16960
         DO 240 J3 = 1,NPARM+1                                          THE16970
            DO 240 J4 = J3,NPARM+1                                      THE16980
240            NORMAL(J4,J3) = NORMAL(J4,J3) + WORK(J3) * WORK(J4)      THE16990
        RETURN                                                          THE17000
        END                                                             THE17010
                                                                        THE17020
                                                                        THE17030




      SUBROUTINE READIN(SPEC,INPUT,HF298,S298,CPL,T,REF,REF3,IE,NE,  
     +     PHASE,NROT,NT,FUNITS,CFLAG,ICOUNT,COMMAN)    
      IMPLICIT REAL*8(A-H,O-Z)                                       
      CHARACTER*16 SPEC,INPUT,FUNITS                                 
      DIMENSION CPL(30),T(30)                                        
      CHARACTER*1 PHASE,CFLAG,BLANK,COMMA                            
      CHARACTER*6 REF3                                               
      CHARACTER*2 IE(4)                                              
      CHARACTER*4 REF(2),INAME                                       
      CHARACTER*70 TEXT
    CHARACTER*(*) COMMAN                                              
      REAL*4 RN                                                      
      INTEGER NE(4),IND(70)                                          
      
      BLANK=' '                                                      
      COMMA=','                                                      
      
      DO 1 I=1,4                                                     
         NE(I)=0                                                        
         IE(I)=' '                                                      
 1    CONTINUE                                                       
      
C**********************************************************************
C***********LST FORMAT INPUT FILE  **********************************
C**********************************************************************
      
      IF(INPUT.EQ.'LST')THEN                                         
         T(1)=300.                                              
         T(2)=400.                                              
         T(3)=500.                                              
         T(4)=600.                                              
         T(5)=800.                                              
         T(6)=1000.                                             
         T(7)=1500.                                             
c**** Changed this wretched formatted read to be XMG-friendly DMM
c     This whole program is an example of what not to do 
c   Added COMMAN (WILHOI,EXP,or HOE) to list of inputs - PEY 

       READ(1,*) COMMAN
         READ(1,168) SPEC(:16)
 168     format(A16)
         READ(1,*) HF298
         READ(1,*) S298
         READ(1,*) CPL(1)
         READ(1,*) CPL(2)
         READ(1,*) CPL(3)
         READ(1,*) CPL(4)
         READ(1,*) CPL(5)
         READ(1,*) CPL(6)
         READ(1,*) CPL(7)
         
         DO 170 I=1,4
            READ(1,169) IE(I)
 169        format(A2)
            READ(1,*) NE(I)
 170     continue

         READ(1,*) PHASE
         READ(1,*) NROT

c 170        FORMAT(A1,A10,7X,F8.2,1X,F8.2,3X,7(F7.3,2X),1X,2A4,A6,7X
c     $           ,4(A2,I3,1X),A1,1x,I2,A2,I3,1x,A1,I2)
         
c     IF(INAME(:3).EQ.'END') THEN                            
c     INPUT='DONE'                                   
c     RETURN                                         

         IF(CPL(7).EQ.0.)THEN                                   
            NT=6                                           
         ELSE                                                   
            NT=7                                           
         ENDIF                                                  
C---------------------------------------------------INTERACTIVE INPUT 
      ELSEIF(INPUT.EQ.'NONLST')THEN                                  
         FUNITS='KCAL'                                          
         read(1,4321,end=901) cflag, spec(:16)
 4321    format(a1,a16)
c*********terminate input file with END starting in column 1
         if(cflag.eq.'E') then
            input='DONE'
            write(*,*)  'reached end of input file'
            return
         endif
         read(1,4324) (ie(i),ne(i), i=1,4)
 4324    format(4(a2,i3,1x))
         read(1,4322) ref(1), ref(2), ref3
 4322    format(2a4,a6)
         read(1,4323) phase,nrot
 4323    format(a1,i2)
         read(1,*) hf298,s298,nt
         read(1,*) (t(i),cpl(i), i=1,nt)
 2112    IF(Nt.LT.5)THEN                                         
            WRITE(*,*)' ERROR : TOO FEW DATA POINTS!'      
            WRITE(*,*)'SPECIES NOT PROCESSED'              
            stop                                      
         elseif(nt.gt.29) then
            write(*,*)'  error:  too many data points'  
            write(*,*)'  29 max allowed'    
            WRITE(*,*)'SPECIES NOT PROCESSED'              
            stop                                       
         ENDIF                                                   
      else 
         write(*,*) 'input option not specified'
         stop
      ENDIF                                                          
      
C**************************************ALL INPUT                     
      IF(CFLAG.EQ.' ')THEN                                           
         WRITE(25,6321) SPEC
 6321    format(//,A16)                            
      elseif(cflag.eq.'*') then
         WRITE(25,6322)SPEC
 6322    format(//,A16,'fit skipped, NASA polynomials generated')
     $        
      elseif(cflag.eq.'@') then
         WRITE(25,6324)SPEC
 6324    format(//,A16,'  no NASA polynomials generated')
     $        
      ENDIF                                                          
C---  CONVERTING HF298 TO CAL/MOLE OR J/MOLE
c     
      IF(FUNITS.EQ.'KCAL')THEN                                       
         HF298=HF298*1000.                                        
      ELSEIF(FUNITS.EQ.'KJ')THEN                                     
         HF298=HF298*1000.0/4.184                               
         S298=S298/4.184                                        
      ENDIF                                                          
      RETURN                                                         
 901  CONTINUE                                                       
 900  INPUT='DONE'                                                   
      RETURN                                                         
      END                                                            
                                                                      



                                                                        THE20980
                                                                        THE20990
        SUBROUTINE XCHEX(SPEC,LNG,SS,NX,ICNT,IADD,ISAME)                THE21000
        CHARACTER*(*) SS(*),SPEC                                        THE21010
        INTEGER NX(*),IADD,ISAME,ICNT                                   THE21020
        IADD=1                                                          THE21030
        DO 100 I=1,ICNT                                                 THE21040
        IF(SPEC(:LNG).EQ.SS(I)(:LNG))THEN                               THE21050
                ISAME=1                                                 THE21060
                IADD=I                                                  THE21070
                RETURN                                                  THE21080
        ENDIF                                                           THE21090
100     CONTINUE                                                        THE21100
        RETURN                                                          THE21110
        END
          SUBROUTINE DBLNK(S1,LST,CCSW)                                 THS00010
C ********************************************************************* THS00020
C *****   5/25/90 MODIFIED TO RUN ON FLIC.  ALL OPEN AND CLOSES   ***** THS00030
C *****   WERE REMOVED AND REPLACED BY FILEDEFS.  IN ADDITION,    ***** THS00040
C *****   ROUTINES TO CAPITOLIZE STRINGS AND MANIPULATE FILENAMES ***** THS00050
C *****   WERE REMOVED, AND THE CLS SUBROUTINE WAS DELETED. (JCB) ***** THS00060
C ********************************************************************* THS00070
          CHARACTER*(*) S1                                              THS00080
          INTEGER LS1,CCSW,FLG                                          THS00090
          FLG=0                                                         THS00100
          LST=LEN(S1)                                                   THS00110
          NNB=0                                                         THS00120
          KOUNT=0                                                       THS00130
C         WRITE(*,*)'ENTER STRING'                                      THS00140
C         READ(*,3)S1                                                   THS00150
3         FORMAT(A70)                                                   THS00160
          DO 100 I=1,LST                                                THS00170
          IF(S1(I:I).NE.' ')NNB=NNB+1                                   THS00180
          IF(S1(I:LST).EQ.' ')THEN                                      THS00190
                LST=I-1                                                 THS00200
                GOTO 450                                                THS00210
          ENDIF                                                         THS00220
C         WRITE(*,*)'NNB= ',NNB                                         THS00230
100       CONTINUE                                                      THS00240
450       CONTINUE                                                      THS00250
400       I=1                                                           THS00260
          KOUNT=KOUNT+1                                                 THS00270
300       CONTINUE                                                      THS00280
          IF(I.GT.LST)GOTO 200                                          THS00290
          IF(S1(I:I).EQ.' '.AND.(CCSW.EQ.0.OR.I.EQ.1))THEN              THS00300
                S1(I:LST-1)=S1(I+1:LST)                                 THS00310
                S1(LST:LST)=' '                                         THS00320
                LST=LST-1                                               THS00330
          ENDIF                                                         THS00340
          IP1=I+1                                                       THS00350
          IF(S1(I:I).EQ.' '.AND.S1(IP1:IP1).EQ.' '.AND.CCSW.EQ.1)THEN   THS00360
                S1(I:LST-1)=S1(I+1:LST)                                 THS00370
                S1(LST:LST)=' '                                         THS00380
                LST=LST-1                                               THS00390
          ENDIF                                                         THS00400
          I=I+1                                                         THS00410
          GOTO 300                                                      THS00420
200       CONTINUE                                                      THS00430
          IF(KOUNT.GT.LST)GOTO 500                                      THS00440
          IF(LST.GT.NNB)GOTO 400                                        THS00450
500       CONTINUE                                                      THS00460
C         WRITE(*,*)S1                                                  THS00470
C         WRITE(*,*)'LST= ',LST                                         THS00480
          RETURN                                                        THS00490
          END                                                           THS00500
                                                                        THS00510
        SUBROUTINE INDXST(S1,S2,N,INDX,LS1,LS2,FLG)                     THS00750
        CHARACTER S1*(*),S2*(*)                                         THS00760
        INTEGER N,INDX(*),LS1,LS2,FLG                                   THS00770
C***********************************************************************THS00780
C                                                                       THS00790
C  SUBROUTINE INDXST  :   INDeX STring                                  THS00800
C                                                                       THS00810
C                                                                       THS00820
C     The purpose of this subroutine is to search string S1  (arbitrary THS00830
C    LENGTH) FOR ALL OCCURRENCES OF SUBSTRING S2.   S2 MUST BE SHORTER  THS00840
C    than S1.  Returned are the number of occurrences of S2 in S1;  the THS00850
C    vector of starting addresses of each occurrence; string lengths forTHS00860
C    both S1 and S2 if not supplied; and an error flag.                 THS00870
C                                                                       THS00880
C    Variable definitions:                                              THS00890
C      S1 - input string to be searched (character length defined by    THS00900
C           calling routine).                                           THS00910
C      S2 - target substring; length determined by calling routine.     THS00920
C      N  - number of occurences of substring S2 in S1.                 THS00930
C  INDX(*)- vector giving the starting address for each occurrence of   THS00940
C           S2 in S1.                                                   THS00950
C      LS1- length (non-blank) of S1; if not supplied it is determined  THS00960
C           and returned to the calling routine.                        THS00970
C      LS2- length ( = 1 for a single blank ) of target substring S2;   THS00980
C           if not supplied it is determined and returned to the callingTHS00990
C           routine.                                                    THS01000
C      FLG- error flag.  FLG = 0 if no error encountered                THS01010
C                        FLG < 0 if LS2 > LS1 or other error occurrs    THS01020
C                                                                       THS01030
C                                                                       THS01040
C***********************************************************************THS01050
        LMAX=LEN(S1)                                                    THS01060
        DO 10 I=1,LMAX                                                  THS01070
10      INDX(I)=0                                                       THS01080
        N=0                                                             THS01090
        FLG=0                                                           THS01100
C       DETERMINE STRING LENGTH IF NOT PASSED TO THIS ROUTINE           THS01110
        IF(LS1.EQ.0.AND.LS2.EQ.0)THEN                                   THS01120
                DO 100 I=1,LMAX                                         THS01130
                IP1=I+1                                                 THS01140
                IF(I.GT.1)THEN                                          THS01150
                        IM1=I-1                                         THS01160
                ELSE                                                    THS01170
                        IM1=1                                           THS01180
                ENDIF                                                   THS01190
                IF(LS1.EQ.0.AND.S1(I:I).EQ.' '.AND.S1(IP1:IP1).EQ.' ')  THS01200
     +          LS1=IM1                                                 THS01210
                IF(LS2.EQ.0.AND.S2(I:I).EQ.' '.AND.S2(IP1:IP1).EQ.' ')  THS01220
     +          LS2=IM1                                                 THS01230
100             CONTINUE                                                THS01240
        ENDIF                                                           THS01250
150     CONTINUE                                                        THS01260
        IF(LS1.LT.LS2)THEN                                              THS01270
                FLG=-10                                                 THS01280
                RETURN                                                  THS01290
        ENDIF                                                           THS01300
        IFOUND=0                                                        THS01310
        DO 200 I=1,LS1                                                  THS01320
        I2=I+LS2-1                                                      THS01330
        IF((I2-I).GT.LS2)GOTO 200                                       THS01340
        IF(S1(I:I2).EQ.S2(:LS2))THEN                                    THS01350
                IFOUND=IFOUND+1                                         THS01360
                INDX(IFOUND)=I                                          THS01370
        ENDIF                                                           THS01380
200     CONTINUE                                                        THS01390
        N=IFOUND                                                        THS01400
        RETURN                                                          THS01410
        END                                                             THS01420
                                                                        THS01430
        SUBROUTINE INQNBK(S1,IANS,NUMLEN,LEN,IBL)                       THS01440
C********************************************************************   THS01450
C                                                                   *   THS01460
C               written by Ed Ritter : 11/20/88                     *   THS01470
C                                                                   *   THS01480
C********************************************************************   THS01490
C                                                                   *   THS01500
C       THIS SUBROUTINE SCANS STRING S1 FOR A LEADING NUMBER        *   THS01510
C       The search is carried out from left to right.               *   THS01520
C       The first non-numeric character encountered ends            *   THS01530
C        the numeric substring.                                     *   THS01540
C                                                                   *   THS01550
C       return codes :IANS                                          *   THS01560
C                                                                   *   THS01570
C       IANS = 0  : the string is non-numeric                       *   THS01580
C            = 1  : string S1 is an integer                         *   THS01590
C            = 2  : S1 is a decimal                                 *   THS01600
C            = 3  : S1 has a leading integer of length NUMLEN       *   THS01610
C            = 4  : S1 has a leading decimal of length NUMLEN.      *   THS01620
                                                                        THS01630
C                                                                   *   THS01640
C       NUMLEN    : the length of the leading numeric portion of S1 *   THS01650
C       LEN       : the length of the S1 without blanks             *   THS01660
C       IBL       : allows the calling routine to control blank     *   THS01670
C                   removal in subroutine DEBLNK                    *   THS01680
C                   IBL=0 : leave no blanks in S1                   *   THS01690
C                   IBL=1 : leave one remaining blank in each blank *   THS01700
C                           field                                   *   THS01710
C                                                                   *   THS01720
C********************************************************************   THS01730
        CHARACTER*(*) S1                                                THS01740
        INTEGER IANS,NUMLEN,LEN,ICTRL,IDEC,IBL                          THS01750
C                                                                       THS01760
C  ICTRL <> 0, when a non-numeric character has been encountered in S1. THS01770
C  IDEC <> 0, if a decimal point is encountered in a numeric field;     THS01780
C             in this case IDEC = address of decimal point.             THS01790
                                                                        THS01800
10      CONTINUE                                                        THS01810
                                                                        THS01820
        NUMLEN=0                                                        THS01830
        ICTRL=0                                                         THS01840
        IDEC=0                                                          THS01850
        IANS=0                                                          THS01860
                                                                        THS01870
                                                                        THS01880
C  inserted for test                                                    THS01890
C       write(*,*)'enter string s1'                                     THS01900
C       read(*,'(a70)')s1                                               THS01910
C       ibl=0                                                           THS01920
C  END OF INSERT#1 TEST                                                 THS01930
                                                                        THS01940
                                                                        THS01950
        CALL DBLNK (S1,LEN,IBL)                                         THS01960
        DO 100 J=1,LEN                                                  THS01970
        I=LEN+1-J                                                       THS01980
        IF((ICHAR(S1(I:I)).GE.48.AND.ICHAR(S1(I:I)).LE.57).OR.          THS01990
     +(S1(I:I).EQ.'.'.AND.IDEC.EQ.0))THEN                               THS02000
                IF(ICTRL.EQ.0)THEN                                      THS02010
                        NUMLEN=NUMLEN+1                                 THS02020
                        IF(S1(I:I).EQ.'.')IDEC=I                        THS02030
                ENDIF                                                   THS02040
        ELSE                                                            THS02050
                ICTRL=1                                                 THS02060
        ENDIF                                                           THS02070
100     CONTINUE                                                        THS02080
        IF(NUMLEN.EQ.0)THEN                                             THS02090
                IANS=0                                                  THS02100
        ELSEIF(NUMLEN.EQ.LEN)THEN                                       THS02110
                IF(IDEC.EQ.0)THEN                                       THS02120
                        IANS=1                                          THS02130
                ELSEIF(IDEC.NE.0)THEN                                   THS02140
                        IANS=2                                          THS02150
                ENDIF                                                   THS02160
        ELSEIF(NUMLEN.LT.LEN)THEN                                       THS02170
                IF(IDEC.EQ.0)THEN                                       THS02180
                        IANS=3                                          THS02190
                ELSEIF(IDEC.NE.0)THEN                                   THS02200
                        IANS=4                                          THS02210
                ENDIF                                                   THS02220
        ENDIF                                                           THS02230
                                                                        THS02240
C       insert for test                                                 THS02250
C       write(*,*)'S1 = ',s1                                            THS02260
C       write(*,*)'IANS = ',ians                                        THS02270
C       write(*,*)'NUMLEN = ',numlen                                    THS02280
C       write(*,*)'IDEC = ',idec                                        THS02290
C       write(*,*)'LEN = ',len                                          THS02300
C       pause ' '                                                       THS02310
C       goto 10                                                         THS02320
C       nd of insert#2                                                  THS02330
                                                                        THS02340
99999   CONTINUE                                                        THS02350
        RETURN                                                          THS02360
C       STOP ' '                                                        THS02370
        END                                                             THS02380
                                                                        THS02390
                                                                        THS02400
      SUBROUTINE LINES(D,I)                                             THS02410
      INTEGER I,D                                                       THS02420
      IF(D.EQ.0) THEN                                                   THS02430
         DO 2020 JJ=1,I                                                 THS02440
2020     WRITE (*,*)' '                                                 THS02450
         RETURN                                                         THS02460
      ENDIF                                                             THS02470
         DO 2021 JJ=1,I                                                 THS02480
2021     WRITE (D,*)' '                                                 THS02490
      RETURN                                                            THS02500
      END                                                               THS02510
                                                                        THS02520
      SUBROUTINE LRDCHR(LS1)                                            THS02530
C THIS SUBROUTINE RETURNS THE SAME STRING IN CAPITOL LETTERS            THS02540
C NUMERIC EXPRESSIONS REMAIN THE SAME                                   THS02550
      CHARACTER*70 LS                                                   THS02560
        CHARACTER*70 LS1,DUMMY                                          THS02570
C       LS1=' '                                                         THS02580
C       LS1=LS(:14)                                                     THS02590
        CALL DEBLNK(LS1,DUMMY,LEN,1)                                    THS02600
C       LS(:14)=LS1(:14)                                                THS02610
C      WRITE(*,*)'READCRin ',LS                                         THS02620
      I12=1                                                             THS02630
1400  IF(I12.GT.LEN)GOTO 1450                                           THS02640
      IF(ICHAR(LS1(I12:I12)).GE.97)THEN                                 THS02650
        N28=ICHAR(LS1(I12:I12))-32                                      THS02660
        LS1(I12:I12)=CHAR(N28)                                          THS02670
      ENDIF                                                             THS02680
      I12=I12+1                                                         THS02690
      IF (I12.GE.(LEN+1))GOTO 1450                                      THS02700
                                                                        THS02710
                                                                        THS02720
      GOTO 1400                                                         THS02730
1450  CONTINUE                                                          THS02740
C      WRITE(*,*)'READCROUT ',LS                                        THS02750
      RETURN                                                            THS02760
      END                                                               THS02770
                                                                        THS02780
                                                                        THS02790
C      SUBROUTINE READCR(LS)                                            THS02800
CC THIS SUBROUTINE RETURNS THE SAME STRING IN CAPITOL LETTERS           THS02810
CC NUMERIC EXPRESSIONS REMAIN THE SAME                                  THS02820
C      CHARACTER*(*) LS                                                 THS02830
C      CALL DBLNK(LS,LEN,1)                                             THS02840
CC      WRITE(*,*)'READCRIN ',LS                                        THS02850
C      I12=1                                                            THS02860
C1400  IF(I12.GT.LEN)GOTO 1450                                          THS02870
C      IF(ICHAR(LS(I12:I12)).GE.97)THEN                                 THS02880
C        N28=ICHAR(LS(I12:I12))-32                                      THS02890
C        LS(I12:I12)=CHAR(N28)                                          THS02900
C      ENDIF                                                            THS02910
C      I12=I12+1                                                        THS02920
C      IF (I12.GT.LEN)GOTO 1450                                         THS02930
C      GOTO 1400                                                        THS02940
C1450  CONTINUE                                                         THS02950
CC      WRITE(*,*)'READCROUT ',LS                                       THS02960
C      RETURN                                                           THS02970
C      END                                                              THS02980
                                                                        THS02990
        SUBROUTINE FORMIN(FORM,LEN1,IE,NE,NATOMT,IENTRY)                THS03000
        CHARACTER*2 ELEM(50)                                            THS03010
        CHARACTER*70 FORM,DUMMY                                         THS03020
        CHARACTER*4 IE(4),IES(4),IE1(4,400),B,C,TMP                     THS03030
        INTEGER NE(4),M,M1,N(4),L,ITEMP,LEN1,ID(70),NOC1(70),           THS03040
     +  NOC2(70),IN2C(35,4),IN1C(70,4),I2CH(35),I1CH(50),N2CH,N1CH,     THS03050
     +  IFM(4),IG(35),NG,INUM(70),NUMBER(70),NNUM,NES(4),NLP,IDLP(70),  THS03060
     +  NRP,IDRP(70),NSLASH,IDSLSH(70),IFOUND(4),NE1(4,400),NNZERO,     THS03070
     +  NELT(20),IELEM(10),NEL(20),IEL(4),ISSQ(20)                      THS03080
        COMMON/CONFG1/ELEM                                              THS03090
        COMMON/CONFG2/NELEM,N1CH,N2CH,I1CH,I2CH                         THS03100
        DO 234 I=1,NELEM                                                THS03110
234     NELT(I)=0                                                       THS03120
        DO 123 I=1,4                                                    THS03130
        NES(I)=0                                                        THS03140
        NE(I)=0                                                         THS03150
        IES(I)=' '                                                      THS03160
123     IE(I)=' '                                                       THS03170
        DO 345 I=1,LEN1                                                 THS03180
        IF(FORM(I:I).EQ.''.OR.FORM(I:I).EQ.'{')FORM(I:I)='('           THS03190
        IF(FORM(I:I).EQ.'|'.OR.FORM(I:I).EQ.'}')FORM(I:I)=')'           THS03200
345     CONTINUE                                                        THS03210
        NSUBS=0                                                         THS03220
        NATOMT=0                                                        THS03230
        NATOM1=0                                                        THS03240
        IF(LEN1.EQ.1)THEN                                               THS03250
                IE(1)(:1)=FORM(:1)                                      THS03260
                NE(1)=1                                                 THS03270
                                                                        THS03280
                NATOM=1                                                 THS03290
                GOTO 71                                                 THS03300
C*                      XFER TO OUTPUT/RETURN BLOCK                     THS03310
        ENDIF                                                           THS03320
C********************************************************************** THS03330
C*              SEARCH FOR NUMERIC CHARACTERS                        *  THS03340
C********************************************************************** THS03350
        CALL INDXNU(FORM,LEN1,NNUM,NUMBER,INUM,KFLG)                    THS03360
C********************************************************************** THS03370
C*                                                                      THS03380
C********************************************************************** THS03390
C*    SEARCH FOR PARENTHESES                                            THS03400
C****************************                                           THS03410
        CALL INDXST(FORM,'(',NLP,IDLP,LEN1,1,KFLG)                      THS03420
        CALL INDXST(FORM,')',NRP,IDRP,LEN1,1,KFLG)                      THS03430
        CALL INDXST(FORM,'/',NSLASH,IDSLSH,LEN1,1,KFLG)                 THS03440
C********************************************************************** THS03450
        IF(NLP.EQ.0.AND.NRP.EQ.0.AND.NSLASH.EQ.0)THEN                   THS03460
                CALL FINTER(FORM,LEN1,IE,NE,NATOMT,IENTRY)              THS03470
                GOTO 71                                                 THS03480
C*                      XFER TO OUTPUT/RETURN BLOCK                     THS03490
        ELSE IF(NSLASH.NE.0.AND.(NLP.NE.0.OR.NRP.NE.0))THEN             THS03500
C*              CAN'T MIX STRING OPTIONS                                THS03510
                IENTRY=-3                                               THS03520
                GOTO 71                                                 THS03530
C*                      XFER TO OUTPUT/RETURN BLOCK                     THS03540
        ELSE IF(NLP.NE.NRP)THEN                                         THS03550
C*               UNBALANCED PARENTHESES                                 THS03560
                IENTRY=-2                                               THS03570
                GOTO 71                                                 THS03580
C*                      XFER TO OUTPUT/RETURN BLOCK                     THS03590
        ELSE IF(NSLASH.NE.0)THEN                                        THS03600
C*               DECOMPOSE SUBSTRINGS AND PASS PIECEWISE                THS03610
                IF(NSLASH.EQ.1.AND.IDSLSH(1).EQ.1)THEN                  THS03620
                        IENTRY=-4                                       THS03630
                        GOTO 71                                         THS03640
C*                      XFER TO OUTPUT/RETURN BLOCK                     THS03650
                ENDIF                                                   THS03660
                IF(IDSLSH(1).NE.1)THEN                                  THS03670
                        DUMMY=' '                                       THS03680
                        DUMMY=FORM(:IDSLSH(1)-1)                        THS03690
                        CALL SLEN(DUMMY,LENS)                           THS03700
                CALL FINTER(DUMMY,LENS,IES,NES,NATOM1,IENTRY)           THS03710
                        NATOMT=NATOMT+NATOM1                            THS03720
                        NSUBS=1                                         THS03730
                        DO 1001 IK=1,4                                  THS03740
                        IE1(IK,1)=IES(IK)                               THS03750
1001                    NE1(IK,1)=NES(IK)                               THS03760
                ENDIF                                                   THS03770
                DO 2001 JK=1,NSLASH-1                                   THS03780
                DUMMY=' '                                               THS03790
                IF(IDSLSH(JK)+1.GT.IDSLSH(JK+1)-1)GOTO 2001             THS03800
                NSUBS=NSUBS+1                                           THS03810
                DUMMY=FORM(IDSLSH(JK)+1:IDSLSH(JK+1)-1)                 THS03820
                CALL SLEN(DUMMY,LENS)                                   THS03830
                                                                        THS03840
                                                                        THS03850
                CALL FINTER(DUMMY,LENS,IES,NES,NATOM1,IENTRY)           THS03860
                NATOMT=NATOM1+NATOMT                                    THS03870
                DO 1002 IK1=1,4                                         THS03880
                IE1(IK1,1+JK)=IES(IK1)                                  THS03890
1002            NE1(IK1,1+JK)=NES(IK1)                                  THS03900
2001            CONTINUE                                                THS03910
                IF(LEN1.NE.IDSLSH(NSLASH))THEN                          THS03920
                DUMMY=' '                                               THS03930
                NSUBS=NSUBS+1                                           THS03940
                DUMMY=FORM(IDSLSH(NSLASH)+1:LEN1)                       THS03950
                CALL SLEN(DUMMY,LENS)                                   THS03960
                CALL FINTER(DUMMY,LENS,IES,NES,NATOM1,IENTRY)           THS03970
                NATOMT=NATOM1+NATOMT                                    THS03980
                DO 1003 IK1=1,4                                         THS03990
                IE1(IK1,NSUBS)=IES(IK1)                                 THS04000
1003            NE1(IK1,NSUBS)=NES(IK1)                                 THS04010
                ENDIF                                                   THS04020
                GOTO 81                                                 THS04030
C*               SUM ELEMENTS OVER ALL SUBSTRINGS                       THS04040
        ELSE IF(NLP.EQ.NRP)THEN                                         THS04050
C*               DECOMPOSE PARENTHETICAL SUBSTRINGS AND PASS            THS04060
C*              CALL FINTER(DUMMY,LENS,IES,NES,NATOM1,IENTRY)           THS04070
                DO 5101 J=1,NRP                                         THS04080
                DO 5201 I=J,NNUM                                        THS04090
                ISSQ(J)=1                                               THS04100
                IF(INUM(I).EQ.IDRP(J)+1)THEN                            THS04110
                        ISSQ(J)=NUMBER(I)                               THS04120
                        GOTO 5101                                       THS04130
                ENDIF                                                   THS04140
5201            CONTINUE                                                THS04150
5101            CONTINUE                                                THS04160
                DO 6301 I=1,NRP-1                                       THS04170
                IF(IDRP(I).GT.IDLP(I+1))THEN                            THS04180
                        WRITE(*,*)' TOO MANY LEVELS OF PARENTHESES'     THS04190
                        PAUSE ' hit  return to continue'                THS04200
                        IENTRY=-6                                       THS04210
                        GOTO 71                                         THS04220
                ENDIF                                                   THS04230
6301            CONTINUE                                                THS04240
                IF(IDRP(NRP).LT.IDLP(NLP))THEN                          THS04250
                        WRITE(*,*)' TOO MANY LEVELS OF PARENTHESES'     THS04260
                        PAUSE ' hit  return to continue'                THS04270
                        IENTRY=-6                                       THS04280
                        GOTO 71                                         THS04290
                ENDIF                                                   THS04300
                IF(IDLP(1).NE.1)THEN                                    THS04310
                        DUMMY=' '                                       THS04320
                        DUMMY=FORM(:IDLP(1)-1)                          THS04330
                        NSUBS=NSUBS+1                                   THS04340
                        CALL SLEN(DUMMY,LENS)                           THS04350
                CALL FINTER(DUMMY,LENS,IES,NES,NATOM1,IENTRY)           THS04360
                        NATOMT=NATOMT+NATOM1                            THS04370
                        DO 5001 IK=1,4                                  THS04380
                        IE1(IK,NSUBS)=IES(IK)                           THS04390
5001                    NE1(IK,NSUBS)=NES(IK)                           THS04400
                                                                        THS04410
                                                                        THS04420
C*      write(*,*)' substring # ',nsubs                                 THS04430
C*      write(*,*)dummy                                                 THS04440
C*      write(*,787)(ies(ik),nes(ik),ik=1,4)                            THS04450
787     FORMAT(1X,4(4X,A4,I3))                                          THS04460
                ENDIF                                                   THS04470
                DO 6201 J=1,NRP                                         THS04480
                NCHARN=1                                                THS04490
                DUMMY=' '                                               THS04500
                IF(ISSQ(J).GT.9)NCHARN=2                                THS04510
                IF(ISSQ(J).GT.99)NCHARN=3                               THS04520
                IF(ISSQ(J).GT.999)NCHARN=4                              THS04530
                IF(ISSQ(J).GT.9999)NCHARN=5                             THS04540
        IF(J.LT.NRP)THEN                                                THS04550
                IF(IDLP(J+1)-1.EQ.IDRP(J))GOTO 6201                     THS04560
                IF(IDRP(J)+NCHARN+1.LT.IDLP(J+1))THEN                   THS04570
                DUMMY=FORM(IDRP(J)+NCHARN+1:IDLP(J+1)-1)                THS04580
                GOTO 91                                                 THS04590
                ENDIF                                                   THS04600
        ENDIF                                                           THS04610
                                                                        THS04620
        GOTO 6201                                                       THS04630
91      CONTINUE                                                        THS04640
                CALL SLEN(DUMMY,LENS)                                   THS04650
                CALL FINTER(DUMMY,LENS,IES,NES,NATOM1,IENTRY)           THS04660
                NATOMT=NATOM1+NATOMT                                    THS04670
                NSUBS=NSUBS+1                                           THS04680
        DO 5301 IK=1,4                                                  THS04690
        IE1(IK,NSUBS)=IES(IK)                                           THS04700
5301    NE1(IK,NSUBS)=NES(IK)                                           THS04710
C*      write(*,*)' substring # ',nsubs                                 THS04720
C*      write(*,*)dummy                                                 THS04730
C*      write(*,787)(ies(ik),nes(ik),ik=1,4)                            THS04740
6201    CONTINUE                                                        THS04750
                IF(IDRP(NRP)+NCHARN+1.LE.LEN1)THEN                      THS04760
                        DUMMY=FORM(IDRP(NRP)+NCHARN+1:LEN1)             THS04770
                        CALL SLEN(DUMMY,LENS)                           THS04780
                CALL FINTER(DUMMY,LENS,IES,NES,NATOM1,IENTRY)           THS04790
                        NATOMT=NATOM1+NATOMT                            THS04800
                        NSUBS=NSUBS+1                                   THS04810
                        DO 5321 IK=1,4                                  THS04820
                        IE1(IK,NSUBS)=IES(IK)                           THS04830
5321                    NE1(IK,NSUBS)=NES(IK)                           THS04840
C*      write(*,*)' substring # ',nsubs                                 THS04850
C*      write(*,*)dummy                                                 THS04860
C*      write(*,787)(ies(ik),nes(ik),ik=1,4)                            THS04870
                                                                        THS04880
C               WRITE(*,*)' SUBSTRING MUST BE ENCLOSED IN PARENTHESES'  THS04890
C               PAUSE 'HIT RETURN TO CONTINUE'                          THS04900
C               IENTRY=-1                                               THS04910
C               GOTO 71                                                 THS04920
                ENDIF                                                   THS04930
                DO 6001 J=1,NRP                                         THS04940
                DUMMY=' '                                               THS04950
                DUMMY=FORM(IDLP(J)+1:IDRP(J)-1)                         THS04960
                CALL SLEN(DUMMY,LENS)                                   THS04970
                                                                        THS04980
                                                                        THS04990
                CALL FINTER(DUMMY,LENS,IES,NES,NATOM1,IENTRY)           THS05000
                NATOMT=(NATOM1*ISSQ(J))+NATOMT                          THS05010
                NSUBS=NSUBS+1                                           THS05020
                DO 5401 IK=1,4                                          THS05030
                IE1(IK,NSUBS)=IES(IK)                                   THS05040
5401            NE1(IK,NSUBS)=(NES(IK)*ISSQ(J))                         THS05050
C*      write(*,*)' substring # ',nsubs                                 THS05060
C*      write(*,*)dummy                                                 THS05070
C*      write(*,787)(ies(ik),nes(ik),ik=1,4)                            THS05080
                                                                        THS05090
6001            CONTINUE                                                THS05100
                GOTO 81                                                 THS05110
C*               SUM ELEMENTS OVER ALL SUBSTRINGS                       THS05120
        ENDIF                                                           THS05130
C******************************************************************     THS05140
C*       SUBSTRING ELEMENT SUMMATION BLOCK                              THS05150
C******************************************************************     THS05160
81              CONTINUE                                                THS05170
                DO 850 IK2=1,NSUBS                                      THS05180
C*      write(*,*)'subr natoms output for substr# ',ik2                 THS05190
                CALL NATOMS(IK2,IE1,NE1,IEL,NEL,NATOM)                  THS05200
                DO 950 I1=1,4                                           THS05210
950             NELT(IEL(I1))=NELT(IEL(I1))+NEL(IEL(I1))                THS05220
C*      WRITE(*,787)(ELEM(IEL(I1)),NELT(IEL(I1)),I1=1,4)                THS05230
850             CONTINUE                                                THS05240
                NNZERO=0                                                THS05250
                DO 750 I1=1,NELEM                                       THS05260
                IF(NELT(I1).NE.0)THEN                                   THS05270
                        NNZERO=NNZERO+1                                 THS05280
                        IELEM(NNZERO)=I1                                THS05290
                ENDIF                                                   THS05300
750             CONTINUE                                                THS05310
                IF(NNZERO.GT.4)THEN                                     THS05320
                        IENTRY=-1                                       THS05330
                        GOTO 71                                         THS05340
C*                      XFER TO OUTPUT/RETURN BLOCK                     THS05350
                ENDIF                                                   THS05360
                DO 650 I1=1,NNZERO                                      THS05370
                NE(I1)=NELT(IELEM(I1))                                  THS05380
                IE(I1)=ELEM(IELEM(I1))                                  THS05390
650             CONTINUE                                                THS05400
C******************************************************************     THS05410
C*      OUTPUT/RETURN BLOCK                                             THS05420
C******************************************************************     THS05430
71      CONTINUE                                                        THS05440
        CALL ORDERL(NE,IE)                                              THS05450
        IF(IENTRY.EQ.0)THEN                                             THS05460
                WRITE(*,*)'FORMULA ENTERED :'                           THS05470
                WRITE(*,*)FORM                                          THS05480
                CALL LINES(0,2)                                         THS05490
                WRITE(*,*)'  ELEMENTAL COMPOSITION:'                    THS05500
                CALL LINES(0,1)                                         THS05510
                DO 23 I=1,4                                             THS05520
23              WRITE(*,*)IE(I),NE(I)                                   THS05530
        ENDIF                                                           THS05540
        RETURN                                                          THS05550
        END                                                             THS05560
        SUBROUTINE INDXNU(S1,LEN1,N,NUM,IND,KFLG)                       THS05570
        CHARACTER*(*) S1,S2*14                                          THS05580
        INTEGER N,NUM(132),IND(132),KFLG,LEN1,ITMP(132),NOC             THS05590
70      FORMAT(A)                                                       THS05600
        KFLG=0                                                          THS05610
        N=0                                                             THS05620
        LSTR=LEN(S1)                                                    THS05630
        DO 10 I=1,LSTR                                                  THS05640
10      IND(I)=0                                                        THS05650
        NOC=0                                                           THS05660
        DO 200 I=1,LEN1                                                 THS05670
        IF(ICHAR(S1(I:I)).LE.57.AND.ICHAR(S1(I:I)).GE.48)THEN           THS05680
                NOC=NOC+1                                               THS05690
                ITMP(NOC)=I                                             THS05700
        ENDIF                                                           THS05710
200     CONTINUE                                                        THS05720
        IF(NOC.EQ.0)THEN                                                THS05730
                KFLG=1                                                  THS05740
                RETURN                                                  THS05750
        ENDIF                                                           THS05760
        IF(NOC.EQ.1)THEN                                                THS05770
                N=1                                                     THS05780
                S2=S1(ITMP(1):ITMP(1))                                  THS05790
                CALL READNU(S2,RN,NUM(1),KFLG)                          THS05800
C*              CALL READNU(S1(ITMP(1):ITMP(1)),RN,NUM(1),KFLG)         THS05810
                IND(1)=ITMP(1)                                          THS05820
                RETURN                                                  THS05830
        ENDIF                                                           THS05840
        I=0                                                             THS05850
3000    I=I+1                                                           THS05860
        IF(ITMP(I).NE.ITMP(I+1)-1)THEN                                  THS05870
                N=N+1                                                   THS05880
                S2=S1(ITMP(I):ITMP(I))                                  THS05890
                CALL READNU(S2,RN,NUM(N),KFLG)                          THS05900
C*              CALL READNU(S1(ITMP(I):ITMP(I)),RN,NUM(N),KFLG)         THS05910
                IND(N)=ITMP(I)                                          THS05920
        ELSE IF(ITMP(I).EQ.ITMP(I+1)-1)THEN                             THS05930
                N=N+1                                                   THS05940
                S2=S1(ITMP(I):ITMP(I+1))                                THS05950
                CALL READNU(S2,RN,NUM(N),KFLG)                          THS05960
C*              CALL READNU(S1(ITMP(I):ITMP(I+1)),RN,NUM(N),KFLG)       THS05970
                IND(N)=ITMP(I)                                          THS05980
                I=I+1                                                   THS05990
        ENDIF                                                           THS06000
        IF(I.LT.NOC-1)GOTO 3000                                         THS06010
        IF(ITMP(NOC-1)+1.NE.ITMP(NOC))THEN                              THS06020
                N=N+1                                                   THS06030
                S2=S1(ITMP(NOC):ITMP(NOC))                              THS06040
                CALL READNU(S2,RN,NUM(N),KFLG)                          THS06050
C*              CALL READNU(S1(ITMP(NOC):ITMP(NOC)),RN,NUM(N),KFLG)     THS06060
                                                                        THS06070
                IND(N)=ITMP(NOC)                                        THS06080
        ENDIF                                                           THS06090
        RETURN                                                          THS06100
        END                                                             THS06110
                                                                        THS06120
        SUBROUTINE FINTER(FORM,LEN1,IE,NE,NATOM,IENTRY)                 THS06130
        CHARACTER*2 ELEM(50)                                            THS06140
        CHARACTER*70 FORM,DUMMY                                         THS06150
        CHARACTER*4 IE(4),B,C,IE1(70),TMP                               THS06160
        INTEGER NE(4),M,M1,N(4),L,NE1(70),ITEMP,LEN1,ID(70),NOC1(70),   THS06170
     +  NOC2(70),IN2C(35,4),IN1C(70,4),I2CH(35),I1CH(50),N2CH,N1CH,     THS06180
     +  IFM(4),IG(35),NG,INUM(70),NUMBER(70),NNUM                       THS06190
        COMMON/CONFG1/ELEM                                              THS06200
        COMMON/CONFG2/NELEM,N1CH,N2CH,I1CH,I2CH                         THS06210
C *   INITALIZE ARRAYS                                                  THS06220
        DO 1234 I=1,4                                                   THS06230
        IE(I)=' '                                                       THS06240
1234    NE(I)=0                                                         THS06250
C       IF(LEN1.EQ.1)GOTO 1000                                          THS06260
C*********************************************************************  THS06270
C               SEARCH FOR NUMERIC CHARACTERS                        *  THS06280
C*********************************************************************  THS06290
        CALL INDXNU(FORM,LEN1,NNUM,NUMBER,INUM,KFLG)                    THS06300
C*********************************************************************  THS06310
C               BEGIN SEARCH FOR 2 CHARACTER ELEMENTS                *  THS06320
C*********************************************************************  THS06330
        I1=0                                                            THS06340
        DO 1000 I=1,N2CH                                                THS06350
        DUMMY=ELEM(I2CH(I))(:2)                                         THS06360
        CALL INDXST(FORM,DUMMY,NOC,ID,LEN1,2,KFLG)                      THS06370
        IF(NOC.NE.0)THEN                                                THS06380
                I1=I1+1                                                 THS06390
                NOC2(I1)=NOC                                            THS06400
                                                                        THS06410
        IF(I1.GT.1)THEN                                                 THS06420
                NNC=0                                                   THS06430
                DO 9100 K=1,NOC                                         THS06440
                DO 9200 IJ=1, I1-1                                      THS06450
                DO 9300 IJK = 1, NOC2(I1-1)                             THS06460
                IF(IN2C(IJK,IJ).EQ.(ID(K)-1))THEN                       THS06470
                        ID(K) = -1                                      THS06480
                        NNC = NNC+1                                     THS06490
                ENDIF                                                   THS06500
9300            CONTINUE                                                THS06510
9200            CONTINUE                                                THS06520
9100            CONTINUE                                                THS06530
                IF(NOC-NNC.GT.0)THEN                                    THS06540
                        NOC2(I1)=NOC-NNC                                THS06550
                ELSE                                                    THS06560
                        I1 = I1 - 1                                     THS06570
                        GOTO 1000                                       THS06580
                ENDIF                                                   THS06590
        ENDIF                                                           THS06600
                IFM(I1)=I2CH(I)                                         THS06610
                J1=0                                                    THS06620
                DO 100 J=1,NOC                                          THS06630
                IF(ID(J).EQ.-1)GOTO 100                                 THS06640
                J1=J1+1                                                 THS06650
                DO 150 K1=1,NNUM                                        THS06660
                IF(ID(J)+2.EQ.INUM(K1))THEN                             THS06670
                        NOC2(I1)=NOC2(I1)+NUMBER(K1)-1                  THS06680
                ENDIF                                                   THS06690
150             CONTINUE                                                THS06700
                IN2C(J1,I1)=ID(J)                                       THS06710
100             CONTINUE                                                THS06720
        ENDIF                                                           THS06730
        NEL2=I1                                                         THS06740
        NELTOT=I1                                                       THS06750
        IF(NELTOT.GT.4)THEN                                             THS06760
C                CALL CLRSCR                                            THS06770
                WRITE(*,*)'  ERROR : FORMULA APPEARS TO CONTAIN'        THS06780
                WRITE(*,*)'  MORE THAN 4 ELEMENTS and CANNOT BE'        THS06790
                WRITE(*,*)' PUT INTO NASA FORMAT FOR USE WITH CHEMKIN.' THS06800
                READ(*,3)B                                              THS06810
3               FORMAT(A4)                                              THS06820
                IENTRY=-1                                               THS06830
                RETURN                                                  THS06840
        ENDIF                                                           THS06850
1000    CONTINUE                                                        THS06860
C                                                                       THS06870
C*********************************************************************  THS06880
C               BEGIN SEARCH FOR 1 CHARACTER ELEMENTS                *  THS06890
C*********************************************************************  THS06900
        I1=0                                                            THS06910
        N1C=0                                                           THS06920
        DO 2000 I=1,N1CH                                                THS06930
        DUMMY=ELEM(I1CH(I))(:1)                                         THS06940
        N1CS=N1C                                                        THS06950
        CALL INDXST(FORM,DUMMY,NOC,ID,LEN1,1,KFLG)                      THS06960
        IF(NOC.NE.0)THEN                                                THS06970
                ISAME=0                                                 THS06980
                I1=I1+1                                                 THS06990
        DO 30 J=1,NEL2                                                  THS07000
        IF(ELEM(I1CH(I))(:1).EQ.ELEM(IFM(J))(:1)                        THS07010
     +.OR.ELEM(I1CH(I))(:1).EQ.ELEM(IFM(J))(2:2))THEN                   THS07020
                ISAME=1                                                 THS07030
                NO1CA=0                                                 THS07040
                NG=0                                                    THS07050
        DO 41 I3=1,NOC                                                  THS07060
        DO 40 I2=1,NOC2(J)                                              THS07070
                IF(ID(I3).EQ.IN2C(I2,J).OR.ID(I3).EQ.(IN2C(I2,J)+1))THENTHS07080
                        NG=NG+1                                         THS07090
                        IG(NG)=I3                                       THS07100
                ENDIF                                                   THS07110
40      CONTINUE                                                        THS07120
41      CONTINUE                                                        THS07130
        DO 45 I2=1,NOC                                                  THS07140
                                                                        THS07150
                                                                        THS07160
        IST=0                                                           THS07170
                DO 46 I3=1,NG                                           THS07180
46              IF(I2.EQ.IG(I3))IST=1                                   THS07190
                IF(IST.EQ.0)THEN                                        THS07200
                        NO1CA=NO1CA+1                                   THS07210
                DO 350 K1=1,NNUM                                        THS07220
                IF(ID(I2)+1.EQ.INUM(K1))THEN                            THS07230
                        NO1CA=NO1CA+NUMBER(K1)-1                        THS07240
                ENDIF                                                   THS07250
350             CONTINUE                                                THS07260
                        IN1C(I2,I1)=ID(I2)                              THS07270
                ENDIF                                                   THS07280
45      CONTINUE                                                        THS07290
        ENDIF                                                           THS07300
                                                                        THS07310
30      CONTINUE                                                        THS07320
        IF(ISAME.EQ.0)THEN                                              THS07330
                NOC1(I1)=NOC                                            THS07340
                IFM(NEL2+I1)=I1CH(I)                                    THS07350
                DO 50 I2=1,NOC                                          THS07360
                DO 250 K1=1,NNUM                                        THS07370
                IF(ID(I2)+1.EQ.INUM(K1))THEN                            THS07380
                        NOC1(I1)=NOC1(I1)+NUMBER(K1)-1                  THS07390
                ENDIF                                                   THS07400
250             CONTINUE                                                THS07410
50              IN1C(I2,I1)=ID(I2)                                      THS07420
                NIC=I1                                                  THS07430
        ELSE IF(ISAME.NE.0)THEN                                         THS07440
                IF(NO1CA.EQ.0)THEN                                      THS07450
                        I1=I1-1                                         THS07460
                ELSE IF(NO1CA.NE.0)THEN                                 THS07470
                        IFM(NEL2+I1)=I1CH(I)                            THS07480
                        NOC1(I1)=NO1CA                                  THS07490
                ENDIF                                                   THS07500
                NIC=I1                                                  THS07510
        ENDIF                                                           THS07520
        IF(NEL2+I1.GT.4)THEN                                            THS07530
C                CALL CLRSCR                                            THS07540
                WRITE(*,*)'  ERROR : FORMULA APPEARS TO CONTAIN'        THS07550
                WRITE(*,*)'  MORE THAN 4 ELEMENTS and CANNOT BE'        THS07560
                WRITE(*,*)' PUT INTO NASA FORMAT FOR USE WITH CHEMKIN.' THS07570
                READ(*,3)B                                              THS07580
                IENTRY=-1                                               THS07590
                RETURN                                                  THS07600
        ENDIF                                                           THS07610
        ENDIF                                                           THS07620
2000    CONTINUE                                                        THS07630
        NEL1=I1                                                         THS07640
        NATOM=0                                                         THS07650
        DO 324 J2=1,NEL1                                                THS07660
324             NATOM=NATOM+NOC1(J2)                                    THS07670
        DO 325 J2=1,NEL2                                                THS07680
325             NATOM=NATOM+NOC2(J2)                                    THS07690
        DO 60 I=1,NEL2                                                  THS07700
        IE(I)(:2)=ELEM(IFM(I))                                          THS07710
                                                                        THS07720
                                                                        THS07730
60      NE(I)=NOC2(I)                                                   THS07740
        DO 70 I=1,NEL1                                                  THS07750
        IE(I+NEL2)(:2)=ELEM(IFM(I+NEL2))                                THS07760
70      NE(I+NEL2)=NOC1(I)                                              THS07770
                                                                        THS07780
71      CONTINUE                                                        THS07790
        RETURN                                                          THS07800
        END                                                             THS07810
                                                                        THS07820
        SUBROUTINE SLEN (S1,LSTR)                                       THS07830
        CHARACTER*(*) S1                                                THS07840
        INTEGER LSTR, I                                                 THS07850
        LSTR=LEN(S1)                                                    THS07860
        IF(S1.EQ.' ')THEN                                               THS07870
                LSTR=0                                                  THS07880
                GOTO 2000                                               THS07890
        ENDIF                                                           THS07900
        I=0                                                             THS07910
1000    I=I+1                                                           THS07920
        IF(S1(I+1:LSTR).EQ.' ')THEN                                     THS07930
                LSTR=I                                                  THS07940
        ELSEIF(I.LE.LSTR)THEN                                           THS07950
                GOTO 1000                                               THS07960
        ENDIF                                                           THS07970
2000    CONTINUE                                                        THS07980
        RETURN                                                          THS07990
        END                                                             THS08000
                                                                        THS08010
        SUBROUTINE SUBSTI(S1,S2,S3,LEN1,LEN2)                           THS08020
        CHARACTER*(*) S1,S2,S3                                          THS08030
        INTEGER LEN1,LEN2                                               THS08040
        DO 100 I=1,(LEN1+1-LEN2)                                        THS08050
        IF(S1(I:I+LEN2-1).EQ.S2(:LEN2))THEN                             THS08060
                S1(I:I+LEN2-1)=S3(:LEN2)                                THS08070
        ENDIF                                                           THS08080
100     CONTINUE                                                        THS08090
        RETURN                                                          THS08100
        END                                                             THS08110
                                                                        THS08120
        SUBROUTINE TBREAK(TMID,TINIT)                                   THS08130
C       FUNCTION TBREAK(T,F)                                            THS08140
        IMPLICIT REAL*8(A-H,O-Z)                                        THS08150
        REAL*8 TMID                                                     THS08160
        DIMENSION T(30),F(14)                                           THS08170
        COMMON/TEMFIT/T,F                                               THS08180
        COMMON/TOL/TOL,TSTEP                                            THS08190
        EXTERNAL FX                                                     THS08200
                                                                        THS08210
        CALL TMIDOP(TMID,F)                                             THS08220
        DEL=ABS(CPT(0,TMID,F)-CPT(7,TMID,F))                            THS08230
        DEL2=ABS(DCPT(0,TMID,F)-DCPT(7,TMID,F))                         THS08240
        IF(DEL.LE.TOL.AND.DEL2.LE.(2.*TOL))THEN                         THS08250
                GOTO 200                                                THS08260
        ENDIF                                                           THS08270
C       IF(DEL.GT.(10.*TOL).OR.DEL2.GT.(20.*TOL))THEN                   THS08280
C               write(*,*)' Tmid convergence failed ! '                 THS08290
C               write(*,*)'Tmid = ',tmid                                THS08300
C               write(*,*)'error cp = ',del                             THS08310
C               write(*,*)'error dcp = ',del2                           THS08320
C       ENDIF                                                           THS08330
200     CONTINUE                                                        THS08340
        RETURN                                                          THS08350
        END                                                             THS08360
                                                                        THS08370
                                                                        THS08380
        SUBROUTINE TMIDOP(T,F)                                          THS08390
C***********************************************************************THS08400
C  THIS SUBROUTINE SOLVES FOR THE BEST MID-POINT TEMPERATURE BETWEEN    THS08410
C  TWO OVERLAPPED POLYNOMIALS.  TMID TRIES TO MATCH BOTH Cp AND dCp/dT  THS08420
C  AT THE BREAK POINT TEMPERATURE                                       THS08430
C                                                                       THS08440
C      WRITTEN BY:   EDWARD RITTER      1/24/89                         THS08450
C***********************************************************************THS08460
        IMPLICIT REAL*8 (A-H,O-Z)                                       THS08470
        DIMENSION X0(2),WX(2),S(2)                                      THS08480
        EXTERNAL FX                                                     THS08490
        WX(1)=0.0                                                       THS08500
        WX(2)=0.0                                                       THS08510
        X0(2)=0.0                                                       THS08520
        S(2)=0.0                                                        THS08530
        A0=0.                                                           THS08540
        A1=0.                                                           THS08550
        A2=0.                                                           THS08560
        F0=0.                                                           THS08570
        F1=0.                                                           THS08580
        F2=0.                                                           THS08590
        T1=T                                                            THS08600
        AS=0.                                                           THS08610
        FS=0.                                                           THS08620
C       OPEN(9,FILE='NUL',STATUS='UNKNOWN')                             THS08630
        S(1)=-1.0                                                       THS08640
        X0(1)=1500.0                                                    THS08650
        EPS=1.0E-6                                                      THS08660
        N=1                                                             THS08670
        CALL BOUND(X0,WX,FX,N,A0,A1,A2,F0,F1,F2,S)                      THS08680
        CALL GOLD(X0,WX,FX,N,A0,A1,A2,F0,F1,F2,EPS,AS,FS,S)             THS08690
        T=X0(1)+AS                                                      THS08700
C       ITMID=INT((T*100.)+0.5)                                         THS08710
C       T=FLOAT(ITMID)/100.0                                            THS08720
        T=FLOAT(INT(T+0.5))                                             THS08730
        RETURN                                                          THS08740
        END                                                             THS08750
                                                                        THS08760
        FUNCTION FX(X)                                                  THS08770
C       FUNCTION FOR DETERMINATION OF BREAK POINT TEMPERATURE           THS08780
C       RETURNS SSR DEVIATION Cph:Cpl                                   THS08790
        IMPLICIT REAL*8 (A-H,O-Z)                                       THS08800
        DIMENSION X(1),T(30),F(14)                                      THS08810
        COMMON/TEMFIT/T,F                                               THS08820
        TMID=X(1)                                                       THS08830
        IF(TMID.LT.500.0.OR.TMID.GT.3000.0)THEN                         THS08840
                FX=1.0E+6                                               THS08850
                RETURN                                                  THS08860
        ENDIF                                                           THS08870
        DEL=ABS(CPT(0,TMID,F)-CPT(7,TMID,F))                            THS08880
        DEL2=ABS(DCPT(0,TMID,F)-DCPT(7,TMID,F))                         THS08890
        SSR=0.0                                                         THS08900
        SSR=(DEL**2)+(DEL2**2)                                          THS08910
100     CONTINUE                                                        THS08920
        FX=SSR                                                          THS08930
        RETURN                                                          THS08940
        END                                                             THS08950
                                                                        THS08960
C     ******************************************************************THS08970
C     *                                                                *THS08980
C     *                                                                *THS08990
C     *      SUBROUTINE NAME: BOUND                                    *THS09000
C     *                                                                *THS09010
C     *      PURPOSE:  THIS SUBROUTINE LOCATES 3 STEP SIZES A0,A1,A2   *THS09020
C     *                THAT BOUND THE MINIMUM OF A UNIMODAL OBJECTIVE  *THS09030
C     *                FUNCTION F(X) ALONG A DIRECTION DEFINED BY S(N).*THS09040
C     *                                                                *THS09050
C     *      CALLED BY SUBROUTINES: $MAIN,GBASE,PCD                    *THS09060
C     *                                                                *THS09070
C     *      CALLING SEQUENCE:                                         *THS09080
C     *                                                                *THS09090
C     *          CALL BOUND (XO,WX,F,N,A0,A1,A2,F0,F1,F2,S)            *THS09100
C     *                                                                *THS09110
C     *      ARGUMENTS:                                                *THS09120
C     *                                                                *THS09130
C     *          INPUT:XO=STARTING POINT                               *THS09140
C     *                F=NAME OF SUBPROGRAM CONTAINING F(X)            *THS09150
C     *                N=NUMBER OF INDEPENDENT VARIABLES               *THS09160
C     *                S=SEARCH DIRECTION VECTOR                       *THS09170
C     *                                                                *THS09180
C     *          OUTPUT:WX=WORKING ARRAY                               *THS09190
C     *                 A0,A1,A2=STEP SIZES                            *THS09200
C     *                 F0,F1,F2=LAST VALUES OF F(X)                   *THS09210
C                                                                       THS09220
C   Added a catch to prevent infinite loops caused by bad Cp data.
C     <pey 2-jul-04>
C                                                                       THS09230
C     ******************************************************************THS09240
      SUBROUTINE BOUND(XO,WX,F,N,A0,A1,A2,F0,F1,F2,S)                   THS09250
      IMPLICIT REAL*8 (A-H,O-Z)                                         THS09260
      DIMENSION XO(N),S(N),WX(N)
      external F

    iter = 1
C     WRITE(9,*)'CALL TO SUBROUTINE BOUND'                              THS09280
      SN=0.0                                                            THS09290
      DO 90 I=1,N                                                       THS09300
   90 SN=SN+S(I)*S(I)                                                   THS09310
      STEP=.25/SQRT(SN)                                                 THS09320
      A0=0.0                                                            THS09330
      A1=STEP                                                           THS09340
C     WRITE(9,*)'FIRST CALL TO FOFX FROM THIS CALL TO BOUND'            THS09350
      F0=F(XO)                                                          THS09360
      DO 100 I=1,N                                                      THS09370
      WX(I)=XO(I)-A1*S(I)                                               THS09380
  100 CONTINUE                                                          THS09390
C     WRITE(9,*)'ANOTHER CALL TO FOFX FROM BOUND:CALC F1'               THS09400
      F1=F(WX)                                                          THS09410
      IF(F0.LT.F1) GO TO 104                                            THS09420
  101 STEP=STEP*2.                                                      THS09430
      A2=STEP+A1                                                        THS09440
      DO 102 I=1,N                                                      THS09450
      WX(I)=XO(I)-A2*S(I)                                               THS09460
  102 CONTINUE                                                          THS09470
C     WRITE(9,*)'ANOTHER CALL TO FOFX FROM BOUND:CALC F2'               THS09480
      F2=F(WX)                                                          THS09490
      IF(F2.GT.F1) RETURN                                               THS09500
      A0=A1                                                             THS09510
      A1=A2                                                             THS09520
      F0=F1                                                             THS09530
      F1=F2                                                             THS09540
c   begin pey
    iter=iter+1
    if(iter.gt.500) stop 'Therfit Error: Failure in subroutine BOUND'
c   end pey
      GO TO 101                                                         THS09550
  104 A2=A1                                                             THS09560
      A1=A0                                                             THS09570
      F2=F1                                                             THS09580
      F1=F0                                                             THS09590
      A0=-STEP+A1                                                       THS09600
      DO 105 I=1,N                                                      THS09610
      WX(I)=XO(I)-A0*S(I)                                               THS09620
  105 CONTINUE                                                          THS09630
C     WRITE(9,*)'ANOTHER CALL TO FOFX FROM BOUND:DET NEW F0(105)'       THS09640
      F0=F(WX)                                                          THS09650
C     IF(F0.GT.F1) WRITE(9,*)'CALL TO BOUND COMPLETED'                  THS09660
      IF(F0.GT.F1) RETURN                                               THS09670
      STEP=2.*STEP                                                      THS09680
C     WRITE(9,*)'STEP SIZE IS DOUBLED FOR NEW ATTEMPT AT F0'            THS09690
c   begin pey
    iter=iter+1
    if(iter.gt.500) stop 'Therfit Error: Failure in subroutine BOUND'
c   end pey
      GO TO 104                                                         THS09700
C                                                                       THS09710
      END                                                               THS09720
                                                                        THS09730
                                                                        THS09740
C     ******************************************************************THS09750
C     *                                                                *THS09760
C     *                                                                *THS09770
C     *      SUBROUTINE NAME: GOLD                                     *THS09780
C     *                                                                *THS09790
C     *      PURPOSE: THIS SUBROUTINE LOCATES THE MINIMUM OF A UNIMODAL*THS09800
C     *               FUNCTION F(X), BY THE METHOD OF GOLDEN SECTION   *THS09810
C     *               THE MINIMUM IS BOUNDED BY SUBROUTINE BOUND.      *THS09820
C     *                                                                *THS09830
C     *      CALLED BY SUBROUTINE: GBASE,PCD                           *THS09840
C     *                                                                *THS09850
C     *      CALLING SEQUENCE:                                         *THS09860
C     *                                                                *THS09870
C     *          CALL GOLD (XO,WX,F,N,A0,A1,A2,F0,F1,F2,EPS,AS,FS,S)   *THS09880
C     *                                                                *THS09890
C     *      ARGUMENTS:                                                *THS09900
C     *                                                                *THS09910
C     *          INPUT: XO   - STARTING POINT                          *THS09920
C     *                 F    - NAME OF FUNCTION SUBPROGRAM DESCRIBING  *THS09930
C     *                        F(X)                                    *THS09940
C     *                 N    - NUMBER OF INDEPENDENT VARIABLES         *THS09950
C     *                 A0,A1,A2 - BOUNDING STEP LENGTHS               *THS09960
C     *                 F0,F1,F2 - CORESPONDING VALUES OF F(X)         *THS09970
C     *                 EPS  - CONVERGENCE CRITERIA                    *THS09980
C     *                 S    - SEARCH DIRECTION VECTOR                 *THS09990
C     *                                                                *THS10000
C     *          OUTPUT: AS - OPTIMAL STEP LENGTH                      *THS10010
C     *                  FS - OPTIMAL VALUE OF F(X)                    *THS10020
C     *                  A0,A1,A2 - BOUNDING STEP LENGTHS              *THS10030
C     *                  F0,F1,F2 - CORSPONDING VALUES OF F(X)         *THS10040
C     *                                                                *THS10050
C     *                                                                *THS10060
C     *                                                                *THS10070
C     ******************************************************************THS10080
      SUBROUTINE GOLD(XO,WX,F,N,A0,A1,A2,F0,F1,F2,EPS,AS,FS,S)          THS10090
      IMPLICIT REAL*8 (A-H,O-Z)                                         THS10100
      DIMENSION XO(N),WX(N),S(N)
      external F
C     WRITE(9,*)'CALL TO SUBROUTINE GOLD'                               THS10120
C      FIBONACC1 NUMBER                                                 THS10130
      FIB=(3.-SQRT(5.))/2.                                              THS10140
      H=ABS(A2-A0)                                                      THS10150
C     GOLDEN SECTION POINT A3 AND A4                                    THS10160
      A3=A0+H*FIB                                                       THS10170
      A4=A2-H*FIB                                                       THS10180
      DO 100 I=1,N                                                      THS10190
      WX(I)=XO(I)-A3*S(I)                                               THS10200
  100 CONTINUE                                                          THS10210
C     WRITE(9,*)'FIRST CALL TO FOFX FROM GOLD :F3'                      THS10220
      F3=F(WX)                                                          THS10230
      DO 101 I=1,N                                                      THS10240
      WX(I)=XO(I)-A4*S(I)                                               THS10250
  101 CONTINUE                                                          THS10260
C     WRITE(9,*)'SECOND CALL TO FOFX FROM GOLD :F4'                     THS10270
      F4=F(WX)                                                          THS10280
C     CONVERGENCE CHECK                                                 THS10290
  102 IF(H.LE.EPS) GO TO 107                                            THS10300
C     MAINTAIN BOUND ON MINIMUM                                         THS10310
      IF(F3.LT.F4) GO TO 104                                            THS10320
C     F4 .LT. F3                                                        THS10330
      A0=A3                                                             THS10340
      F0=F3                                                             THS10350
      A3=A4                                                             THS10360
      F3=F4                                                             THS10370
      A4=A2-(A2-A0)*FIB                                                 THS10380
      DO 103 I=1,N                                                      THS10390
      WX(I)=XO(I)-A4*S(I)                                               THS10400
  103 CONTINUE                                                          THS10410
C     WRITE(9,*)'ANOTHER CALL TO FOFX FOR F4'                           THS10420
      F4=F(WX)                                                          THS10430
      GO TO 106                                                         THS10440
C     F3 .LT. F4                                                        THS10450
  104 A2=A4                                                             THS10460
      F2=F4                                                             THS10470
      A4=A3                                                             THS10480
      F4=F3                                                             THS10490
      A3=A0+(A2-A0)*FIB                                                 THS10500
      DO 105 I=1,N                                                      THS10510
      WX(I)=XO(I)-A3*S(I)                                               THS10520
  105 CONTINUE                                                          THS10530
C     WRITE(9,*)'ANOTHER CALL TO FOFX FOR F3'                           THS10540
      F3=F(WX)                                                          THS10550
  106 H=A2-A0                                                           THS10560
      GO TO 102                                                         THS10570
  107 IF(F1.LT.F3.AND.F1.LT.F4) GO TO 109                               THS10580
      IF(F3.LT.F4) GO TO 108                                            THS10590
      AS=A4                                                             THS10600
      FS=F4                                                             THS10610
                                                                        THS10620
                                                                        THS10630
C     WRITE(9,*) 'CALL TO GOLD COMPLETED'                               THS10640
      RETURN                                                            THS10650
  108 AS=A3                                                             THS10660
      FS=F3                                                             THS10670
C     WRITE(9,*)'CALL TO GOLD COMPLETED'                                THS10680
      RETURN                                                            THS10690
  109 AS=A1                                                             THS10700
      FS=F1                                                             THS10710
C     WRITE(9,*) 'CALL TO GOLD COMPLETED'                               THS10720
      RETURN                                                            THS10730
      END                                                               THS10740
                                                                        THS10750
        FUNCTION HF(I,T,F)                                              THS10760
C                                                                       THS10770
        IMPLICIT REAL*8 (A-H,O-Z)                                       THS10780
        DIMENSION F(*)                                                  THS10790
        HF = F(1+I)+(F(6+I)/T)+((F(2+I)/2.)*T)+((F(3+I)/3.)*(T**2))     THS10800
     + +((F(4+I)/4.)*(T**3))+((F(5+I)/5.)*(T**4))                       THS10810
        RETURN                                                          THS10820
        END                                                             THS10830
                                                                        THS10840
        FUNCTION DHF(I,T,F)                                             THS10850
C                                                                       THS10860
        IMPLICIT REAL*8 (A-H,O-Z)                                       THS10870
        DIMENSION F(*)                                                  THS10880
        DHF=(F(2+I)/2.)+(-1.*F(6+I)/(T**2))+((2./3.)*F(3+I)*T)+         THS10890
     + ((3./4.)*F(4+I)*(T**2))+((4./5.)*F(5+I)*(T**3))                  THS10900
        RETURN                                                          THS10910
        END                                                             THS10920
                                                                        THS10930
                                                                        THS10940
        FUNCTION CPT(I,T,F)                                             THS10950
C                                                                       THS10960
        IMPLICIT REAL*8 (A-H,O-Z)                                       THS10970
        DIMENSION F(*)                                                  THS10980
        CPT=F(1+I)+(F(2+I)*T)+(F(3+I)*(T**2))+(F(4+I)*(T**3))+          THS10990
     +  (F(5+I)*(T**4))                                                 THS11000
        RETURN                                                          THS11010
        END                                                             THS11020
                                                                        THS11030
        FUNCTION DCPT(I,T,F)                                            THS11040
C                                                                       THS11050
        IMPLICIT REAL*8 (A-H,O-Z)                                       THS11060
        DIMENSION F(*)                                                  THS11070
        DCPT=F(2+I)+(2.*F(3+I)*T)+(3.*F(4+I)*(T**2))+(4.*F(5+I)*(T**3)) THS11080
        RETURN                                                          THS11090
        END                                                             THS11100
                                                                        THS11110
        FUNCTION ST(I,T,F)                                              THS11120
C                                                                       THS11130
        IMPLICIT REAL*8 (A-H,O-Z)                                       THS11140
        DIMENSION F(*)                                                  THS11150
        ST=F(7+I)+(F(1+I)*LOG(T))+(F(2+I)*T)+((F(3+I)/2.)*(T**2))+      THS11160
     +  ((F(4+I)/3.)*(T**3))+((F(5+I)/4.)*(T**4))                       THS11170
        RETURN                                                          THS11180
        END                                                             THS11190
                                                                        THS11200
        FUNCTION DST(I,T,F)                                             THS11210
C                                                                       THS11220
      IMPLICIT REAL*8 (A-H,O-Z)                                         THS11230
      DIMENSION F(*)                                                    THS11240
      DST=(F(1+I)/T)+F(2+I)+(F(3+I)*T)+(F(4+I)*(T**2))+(F(5+I)*(T**3))  THS11250
      RETURN                                                            THS11260
      END                                                               THS11270
                                                                        THS11280
        FUNCTION F6(I,HOFT,T,F)                                         THS11290
C                                                                       THS11300
        IMPLICIT REAL*8(A-H,O-Z)                                        THS11310
        DIMENSION F(*)                                                  THS11320
      X=HOFT-((T*F(1+I))+((F(2+I)/2.)*(T**2))+((F(3+I)/3.)*             THS11330
     #(T**3))+((F(4+I)/4.)*(T**4))+((F(5+I)/5.)*(T**5)))                THS11340
        F6=X                                                            THS11350
        RETURN                                                          THS11360
        END                                                             THS11370
                                                                        THS11380
                                                                        THS11390
        FUNCTION F7(I,SOFT,T,F)                                         THS11400
C                                                                       THS11410
        IMPLICIT REAL*8(A-H,O-Z)                                        THS11420
        REAL*8 T,F7                                                     THS11430
        DIMENSION F(*)                                                  THS11440
      X=SOFT-((F(1+I)*DLOG(T))+(F(2+I)*T)+((F(3)/2.)*(T**2))            THS11450
     #+((F(4+I)/3.)*(T**3))+((F(5+I)/4.)*(T**4)))                       THS11460
        F7=X                                                            THS11470
        RETURN                                                          THS11480
        END                                                             THS11490
                                                                        THS11500
                                                                        THS11510
        SUBROUTINE CVVIB(FREQ,T,CV,THETA)                               THS11520
C                                                                       THS11530
        IMPLICIT REAL*8(A-H,O-Z)                                        THS11540
        REAL*8 HP,KB,C,X,THETA,TEST                                     THS11550
        HP=6.62517E-34                                                  THS11560
        KB=1.3804E-23                                                   THS11570
        C=2.998E+10                                                     THS11580
        IF(FREQ.GT.25000.0.OR.FREQ.LT.1.)THEN                           THS11590
                CV=1000.                                                THS11600
                RETURN                                                  THS11610
        ENDIF                                                           THS11620
        IF(FREQ.NE.0.0)THEN                                             THS11630
                THETA=((HP*C)/KB)*FREQ                                  THS11640
        ELSE                                                            THS11650
                IF(THETA.EQ.0.)THEN                                     THS11660
                        CV=1000.                                        THS11670
                        RETURN                                          THS11680
                                                                        THS11690
                                                                        THS11700
                ENDIF                                                   THS11710
        ENDIF                                                           THS11720
        X=THETA/T                                                       THS11730
C       WRITE(*,*)'X= ',X                                               THS11740
C       WRITE(*,*)'THETA= ',THETA                                       THS11750
C       WRITE(*,*)'FREQ= ',FREQ                                         THS11760
C       WRITE(*,*)'T= ',T                                               THS11770
C       IF(X.LT.10-6)THEN                                               THS11780
C               CV=100.                                                 THS11790
C               RETURN                                                  THS11800
C       ENDIF                                                           THS11810
        IF(X.GT.230.)THEN                                               THS11820
                CV=1000.                                                THS11830
                RETURN                                                  THS11840
        ENDIF                                                           THS11850
        TEST=DEXP(X)                                                    THS11860
        IF(TEST-1.0.LT.10.0E-6)THEN                                     THS11870
                CV=100.                                                 THS11880
                RETURN                                                  THS11890
        ENDIF                                                           THS11900
        TOP=EXP(X)                                                      THS11910
        BOTTOM=(TOP-1.0)**2                                             THS11920
        CV=(X**2)*(TOP/BOTTOM)                                          THS11930
        RETURN                                                          THS11940
        END                                                             THS11950
                                                                        THS11960
                                                                        THS11970
                                                                        THS11980
        SUBROUTINE SIGMCP(B)                                            THS11990
        IMPLICIT REAL*8(A-H,O-Z)                                        THS12000
        CHARACTER*70 FILPRN                                             THS12010
        DIMENSION B(*),CONST(3)                                         THS12020
        COMMON/CONST/CONST,NTERMS                                       THS12030
        WRITE(*,*)'ENTER PRN FILENAME'                                  THS12040
        READ(*,70)FILPRN                                                THS12050
70      FORMAT(A70)                                                     THS12060
C       OPEN(12,FILE=FILPRN,STATUS='UNKNOWN')                           THS12070
        T=100.                                                          THS12080
        DO 100 I=1,60                                                   THS12090
        CP=0.                                                           THS12100
        DO 200 J=1,NTERMS                                               THS12110
        THETA=0.0                                                       THS12120
        IF(NTERMS.EQ.1)THEN                                             THS12130
                CALL CVVIB(B(2),T,CV,THETA)                             THS12140
        ELSE                                                            THS12150
                CALL CVVIB(B(NTERMS+J-1),T,CV,THETA)                    THS12160
        ENDIF                                                           THS12170
        IF(J.NE.NTERMS.OR.NTERMS.EQ.1)THEN                              THS12180
                TERM=B(J)*CV                                            THS12190
        ELSE IF(J.EQ.NTERMS.AND.NTERMS.NE.1)THEN                        THS12200
                                                                        THS12210
                                                                        THS12220
                SUM=0.0                                                 THS12230
                DO 234 IK=1,NTERMS-1                                    THS12240
234             SUM=SUM+B(IK)                                           THS12250
                TERM=(CONST(2)-SUM)*CV                                  THS12260
        ENDIF                                                           THS12270
        CP=CP+TERM                                                      THS12280
200     CONTINUE                                                        THS12290
        CP=(CP+4.0)*1.987*CONST(1)                                      THS12300
        WRITE(12,*)T,CP                                                 THS12310
        T=T+100.                                                        THS12320
100     CONTINUE                                                        THS12330
C       CLOSE(12,STATUS='KEEP')                                         THS12340
        RETURN                                                          THS12350
        END                                                             THS12360
                                                                        THS12370
                                                                        THS12380
        FUNCTION CPSERI(B,T)                                            THS12390
        IMPLICIT REAL*8(A-H,O-Z)                                        THS12400
        DIMENSION B(6),CONST(3)                                         THS12410
        COMMON/CONST/CONST,NTERMS                                       THS12420
        CP=0.0                                                          THS12430
        CV=0.0                                                          THS12440
        DO 100 I=1,3                                                    THS12450
                CALL CVVIB(B(I+3),T,CV,THETA)                           THS12460
                CP=CP+(B(I)*CV)                                         THS12470
100     CONTINUE                                                        THS12480
        IF(CONST(3).EQ.0.)THEN                                          THS12490
C               NON-LINEAR                                              THS12500
C               CP=(CP+4.0)*1.987*CONST(1)                              THS12510
                CP=(CP+4.0)                                             THS12520
        ELSE                                                            THS12530
C               LINEAR                                                  THS12540
C               CP=(CP+(7./2.))*1.987*CONST(1)                          THS12550
                CP=(CP+(7./2.))                                         THS12560
        ENDIF                                                           THS12570
        CPSERI=CP                                                       THS12580
        RETURN                                                          THS12590
        END                                                             THS12600
                                                                        THS12610
                                                                        THS12620
        SUBROUTINE READNU(STRING,RNUM,NUM,IERR)                         THS12630
C*--------------------------------------------------------              THS12640
C*    this subroutine is an interface to the revised                    THS12650
C*    string operation functions written 12/20/89                       THS12660
C*    by Ed Ritter.  This subroutine replaces the                       THS12670
C*    previous version of READNU                                        THS12680
C*--------------------------------------------------------              THS12690
        REAL*4 RNUM, REALN                                              THS12700
        CHARACTER*(*) STRING                                            THS12710
        INTEGER NUM, INTNUM, IERR                                       THS12720
        CALL DBLNK(STRING,LST,0)                                        THS12730
        RNUM = REALN(STRING,IERR)                                       THS12740
        IF(IERR.NE.0)GOTO 10                                            THS12750
        NUM = INTNUM (STRING,IERR)                                      THS12760
10      CONTINUE                                                        THS12770
                                                                        THS12780
        RETURN                                                          THS12790
        END                                                             THS12800
C*======================================================================THS12810
        SUBROUTINE DRDNUM(STRING,DNUM,IERR)                             THS12820
C*--------------------------------------------------------              THS12830
C*    this subroutine is an interface to the revised                    THS12840
C*    string operation functions written 12/20/89                       THS12850
C*    by Ed Ritter.  This subroutine replaces the                       THS12860
C*    previous version of DRDNUM                                        THS12870
C*--------------------------------------------------------              THS12880
        REAL*8 DNUM, DREALN                                             THS12890
        CHARACTER*(*) STRING                                            THS12900
        INTEGER IERR                                                    THS12910
        DNUM = DREALN(STRING,IERR)                                      THS12920
        RETURN                                                          THS12930
        END                                                             THS12940
C*======================================================================THS12950
        FUNCTION STLEN(STRING)                                          THS12960
C*------------------------------------------------------                THS12970
C*                                                                      THS12980
C*    FUNCTION STLEN                                                    THS12990
C*                                                                      THS13000
C* written by: Ed Ritter                                                THS13010
C*                                                                      THS13020
C*    this function returns the length of a character string            THS13030
C*    ignoring trailing blanks                                          THS13040
C*                                                                      THS13050
C*  NOTE:  the calling routine MUST declare STLEN an integer            THS13060
C*                                                                      THS13070
C*------------------------------------------------------                THS13080
        CHARACTER*(*) STRING                                            THS13090
        INTEGER STLEN,LST                                               THS13100
        CALL SLEN(STRING,LST)                                           THS13110
        STLEN = LST                                                     THS13120
        RETURN                                                          THS13130
        END                                                             THS13140
C*======================================================================THS13150
        FUNCTION DREALN(STRING,IERR)                                    THS13160
C*-----------------------------------------------------                 THS13170
C*                                                                      THS13180
C*  FUNCTION DREALN                                                     THS13190
C*                                                                      THS13200
C*   written by: Ed Ritter                                              THS13210
C*                                                                      THS13220
C*                                                                      THS13230
C*   this function converts a numeric character string                  THS13240
C*   into a double precision number.  The error flag                    THS13250
C*   IERR is set to 0 upon entering this function                       THS13260
C*   and has the following return values:                               THS13270
C*                                                                      THS13280
C*        ierr =  0 : normal return                                     THS13290
C*             = -1 : more than 1 decimal point detected                THS13300
C*             = -2 : string contains non-numeric characters            THS13310
C*             = -3 : function INUM error (this should not occur)       THS13320
C*             = -4 : integer expression out of range                   THS13330
                                                                        THS13340
C*                    (should not occur)                                THS13350
C*             = -5 : double precision number out of range              THS13360
C*             = -6 : single precision number out of range              THS13370
C*                    (should not occur)                                THS13380
C*------------------------------------------------------                THS13390
        REAL*8 DREALN, PROD, SUM, TERM, DCHECK                          THS13400
        CHARACTER*(*) STRING                                            THS13410
        CHARACTER EXPONE                                                THS13420
        INTEGER STLEN, LST, IERR, NCH, J1, ID, J, I                     THS13430
        LST = STLEN(STRING)                                             THS13440
        IERR = 0                                                        THS13450
        ID = 0                                                          THS13460
        IEXP = 0                                                        THS13470
        NEXP = 0                                                        THS13480
C*------------------------------------------------------------          THS13490
C*         verify string is numeric / return error code if not          THS13500
C*------------------------------------------------------------          THS13510
        DO 100 I = 1,LST                                                THS13520
        NCH = ICHAR(STRING(I:I))                                        THS13530
        IF(STRING(I:I).EQ.'E'.OR.STRING(I:I).EQ.'D'                     THS13540
     +.OR.STRING(I:I).EQ.'D')THEN                                       THS13550
                STRING(I:I) = 'E'                                       THS13560
        ENDIF                                                           THS13570
        IF(STRING(I:I).EQ.'.'.AND.ID.EQ.0)THEN                          THS13580
                ID = I                                                  THS13590
        ELSEIF(STRING(I:I).EQ.'-'.AND.I.EQ.1)THEN                       THS13600
C*              no action                                               THS13610
        ELSEIF(STRING(I:I).EQ.'E'.AND.IEXP.EQ.0)THEN                    THS13620
                IEXP = I                                                THS13630
                LEXP = LST - (IEXP + 1)                                 THS13640
                LSTS = LST                                              THS13650
        ELSEIF(IEXP.NE.0.AND.I.EQ.IEXP+1)THEN                           THS13660
                IF(STRING(I:I).EQ.'+')THEN                              THS13670
                        EXPONE = '+'                                    THS13680
                        LEXP = LEXP - 1                                 THS13690
                ELSEIF(STRING(I:I).EQ.'-')THEN                          THS13700
                        EXPONE = '-'                                    THS13710
                        LEXP = LEXP - 1                                 THS13720
                ELSEIF((NCH.GE.48.AND.NCH.LE.57).AND.(STRING(I:I).NE.'-'THS13730
     +.AND.STRING(I:I).NE.'+'))THEN                                     THS13740
                        EXPONE = '+'                                    THS13750
                ELSE                                                    THS13760
                        DREALN = 1.                                     THS13770
                        IERR = -2                                       THS13780
                        RETURN                                          THS13790
                ENDIF                                                   THS13800
        ELSEIF(STRING(I:I).EQ.'.'.AND.ID.NE.0)THEN                      THS13810
                DREALN = 1.                                             THS13820
                IERR = -1                                               THS13830
                RETURN                                                  THS13840
        ELSEIF(NCH.LT.48.OR.NCH.GT.57)THEN                              THS13850
                DREALN = 1.                                             THS13860
                IERR = -2                                               THS13870
                RETURN                                                  THS13880
        ENDIF                                                           THS13890
100     CONTINUE                                                        THS13900
        IF(IEXP.NE.0)THEN                                               THS13910
                LST = IEXP - 1                                          THS13920
        ENDIF                                                           THS13930
C*----------------------------------------------------------------      THS13940
C*    if the string is an E format number, determine the exponent       THS13950
C*----------------------------------------------------------------      THS13960
        IF(LEXP.GT.3)THEN                                               THS13970
                DREALN = 1.                                             THS13980
                IERR = -5                                               THS13990
                RETURN                                                  THS14000
        ENDIF                                                           THS14010
        J1 = LSTS + 1                                                   THS14020
        J2 = 0                                                          THS14030
        DO 10 K = 1, LEXP+1                                             THS14040
                J1 = J1 - 1                                             THS14050
                J2 = J2 + 1                                             THS14060
                NUM = INUM(STRING(J1:J1))                               THS14070
                NTERM = (10**(J2-1))*NUM                                THS14080
                NEXP = NEXP + NTERM                                     THS14090
10      CONTINUE                                                        THS14100
        IF(NEXP.GT.308)THEN                                             THS14110
                DREALN = 1.                                             THS14120
                IERR = -5                                               THS14130
                RETURN                                                  THS14140
        ENDIF                                                           THS14150
        IF(EXPONE.EQ.'-')THEN                                           THS14160
                NEXP = -1 * NEXP                                        THS14170
        ENDIF                                                           THS14180
C*----------------------------------------------------------------      THS14190
C*      determine the value of the number                               THS14200
C*----------------------------------------------------------------      THS14210
        J1 = 0                                                          THS14220
        SUM = 0.0                                                       THS14230
        IF(ID.EQ.0)THEN                                                 THS14240
                LD = LST                                                THS14250
        ELSE                                                            THS14260
                LD = ID - 1                                             THS14270
        ENDIF                                                           THS14280
C*                                                                      THS14290
      DO 2000 J = 1,LST                                                 THS14300
C*                                                                      THS14310
        IF(J.EQ.ID)THEN                                                 THS14320
                GOTO 2000                                               THS14330
        ELSE                                                            THS14340
                J1 = J1 + 1                                             THS14350
                I = LD - J1                                             THS14360
                IF(STRING(J:J).EQ.'-')GOTO 2000                         THS14370
        ENDIF                                                           THS14380
        NUM = INUM(STRING(J:J))                                         THS14390
        IF(NUM.LT.0)THEN                                                THS14400
                DREALN = 1.                                             THS14410
                IERR = -3                                               THS14420
                RETURN                                                  THS14430
        ENDIF                                                           THS14440
        IA=ABS(I)                                                       THS14450
        IF(STRING(1:1).EQ.'-')THEN                                      THS14460
                PROD=-1.                                                THS14470
        ELSE                                                            THS14480
                PROD=1.                                                 THS14490
        ENDIF                                                           THS14500
        IF(I.EQ.0)THEN                                                  THS14510
                GOTO 101                                                THS14520
        ENDIF                                                           THS14530
        DO 200 II=1,IA                                                  THS14540
        IF(I.GT.0)THEN                                                  THS14550
                PROD=PROD*10.                                           THS14560
        ELSEIF(I.LT.0)THEN                                              THS14570
                PROD=PROD/10.                                           THS14580
        ENDIF                                                           THS14590
200     CONTINUE                                                        THS14600
101     CONTINUE                                                        THS14610
        TERM = PROD * FLOAT(NUM)                                        THS14620
        SUM = SUM + TERM                                                THS14630
C*                                                                      THS14640
2000    CONTINUE                                                        THS14650
C*                                                                      THS14660
        IF(IEXP.NE.0)THEN                                               THS14670
                DCHECK = (LOG(ABS(SUM))/2.302585) + ABS(DBLE(NEXP))     THS14680
                IF(DCHECK.GT.308.2)THEN                                 THS14690
                        DREALN = 1.                                     THS14700
                        IERR = -5                                       THS14710
                        RETURN                                          THS14720
                ENDIF                                                   THS14730
                DREALN = SUM * (10.**NEXP)                              THS14740
        ELSE                                                            THS14750
                DREALN = SUM                                            THS14760
        ENDIF                                                           THS14770
        RETURN                                                          THS14780
        END                                                             THS14790
C*======================================================================THS14800
        FUNCTION REALN(STRING,IERR)                                     THS14810
C*-----------------------------------------------------                 THS14820
C*                                                                      THS14830
C*   this function converts a numeric character string                  THS14840
C*   into a single precision number.  The error flag                    THS14850
C*   IERR is set to 0 upon entering this function                       THS14860
C*   and has the following return values:                               THS14870
C*                                                                      THS14880
C*        ierr =  0 : normal return                                     THS14890
C*             = -1 : more than 1 decimal point detected                THS14900
C*             = -2 : string contains non-numeric characters            THS14910
C*             = -3 : function INUM error (this should not occur)       THS14920
C*             = -4 : integer expression out of range                   THS14930
C*                    (should not occur)                                THS14940
C*             = -5 : double precision number out of range              THS14950
C*             = -6 : single precision number out of range              THS14960
C*                                                                      THS14970
C*------------------------------------------------------                THS14980
        REAL*4 REALN                                                    THS14990
        REAL*8 DREALN, DRNUM                                            THS15000
        CHARACTER*(*) STRING                                            THS15010
        IERR = 0                                                        THS15020
        DRNUM = DREALN(STRING,IERR)                                     THS15030
        IF(DRNUM.EQ.0.)THEN                                             THS15040
                REALN=0.                                                THS15050
                RETURN                                                  THS15060
        ELSEIF(ABS(DRNUM).GT.3.4E+38.OR.ABS(DRNUM).LT.1.2E-38)THEN      THS15070
                REALN = 1.                                              THS15080
                IERR = -6                                               THS15090
                RETURN                                                  THS15100
        ELSE                                                            THS15110
                REALN = SNGL(DRNUM)                                     THS15120
        ENDIF                                                           THS15130
        RETURN                                                          THS15140
        END                                                             THS15150
C*======================================================================THS15160
        FUNCTION INTNUM(STRING,IERR)                                    THS15170
C*-----------------------------------------------------                 THS15180
C*                                                                      THS15190
C*   this function converts a numeric character string                  THS15200
C*   into an integer*4 number.  The error flag                          THS15210
C*   IERR is set to 0 upon entering this function                       THS15220
C*   and has the following return values:                               THS15230
C*                                                                      THS15240
C*        ierr =  0 : normal return                                     THS15250
C*             = -1 : more than 1 decimal point detected                THS15260
C*             = -2 : string contains non-numeric characters            THS15270
C*             = -3 : function INUM error (this should not occur)       THS15280
C*             = -4 : integer expression out of range                   THS15290
C*                    num > 2,147,483,640                               THS15300
C*             = -5 : double precision number out of range              THS15310
C*             = -6 : single precision number out of range              THS15320
C*                                                                      THS15330
C*------------------------------------------------------                THS15340
        REAL*8 RNUM, DREALN                                             THS15350
        CHARACTER*(*) STRING                                            THS15360
        INTEGER INTNUM, IERR                                            THS15370
        IERR = 0                                                        THS15380
        RNUM = DREALN(STRING,IERR)                                      THS15390
        IF(ABS(RNUM).GT.2147483640)THEN                                 THS15400
                INTNUM = 1                                              THS15410
                IERR = -4                                               THS15420
                RETURN                                                  THS15430
        ENDIF                                                           THS15440
        IF(RNUM.GE.0.)THEN                                              THS15450
                INTNUM = INT(RNUM + 0.5)                                THS15460
        ELSE                                                            THS15470
                INTNUM = INT(RNUM - 0.5)                                THS15480
        ENDIF                                                           THS15490
        RETURN                                                          THS15500
        END                                                             THS15510
        FUNCTION INUM(CH)                                               THS15520
        CHARACTER CH                                                    THS15530
                                                                        THS15540
        IF(ICHAR(CH).EQ.64) INUM=0                                      THS15550
        IF(ICHAR(CH).EQ.64) GOTO 100                                    THS15560
        IF(ICHAR(CH).GE.240.AND.ICHAR(CH).LE.249)THEN                   THS15570
                INUM=ICHAR(CH)-240                                      THS15580
        ELSE                                                            THS15590
C*              ERROR                                                   THS15600
                INUM=-1                                                 THS15610
        ENDIF                                                           THS15620
100     CONTINUE                                                        THS15630
        RETURN                                                          THS15640
        END                                                             THS15650
                                                                        THS15660
                                                                        THS15670
        SUBROUTINE CSTRIN(S,N)                                          THS15680
        CHARACTER*(*) S(*)                                              THS15690
100     CONTINUE                                                        THS15700
        DO 200 I=1,N                                                    THS15710
        IF(S(I).EQ.' ')THEN                                             THS15720
                DO 300 I1=I+1,N                                         THS15730
300             S(I1-1)=S(I1)                                           THS15740
                N=N-1                                                   THS15750
                GOTO 100                                                THS15760
        ENDIF                                                           THS15770
200     CONTINUE                                                        THS15780
        RETURN                                                          THS15790
        END                                                             THS15800
                                                                        THS15810
                                                                        THS15820
        SUBROUTINE DCREAL(X,N)                                          THS15830
        IMPLICIT REAL*8(A-H,O-Z)                                        THS15840
        DIMENSION X(*)                                                  THS15850
100     CONTINUE                                                        THS15860
        DO 200 I=1,N                                                    THS15870
        IF(X(I).EQ.0.0)THEN                                             THS15880
                DO 300 I1=I+1,N                                         THS15890
300             X(I1-1)=X(I1)                                           THS15900
                N=N-1                                                   THS15910
                GOTO 100                                                THS15920
        ENDIF                                                           THS15930
200     CONTINUE                                                        THS15940
        RETURN                                                          THS15950
        END                                                             THS15960
                                                                        THS15970
        SUBROUTINE ORDERL(NE,IE)                                        THS15980
        INTEGER NE(4),NEL(4), NTMP                                      THS15990
        CHARACTER*4 IE(4), PASS, ITMP                                   THS16000
        DO 100 I=1,4                                                    THS16010
        IF(NE(I).EQ.0)THEN                                              THS16020
                NEL(I)=999                                              THS16030
        ELSE                                                            THS16040
                NEL(I)=NELEMN(IE(I))                                    THS16050
        ENDIF                                                           THS16060
100     CONTINUE                                                        THS16070
        DO 200 I=1,3                                                    THS16080
                DO 300 J=1,4-I                                          THS16090
                IF(NEL(J).GT.NEL(J+1))THEN                              THS16100
                        NTMP = NEL(J)                                   THS16110
                        NEL(J) = NEL(J+1)                               THS16120
                        NEL(J+1) = NTMP                                 THS16130
                        NTMP = NE(J)                                    THS16140
                        ITMP = IE(J)                                    THS16150
                        NE(J) = NE(J+1)                                 THS16160
                        IE(J) = IE(J+1)                                 THS16170
                        NE(J+1) = NTMP                                  THS16180
                        IE(J+1) = ITMP                                  THS16190
                ENDIF                                                   THS16200
300             CONTINUE                                                THS16210
200     CONTINUE                                                        THS16220
        RETURN                                                          THS16230
        END                                                             THS16240
        FUNCTION NELEMN(IE)                                             THS16250
        INTEGER NELEM,N1CH,N2CH,I1CH(50),I2CH(35)                       THS16260
        CHARACTER*2 ELEM(50)                                            THS16270
        CHARACTER*4 IE, PASS                                            THS16280
        COMMON/CONFG1/ELEM                                              THS16290
        COMMON/CONFG2/NELEM,N1CH,N2CH,I1CH,I2CH                         THS16300
        DO 100 I = 1, NELEM                                             THS16310
                IF(IE(:2).EQ.ELEM(I)(:2))THEN                           THS16320
                        NELEMN = I                                      THS16330
                        RETURN                                          THS16340
                ENDIF                                                   THS16350
100     CONTINUE                                                        THS16360
        NELEMN = 999                                                    THS16370
        RETURN                                                          THS16380
        END                                                             THS16390
                                                                        THS16400
        SUBROUTINE ATOMSO(ISORT,NE1,IE1,ICNT,KEY)                       THS16410
        INTEGER ISORT(*),NE1(4,400),ICNT,LEVEL,KEY(10),IL1(2),          THS16420
     +IL2(2),INC,ICTRL,IK,J,I,K,ITMP,NELEM,N1CH,N2CH,I1CH(50),          THS16430
     +I2CH(35)                                                          THS16440
        CHARACTER*2 ELEM(50)                                            THS16450
                                                                        THS16460
        CHARACTER*4 IE1(4,400), PASS                                    THS16470
        COMMON/CONFG1/ELEM                                              THS16480
        COMMON/CONFG2/NELEM,N1CH,N2CH,I1CH,I2CH                         THS16490
        COMMON/CONFG5/NKEYS                                             THS16500
C       INITIALIZE INDIRECT ADDERSS VECTOR: ISORT                       THS16510
        DO 10 I=1,ICNT                                                  THS16520
10      ISORT(I)=I                                                      THS16530
C       INITIALIZE DEFAULT SORT KEY VECTOR. KEY IS DEFINED AS THE FIRST THS16540
C       6 ELEMENTS LISTED IN FILE : "THERM.CFG"                         THS16550
                DO 20 I=1,NKEYS                                         THS16560
20              KEY(I)=I                                                THS16570
        DO 10000 LEVEL=1,NKEYS                                          THS16580
C*      level = 1                                                       THS16590
        ICTRL=0                                                         THS16600
        JLAST = 1                                                       THS16610
        DO 1000 I=1,ICNT-1                                              THS16620
                IF(ICTRL.EQ.1.AND.LEVEL.EQ.NKEYS)THEN                   THS16630
                        RETURN                                          THS16640
                ELSEIF(ICTRL.EQ.1.AND.LEVEL.LT.NKEYS)THEN               THS16650
                        GOTO 1000                                       THS16660
                ENDIF                                                   THS16670
                ICTRL=1                                                 THS16680
                IF((FLOAT(I)/2.).EQ.FLOAT(I/2))THEN                     THS16690
                        PASS = 'EVEN'                                   THS16700
                ELSE                                                    THS16710
                        PASS = 'ODD'                                    THS16720
                ENDIF                                                   THS16730
                IF(PASS.EQ.'ODD')THEN                                   THS16740
                        INC=1                                           THS16750
                        J=JLAST - INC                                   THS16760
                ELSE                                                    THS16770
                        INC = -1                                        THS16780
                        J = JLAST - INC                                 THS16790
                ENDIF                                                   THS16800
                DO 2000 IK=1,ICNT-I                                     THS16810
                        J=J+INC                                         THS16820
                        JLAST = J                                       THS16830
                        IF((J+INC.GT.ICNT.OR.J+INC.LT.0)                THS16840
     +                  .OR.(J.GT.ICNT.OR.J.LT.1))GOTO 2000             THS16850
C*                      FUNCTION NUMELM(J,K,IE1,NE1)                    THS16860
                        IL1(1) = NUMELM(LEVEL,ISORT(J),IE1,NE1)         THS16870
                        IL2(1) = NUMELM(LEVEL,ISORT(J+INC),IE1,NE1)     THS16880
                        IF(LEVEL.GT.1)THEN                              THS16890
                        DO 1010 IJK=1,LEVEL-1                           THS16900
                           IL1(2) = NUMELM(IJK,ISORT(J),IE1,NE1)        THS16910
                           IL2(2) = NUMELM(IJK,ISORT(J+INC),IE1,NE1)    THS16920
                           IF(IL1(2).NE.IL2(2))GOTO 2000                THS16930
1010                    CONTINUE                                        THS16940
                        ENDIF                                           THS16950
                        IF(PASS.EQ.'ODD')THEN                           THS16960
                                IF(IL2(1).GE.IL1(1))GOTO 2000           THS16970
                                GOTO 150                                THS16980
                        ELSE                                            THS16990
                                IF(IL2(1).LE.IL1(1))GOTO 2000           THS17000
                                GOTO 150                                THS17010
                        ENDIF                                           THS17020
150                     CONTINUE                                        THS17030
                        ICTRL=0                                         THS17040
                        ITMP=ISORT(J)                                   THS17050
                        ISORT(J)=ISORT(J+INC)                           THS17060
                        ISORT(J+INC)=ITMP                               THS17070
2000            CONTINUE                                                THS17080
1000    CONTINUE                                                        THS17090
10000   CONTINUE                                                        THS17100
3       FORMAT(A14)                                                     THS17110
C        CALL CLRSCR                                                    THS17120
        RETURN                                                          THS17130
        END                                                             THS17140
                                                                        THS17150
        FUNCTION NUMELM(J,K,IE1,NE1)                                    THS17160
        INTEGER NE1(4,400)                                              THS17170
        CHARACTER*2 ELEM(50)                                            THS17180
        CHARACTER*4 IE1(4,400)                                          THS17190
        COMMON/CONFG1/ELEM                                              THS17200
        DO 10 I=1, 4                                                    THS17210
        IF(IE1(I,K)(:2).EQ.ELEM(J)(:2))THEN                             THS17220
                NUMELM = NE1(I,K)                                       THS17230
                RETURN                                                  THS17240
        ENDIF                                                           THS17250
10      CONTINUE                                                        THS17260
        NUMELM = 0                                                      THS17270
        RETURN                                                          THS17280
        END                                                             THS17290
                                                                        THS17300
        SUBROUTINE INQNUM(S1,IANS,NUMLEN,LEN,IBL)                       THS17310
C*********************************************************************  THS17320
C*                                                                  *   THS17330
C*              written by Ed Ritter : 11/20/88                     *   THS17340
C*                                                                  *   THS17350
C*********************************************************************  THS17360
C*                                                                  *   THS17370
C*      THIS SUBROUTINE SCANS STRING S1 FOR A LEADING NUMBER        *   THS17380
C*       The search is carried out from left to right.               *  THS17390
C*       The first non-numeric character encountered ends            *  THS17400
C*        the numeric substring.                                        THS17410
C*                                                                  *   THS17420
C*      return codes :IANS                                          *   THS17430
C*                                                                  *   THS17440
C*      IANS = 0  : the string is non-numeric                       *   THS17450
C*            = 1  : string S1 is an integer                        *   THS17460
C*            = 2  : S1 is a decimal                                *   THS17470
C*            = 3  : S1 has a leading integer of length NUMLEN       *  THS17480
C*            = 4  : S1 has a leading decimal of length NUMLEN.      *  THS17490
C*                                                                  *   THS17500
C*       NUMLEN    : the length of the leading numeric portion of S1 *  THS17510
C*       LEN       : the length of the S1 without blanks             *  THS17520
C*       IBL       : allows the calling routine to control blank     *  THS17530
                                                                        THS17540
                                                                        THS17550
C*                   removal in subroutine DEBLNK                    *  THS17560
C*                   IBL=0 : leave no blanks in S1                   *  THS17570
C*                  IBL=1 : leave one remaining blank in each blank *   THS17580
C*                          field                                   *   THS17590
C*                                                                  *   THS17600
C*********************************************************************  THS17610
        CHARACTER*(*) S1                                                THS17620
        INTEGER IANS,NUMLEN,LEN,ICTRL,IDEC,IBL                          THS17630
C*                                                                      THS17640
C*  ICTRL <> 0, when a non-numeric character has been encountered in S1.THS17650
C*  IDEC <> 0, if a decimal point is encountered in a numeric field;    THS17660
C*             in this case IDEC = address of decimal point.            THS17670
                                                                        THS17680
10      CONTINUE                                                        THS17690
                                                                        THS17700
        NUMLEN=0                                                        THS17710
        ICTRL=0                                                         THS17720
        IDEC=0                                                          THS17730
        IANS=0                                                          THS17740
                                                                        THS17750
                                                                        THS17760
C*  inserted for test                                                   THS17770
C*      write(*,*)'enter string s1'                                     THS17780
C*      read(*,'(a70)')s1                                               THS17790
C*      ibl=0                                                           THS17800
C*  end of insert#1                                                     THS17810
                                                                        THS17820
                                                                        THS17830
        CALL DBLNK (S1,LEN,IBL)                                         THS17840
        DO 100 I=1,LEN                                                  THS17850
        IF((ICHAR(S1(I:I)).GE.48.AND.ICHAR(S1(I:I)).LE.57).OR.          THS17860
     +(S1(I:I).EQ.'.'.AND.IDEC.EQ.0))THEN                               THS17870
                IF(ICTRL.EQ.0)THEN                                      THS17880
                        NUMLEN=NUMLEN+1                                 THS17890
                        IF(S1(I:I).EQ.'.')IDEC=I                        THS17900
                ENDIF                                                   THS17910
        ELSE                                                            THS17920
                ICTRL=1                                                 THS17930
        ENDIF                                                           THS17940
100     CONTINUE                                                        THS17950
        IF(NUMLEN.EQ.0)THEN                                             THS17960
                IANS=0                                                  THS17970
        ELSEIF(NUMLEN.EQ.LEN)THEN                                       THS17980
                IF(IDEC.EQ.0)THEN                                       THS17990
                        IANS=1                                          THS18000
                ELSEIF(IDEC.NE.0)THEN                                   THS18010
                        IANS=2                                          THS18020
                ENDIF                                                   THS18030
        ELSEIF(NUMLEN.LT.LEN)THEN                                       THS18040
                IF(IDEC.EQ.0)THEN                                       THS18050
                        IANS=3                                          THS18060
                ELSEIF(IDEC.NE.0)THEN                                   THS18070
                        IANS=4                                          THS18080
                ENDIF                                                   THS18090
        ENDIF                                                           THS18100
                                                                        THS18110
                                                                        THS18120
                                                                        THS18130
C*      insert for test                                                 THS18140
C*      write(*,*)'S1 = ',s1                                            THS18150
C*      write(*,*)'IANS = ',ians                                        THS18160
C*      write(*,*)'NUMLEN = ',numlen                                    THS18170
C*      write(*,*)'IDEC = ',idec                                        THS18180
C*      write(*,*)'LEN = ',len                                          THS18190
C*      pause ' '                                                       THS18200
C*      goto 10                                                         THS18210
C*      end of insert#2                                                 THS18220
                                                                        THS18230
99999   CONTINUE                                                        THS18240
        RETURN                                                          THS18250
C*      STOP ' '                                                        THS18260
        END                                                             THS18270
                                                                        THS18280
        SUBROUTINE CONFIG(NFILES,KS,ISTAT,OUTPUT)                       THS18290
        CHARACTER*70 NUL,DUMMY,KS(1),HLP(10),ERRORM                     THS18300
        CHARACTER*14 NUMIN,OUTPUT,UNITS,TRANGE                          THS18310
        CHARACTER*2 ELEM(50),DRIVE1,DRIVE2,DRIVE3,ERRORC                THS18320
        INTEGER I,LEN1,NELEM,NFILES,ISTAT,I1CH(50),I2CH(35),
     +     N1CH,N2CH
        COMMON/CONFG1/ELEM                                              THS18340
        COMMON/CONFG2/NELEM,N1CH,N2CH,I1CH,I2CH                         THS18350
        COMMON/CONFG3/HLP                                               THS18360
        COMMON/CONFG4/DRIVE1,DRIVE2,DRIVE3                              THS18370
        COMMON/CONFG5/NKEYS                                             THS18380
        COMMON/CONFG6/UNITS                                             THS18390
        COMMON/CONFG7/ERRORC                                            THS18400
        COMMON/CONFG8/TRANGE                                            THS18410
        TRANGE='2000'                                                   THS18420
        UNITS='KCAL'                                                    THS18430
        OUTPUT='QUICK'                                                  THS18440
        DRIVE1=' '                                                      THS18450
        DRIVE2=' '                                                      THS18460
        DRIVE3=' '                                                      THS18470
        NKEYS=3                                                         THS18480
        ERRORC=' '                                                      THS18490
2       FORMAT(A2)                                                      THS18500
14      FORMAT(A14)                                                     THS18510
70      FORMAT(A70)                                                     THS18520
        ISTAT=0                                                         THS18530
C       OPEN(96,FILE='THERM ',STATUS='OLD',ERR=9999)                    THS18540
        ERRORM=' KEYWOR " #FILES " MISSING OR NOT FIRST LINE'           THS18550
        READ(96,70,END=9999)NUL                                         THS18560
        CALL LRDCHR(NUL)                                                THS18570
        CALL DEBLNK(NUL,DUMMY,LEN1,0)                                   THS18580
        IF(NUL(:6).NE.'#FILES')GOTO 9999                                THS18590
        NUMIN=NUL(7:LEN1)                                               THS18600
        CALL READNU(NUMIN,RNUM,NFILES,KFLG)                             THS18610
        ERRORM=' FEWER FILES FOUND THAN SPECIFIED IN " #FILES " '       THS18620
        DO 10 I=1,NFILES                                                THS18630
        READ(96,70,END=9999)KS(I)                                       THS18640
        CALL DBLNK(KS(I),LEN1,0)                                        THS18650
                                                                        THS18660
                                                                        THS18670
        IF(KS(I).EQ.'#ELEMENTS')GOTO 9999                               THS18680
10      CONTINUE                                                        THS18690
        ERRORM=' KEYWOR " #ELEMENTS " IS MISSING'                       THS18700
        READ(96,70,END=9999)NUL                                         THS18710
        CALL LRDCHR(NUL)                                                THS18720
        CALL DEBLNK(NUL,DUMMY,LEN1,0)                                   THS18730
        IF(NUL(:9).NE.'#ELEMENTS')GOTO 9999                             THS18740
        NUMIN=NUL(10:LEN1)                                              THS18750
        CALL READNU(NUMIN,RNUM,NELEM,KFLG)                              THS18760
        N1CH=0                                                          THS18770
        N2CH=0                                                          THS18780
        ERRORM=' FEWER ELEMENTS FOUND THAN SPECIIED IN " #ELEMENTS"'    THS18790
        DO 20 I=1,NELEM                                                 THS18800
        READ(96,70,END=9999)NUL                                         THS18810
        CALL DEBLNK(NUL,DUMMY,LEN1,0)                                   THS18820
        ELEM(I)=NUL(:LEN1)                                              THS18830
C       WRITE(*,*)NUL                                                   THS18840
C       WRITE(*,*)'LENGTH=',LEN1                                        THS18850
C       WRITE(*,*)'ELEM(I)=',ELEM(I)                                    THS18860
        IF(LEN1.EQ.1)THEN                                               THS18870
                N1CH=N1CH+1                                             THS18880
                I1CH(N1CH)=I                                            THS18890
                GOTO 20                                                 THS18900
        ENDIF                                                           THS18910
        IF(LEN1.EQ.2)THEN                                               THS18920
                N2CH=N2CH+1                                             THS18930
                I2CH(N2CH)=I                                            THS18940
        ENDIF                                                           THS18950
        IF(LEN1.GT.2)THEN                                               THS18960
C*              CALL CLS                                                THS18970
C*              WRITE(*,*)'ERROR IN CONFIG FILE'                        THS18980
C*              WRITE(*,*)NUL                                           THS18990
C*              WRITE(*,*)'CONTAINS MORE THAN TWO CHARACTERS'           THS19000
C*              WRITE(*,*)'HIT RETURN TO CONTINUE'                      THS19010
C*              READ(*,70)DUMMY                                         THS19020
                GOTO 9999                                               THS19030
        ENDIF                                                           THS19040
20      CONTINUE                                                        THS19050
        ERRORM=' KEYWOR " #HELP " MISSING'                              THS19060
        READ(96,70,END=9999)NUL                                         THS19070
        CALL LRDCHR(NUL)                                                THS19080
        CALL DEBLNK(NUL,DUMMY,LEN1,0)                                   THS19090
        IF(NUL(:5).NE.'#HELP')GOTO 9999                                 THS19100
        NUMIN=NUL(6:LEN1)                                               THS19110
        CALL READNU(NUMIN,RNUM,NHELP,KFLG)                              THS19120
        ERRORM=' FEWER HELP FILES FOUND THAN SPECIFIED IN " #HELP "'    THS19130
        DO 60 I=1,NHELP                                                 THS19140
        READ(96,70,END=9999)HLP(I)                                      THS19150
        CALL DBLNK(HLP(I),LEN1,0)                                       THS19160
        IF(HLP(I)(:6).EQ.'OPTION')GOTO 9999                             THS19170
60      CONTINUE                                                        THS19180
                                                                        THS19190
C*                                                                      THS19200
C*    OPTIONS ANY ORDER IS ACCEPTABLE                                   THS19210
C*                                                                      THS19220
                                                                        THS19230
                                                                        THS19240
C*     DEFAULTS ARE DEFINED: ONLY NECESSARY IF WANT TO SPECIFY          THS19250
C*      OPTION DIFFERENT FROM DEFAULT                                   THS19260
C*                                                                      THS19270
10000   CONTINUE                                                        THS19280
        READ(96,70,END=9777)NUL                                         THS19290
        CALL DEBLNK(NUL,DUMMY,LEN1,0)                                   THS19300
        IF(NUL(:12).EQ.'THERMFIT.EXE')THEN                              THS19310
                IF(NUL(13:14).NE.' ')DRIVE1=NUL(13:14)                  THS19320
        ELSE IF(NUL(:12).EQ.'THERMLST.EXE')THEN                         THS19330
                IF(NUL(13:14).NE.' ')DRIVE2=NUL(13:14)                  THS19340
        ELSE IF(NUL(:12).EQ.'THERMRXN.EXE')THEN                         THS19350
                IF(NUL(13:14).NE.' ')DRIVE3=NUL(13:14)                  THS19360
        ELSEIF(NUL(:5).EQ.'#KEYS')THEN                                  THS19370
                NUMIN=NUL(6:LEN1)                                       THS19380
                CALL READNU(NUMIN,RNUM,NKEYS,KFLG)                      THS19390
        ELSEIF(NUL(:11).EQ.'ERRORC=OFF')THEN                            THS19400
                ERRORC='NO'                                             THS19410
        ELSEIF(NUL(:9).EQ.'TRANGE=EX'.OR.NUL(:8).EQ.'TRANGE=5')THEN     THS19420
                TRANGE='5000'                                           THS19430
        ELSEIF(NUL(:6).EQ.'UNITS:')THEN                                 THS19440
                IF(NUL(7:8).EQ.'KJ')UNITS='KJ'                          THS19450
        ELSEIF(NUL(:7).EQ.'OUTPUT:')THEN                                THS19460
                IF(NUL(8:LEN1).EQ.'CLOSE')OUTPUT='CLOSE'                THS19470
        ENDIF                                                           THS19480
        GOTO 10000                                                      THS19490
C9777   CLOSE(96,STATUS='KEEP')                                         THS19500
 9777   CONTINUE                                                        THS19510
C*      WRITE(*,*)'TRANGE = ',TRANGE(:4)                                THS19520
        ISTAT=0                                                         THS19530
        RETURN                                                          THS19540
9999    ISTAT=1                                                         THS19550
        CALL LINES(0,1)                                                 THS19560
        WRITE(*,'(1X,A70)')ERRORM                                       THS19570
        CALL LINES(0,1)                                                 THS19580
        WRITE(*,*)'   ERROR ENCOUNTERED IN THERM.CFG'                   THS19590
        WRITE(*,*)'   REQUIRED CONFIG DATA IS MISSING'                  THS19600
C       CLOSE(96,STATUS='KEEP')                                         THS19610
        RETURN                                                          THS19620
        END                                                             THS19630
                                                                        THS19640
        SUBROUTINE NATOMS(J,IE1,NE1,IEL,NEL,NATOM)                      THS19650
        INTEGER NE1(4,400),NELEM,N1CH,N2CH,I1CH(50),I2CH(35),IEL(*),    THS19660
     +NEL(*),NATOM                                                      THS19670
        CHARACTER*2 ELEM(50)                                            THS19680
        CHARACTER*4 IE1(4,400)                                          THS19690
        CHARACTER*14 NUL                                                THS19700
        COMMON/CONFG1/ELEM                                              THS19710
        COMMON/CONFG2/NELEM,N1CH,N2CH,I1CH,I2CH                         THS19720
        DO 10 I=1,4                                                     THS19730
10      IEL(I)=0                                                        THS19740
        DO 20 I=1,NELEM                                                 THS19750
20      NEL(I)=0                                                        THS19760
                                                                        THS19770
                                                                        THS19780
        NATOM=0                                                         THS19790
        I1=0                                                            THS19800
        DO 1000 I=1,4                                                   THS19810
        DO 1050 K=1,NELEM                                               THS19820
        IF(ELEM(K)(:2).EQ.IE1(I,J)(:2))THEN                             THS19830
                I1=I1+1                                                 THS19840
                IEL(I1)=K                                               THS19850
                NEL(K)=NEL(K)+NE1(I,J)                                  THS19860
        ENDIF                                                           THS19870
1050    CONTINUE                                                        THS19880
1000    CONTINUE                                                        THS19890
        RETURN                                                          THS19900
        END                                                             THS19910
                                                                        THS19920
        SUBROUTINE MATBAL(NLHS,LX,NRHS,RX,IE1,NE1,COEF,IERR)            THS19930
        INTEGER LX(*),RX(*),NE1(4,400),IELL(4),NELL(50),IELR(4),        THS19940
     +NELR(20),NELLT(20),NELRT(20),NELEM,I1CH(50),I2CH(35),LNATOM,      THS19950
     +RNATOM                                                            THS19960
        REAL*8 COEF(*),DNELL(50),DNELLT(50),DNELR(50),DNELRT(50)        THS19970
        CHARACTER*4 IE1(4,400)                                          THS19980
        CHARACTER*2 ELEM(50)                                            THS19990
        COMMON/CONFG1/ELEM                                              THS20000
        COMMON/CONFG2/NELEM,N1CH,N2CH,I1CH,I2CH                         THS20010
C*      CALL CONFIG(ELEM,NFILES,KS,ISTAT,OUTPUT)                        THS20020
C*      WRITE(*,*)'   CALL TO MATBAL'                                   THS20030
C*      WRITE(*,*)' LEFT HAND SIDE COMPOSITION'                         THS20040
C*      DO 11 I=1,NLHS                                                  THS20050
C*11    WRITE(*,'(1X,4(1X,A2,1X,I3,3X))')(IE1(P,LX(I)),NE1(P,LX(I)),P=1,THS20060
C*      WRITE(*,*)' RIGHT HAND SIDE COMPOSITION'                        THS20070
C*      DO 12 I=1,NRHS                                                  THS20080
C*12    WRITE(*,'(1X,4(1X,A2,1X,I3,3X))')(IE1(P,RX(I)),NE1(P,RX(I)),P=1,THS20090
        DO 10 I=1,NELEM                                                 THS20100
        DNELLT(I)=0.0                                                   THS20110
        DNELL(I)=0.0                                                    THS20120
        NELLT(I)=0                                                      THS20130
        DNELR(I)=0.0                                                    THS20140
        DNELRT(I)=0.0                                                   THS20150
10      NELRT(I)=0                                                      THS20160
        LNATOM=0                                                        THS20170
        RNATOM=0                                                        THS20180
        DO 800 I=1,NLHS                                                 THS20190
        CALL NATOMS(LX(I),IE1,NE1,IELL,NELL,NATOM)                      THS20200
        DO 850 I1=1,NELEM                                               THS20210
        DNELL(I1)=COEF(I)*FLOAT(NELL(I1))                               THS20220
C*850   NELLT(I1)=NELLT(I1)+NELL(I1)                                    THS20230
C*      WRITE(*,*)'ELEM = ',ELEM(I1)                                    THS20240
        DNELLT(I1)=DNELLT(I1)+DNELL(I1)                                 THS20250
C*      WRITE(*,*)'DNELLT(I)= ',DNELLT(I1)                              THS20260
850     CONTINUE                                                        THS20270
        LNATOM=NATOM                                                    THS20280
800     CONTINUE                                                        THS20290
                                                                        THS20300
                                                                        THS20310
        DO 900 I=1,NRHS                                                 THS20320
        CALL NATOMS(RX(I),IE1,NE1,IELL,NELL,NATOM)                      THS20330
        DO 950 I1=1,NELEM                                               THS20340
        DNELR(I1)=COEF(I+NLHS)*FLOAT(NELL(I1))                          THS20350
C*      NELL(I1)=REAL(COEF(I+NLHS))*NELL(I1)                            THS20360
C*950   NELRT(I1)=NELRT(I1)+NELL(I1)                                    THS20370
        DNELRT(I1)=DNELRT(I1)+DNELR(I1)                                 THS20380
C*      WRITE(*,*)'ELEM = ',ELEM(I1)                                    THS20390
C*      WRITE(*,*)'DNELRT(I)= ',DNELRT(I1)                              THS20400
950     CONTINUE                                                        THS20410
        RNATOM=NATOM                                                    THS20420
900     CONTINUE                                                        THS20430
        DO 1000 I=1,NELEM                                               THS20440
C*      WRITE(*,'(1X,A2,2X,F10.3,2X,F10.3)')ELEM(I),DNELRT(I),DNELLT(I) THS20450
        IF(INT(100.*DNELRT(I)).NE.INT(100.*DNELLT(I)))THEN              THS20460
                WRITE(*,*)' REACTION DOES NOT BALANCE!'                 THS20470
                WRITE(*,7000)ELEM(I),DNELLT(I)                          THS20480
                WRITE(*,7001)ELEM(I),DNELRT(I)                          THS20490
7000    FORMAT(2X,'# of ',A,' on left side = ',F5.0)                    THS20500
7001    FORMAT(2X,'# of ',A,' on right side = ',F5.0)                   THS20510
                IERR=1                                                  THS20520
                RETURN                                                  THS20530
        ENDIF                                                           THS20540
1000    CONTINUE                                                        THS20550
        IERR=0                                                          THS20560
        RETURN                                                          THS20570
        END                                                             THS20580
                                                                        THS20590
                                                                        THS20600
        SUBROUTINE IND127(S1,S2,N,INDX,LS1,LS2,FLG)                     THS20610
        CHARACTER*70 S2,NUL                                             THS20620
        CHARACTER*127 S1                                                THS20630
        INTEGER N,INDX(70),LS1,LS2,FLG                                  THS20640
3       FORMAT(A70)                                                     THS20650
        DO 10 I=1,70                                                    THS20660
10      INDX(I)=0                                                       THS20670
        N=0                                                             THS20680
        FLG=0                                                           THS20690
        IF(LS1.LT.LS2)RETURN                                            THS20700
        IF(LS1.NE.0.AND.LS2.NE.0)GOTO 150                               THS20710
        LS1=0                                                           THS20720
        LS2=0                                                           THS20730
        DO 100 I=1,127                                                  THS20740
        IF(LS1.NE.0.AND.LS2.NE.0)GOTO 100                               THS20750
        IP1=I+1                                                         THS20760
        IM1=I-1                                                         THS20770
        IF(IM1.LE.0)IM1=1                                               THS20780
        IF(LS1.EQ.0.AND.S1(I:I).EQ.' '.AND.S1(IP1:IP1).EQ.' ')LS1=IM1   THS20790
        IF(LS2.EQ.0.AND.S2(I:I).EQ.' '.AND.S2(IP1:IP1).EQ.' ')LS2=IM1   THS20800
100     CONTINUE                                                        THS20810
150     CONTINUE                                                        THS20820
                                                                        THS20830
                                                                        THS20840
C   NOW WE KNOW THE LENGTH OF BOTH STRINGS                              THS20850
C        IF(LS2.GT.LS1.AND.S2.NE.S1)THEN                                THS20860
C                WRITE(*,*)'   ERROR ENCOUNTERED IN SUBROUTINE INDXST'  THS20870
C                WRITE(*,*)' '                                          THS20880
C                WRITE(*,*)'         { hit return to continue }'        THS20890
C                READ(*,3)NUL                                           THS20900
C                FLG=0                                                  THS20910
C                RETURN                                                 THS20920
C        ENDIF                                                          THS20930
        IFOUND=0                                                        THS20940
        DO 200 I=1,LS1                                                  THS20950
        I2=I+LS2-1                                                      THS20960
        IF((I2-I).GT.LS2)GOTO 200                                       THS20970
        IF(S1(I:I2).EQ.S2(:LS2))THEN                                    THS20980
                IFOUND=IFOUND+1                                         THS20990
                IF(IFOUND.GT.70)THEN                                    THS21000
                        FLG=1                                           THS21010
                        RETURN                                          THS21020
                ENDIF                                                   THS21030
                INDX(IFOUND)=I                                          THS21040
        ENDIF                                                           THS21050
200     CONTINUE                                                        THS21060
        N=IFOUND                                                        THS21070
        RETURN                                                          THS21080
        END                                                             THS21090
                                                                        THS21100
          SUBROUTINE DEBLNK(S1,S2,LST,CCSW)                             THS21110
          CHARACTER*(*) S1,S2                                           THS21120
          INTEGER LS1,CCSW,FLG                                          THS21130
          FLG=0                                                         THS21140
          LST=LEN(S1)                                                   THS21150
          NNB=0                                                         THS21160
          KOUNT=0                                                       THS21170
C         WRITE(*,*)'ENTER STRING'                                      THS21180
C         READ(*,3)S1                                                   THS21190
3         FORMAT(A70)                                                   THS21200
          S2=S1                                                         THS21210
          DO 100 I=1,LST                                                THS21220
          IF(S1(I:I).NE.' ')NNB=NNB+1                                   THS21230
          IF(S1(I:LST).EQ.' ')THEN                                      THS21240
                LST=I-1                                                 THS21250
                GOTO 450                                                THS21260
          ENDIF                                                         THS21270
C         WRITE(*,*)'NNB= ',NNB                                         THS21280
100       CONTINUE                                                      THS21290
450       CONTINUE                                                      THS21300
400       I=1                                                           THS21310
          KOUNT=KOUNT+1                                                 THS21320
300       CONTINUE                                                      THS21330
          IF(I.GT.LST)GOTO 200                                          THS21340
          IF(S1(I:I).EQ.' '.AND.(CCSW.EQ.0.OR.I.EQ.1))THEN              THS21350
                                                                        THS21360
                                                                        THS21370
                S1(I:LST-1)=S1(I+1:LST)                                 THS21380
                S1(LST:LST)=' '                                         THS21390
                LST=LST-1                                               THS21400
          ENDIF                                                         THS21410
          IP1=I+1                                                       THS21420
          IF(S1(I:I).EQ.' '.AND.S1(IP1:IP1).EQ.' '.AND.CCSW.EQ.1)THEN   THS21430
                S1(I:LST-1)=S1(I+1:LST)                                 THS21440
                S1(LST:LST)=' '                                         THS21450
                LST=LST-1                                               THS21460
          ENDIF                                                         THS21470
          I=I+1                                                         THS21480
          GOTO 300                                                      THS21490
200       CONTINUE                                                      THS21500
          IF(KOUNT.GT.LST)GOTO 500                                      THS21510
          IF(LST.GT.NNB)GOTO 400                                        THS21520
500       CONTINUE                                                      THS21530
C         WRITE(*,*)S2                                                  THS21540
C         WRITE(*,*)S1                                                  THS21550
C         WRITE(*,*)'LST= ',LST                                         THS21560
          RETURN                                                        THS21570
          END                                                           THS21580
                                                                        THS21590
        SUBROUTINE CHCOEF(LS,NL,RS,NR,C,S)                              THS21600
        INTEGER NL, NR, N, IC(20), STLEN                                THS21610
        REAL*8 C(*)                                                     THS21620
        CHARACTER*5 CC                                                  THS21630
        CHARACTER*(*) LS(*), RS(*), S(*)                                THS21640
C*                                                                      THS21650
        N = NL + NR                                                     THS21660
        DO 1 I = 1, N                                                   THS21670
1       IC(I) = IDINT(C(I))                                             THS21680
        DO 2 I = 1, NL                                                  THS21690
        CALL NTS(IC(I),CC,LCC,IERR)                                     THS21700
        LST = STLEN(S(I))                                               THS21710
        LS(I) = ' '                                                     THS21720
        IF(IC(I).GT.1)THEN                                              THS21730
                LS(I)(:LCC) = CC(:LCC)                                  THS21740
                LS(I)(LCC+2:LCC+LST+2) = S(I)(:LST)                     THS21750
        ELSE                                                            THS21760
                LS(I)(:LST) = S(I)(:LST)                                THS21770
        ENDIF                                                           THS21780
2       CONTINUE                                                        THS21790
        DO 3 I = 1, NR                                                  THS21800
        CALL NTS(IC(I+NL), CC, LCC, IERR)                               THS21810
        LST = STLEN(S(I+NL))                                            THS21820
        RS(I) = ' '                                                     THS21830
        IF(IC(I+NL).GT.1)THEN                                           THS21840
                RS(I)(:LCC) = CC(:LCC)                                  THS21850
                RS(I)(LCC+2:LCC+2+LST) = S(I+NL)(:LST)                  THS21860
        ELSE                                                            THS21870
                                                                        THS21880
                                                                        THS21890
                RS(I)(:LST) = S(I+NL)(:LST)                             THS21900
        ENDIF                                                           THS21910
3       CONTINUE                                                        THS21920
        RETURN                                                          THS21930
        END                                                             THS21940
                                                                        THS21950
        SUBROUTINE NTS ( N, S, LST, IERR )                              THS21960
C*                                                                      THS21970
C*    This subroutine converts an integer into a character string       THS21980
C*                                                                      THS21990
C*    Written by:  Ed Ritter                                            THS22000
C*                 12/31/89                                             THS22010
C*                                                                      THS22020
C*    N - the integer to be converted                                   THS22030
C*    S - the resulting character string                                THS22040
C*    LST - length of the character string (numeric portion)            THS22050
C*    IERR - an error flag                                              THS22060
C*           IERR = 0 : successful  conversion                          THS22070
C*                = -1 : integer longer than character string dimension THS22080
C*                                                                      THS22090
C*                                                                      THS22100
        INTEGER N, LST, LMX, IERR, NTEST, NOFI(20)                      THS22110
        CHARACTER*(*) S                                                 THS22120
        REAL*8 TMP1, TMP2, RNUM                                         THS22130
        IERR = 0                                                        THS22140
        LMX = LEN(S)                                                    THS22150
        LST = 0                                                         THS22160
        NTEST = N                                                       THS22170
        DO 101 I = 1,LMX                                                THS22180
        IF(NTEST.GT.10)THEN                                             THS22190
                NTEST = NTEST / 10                                      THS22200
        ELSE                                                            THS22210
                NORDER = I                                              THS22220
                GOTO 8                                                  THS22230
        ENDIF                                                           THS22240
101     CONTINUE                                                        THS22250
8       CONTINUE                                                        THS22260
        TMP1 = FLOAT(N)                                                 THS22270
        DO 102 I = 2,NORDER+1                                           THS22280
        TMP2 = 10.0**(I-1)                                              THS22290
        RNUM = MOD(TMP1,TMP2)                                           THS22300
        NOFI(I-1)=IDINT((10.0*RNUM)/TMP2)                               THS22310
102     CONTINUE                                                        THS22320
        DO 104 I = 1, NORDER                                            THS22330
        IND = NORDER - I + 1                                            THS22340
        LST = LST + 1                                                   THS22350
        IF(LST.GT.LMX)THEN                                              THS22360
                IERR = -1                                               THS22370
                S = ' '                                                 THS22380
                LST = 1                                                 THS22390
                RETURN                                                  THS22400
                                                                        THS22410
                                                                        THS22420
        ENDIF                                                           THS22430
        S(LST:LST) = CHAR (NOFI(IND)+48)                                THS22440
104     CONTINUE                                                        THS22450
        RETURN                                                          THS22460
        END                                                             THS22470
                                                                        THS22480
        SUBROUTINE SAMESP(S,C,NL,NR)                                    THS22490
        IMPLICIT REAL*8 (A-H,O-Z)                                       THS22500
        CHARACTER*14 SL(5), SR(5)                                       THS22510
        CHARACTER*(*) S(*)                                              THS22520
        REAL*8 C(*)                                                     THS22530
        INTEGER NL, NR, NSL, NSR, N                                     THS22540
        DO 101  I = 1, NL                                               THS22550
101     SL(I) = S(I)                                                    THS22560
        DO 102 I = 1, NR                                                THS22570
102     SR(I) = S(I+NL)                                                 THS22580
1       CONTINUE                                                        THS22590
        DO 8 I = 1, NL                                                  THS22600
                DO 9 J = I, NL                                          THS22610
                IF(SL(I).EQ.SL(J).AND.I.NE.J)THEN                       THS22620
                        C(I) = C(I) + C(J)                              THS22630
                        C(J) = 0.0                                      THS22640
                        SL(J) = ' '                                     THS22650
                ENDIF                                                   THS22660
9               CONTINUE                                                THS22670
8       CONTINUE                                                        THS22680
2       CONTINUE                                                        THS22690
        DO 10 I = 1, NR                                                 THS22700
                DO 11 J = I, NR                                         THS22710
                IF(SR(I).EQ.SR(J).AND.I.NE.J)THEN                       THS22720
                        C(NL+I) = C(NL+I) + C(NL+J)                     THS22730
                        C(NL+J) = 0.0                                   THS22740
                        SR(J) = ' '                                     THS22750
                ENDIF                                                   THS22760
11              CONTINUE                                                THS22770
10      CONTINUE                                                        THS22780
3       CONTINUE                                                        THS22790
        DO 12 I = 1, NL                                                 THS22800
                DO 13 J = 1, NR                                         THS22810
                IF(SL(I).EQ.SR(J).AND.C(I).EQ.C(J+NL))THEN              THS22820
                        C(I) = 0.0                                      THS22830
                        C(J+NL) = 0.0                                   THS22840
                        SL(I) = ' '                                     THS22850
                        SR(I) = ' '                                     THS22860
                ELSEIF(SL(I).EQ.SR(J).AND.C(I).NE.C(J+NL))THEN          THS22870
                        DEL = C(J+NL) - C(I)                            THS22880
                        IF(DEL.GT.0.0)THEN                              THS22890
                                C(I) = 0.0                              THS22900
                                C(J+NL) = DEL                           THS22910
                                SL(I) = ' '                             THS22920
                        ELSEIF(DEL.LT.0.0)THEN                          THS22930
                                C(J+NL) = 0.0                           THS22940
                                                                        THS22950
                                                                        THS22960
                                C(I) = -1.0*DEL                         THS22970
                                SR(J) = ' '                             THS22980
                        ENDIF                                           THS22990
                ENDIF                                                   THS23000
13              CONTINUE                                                THS23010
12      CONTINUE                                                        THS23020
                        N = NL + NR                                     THS23030
                        CALL CSTRIN ( SL, NL )                          THS23040
                        CALL CSTRIN ( SR, NR )                          THS23050
                        CALL DCREAL ( C, N )                            THS23060
                                                                        THS23070
        DO 100 I = 1,NL                                                 THS23080
100     S(I) = SL(I)                                                    THS23090
        DO 200 I = 1,NR                                                 THS23100
200     S(I+NL) = SR(I)                                                 THS23110
C*      WRITE(*,'(1X,5(F5.0,1X,A,1X))')(C(I),S(I),I=1,NL)               THS23120
C*      WRITE(*,'(1X,5(F5.0,1X,A,1X))')(C(I),S(I),I=NL+1,NR+NL)         THS23130
C*      PAUSE                                                           THS23140
        RETURN                                                          THS23150
        END                                                             THS23160
                                                                        THS23170
                                                                        THS23180
      SUBROUTINE CPEXFT(NT,CPIN,T,CPZERO,CPINF,A)                       THS23190
      IMPLICIT REAL*8 (A-H,O-Z)                                         THS23200
      DIMENSION A(*),Y(30),B1(5),T(*),CPIN(*)                           THS23210
                                                                        THS23220
        T(NT+1)=6000.0                                                  THS23230
        CPIN(NT+1)=0.995*CPINF                                          THS23240
                                                                        THS23250
      DO 200 I=1,NT+1                                                   THS23260
      Y(I) = EYOFT(T(I),CPZERO,CPINF,CPIN(I))                           THS23270
200   CONTINUE                                                          THS23280
                                                                        THS23290
      CALL LLSQ(T,Y,B1,4,NT+1,SSR,AVGERR)                               THS23300
                                                                        THS23310
      DO 202 I=1,4                                                      THS23320
202   A(I)=B1(I)                             

        errmax=0.
        DO 100 I=1,NT
           error = abs(cpexp(a,t(i),cpzero,cpinf)-cpin(i))/cpin(i)
           if (error.gt.errmax) errmax = error
100     continue
        if (errmax.gt.0.045) write(25,150) errmax*100.
150     format(' WARNING:  EXP max fitting error of ',f4.1,'%')
    
*100    WRITE(*,'(1X,2(2X,F9.3))')T(I),CPIN(I)                          THS23350
*       DO 101 I=1,30                                                   THS23360
*       TIN=100.*I                                                      THS23370
*       CPOUT=CPEXP(A,TIN,CPZERO,CPINF)                                 THS23380
*       WRITE(*,'(1X,2(2X,F9.3))')TIN,CPOUT                             THS23390
*101    CONTINUE                                                        THS23400
      RETURN                                                            THS23410
      END                                                               THS23420
*--------------------------------------------------------               THS23430
                                                                        THS23440
      FUNCTION CPEXP (A,T,CPZERO,CPINF)                                 THS23450
*                                                                       THS23460
*   THIS FUNCTION RETURNS Cp OF A MOLECULE                              THS23470
*   AS CALCULATED USING THE EXP POLYNOMIALS                             THS23480
*                                                                       THS23490
      IMPLICIT REAL*8 (A-H,O-Z)                                         THS23500
      DIMENSION A(*)                                                    THS23510
C        KUF = O                                                        THS23520
C        CALL XUFLOW(KUF)                                               THS23530
      TERM=EXP(-1. * EYCALC(A,T))                                       THS23540
      CP = CPZERO + ((CPINF-CPZERO)*(1.0-TERM))                         THS23550
      CPEXP = CP                                                        THS23560
      RETURN                                                            THS23570
      END                                                               THS23580
*--------------------------------------------------------               THS23590
                                                                        THS23600
      FUNCTION EYOFT (T,CPZERO,CPINF,CPOFT)                             THS23610
*                                                                       THS23620
*    THIS FUNCTION TRANSFORMS A Cp DATA POINT INTO                      THS23630
*    THE APPROPRIATE VALYE OF THE DEPENDENT VARIABLE                    THS23640
*    FOR THE EXP POLYNOMIAL FORM.                                       THS23650
*                                                                       THS23660
      IMPLICIT REAL*8 (A-H,O-Z)                                         THS23670
        IF(CPOFT.GE.CPINF)THEN                                          THS23680
                CPINF=1.01*CPOFT                                        THS23690
        ENDIF                                                           THS23700
      Y = -1.0*LOG(1.0-((CPOFT-CPZERO)/(CPINF-CPZERO)))                 THS23710
      EYOFT=Y                                                           THS23720
      RETURN                                                            THS23730
      END                                                               THS23740
                                                                        THS23750
*---------------------------------------------------------              THS23760
      FUNCTION EYCALC (A,T)                                             THS23770
*                                                                       THS23780
*     THIS FUNCTION COMPUTES THE CURRENT VALUE OF THE                   THS23790
*     DEPENDENT VARIABLE IN TERMS OF EXP TRANSFORMATION                 THS23800
*                                                                       THS23810
      IMPLICIT REAL*8(A-H,O-Z)                                          THS23820
      DIMENSION A(*)                                                    THS23830
        COMMON/EXPAR/B                                                  THS23840
      X=T                                                               THS23850
      SUM = 0.0                                                         THS23860
        SLOPE=0.0                                                       THS23870
      DO 100 I=1,4                                                      THS23880
        IF(I.GT.1)SLOPE=SLOPE+DBLE(I-1)*(A(I)*(X**(I-2)))               THS23890
100   SUM=SUM+(A(I)*(X**(I-1)))
c      the following is unnecessary   <ayc 2/93>
c       IF(SLOPE.GT.0)THEN                                              THS23910
c               EYCALC=SUM                                              THS23920
c               B=SUM                                                   THS23930
c       ELSE                                                            THS23940
c               EYCALC=B*X                                              THS23950
c       ENDIF                                                           THS23960
      eycalc=sum
      RETURN                                                            THS23970
      END                                                               THS23980
                                                                        THS23990
                                                                        THS24000
*******************************************************************     THS24010
*******************************************************************     THS24020
                                                                        THS24030
                                                                        THS24040
      SUBROUTINE WILHOI(NT,CPIN,T,CPZERO,CPINF,A,B)                     THS24050
                                                                        THS24060
                                                                        THS24070
      IMPLICIT REAL*8 (A-H,O-Z)                                         THS24080
      DIMENSION A(*),Y(30),X(30), B1(5),T(*),CPIN(*)                    THS24090
                                                                        THS24100
        B=500.                                                          THS24110
*       WRITE(*,*)'ENTER A GUESS FOR B'                                 THS24120
*       READ(*,*)B                                                      THS24130
      DO 200 I=1,NT                                                     THS24140
      Y(I) = WYOFT(B,T(I),CPZERO,CPINF,CPIN(I))
    X(I) = T(I)/(T(I)+B)
200   CONTINUE                                                          THS24160
c   begin pey (25may04) -- fixed bug in following call to LLSQ
c      CALL LLSQ(T,Y,B1,4,NT,SSR,AVGERR)                                 THS24170
    CALL LLSQ(X,Y,B1,4,NT,SSR,AVGERR)
c   end pey
      DO 202 I=1,4                                                      THS24180
202   A(I)=B1(I)                                                        THS24190
c       DO 100 I=1,NT                                                    THS24200
c100    WRITE(*,'(1X,2(2X,F9.3))')T(I),CPIN(I)                           THS24210
c       DO 101 I=1,30                                                    THS24220
c       TIN=100.*I                                                       THS24230
c       CPOUT=CPWH(A,B,TIN,CPZERO,CPINF)                                 THS24240
c       WRITE(*,'(1X,2(2X,F9.3))')TIN,CPOUT                              THS24250
c101    CONTINUE                                                         THS24260
      RETURN                                                            THS24270
      END                                                               THS24280
*--------------------------------------------------------               THS24290
                                                                        THS24300
      FUNCTION CPWH (A,B,T,CPZERO,CPINF)                                THS24310
*                                                                       THS24320
*   THIS FUNCTION RETURNS Cp OF A MOLECULE                              THS24330
*   AS CALCULATED USING THE WILHOIT POLYNOMIALS                         THS24340
*                                                                       THS24350
      IMPLICIT REAL*8 (A-H,O-Z)                                         THS24360
      DIMENSION A(*)                                                    THS24370
      X = (T/(T+B))                                                     THS24380
      SUM = WYCALC(A,B,T)                                               THS24390
      TERM1 = (1.0+ ((X-1.0)*SUM))                                      THS24400
      TERM2 = CPINF-CPZERO                                              THS24410
      CP = CPZERO + (TERM2*(X**2)*TERM1)                                THS24420
      CPWH = CP                                                         THS24430
      RETURN                                                            THS24440
      END                                                               THS24450
                                                                        THS24460
                                                                        THS24470
*--------------------------------------------------------               THS24480
                                                                        THS24490
      FUNCTION WYOFT (B,T,CPZERO,CPINF,CPOFT)                           THS24500
*                                                                       THS24510
*    THIS FUNCTION TRANSFORMS A Cp DATA POINT INTO                      THS24520
*    THE APPROPRIATE VALUE OF THE DEPENDENT VARIABLE                    THS24530
*    FOR THE WILHOIT POLYNOMIAL FORM.                                   THS24540
*                                                                       THS24550
      IMPLICIT REAL*8 (A-H,O-Z)                                         THS24560
      X = T/(T+B)                                                       THS24570
      TOP = CPOFT - CPZERO - ((X**2)*(CPINF-CPZERO))                    THS24580
      BOTTOM = (CPINF-CPZERO)*(X**2)*(X-1.0)                            THS24590
      WYOFT=TOP/BOTTOM                                                  THS24600
      RETURN                                                            THS24610
      END                                                               THS24620
                                                                        THS24630
*--------------------------------------------------------               THS24660
                                                                        THS24670
      FUNCTION WYCALC (A,B,T)                                           THS24680
*                                                                       THS24690
*     THIS FUNCTION COMPUTES THE CURRENT VALUE OF THE                   THS24700
*     DEPENDENT VARIABLE IN TERMS OF WILHOIT'S TRANSFORMATION           THS24710
*                                                                       THS24720
      IMPLICIT REAL*8(A-H,O-Z)                                          THS24730
      DIMENSION A(*)                                                    THS24740
      X=T/(T+B)                                                         THS24750
      SUM = 0.0                                                         THS24760
      DO 100 I=1,4                                                      THS24770
100   SUM=SUM+(A(I)*(X**(I-1)))                                         THS24780
      WYCALC=SUM                                                        THS24790
      RETURN                                                            THS24800
      END                                                               THS24810
                                                                        THS24820
*---------------------------------------------------------              THS24830

      SUBROUTINE LLSQ (XIN,YIN,B,NPAR,NDATA,SSR,AVGERR)                 THS24840
      IMPLICIT REAL*8 (A-H,O-Z)                                         THS24850
        DIMENSION XIN(*),YIN(*)                                         THS24860
      DIMENSION XT(5,30),X(30,5),XTX(5,5)                               THS24870
      DIMENSION XINV(5,5),WORK(5,5),DELPCT(30)                          THS24880
      DIMENSION Y(30,1),B(5),XTY(5,1)                                   THS24890
C                                                                       THS24900
C   fill the X matrix and X transpose                                   THS24910
        DO 98 J=1,NPAR                                                  THS24920
      DO 99 K=1,NDATA                                                   THS24930
      X(K,J)=XIN(K)**(J-1)                                              THS24940
99      XT(J,K)=XIN(K)**(J-1)                                           THS24950
98      CONTINUE                                                        THS24960
                                                                        THS24970
C   COMPUTE THE MATRIX PRODUCT XTX                                      THS24980
      DO 1103 M=1,NPAR                                                  THS24990
      DO 1102 N=1,NPAR                                                  THS25000
      TSUM=0.                                                           THS25010
      DO 1101 L=1,NDATA                                                 THS25020
 1101  TSUM=TSUM+(XT(M,L)*X(L,N))                                       THS25030
 1102  XTX(M,N)=TSUM                                                    THS25040
 1103  CONTINUE                                                         THS25050
C                                                                       THS25060
C   COMPUTE THE INVERSE OF XTX                                          THS25070
      CALL INVERT(NPAR,XTX,XINV,WORK,1)                                 THS25080
C                                                                       THS25090
C   DEFINE THE VECTOR YP(I,1,1): THESE ARE THE FUNCTION VALUES          THS25100
      DO 335 I=1,NDATA                                                  THS25110
  335 Y(I,1)=YIN(I)                                                     THS25120
C                                                                       THS25130
C   COMPUTE THE MATRIX VECTOR PRODUCT XTY                               THS25140
      DO 2103 M=1,NPAR                                                  THS25150
      TSUM=0.                                                           THS25160
      DO 2101 L=1,NDATA                                                 THS25170
 2101  TSUM=TSUM+(XT(M,L)*Y(L,1))                                       THS25180
                                                                        THS25190
                                                                        THS25200
 2102  XTY(M,1)=TSUM                                                    THS25210
 2103  CONTINUE                                                         THS25220
                                                                        THS25230
C                                                                       THS25240
C   COMPUTE THE VECTOR OF LINEAR COEFFICIENTS                           THS25250
C                                                                       THS25260
      DO 3103 M=1,NPAR                                                  THS25270
      TSUM=0.                                                           THS25280
      DO 3101 L=1,NPAR                                                  THS25290
 3101  TSUM=TSUM+(XINV(M,L)*XTY(L,1))                                   THS25300
 3102  B(M)=TSUM                                                        THS25310
 3103  CONTINUE                                                         THS25320
                                                                        THS25330
      SSR=0.                                                            THS25340
      DELSUM=0.                                                         THS25350
      AVGERR=0.                                                         THS25360
C                                                                       THS25370
C   CALCULATE PREDICTED RESPONSE, AVG ABS ERROR, AND SUM OF SQUARES     THS25380
      YCALC=0.                                                          THS25390
      DO 333 I=1,NPAR                                                   THS25400
  333 YCALC=YCALC+(B(I)*X(J,I))                                         THS25410
      DEL= YIN(J)-YCALC                                                 THS25420
      IF(YIN(J).EQ.0.0)GO TO 334                                        THS25430
      DELPCT(J)=(DEL/YIN(J))*100.                                       THS25440
  334 DELSUM=DELSUM+DABS(DELPCT(J))                                     THS25450
  444 SSR=SSR+(DEL**2)                                                  THS25460
      AVGERR=DELSUM/NDATA                                               THS25470
      RETURN                                                            THS25480
      END                                                               THS25490
                                                                        THS25500
                                                                        THS25510
      SUBROUTINE INVERT(NCON,C,CIN,B,N)                                 THS25520
      IMPLICIT REAL*8 (A-H,O-Z)                                         THS25530
      DIMENSION C(5,5),CIN(5,5),B(5,5)                                  THS25540
C  CREATE THE IDENTITY MATRIX B                                         THS25550
      DO 10 I=1,NCON                                                    THS25560
      DO 20 J=1,NCON                                                    THS25570
        IF(J.EQ.I) THEN                                                 THS25580
                B(I,J)=1.0                                              THS25590
        ELSE                                                            THS25600
                B(I,J)=0.0                                              THS25610
        ENDIF                                                           THS25620
   20 CONTINUE                                                          THS25630
   10 CONTINUE                                                          THS25640
      NR=NCON                                                           THS25650
      CALL GJELIM(NCON,NR,N,C,B,CIN)                                    THS25660
      RETURN                                                            THS25670
      END                                                               THS25680
                                                                        THS25690
                                                                        THS25700
                                                                        THS25710
      SUBROUTINE GJELIM(NCON,NR,NA,C,B,CIN)                             THS25720
      IMPLICIT REAL*8 (A-H,O-Z)                                         THS25730
      DIMENSION C(5,5),B(5,5),CIN(5,5)                                  THS25740
C                                                                       THS25750
C     ARRAYS EXTENDED TO INTERFACE LITH LINEAR LEAST SQUARE ROUTINE     THS25760
C     EXTENDED TO DOUBLE PRECISION                                      THS25770
C                                                                       THS25780
      IF(NCON.EQ.1) GO TO 1120                                          THS25790
      LAST=NCON-1                                                       THS25800
      DO 1700 I=1,LAST                                                  THS25810
      M=I                                                               THS25820
      ITEMP=I+1                                                         THS25830
      DO 1100 J=ITEMP,NCON                                              THS25840
      IF(ABS(C(M,I)).GT.ABS(C(J,I))) GO TO 1100                         THS25850
      M=J                                                               THS25860
 1100 CONTINUE                                                          THS25870
      IF(ABS(C(M,I)).LT.1.E-15) GO TO 1130                              THS25880
      IF(M.EQ.I) GO TO 1400                                             THS25890
      DO 1200 J=I,NCON                                                  THS25900
      TEMP=C(M,J)                                                       THS25910
      C(M,J)=C(I,J)                                                     THS25920
 1200 C(I,J)=TEMP                                                       THS25930
      DO 1300 J=1,NR                                                    THS25940
      TEMP=B(M,J)                                                       THS25950
      B(M,J)=B(I,J)                                                     THS25960
 1300 B(I,J)=TEMP                                                       THS25970
 1400 DO 1600 J=ITEMP,NCON                                              THS25980
      TEMP=C(J,I)/C(I,I)                                                THS25990
      C(J,I)=0.0                                                        THS26000
      DO 1500 K=ITEMP,NCON                                              THS26010
 1500 C(J,K)=C(J,K)-TEMP*C(I,K)                                         THS26020
      DO 1600 K=1,NR                                                    THS26030
 1600 B(J,K)=B(J,K)-TEMP*B(I,K)                                         THS26040
 1700 CONTINUE                                                          THS26050
      IF(ABS(C(NCON,NCON)).LT.1.E-15) GO TO 1130                        THS26060
      DO 1800 J=1,NR                                                    THS26070
 1800 CIN(NCON,J)=B(NCON,J)/C(NCON,NCON)                                THS26080
      K=NCON-1                                                          THS26090
 1900 DO 1111 J=1,NR                                                    THS26100
      IT=K+1                                                            THS26110
      TEMP=0.0                                                          THS26120
      DO 1110 I=IT,NCON                                                 THS26130
 1110 TEMP=TEMP+C(K,I)*CIN(I,J)                                         THS26140
 1111 CIN(K,J)=(B(K,J)-TEMP)/C(K,K)                                     THS26150
      K=K-1                                                             THS26160
      IF(K.NE.0) GO TO 1900                                             THS26170
      RETURN                                                            THS26180
 1120 CIN(NCON,NR)=B(NCON,NR)/C(NCON,NCON)                              THS26190
      RETURN                                                            THS26200
 1130 CONTINUE                                                          THS26210
                                                                        THS26220
      WRITE(*,120)                                                      THS26230
  120 FORMAT(//10X,'*** MATRIX IS SINGULAR  ***'   )                    THS26240
      WRITE(*,*)'ERROR IN SUBROUTINE GJELIM: MATRIX FOLLOWS'            THS26250
      DO 1066 I=1,NCON                                                  THS26260
      WRITE(*,*)(C(I,J),J=1,NCON)                                       THS26270
 1066 CONTINUE                                                          THS26280
      STOP                                                              THS26290
      END                                                               THS26300
c  dmmatheu revisions
c
c $Id$
c $Author$
c $Date$
c $Source$
c $Revision$
c $Log$
c Revision 1.2  2009-01-20 18:22:49  rwest
c changed end of therfit output to "Therfit Job Complete". (must update fortran and java code simultaneously)
c
c Revision 1.1.1.1  2007/02/20 23:10:23  sandeeps
c Fortran_software contain
c 1. therfir
c 2. chemdis
c 3. reactorModel (chemkin, as modified by Paul)
c 4. ODESolver (the fortran code which RMG uses to call daspk or dassl)
c 5. fit3pbnd
c
c Revision 1.2  2002/04/18 21:29:11  dmmatheu
c version for working with XMG, courtesy of amamo, modified by dmm.
c Uses different input format.  This whole dumb code should be
c completely rewritten by another mind in another language.  It's too
c hard to maintain.
c -- dmm.
c
c Revision 1.1  2002/04/18 21:27:21  dmmatheu
c Initial revision
c
c Revision 1.4  2000/06/02 16:15:48  dmmatheu
c Working integer version checks out ok on all of Henning's molecules.
c Before adding user option for integer degeneracy fits ...
c
c Revision 1.3  2000/05/31 21:23:57  dmmatheu
c Just before swapping out THCPFI with IDCPFIT ...
c
c Revision 1.2  2000/05/30 22:58:53  dmmatheu
c Before changing definition of X in FOFX to combat bad coding style
c
c Revision 1.1  2000/05/23 15:38:30  dmmatheu
c Initial revision
c




