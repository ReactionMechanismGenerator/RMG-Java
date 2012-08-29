C Revised by John Z. WEN at MIT, Oct 2005
C Add the minumim and maximum temperatures:
C     The code calculate from 200K to 6000K
C
C Add the number of elements up to 4
C Change the species name up to 16 characters
C
C----------------------------------------------------------------------------------
C Revised by John Z. WEN at MIT, August 2005
C Change the input enthalpies in 'Kcal/mole'
C
C Modified by Richard H. West, June 2009, to read from STDIN rather than a file
C  and by Ramanan Sankaran and Richard H West (June 2012) to read continuously.
C
C**********************************************************************************
C Programmed by John Z. WEN at MIT, June 2005
C Function: to convert the thermo data taken from the group additivity calculation
C           into thermo data sets in different formats: CHEMKIN, NASA, WILH, NIST, 
C           etc.
C
C Unit: Cp: cal/mole.K
C       H : kcal/mole
C       S : cal/mole.K
C adjustment:
C: M_in: the number of input thermo data, use 'DATAGROUP' to read the input data
C***********************************************************************************
      PROGRAM LSAP
C
      implicit none
C
      INTEGER LIN,LOUT
      INTEGER I,J,J_elem,IJ
C max atoms number allowed in a species is 999
C max number of different atoms is 5 
      integer indc_elno(5,3)
      PARAMETER (LIN=5,LOUT=6) ! 5 is STDIN and 6 is STDOUT
C For input file
      CHARACTER(LEN=4) MARK, DATATYPE
      CHARACTER(LEN=48) TEXT, SNAM ! beware truncation!
      CHARACTER(LEN=9) STRUC_MOL
      CHARACTER(LEN=3) ELNO(5)
      CHARACTER(LEN=2) ENAM(5)
      DOUBLE PRECISION ROL,H_298,DLTH,S_298,MW,T_int,TINT
      DOUBLE PRECISION T_min,T_max,T_int1,TINT1,ATOMS,ROTORS
C
C: 'M' is the size of sample dataset
      INTEGER M_in
C     PARAMETER (M_in=59)
      PARAMETER (M_in=7)
      DOUBLE PRECISION TEMP,CPT,CH,CS,TEMF,CPF,CHF,CSF
      DIMENSION TEMP(M_in),CPT(M_in),CH(M_in),CS(M_in)
      DIMENSION TEMF(101),CPF(101),CHF(101),CSF(101)
C
C: THERMs contain the coefficients in different temperature ranges 
      DOUBLE PRECISION THERM1,THERM2,THERM3
      DIMENSION THERM1(20),THERM2(20),THERM3(20)
C
      COMMON /GAS/ROL,H_298,S_298,DLTH
      COMMON /INDICTOR/INDC_ELNO
      DATA ROL/1.9872D0/
C
C: For the format of input file, please ref. a sample file. 
C  This now commented out because LIN=5 is standard input 
C  and LOUT=6 is standard output, defined above.
C      OPEN (LIN, FORM='FORMATTED', STATUS='OLD',
C     1 FILE='INPUT.txt')
C      OPEN (LOUT, FORM='FORMATTED', STATUS='UNKNOWN',
C     1 FILE='OUTPUT.txt')
C
34501 CONTINUE
      MARK=''
      DATATYPE=''
      TEXT=''
      STRUC_MOL=''
      DO I=1,5
         ENAM(I)=''
         ELNO(I)=''
      ENDDO
C
      IJ=0
      J_elem=0
      do i=1,5
         do j=1,3
            indc_elno(i,j)=0
         enddo
      enddo
C 
      DO I=1,8
        IJ=IJ+1
        READ (LIN, 100, END=34511) MARK, TEXT
        IF (MARK .EQ. 'SPEC') THEN
C        WRITE(*,*) MARK
C        WRITE(*,*) TEXT
         SNAM=TEXT(1:16)
C        WRITE(*,103) MARK,SNAM
        ELSE IF (MARK .EQ. 'ELEM') THEN
         J_elem=J_elem+1
C        WRITE(*,*) TEXT
         ENAM(J_elem)=TEXT(1:(INDEX(TEXT,' ')-1))
         do J=1,3
           if (TEXT(INDEX(TEXT,' ')+J:INDEX(TEXT,' ')+J) .EQ. ' ') then
               indc_elno(j_elem,J)=1
            endif
         enddo
C 
         ELNO(J_elem)=TEXT(INDEX(TEXT,' ')+1:INDEX(TEXT,' ')+4)
C        WRITE(*,102) MARK,ENAM(J_ELEM),ELNO(J_ELEM)
        ELSE IF (MARK .EQ. 'H298') THEN 
C       We've read beyond the end of the elements, so read TEXT 
C       as a double into H_298 and exit the loop
C    WARNING: if the number in the input spans beyond the end of TEXT
C             then it will be truncated, eg: 1.345678901234567E+012
C                                            ---------------#######
C             this number could be read in as 10^12 too small!
          IF (TEXT(LEN(TEXT):LEN(TEXT)) .ne. ' ') THEN
            WRITE(*,*) "WARNING! The H_298 value is longer than ",
     1        "the string being used to read it.  ",
     2        "Truncating an exponent could give a very wrong value."
            STOP
          ENDIF
          READ(TEXT,*) H_298
          GOTO 200
        ENDIF
      ENDDO
 100  FORMAT (A4,1X,A48)  ! beware truncation
 200  CONTINUE
C
      DO I=1,5
         IF (ELNO(I) .EQ. '') THEN
            ELNO(I)='0'
         ENDIF
      ENDDO
C

C   We have already read H_298 line, so continue with S_298
C     READ (LIN, *) MARK, H_298
C     WRITE(*,*) MARK, H_298
      READ (LIN, *) MARK, S_298
C     WRITE(*,*) MARK, S_298
      READ (LIN, *) MARK, DLTH
C     WRITE(*,*) 'HEAT OF FORMATION at 298K', DLTH
      READ (LIN, *) MARK, MW
C     WRITE(*,*) 'MOLAR WEIGHT', MW
      READ (LIN, *) MARK, T_int
C     WRITE(*,*) MARK, T_int
      READ (LIN, *) MARK, T_min
C     WRITE(*,*) MARK, T_min
      READ (LIN, *) MARK, T_max
C     WRITE(*,*) MARK, T_max
      if (T_max .gt. 6000.0) then
      write(*,*) 'Warning!!!'
      write(*,*) 'The maximum T should be no greater than 6000K!'
      endif
 101  FORMAT (A4,1X,F8.1)
 102  FORMAT (A5,1X,A2,1X,A4)
 103  FORMAT (A5,1X,A8)
      READ (LIN, *) DATATYPE
C     WRITE(*,*) 'DATA FORMAT:  ', DATATYPE
C
      IF (DATATYPE .EQ. 'NASA') THEN
         READ (LIN, *) MARK, T_int1
C        WRITE(*,*) MARK, T_int1
      ELSE
         READ (LIN, 101)
         T_int1=0.D0
      ENDIF
C
      READ (LIN, *) STRUC_MOL
C     WRITE(*,*) STRUC_MOL
      READ (LIN, *) ATOMS
C     WRITE(*,*) 'NO. OF ATOMS', ATOMS
      READ (LIN, *) ROTORS
C     WRITE(*,*) 'NO. OF ROTORS', ROTORS
 104  FORMAT(F8.1)
C
C: Keywords input over
C
C: INPUT the Cp, H, and S sample dataset
C:
C
      DO I=1,M_in
         TEMP(I)=0.D0
         CPT(I) =0.D0
         CH(I)  =0.D0
         CS(I)  =0.D0
      ENDDO
C
C      CALL DATAINPUT (M_in,TEMP,CPT,CH,CS)
      CALL DATAGROUP (M_in,TEMP,CPT,CH,CS,LIN)
C
      TINT  = 0.D0
      TINT1 = 0.D0
C
C SWILT uses the enthalpy in: cal/mol
      DO I=1,M_in
         CH(I) = CH(I)*1.D3
      ENDDO
C
C First fit SWILT because it takes into account limits at 0 and infinite T
         CALL SWILT(SNAM,ENAM,ELNO,STRUC_MOL,ATOMS,ROTORS,THERM1,
     &                TEMP,CPT,CH,CS,M_in)
C
      IF (DATATYPE .EQ. 'WILH') GOTO 3000
C
C: INTERPOLATE ACCORDING TO THE SWILHOIT FORMAT
C
      DO I=1,11
C     10 points between 1 and 298
        TEMF(I) = 29.815D0*DFLOAT(I-1)
        IF (TEMF(I) .EQ. 0.D0) TEMF(I)=1.0D0
        CALL DATAFIND(THERM1,TEMF(I),CPF(I),CHF(I),CSF(I))
        CHF(I) = CHF(I)/1.D3
      ENDDO
C
      DO I=12,101
C     90 points between 298 and 6000
        TEMF(I)=298.15D0+(6000.D0-298.15D0)*DFLOAT(I-11)/90.D0
        CALL DATAFIND(THERM1,TEMF(I),CPF(I),CHF(I),CSF(I))
        CHF(I)=CHF(I)/1.D3
      ENDDO
C
C Then re-fit the other formats to these data points
C
      IF (DATATYPE .EQ. 'CHEM') THEN
         CALL CHEM(T_int,TINT,THERM1,THERM2,TEMF,CPF,CHF,CSF)
         IF (TINT .NE. 0.D0 ) T_int =TINT
C
      ELSEIF (DATATYPE .EQ. 'NASA') THEN
         CALL NASA(T_int,TINT,T_int1,TINT1,THERM1,THERM2,
     &               THERM3,TEMF,CPF,CHF,CSF)
         IF (TINT  .NE. 0.D0 ) T_int =TINT
         IF (TINT1 .NE. 0.D0 ) T_int1=TINT1
C
      ELSEIF (DATATYPE .EQ. 'NIST') THEN
         CALL NIST(T_int,TINT,THERM1,THERM2,TEMF,CPF,CHF,CSF)
         IF (TINT .NE. 0.D0 ) T_int =TINT
C
      ENDIF
C
 3000 CONTINUE
C
      CALL PRINT_THERM(T_int,T_int1,SNAM,ENAM,ELNO,MW,H_298,DLTH,
     &                   DATATYPE,LOUT,THERM1,THERM2,THERM3,J_elem,
     &                 T_min,T_max)
C
C     CLOSE(LOUT)
      WRITE(LOUT,*) 'GATPFIT_HAS_FINISHED_ONE_INPUT'
      FLUSH(LOUT)
      GOTO 34501

34511 CONTINUE
C
      END
      

C***********************************************************************
C
      SUBROUTINE DATAINPUT(M,TEMP,CPT,CH,CS)
C
      INTEGER M,I,LIN,IJ,IK,IM
      CHARACTER(LEN=4) MARK
      CHARACTER(LEN=40) TEXT
      DOUBLE PRECISION TEMP,CPT,CH,CS
      DIMENSION TEMP(*),CPT(*),CH(*),CS(*)
      DOUBLE PRECISION TEMP1,C1,CH1,CS1
      DIMENSION TEMP1(7),C1(7),CH1(7),CS1(7)
C
      DATA TEMP1/290.262D0, 393.258D0, 496.255D0, 589.888D0, 795.88D0,
     &           1001.87D0, 1488.76D0/
C
      DATA  C1/17.7725D0, 22.4645D0, 27.0142D0, 30.7109D0, 36.9668D0,
     &             41.8009D0, 49.3365D0/
C
      DATA CH1/0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0/
C
      DATA CS1/0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0/
C
      DO I=1,M
         TEMP(I)= TEMP1(I)
         CPT(I) = C1(I)
         CH(I)  = CH1(I)
         CS(I)  = CS1(I)
      ENDDO
      RETURN
C
      END
C***********************************************************************
C
      SUBROUTINE PRINT_THERM(Tint,Tint1,SNAM,ENAM,ELNO,MW,H298,DLTH,
     &                       DATATYPE,LOUT,THERM1,THERM2,THERM3,J_elem,
     &                       T_min,T_max)
C
      implicit none
C
      INTEGER LOUT,J_elem
      integer indc_elno(5,3)
      CHARACTER(LEN=4) DATATYPE
      CHARACTER(LEN=16) SNAM
      CHARACTER(LEN=3) ELNO(5)
      CHARACTER(LEN=2) ENAM(5)
      DOUBLE PRECISION THERM1,THERM2,THERM3,Tint,Tint1,H298,DLTH,MW
      DOUBLE PRECISION DLTH298,T_min,T_max
      DIMENSION THERM1(*),THERM2(*),THERM3(*)
      EXTERNAL DLTH298
      COMMON /INDICTOR/INDC_ELNO
C
      IF (DATATYPE .EQ. 'WILH') THEN
C       write(*,*) 'THE COEFFICIENTS OF WILHOIT FORM'
C       write(*,*) THERM1(1:6)
        WRITE(LOUT,*) 'The Wilhoit polynomial coefficients calculated:'
        WRITE(LOUT,200)
        WRITE(LOUT,201) THERM1(1:6)
 200  FORMAT(10X,'a1',15X,'a2',15X,'a3',15X,'a4',15X,'I',15X,'J')
 201  FORMAT(1X,E15.6,2X,E15.6,2X,E15.6,2X,E15.6,2X,E15.6,2X,E15.6)
C
      ELSEIF (DATATYPE .EQ. 'CHEM') THEN
C      write(*,*) SNAM
C
       WRITE(LOUT,*) 'The Chemkin polynomial coefficients calculated:'
      if (J_elem .eq. 2) then
       IF (INDC_ELNO(1,2) .EQ. 1) THEN
          IF (INDC_ELNO(2,2) .EQ. 1) THEN
               WRITE(LOUT,330) SNAM,ENAM(1),ELNO(1),ENAM(2),ELNO(2),
     &                        'G',T_min,T_max,Tint
          ELSEIF (INDC_ELNO(2,3) .EQ. 1) THEN
               WRITE(LOUT,331) SNAM,ENAM(1),ELNO(1),ENAM(2),ELNO(2),
     &                        'G',T_min,T_max,Tint
          ELSE
               WRITE(LOUT,332) SNAM,ENAM(1),ELNO(1),ENAM(2),ELNO(2),
     &                        'G',T_min,T_max,Tint
          ENDIF
        ELSEIF (INDC_ELNO(1,3) .EQ. 1) THEN
          IF (INDC_ELNO(2,2) .EQ. 1) THEN
               WRITE(LOUT,333) SNAM,ENAM(1),ELNO(1),ENAM(2),ELNO(2),
     &                        'G',T_min,T_max,Tint
          ELSEIF (INDC_ELNO(2,3) .EQ. 1) THEN
               WRITE(LOUT,334) SNAM,ENAM(1),ELNO(1),ENAM(2),ELNO(2),
     &                        'G',T_min,T_max,Tint
          ELSE
               WRITE(LOUT,335) SNAM,ENAM(1),ELNO(1),ENAM(2),ELNO(2),
     &                        'G',T_min,T_max,Tint
          ENDIF
        ELSE
          IF (INDC_ELNO(2,2) .EQ. 1) THEN
               WRITE(LOUT,336) SNAM,ENAM(1),ELNO(1),ENAM(2),ELNO(2),
     &                        'G',T_min,T_max,Tint
          ELSEIF (INDC_ELNO(2,3) .EQ. 1) THEN
               WRITE(LOUT,337) SNAM,ENAM(1),ELNO(1),ENAM(2),ELNO(2),
     &                        'G',T_min,T_max,Tint
          ELSE
               WRITE(LOUT,338) SNAM,ENAM(1),ELNO(1),ENAM(2),ELNO(2),
     &                        'G',T_min,T_max,Tint
          ENDIF
        ENDIF
      ENDIF
C
      if (J_elem .eq. 3) then
       IF (INDC_ELNO(1,2) .EQ. 1) THEN
          IF (INDC_ELNO(2,2) .EQ. 1) THEN
             IF (INDC_ELNO(3,2) .EQ. 1) THEN
                 WRITE(LOUT,340) SNAM,ENAM(1),ELNO(1),ENAM(2),ELNO(2),
     &                        ENAM(3),ELNO(3),'G',T_min,T_max,Tint
             ELSEIF (INDC_ELNO(3,3) .EQ. 1) THEN
                 WRITE(LOUT,341) SNAM,ENAM(1),ELNO(1),ENAM(2),ELNO(2),
     &                        ENAM(3),ELNO(3),'G',T_min,T_max,Tint
             ELSE
                 WRITE(LOUT,342) SNAM,ENAM(1),ELNO(1),ENAM(2),ELNO(2),
     &                        ENAM(3),ELNO(3),'G',T_min,T_max,Tint
             ENDIF
          ELSEIF (INDC_ELNO(2,3) .EQ. 1) THEN
             IF (INDC_ELNO(3,2) .EQ. 1) THEN
                 WRITE(LOUT,343) SNAM,ENAM(1),ELNO(1),ENAM(2),ELNO(2),
     &                        ENAM(3),ELNO(3),'G',T_min,T_max,Tint
             ELSEIF (INDC_ELNO(3,3) .EQ. 1) THEN
                 WRITE(LOUT,344) SNAM,ENAM(1),ELNO(1),ENAM(2),ELNO(2),
     &                        ENAM(3),ELNO(3),'G',T_min,T_max,Tint
             ELSE
                 WRITE(LOUT,345) SNAM,ENAM(1),ELNO(1),ENAM(2),ELNO(2),
     &                        ENAM(3),ELNO(3),'G',T_min,T_max,Tint
             ENDIF
          ELSE
             IF (INDC_ELNO(3,2) .EQ. 1) THEN
                 WRITE(LOUT,346) SNAM,ENAM(1),ELNO(1),ENAM(2),ELNO(2),
     &                        ENAM(3),ELNO(3),'G',T_min,T_max,Tint
             ELSEIF (INDC_ELNO(3,3) .EQ. 1) THEN
                 WRITE(LOUT,347) SNAM,ENAM(1),ELNO(1),ENAM(2),ELNO(2),
     &                        ENAM(3),ELNO(3),'G',T_min,T_max,Tint
             ELSE
                 WRITE(LOUT,348) SNAM,ENAM(1),ELNO(1),ENAM(2),ELNO(2),
     &                        ENAM(3),ELNO(3),'G',T_min,T_max,Tint
             ENDIF
          ENDIF
       ELSEIF (INDC_ELNO(1,3) .EQ. 1) THEN
          IF (INDC_ELNO(2,2) .EQ. 1) THEN
             IF (INDC_ELNO(3,2) .EQ. 1) THEN
                 WRITE(LOUT,350) SNAM,ENAM(1),ELNO(1),ENAM(2),ELNO(2),
     &                        ENAM(3),ELNO(3),'G',T_min,T_max,Tint
             ELSEIF (INDC_ELNO(3,3) .EQ. 1) THEN
                 WRITE(LOUT,351) SNAM,ENAM(1),ELNO(1),ENAM(2),ELNO(2),
     &                        ENAM(3),ELNO(3),'G',T_min,T_max,Tint
             ELSE
                 WRITE(LOUT,352) SNAM,ENAM(1),ELNO(1),ENAM(2),ELNO(2),
     &                        ENAM(3),ELNO(3),'G',T_min,T_max,Tint
             ENDIF
          ELSEIF (INDC_ELNO(2,3) .EQ. 1) THEN
             IF (INDC_ELNO(3,2) .EQ. 1) THEN
                 WRITE(LOUT,353) SNAM,ENAM(1),ELNO(1),ENAM(2),ELNO(2),
     &                        ENAM(3),ELNO(3),'G',T_min,T_max,Tint
             ELSEIF (INDC_ELNO(3,3) .EQ. 1) THEN
                 WRITE(LOUT,354) SNAM,ENAM(1),ELNO(1),ENAM(2),ELNO(2),
     &                        ENAM(3),ELNO(3),'G',T_min,T_max,Tint
             ELSE
                 WRITE(LOUT,355) SNAM,ENAM(1),ELNO(1),ENAM(2),ELNO(2),
     &                        ENAM(3),ELNO(3),'G',T_min,T_max,Tint
             ENDIF
          ELSE
             IF (INDC_ELNO(3,2) .EQ. 1) THEN
                 WRITE(LOUT,356) SNAM,ENAM(1),ELNO(1),ENAM(2),ELNO(2),
     &                        ENAM(3),ELNO(3),'G',T_min,T_max,Tint
             ELSEIF (INDC_ELNO(3,3) .EQ. 1) THEN
                 WRITE(LOUT,357) SNAM,ENAM(1),ELNO(1),ENAM(2),ELNO(2),
     &                        ENAM(3),ELNO(3),'G',T_min,T_max,Tint
             ELSE
                 WRITE(LOUT,358) SNAM,ENAM(1),ELNO(1),ENAM(2),ELNO(2),
     &                        ENAM(3),ELNO(3),'G',T_min,T_max,Tint
             ENDIF
          ENDIF
       ELSE
          IF (INDC_ELNO(2,2) .EQ. 1) THEN
             IF (INDC_ELNO(3,2) .EQ. 1) THEN
                 WRITE(LOUT,360) SNAM,ENAM(1),ELNO(1),ENAM(2),ELNO(2),
     &                        ENAM(3),ELNO(3),'G',T_min,T_max,Tint
             ELSEIF (INDC_ELNO(3,3) .EQ. 1) THEN
                 WRITE(LOUT,361) SNAM,ENAM(1),ELNO(1),ENAM(2),ELNO(2),
     &                        ENAM(3),ELNO(3),'G',T_min,T_max,Tint
             ELSE
                 WRITE(LOUT,362) SNAM,ENAM(1),ELNO(1),ENAM(2),ELNO(2),
     &                        ENAM(3),ELNO(3),'G',T_min,T_max,Tint
             ENDIF
          ELSEIF (INDC_ELNO(2,3) .EQ. 1) THEN
             IF (INDC_ELNO(3,2) .EQ. 1) THEN
                 WRITE(LOUT,363) SNAM,ENAM(1),ELNO(1),ENAM(2),ELNO(2),
     &                        ENAM(3),ELNO(3),'G',T_min,T_max,Tint
             ELSEIF (INDC_ELNO(3,3) .EQ. 1) THEN
                 WRITE(LOUT,364) SNAM,ENAM(1),ELNO(1),ENAM(2),ELNO(2),
     &                        ENAM(3),ELNO(3),'G',T_min,T_max,Tint
             ELSE
                 WRITE(LOUT,365) SNAM,ENAM(1),ELNO(1),ENAM(2),ELNO(2),
     &                        ENAM(3),ELNO(3),'G',T_min,T_max,Tint
             ENDIF
          ELSE
             IF (INDC_ELNO(3,2) .EQ. 1) THEN
                 WRITE(LOUT,366) SNAM,ENAM(1),ELNO(1),ENAM(2),ELNO(2),
     &                        ENAM(3),ELNO(3),'G',T_min,T_max,Tint
             ELSEIF (INDC_ELNO(3,3) .EQ. 1) THEN
                 WRITE(LOUT,367) SNAM,ENAM(1),ELNO(1),ENAM(2),ELNO(2),
     &                        ENAM(3),ELNO(3),'G',T_min,T_max,Tint
             ELSE
                 WRITE(LOUT,368) SNAM,ENAM(1),ELNO(1),ENAM(2),ELNO(2),
     &                        ENAM(3),ELNO(3),'G',T_min,T_max,Tint
             ENDIF
          ENDIF
       ENDIF
      ENDIF
      if (J_elem .eq. 4) then
       IF (INDC_ELNO(1,2) .EQ. 1) THEN
          IF (INDC_ELNO(2,2) .EQ. 1) THEN
             IF (INDC_ELNO(3,2) .EQ. 1) THEN
               IF (INDC_ELNO(4,2) .EQ. 1) THEN
                   WRITE(LOUT,370) SNAM,ENAM(1),ELNO(1),ENAM(2),ELNO(2),
     &             ENAM(3),ELNO(3),ENAM(4),ELNO(4),'G',T_min,T_max,Tint
               ELSEIF (INDC_ELNO(4,3) .EQ. 1) THEN
                   WRITE(LOUT,371) SNAM,ENAM(1),ELNO(1),ENAM(2),ELNO(2),
     &             ENAM(3),ELNO(3),ENAM(4),ELNO(4),'G',T_min,T_max,Tint
               ELSE
                   WRITE(LOUT,372) SNAM,ENAM(1),ELNO(1),ENAM(2),ELNO(2),
     &             ENAM(3),ELNO(3),ENAM(4),ELNO(4),'G',T_min,T_max,Tint
               ENDIF
             ELSEIF (INDC_ELNO(3,3) .EQ. 1) THEN
               IF (INDC_ELNO(4,2) .EQ. 1) THEN
                   WRITE(LOUT,373) SNAM,ENAM(1),ELNO(1),ENAM(2),ELNO(2),
     &             ENAM(3),ELNO(3),ENAM(4),ELNO(4),'G',T_min,T_max,Tint
               ELSEIF (INDC_ELNO(4,3) .EQ. 1) THEN
                   WRITE(LOUT,374) SNAM,ENAM(1),ELNO(1),ENAM(2),ELNO(2),
     &             ENAM(3),ELNO(3),ENAM(4),ELNO(4),'G',T_min,T_max,Tint
               ELSE
                   WRITE(LOUT,375) SNAM,ENAM(1),ELNO(1),ENAM(2),ELNO(2),
     &             ENAM(3),ELNO(3),ENAM(4),ELNO(4),'G',T_min,T_max,Tint
               ENDIF
             ELSE
               IF (INDC_ELNO(4,2) .EQ. 1) THEN
                   WRITE(LOUT,376) SNAM,ENAM(1),ELNO(1),ENAM(2),ELNO(2),
     &             ENAM(3),ELNO(3),ENAM(4),ELNO(4),'G',T_min,T_max,Tint
               ELSEIF (INDC_ELNO(4,3) .EQ. 1) THEN
                   WRITE(LOUT,377) SNAM,ENAM(1),ELNO(1),ENAM(2),ELNO(2),
     &             ENAM(3),ELNO(3),ENAM(4),ELNO(4),'G',T_min,T_max,Tint
               ELSE
                   WRITE(LOUT,378) SNAM,ENAM(1),ELNO(1),ENAM(2),ELNO(2),
     &             ENAM(3),ELNO(3),ENAM(4),ELNO(4),'G',T_min,T_max,Tint
               ENDIF
             ENDIF
          ELSEIF (INDC_ELNO(2,3) .EQ. 1) THEN
             IF (INDC_ELNO(3,2) .EQ. 1) THEN
               IF (INDC_ELNO(4,2) .EQ. 1) THEN
                   WRITE(LOUT,380) SNAM,ENAM(1),ELNO(1),ENAM(2),ELNO(2),
     &             ENAM(3),ELNO(3),ENAM(4),ELNO(4),'G',T_min,T_max,Tint
               ELSEIF (INDC_ELNO(4,3) .EQ. 1) THEN
                   WRITE(LOUT,381) SNAM,ENAM(1),ELNO(1),ENAM(2),ELNO(2),
     &             ENAM(3),ELNO(3),ENAM(4),ELNO(4),'G',T_min,T_max,Tint
               ELSE
                   WRITE(LOUT,382) SNAM,ENAM(1),ELNO(1),ENAM(2),ELNO(2),
     &             ENAM(3),ELNO(3),ENAM(4),ELNO(4),'G',T_min,T_max,Tint
               ENDIF
             ELSEIF (INDC_ELNO(3,3) .EQ. 1) THEN
               IF (INDC_ELNO(4,2) .EQ. 1) THEN
                   WRITE(LOUT,383) SNAM,ENAM(1),ELNO(1),ENAM(2),ELNO(2),
     &             ENAM(3),ELNO(3),ENAM(4),ELNO(4),'G',T_min,T_max,Tint
               ELSEIF (INDC_ELNO(4,3) .EQ. 1) THEN
                   WRITE(LOUT,384) SNAM,ENAM(1),ELNO(1),ENAM(2),ELNO(2),
     &             ENAM(3),ELNO(3),ENAM(4),ELNO(4),'G',T_min,T_max,Tint
               ELSE
                   WRITE(LOUT,385) SNAM,ENAM(1),ELNO(1),ENAM(2),ELNO(2),
     &             ENAM(3),ELNO(3),ENAM(4),ELNO(4),'G',T_min,T_max,Tint
               ENDIF
             ELSE
               IF (INDC_ELNO(4,2) .EQ. 1) THEN
                   WRITE(LOUT,386) SNAM,ENAM(1),ELNO(1),ENAM(2),ELNO(2),
     &             ENAM(3),ELNO(3),ENAM(4),ELNO(4),'G',T_min,T_max,Tint
               ELSEIF (INDC_ELNO(4,3) .EQ. 1) THEN
                   WRITE(LOUT,387) SNAM,ENAM(1),ELNO(1),ENAM(2),ELNO(2),
     &             ENAM(3),ELNO(3),ENAM(4),ELNO(4),'G',T_min,T_max,Tint
               ELSE
                   WRITE(LOUT,388) SNAM,ENAM(1),ELNO(1),ENAM(2),ELNO(2),
     &             ENAM(3),ELNO(3),ENAM(4),ELNO(4),'G',T_min,T_max,Tint
               ENDIF
             ENDIF
          ELSE
             IF (INDC_ELNO(3,2) .EQ. 1) THEN
               IF (INDC_ELNO(4,2) .EQ. 1) THEN
                   WRITE(LOUT,390) SNAM,ENAM(1),ELNO(1),ENAM(2),ELNO(2),
     &             ENAM(3),ELNO(3),ENAM(4),ELNO(4),'G',T_min,T_max,Tint
               ELSEIF (INDC_ELNO(4,3) .EQ. 1) THEN
                   WRITE(LOUT,391) SNAM,ENAM(1),ELNO(1),ENAM(2),ELNO(2),
     &             ENAM(3),ELNO(3),ENAM(4),ELNO(4),'G',T_min,T_max,Tint
               ELSE
                   WRITE(LOUT,392) SNAM,ENAM(1),ELNO(1),ENAM(2),ELNO(2),
     &             ENAM(3),ELNO(3),ENAM(4),ELNO(4),'G',T_min,T_max,Tint
               ENDIF
             ELSEIF (INDC_ELNO(3,3) .EQ. 1) THEN
               IF (INDC_ELNO(4,2) .EQ. 1) THEN
                   WRITE(LOUT,393) SNAM,ENAM(1),ELNO(1),ENAM(2),ELNO(2),
     &             ENAM(3),ELNO(3),ENAM(4),ELNO(4),'G',T_min,T_max,Tint
               ELSEIF (INDC_ELNO(4,3) .EQ. 1) THEN
                   WRITE(LOUT,394) SNAM,ENAM(1),ELNO(1),ENAM(2),ELNO(2),
     &             ENAM(3),ELNO(3),ENAM(4),ELNO(4),'G',T_min,T_max,Tint
               ELSE
                   WRITE(LOUT,395) SNAM,ENAM(1),ELNO(1),ENAM(2),ELNO(2),
     &             ENAM(3),ELNO(3),ENAM(4),ELNO(4),'G',T_min,T_max,Tint
               ENDIF
             ELSE
               IF (INDC_ELNO(4,2) .EQ. 1) THEN
                   WRITE(LOUT,396) SNAM,ENAM(1),ELNO(1),ENAM(2),ELNO(2),
     &             ENAM(3),ELNO(3),ENAM(4),ELNO(4),'G',T_min,T_max,Tint
               ELSEIF (INDC_ELNO(4,3) .EQ. 1) THEN
                   WRITE(LOUT,397) SNAM,ENAM(1),ELNO(1),ENAM(2),ELNO(2),
     &             ENAM(3),ELNO(3),ENAM(4),ELNO(4),'G',T_min,T_max,Tint
               ELSE
                   WRITE(LOUT,398) SNAM,ENAM(1),ELNO(1),ENAM(2),ELNO(2),
     &             ENAM(3),ELNO(3),ENAM(4),ELNO(4),'G',T_min,T_max,Tint
               ENDIF
             ENDIF
          ENDIF
       ELSEIF (INDC_ELNO(1,3) .EQ. 1) THEN
          IF (INDC_ELNO(2,2) .EQ. 1) THEN
             IF (INDC_ELNO(3,2) .EQ. 1) THEN
               IF (INDC_ELNO(4,2) .EQ. 1) THEN
                   WRITE(LOUT,470) SNAM,ENAM(1),ELNO(1),ENAM(2),ELNO(2),
     &             ENAM(3),ELNO(3),ENAM(4),ELNO(4),'G',T_min,T_max,Tint
               ELSEIF (INDC_ELNO(4,3) .EQ. 1) THEN
                   WRITE(LOUT,471) SNAM,ENAM(1),ELNO(1),ENAM(2),ELNO(2),
     &             ENAM(3),ELNO(3),ENAM(4),ELNO(4),'G',T_min,T_max,Tint
               ELSE
                   WRITE(LOUT,472) SNAM,ENAM(1),ELNO(1),ENAM(2),ELNO(2),
     &             ENAM(3),ELNO(3),ENAM(4),ELNO(4),'G',T_min,T_max,Tint
               ENDIF
             ELSEIF (INDC_ELNO(3,3) .EQ. 1) THEN
               IF (INDC_ELNO(4,2) .EQ. 1) THEN
                   WRITE(LOUT,473) SNAM,ENAM(1),ELNO(1),ENAM(2),ELNO(2),
     &             ENAM(3),ELNO(3),ENAM(4),ELNO(4),'G',T_min,T_max,Tint
               ELSEIF (INDC_ELNO(4,3) .EQ. 1) THEN
                   WRITE(LOUT,474) SNAM,ENAM(1),ELNO(1),ENAM(2),ELNO(2),
     &             ENAM(3),ELNO(3),ENAM(4),ELNO(4),'G',T_min,T_max,Tint
               ELSE
                   WRITE(LOUT,475) SNAM,ENAM(1),ELNO(1),ENAM(2),ELNO(2),
     &             ENAM(3),ELNO(3),ENAM(4),ELNO(4),'G',T_min,T_max,Tint
               ENDIF
             ELSE
               IF (INDC_ELNO(4,2) .EQ. 1) THEN
                   WRITE(LOUT,476) SNAM,ENAM(1),ELNO(1),ENAM(2),ELNO(2),
     &             ENAM(3),ELNO(3),ENAM(4),ELNO(4),'G',T_min,T_max,Tint
               ELSEIF (INDC_ELNO(4,3) .EQ. 1) THEN
                   WRITE(LOUT,477) SNAM,ENAM(1),ELNO(1),ENAM(2),ELNO(2),
     &             ENAM(3),ELNO(3),ENAM(4),ELNO(4),'G',T_min,T_max,Tint
               ELSE
                   WRITE(LOUT,478) SNAM,ENAM(1),ELNO(1),ENAM(2),ELNO(2),
     &             ENAM(3),ELNO(3),ENAM(4),ELNO(4),'G',T_min,T_max,Tint
               ENDIF
             ENDIF
          ELSEIF (INDC_ELNO(2,3) .EQ. 1) THEN
             IF (INDC_ELNO(3,2) .EQ. 1) THEN
               IF (INDC_ELNO(4,2) .EQ. 1) THEN
                   WRITE(LOUT,480) SNAM,ENAM(1),ELNO(1),ENAM(2),ELNO(2),
     &             ENAM(3),ELNO(3),ENAM(4),ELNO(4),'G',T_min,T_max,Tint
               ELSEIF (INDC_ELNO(4,3) .EQ. 1) THEN
                   WRITE(LOUT,481) SNAM,ENAM(1),ELNO(1),ENAM(2),ELNO(2),
     &             ENAM(3),ELNO(3),ENAM(4),ELNO(4),'G',T_min,T_max,Tint
               ELSE
                   WRITE(LOUT,482) SNAM,ENAM(1),ELNO(1),ENAM(2),ELNO(2),
     &             ENAM(3),ELNO(3),ENAM(4),ELNO(4),'G',T_min,T_max,Tint
               ENDIF
             ELSEIF (INDC_ELNO(3,3) .EQ. 1) THEN
               IF (INDC_ELNO(4,2) .EQ. 1) THEN
                   WRITE(LOUT,483) SNAM,ENAM(1),ELNO(1),ENAM(2),ELNO(2),
     &             ENAM(3),ELNO(3),ENAM(4),ELNO(4),'G',T_min,T_max,Tint
               ELSEIF (INDC_ELNO(4,3) .EQ. 1) THEN
                   WRITE(LOUT,484) SNAM,ENAM(1),ELNO(1),ENAM(2),ELNO(2),
     &             ENAM(3),ELNO(3),ENAM(4),ELNO(4),'G',T_min,T_max,Tint
               ELSE
                   WRITE(LOUT,485) SNAM,ENAM(1),ELNO(1),ENAM(2),ELNO(2),
     &             ENAM(3),ELNO(3),ENAM(4),ELNO(4),'G',T_min,T_max,Tint
               ENDIF
             ELSE
               IF (INDC_ELNO(4,2) .EQ. 1) THEN
                   WRITE(LOUT,486) SNAM,ENAM(1),ELNO(1),ENAM(2),ELNO(2),
     &             ENAM(3),ELNO(3),ENAM(4),ELNO(4),'G',T_min,T_max,Tint
               ELSEIF (INDC_ELNO(4,3) .EQ. 1) THEN
                   WRITE(LOUT,487) SNAM,ENAM(1),ELNO(1),ENAM(2),ELNO(2),
     &             ENAM(3),ELNO(3),ENAM(4),ELNO(4),'G',T_min,T_max,Tint
               ELSE
                   WRITE(LOUT,488) SNAM,ENAM(1),ELNO(1),ENAM(2),ELNO(2),
     &             ENAM(3),ELNO(3),ENAM(4),ELNO(4),'G',T_min,T_max,Tint
               ENDIF
             ENDIF
          ELSE
             IF (INDC_ELNO(3,2) .EQ. 1) THEN
               IF (INDC_ELNO(4,2) .EQ. 1) THEN
                   WRITE(LOUT,490) SNAM,ENAM(1),ELNO(1),ENAM(2),ELNO(2),
     &             ENAM(3),ELNO(3),ENAM(4),ELNO(4),'G',T_min,T_max,Tint
               ELSEIF (INDC_ELNO(4,3) .EQ. 1) THEN
                   WRITE(LOUT,491) SNAM,ENAM(1),ELNO(1),ENAM(2),ELNO(2),
     &             ENAM(3),ELNO(3),ENAM(4),ELNO(4),'G',T_min,T_max,Tint
               ELSE
                   WRITE(LOUT,492) SNAM,ENAM(1),ELNO(1),ENAM(2),ELNO(2),
     &             ENAM(3),ELNO(3),ENAM(4),ELNO(4),'G',T_min,T_max,Tint
               ENDIF
             ELSEIF (INDC_ELNO(3,3) .EQ. 1) THEN
               IF (INDC_ELNO(4,2) .EQ. 1) THEN
                   WRITE(LOUT,493) SNAM,ENAM(1),ELNO(1),ENAM(2),ELNO(2),
     &             ENAM(3),ELNO(3),ENAM(4),ELNO(4),'G',T_min,T_max,Tint
               ELSEIF (INDC_ELNO(4,3) .EQ. 1) THEN
                   WRITE(LOUT,494) SNAM,ENAM(1),ELNO(1),ENAM(2),ELNO(2),
     &             ENAM(3),ELNO(3),ENAM(4),ELNO(4),'G',T_min,T_max,Tint
               ELSE
                   WRITE(LOUT,495) SNAM,ENAM(1),ELNO(1),ENAM(2),ELNO(2),
     &             ENAM(3),ELNO(3),ENAM(4),ELNO(4),'G',T_min,T_max,Tint
               ENDIF
             ELSE
               IF (INDC_ELNO(4,2) .EQ. 1) THEN
                   WRITE(LOUT,496) SNAM,ENAM(1),ELNO(1),ENAM(2),ELNO(2),
     &             ENAM(3),ELNO(3),ENAM(4),ELNO(4),'G',T_min,T_max,Tint
               ELSEIF (INDC_ELNO(4,3) .EQ. 1) THEN
                   WRITE(LOUT,497) SNAM,ENAM(1),ELNO(1),ENAM(2),ELNO(2),
     &             ENAM(3),ELNO(3),ENAM(4),ELNO(4),'G',T_min,T_max,Tint
               ELSE
                   WRITE(LOUT,498) SNAM,ENAM(1),ELNO(1),ENAM(2),ELNO(2),
     &             ENAM(3),ELNO(3),ENAM(4),ELNO(4),'G',T_min,T_max,Tint
               ENDIF
             ENDIF
          ENDIF
       ELSE
          IF (INDC_ELNO(2,2) .EQ. 1) THEN
             IF (INDC_ELNO(3,2) .EQ. 1) THEN
               IF (INDC_ELNO(4,2) .EQ. 1) THEN
                   WRITE(LOUT,570) SNAM,ENAM(1),ELNO(1),ENAM(2),ELNO(2),
     &             ENAM(3),ELNO(3),ENAM(4),ELNO(4),'G',T_min,T_max,Tint
               ELSEIF (INDC_ELNO(4,3) .EQ. 1) THEN
                   WRITE(LOUT,571) SNAM,ENAM(1),ELNO(1),ENAM(2),ELNO(2),
     &             ENAM(3),ELNO(3),ENAM(4),ELNO(4),'G',T_min,T_max,Tint
               ELSE
                   WRITE(LOUT,572) SNAM,ENAM(1),ELNO(1),ENAM(2),ELNO(2),
     &             ENAM(3),ELNO(3),ENAM(4),ELNO(4),'G',T_min,T_max,Tint
               ENDIF
             ELSEIF (INDC_ELNO(3,3) .EQ. 1) THEN
               IF (INDC_ELNO(4,2) .EQ. 1) THEN
                   WRITE(LOUT,573) SNAM,ENAM(1),ELNO(1),ENAM(2),ELNO(2),
     &             ENAM(3),ELNO(3),ENAM(4),ELNO(4),'G',T_min,T_max,Tint
               ELSEIF (INDC_ELNO(4,3) .EQ. 1) THEN
                   WRITE(LOUT,574) SNAM,ENAM(1),ELNO(1),ENAM(2),ELNO(2),
     &             ENAM(3),ELNO(3),ENAM(4),ELNO(4),'G',T_min,T_max,Tint
               ELSE
                   WRITE(LOUT,575) SNAM,ENAM(1),ELNO(1),ENAM(2),ELNO(2),
     &             ENAM(3),ELNO(3),ENAM(4),ELNO(4),'G',T_min,T_max,Tint
               ENDIF
             ELSE
               IF (INDC_ELNO(4,2) .EQ. 1) THEN
                   WRITE(LOUT,576) SNAM,ENAM(1),ELNO(1),ENAM(2),ELNO(2),
     &             ENAM(3),ELNO(3),ENAM(4),ELNO(4),'G',T_min,T_max,Tint
               ELSEIF (INDC_ELNO(4,3) .EQ. 1) THEN
                   WRITE(LOUT,577) SNAM,ENAM(1),ELNO(1),ENAM(2),ELNO(2),
     &             ENAM(3),ELNO(3),ENAM(4),ELNO(4),'G',T_min,T_max,Tint
               ELSE
                   WRITE(LOUT,578) SNAM,ENAM(1),ELNO(1),ENAM(2),ELNO(2),
     &             ENAM(3),ELNO(3),ENAM(4),ELNO(4),'G',T_min,T_max,Tint
               ENDIF
             ENDIF
          ELSEIF (INDC_ELNO(2,3) .EQ. 1) THEN
             IF (INDC_ELNO(3,2) .EQ. 1) THEN
               IF (INDC_ELNO(4,2) .EQ. 1) THEN
                   WRITE(LOUT,580) SNAM,ENAM(1),ELNO(1),ENAM(2),ELNO(2),
     &             ENAM(3),ELNO(3),ENAM(4),ELNO(4),'G',T_min,T_max,Tint
               ELSEIF (INDC_ELNO(4,3) .EQ. 1) THEN
                   WRITE(LOUT,581) SNAM,ENAM(1),ELNO(1),ENAM(2),ELNO(2),
     &             ENAM(3),ELNO(3),ENAM(4),ELNO(4),'G',T_min,T_max,Tint
               ELSE
                   WRITE(LOUT,582) SNAM,ENAM(1),ELNO(1),ENAM(2),ELNO(2),
     &             ENAM(3),ELNO(3),ENAM(4),ELNO(4),'G',T_min,T_max,Tint
               ENDIF
             ELSEIF (INDC_ELNO(3,3) .EQ. 1) THEN
               IF (INDC_ELNO(4,2) .EQ. 1) THEN
                   WRITE(LOUT,583) SNAM,ENAM(1),ELNO(1),ENAM(2),ELNO(2),
     &             ENAM(3),ELNO(3),ENAM(4),ELNO(4),'G',T_min,T_max,Tint
               ELSEIF (INDC_ELNO(4,3) .EQ. 1) THEN
                   WRITE(LOUT,584) SNAM,ENAM(1),ELNO(1),ENAM(2),ELNO(2),
     &             ENAM(3),ELNO(3),ENAM(4),ELNO(4),'G',T_min,T_max,Tint
               ELSE
                   WRITE(LOUT,585) SNAM,ENAM(1),ELNO(1),ENAM(2),ELNO(2),
     &             ENAM(3),ELNO(3),ENAM(4),ELNO(4),'G',T_min,T_max,Tint
               ENDIF
             ELSE
               IF (INDC_ELNO(4,2) .EQ. 1) THEN
                   WRITE(LOUT,586) SNAM,ENAM(1),ELNO(1),ENAM(2),ELNO(2),
     &             ENAM(3),ELNO(3),ENAM(4),ELNO(4),'G',T_min,T_max,Tint
               ELSEIF (INDC_ELNO(4,3) .EQ. 1) THEN
                   WRITE(LOUT,587) SNAM,ENAM(1),ELNO(1),ENAM(2),ELNO(2),
     &             ENAM(3),ELNO(3),ENAM(4),ELNO(4),'G',T_min,T_max,Tint
               ELSE
                   WRITE(LOUT,588) SNAM,ENAM(1),ELNO(1),ENAM(2),ELNO(2),
     &             ENAM(3),ELNO(3),ENAM(4),ELNO(4),'G',T_min,T_max,Tint
               ENDIF
             ENDIF
          ELSE
             IF (INDC_ELNO(3,2) .EQ. 1) THEN
               IF (INDC_ELNO(4,2) .EQ. 1) THEN
                   WRITE(LOUT,590) SNAM,ENAM(1),ELNO(1),ENAM(2),ELNO(2),
     &             ENAM(3),ELNO(3),ENAM(4),ELNO(4),'G',T_min,T_max,Tint
               ELSEIF (INDC_ELNO(4,3) .EQ. 1) THEN
                   WRITE(LOUT,591) SNAM,ENAM(1),ELNO(1),ENAM(2),ELNO(2),
     &             ENAM(3),ELNO(3),ENAM(4),ELNO(4),'G',T_min,T_max,Tint
               ELSE
                   WRITE(LOUT,592) SNAM,ENAM(1),ELNO(1),ENAM(2),ELNO(2),
     &             ENAM(3),ELNO(3),ENAM(4),ELNO(4),'G',T_min,T_max,Tint
               ENDIF
             ELSEIF (INDC_ELNO(3,3) .EQ. 1) THEN
               IF (INDC_ELNO(4,2) .EQ. 1) THEN
                   WRITE(LOUT,593) SNAM,ENAM(1),ELNO(1),ENAM(2),ELNO(2),
     &             ENAM(3),ELNO(3),ENAM(4),ELNO(4),'G',T_min,T_max,Tint
               ELSEIF (INDC_ELNO(4,3) .EQ. 1) THEN
                   WRITE(LOUT,594) SNAM,ENAM(1),ELNO(1),ENAM(2),ELNO(2),
     &             ENAM(3),ELNO(3),ENAM(4),ELNO(4),'G',T_min,T_max,Tint
               ELSE
                   WRITE(LOUT,595) SNAM,ENAM(1),ELNO(1),ENAM(2),ELNO(2),
     &             ENAM(3),ELNO(3),ENAM(4),ELNO(4),'G',T_min,T_max,Tint
               ENDIF
             ELSE
               IF (INDC_ELNO(4,2) .EQ. 1) THEN
                   WRITE(LOUT,596) SNAM,ENAM(1),ELNO(1),ENAM(2),ELNO(2),
     &             ENAM(3),ELNO(3),ENAM(4),ELNO(4),'G',T_min,T_max,Tint
               ELSEIF (INDC_ELNO(4,3) .EQ. 1) THEN
                   WRITE(LOUT,597) SNAM,ENAM(1),ELNO(1),ENAM(2),ELNO(2),
     &             ENAM(3),ELNO(3),ENAM(4),ELNO(4),'G',T_min,T_max,Tint
               ELSE
                   WRITE(LOUT,598) SNAM,ENAM(1),ELNO(1),ENAM(2),ELNO(2),
     &             ENAM(3),ELNO(3),ENAM(4),ELNO(4),'G',T_min,T_max,Tint
               ENDIF
             ENDIF
          ENDIF
       ENDIF
      ENDIF
      IF (J_elem .eq. 5) THEN
            WRITE(LOUT,601) SNAM,T_min,T_max,Tint
            WRITE(LOUT,602) ENAM(1),ELNO(1),ENAM(2),ELNO(2),
     &      ENAM(3),ELNO(3),ENAM(4),ELNO(4),ENAM(5),ELNO(5)
      ENDIF
      WRITE(LOUT,301) THERM2(1),THERM2(2),THERM2(3),THERM2(4),THERM2(5)
       WRITE(LOUT,302) THERM2(6),THERM2(7),THERM1(1),THERM1(2),THERM1(3)
       WRITE(LOUT,303) THERM1(4),THERM1(5),THERM1(6),THERM1(7)
C
 330  FORMAT(A16,8X,A2,2X,A1,A2,2X,A1,10X,A1,3X,F7.3,2X,F8.3,2X,F8.3,4X,
     &      '1')
 331  FORMAT(A16,8X,A2,2X,A1,A2,1X,A2,10X,A1,3X,F7.3,2X,F8.3,2X,F8.3,4X,
     &      '1')
 332  FORMAT(A16,8X,A2,2X,A1,A3,A2,10X,A1,3X,F7.3,2X,F8.3,2X,F8.3,4X,
     &      '1')
 333  FORMAT(A16,8X,A2,1X,A2,A2,2X,A1,10X,A1,3X,F7.3,2X,F8.3,2X,F8.3,4X,
     &      '1')
 334  FORMAT(A16,8X,A2,1X,A2,A2,1X,A2,10X,A1,3X,F7.3,2X,F8.3,2X,F8.3,4X,
     &      '1')
 335  FORMAT(A16,8X,A2,1X,A2,A3,A2,10X,A1,3X,F7.3,2X,F8.3,2X,F8.3,4X,
     &      '1')
 336  FORMAT(A16,8X,A2,A3,A2,2X,A1,10X,A1,3X,F7.3,2X,F8.3,2X,F8.3,4X,
     &      '1')
 337  FORMAT(A16,8X,A2,A3,A2,1X,A2,10X,A1,3X,F7.3,2X,F8.3,2X,F8.3,4X,
     &      '1')
 338  FORMAT(A16,8X,A2,A3,A3,A2,10X,A1,3X,F7.3,2X,F8.3,2X,F8.3,4X,
     &      '1')
 340  FORMAT(A16,8X,A2,2X,A1,A2,2X,A1,A2,2X,A1,5X,A1,3X,F7.3,2X,F8.3,2X,
     &      F8.3,4X,'1')
 341  FORMAT(A16,8X,A2,2X,A1,A2,2X,A1,A2,1X,A2,5X,A1,3X,F7.3,2X,F8.3,2X,
     &      F8.3,4X,'1')
 342  FORMAT(A16,8X,A2,2X,A1,A2,2X,A1,A2,A3,5X,A1,3X,F7.3,2X,F8.3,2X,
     &      F8.3,4X,'1')
 343  FORMAT(A16,8X,A2,2X,A1,A2,1X,A2,A2,2X,A1,5X,A1,3X,F7.3,2X,F8.3,2X,
     &      F8.3,4X,'1')
 344  FORMAT(A16,8X,A2,2X,A1,A2,1X,A2,A2,1X,A2,5X,A1,3X,F7.3,2X,F8.3,2X,
     &      F8.3,4X,'1')
 345  FORMAT(A16,8X,A2,2X,A1,A2,1X,A2,A2,A3,5X,A1,3X,F7.3,2X,F8.3,2X,
     &      F8.3,4X,'1')
 346  FORMAT(A16,8X,A2,2X,A1,A2,A3,A2,2X,A1,5X,A1,3X,F7.3,2X,F8.3,2X,
     &      F8.3,4X,'1')
 347  FORMAT(A16,8X,A2,2X,A1,A2,A3,A2,1X,A2,5X,A1,3X,F7.3,2X,F8.3,2X,
     &      F8.3,4X,'1')
 348  FORMAT(A16,8X,A2,2X,A1,A2,A3,A2,A3,5X,A1,3X,F7.3,2X,F8.3,2X,
     &      F8.3,4X,'1')
 350  FORMAT(A16,8X,A2,1X,A2,A2,2X,A1,A2,2X,A1,5X,A1,3X,F7.3,2X,F8.3,2X,
     &      F8.3,4X,'1')
 351  FORMAT(A16,8X,A2,1X,A2,A2,2X,A1,A2,1X,A2,5X,A1,3X,F7.3,2X,F8.3,2X,
     &      F8.3,4X,'1')
 352  FORMAT(A16,8X,A2,1X,A2,A2,2X,A1,A2,A3,5X,A1,3X,F7.3,2X,F8.3,2X,
     &      F8.3,4X,'1')
 353  FORMAT(A16,8X,A2,1X,A2,A2,1X,A2,A2,2X,A1,5X,A1,3X,F7.3,2X,F8.3,2X,
     &      F8.3,4X,'1')
 354  FORMAT(A16,8X,A2,1X,A2,A2,1X,A2,A2,1X,A2,5X,A1,3X,F7.3,2X,F8.3,2X,
     &      F8.3,4X,'1')
 355  FORMAT(A16,8X,A2,1X,A2,A2,1X,A2,A2,A3,5X,A1,3X,F7.3,2X,F8.3,2X,
     &      F8.3,4X,'1')
 356  FORMAT(A16,8X,A2,1X,A2,A2,A3,A2,2X,A1,5X,A1,3X,F7.3,2X,F8.3,2X,
     &      F8.3,4X,'1')
 357  FORMAT(A16,8X,A2,1X,A2,A2,A3,A2,1X,A2,5X,A1,3X,F7.3,2X,F8.3,2X,
     &      F8.3,4X,'1')
 358  FORMAT(A16,8X,A2,1X,A2,A2,A3,A2,A3,5X,A1,3X,F7.3,2X,F8.3,2X,
     &      F8.3,4X,'1')
 360  FORMAT(A16,8X,A2,A3,A2,2X,A1,A2,2X,A1,5X,A1,3X,F7.3,2X,F8.3,2X,
     &      F8.3,4X,'1')
 361  FORMAT(A16,8X,A2,A3,A2,2X,A1,A2,1X,A2,5X,A1,3X,F7.3,2X,F8.3,2X,
     &      F8.3,4X,'1')
 362  FORMAT(A16,8X,A2,A3,A2,2X,A1,A2,A3,5X,A1,3X,F7.3,2X,F8.3,2X,
     &      F8.3,4X,'1')
 363  FORMAT(A16,8X,A2,A3,A2,1X,A2,A2,2X,A1,5X,A1,3X,F7.3,2X,F8.3,2X,
     &      F8.3,4X,'1')
 364  FORMAT(A16,8X,A2,A3,A2,1X,A2,A2,1X,A2,5X,A1,3X,F7.3,2X,F8.3,2X,
     &      F8.3,4X,'1')
 365  FORMAT(A16,8X,A2,A3,A2,1X,A2,A2,A3,5X,A1,3X,F7.3,2X,F8.3,2X,
     &      F8.3,4X,'1')
 366  FORMAT(A16,8X,A2,A3,A2,A3,A2,2X,A1,5X,A1,3X,F7.3,2X,F8.3,2X,
     &      F8.3,4X,'1')
 367  FORMAT(A16,8X,A2,A3,A2,A3,A2,1X,A2,5X,A1,3X,F7.3,2X,F8.3,2X,
     &      F8.3,4X,'1')
 368  FORMAT(A16,8X,A2,A3,A2,A3,A2,A3,5X,A1,3X,F7.3,2X,F8.3,2X,
     &      F8.3,4X,'1')
 370  FORMAT(A16,8X,A2,2X,A1,A2,2X,A1,A2,2X,A1,A2,2X,A1,A1,3X,F7.3,2X,
     &      F8.3,2X,F8.3,4X,'1')
 371  FORMAT(A16,8X,A2,2X,A1,A2,2X,A1,A2,2X,A1,A2,1X,A2,A1,3X,F7.3,2X,
     &      F8.3,2X,F8.3,4X,'1')
 372  FORMAT(A16,8X,A2,2X,A1,A2,2X,A1,A2,2X,A1,A2,A3,A1,3X,F7.3,2X,F8.3,
     &      2X,F8.3,4X,'1')
 373  FORMAT(A16,8X,A2,2X,A1,A2,2X,A1,A2,1X,A2,A2,2X,A1,A1,3X,F7.3,2X,
     &      F8.3,2X,F8.3,4X,'1')
 374  FORMAT(A16,8X,A2,2X,A1,A2,2X,A1,A2,1X,A2,A2,1X,A2,A1,3X,F7.3,2X,
     &      F8.3,2X,F8.3,4X,'1')
 375  FORMAT(A16,8X,A2,2X,A1,A2,2X,A1,A2,1X,A2,A2,A3,A1,3X,F7.3,2X,F8.3,
     &      2X,F8.3,4X,'1')
 376  FORMAT(A16,8X,A2,2X,A1,A2,2X,A1,A2,A3,A2,2X,A1,A1,3X,F7.3,2X,
     &      F8.3,2X,F8.3,4X,'1')
 377  FORMAT(A16,8X,A2,2X,A1,A2,2X,A1,A2,A3,A2,1X,A2,A1,3X,F7.3,2X,
     &      F8.3,2X,F8.3,4X,'1')
 378  FORMAT(A16,8X,A2,2X,A1,A2,2X,A1,A2,A3,A2,A3,A1,3X,F7.3,2X,F8.3,
     &      2X,F8.3,4X,'1')
 380  FORMAT(A16,8X,A2,2X,A1,A2,1X,A2,A2,2X,A1,A2,2X,A1,A1,3X,F7.3,2X,
     &      F8.3,2X,F8.3,4X,'1')
 381  FORMAT(A16,8X,A2,2X,A1,A2,1X,A2,A2,2X,A1,A2,1X,A2,A1,3X,F7.3,2X,
     &      F8.3,2X,F8.3,4X,'1')
 382  FORMAT(A16,8X,A2,2X,A1,A2,1X,A2,A2,2X,A1,A2,A3,A1,3X,F7.3,2X,F8.3,
     &      2X,F8.3,4X,'1')
 383  FORMAT(A16,8X,A2,2X,A1,A2,1X,A2,A2,1X,A2,A2,2X,A1,A1,3X,F7.3,2X,
     &      F8.3,2X,F8.3,4X,'1')
 384  FORMAT(A16,8X,A2,2X,A1,A2,1X,A2,A2,1X,A2,A2,1X,A2,A1,3X,F7.3,2X,
     &      F8.3,2X,F8.3,4X,'1')
 385  FORMAT(A16,8X,A2,2X,A1,A2,1X,A2,A2,1X,A2,A2,A3,A1,3X,F7.3,2X,F8.3,
     &      2X,F8.3,4X,'1')
 386  FORMAT(A16,8X,A2,2X,A1,A2,1X,A2,A2,A3,A2,2X,A1,A1,3X,F7.3,2X,
     &      F8.3,2X,F8.3,4X,'1')
 387  FORMAT(A16,8X,A2,2X,A1,A2,1X,A2,A2,A3,A2,1X,A2,A1,3X,F7.3,2X,
     &      F8.3,2X,F8.3,4X,'1')
 388  FORMAT(A16,8X,A2,2X,A1,A2,1X,A2,A2,A3,A2,A3,A1,3X,F7.3,2X,F8.3,
     &      2X,F8.3,4X,'1')
 390  FORMAT(A16,8X,A2,2X,A1,A2,A3,A2,2X,A1,A2,2X,A1,A1,3X,F7.3,2X,
     &      F8.3,2X,F8.3,4X,'1')
 391  FORMAT(A16,8X,A2,2X,A1,A2,A3,A2,2X,A1,A2,1X,A2,A1,3X,F7.3,2X,
     &      F8.3,2X,F8.3,4X,'1')
 392  FORMAT(A16,8X,A2,2X,A1,A2,A3,A2,2X,A1,A2,A3,A1,3X,F7.3,2X,F8.3,
     &      2X,F8.3,4X,'1')
 393  FORMAT(A16,8X,A2,2X,A1,A2,A3,A2,1X,A2,A2,2X,A1,A1,3X,F7.3,2X,
     &      F8.3,2X,F8.3,4X,'1')
 394  FORMAT(A16,8X,A2,2X,A1,A2,A3,A2,1X,A2,A2,1X,A2,A1,3X,F7.3,2X,
     &      F8.3,2X,F8.3,4X,'1')
 395  FORMAT(A16,8X,A2,2X,A1,A2,A3,A2,1X,A2,A2,A3,A1,3X,F7.3,2X,F8.3,
     &      2X,F8.3,4X,'1')
 396  FORMAT(A16,8X,A2,2X,A1,A2,A3,A2,A3,A2,2X,A1,A1,3X,F7.3,2X,
     &      F8.3,2X,F8.3,4X,'1')
 397  FORMAT(A16,8X,A2,2X,A1,A2,A3,A2,A3,A2,1X,A2,A1,3X,F7.3,2X,
     &      F8.3,2X,F8.3,4X,'1')
 398  FORMAT(A16,8X,A2,2X,A1,A2,A3,A2,A3,A2,A3,A1,3X,F7.3,2X,F8.3,
     &      2X,F8.3,4X,'1')
 470  FORMAT(A16,8X,A2,1X,A2,A2,2X,A1,A2,2X,A1,A2,2X,A1,A1,3X,F7.3,2X,
     &      F8.3,2X,F8.3,4X,'1')
 471  FORMAT(A16,8X,A2,1X,A2,A2,2X,A1,A2,2X,A1,A2,1X,A2,A1,3X,F7.3,2X,
     &      F8.3,2X,F8.3,4X,'1')
 472  FORMAT(A16,8X,A2,1X,A2,A2,2X,A1,A2,2X,A1,A2,A3,A1,3X,F7.3,2X,F8.3,
     &      2X,F8.3,4X,'1')
 473  FORMAT(A16,8X,A2,1X,A2,A2,2X,A1,A2,1X,A2,A2,2X,A1,A1,3X,F7.3,2X,
     &      F8.3,2X,F8.3,4X,'1')
 474  FORMAT(A16,8X,A2,1X,A2,A2,2X,A1,A2,1X,A2,A2,1X,A2,A1,3X,F7.3,2X,
     &      F8.3,2X,F8.3,4X,'1')
 475  FORMAT(A16,8X,A2,1X,A2,A2,2X,A1,A2,1X,A2,A2,A3,A1,3X,F7.3,2X,F8.3,
     &      2X,F8.3,4X,'1')
 476  FORMAT(A16,8X,A2,1X,A2,A2,2X,A1,A2,A3,A2,2X,A1,A1,3X,F7.3,2X,
     &      F8.3,2X,F8.3,4X,'1')
 477  FORMAT(A16,8X,A2,1X,A2,A2,2X,A1,A2,A3,A2,1X,A2,A1,3X,F7.3,2X,
     &      F8.3,2X,F8.3,4X,'1')
 478  FORMAT(A16,8X,A2,1X,A2,A2,2X,A1,A2,A3,A2,A3,A1,3X,F7.3,2X,F8.3,
     &      2X,F8.3,4X,'1')
 480  FORMAT(A16,8X,A2,1X,A2,A2,1X,A2,A2,2X,A1,A2,2X,A1,A1,3X,F7.3,2X,
     &      F8.3,2X,F8.3,4X,'1')
 481  FORMAT(A16,8X,A2,1X,A2,A2,1X,A2,A2,2X,A1,A2,1X,A2,A1,3X,F7.3,2X,
     &      F8.3,2X,F8.3,4X,'1')
 482  FORMAT(A16,8X,A2,1X,A2,A2,1X,A2,A2,2X,A1,A2,A3,A1,3X,F7.3,2X,F8.3,
     &      2X,F8.3,4X,'1')
 483  FORMAT(A16,8X,A2,1X,A2,A2,1X,A2,A2,1X,A2,A2,2X,A1,A1,3X,F7.3,2X,
     &      F8.3,2X,F8.3,4X,'1')
 484  FORMAT(A16,8X,A2,1X,A2,A2,1X,A2,A2,1X,A2,A2,1X,A2,A1,3X,F7.3,2X,
     &      F8.3,2X,F8.3,4X,'1')
 485  FORMAT(A16,8X,A2,1X,A2,A2,1X,A2,A2,1X,A2,A2,A3,A1,3X,F7.3,2X,F8.3,
     &      2X,F8.3,4X,'1')
 486  FORMAT(A16,8X,A2,1X,A2,A2,1X,A2,A2,A3,A2,2X,A1,A1,3X,F7.3,2X,
     &      F8.3,2X,F8.3,4X,'1')
 487  FORMAT(A16,8X,A2,1X,A2,A2,1X,A2,A2,A3,A2,1X,A2,A1,3X,F7.3,2X,
     &      F8.3,2X,F8.3,4X,'1')
 488  FORMAT(A16,8X,A2,1X,A2,A2,1X,A2,A2,A3,A2,A3,A1,3X,F7.3,2X,F8.3,
     &      2X,F8.3,4X,'1')
 490  FORMAT(A16,8X,A2,1X,A2,A2,A3,A2,2X,A1,A2,2X,A1,A1,3X,F7.3,2X,
     &      F8.3,2X,F8.3,4X,'1')
 491  FORMAT(A16,8X,A2,1X,A2,A2,A3,A2,2X,A1,A2,1X,A2,A1,3X,F7.3,2X,
     &      F8.3,2X,F8.3,4X,'1')
 492  FORMAT(A16,8X,A2,1X,A2,A2,A3,A2,2X,A1,A2,A3,A1,3X,F7.3,2X,F8.3,
     &      2X,F8.3,4X,'1')
 493  FORMAT(A16,8X,A2,1X,A2,A2,A3,A2,1X,A2,A2,2X,A1,A1,3X,F7.3,2X,
     &      F8.3,2X,F8.3,4X,'1')
 494  FORMAT(A16,8X,A2,1X,A2,A2,A3,A2,1X,A2,A2,1X,A2,A1,3X,F7.3,2X,
     &      F8.3,2X,F8.3,4X,'1')
 495  FORMAT(A16,8X,A2,1X,A2,A2,A3,A2,1X,A2,A2,A3,A1,3X,F7.3,2X,F8.3,
     &      2X,F8.3,4X,'1')
 496  FORMAT(A16,8X,A2,1X,A2,A2,A3,A2,A3,A2,2X,A1,A1,3X,F7.3,2X,
     &      F8.3,2X,F8.3,4X,'1')
 497  FORMAT(A16,8X,A2,1X,A2,A2,A3,A2,A3,A2,1X,A2,A1,3X,F7.3,2X,
     &      F8.3,2X,F8.3,4X,'1')
 498  FORMAT(A16,8X,A2,1X,A2,A2,A3,A2,A3,A2,A3,A1,3X,F7.3,2X,F8.3,
     &      2X,F8.3,4X,'1')
 570  FORMAT(A16,8X,A2,A3,A2,2X,A1,A2,2X,A1,A2,2X,A1,A1,3X,F7.3,2X,
     &      F8.3,2X,F8.3,4X,'1')
 571  FORMAT(A16,8X,A2,A3,A2,2X,A1,A2,2X,A1,A2,1X,A2,A1,3X,F7.3,2X,
     &      F8.3,2X,F8.3,4X,'1')
 572  FORMAT(A16,8X,A2,A3,A2,2X,A1,A2,2X,A1,A2,A3,A1,3X,F7.3,2X,F8.3,
     &      2X,F8.3,4X,'1')
 573  FORMAT(A16,8X,A2,A3,A2,2X,A1,A2,1X,A2,A2,2X,A1,A1,3X,F7.3,2X,
     &      F8.3,2X,F8.3,4X,'1')
 574  FORMAT(A16,8X,A2,A3,A2,2X,A1,A2,1X,A2,A2,1X,A2,A1,3X,F7.3,2X,
     &      F8.3,2X,F8.3,4X,'1')
 575  FORMAT(A16,8X,A2,A3,A2,2X,A1,A2,1X,A2,A2,A3,A1,3X,F7.3,2X,F8.3,
     &      2X,F8.3,4X,'1')
 576  FORMAT(A16,8X,A2,A3,A2,2X,A1,A2,A3,A2,2X,A1,A1,3X,F7.3,2X,
     &      F8.3,2X,F8.3,4X,'1')
 577  FORMAT(A16,8X,A2,A3,A2,2X,A1,A2,A3,A2,1X,A2,A1,3X,F7.3,2X,
     &      F8.3,2X,F8.3,4X,'1')
 578  FORMAT(A16,8X,A2,A3,A2,2X,A1,A2,A3,A2,A3,A1,3X,F7.3,2X,F8.3,
     &      2X,F8.3,4X,'1')
 580  FORMAT(A16,8X,A2,A3,A2,1X,A2,A2,2X,A1,A2,2X,A1,A1,3X,F7.3,2X,
     &      F8.3,2X,F8.3,4X,'1')
 581  FORMAT(A16,8X,A2,A3,A2,1X,A2,A2,2X,A1,A2,1X,A2,A1,3X,F7.3,2X,
     &      F8.3,2X,F8.3,4X,'1')
 582  FORMAT(A16,8X,A2,A3,A2,1X,A2,A2,2X,A1,A2,A3,A1,3X,F7.3,2X,F8.3,
     &      2X,F8.3,4X,'1')
 583  FORMAT(A16,8X,A2,A3,A2,1X,A2,A2,1X,A2,A2,2X,A1,A1,3X,F7.3,2X,
     &      F8.3,2X,F8.3,4X,'1')
 584  FORMAT(A16,8X,A2,A3,A2,1X,A2,A2,1X,A2,A2,1X,A2,A1,3X,F7.3,2X,
     &      F8.3,2X,F8.3,4X,'1')
 585  FORMAT(A16,8X,A2,A3,A2,1X,A2,A2,1X,A2,A2,A3,A1,3X,F7.3,2X,F8.3,
     &      2X,F8.3,4X,'1')
 586  FORMAT(A16,8X,A2,A3,A2,1X,A2,A2,A3,A2,2X,A1,A1,3X,F7.3,2X,
     &      F8.3,2X,F8.3,4X,'1')
 587  FORMAT(A16,8X,A2,A3,A2,1X,A2,A2,A3,A2,1X,A2,A1,3X,F7.3,2X,
     &      F8.3,2X,F8.3,4X,'1')
 588  FORMAT(A16,8X,A2,A3,A2,1X,A2,A2,A3,A2,A3,A1,3X,F7.3,2X,F8.3,
     &      2X,F8.3,4X,'1')
 590  FORMAT(A16,8X,A2,A3,A2,A3,A2,2X,A1,A2,2X,A1,A1,3X,F7.3,2X,
     &      F8.3,2X,F8.3,4X,'1')
 591  FORMAT(A16,8X,A2,A3,A2,A3,A2,2X,A1,A2,1X,A2,A1,3X,F7.3,2X,
     &      F8.3,2X,F8.3,4X,'1')
 592  FORMAT(A16,8X,A2,A3,A2,A3,A2,2X,A1,A2,A3,A1,3X,F7.3,2X,F8.3,
     &      2X,F8.3,4X,'1')
 593  FORMAT(A16,8X,A2,A3,A2,A3,A2,1X,A2,A2,2X,A1,A1,3X,F7.3,2X,
     &      F8.3,2X,F8.3,4X,'1')
 594  FORMAT(A16,8X,A2,A3,A2,A3,A2,1X,A2,A2,1X,A2,A1,3X,F7.3,2X,
     &      F8.3,2X,F8.3,4X,'1')
 595  FORMAT(A16,8X,A2,A3,A2,A3,A2,1X,A2,A2,A3,A1,3X,F7.3,2X,F8.3,
     &      2X,F8.3,4X,'1')
 596  FORMAT(A16,8X,A2,A3,A2,A3,A2,A3,A2,2X,A1,A1,3X,F7.3,2X,
     &      F8.3,2X,F8.3,4X,'1')
 597  FORMAT(A16,8X,A2,A3,A2,A3,A2,A3,A2,1X,A2,A1,3X,F7.3,2X,
     &      F8.3,2X,F8.3,4X,'1')
 598  FORMAT(A16,8X,A2,A3,A2,A3,A2,A3,A2,A3,A1,3X,F7.3,2X,F8.3,
     &      2X,F8.3,4X,'1')
 601  FORMAT(A16,28X,'G',3X,F7.3,2X,F8.3,2X,F8.3,4X,'1&')
C'
 602  FORMAT(A2,1X,A3,A2,1X,A3,A2,1X,A3,A2,1X,A3,A2,1X,A3)
 300  FORMAT(A16,8X,A4,A1,A4,A1,10X,A1,3X,F7.3,2X,F8.3,2X,F8.3,4X,'1')
 301  FORMAT(ES15.8,ES15.8,ES15.8,ES15.8,ES15.8,4X,'2')
 302  FORMAT(ES15.8,ES15.8,ES15.8,ES15.8,ES15.8,4X,'3')
 303  FORMAT(ES15.8,ES15.8,ES15.8,ES15.8,15X,4X,'4')
C
      ELSEIF (DATATYPE .EQ. 'NIST') THEN
C      write(*,*) 'THE COEFFICIENTS OF NIST FORM'
C      write(*,*) THERM1(1:7)
C      write(*,*) THERM2(1:7)
C
       WRITE(LOUT,*) 'The NIST polynomial coefficients calculated:'
       WRITE(LOUT,*) 'H0-H0_298=A*t+B*t*t/2+C*t*t*t/3+D*t*t*t*t/4'
       WRITE(LOUT,*) '         -E/t+F-H'
       WRITE(LOUT,*) 'H=DLTH_formation=',DLTH
       WRITE(LOUT,300) SNAM,ENAM(1),ELNO(1),ENAM(2),ELNO(2),'G',
     &                   T_min,T_max,Tint
       WRITE(LOUT,301) THERM1(1),THERM1(2),THERM1(3),THERM1(4),THERM1(5)
       WRITE(LOUT,302) THERM1(6),THERM1(7),THERM2(1),THERM2(2),THERM2(3)
       WRITE(LOUT,303) THERM2(4),THERM2(5),THERM2(6),THERM2(7)
C
      ELSE IF (DATATYPE .EQ. 'NASA') THEN
C      write(*,*) 'THE COEFFICIENTS OF NASA FORM'
C      write(*,*) THERM1(1:9)
C      write(*,*) THERM2(1:9)
C      write(*,*) THERM3(1:9)
C
       WRITE(LOUT,*) 'The NASA polynomial coefficients calculated:'
       WRITE(LOUT,399) SNAM
C HF_298 in 'J/mol' in NASA format: 1 cal = 4.184 J
C      WRITE(*,*) MW
C      WRITE(*,*) DLTH*4.184D0
       WRITE(LOUT,400) 3,ENAM(1),ELNO(1),'.00',ENAM(2),ELNO(2),'.00',
     &                   ENAM(3),ELNO(3),'.00',ENAM(4),ELNO(4),'.00',
     &                 ENAM(5),ELNO(5),'.00',MW,DLTH*4.184D0
       WRITE(LOUT,401) T_min,Tint,7,'-2.0','-1.0','0.0','1.0','2.0',
     &                   '3.0','4.0','0.0',DLTH298(THERM1(1:7))*4.184D0
       WRITE(LOUT,404) THERM1(1),THERM1(2),THERM1(3),THERM1(4),THERM1(5)
       WRITE(LOUT,404) THERM1(6),THERM1(7),0.D0,THERM1(8),THERM1(9)
       WRITE(LOUT,402) Tint,Tint1,7,'-2.0','-1.0','0.0','1.0','2.0',
     &                   '3.0','4.0','0.0',DLTH298(THERM1(1:7))*4.184D0
       WRITE(LOUT,404) THERM2(1),THERM2(2),THERM2(3),THERM2(4),THERM2(5)
       WRITE(LOUT,404) THERM2(6),THERM2(7),0.D0,THERM2(8),THERM2(9)
       WRITE(LOUT,403) Tint1,T_max,7,'-2.0','-1.0','0.0','1.0','2.0',
     &                   '3.0','4.0','0.0',DLTH298(THERM1(1:7))*4.184D0
       WRITE(LOUT,404) THERM3(1),THERM3(2),THERM3(3),THERM3(4),THERM3(5)
       WRITE(LOUT,404) THERM3(6),THERM3(7),0.D0,THERM3(8),THERM3(9)
C
      ENDIF
C
 399  FORMAT(A16)
 400  FORMAT(1X,I1,8X,A4,A1,A3,A4,A1,A3,A4,A1,A3,A4,A1,A3,A4,A1,A3,
     &       1X,'0',2X,F11.7,4X,F11.3)
 401  FORMAT(4X,F7.3,3X,F7.3,1X,I1,1X,A4,1X,A4,2X,A3,2X,A3,2X,A3,2X,A3,
     &       2X,A3,2X,A3,6X,F11.3)
 402  FORMAT(4X,F7.3,2X,F8.3,1X,I1,1X,A4,1X,A4,2X,A3,2X,A3,2X,A3,2X,A3,
     &       2X,A3,2X,A3,6X,F11.3)
 403  FORMAT(3X,F8.3,2X,F8.3,1X,I1,1X,A4,1X,A4,2X,A3,2X,A3,2X,A3,2X,A3,
     &       2X,A3,2X,A3,6X,F11.3)
 404  FORMAT(5ES16.9)
C
      RETURN
C
      END
C***********************************************************************
C
      SUBROUTINE DATAGROUP(M_in,TEMP,CPT,CH,CS,LIN)
C
      INTEGER M_in,I,LIN,IJ,IK,IM
      CHARACTER(LEN=4) MARK
      CHARACTER(LEN=40) TEXT
      DOUBLE PRECISION TEMP,CPT,CH,CS
      DIMENSION TEMP(*),CPT(*),CH(*),CS(*)
      DOUBLE PRECISION TEMP1,C1,CH1,CS1
C     DIMENSION TEMP1(50),C1(50),CH1(50),CS1(50)
C
C     DATA TEMP1/290.262D0, 393.258D0, 496.255D0, 589.888D0, 795.88D0,
C     &           1001.87D0, 1488.76D0/
C
C     DATA  C1/17.7725D0, 22.4645D0, 27.0142D0, 30.7109D0, 36.9668D0,
C     &            41.8009D0, 49.3365D0/
C
C     DATA CH1/0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0/
C
C     DATA CS1/0.D0,0.D0,0.D0,0.D0,0.D0,0.D0,0.D0/
C
      MARK=''
      TEXT=''
 100  FORMAT (A4,1X,A40)
 101  FORMAT (A4,1X,F4.0,1X,F8.1)
C
C
      DO I=1,M_in
         READ(LIN, *, END=9998) MARK, TEMP(I), CPT(I)
         CH(I)  = 0.0D0
         CS(I)  = 0.0D0
         IF (MARK .NE. 'TECP') THEN
              write(*,*) 'Error. Need more TECP'
            GOTO 250
         ENDIF
      ENDDO
C
 250  CONTINUE
      RETURN
 9998 CONTINUE
      WRITE(*,*) 'Error. Input ended sooner than expected!'
      STOP
C
      END
