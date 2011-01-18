C: Add the constrains to the first and second derivatives
C: Therefore, there are three constrains: Cp, DCp, and D2Cp
C: at each intermediate temperature
C: --------------------------------------------------------------
C: For NASA
C: Use a single matrix to calculate all 21 coefficients
C: [A]: M x 21
C: [X]: 21 x 1
C:     [alpha1_l,alpha2_l,alpha3_l,alpha4_l,alpha5_l,alpha6_l,alpha7_l,
C:     [alpha1_m,alpha2_m,alpha3_m,alpha4_m,alpha5_m,alpha6_m,alpha7_m,
C:     [alpha1_h,alpha2_h,alpha3_h,alpha4_h,alpha5_h,alpha6_h,alpha7_h]
C: [C]: M x 1
C: Contrains at both temperatures:
C: [B]: 6 x 21
C: [D]: 6 x 1
C----------------------------------------------------------------
C FOR NASA CHEMTHERMO DATA FORMAT
C number of coefficients: N=7 x 3
C number of discrete data point: M>=N, i.e., M=103
C
C use the DGGLSE to solve the linear least squares equality-constrained
C (LLSE) problem.
!****************************************************************
! Programmed by John Z. Wen, MIT, June 2005
!****************************************************************
C
      SUBROUTINE NASA(Tint,TNEW,Tint1,T1NEW,
     &      THERM1,THERM2,THERM3,TEMP1,C1,CH1,CS1)
C
      implicit none
C
      DOUBLE PRECISION Tint,Tint1,THERM1,THERM2,THERM3,TNEW,T1NEW
      DOUBLE PRECISION H_298,S_298,DLTH,ROL
C
      INTEGER M,N,P,LDA,LDB,LWORK,INFO
      INTEGER MH,NH,PH,LDAH,LDBH,LWORKH
      INTEGER MS,NS,PS,LDAS,LDBS,LWORKS
      DOUBLE PRECISION A1,B1,Bk,C1,D1,X1,WORK1
      DOUBLE PRECISION TEMP1,PHI_NASA1,PHI_NASA2,PHI_NASA3,
     &                   PHI_NASA4,PHI_NASA5,PHI_NASA6,PHI_NASA7
      DOUBLE PRECISION PHI_NASA_DCP1,PHI_NASA_DCP2,PHI_NASA_DCP3,
     &       PHI_NASA_DCP4,PHI_NASA_DCP5,PHI_NASA_DCP6,PHI_NASA_DCP7,
     &       PHI_NASA_D2CP1,PHI_NASA_D2CP2,PHI_NASA_D2CP3,
     &       PHI_NASA_D2CP4,PHI_NASA_D2CP5,PHI_NASA_D2CP6,
     &       PHI_NASA_D2CP7
      DOUBLE PRECISION A_H1,A_S1,A_H2,A_S2,A_H3,A_S3
      DOUBLE PRECISION AH1,BH1,CH1,DH1,XH1,WORKH,BHk,CH,CH2
      DOUBLE PRECISION AS1,BS1,CS1,DS1,XS1,WORKS,BSk,CS,CS2
      DOUBLE PRECISION PHI_NASA_H,PHI_NASA_S
      INTEGER MARK1,MARK2
      DOUBLE PRECISION TEMP,C
C: Required by llse solver:
! LDA >= MAX(1,M)
! LDB >= MAX(1,P)
! LWORK >= max(1,M+N+P)
! P=1 for only one constrain at the intermediate temp.
! P=3 for three constrains at the intermediate temp.
      PARAMETER (M=103)
      PARAMETER (LDB=6,P=6)
      PARAMETER (N=21)
      DIMENSION A1(M,N),B1(LDB,N),C1(M-2),D1(P),X1(N),WORK1(M+N+P),
     &            Bk(LDB,N)
      PARAMETER (LDBH=3,PH=3)
      PARAMETER (NH=3)
      PARAMETER (LDBS=3,PS=3)
      PARAMETER (NS=3)
      DIMENSION AH1(M,NH),BH1(LDBH,NH),CH1(M-2),DH1(PH),XH1(NH),
     &            WORKH(M+NH+PH),BHk(LDBH,NH),CH(M),CH2(M)
      DIMENSION AS1(M,NS),BS1(LDBS,NS),CS1(M-2),DS1(PS),XS1(NS),
     &            WORKS(M+NS+PS),BSk(LDBS,NS),CS(M),CS2(M)
      DIMENSION TEMP1(M-2),THERM1(20),THERM2(20),THERM3(20)
      DIMENSION TEMP(M),C(M)
      INTEGER I,J
C
! GAS CONSTANT
      COMMON /GAS/ROL,H_298,S_298,DLTH
C
      EXTERNAL PHI_NASA1,PHI_NASA2,PHI_NASA3,PHI_NASA4,PHI_NASA5,
     &         PHI_NASA6,PHI_NASA7,
     &         PHI_NASA_DCP1,PHI_NASA_DCP2,PHI_NASA_DCP3,
     &         PHI_NASA_DCP4,PHI_NASA_DCP5,PHI_NASA_DCP6,
     &         PHI_NASA_DCP7,
     &         PHI_NASA_D2CP1,PHI_NASA_D2CP2,PHI_NASA_D2CP3,
     &         PHI_NASA_D2CP4,PHI_NASA_D2CP5,PHI_NASA_D2CP6,
     &         PHI_NASA_D2CP7,
     &         PHI_NASA_H,PHI_NASA_S
C
!     WRITE(*,*) 'CALLED NASA'
C
      LDA   = M
      LWORK = M+N+P
      MH    = M
      LDAH  = MH
      LWORKH= MH+NH+PH
      MS    = M
      LDAS  = MS
      LWORKS= MS+NS+PS
C
C FIND 'Tint' and 'Tint1' in sample data
C
      MARK1 = 0
      DO I=1,M-2
        IF (TEMP1(I) .GT. Tint) THEN
           MARK1=I-1
           IF (Tint .NE. TEMP1(I-1)) THEN
                Tint=TEMP1(I-1)
              TNEW=Tint
!             WRITE(*,*) 'Tint has been replaced!'
           ENDIF
           GOTO 995
        ENDIF
      ENDDO
 995  CONTINUE
!     WRITE(*,*) 'MARK1=',MARK1
C
      MARK2 = 0
      DO I=1,M-2
        IF (TEMP1(I) .GT. Tint1) THEN
           MARK2=I-1
           IF (Tint1 .NE. TEMP1(I-1)) THEN
                Tint1=TEMP1(I-1)
              T1NEW=Tint1
!             WRITE(*,*) 'Tint1 has been replaced!'
           ENDIF
           GOTO 996
        ENDIF
      ENDDO
 996  CONTINUE
!     WRITE(*,*) 'MARK2=',MARK2
C
C Rearrange the matrics
C
      DO I=1,M
         TEMP(I)=0.D0
         C   (I)=0.D0
         CH  (I)=0.D0
         CS  (I)=0.D0
         CH2 (I)=0.D0
         CS2 (I)=0.D0
      ENDDO
C
      DO I=1,MARK1
         TEMP(I)=TEMP1(I)
         C   (I)=C1   (I)
         CH  (I)=CH1  (I)
         CS  (I)=CS1  (I)
      ENDDO
      TEMP(MARK1+1)=TEMP(MARK1)
      C   (MARK1+1)=C   (MARK1)
      CH  (MARK1+1)=CH  (MARK1)
      CS  (MARK1+1)=CS  (MARK1)
      DO I=MARK1+2,MARK2+1
         TEMP(I)=TEMP1(I-1)
         C   (I)=C1   (I-1)
         CH  (I)=CH1  (I-1)
         CS  (I)=CS1  (I-1)
      ENDDO
      TEMP(MARK2+2)=TEMP(MARK2+1)
      C   (MARK2+2)=C   (MARK2+1)
      CH  (MARK2+2)=CH  (MARK2+1)
      CS  (MARK2+2)=CS  (MARK2+1)
      DO I=MARK2+3,M
         TEMP(I)=TEMP1(I-2)
         C   (I)=C1   (I-2)
         CH  (I)=CH1  (I-2)
         CS  (I)=CS1  (I-2)
      ENDDO
!     WRITE(*,*) 'Tint =', Tint
!     WRITE(*,*) 'Tint1=', Tint1
!     WRITE(*,*) TEMP
!     WRITE(*,*) C
!     WRITE(*,*) CH
!     WRITE(*,*) CS
C
C: Redefine the MARKS
C
      DO I=1,M
         IF (TEMP(I) .EQ. Tint ) THEN
             MARK1=I
             GOTO 997
         ENDIF
      ENDDO
 997  CONTINUE
      DO I=1,M
         IF (TEMP(I) .EQ. Tint1) THEN
             MARK2=I
             GOTO 998
         ENDIF
      ENDDO
 998  CONTINUE
!     WRITE(*,*) 'MARK1 =', MARK1
!     WRITE(*,*) 'MARK2 =', MARK2
!     PAUSE
C
      DO I=1,N
         X1(I)=0.D0
      ENDDO
      DO I=1,N
         DO J=1,LDA
           A1(J,I)=0.D0
         ENDDO
         DO J=1,LDB
           B1(J,I)=0.D0
           Bk(J,I)=0.D0
         ENDDO
      ENDDO
      DO I=1,LWORK
         WORK1(I)=0.D0
      ENDDO
C
C For PHI matrix in FIRST temp range
C
        DO I=1,MARK1
           A1(I,1)=PHI_NASA1(TEMP(I))
           A1(I,2)=PHI_NASA2(TEMP(I))
           A1(I,3)=PHI_NASA3(TEMP(I))
           A1(I,4)=PHI_NASA4(TEMP(I))
           A1(I,5)=PHI_NASA5(TEMP(I))
           A1(I,6)=PHI_NASA6(TEMP(I))
           A1(I,7)=PHI_NASA7(TEMP(I))
        ENDDO
C
C For PHI matrix in SECOND temp range
C
        DO I=MARK1+1,MARK2
           A1(I,8) =PHI_NASA1(TEMP(I))
           A1(I,9) =PHI_NASA2(TEMP(I))
           A1(I,10)=PHI_NASA3(TEMP(I))
           A1(I,11)=PHI_NASA4(TEMP(I))
           A1(I,12)=PHI_NASA5(TEMP(I))
           A1(I,13)=PHI_NASA6(TEMP(I))
           A1(I,14)=PHI_NASA7(TEMP(I))
        ENDDO
C
C For PHI matrix in THIRD temp range
C
        DO I=MARK2+1,M
           A1(I,15)=PHI_NASA1(TEMP(I))
           A1(I,16)=PHI_NASA2(TEMP(I))
           A1(I,17)=PHI_NASA3(TEMP(I))
           A1(I,18)=PHI_NASA4(TEMP(I))
           A1(I,19)=PHI_NASA5(TEMP(I))
           A1(I,20)=PHI_NASA6(TEMP(I))
           A1(I,21)=PHI_NASA7(TEMP(I))
        ENDDO
C
C For Constrains at Tint-------------------------
C
C For CP values
        Bk(1,1) = PHI_NASA1(Tint)
        Bk(1,2) = PHI_NASA2(Tint)
        Bk(1,3) = PHI_NASA3(Tint)
        Bk(1,4) = PHI_NASA4(Tint)
        Bk(1,5) = PHI_NASA5(Tint)
        Bk(1,6) = PHI_NASA6(Tint)
        Bk(1,7) = PHI_NASA7(Tint)
        Bk(1,8) =-PHI_NASA1(Tint)
        Bk(1,9) =-PHI_NASA2(Tint)
        Bk(1,10)=-PHI_NASA3(Tint)
        Bk(1,11)=-PHI_NASA4(Tint)
        Bk(1,12)=-PHI_NASA5(Tint)
        Bk(1,13)=-PHI_NASA6(Tint)
        Bk(1,14)=-PHI_NASA7(Tint)
        D1(1)=0.D0
C For DCp values
        Bk(2,1) = PHI_NASA_DCP1(Tint)
        Bk(2,2) = PHI_NASA_DCP2(Tint)
        Bk(2,3) = PHI_NASA_DCP3(Tint)
        Bk(2,4) = PHI_NASA_DCP4(Tint)
        Bk(2,5) = PHI_NASA_DCP5(Tint)
        Bk(2,6) = PHI_NASA_DCP6(Tint)
        Bk(2,7) = PHI_NASA_DCP7(Tint)
        Bk(2,8) =-PHI_NASA_DCP1(Tint)
        Bk(2,9) =-PHI_NASA_DCP2(Tint)
        Bk(2,10)=-PHI_NASA_DCP3(Tint)
        Bk(2,11)=-PHI_NASA_DCP4(Tint)
        Bk(2,12)=-PHI_NASA_DCP5(Tint)
        Bk(2,13)=-PHI_NASA_DCP6(Tint)
        Bk(2,14)=-PHI_NASA_DCP7(Tint)
        D1(2)=0.D0
C For D2Cp values
        Bk(3,1) = PHI_NASA_D2CP1(Tint)
        Bk(3,2) = PHI_NASA_D2CP2(Tint)
        Bk(3,3) = PHI_NASA_D2CP3(Tint)
        Bk(3,4) = PHI_NASA_D2CP4(Tint)
        Bk(3,5) = PHI_NASA_D2CP5(Tint)
        Bk(3,6) = PHI_NASA_D2CP6(Tint)
        Bk(3,7) = PHI_NASA_D2CP7(Tint)
        Bk(3,8) =-PHI_NASA_D2CP1(Tint)
        Bk(3,9) =-PHI_NASA_D2CP2(Tint)
        Bk(3,10)=-PHI_NASA_D2CP3(Tint)
        Bk(3,11)=-PHI_NASA_D2CP4(Tint)
        Bk(3,12)=-PHI_NASA_D2CP5(Tint)
        Bk(3,13)=-PHI_NASA_D2CP6(Tint)
        Bk(3,14)=-PHI_NASA_D2CP7(Tint)
        D1(3)=0.D0
C
C For Constrains at Tint1----------------------------
C
C For CP values
        Bk(4,8) = PHI_NASA1(Tint1)
        Bk(4,9) = PHI_NASA2(Tint1)
        Bk(4,10)= PHI_NASA3(Tint1)
        Bk(4,11)= PHI_NASA4(Tint1)
        Bk(4,12)= PHI_NASA5(Tint1)
        Bk(4,13)= PHI_NASA6(Tint1)
        Bk(4,14)= PHI_NASA7(Tint1)
        Bk(4,15)=-PHI_NASA1(Tint1)
        Bk(4,16)=-PHI_NASA2(Tint1)
        Bk(4,17)=-PHI_NASA3(Tint1)
        Bk(4,18)=-PHI_NASA4(Tint1)
        Bk(4,19)=-PHI_NASA5(Tint1)
        Bk(4,20)=-PHI_NASA6(Tint1)
        Bk(4,21)=-PHI_NASA7(Tint1)
        D1(4)=0.D0
C For DCp values
        Bk(5,8) = PHI_NASA_DCP1(Tint1)
        Bk(5,9) = PHI_NASA_DCP2(Tint1)
        Bk(5,10)= PHI_NASA_DCP3(Tint1)
        Bk(5,11)= PHI_NASA_DCP4(Tint1)
        Bk(5,12)= PHI_NASA_DCP5(Tint1)
        Bk(5,13)= PHI_NASA_DCP6(Tint1)
        Bk(5,14)= PHI_NASA_DCP7(Tint1)
        Bk(5,15)=-PHI_NASA_DCP1(Tint1)
        Bk(5,16)=-PHI_NASA_DCP2(Tint1)
        Bk(5,17)=-PHI_NASA_DCP3(Tint1)
        Bk(5,18)=-PHI_NASA_DCP4(Tint1)
        Bk(5,19)=-PHI_NASA_DCP5(Tint1)
        Bk(5,20)=-PHI_NASA_DCP6(Tint1)
        Bk(5,21)=-PHI_NASA_DCP7(Tint1)
        D1(5)=0.D0
C For D2Cp values
        Bk(6,8) = PHI_NASA_D2CP1(Tint1)
        Bk(6,9) = PHI_NASA_D2CP2(Tint1)
        Bk(6,10)= PHI_NASA_D2CP3(Tint1)
        Bk(6,11)= PHI_NASA_D2CP4(Tint1)
        Bk(6,12)= PHI_NASA_D2CP5(Tint1)
        Bk(6,13)= PHI_NASA_D2CP6(Tint1)
        Bk(6,14)= PHI_NASA_D2CP7(Tint1)
        Bk(6,15)=-PHI_NASA_D2CP1(Tint1)
        Bk(6,16)=-PHI_NASA_D2CP2(Tint1)
        Bk(6,17)=-PHI_NASA_D2CP3(Tint1)
        Bk(6,18)=-PHI_NASA_D2CP4(Tint1)
        Bk(6,19)=-PHI_NASA_D2CP5(Tint1)
        Bk(6,20)=-PHI_NASA_D2CP6(Tint1)
        Bk(6,21)=-PHI_NASA_D2CP7(Tint1)
        D1(6)=0.D0
C
      DO I=1,N
         DO J=1,LDB
            B1(J,I)=Bk(J,I)
         ENDDO
      ENDDO
C
! NASA FORMAT
C
      call DGGLSE(M,N,P,A1,LDA,B1,LDB,C,D1,X1,WORK1,LWORK,INFO)
!     WRITE(*,*) 'X1=',X1(1:21)
!     pause
C
C For enthalpy, NEED define 2 coeffs, N=3
      DO I=1,NH
         XH1(I)=0.D0
      ENDDO
      DO I=1,NH
         DO J=1,LDAH
           AH1(J,I)=0.D0
         ENDDO
         DO J=1,LDBH
           BH1(J,I)=0.D0
           BHk(J,I)=0.D0
         ENDDO
      ENDDO
      DO I=1,LWORKH
         WORKH(I)=0.D0
      ENDDO
C
      DO I=1,MARK1
         CH2(I)=CH(I) - PHI_NASA_H(TEMP(I),X1(1:7))
      ENDDO
C
      DO I=MARK1+1,MARK2
         CH2(I)=CH(I) - PHI_NASA_H(TEMP(I),X1(8:14))
      ENDDO
C
      DO I=MARK2+1,MH
         CH2(I)=CH(I) - PHI_NASA_H(TEMP(I),X1(15:21))
      ENDDO
C
      DO I=1,MARK1
         AH1(I,1)=ROL
      ENDDO
      DO I=MARK1+1,MARK2
         AH1(I,2)=ROL
      ENDDO
      DO I=MARK2+1,MH
         AH1(I,3)=ROL
      ENDDO
C
C Contrain for H at Tint
      BHk(1,1)=ROL
      BHk(1,2)=-ROL
      DH1(1)=PHI_NASA_H(Tint,X1(8:14))
     &      -PHI_NASA_H(Tint,X1(1:7))
C
C Contrain for H at 298.15K
      BHk(2,1)=ROL
      BHk(2,2)=0.D0
!     DH1(2)=CH(1) - PHI_NASA_H(298.D0,X1(1:7))
      DH1(2)=DLTH - PHI_NASA_H(298.15D0,X1(1:7))
C
C Contrain for H at Tint1
      BHk(3,2)=ROL
      BHk(3,3)=-ROL
      DH1(3)=PHI_NASA_H(Tint1,X1(15:21))
     &      -PHI_NASA_H(Tint1,X1(8:14))
C
      DO I=1,NH
         DO J=1,LDBH
            BH1(J,I)=BHk(J,I)
         ENDDO
      ENDDO
C
      call DGGLSE(MH,NH,PH,AH1,LDAH,BH1,LDBH,CH2,DH1,XH1,WORKH,
     &              LWORKH,INFO)
!     WRITE(*,*) 'H_CONST1=',XH1(1)
!     WRITE(*,*) 'H_CONST2=',XH1(2)
!     WRITE(*,*) 'H_CONST3=',XH1(3)
C
      A_H1=XH1(1)
      A_H2=XH1(2)
      A_H3=XH1(3)
C
C For entropy, NEED define 2 coeffs, N=2
      DO I=1,NS
         XS1(I)=0.D0
      ENDDO
      DO I=1,NS
         DO J=1,LDAS
           AS1(J,I)=0.D0
         ENDDO
         DO J=1,LDBS
           BS1(J,I)=0.D0
           BSk(J,I)=0.D0
         ENDDO
      ENDDO
      DO I=1,LWORKS
         WORKS(I)=0.D0
      ENDDO
C
      DO I=1,MARK1
         CS2(I)=CS(I) - PHI_NASA_S(TEMP(I),X1(1:7))
      ENDDO
      DO I=MARK1+1,MARK2
         CS2(I)=CS(I) - PHI_NASA_S(TEMP(I),X1(8:14))
      ENDDO
      DO I=MARK2+1,MS
         CS2(I)=CS(I) - PHI_NASA_S(TEMP(I),X1(15:21))
      ENDDO
C
      DO I=1,MARK1
         AS1(I,1)=ROL
      ENDDO
      DO I=MARK1+1,MARK2
         AS1(I,2)=ROL
      ENDDO
      DO I=MARK2+1,MS
         AS1(I,3)=ROL
      ENDDO
C
C Contrain for S at Tint
      BSk(1,1)=ROL
      BSk(1,2)=-ROL
      DS1(1)=PHI_NASA_S(Tint,X1(8:14))
     &      -PHI_NASA_S(Tint,X1(1:7))
C
C Contrain for S at 298.15K
      BSk(2,1)=ROL
      BSk(2,2)=0.D0
!     DS1(2)=CS(1) - PHI_NASA_S(298.D0,X1(1:7))
      DS1(2)=S_298 - PHI_NASA_S(298.15D0,X1(1:7))
C
C Contrain for S at Tint1
      BSk(3,2)=ROL
      BSk(3,3)=-ROL
      DS1(3)=PHI_NASA_S(Tint1,X1(15:21))
     &      -PHI_NASA_S(Tint1,X1(8:14))
C
      DO I=1,NS
         DO J=1,LDBS
            BS1(J,I)=BSk(J,I)
         ENDDO
      ENDDO
C
      call DGGLSE(MS,NS,PS,AS1,LDAS,BS1,LDBS,CS2,DS1,XS1,WORKS,
     &              LWORKS,INFO)
!     WRITE(*,*) XS1
!     WRITE(*,*) 'S_CONST1=',XS1(1)
!     WRITE(*,*) 'S_CONST2=',XS1(2)
!     WRITE(*,*) 'S_CONST3=',XS1(3)
C
      A_S1=XS1(1)
      A_S2=XS1(2)
      A_S3=XS1(3)
C
C
C
        DO I=1,7
           THERM1(I)=X1(I)
        ENDDO
        THERM1(7+1)=A_H1
        THERM1(7+2)=A_S1
        DO I=1,7
           THERM2(I)=X1(I+7)
        ENDDO
        THERM2(7+1)=A_H2
        THERM2(7+2)=A_S2
        DO I=1,7
           THERM3(I)=X1(I+7*2)
        ENDDO
        THERM3(7+1)=A_H3
        THERM3(7+2)=A_S3
!     write(*,*) therm1(1:9)
!     write(*,*) therm2(1:9)
!     write(*,*) therm3(1:9)
C
!     write(*,*) 'NASA OVER!'
C
      RETURN
C
      END
C
C
      DOUBLE PRECISION FUNCTION PHI_NASA1(TEM)
C
      implicit none
C
      DOUBLE PRECISION TEM,ROL,H_298,S_298,DLTH
      COMMON /GAS/ROL,H_298,S_298,DLTH
C
      PHI_NASA1=ROL/(TEM*TEM)
C
      RETURN
C
      END
C
      DOUBLE PRECISION FUNCTION PHI_NASA2(TEM)
C
      implicit none
C
      DOUBLE PRECISION TEM,ROL,H_298,S_298,DLTH
      COMMON /GAS/ROL,H_298,S_298,DLTH
C
      PHI_NASA2=ROL/TEM
C
      RETURN
C
      END
C
      DOUBLE PRECISION FUNCTION PHI_NASA3(TEM)
C
      implicit none
C
      DOUBLE PRECISION TEM,ROL,H_298,S_298,DLTH
      COMMON /GAS/ROL,H_298,S_298,DLTH
C
      PHI_NASA3=ROL
C
      RETURN
C
      END
C
C
      DOUBLE PRECISION FUNCTION PHI_NASA4(TEM)
C
      implicit none
C
      DOUBLE PRECISION TEM,ROL,H_298,S_298,DLTH
      COMMON /GAS/ROL,H_298,S_298,DLTH
C
      PHI_NASA4=ROL*TEM
C
      RETURN
C
      END
C
      DOUBLE PRECISION FUNCTION PHI_NASA5(TEM)
C
      implicit none
C
      DOUBLE PRECISION TEM,ROL,H_298,S_298,DLTH
      COMMON /GAS/ROL,H_298,S_298,DLTH
C
      PHI_NASA5=ROL*TEM*TEM
C
      RETURN
C
      END
C
      DOUBLE PRECISION FUNCTION PHI_NASA6(TEM)
C
      implicit none
C
      DOUBLE PRECISION TEM,ROL,H_298,S_298,DLTH
      COMMON /GAS/ROL,H_298,S_298,DLTH
C
      PHI_NASA6=ROL*TEM*TEM*TEM
C
      RETURN
C
      END
C
      DOUBLE PRECISION FUNCTION PHI_NASA7(TEM)
C
      implicit none
C
      DOUBLE PRECISION TEM,ROL,H_298,S_298,DLTH
      COMMON /GAS/ROL,H_298,S_298,DLTH
C
      PHI_NASA7=ROL*TEM*TEM*TEM*TEM
C
      RETURN
C
      END
C
      DOUBLE PRECISION FUNCTION PHI_NASA_DCP1(TEM)
C
      implicit none
C
      DOUBLE PRECISION TEM,ROL,H_298,S_298,DLTH
      COMMON /GAS/ROL,H_298,S_298,DLTH
C
      PHI_NASA_DCP1=-2.D0*ROL/(TEM*TEM*TEM)
C
      RETURN
C
      END
C
      DOUBLE PRECISION FUNCTION PHI_NASA_DCP2(TEM)
C
      implicit none
C
      DOUBLE PRECISION TEM,ROL,H_298,S_298,DLTH
      COMMON /GAS/ROL,H_298,S_298,DLTH
C
      PHI_NASA_DCP2=-ROL/(TEM*TEM)
C
      RETURN
C
      END
C
      DOUBLE PRECISION FUNCTION PHI_NASA_DCP3(TEM)
C
      implicit none
C
      DOUBLE PRECISION TEM,ROL,H_298,S_298,DLTH
      COMMON /GAS/ROL,H_298,S_298,DLTH
C
      PHI_NASA_DCP3=0.D0
C
      RETURN
C
      END
C
C
      DOUBLE PRECISION FUNCTION PHI_NASA_DCP4(TEM)
C
      implicit none
C
      DOUBLE PRECISION TEM,ROL,H_298,S_298,DLTH
      COMMON /GAS/ROL,H_298,S_298,DLTH
C
      PHI_NASA_DCP4=ROL
C
      RETURN
C
      END
C
      DOUBLE PRECISION FUNCTION PHI_NASA_DCP5(TEM)
C
      implicit none
C
      DOUBLE PRECISION TEM,ROL,H_298,S_298,DLTH
      COMMON /GAS/ROL,H_298,S_298,DLTH
C
      PHI_NASA_DCP5=2.D0*ROL*TEM
C
      RETURN
C
      END
C
      DOUBLE PRECISION FUNCTION PHI_NASA_DCP6(TEM)
C
      implicit none
C
      DOUBLE PRECISION TEM,ROL,H_298,S_298,DLTH
      COMMON /GAS/ROL,H_298,S_298,DLTH
C
      PHI_NASA_DCP6=3.D0*ROL*TEM*TEM
C
      RETURN
C
      END
C
      DOUBLE PRECISION FUNCTION PHI_NASA_DCP7(TEM)
C
      implicit none
C
      DOUBLE PRECISION TEM,ROL,H_298,S_298,DLTH
      COMMON /GAS/ROL,H_298,S_298,DLTH
C
      PHI_NASA_DCP7=4.D0*ROL*TEM*TEM*TEM
C
      RETURN
C
      END
C
      DOUBLE PRECISION FUNCTION PHI_NASA_D2CP1(TEM)
C
      implicit none
C
      DOUBLE PRECISION TEM,ROL,H_298,S_298,DLTH
      COMMON /GAS/ROL,H_298,S_298,DLTH
C
      PHI_NASA_D2CP1=6.D0*ROL/(TEM*TEM*TEM*TEM)
C
      RETURN
C
      END
C
      DOUBLE PRECISION FUNCTION PHI_NASA_D2CP2(TEM)
C
      implicit none
C
      DOUBLE PRECISION TEM,ROL,H_298,S_298,DLTH
      COMMON /GAS/ROL,H_298,S_298,DLTH
C
      PHI_NASA_D2CP2=2.D0*ROL/(TEM*TEM*TEM)
C
      RETURN
C
      END
C
      DOUBLE PRECISION FUNCTION PHI_NASA_D2CP3(TEM)
C
      implicit none
C
      DOUBLE PRECISION TEM,ROL,H_298,S_298,DLTH
      COMMON /GAS/ROL,H_298,S_298,DLTH
C
      PHI_NASA_D2CP3=0.D0
C
      RETURN
C
      END
C
C
      DOUBLE PRECISION FUNCTION PHI_NASA_D2CP4(TEM)
C
      implicit none
C
      DOUBLE PRECISION TEM,ROL,H_298,S_298,DLTH
      COMMON /GAS/ROL,H_298,S_298,DLTH
C
      PHI_NASA_D2CP4=0.D0
C
      RETURN
C
      END
C
      DOUBLE PRECISION FUNCTION PHI_NASA_D2CP5(TEM)
C
      implicit none
C
      DOUBLE PRECISION TEM,ROL,H_298,S_298,DLTH
      COMMON /GAS/ROL,H_298,S_298,DLTH
C
      PHI_NASA_D2CP5=2.D0*ROL
C
      RETURN
C
      END
C
      DOUBLE PRECISION FUNCTION PHI_NASA_D2CP6(TEM)
C
      implicit none
C
      DOUBLE PRECISION TEM,ROL,H_298,S_298,DLTH
      COMMON /GAS/ROL,H_298,S_298,DLTH
C
      PHI_NASA_D2CP6=6.D0*ROL*TEM
C
      RETURN
C
      END
C
      DOUBLE PRECISION FUNCTION PHI_NASA_D2CP7(TEM)
C
      implicit none
C
      DOUBLE PRECISION TEM,ROL,H_298,S_298,DLTH
      COMMON /GAS/ROL,H_298,S_298,DLTH
C
      PHI_NASA_D2CP7=12.D0*ROL*TEM*TEM
C
      RETURN
C
      END
C
C
      DOUBLE PRECISION FUNCTION PHI_NASA_H(TEM,alpha)
C
      implicit none
C
      DOUBLE PRECISION TEM,ROL,alpha,H_298,S_298,DLTH
      DIMENSION alpha(7)
      COMMON /GAS/ROL,H_298,S_298,DLTH
C
!     WRITE(*,*) ALPHA
!     WRITE(*,*) TEM
!     WRITE(*,*) PHI_NASA_H
      PHI_NASA_H=((-1.0D0)*alpha(1)/TEM+alpha(2)*DLOG(TEM)+alpha(3)*TEM
     &          +alpha(4)*TEM*TEM/2.D0+alpha(5)*TEM*TEM*TEM/3.D0
     &          +alpha(6)*TEM*TEM*TEM*TEM/4.D0
     &          +alpha(7)*TEM*TEM*TEM*TEM*TEM/5.D0)*ROL
C
      RETURN
C
      END
C
      DOUBLE PRECISION FUNCTION PHI_NASA_S(TEM,alpha)
C
      implicit none
C
      DOUBLE PRECISION TEM,ROL,alpha,H_298,S_298,DLTH
      DIMENSION alpha(7)
      COMMON /GAS/ROL,H_298,S_298,DLTH
C
!     WRITE(*,*) ALPHA(1:7)
!     WRITE(*,*) TEM
      PHI_NASA_S=(-alpha(1)/TEM/TEM/2.D0-alpha(2)/TEM
     &          +alpha(3)*DLOG(TEM)+alpha(4)*TEM
     &          +alpha(5)*TEM*TEM/2.D0+alpha(6)*TEM*TEM*TEM/3.D0
     &          +alpha(7)*TEM*TEM*TEM*TEM/4.D0)*ROL
C
      RETURN
C
      END
C
      DOUBLE PRECISION FUNCTION DLTH298(alpha)
C
      implicit none
C
      DOUBLE PRECISION TEM,ROL,alpha,H_298,S_298,DLTH
      DIMENSION alpha(7)
      COMMON /GAS/ROL,H_298,S_298,DLTH
C
      TEM=298.15D0
      DLTH298=((-1.0D0)*alpha(1)/TEM+alpha(2)*DLOG(TEM)+alpha(3)*TEM
     &          +alpha(4)*TEM*TEM/2.D0+alpha(5)*TEM*TEM*TEM/3.D0
     &          +alpha(6)*TEM*TEM*TEM*TEM/4.D0
     &          +alpha(7)*TEM*TEM*TEM*TEM*TEM/5.D0)*ROL
C
      RETURN
C
      END