C: Add the constrains to the first and second derivatives 
C: Therefore, there are three constrains: Cp, DCp, and D2Cp
C: --------------------------------------------------------------
C: For CHEMKIN
C: Use a single matrix to calculate all 10 coefficients
C: [A]: M x 10 
C: [X]: 10 x 1 (alpha1_low,alpha2_low,alpha3_low,alpha4_low,alpha5_low,
C:              alpha1_hig,alpha2_hig,alpha3_hig,alpha4_hig,alpha5_hig)
C: [C]: M x 1
C: Contrains: CP, dCP/dT, and d2CP/d2T
C:            at the intermediate temperature Tint
C: [B]: 3 x 10
C: [D]: 3 x 1
C----------------------------------------------------------------
C FOR CHEMKIN CHEMTHERMO DATA FORMAT
C number of coefficients: N=5 x 2
C number of discrete data point: M>=N, i.e., M=102
C 
C use the DGGLSE to solve the linear least squares equality-constrained  
C (LlSE) problem.
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C: Need to change the value of M according to the size of sample data
!********************************************************************
! Programmed by John Z. Wen, MIT, June 2005
!********************************************************************
C
      SUBROUTINE CHEM(Tint,TNEW,THERM1,THERM2,TEMP1,C1,CH1,CS1)
C
      implicit none
C
      DOUBLE PRECISION Tint,TNEW,THERM1,THERM2
      DOUBLE PRECISION H_298,S_298,DLTH
C
      INTEGER M,N,P,LDA,LDB,LWORK,INFO,M1
      INTEGER MH,NH,PH,LDAH,LDBH,LWORKH
      INTEGER MS,NS,PS,LDAS,LDBS,LWORKS
      DOUBLE PRECISION A1,B1,Bk,C1,D1,X1,WORK1
      DOUBLE PRECISION TEMP1,PHI_CHEM1,PHI_CHEM2,PHI_CHEM3,
     &                   PHI_CHEM4,PHI_CHEM5,ROL
      DOUBLE PRECISION PHI_CHEM_DCP1,PHI_CHEM_DCP2,PHI_CHEM_DCP3,
     &           PHI_CHEM_DCP4,PHI_CHEM_DCP5,PHI_CHEM_D2CP1,
     &             PHI_CHEM_D2CP2,PHI_CHEM_D2CP3,PHI_CHEM_D2CP4,
     &         PHI_CHEM_D2CP5
      DOUBLE PRECISION A_H1,A_S1,A_H2,A_S2
      DOUBLE PRECISION AH1,BH1,CH1,DH1,XH1,WORKH,BHk,CH,CH2
      DOUBLE PRECISION AS1,BS1,CS1,DS1,XS1,WORKS,BSk,CS,CS2
      DOUBLE PRECISION PHI_CHEM_H,PHI_CHEM_S
      INTEGER MARK1
      DOUBLE PRECISION TEMP,C
C
C Required by LSEC solver:
! LDA >= MAX(1,M)
! LDB >= MAX(1,P)
! LWORK >= max(1,M+N+P)
! P=1 for only one constrain at the intermediate temp.
! P=3 for three constrains at the intermediate temp.
C
C CHANGE 'LWORK' if NECESSARY
      PARAMETER (M=102)
      PARAMETER (LDB=3,P=3)
      PARAMETER (N=10)
      DIMENSION A1(M,N),B1(LDB,N),C1(M-1),D1(P),X1(N),WORK1(M+N+P),
     &            Bk(LDB,N)
      PARAMETER (LDBH=2,PH=2)
      PARAMETER (NH=2)
      PARAMETER (LDBS=2,PS=2)
      PARAMETER (NS=2)
      DIMENSION AH1(M,NH),BH1(LDBH,NH),CH1(M-1),DH1(PH),XH1(NH),
     &            WORKH(M+NH+PH),BHk(LDBH,NH),CH(M),CH2(M)
      DIMENSION AS1(M,NS),BS1(LDBS,NS),CS1(M-1),DS1(PS),XS1(NS),
     &            WORKS(M+NS+PS),BSk(LDBS,NS),CS(M),CS2(M)
      DIMENSION TEMP1(M-1),THERM1(7),THERM2(7)
      DIMENSION TEMP(M),C(M)
      INTEGER I,J
C
! GAS CONSTANT
      COMMON /GAS/ROL,H_298,S_298,DLTH
C
      EXTERNAL PHI_CHEM1,PHI_CHEM2,PHI_CHEM3,PHI_CHEM4,PHI_CHEM5,
     &         PHI_CHEM_DCP1,PHI_CHEM_DCP2,PHI_CHEM_DCP3,
     &         PHI_CHEM_DCP4,PHI_CHEM_DCP5,PHI_CHEM_D2CP1,
     &             PHI_CHEM_D2CP2,PHI_CHEM_D2CP3,PHI_CHEM_D2CP4,
     &         PHI_CHEM_D2CP5,PHI_CHEM_H,PHI_CHEM_S
C
!     WRITE(*,*) 'CALLED CHEM'
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
      DO I=1,M-1
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
         C (I)=C1(I)
         CH(I)=CH1(I)
         CS(I)=CS1(I)
      ENDDO
      TEMP(MARK1+1)=TEMP(MARK1)
      C (MARK1+1)=C (MARK1)
      CH(MARK1+1)=CH(MARK1)
      CS(MARK1+1)=CS(MARK1)
      DO I=MARK1+2,M
         TEMP(I)=TEMP1(I-1)
         C (I)=C1 (I-1)
         CH(I)=CH1(I-1)
         CS(I)=CS1(I-1)
      ENDDO
C
!     WRITE(*,*) 'Tint=', Tint
!     WRITE(*,*) TEMP
!     WRITE(*,*) C
!     WRITE(*,*) CH
!     WRITE(*,*) CS
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
! For PHI matrix, first temperature range
        DO I=1,MARK1
           A1(I,1)=PHI_CHEM1(TEMP(I))
           A1(I,2)=PHI_CHEM2(TEMP(I))
           A1(I,3)=PHI_CHEM3(TEMP(I))
           A1(I,4)=PHI_CHEM4(TEMP(I))
           A1(I,5)=PHI_CHEM5(TEMP(I))
        ENDDO
! For PHI matrix, second temperature range 
        DO I=MARK1+1,M
           A1(I,6) =PHI_CHEM1(TEMP(I))
           A1(I,7) =PHI_CHEM2(TEMP(I))
           A1(I,8) =PHI_CHEM3(TEMP(I))
           A1(I,9) =PHI_CHEM4(TEMP(I))
           A1(I,10)=PHI_CHEM5(TEMP(I))
        ENDDO
C
!     WRITE(*,*) TEMP1
!     WRITE(*,*) C1
!     PAUSE
C
C For Constrains at Tint ------------------------
C
! For CP values
        Bk(1,1) = PHI_CHEM1(Tint)
        Bk(1,2) = PHI_CHEM2(Tint)
        Bk(1,3) = PHI_CHEM3(Tint)
        Bk(1,4) = PHI_CHEM4(Tint)
        Bk(1,5) = PHI_CHEM5(Tint)
        Bk(1,6) =-PHI_CHEM1(Tint)
        Bk(1,7) =-PHI_CHEM2(Tint)
        Bk(1,8) =-PHI_CHEM3(Tint)
        Bk(1,9) =-PHI_CHEM4(Tint)
        Bk(1,10)=-PHI_CHEM5(Tint)
        D1(1)=0.D0
C For DCp values
        Bk(2,1) = PHI_CHEM_DCP1(Tint)
        Bk(2,2) = PHI_CHEM_DCP2(Tint)
        Bk(2,3) = PHI_CHEM_DCP3(Tint)
        Bk(2,4) = PHI_CHEM_DCP4(Tint)
        Bk(2,5) = PHI_CHEM_DCP5(Tint)
        Bk(2,6) =-PHI_CHEM_DCP1(Tint)
        Bk(2,7) =-PHI_CHEM_DCP2(Tint)
        Bk(2,8) =-PHI_CHEM_DCP3(Tint)
        Bk(2,9) =-PHI_CHEM_DCP4(Tint)
        Bk(2,10)=-PHI_CHEM_DCP5(Tint)
        D1(2)=0.D0
C For D2Cp values
        Bk(3,1) = PHI_CHEM_D2CP1(Tint)
        Bk(3,2) = PHI_CHEM_D2CP2(Tint)
        Bk(3,3) = PHI_CHEM_D2CP3(Tint)
        Bk(3,4) = PHI_CHEM_D2CP4(Tint)
        Bk(3,5) = PHI_CHEM_D2CP5(Tint)
        Bk(3,6) =-PHI_CHEM_D2CP1(Tint)
        Bk(3,7) =-PHI_CHEM_D2CP2(Tint)
        Bk(3,8) =-PHI_CHEM_D2CP3(Tint)
        Bk(3,9) =-PHI_CHEM_D2CP4(Tint)
        Bk(3,10)=-PHI_CHEM_D2CP5(Tint)
        D1(3)=0.D0
C
      DO I=1,N
         DO J=1,LDB
            B1(J,I)=Bk(J,I)
         ENDDO
      ENDDO
C
! CHEM FORMAT-----------------------------------------
C
      call DGGLSE(M,N,P,A1,LDA,B1,LDB,C,D1,X1,WORK1,LWORK,INFO)
!     WRITE(*,*) 'X1=',X1(1:10)
!     pause
C
C For enthalpy, NEED define 2 coeffs, N=2 
C
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
         CH2(I)=CH(I) - PHI_CHEM_H(TEMP(I),X1(1:5))
      ENDDO
      DO I=MARK1+1,MH
         CH2(I)=CH(I) - PHI_CHEM_H(TEMP(I),X1(6:10))
      ENDDO
C
      DO I=1,MARK1
         AH1(I,1)=ROL/1000.d0
      ENDDO
      DO I=MARK1+1,MH
         AH1(I,2)=ROL/1000.d0
      ENDDO
C
C Contrain for H at Tint
C
      BHk(1,1)=ROL/1000.d0
      BHk(1,2)=-ROL/1000.d0
      DH1(1)=PHI_CHEM_H(Tint,X1(6:10))
     &      -PHI_CHEM_H(Tint,X1(1:5))
C
C Contrain for H at 298.15K
C
      BHk(2,1)=ROL/1000.d0
      BHk(2,2)=0.D0
      DH1(2)=DLTH - PHI_CHEM_H(298.15D0,X1(1:5))
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
C
      A_H1=XH1(1)
      A_H2=XH1(2)
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
         CS2(I)=CS(I) - PHI_CHEM_S(TEMP(I),X1(1:5))
      ENDDO
      DO I=MARK1+1,MS
         CS2(I)=CS(I) - PHI_CHEM_S(TEMP(I),X1(6:10))
      ENDDO
      DO I=1,MARK1
         AS1(I,1)=ROL
      ENDDO
      DO I=MARK1+1,MS
         AS1(I,2)=ROL
      ENDDO
C
C Contrain for S at Tint
C
      BSk(1,1)=ROL
      BSk(1,2)=-ROL
      DS1(1)=PHI_CHEM_S(Tint,X1(6:10))
     &      -PHI_CHEM_S(Tint,X1(1:5))
C
C Contrain for S at 298.15K
C
      BSk(2,1)=ROL
      BSk(2,2)=0.D0
      DS1(2)=S_298 - PHI_CHEM_S(298.15D0,X1(1:5))
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
C
      A_S1=XS1(1)
      A_S2=XS1(2)
C
C
        DO I=1,5
           THERM1(I)=X1(I)
        ENDDO
        THERM1(5+1)=A_H1
        THERM1(5+2)=A_S1
        DO I=1,5
           THERM2(I)=X1(I+5)
        ENDDO
        THERM2(5+1)=A_H2
        THERM2(5+2)=A_S2
C
!     WRITE(*,*) 'CHEM OVER'
C
      RETURN
C
      END
C
      DOUBLE PRECISION FUNCTION PHI_CHEM1(TEM)
C
      implicit none
C
      DOUBLE PRECISION TEM,ROL,H_298,S_298,DLTH
      COMMON /GAS/ROL,H_298,S_298,DLTH
C
      PHI_CHEM1=ROL
C
      RETURN
C
      END
C
      DOUBLE PRECISION FUNCTION PHI_CHEM2(TEM)
C
      implicit none
C
      DOUBLE PRECISION TEM,ROL,H_298,S_298,DLTH
      COMMON /GAS/ROL,H_298,S_298,DLTH
C
      PHI_CHEM2=ROL*TEM
C
      RETURN
C
      END
C
      DOUBLE PRECISION FUNCTION PHI_CHEM3(TEM)
C
      implicit none
C
      DOUBLE PRECISION TEM,ROL,H_298,S_298,DLTH
      COMMON /GAS/ROL,H_298,S_298,DLTH
C
      PHI_CHEM3=ROL*TEM*TEM
C
      RETURN
C
      END
C
C
      DOUBLE PRECISION FUNCTION PHI_CHEM4(TEM)
C
      implicit none
C
      DOUBLE PRECISION TEM,ROL,H_298,S_298,DLTH
      COMMON /GAS/ROL,H_298,S_298,DLTH
C
      PHI_CHEM4=ROL*TEM*TEM*TEM
C
      RETURN
C
      END
C
      DOUBLE PRECISION FUNCTION PHI_CHEM5(TEM)
C
      implicit none
C
      DOUBLE PRECISION TEM,ROL,H_298,S_298,DLTH
      COMMON /GAS/ROL,H_298,S_298,DLTH
C
      PHI_CHEM5=ROL*TEM*TEM*TEM*TEM
C
      RETURN
C
      END
C
      DOUBLE PRECISION FUNCTION PHI_CHEM_DCP1(TEM)
C
      implicit none
C
      DOUBLE PRECISION TEM,ROL,H_298,S_298,DLTH
      COMMON /GAS/ROL,H_298,S_298,DLTH
C
      PHI_CHEM_DCP1=0.D0
C
      RETURN
C
      END
C
      DOUBLE PRECISION FUNCTION PHI_CHEM_DCP2(TEM)
C
      implicit none
C
      DOUBLE PRECISION TEM,ROL,H_298,S_298,DLTH
      COMMON /GAS/ROL,H_298,S_298,DLTH
C
      PHI_CHEM_DCP2=ROL
C
      RETURN
C
      END
C
      DOUBLE PRECISION FUNCTION PHI_CHEM_DCP3(TEM)
C
      implicit none
C
      DOUBLE PRECISION TEM,ROL,H_298,S_298,DLTH
      COMMON /GAS/ROL,H_298,S_298,DLTH
C
      PHI_CHEM_DCP3=2.D0*ROL*TEM
C
      RETURN
C
      END
C
      DOUBLE PRECISION FUNCTION PHI_CHEM_DCP4(TEM)
C
      implicit none
C
      DOUBLE PRECISION TEM,ROL,H_298,S_298,DLTH
      COMMON /GAS/ROL,H_298,S_298,DLTH
C
      PHI_CHEM_DCP4=3.D0*ROL*TEM*TEM
C
      RETURN
C
      END
C
      DOUBLE PRECISION FUNCTION PHI_CHEM_DCP5(TEM)
C
      implicit none
C
      DOUBLE PRECISION TEM,ROL,H_298,S_298,DLTH
      COMMON /GAS/ROL,H_298,S_298,DLTH
C
      PHI_CHEM_DCP5=4.D0*ROL*TEM*TEM*TEM
C
      RETURN
C
      END
C
      DOUBLE PRECISION FUNCTION PHI_CHEM_D2CP1(TEM)
C
      implicit none
C
      DOUBLE PRECISION TEM,ROL,H_298,S_298,DLTH
      COMMON /GAS/ROL,H_298,S_298,DLTH
C
      PHI_CHEM_D2CP1=0.D0
C
      RETURN
C
      END
C
      DOUBLE PRECISION FUNCTION PHI_CHEM_D2CP2(TEM)
C
      implicit none
C
      DOUBLE PRECISION TEM,ROL,H_298,S_298,DLTH
      COMMON /GAS/ROL,H_298,S_298,DLTH
C
      PHI_CHEM_D2CP2=0.D0
C
      RETURN
C
      END
C
      DOUBLE PRECISION FUNCTION PHI_CHEM_D2CP3(TEM)
C
      implicit none
C
      DOUBLE PRECISION TEM,ROL,H_298,S_298,DLTH
      COMMON /GAS/ROL,H_298,S_298,DLTH
C
      PHI_CHEM_D2CP3=2.D0*ROL
C
      RETURN
C
      END
C
      DOUBLE PRECISION FUNCTION PHI_CHEM_D2CP4(TEM)
C
      implicit none
C
      DOUBLE PRECISION TEM,ROL,H_298,S_298,DLTH
      COMMON /GAS/ROL,H_298,S_298,DLTH
C
      PHI_CHEM_D2CP4=6.D0*ROL*TEM
C
      RETURN
C
      END
C
      DOUBLE PRECISION FUNCTION PHI_CHEM_D2CP5(TEM)
C
      implicit none
C
      DOUBLE PRECISION TEM,ROL,H_298,S_298,DLTH
      COMMON /GAS/ROL,H_298,S_298,DLTH
C
      PHI_CHEM_D2CP5=12.D0*ROL*TEM*TEM
C
      RETURN
C
      END
C
      DOUBLE PRECISION FUNCTION PHI_CHEM_H(TEM,alpha)
C
      implicit none
C
      DOUBLE PRECISION TEM,ROL,alpha,H_298,S_298,DLTH
      DIMENSION alpha(5)
      COMMON /GAS/ROL,H_298,S_298,DLTH
C
      PHI_CHEM_H=(alpha(1)*TEM+alpha(2)*TEM*TEM/2.D0
     &          +alpha(3)*TEM*TEM*TEM/3.D0+alpha(4)*TEM*TEM*TEM*TEM/4.D0
     &          +alpha(5)*TEM*TEM*TEM*TEM*TEM/5.D0)*ROL/1000.d0
C
      RETURN
C
      END
C
      DOUBLE PRECISION FUNCTION PHI_CHEM_S(TEM,alpha)
C
      implicit none
C
      DOUBLE PRECISION TEM,ROL,alpha,H_298,S_298,DLTH
      DIMENSION alpha(5)
      COMMON /GAS/ROL,H_298,S_298,DLTH
C
!     WRITE(*,*) ALPHA(1:5)
!     WRITE(*,*) TEM
      PHI_CHEM_S=(alpha(1)*DLOG(TEM)+alpha(2)*TEM
     &          +alpha(3)*TEM*TEM/2.D0+alpha(4)*TEM*TEM*TEM/3.D0
     &          +alpha(5)*TEM*TEM*TEM*TEM/4.D0)*ROL
C
      RETURN
C
      END
