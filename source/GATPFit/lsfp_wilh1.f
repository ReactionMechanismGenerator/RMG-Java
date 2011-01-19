C WEN, August 2005,
C add the Cp for 0K
C**********************************************************************
C FOR WILHOID CHEMTHERMO DATE FORMAT, N=4, M>=4
C
      SUBROUTINE SWILT(SNAM,ENAM,ELNO,STRUC_MOL,ATOMS,
     &                   ROTORS,THERM1,TEMP,CPT,CH,CS,M)
C
      implicit none
C
      CHARACTER(LEN=4) ENAM(5)
      CHARACTER(LEN=16) SNAM
      CHARACTER(LEN=9) STRUC_MOL
      CHARACTER(LEN=1) ELNO(5)
      DOUBLE PRECISION ROL,H_298,S_298,DLTH
      DOUBLE PRECISION TEMP,CPT,CH,CS,PHI1,PHI2,PHI3,PHI4,B_WILH,
     &                   CP_INF,CP_0,ATOMS,ROTORS
      DOUBLE PRECISION CPT2,TEM2
      DOUBLE PRECISION Y_H,A_H,A_S,FUN_H,FUN_S
      DOUBLE PRECISION THERM1
C
      CHARACTER TRANS
      INTEGER M,N,NRHS,LDA,LDB,LWORK,INFO
      DOUBLE PRECISION A,B,WORK
      PARAMETER (N=4, NRHS=1)
C
C REQUIRED BY dgels:
! LDA >= MAX(1,M)
! LDB >= MAX(1,M,N)
! LWORK >= max( 1, MN + max( MN, NRHS ) )
! MN= MIN( M, N )
C
      PARAMETER (LWORK=8)
      DIMENSION A(M+1,N),B(M+1,NRHS),WORK(LWORK)
      DIMENSION TEMP(M),CPT(M),CH(M),CS(M),THERM1(6)
      DIMENSION TEM2(M+1),CPT2(M+1)
C
      INTEGER I,J
C
      COMMON /GAS/ROL,H_298,S_298,DLTH
      COMMON /CP/CP_INF,CP_0
C
      EXTERNAL PHI1,PHI2,PHI3,PHI4,B_WILH,FUN_H
C
C Input 100 temperature data
! TEMP(1:M) from group additivity
!     DATA TEMP/300.D0,400.D0,500.D0,600.D0,
!     &           800.D0,1000.D0,1500.D0/
!     DATA CPT/17.7725D0,22.4645D0,27.0142D0,30.7109D0,
!     &         36.9668D0,41.8009D0,49.3365D0/
C
! DEFINE THE TOTAL NUMBER OF DISCRETE DATA******************
!
!     WRITE(*,*) 'CALLED WILH'
C
      LDA=M+1
      LDB=M+1
      TRANS ='N'
C
      DO I=1,N
         DO J=1,LDA
            A(J,I)=0.D0
         ENDDO
      ENDDO
C
      DO I=1,NRHS
         DO J=1,LDB
            B(J,I)=0.D0
         ENDDO
      ENDDO
C
      DO I=1,LWORK
         WORK(I)=0.D0
      ENDDO
C
C Calculate the Cp_0 and Cp_int, here *******************
C gmagoon 10/2/09: added monoatomic case; this needs to be considered separately
C and was causing problems for H radical (presumably similar for O atom)
        IF (ATOMS .EQ. 1) THEN
            CP_0 = 2.5D0*1.9872D0
            CP_INF = 2.5D0*1.9872D0
      ELSE IF (STRUC_MOL .EQ. 'NONLINEAR') THEN
         CP_0  =4.0D0*1.9872D0
         CP_INF=(3.0D0*ATOMS-(2.D0+0.5D0*ROTORS))*1.9872D0
      ELSE
         CP_0  =3.5D0*1.9872D0
         CP_INF=(3.0D0*ATOMS-1.5D0)*1.9872D0
      ENDIF
!     WRITE(*,*) 'CP_INF=', CP_INF
!     WRITE(*,*) 'CP_0  =', CP_0
!     PAUSE
C********************************************************
C RESET 'TEMP' AND 'CPT'
      DO I=1,M+1
         TEM2(I)=0.D0
         CPT2(I)=0.D0
      ENDDO
!
      TEM2(1)  = 0.D0
      CPT2(1)  = CP_0
!
      DO I=1,M
         TEM2(I+1)=TEMP(I)
         CPT2(I+1)=CPT(I)
      ENDDO
! A(1:N,1:M+1)
      DO I=1,M+1
         A(I,1)=PHI1(TEM2(I))
         A(I,2)=PHI2(TEM2(I))
         A(I,3)=PHI3(TEM2(I))
         A(I,4)=PHI4(TEM2(I))
      ENDDO
! C(1:M+1), C(1:M+1,1)=CP(T)-CP_0-(CP_INF-CP_0)*Y*Y
      B(1,1)   = 0.D0
      DO I=1,M
        B(I+1,1)=B_WILH(TEM2(I+1),CPT2(I+1))
      ENDDO
!
C Contrain for Cp at 0K
!
!     B(1,1) = PHI1(0.D0)
!     B(1,2) = PHI2(0.D0)
!     B(1,3) = PHI3(0.D0)
!     B(1,4) = PHI4(0.D0)
!     WRITE(*,*) 'B=', B
!     PAUSE
!
!     D(1) = 0.D0
!
C: Call LS solver wihtout contrains
C: Solver   A*X=B
C
      call dgels(trans,m+1,n,nrhs,a,lda,b,ldb,work,lwork,info)
!     call DGGLSE(M+1,N,P,A,LDA,B,LDB,C,D,X1,WORK,
!     &             LWORK,INFO)
C
C The alpha coefficients are stored in b(1:4,1)
!     WRITE(*,*) B(1:4)
C
C Calculate the constant 'I' for Enthalpy
C
      Y_H=298.15D0/(298.15D0+500.D0)
C
      A_H=DLTH*1000.d0-(CP_0*298.15D0)+(CP_INF-CP_0)*298.15D0
     &   *FUN_H(298.15D0,B(1:4,1))
C
C Calculate the constant 'J' for Entropy
C
      A_S=S_298 - CP_INF*DLOG(298.15D0) + (CP_INF-CP_0)*
     &      FUN_S(298.15D0,B(1:4,1))
C
C
      DO I=1,4
         THERM1(I)=B(I,1)
      ENDDO
      THERM1(4+1)=A_H
      THERM1(4+2)=A_S
!     WRITE(*,*) THERM1(1:6)
C
!     WRITE(*,*) 'WILH OVER'
      RETURN
C
      END
C
C************************************************************************
      DOUBLE PRECISION FUNCTION PHI1(TEM)
C
      implicit none
C
      DOUBLE PRECISION TEM,Y,CP_INF,CP_0
      COMMON /CP/CP_INF,CP_0
C
      Y=TEM/(TEM+500.D0)
      PHI1=1.0D0*(CP_INF-CP_0)*Y*Y*(Y-1)
C
      RETURN
C
      END
C
      DOUBLE PRECISION FUNCTION PHI2(TEM)
C
      implicit none
C
      DOUBLE PRECISION TEM,Y,CP_INF,CP_0
      COMMON /CP/CP_INF,CP_0
C
      Y=TEM/(TEM+500.D0)
      PHI2=Y*(CP_INF-CP_0)*Y*Y*(Y-1)
C
      RETURN
C
      END
C
      DOUBLE PRECISION FUNCTION PHI3(TEM)
C
      implicit none
C
      DOUBLE PRECISION TEM,Y,CP_INF,CP_0
      COMMON /CP/CP_INF,CP_0
C
      Y=TEM/(TEM+500.D0)
      PHI3=Y*Y*(CP_INF-CP_0)*Y*Y*(Y-1)
C
      RETURN
C
      END
C
C
      DOUBLE PRECISION FUNCTION PHI4(TEM)
C
      implicit none
C
      DOUBLE PRECISION TEM,Y,CP_INF,CP_0
      COMMON /CP/CP_INF,CP_0
C
      Y=TEM/(TEM+500.D0)
      PHI4=Y*Y*Y*(CP_INF-CP_0)*Y*Y*(Y-1)
C
      RETURN
C
      END
C
      DOUBLE PRECISION FUNCTION B_WILH(TEM,CP)
C
      implicit none
C
      DOUBLE PRECISION TEM,Y,CP_INF,CP_0,CP
      COMMON /CP/CP_INF,CP_0
C
      Y=TEM/(TEM+500.D0)
      B_WILH=CP-CP_0-Y*Y*(CP_INF-CP_0)
C
      RETURN
C
      END
C
      DOUBLE PRECISION FUNCTION FUN_H(TEM,ALPHA)
C
      implicit none
C
      DOUBLE PRECISION TEM,ALPHA,Y
      DOUBLE PRECISION FUN1,FUN2,FUN3,FUN4,FUN5
      DIMENSION alpha(4)
C
      Y=TEM/(TEM+500.D0)
      IF (Y .EQ. 0.D0) Y=1.0D0
      FUN1=(2.D0+ALPHA(1)+ALPHA(2)+ALPHA(3)+ALPHA(4))*
     &       (Y/2.D0-1+(1.D0/Y-1)*DLOG(TEM+500.D0))
      FUN2=(1.D0/6.D0)*(3.D0*ALPHA(1)+ALPHA(2)+ALPHA(3)+ALPHA(4))
      FUN3=(Y/12.D0)*(4.D0*ALPHA(2)+ALPHA(3)+ALPHA(4))
      FUN4=(Y*Y/20.D0)*(5.D0*ALPHA(3)+ALPHA(4))
      FUN5=(Y*Y*Y/30.D0)*6.D0*ALPHA(4)
C
      FUN_H=FUN1+Y*Y*(FUN2+FUN3+FUN4+FUN5)
C
      RETURN
C
      END
      DOUBLE PRECISION FUNCTION FUN_S(TEM,ALPHA)
C
      implicit none
C
      DOUBLE PRECISION TEM,ALPHA,Y
      DIMENSION alpha(4)
C
      Y=TEM/(TEM+500.D0)
      IF (Y .EQ. 0.D0) Y=1.0D0
C
      FUN_S=DLOG(Y)+
     &     (1.D0+Y*(ALPHA(1)/2.D0+ALPHA(2)*Y/3.D0+ALPHA(3)*Y*Y/4.D0
     &     +ALPHA(4)*Y*Y*Y/5.D0))*Y
C
      RETURN
C
      END
C
      DOUBLE PRECISION FUNCTION CP_WILH(TEM,BETA)
C
      implicit none
C
      DOUBLE PRECISION TEM,BETA,Y,CP_INF,CP_0
      DIMENSION BETA(6)
      COMMON /CP/CP_INF,CP_0
C
      Y=TEM/(TEM+500.D0)
      CP_WILH=CP_0+(CP_INF-CP_0)*Y*Y
     &         *(1.D0+(Y-1.D0)*(BETA(1)+BETA(2)*Y+BETA(3)*Y*Y
     &       +BETA(4)*Y*Y*Y))
C
      RETURN
C
      END
C
      DOUBLE PRECISION FUNCTION H_WILH(TEM,BETA)
C
      implicit none
C
      DOUBLE PRECISION TEM,BETA,Y,CP_INF,CP_0,FUN_H
      DIMENSION BETA(6)
      EXTERNAL FUN_H
      COMMON /CP/CP_INF,CP_0
C
      Y=TEM/(TEM+500.D0)
      H_WILH=BETA(5)
     &        +(CP_0*TEM)-(CP_INF-CP_0)*TEM*
     &      FUN_H(TEM,BETA(1:4))
C
      RETURN
C
      END
C
      DOUBLE PRECISION FUNCTION S_WILH(TEM,BETA)
C
      implicit none
C
      DOUBLE PRECISION TEM,BETA,Y,CP_INF,CP_0,FUN_S
      DIMENSION BETA(6)
      COMMON /CP/CP_INF,CP_0
      EXTERNAL FUN_S
C
      Y=TEM/(TEM+500.D0)
      IF (TEM .EQ. 0.D0) TEM = 1.0D0
      S_WILH=BETA(6)
     &      +CP_INF*DLOG(TEM) - (CP_INF-CP_0)*
     &      FUN_S(TEM,BETA(1:4))
C
      RETURN
C
      END
      SUBROUTINE DATAFIND(BETA,TEM,CP,CH,CS)
C
      implicit none
C
      DOUBLE PRECISION BETA,TEM,CP,CH,CS
      DOUBLE PRECISION CP_WILH,H_WILH,S_WILH
      DIMENSION BETA(6)
C
      EXTERNAL CP_WILH,H_WILH,S_WILH
C
      CP=CP_WILH(TEM,BETA)
      CH=H_WILH(TEM,BETA)
      CS=S_WILH(TEM,BETA)
      RETURN
C
      END