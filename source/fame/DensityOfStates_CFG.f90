MODULE DensityOfStates

CONTAINS

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! This module contains all the subroutines necessary to calculate the density
! of states.
!
! The first subroutine is CElliptic.  This routine was taken by Argonne (??);
! it calculates the complete elliptical integral of the first and second kind,
! which is needed to calculated the density of states for a hindered rotor.
!
! The second subroutine is calc_hr_rho.  This subroutine calculates the
! density of states for a hindered rotor.  It needs more commenting, which
! CFG will do upon returning to MIT.
!
! The third subroutine is convolution_rho.  This subroutine performs a  
! convolution intregral for two density of states.  It is used to 
! convolute the density of states for multiple hindered rotor density of states
! as well as convoluting the density of states for the active K-rotor and the 
! torsional rotors.
!
! The forth subroutine is beyer.  This subroutine performs the 
! Beyer-Swinehart algorithm to calculate the vibrational density of states.
! Currently, it is used in the method described by Astholz; in this method
! the B.S. algorithm is initialized with the density of states for all the 
! other active degrees of freedom (i.e. active K-rotor, hindered rotors, etc).
! The B.S. algorithm then convolutes the active d.o.s. with the vibrational
! modes.
!
! The fifth subroutine is calculate_rho.  This subroutine is the heart of the
! the module.  It takes as input a filename, the maximum energy value, the 
! bin spacing, and an empty array.  It return the array with the 
! density of states.
!
! IMPORTANT:  Currently this code is set to ignore all multiplicative constants
! That would otherwise be in front of the density of states, such as rotational
! constants, external symmetry numbers, etc.
! Why?  1) RMG does not have these numbers, and 
! 2) As long as we do Inverse Laplace, they will cancel out.
! One area where they do no cancel out is in the detailed balance calculations.
! Here we must assume that the rotational constants, etc, are equivalent
! between the isomers.  Not good, but the best we can do.
!------------------------------------------------------------------------------

!******************************************************
!* Complete elliptic integral of the first and second *
!* kind. The input parameter is xk, which should be   *
!* between 0 and 1. Technique uses Gauss' formula for *
!* the arithmogeometrical mean. e is a measure of the * 
!* convergence accuracy. The returned values are e1,  *
!* the elliptic integral of the first kind, and e2,   *
!* the elliptic integral of the second kind.          *
!* -------------------------------------------------- *
!* Reference: Ball, algorithms for RPN calculators.   *
!******************************************************
SUBROUTINE CElliptic(e,m,K,n)  
! Label: 100
  REAL(8) e,m,K,pi
  REAL(8)  A(0:99), B(0:99)
  INTEGER j,l,n
  pi = 4.d0*datan(1.d0)
  A(0)=1.d0+m ; B(0)=1.d0-m
  n=0
  if (m < 0.d0) return
  if (m > 1.d0) return
  if (e <= 0.d0) return
100  n = n + 1
  ! Generate improved values
  A(n)=(A(n-1)+B(n-1))/2.d0
  B(n)=dsqrt(A(n-1)*B(n-1))
  if (dabs(A(n)-B(n)) > e) goto 100
  K=pi/2.d0/A(n)
  l=1
  
  return
END SUBROUTINE CElliptic
!------------------------------------------------------------------------------
! This subroutine calculates the density of state for a hindered internal rotor
! Using the formulation of Forst.  See my notes for more details.
SUBROUTINE calc_HR_rho(energy, nu, V0, rho)
IMPLICIT NONE

REAL(8), DIMENSION(:), INTENT(IN) :: energy
REAL(8), INTENT(IN) :: nu
REAL(8), INTENT(IN) :: V0
REAL(8), DIMENSION(:), INTENT(OUT) :: rho

REAL(8) ::m, K, e
INTEGER :: n, i

!CFG:  WILL ADD SOME COMMENTS WHEN BACK AT MIT
DO i = 1, size(energy)

   IF (energy(i) < nu/2.0) THEN
   !   rho(i) = 0.0
      m = 0
      e=1.d-7
      call CElliptic(e,m,K,n)
      rho(i) = K/SQRT(V0)


   ELSEIF ( (energy(i) - nu/2) < V0) THEN
      m = SQRT( ( energy(i) - nu/2.0)/V0 )
      e=1.d-7
      call CElliptic(e,m,K,n)
      rho(i) = K/SQRT(V0)

   ELSEIF  ( (energy(i) - nu/2) == V0) THEN
      m = 1.0 - 1.e-7
      e=1.d-7
      call CElliptic(e,m,K,n)
      rho(i) = K/SQRT(V0)
      m = 1.0 - 1e-7
      e=1.d-7
      call CElliptic(e,m,K,n)
      rho(i) =( rho(i) +  K/SQRT( energy(i) - nu/2.0 ) )/2.0

   ELSE
      m = SQRT( V0 / ( energy(i) - nu/2.0) )
      e=1.d-7
      call CElliptic(e,m,K,n)
      rho(i) = K/SQRT( energy(i) - nu/2.0 )
   
   ENDIF

ENDDO

END SUBROUTINE calc_HR_rho
!------------------------------------------------------------------------------
! This subroutine takes two density of states and convolutes them.
SUBROUTINE convolution_rho(rho_1, rho_2, delta_nu, rho_out)
IMPLICIT NONE

REAL(8), DIMENSION(:), INTENT(IN) :: rho_1, rho_2
REAL(8), DIMENSION(size(rho_1)) :: rho_mid
REAL(8), DIMENSION(:), INTENT(OUT) :: rho_out
INTEGER, INTENT(IN) :: delta_nu
INTEGER i, j, k
REAL(8), ALLOCATABLE, DIMENSION(:) :: f

rho_out = rho_2
DO i = 1, size(rho_1)
   ALLOCATE(f(i))  
   DO j = 1,i
      f(j) = rho_2(i - j +1) * rho_1(j) * REAL(delta_nu)
   ENDDO
   rho_out(i) = SUM(f)
   DEALLOCATE(f)
ENDDO

END SUBROUTINE convolution_rho
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! This subroutine is the Beyer-Swinehart algorithm.
SUBROUTINE beyer(nu, Delta_nu, E_max, T)
IMPLICIT NONE

REAL(8), DIMENSION(:), INTENT(IN) :: nu 
INTEGER, INTENT(IN) :: Delta_nu 
INTEGER,INTENT(IN) :: E_Max
REAL(8),  DIMENSION(:), INTENT(INOUT) :: T
INTEGER, DIMENSION(size(nu)) :: R
INTEGER i, j
INTEGER M

! Convert the incoming frequencies into integers divided by the bin size
DO i = 1, SIZE(nu)
   R(i) = NINT( nu(i) / Delta_nu )
ENDDO

! M is the number of elements in the T array
M = E_max / Delta_nu +1

! The Beyer-Swinehart Algorithm to calculate rho
DO i = 1, SIZE(nu)
   DO j = 1, ( M - R(i) )
      T(R(i) + j) = T(R(i) + j) + T(j)
   ENDDO
ENDDO

END SUBROUTINE beyer

!------------------------------------------------------------------------------
! This is the main subroutine to calculate the total density of states 
! for all the active degrees of freedom.
SUBROUTINE Calculate_rho(filename, E_max, delta_nu, rho)

IMPLICIT NONE

! These variables are specific to the molecule
  INTEGER :: N_atoms
  INTEGER :: linearity
  INTEGER :: N_RRHO
  INTEGER :: N_rotors
  REAL(8), ALLOCATABLE, DIMENSION(:) :: Total_harm_osc_freq
  REAL(8), ALLOCATABLE, DIMENSION(:) :: Hind_rot_nu
  REAL(8), ALLOCATABLE, DIMENSION(:) :: Hind_rot_V
  REAL(8) :: K_rot_const
  INTEGER :: i, j

 ! These variables are used to open the file
  INTEGER :: OpenStatus
  CHARACTER(20) :: filename

! These variables are used in the computation of the density of states
  REAL(8), ALLOCATABLE, DIMENSION(:)    :: energy
  REAL(8), ALLOCATABLE, DIMENSION(:)    :: K_rot_rho
  REAL(8), ALLOCATABLE, DIMENSION(:)    :: HR_K_rot_rho
  REAL(8), ALLOCATABLE, DIMENSION(:,:)  :: HR_rho
  REAL(8), DIMENSION(:), INTENT(INOUT)  :: rho
  INTEGER, INTENT(IN) :: Delta_nu  
  INTEGER, INTENT(IN) :: E_Max

! Open the input file containing the data
  open (UNIT = 22, FILE = filename, STATUS = 'OLD', & 
       ACTION = 'READ', IOSTAT = OpenStatus)

! Read the first three lines, which should be the
! The number of atoms, the number of rotors, and the linearity 
  READ (22,*) N_atoms, N_rotors, linearity

! Depending upon the linearity of the molecule:
! Determine whether the total number of rigid-rotor harmonic-oscillator
! (RRHO) vibrational modes is 3N-5 or 3N-6
! (minus the number of rotors, naturally)
  IF (linearity<0) THEN
     WRITE(*,*) 'ERROR!  Linearity is less than zero!'
  ELSE IF (linearity==0) THEN
     WRITE(*,*) 'Linear molecule'
     N_RRHO = 3 * (N_atoms) - 5 - (N_rotors)
    ELSE IF (linearity==1) THEN
     WRITE(*,*) 'Nonlinear molecule'
     N_RRHO = 3 * (N_atoms) - 6 - (N_rotors)
  ELSE
     WRITE(*,*) 'ERROR!  Linearity greater than one!'
  END IF

! Allocate the array containing all the
! rigid-rotor harmonic oscillator frequencies
  ALLOCATE( Total_harm_osc_freq( N_RRHO ) )
! Read the frequencies from the input file
  READ (22,*) Total_harm_osc_freq

! Allocate the array containing all the hindered-rotor frequencies
  ALLOCATE( hind_rot_nu (N_rotors) )
! Allocate the array containing all the hindered-rotor barrier heights
  ALLOCATE( hind_rot_V (N_rotors) )
! Read the hindered-rotor frequencies and barrier heights
  DO i = 1, N_rotors
     READ(22,*) hind_rot_nu(i), hind_rot_V(i)
  ENDDO

! Close the input file
  CLOSE(22)

! Allocate the energy
  ALLOCATE(energy( (E_max/Delta_nu) + 1) )
  
! initialize the energy: energy goes from 0 to E_max in increments of 
! Delta_nu
  energy(1) = 0.0
  DO i = 1, (E_max/Delta_nu)
     energy(i+1) = energy(i) + Delta_nu
  END DO

!----------------------------------------------------------------------------
! BEGIN CALCULATING THE TOTAL DENSITY OF STATES
!----------------------------------------------------------------------------

! CALCULATE THE DENSITY OF STATES FOR THE ACTIVE K-ROTOR
! This assumes that the molecule is a symmetric top in which the unique
! axis exchanges energy with the vibrational modes, but the two
! degenerate rotations are inactive.
! The density of state for a one-dimensional rotor is:
! rho = DOUBLE CHECK THIS!!!

! Allocate the vectors containing the active rotor density of state
  ALLOCATE(K_rot_rho( size(energy) ) )

! Calculate the K-rotor density of states
! Assume that number of states at E = 0 is 1
  K_rot_rho(1) = 1.0
! Set the rotational constant = 1
  K_rot_const = 1.0
  DO i = 2,size(energy)
     K_rot_rho(i) = K_rot_const * SQRT( energy(i) )
  ENDDO

! CALCULATE THE HINDERED-ROTOR DENSITY OF STATES
! The first N_rotor column vectors contain
! the density of states for each hindered rotor.
! The next N_rotor - 2 column vectors contain the successive convolution
! integrals of the density of states.
! The final vector is the total convolution of all the hindered rotors.

! Allocate the vectors containing the hindered-rotor density of states
  ALLOCATE(HR_rho((2*N_rotors-1), size(K_rot_rho) ) )

! Start by making sure that there is at least one hindered rotor.
! If there is no hindered rotor, do nothing.
  IF (n_rotors>0) THEN

     DO i = 1, N_rotors
        ! calculate the density of state for each rotor
        call calc_HR_rho(energy, hind_rot_nu(i), hind_rot_V(i), HR_rho(i,:) )
     ENDDO

     ! Now that we have the individual density of states, begin the convolution
     DO i = 1, (2*N_rotors-3),2
        j = CEILING(REAL(i/2.0)) + N_rotors
        CALL convolution_rho(HR_rho(i,:), HR_rho(i+1,:), delta_nu, HR_rho(j,:) )
     ENDDO
  ENDIF

! Deallocate energy
  DEALLOCATE( energy )

! CALCULATE THE HINDERED-ROTOR, TRANSLATIONAL, AND ROTATIONAL DENSITY OF STATES
! Allocate the vectors containing the active-rotor-hindered-rotor
! density of states
  ALLOCATE(HR_K_rot_rho( size(K_rot_rho) ) )

! If there are no hindered rotors, then the density of state remains unchanged
  IF (N_rotors==0) THEN
     DO i = 1, size(K_rot_rho)
        HR_K_rot_rho(i) = K_rot_rho(i)
     ENDDO
  ! If there are hindered rotors, convolute the translation-rotational
  ! density of state with the total hindered-rotor density of state
  ELSE
     call convolution_rho(HR_rho((2*N_rotors-1), :), K_rot_rho, delta_nu, HR_K_rot_rho)
  ENDIF

! Deallocate the K-rotor and hindered rotor density of states
  DEALLOCATE(K_rot_rho)
  DEALLOCATE(HR_rho)

! CALCULATE THE TOTAL DENSITY OF STATES
! Initialize the total density of state with the total HR density of state
  rho = HR_k_rot_rho

! Deallocate the K-rotor-hindered-rotor density of states
  DEALLOCATE( HR_k_rot_rho )
  
! Call the Beyer-Swinehart Algorithm, which is being initialized with the 
! density of states for the combined active K rotor and all hindered rotors
! in the manner of Astholz.
  call beyer(Total_harm_osc_freq, Delta_nu, E_max, rho)

WRITE(*,*) rho

END SUBROUTINE Calculate_rho
!---------------------------------------------------------------------------

END MODULE
!---------------------------------------------------------------------------

! PROGRAM main
! 
! USE DensityOfStates
! 
!   IMPLICIT NONE
! 
! ! These variables are used in the computation of the density of states and the
! ! Corresponding microcanonical rate equation
!   REAL(8), ALLOCATABLE, DIMENSION(:) :: rho
!   INTEGER :: Delta_nu  
!   INTEGER :: E_Max 
! 
!   ! These variables are used to open the file
!   CHARACTER(20) :: filename
!   
!   Delta_nu = 10 ! Make sure that E_max is divisible by Delta_nu!
!   E_max = 2E4      ! 3.5 cm^-1 = 100 kcal/mol
!   filename = 'rho_input'
! 
! 
! ! The file should have the following structure:
! ! Number of atoms
! ! number of rotors
! ! linearity (0 if linear, 1 otherwise)
! ! list of all the vibrational modes
! ! list of all the frequencies and barrier heights for the hindered rotors.
! ! It is important that everything add up:  
! ! 3*N - 6 - #frequencies - #hindered rotors = 0 (for nonlinear system).
!    
! ! Allocate the density of state:
! ! it should be equal to the total number of energy bins, which is equal to
! ! the maximum energy value divided by the bin size, plus one (for E = 0).
!   ALLOCATE( rho( (E_max/Delta_nu) + 1) )
! 
!   CALL Calculate_rho(filename, E_max, delta_nu, rho)
! 
!   STOP
! END PROGRAM main
! 

