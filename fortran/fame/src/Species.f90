! ==============================================================================
!
!	Species.f90
!
!	Written by Josh Allen (jwallen@mit.edu)
!
! ==============================================================================

module SpeciesModule

	implicit none

	! Struct: Isomer
	! 
	! Contains the information for a (unimolecular or multimolecular) isomer well.
	type Species
		real(8)								:: 	E0				! Ground-state electronic + zero-point energy in kJ/mol
		real(8)								:: 	dHf				! Standard enthalpy of formation in kJ/mol
		real(8)								:: 	dGf				! Standard Gibbs free energy of formation in kJ/mol
		real(8)								:: 	MW				! Molecular weight in g/mol
		real(8)								::	sigma			! Lennard-Jones sigma parameter in m
		real(8)								:: 	eps				! Lennard-Jones epsilon paramter in J
		real(8), dimension(:), allocatable	:: 	vibFreq			! 1D harmonic oscillator frequencies in cm^-1
		real(8), dimension(:), allocatable	:: 	rotFreq			! 1D rigid rotor frequencies in cm^-1
		real(8), dimension(:), allocatable	:: 	hindFreq		! 1D hindered rotor frequencies in cm^-1
		real(8), dimension(:), allocatable	:: 	hindBarrier		! 1D hindered rotor barrier heights in cm^-1
		integer								:: 	symmNum			! Symmetry number
	end type	
	
end module
