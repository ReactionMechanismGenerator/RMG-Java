!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	isomer.f90
!
!	Written by Josh Allen (jwallen@mit.edu)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module IsomerModule

	implicit none
	
	type Isomer	
		integer, dimension(:), allocatable		::	species			! List of indices to species in the well
		real(8), dimension(:), allocatable		::	densStates		! Density of states in mol/J
		real(8)									::	E0				! Ground state energy in J/mol
		real(8)									::	Q				! Partition function
		real(8), dimension(:), allocatable		::	eqDist			! Equilibrium distribution
		real(8)									::	collFreq		! Collision frequency in Hz
		real(8), dimension(:,:), allocatable	::	Mcoll			! Collisional transfer probability matrix
	end type



end module
