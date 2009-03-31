! ==============================================================================
!
!	Simulation.f90
!
!	Written by Josh Allen (jwallen@mit.edu)
!
! ==============================================================================

module SimulationModule

	implicit none
	
	! Struct: BathGas
	! 
	! Contains information about the bath gas.
	type BathGas
		real(8)	MW				! The molecular weight of the bath gas in g/mol
		real(8)	sigma			! The Lennard-Jones sigma parameter of the bath gas in m
		real(8)	eps				! The Lennard-Jones epsilon parameter of the bath gas in J
	end type
	
	! Struct: Simulation
	! 
	! Contains the parameters necessary to complete the master equation 
	! calculation.
	type Simulation
		integer								:: mode		! 1 = modified strong collision, 2 = reservoir state
		real(8), dimension(:), allocatable	:: Tlist	! An array of temperatures to evaluate k(T, P) at
		real(8), dimension(:), allocatable	:: Plist	! An array of pressures to evaluate k(T, P) at
		real(8)	dE				! The grain size in kJ/mol
		real(8) Emin            ! The minimum energy in kJ/mol
        real(8)	Emax			! The maximum energy in kJ/mol
		real(8)	tfinal			! The length of the simulation in s
		real(8)	alpha			! The average energy loss in a deactivating collision in kJ/mol
		integer	nSpecies		! The number of species
		integer	nIsom			! The number of wells
		integer	nRxn			! The number of reactions
		type(BathGas) :: bathGas							! Data for the bath gas
		real(8), dimension(:), allocatable		:: E		! The array of energies in kJ/mol
		integer	nGrains			! The number of energy grains
		integer	nChebT			! The number of temperatures to use to fit Chebyshev polynomials
		integer	nChebP			! The number of pressures to use to fit Chebyshev polynomials
		
	end type

end module
