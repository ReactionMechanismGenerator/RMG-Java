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
		real(8)	T				! The currently-active temperature of the simulation in K
		real(8)	P				! The currently-active pressure of the simulation in bar
		real(8)	dE				! The grain size in kJ/mol
		real(8) Emin            ! The minimum energy in kJ/mol
        real(8)	Emax			! The maximum energy in kJ/mol
		real(8)	tfinal			! The length of the simulation in s
		real(8)	alpha			! The average energy loss in a deactivating collision in kJ/mol
		integer	nUni			! The number of unimolecular isomer wells
		integer	nMulti			! The number of bimolecular source/sink wells
		integer	nRxn			! The number of reactions
		type(BathGas) :: bathGas							! Data for the bath gas
		real(8), dimension(:), allocatable		:: E		! The array of energies in kJ/mol
		integer	nGrains			! The number of energy grains
		integer	nChebT			! The number of temperatures to use to fit Chebyshev polynomials
		integer	nChebP			! The number of pressures to use to fit Chebyshev polynomials
		
	end type

contains

	! --------------------------------------------------------------------------
	!
	! Subroutine: writeSimulationData()
	! 
	! Writes the loaded simulation data to the screen (for debugging).
	!
	! Parameters:
	!		simData - The simulation parameters.
	subroutine writeSimulationData(simData)
		
		! Provide parameter type checking of inputs and outputs
		type(Simulation), intent(in)		:: 	simData

		write (*,*), 'Simulation parameters'
		write (*,*), '\tTemperatures:', simData%Tlist, 'K'
		write (*,*), '\tPressures:', simData%Plist, 'bar'
		write (*,*), '\tGrain size:', simData%dE, 'kJ/mol'
		write (*,*), '\tMaximum grain energy:', simData%Emax, 'kJ/mol'
		write (*,*), '\tExponential down parameter:', simData%alpha, 'kJ/mol'
		write (*,*), '\tBath gas LJ sigma parameter:', simData%bathGas%sigma, 'm'
		write (*,*), '\tBath gas LJ epsilon parameter:', simData%bathGas%eps, 'J'
		write (*,*), '\tBath gas molecular weight:', simData%bathGas%MW, 'g/mol'
		write (*,*), '\tNumber of unimolecular wells:', simData%nUni
		write (*,*), '\tNumber of multimolecular wells:', simData%nMulti
		write (*,*), '\tNumber of reactions:', simData%nRxn
		write (*,*), '\tNumber of grains:', simData%nGrains
		write (*,*), '\tNumber of Chebyshev temperatures:', simData%nChebT
		write (*,*), '\tNumber of Chebyshev pressures:', simData%nChebP
		
	end subroutine

	! --------------------------------------------------------------------------

end module
