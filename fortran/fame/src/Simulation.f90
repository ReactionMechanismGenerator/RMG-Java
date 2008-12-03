! ==============================================================================
!
!	Simulation.f90
!
!	Written by Josh Allen (jwallen@mit.edu)
!
! ==============================================================================

module SimulationModule

	use fileIO

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
		real(8), dimension(:), allocatable	:: Tlist	! An array of temperatures to evaluate k(T, P) at
		real(8), dimension(:), allocatable	:: Plist	! An array of pressures to evaluate k(T, P) at
		real(8)	T				! The currently-active temperature of the simulation in K
		real(8)	P				! The currently-active pressure of the simulation in bar
		real(8)	dE				! The grain size in kJ/mol
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
	! Subroutine: loadSimulationData()
	! 
	! Loads the parameters from the specified file needed to complete the
	! master equation simulation. This function assumes that the input file
	! is already open and set to the appropriate read location.
	!
	! Parameters:
	!		simData - The simulation data object to store the parameters in.
	subroutine loadSimulationData(simData)
		
		! Provide parameter type checking of inputs and outputs
		type(Simulation), intent(inout)		:: 	simData

		character(len=128)  str
		character(len=32) 	unit
		real(8) 			tempreal
		integer				i
		integer				found, ios
		
		ios = 0
		found = 0
		do while (ios == 0 .and. found == 0)
			
			! Read one line from the file
			read (1, fmt='(a128)', iostat=ios), str
			
			! Skip if comment line
			if (index(str(1:1), '#') /= 0) cycle

			! Break if blank line
			if (len(trim(str)) == 0) then
				found = 1
				cycle
			end if

			! Process if not a comment and if a recognized entry
			if (index(str(1:12), 'Temperatures') /= 0) then
				call readInteger(str(13:), i)
				allocate(simData%Tlist(1:i))
				do i = 1, i
					read (1, fmt='(a128)', iostat=ios), str
					call readNumberAndUnit(str, simData%Tlist(i), unit)
				end do
			else if (index(str(1:9), 'Pressures') /= 0) then
				call readInteger(str(10:), i)
				allocate(simData%Plist(1:i))
				do i = 1, i
					read (1, fmt='(a128)', iostat=ios), str
					call readNumberAndUnit(str, simData%Plist(i), unit)
				end do
			else if (index(str(1:10), 'Grain size') /= 0) then
				call readNumberAndUnit(str(11:), simData%dE, unit)
			else if (index(str(1:20), 'Maximum grain energy') /= 0) then
				call readNumberAndUnit(str(21:), simData%Emax, unit)
			else if (index(str(1:15), 'Simulation time') /= 0) then
				call readNumberAndUnit(str(16:), simData%tfinal, unit)
			else if (index(str(1:26), 'Exponential down parameter') /= 0) then
				call readNumberAndUnit(str(27:), simData%alpha, unit)
			else if (index(str(1:28), 'Number of unimolecular wells') /= 0) then
				call readInteger(str(29:), simData%nUni)
			else if (index(str(1:30), 'Number of multimolecular wells') /= 0) then
				call readInteger(str(31:), simData%nMulti)
			else if (index(str(1:19), 'Number of reactions') /= 0) then
				call readInteger(str(20:), simData%nRxn)
			else if (index(str(1:27), 'Bath gas LJ sigma parameter') /= 0) then
				call readNumberAndUnit(str(28:), simData%bathGas%sigma, unit)
			else if (index(str(1:29), 'Bath gas LJ epsilon parameter') /= 0) then
				call readNumberAndUnit(str(30:), simData%bathGas%eps, unit)
			else if (index(str(1:25), 'Bath gas molecular weight') /= 0) then
				call readNumberAndUnit(str(26:), simData%bathGas%MW, unit)
			else if (index(str(1:32), 'Number of Chebyshev temperatures') /= 0) then
				call readInteger(str(33:), simData%nChebT)
			else if (index(str(1:29), 'Number of Chebyshev pressures') /= 0) then
				call readInteger(str(30:), simData%nChebP)
			end if

		end do
		
		! Sets the energy levels of the individual grains
		simData%nGrains = simData%Emax / simData%dE + 1
		allocate (simData%E(1:simData%nGrains))
		do i = 1, simData%nGrains
 			simData%E(i) = simData%dE * (i - 1)
 		end do	
		
	end subroutine

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
