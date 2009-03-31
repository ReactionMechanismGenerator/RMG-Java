!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	Input.f90
!
!	Written by Josh Allen (jwallen@mit.edu)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module InputModule

	use SimulationModule
	use SpeciesModule
	use IsomerModule
	use ReactionModule

	implicit none
	
contains

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!
	! Subroutine: loadNetwork()
	! 
	! Loads the parameters from the specified file needed to complete the
	! master equation simulation.
	!
	! Parameters:
	!		simData - The simulation data object to store the parameters in.
	subroutine loadNetwork(path, simData, speciesList, isomerList, rxnList, verbose)
		
		character(len=*), intent(in)								:: 	path
		type(Simulation), intent(inout)								:: 	simData
		type(Species), dimension(:), allocatable, intent(inout)		:: 	speciesList
		type(Isomer), dimension(:), allocatable, intent(inout)		:: 	isomerList
		type(Reaction), dimension(:), allocatable, intent(inout)	:: 	rxnList
		integer, intent(in)											:: 	verbose
		
		integer ios		! Input file status
		
		integer i, n	! Dummy indices
		
		! Open file for reading; fail if unable to open
		open(1, iostat=ios, file=path, action='read', status='old', access='sequential', recl=128)
		if (ios /= 0) then
			print *, 'Error: Unable to open file "', path, '".'
			return
		end if

		if (verbose >= 1) write (*,*), 'Reading from file "fame_input.txt"...'

		! Load simulation data from input file
		if (verbose >= 2) write (*,*), '\tReading simulation parameters...'
		call loadSimulationData(simData)
				
		! Load unimolecular isomer well data from input file
		if (verbose >= 2) write (*,*), '\tReading species data...'
		allocate (speciesList(1:simData%nSpecies))
		do i = 1, simData%nSpecies
			call loadSpeciesData(speciesList(i))
		end do

		! Load unimolecular isomer well data from input file
		if (verbose >= 2) write (*,*), '\tReading isomer data...'
		allocate (isomerList(1:simData%nIsom))
		do i = 1, simData%nIsom
			call loadIsomerData(isomerList(i), speciesList)
		end do

		! Load reaction data from input file
		if (verbose >= 2) write (*,*), '\tReading reaction data...'
		allocate (rxnList(1:simData%nRxn))
		do i = 1, simData%nRxn
			call loadReactionData(rxnList(i))
		end do
		
		! Close file
		close(1) 
		
	end subroutine
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
			if (index(str(1:4), 'Mode') /= 0) then
				call removeWhitespace(str(5:))
				if (index(str(5:27), 'ModifiedStrongCollision') /= 0) then
					simData%mode = 1
				elseif (index(str(5:18), 'ReservoirState') /= 0) then
					simData%mode = 2
				else
					simData%mode = 0
				end if
			else if (index(str(1:12), 'Temperatures') /= 0) then
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
			else if (index(str(1:20), 'Minimum grain energy') /= 0) then
				call readNumberAndUnit(str(21:), simData%Emin, unit)
			else if (index(str(1:20), 'Maximum grain energy') /= 0) then
                call readNumberAndUnit(str(21:), simData%Emax, unit)
            else if (index(str(1:15), 'Simulation time') /= 0) then
				call readNumberAndUnit(str(16:), simData%tfinal, unit)
			else if (index(str(1:17), 'Number of species') /= 0) then
				call readInteger(str(18:), simData%nSpecies)
			else if (index(str(1:15), 'Number of wells') /= 0) then
				call readInteger(str(16:), simData%nIsom)
			else if (index(str(1:19), 'Number of reactions') /= 0) then
				call readInteger(str(20:), simData%nRxn)
			else if (index(str(1:26), 'Exponential down parameter') /= 0) then
				call readNumberAndUnit(str(27:), simData%alpha, unit)
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
        simData%nGrains = nint((simData%Emax - simData%Emin) / simData%dE)
		allocate (simData%E(1:simData%nGrains))
        do i = 1, simData%nGrains
            simData%E(i) = simData%dE * (i - 0.5) + simData%Emin
        end do
		
! 		write (*,*) simData%mode
! 		write (*,*) simData%Tlist
! 		write (*,*) simData%Plist
! 		write (*,*) simData%dE
! 		write (*,*) simData%Emin
! 		write (*,*) simData%Emax
! 		write (*,*) simData%alpha
! 		write (*,*) simData%nSpecies
! 		write (*,*) simData%nIsom
! 		write (*,*) simData%nRxn
! 		write (*,*) simData%nGrains
! 		write (*,*) simData%nChebT
! 		write (*,*) simData%nChebP
! 		write (*,*) simData%bathGas%MW
! 		write (*,*) simData%bathGas%sigma
! 		write (*,*) simData%bathGas%eps
		
	end subroutine
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!
	! Subroutine: loadSpeciesData()
	! 
	! Loads data for a single species well into the specified object.
	! Assumes that the input file is already open and in the correct position for
	! reading.
	!
	! Parameters:
	!		spec - An Isomer object to load the data into.
	subroutine loadSpeciesData(spec)
		
		! Provide parameter type checking of inputs and outputs
		type(Species), intent(inout)		:: spec
		
		! Local variables
		integer found, ios, i
		character(len=256) str
		character(len=32) unit
		integer number
		
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
			if (index(str(1:19), 'Ground-state energy') /= 0) then
				call readNumberAndUnit(str(21:), spec%E0, unit)
			else if (index(str(1:21), 'Enthalpy of formation') /= 0) then
				call readNumberAndUnit(str(22:), spec%dHf, unit)
			else if (index(str(1:24), 'Free energy of formation') /= 0) then
				call readNumberAndUnit(str(25:), spec%dGf, unit)
			else if (index(str(1:16), 'Molecular weight') /= 0) then
				call readNumberAndUnit(str(25:), spec%MW, unit)
			else if (index(str(1:18), 'LJ sigma parameter') /= 0) then
				call readNumberAndUnit(str(25:), spec%sigma, unit)
			else if (index(str(1:20), 'LJ epsilon parameter') /= 0) then
				call readNumberAndUnit(str(25:), spec%eps, unit)
			else if (index(str(1:20), 'Harmonic oscillators') /= 0) then
				call readInteger(str(21:), number)
				allocate(spec%vibFreq(1:number))
				do i = 1, number
					read (1, fmt='(a128)', iostat=ios), str
					call readNumberAndUnit(str, spec%vibFreq(i), unit)
				end do
			else if (index(str(1:12), 'Rigid rotors') /= 0) then
				call readInteger(str(13:), number)
				allocate(spec%rotFreq(1:number))
				do i = 1, number
					read (1, fmt='(a128)', iostat=ios), str
					call readNumberAndUnit(str, spec%rotFreq(i), unit)
				end do
			else if (index(str(1:15), 'Hindered rotors') /= 0) then
				call readInteger(str(16:), number)
				allocate(spec%hindFreq(1:number))
				allocate(spec%hindBarrier(1:number))
				do i = 1, number
					read (1, fmt='(a128)', iostat=ios), str
					call readNumberAndUnit(str, spec%hindFreq(i), unit)
				end do
				do i = 1, number
					read (1, fmt='(a128)', iostat=ios), str
					call readNumberAndUnit(str, spec%hindBarrier(i), unit)
				end do
			else if (index(str(1:15), 'Symmetry number') /= 0) then
				call readInteger(str(21:), spec%symmNum)
			end if

		end do
		
! 		write (*,*) spec%E0
! 		write (*,*) spec%dHf
! 		write (*,*) spec%dGf
! 		write (*,*) spec%MW
! 		write (*,*) spec%sigma
! 		write (*,*) spec%eps
! 		write (*,*) spec%vibFreq
! 		write (*,*) spec%rotFreq
! 		write (*,*) spec%hindFreq
! 		write (*,*) spec%hindBarrier
! 		write (*,*) spec%symmNum
			
	end subroutine
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!
	! Subroutine: loadIsomerData()
	! 
	! Loads data for a single isomer well into the specified object.
	! Assumes that the input file is already open and in the correct position for
	! reading.
	!
	! Parameters:
	!		isomer - An Isomer object to load the data into.
	subroutine loadIsomerData(isom, speciesList)
		
		! Provide parameter type checking of inputs and outputs
		type(Isomer), intent(inout)					:: isom
		type(Species), dimension(:), intent(in)		:: speciesList
		
		! Local variables
		integer found, ios, i
		character(len=256) str
		character(len=32) unit
		integer number
		
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
			if (index(str(1:7), 'Species') /= 0) then
				call readInteger(str(8:), isom%numSpecies)
				allocate(isom%speciesList(1:isom%numSpecies))
				do i = 1, isom%numSpecies
					read (1, fmt='(a128)', iostat=ios), str
					call readInteger(str, isom%speciesList(i))
				end do
			end if

		end do
		
		isom%E0 = 0
		isom%dHf = 0
		isom%dGf = 0
		isom%MW = 0
		do i = 1, isom%numSpecies
			isom%E0 = isom%E0 + speciesList(isom%speciesList(i))%E0
			isom%dHf = isom%dHf + speciesList(isom%speciesList(i))%dHf
			isom%dGf = isom%dGf + speciesList(isom%speciesList(i))%dGf
			isom%MW = isom%MW + speciesList(isom%speciesList(i))%MW
		end do
		
		isom%omega = 0
		
! 		write (*,*) isom%numSpecies
! 		write (*,*) isom%speciesList
! 		write (*,*) isom%E0
! 		write (*,*) isom%dHf
! 		write (*,*) isom%dGf
! 		write (*,*) isom%MW
! 		write (*,*) isom%omega
		
	end subroutine
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!
	! Subroutine: loadReactionData()
	! 
	! Loads the data from the specified files into the specified isomerData 
	! objects.
	!
	! Parameters:
	!		rxnData - A reaction object to fill with data
	subroutine loadReactionData(rxnData)
		
		! Provide parameter type checking of inputs and outputs
		type(Reaction), intent(inout)		:: rxnData
		
		! Local variables
		integer found, ios
		character(len=128)  str
		character(len=32) 	unit
		integer				r
		real(8)				tempreal

		! Initialize data members
		rxnData%isomerList(1) = 0
		rxnData%isomerList(2) = 0
		rxnData%E0 = 0
		rxnData%arrh_A = 0
		rxnData%arrh_n = 0
		rxnData%arrh_Ea = 0
		
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
			if (index(str(1:19), 'Ground-state energy') /= 0) then
				call readNumberAndUnit(str(20:), rxnData%E0, unit)
			else if (index(str(1:24), 'Arrhenius preexponential') /= 0) then
					call readNumberAndUnit(str(25:), rxnData%arrh_A, unit)
			else if (index(str(1:27), 'Arrhenius activation energy') /= 0) then
					call readNumberAndUnit(str(28:), rxnData%arrh_Ea, unit)
			else if (index(str(1:30), 'Arrhenius temperature exponent') /= 0) then
				call readNumber(str(31:), rxnData%arrh_n)
			else if (index(str(1:15), 'Reactant isomer') /= 0) then
				call readInteger(str(16:), rxnData%isomerList(1))
			else if (index(str(1:14), 'Product isomer') /= 0) then
				call readInteger(str(15:), rxnData%isomerList(2))
			end if

		end do
		
! 		write (*,*) rxnData%isomerList
! 		write (*,*) rxnData%E0
! 		write (*,*) rxnData%arrh_A
! 		write (*,*) rxnData%arrh_n
! 		write (*,*) rxnData%arrh_Ea
		
	end subroutine
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!
	! Subroutine: readNumberAndUnit()
	! 
	! Loads the parameters from the specified file needed to complete the
	! master equation simulation.
	!
	! Parameters:
	!		str - A string containing a number followed by a unit, with any 
	!			number of whitespace before, after, or in between.
	!		number - The number extracted from str.
	!		unit - The unit extracted from str.
	subroutine readNumberAndUnit(str, number, unit)

		! Provide parameter type checking of inputs and outputs
		character(len=*), intent(inout)	:: 	str
		real(8), intent(out)			::	number
		character(len=*), intent(out) 	:: 	unit
		
		integer i

		call readNumber(str, number)

		! Remove leading whitespace
		call removeWhitespace(str)
		
		! Find length of unit
		i = index(str, ' ')
		if (index(str, '\t') /= 0 .and. (i == 0 .or. i > index(str, '\t'))) i = index(str, '\t')
		if (index(str, '\n') /= 0 .and. (i == 0 .or. i > index(str, '\n'))) i = index(str, '\n')
		if (index(str, '\r') /= 0 .and. (i == 0 .or. i > index(str, '\r'))) i = index(str, '\r')
		
		if (i <= 0) then
			unit = ''
			return
		end if

		! Read unit
		unit = str(1:i-1)
		str = str(i:)

	end subroutine

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!
	! Subroutine: readNumber()
	! 
	! Loads the parameters from the specified file needed to complete the
	! master equation simulation.
	!
	! Parameters:
	!		str - A string containing a number, with any number of whitespace 
	!			before or after.
	!		number - The number extracted from str.
	subroutine readNumber(str, number)

		! Provide parameter type checking of inputs and outputs
		character(len=*), intent(inout)		:: 	str
		double precision, intent(out)		::	number
		
		integer		:: i

		! Remove leading whitespace
		call removeWhitespace(str)
		
		! Find length of number
		i = index(str, ' ')
		if (index(str, '\t') /= 0 .and. (i == 0 .or. i > index(str, '\t'))) i = index(str, '\t')
		if (index(str, '\n') /= 0 .and. (i == 0 .or. i > index(str, '\n'))) i = index(str, '\n')
		if (index(str, '\r') /= 0 .and. (i == 0 .or. i > index(str, '\r'))) i = index(str, '\r')
		
		! Read number
		read(str(1:), *), number
		str = str(i:)
		
		
	end subroutine

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!
	! Subroutine: readInteger()
	! 
	! Loads the parameters from the specified file needed to complete the
	! master equation simulation.
	!
	! Parameters:
	!		str - A string containing a number, with any number of whitespace 
	!			before or after.
	!		number - The number extracted from str.
	subroutine readInteger(str, number)

		! Provide parameter type checking of inputs and outputs
		character(len=*), intent(inout)		:: 	str
		integer, intent(out)				::	number
		
		real(8) num
		
		call readNumber(str, num)
		
		! Convert to integer
		number = int(num)
		
	end subroutine

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!
	! Subroutine: removeWhitespace()
	! 
	! Removes leading whitespace from a string.
	!
	! Parameters:
	!		str - The string to remove the whitespace from.
	subroutine removeWhitespace(str)

		! Provide parameter type checking of inputs and outputs
		character(len=*), intent(inout)	:: 	str
		
		! Remove leading whitespace
		do while (index(str(1:1), '\t') /= 0 .or. index(str(1:1), ' ') /= 0 .or. &
			index(str(1:1), '\n') /= 0 .or. index(str(1:1), '\r') /= 0 .and. len(str) > 0)
			str = str(2:)
		end do

	end subroutine

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
end module
