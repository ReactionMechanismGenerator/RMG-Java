! ==============================================================================
!
!	Species.f90
!
!	Written by Josh Allen (jwallen@mit.edu)
!
! ==============================================================================

module SpeciesModule

	use fileIO
	
	! Struct: UniWell
	! 
	! Contains the information for a unimolecular isomer well.
	type UniWell
		real(8)			:: 	E			! Ground-state electronic + zero-point energy in kJ/mol
		real(8)			:: 	H			! Standard enthalpy of formation in kJ/mol
		real(8)			:: 	G			! Standard Gibbs free energy of formation in kJ/mol
		real(8)			:: 	MW			! Molecular weight in g/mol
		real(8)			::	sigma		! Lennard-Jones sigma parameter in m
		real(8)			:: 	eps			! Lennard-Jones epsilon paramter in J
		real(8), dimension(:), allocatable	:: 	vibFreq			! 1D harmonic oscillator frequencies in cm^-1
		real(8), dimension(:), allocatable	:: 	rotFreq			! 1D rigid rotor frequencies in cm^-1
		real(8), dimension(:), allocatable	:: 	hindFreq		! 1D hindered rotor frequencies in cm^-1
		real(8), dimension(:), allocatable	:: 	hindBarrier		! 1D hindered rotor barrier heights in cm^-1
		integer								:: 	symmNum			! Symmetry number for each molecule
		real(8), dimension(:), allocatable	:: 	densStates		! Density of states
	end type
	
	! Struct: MultiWell
	! 
	! Contains the information for a multimolecular isomer source/sink.
	type MultiWell
		integer									:: 	numSpecies		! Number of species in well (usually 2 or 3)
		real(8), dimension(:), allocatable		:: 	E				! Ground-state electronic + zero-point energy in kJ/mol
		real(8), dimension(:), allocatable		:: 	H				! Standard enthalpy of formation in kJ/mol
		real(8), dimension(:), allocatable		:: 	G				! Standard Gibbs free energy of formation in kJ/mol
		real(8), dimension(:,:), allocatable	:: 	vibFreq			! 1D harmonic oscillator frequencies in cm^-1
		real(8), dimension(:,:), allocatable	:: 	rotFreq			! 1D rigid rotor frequencies in cm^-1
		real(8), dimension(:,:), allocatable	:: 	hindFreq		! 1D hindered rotor frequencies in cm^-1
		real(8), dimension(:,:), allocatable	:: 	hindBarrier		! 1D hindered rotor barrier heights in cm^-1
		integer, dimension(:), allocatable		:: 	symmNum			! Symmetry number for each molecule
		real(8), dimension(:), allocatable		:: 	densStates		! Convolved density of states for all species in well
	end type

contains

	! --------------------------------------------------------------------------
	!
	! Subroutine: calcEqDists()
	! 
	! Calculates the equilibrium distribution and partition function from
	! density-of-states data for an isomer.
	!
	! Parameters:
	!		uniData - The unimolecular well data.
	!		multiData - The multimolecular source/sink data.
	!		energies - The grain energies in kJ/mol.
	!		T - The absolute temperature in K.
	!		uniEqDist - The matrix of Boltzmann distributions for each unimolecular isomer.
	!		multiEqDist - The matrix of Boltzmann distributions for each multimolecular isomer.
	subroutine calcEqDists(uniData, multiData, energies, T, uniEqDist, multiEqDist)
		
		type(UniWell), dimension(:), intent(in)		:: uniData
		type(MultiWell), dimension(:), intent(in)	:: multiData
		real(8), dimension(:), intent(in)			:: energies
		real(8), intent(in)							:: T
		real(8), dimension(:,:), intent(out)		:: uniEqDist
		real(8), dimension(:,:), intent(out)		:: multiEqDist
		
		integer nGrains, nUni, nMulti
		real(8)	Q					! Used for normalization of distribution
		real(8)	R 					! Gas constant in J mol^-1 K^-1
		integer	i, n, s				! Dummy indices

		nGrains = size(energies)
		nUni = size(uniData)
		nMulti = size(multiData)

		R = 8.314472

		! Calculate equilibrium distributions for unimolecular wells
		do i = 1, nUni
			
			! Calculate unnormalized eqDist
			do s = 1, nGrains
				uniEqDist(s,i) = uniData(i)%densStates(s) * exp(-energies(s) * 1000 / (R * T))
			end do
			
			! Normalize eqDist
			Q = sum(uniEqDist(1:nGrains,i))			
			uniEqDist(1:nGrains,i) = uniEqDist(1:nGrains,i) / Q

		end do
		
		! Calculate equilibrium distributions for multimolecular sources/sinks
		do n = 1, nMulti
			
			! Calculate unnormalized eqDist
			do s = 1, nGrains
				multiEqDist(s,n) = multiData(n)%densStates(s) * exp(-energies(s) * 1000 / (R * T))
			end do
			
			! Normalize eqDist
			Q = sum(multiEqDist(1:nGrains,n))			
			multiEqDist(1:nGrains,n) = multiEqDist(1:nGrains,n) / Q

		end do

	end subroutine
	
	! --------------------------------------------------------------------------
	!
	! Subroutine: loadUniWellData()
	! 
	! Loads data for a single unimolecular isomer well into the specified object.
	! Assumes that the input file is already open and in the correct position for
	! reading.
	!
	! Parameters:
	!		uniData - A UniWell object to load the data into.
	!		nGrains - The number of grains in each isomer well.
	subroutine loadUniWellData(uniData, nGrains)
		
		! Provide parameter type checking of inputs and outputs
		type(UniWell), intent(inout)		:: uniData
		integer, intent(in) 				:: nGrains
		
		! Local variables
		integer found, ios
		character(len=256) str
		character(len=32) unit
		
		! Initialize data members
		uniData%MW = 0
		uniData%E = 0
		uniData%H = 0
		uniData%G = 0
		uniData%sigma = 0
		uniData%eps = 0
		uniData%symmNum = 1
		allocate (uniData%densStates(1:nGrains))
		
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
			if (index(str(1:16), 'Molecular weight') /= 0) then
				call readNumberAndUnit(str(17:), uniData%MW, unit)
			else if (index(str(1:18), 'LJ sigma parameter') /= 0) then
				call readNumberAndUnit(str(19:), uniData%sigma, unit)
			else if (index(str(1:20), 'LJ epsilon parameter') /= 0) then
				call readNumberAndUnit(str(21:), uniData%eps, unit)
			else if (index(str(1:19), 'Ground-state energy') /= 0) then
				call readNumberAndUnit(str(20:), uniData%E, unit)
			else if (index(str(1:21), 'Enthalpy of formation') /= 0) then
				call readNumberAndUnit(str(22:), uniData%H, unit)
			else if (index(str(1:24), 'Free energy of formation') /= 0) then
				call readNumberAndUnit(str(25:), uniData%G, unit)
			else if (index(str(1:20), 'Harmonic oscillators') /= 0) then
				call readInteger(str(21:), j)
				allocate(uniData%vibFreq(1:j))
				do i = 1, j
					read (1, fmt='(a128)', iostat=ios), str
					call readNumberAndUnit(str, uniData%vibFreq(i), unit)
				end do
			else if (index(str(1:12), 'Rigid rotors') /= 0) then
				call readInteger(str(13:), j)
				if (j > 0) then
					allocate(uniData%rotFreq(1:j))
					do i = 1, j
						read (1, fmt='(a128)', iostat=ios), str
						call readNumberAndUnit(str, uniData%rotFreq(i), unit)
					end do
				end if
			else if (index(str(1:15), 'Hindered rotors') /= 0) then
				call readInteger(str(16:), j)
				if (j > 0) then
					allocate(uniData%hindFreq(1:j))
					allocate(uniData%hindBarrier(1:j))
					do i = 1, j
						read (1, fmt='(a128)', iostat=ios), str
						call readNumberAndUnit(str, uniData%hindFreq(i), unit)
					end do
					do i = 1, j
						read (1, fmt='(a128)', iostat=ios), str
						call readNumberAndUnit(str, uniData%hindBarrier(i), unit)
					end do
				end if
			else if (index(str(1:15), 'Symmetry number') /= 0) then
				call readInteger(str(16:), uniData%symmNum)
			end if

		end do
		
		if (.not. allocated(uniData%rotFreq)) then
			allocate(uniData%rotFreq(1))
			uniData%rotFreq(1) = 0
		end if
		if (.not. allocated(uniData%hindFreq)) then
			allocate(uniData%hindFreq(1))
			uniData%hindFreq(1) = 0
			allocate(uniData%hindBarrier(1))
			uniData%hindBarrier(1) = 0
		end if
		
	end subroutine
	
	! --------------------------------------------------------------------------
	
	! Subroutine: loadMultiWellData()
	! 
	! Loads data for a single multimolecular isomer well into the specified object.
	! Assumes that the input file is already open and in the correct position for
	! reading.
	!
	! Parameters:
	!		multiData - A MultiWell object to load the data into.
	!		nGrains - The number of grains in each isomer well.
	subroutine loadMultiWellData(multiData, nGrains)
		
		! Provide parameter type checking of inputs and outputs
		type(MultiWell), intent(inout)		:: multiData
		integer, intent(in) 				:: nGrains
		
		! Local variables
		integer found, ios
		character(len=256) str
		character(len=32) unit
		integer, dimension(:), allocatable :: number
		
		! Initialize data members
		allocate (multiData%densStates(1:nGrains))
	
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
			if (index(str(1:17), 'Number of species') /= 0) then
				call readInteger(str(18:), multiData%numSpecies)
			else if (index(str(1:19), 'Ground-state energy') /= 0) then
				allocate(multiData%E(1:multiData%numSpecies))
				do i = 1, multiData%numSpecies
					call readNumberAndUnit(str(21:), multiData%E(i), unit)
				end do
			else if (index(str(1:21), 'Enthalpy of formation') /= 0) then
				allocate(multiData%H(1:multiData%numSpecies))
				do i = 1, multiData%numSpecies
					call readNumberAndUnit(str(22:), multiData%H(i), unit)
				end do
			else if (index(str(1:24), 'Free energy of formation') /= 0) then
				allocate(multiData%G(1:multiData%numSpecies))
				do i = 1, multiData%numSpecies
					call readNumberAndUnit(str(25:), multiData%G(i), unit)
				end do
			else if (index(str(1:20), 'Harmonic oscillators') /= 0) then
				allocate(number(1:multiData%numSpecies))
				do i = 1, multiData%numSpecies
					call readInteger(str(21:), number(i))
				end do
				allocate(multiData%vibFreq(1:multiData%numSpecies, 1:maxval(number)))
				multiData%vibFreq = 0 * multiData%vibFreq
				do i = 1, multiData%numSpecies
					do j = 1, number(i)
						read (1, fmt='(a128)', iostat=ios), str
						call readNumberAndUnit(str, multiData%vibFreq(i, j), unit)
					end do
				end do
				deallocate(number)
			else if (index(str(1:12), 'Rigid rotors') /= 0) then
				allocate(number(1:multiData%numSpecies))
				do i = 1, multiData%numSpecies
					call readInteger(str(21:), number(i))
				end do
				if (maxval(number) > 0) then
					allocate(multiData%rotFreq(1:multiData%numSpecies, 1:maxval(number)))
					multiData%rotFreq = 0 * multiData%rotFreq
					do i = 1, multiData%numSpecies
						do j = 1, number(i)
							read (1, fmt='(a128)', iostat=ios), str
							call readNumberAndUnit(str, multiData%rotFreq(i, j), unit)
						end do
					end do
				end if
				deallocate(number)
			else if (index(str(1:15), 'Hindered rotors') /= 0) then
				allocate(number(1:multiData%numSpecies))
				do i = 1, multiData%numSpecies
					call readInteger(str(16:), number(i))
				end do
				if (maxval(number) > 0) then
					allocate(multiData%hindFreq(1:multiData%numSpecies, 1:maxval(number)))
					allocate(multiData%hindBarrier(1:multiData%numSpecies, 1:maxval(number)))
					multiData%hindFreq = 0 * multiData%hindFreq
					multiData%hindBarrier = 0 * multiData%hindBarrier
					do i = 1, multiData%numSpecies
						do j = 1, number(i)
							read (1, fmt='(a128)', iostat=ios), str
							call readNumberAndUnit(str, multiData%hindFreq(i, j), unit)
						end do
						do j = 1, number(i)
							read (1, fmt='(a128)', iostat=ios), str
							call readNumberAndUnit(str, multiData%hindBarrier(i, j), unit)
						end do
					end do
				end if
				deallocate(number)
			else if (index(str(1:15), 'Symmetry number') /= 0) then
				allocate(multiData%symmNum(1:multiData%numSpecies))
				do i = 1, multiData%numSpecies
					call readInteger(str(21:), multiData%symmNum(i))
				end do
			end if

		end do
		
		if (.not. allocated(multiData%rotFreq)) then
			allocate(multiData%rotFreq(1:multiData%numSpecies, 1))
			multiData%rotFreq = 0 * multiData%rotFreq
		end if
		if (.not. allocated(multiData%hindFreq)) then
			allocate(multiData%hindFreq(1:multiData%numSpecies, 1))
			multiData%hindFreq = 0 * multiData%hindFreq
			allocate(multiData%hindBarrier(1:multiData%numSpecies, 1))
			multiData%hindBarrier = 0 * multiData%hindBarrier
		end if
		
	end subroutine
	
	! --------------------------------------------------------------------------
	!
	! Subroutine: writeIsomerData()
	! 
	! Writes the loaded isomer data to the screen (for debugging).
	!
	! Parameters:
	!		uniData - The unimolecular isomer wells.
	!		multiData - The multimolecular isomer sources/sinks.
	subroutine writeIsomerData(uniData, multiData)
	
		type(UniWell), dimension(:), intent(in)		:: uniData
		type(MultiWell), dimension(:), intent(in)	:: multiData
		
		integer i, n
		
		do i = 1, size(uniData)
			write (*,*), 'Unimolecular well', i
			write (*,*), '\tElectronic energy:', uniData(i)%E, 'kJ/mol'
			write (*,*), '\tEnthalpy of formation:', uniData(i)%H, 'kJ/mol'
			write (*,*), '\tFree energy of formation:', uniData(i)%G, 'kJ/mol'
			write (*,*), '\tMolecular weight:', uniData(i)%MW, 'g/mol'
			write (*,*), '\tLJ sigma parameter:', uniData(i)%sigma, 'm'
			write (*,*), '\tLJ epsilon parameter:', uniData(i)%eps, 'J'
			if (size(uniData(i)%vibFreq) > 0) then
				write (*,*), '\tHarmonic oscillators:', uniData(i)%vibFreq, 'cm^-1'
			end if
			if (size(uniData(i)%rotFreq) > 0) then
				write (*,*), '\tRigid rotors:', uniData(i)%rotFreq, 'cm^-1'
			end if
			if (size(uniData(i)%hindFreq) > 0) then
				write (*,*), '\tHindered rotors:', uniData(i)%hindFreq, 'cm^-1', uniData(i)%hindBarrier, 'cm^-1'
			end if
			write (*,*), '\tSymmetry number:', uniData(i)%symmNum
		end do
		
		do n = 1, size(multiData)
			write (*,*), 'Multimolecular source/sink', n
			write (*,*), '\tElectronic energy:', multiData(n)%E, 'kJ/mol'
			write (*,*), '\tEnthalpy of formation:', multiData(n)%H, 'kJ/mol'
			write (*,*), '\tFree energy of formation:', multiData(n)%G, 'kJ/mol'
			write (*,*), '\tNumber of species:', multiData(n)%numSpecies
			do i = 1, multiData(n)%numSpecies
				if (size(multiData(n)%vibFreq(i,:)) > 0) then
					write (*,*), '\tSpecies', i, 'Harmonic oscillators:', multiData(n)%vibFreq(i,:), 'cm^-1'
				end if
				if (size(multiData(n)%rotFreq(i,:)) > 0) then
					write (*,*), '\tSpecies', i, 'Rigid rotors:', multiData(n)%rotFreq(i,:), 'cm^-1'
				end if
				if (size(multiData(n)%hindFreq(i,:)) > 0) then
					write (*,*), '\tSpecies', i, 'Hindered rotors:', multiData(n)%hindFreq(i,:), 'cm^-1', uniData(i)%hindBarrier, 'cm^-1'
				end if
				write (*,*), '\tSpecies', i, 'Symmetry number:', multiData(n)%symmNum(i)
			end do
		end do
		
	end subroutine
	
	! --------------------------------------------------------------------------
	
end module
