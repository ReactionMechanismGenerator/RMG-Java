! ==============================================================================
!
!	Isomer.f90
!
!	Written by Josh Allen (jwallen@mit.edu)
!
! ==============================================================================

module IsomerModule

	use fileIO
	
	! Struct: Isomer
	! 
	! Contains the information for a (unimolecular or multimolecular) isomer well.
	type Isomer
		integer									:: 	numSpecies		! Number of species in well (usually 2 or 3)
		real(8), dimension(:), allocatable		:: 	E				! Ground-state electronic + zero-point energy in kJ/mol
		real(8), dimension(:), allocatable		:: 	H				! Standard enthalpy of formation in kJ/mol
		real(8), dimension(:), allocatable		:: 	G				! Standard Gibbs free energy of formation in kJ/mol
		real(8), dimension(:), allocatable		:: 	MW			! Molecular weight in g/mol
		real(8), dimension(:), allocatable		::	sigma		! Lennard-Jones sigma parameter in m
		real(8), dimension(:), allocatable		:: 	eps			! Lennard-Jones epsilon paramter in J
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
	! Subroutine: calcEqDist()
	! 
	! Calculates the equilibrium distribution and partition function from
	! density-of-states data for an isomer.
	!
	! Parameters:
	!		isom - The isom to calculate the distribution for.
	!		E - The grain energies in kJ/mol.
	!		T - The absolute temperature in K.
	!		eqDist - The vector in which to place the equilibrium distribution.
	subroutine calcEqDist(isom, E, T, eqDist)
		
		type(Isomer), intent(in)				:: isom
		real(8), dimension(:), intent(in)		:: E
		real(8), intent(in)						:: T
		real(8), dimension(:), intent(out)		:: eqDist
		
		integer nGrains
		real(8)	Q					! Used for normalization of distribution
		real(8)	R 					! Gas constant in J mol^-1 K^-1
		integer	i, n, s				! Dummy indices

		nGrains = size(E)
		
		R = 8.314472

		! Calculate unnormalized eqDist
		do s = 1, nGrains
			eqDist(s) = isom%densStates(s) * exp(-E(s) * 1000 / (R * T))
		end do
		
		! Normalize eqDist
		Q = sum(eqDist)			
		eqDist = eqDist / Q

	end subroutine
	
	! Subroutine: loadIsomerData()
	! 
	! Loads data for a single isomer well into the specified object.
	! Assumes that the input file is already open and in the correct position for
	! reading.
	!
	! Parameters:
	!		isom - An Isomer object to load the data into.
	!		nGrains - The number of grains in each isomer well.
	subroutine loadIsomerData(isom, nGrains)
		
		! Provide parameter type checking of inputs and outputs
		type(Isomer), intent(inout)		:: isom
		integer, intent(in) 			:: nGrains
		
		! Local variables
		integer found, ios
		character(len=256) str
		character(len=32) unit
		integer, dimension(:), allocatable :: number
		
		! Initialize data members
		allocate (isom%densStates(1:nGrains))
		isom%numSpecies = 1
		
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
				call readInteger(str(18:), isom%numSpecies)
			else if (index(str(1:19), 'Ground-state energy') /= 0) then
				allocate(isom%E(1:isom%numSpecies))
				do i = 1, isom%numSpecies
					call readNumberAndUnit(str(21:), isom%E(i), unit)
				end do
			else if (index(str(1:21), 'Enthalpy of formation') /= 0) then
				allocate(isom%H(1:isom%numSpecies))
				do i = 1, isom%numSpecies
					call readNumberAndUnit(str(22:), isom%H(i), unit)
				end do
			else if (index(str(1:24), 'Free energy of formation') /= 0) then
				allocate(isom%G(1:isom%numSpecies))
				do i = 1, isom%numSpecies
					call readNumberAndUnit(str(25:), isom%G(i), unit)
				end do
			else if (index(str(1:16), 'Molecular weight') /= 0) then
				allocate(isom%MW(1:isom%numSpecies))
				do i = 1, isom%numSpecies
					call readNumberAndUnit(str(25:), isom%MW(i), unit)
				end do
			else if (index(str(1:18), 'LJ sigma parameter') /= 0) then
				allocate(isom%sigma(1:isom%numSpecies))
				do i = 1, isom%numSpecies
					call readNumberAndUnit(str(25:), isom%sigma(i), unit)
				end do
			else if (index(str(1:20), 'LJ epsilon parameter') /= 0) then
				allocate(isom%eps(1:isom%numSpecies))
				do i = 1, isom%numSpecies
					call readNumberAndUnit(str(25:), isom%eps(i), unit)
				end do
			else if (index(str(1:20), 'Harmonic oscillators') /= 0) then
				allocate(number(1:isom%numSpecies))
				do i = 1, isom%numSpecies
					call readInteger(str(21:), number(i))
				end do
				allocate(isom%vibFreq(1:isom%numSpecies, 1:maxval(number)))
				isom%vibFreq = 0 * isom%vibFreq
				do i = 1, isom%numSpecies
					do j = 1, number(i)
						read (1, fmt='(a128)', iostat=ios), str
						call readNumberAndUnit(str, isom%vibFreq(i, j), unit)
					end do
				end do
				deallocate(number)
			else if (index(str(1:12), 'Rigid rotors') /= 0) then
				allocate(number(1:isom%numSpecies))
				do i = 1, isom%numSpecies
					call readInteger(str(21:), number(i))
				end do
				if (maxval(number) > 0) then
					allocate(isom%rotFreq(1:isom%numSpecies, 1:maxval(number)))
					isom%rotFreq = 0 * isom%rotFreq
					do i = 1, isom%numSpecies
						do j = 1, number(i)
							read (1, fmt='(a128)', iostat=ios), str
							call readNumberAndUnit(str, isom%rotFreq(i, j), unit)
						end do
					end do
				end if
				deallocate(number)
			else if (index(str(1:15), 'Hindered rotors') /= 0) then
				allocate(number(1:isom%numSpecies))
				do i = 1, isom%numSpecies
					call readInteger(str(16:), number(i))
				end do
				if (maxval(number) > 0) then
					allocate(isom%hindFreq(1:isom%numSpecies, 1:maxval(number)))
					allocate(isom%hindBarrier(1:isom%numSpecies, 1:maxval(number)))
					isom%hindFreq = 0 * isom%hindFreq
					isom%hindBarrier = 0 * isom%hindBarrier
					do i = 1, isom%numSpecies
						do j = 1, number(i)
							read (1, fmt='(a128)', iostat=ios), str
							call readNumberAndUnit(str, isom%hindFreq(i, j), unit)
						end do
						do j = 1, number(i)
							read (1, fmt='(a128)', iostat=ios), str
							call readNumberAndUnit(str, isom%hindBarrier(i, j), unit)
						end do
					end do
				end if
				deallocate(number)
			else if (index(str(1:15), 'Symmetry number') /= 0) then
				allocate(isom%symmNum(1:isom%numSpecies))
				do i = 1, isom%numSpecies
					call readInteger(str(21:), isom%symmNum(i))
				end do
			end if

		end do
		
		if (.not. allocated(isom%rotFreq)) then
			allocate(isom%rotFreq(1:isom%numSpecies, 1))
			isom%rotFreq = 0 * isom%rotFreq
		end if
		if (.not. allocated(isom%hindFreq)) then
			allocate(isom%hindFreq(1:isom%numSpecies, 1))
			isom%hindFreq = 0 * isom%hindFreq
			allocate(isom%hindBarrier(1:isom%numSpecies, 1))
			isom%hindBarrier = 0 * isom%hindBarrier
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
	
		type(Isomer), dimension(:), intent(in)		:: uniData
		type(Isomer), dimension(:), intent(in)		:: multiData
		
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
