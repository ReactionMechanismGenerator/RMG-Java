! ==============================================================================
!
!	Isomer.f90
!
!	Written by Josh Allen (jwallen@mit.edu)
!
! ==============================================================================

module IsomerModule

	implicit none

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
		real(8)									::	Q				! Partition function at current temperature
	end type

contains

	! --------------------------------------------------------------------------
	!
	! Function: eqDist()
	! 
	! Calculates the equilibrium distribution and partition function from
	! density-of-states data for an isomer.
	!
	! Parameters:
	!		isom - The isom to calculate the distribution for.
	!		E - The grain energies in kJ/mol.
	!		T - The absolute temperature in K.
	!		eqDist - The vector in which to place the equilibrium distribution.
	subroutine eqDist(isom, E, T, f, Q)
		
		type(Isomer), intent(in)				:: isom
		real(8), dimension(:), intent(in)		:: E
		real(8), intent(in)						:: T
		real(8), dimension(:), intent(out)		:: f
		real(8), intent(out)					:: Q
		
		integer nGrains
		real(8)	R 					! Gas constant in J mol^-1 K^-1
		integer	s					! Dummy index
		
		nGrains = size(E)
		
		R = 8.314472 / 1000.

		! Calculate unnormalized eqDist
		do s = 1, size(E)
			f(s) = isom%densStates(s) * exp(-E(s) / (R * T))
		end do
		
		! Normalize eqDist
		Q = sum(f)			
		f = f / Q
		
	end subroutine
	
	function eqRatio(dGf, dHf, T)

		real(8), intent(in)			:: dGf
		real(8), intent(in)			:: dHf
		real(8), intent(in)			:: T
		real(8)						:: eqRatio
		
		real(8) R, T0
		
		R = 8.314472 / 1000
		T0 = 298

		eqRatio = exp(- dGf / (R * T0) - dHf / R * (1. / T - 1. / T0))

	end function
	
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
