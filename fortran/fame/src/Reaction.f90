! ==============================================================================
!
!	Reaction.f90
!
!	Written by Josh Allen (jwallen@mit.edu)
!
! ==============================================================================

module ReactionModule

	implicit none
	
	! Struct: Reaction
	! 
	! Contains the information for a single reaction.
	type Reaction
		integer			::	isomer(2)	! Isomer wells connected to transition state
		real(8)			:: 	E			! Ground-state electronic + zero-point energy of transition state in kJ/mol
		real(8)			::	arrh_A		! Arrhenius preexponential factor in s^-1
		real(8)			::	arrh_n		! Arrhenius temperature exponent
		real(8)			::	arrh_Ea		! Arrhenius activation energy in kJ/mol
		! real(8), dimension(:), allocatable	:: 	sumStates			! Sum of states
	end type

contains
	
	! --------------------------------------------------------------------------
	!
	! Subroutine: writeReactionData()
	! 
	! Writes the loaded reaction data to the screen (for debugging).
	!
	! Parameters:
	!		uniData - The unimolecular isomer wells.
	!		multiData - The multimolecular isomer sources/sinks.
	subroutine writeReactionData(rxnData)
	
		type(Reaction), dimension(:), intent(in)	:: rxnData
		
		integer r
		
		do r = 1, size(rxnData)
			write (*,*), 'Reaction', r
			write (*,*), '\tIsomers:', rxnData(r)%isomer(1), 'and', rxnData(r)%isomer(2)
			write (*,*), '\tTS electronic energy:', rxnData(r)%E, 'kJ/mol'
			write (*,*), '\tArrhenius preexponential factor:', rxnData(r)%arrh_A, 's^-1'
			write (*,*), '\tArrhenius temperature exponent:', rxnData(r)%arrh_n
			write (*,*), '\tArrhenius activation energy:', rxnData(r)%arrh_Ea, 'kJ/mol'
		end do
		
	end subroutine
	
	! --------------------------------------------------------------------------
	
end module
