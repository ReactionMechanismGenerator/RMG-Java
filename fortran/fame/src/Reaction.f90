! ==============================================================================
!
!	Reaction.f90
!
!	Written by Josh Allen (jwallen@mit.edu)
!
! ==============================================================================

module ReactionModule

	use fileIO

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
	! Subroutine: loadReactionData()
	! 
	! Loads the data from the specified files into the specified isomerData 
	! objects.
	!
	! Parameters:
	!		path - The path to a file containing the reaction data.
	!		rxnData - A vector of objects to load the isomer data into.
	!		iCount - The number of reactions.
	!		nGrains - The number of grains in each isomer well.
	subroutine loadReactionData(rxnData, nGrains)
		
		! Provide parameter type checking of inputs and outputs
		type(Reaction), intent(inout)		:: rxnData
		integer, intent(in) 				:: nGrains
		
		! Local variables
		integer found, ios
		character(len=128)  str
		character(len=32) 	unit
		integer				r
		real(8)				tempreal

		! Initialize data members
		rxnData%isomer(1) = 0
		rxnData%isomer(2) = 0
		rxnData%E = 0
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
				call readNumberAndUnit(str(20:), rxnData%E, unit)
			else if (index(str(1:24), 'Arrhenius preexponential') /= 0) then
					call readNumberAndUnit(str(25:), rxnData%arrh_A, unit)
			else if (index(str(1:27), 'Arrhenius activation energy') /= 0) then
					call readNumberAndUnit(str(28:), rxnData%arrh_Ea, unit)
			else if (index(str(1:30), 'Arrhenius temperature exponent') /= 0) then
				call readNumber(str(31:), rxnData%arrh_n)
			else if (index(str(1:8), 'Isomer 1') /= 0) then
				call readInteger(str(9:), rxnData%isomer(1))
			else if (index(str(1:8), 'Isomer 2') /= 0) then
				call readInteger(str(9:), rxnData%isomer(2))
			end if

		end do
		
	end subroutine
	
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
