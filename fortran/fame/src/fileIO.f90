! ==============================================================================
!
!	fileIO.f90
!
!	Written by Josh Allen (jwallen@mit.edu)
!
! ==============================================================================

module fileIO



contains

	! --------------------------------------------------------------------------

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

	! --------------------------------------------------------------------------

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

	! --------------------------------------------------------------------------

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

	! --------------------------------------------------------------------------

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

	! --------------------------------------------------------------------------

end module
