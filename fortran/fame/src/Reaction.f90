! ==============================================================================
!
!	Reaction.f90
!
!	Written by Josh Allen (jwallen@mit.edu)
!
! ==============================================================================

module ReactionModule

	use IsomerModule

	implicit none
	
	! Struct: ArrheniusKinetics
	! 
	! Contains modified Arrhenius kinetics for a single reaction.
	type ArrheniusKinetics
		real(8)			::	A			! Arrhenius preexponential factor in s^-1
		real(8)			::	n			! Arrhenius temperature exponent
		real(8)			::	Ea			! Arrhenius activation energy in kJ/mol
	end type

	! Struct: Reaction
	! 
	! Contains the information for a single reaction.
	type Reaction
		integer								::	isomerList(2)	! Isomer wells connected to transition state
		real(8)								:: 	E0				! Ground-state electronic + zero-point energy of transition state in kJ/mol
		type(ArrheniusKinetics)				::	kinetics		! Arrhenius kinetics for high pressure limit
		real(8), dimension(:), allocatable 	::	kf				! Microcanonical forward reaction rate in s^-1
		real(8), dimension(:), allocatable 	::	kb				! Microcanonical backward reaction rate in s^-1
	end type

contains
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!
	! Subroutine: microRates()
	!
	! Calculates the microcanonical forward and backward rates for a given 
	! reaction at the given conditions.
	subroutine microRates(T, P, rxn, E, isomerList)
	
		real(8), intent(in) :: T
		real(8), intent(in) :: P
		real(8), dimension(:), intent(in) :: E
		type(Reaction), intent(inout) :: rxn
		type(Isomer), dimension(:), intent(in) :: isomerList
		
		type(Isomer) :: reac, prod
		integer nGrains, r
		real(8) dE, Keq
		
		nGrains = size(E)
		dE = E(2) - E(1)

		! Allocate and zero rate coefficients
		if (.not. allocated(rxn%kf)) then
			allocate( rxn%kf(1:nGrains), rxn%kb(1:nGrains) )
		end if
		rxn%kf = 0 * rxn%kf
		rxn%kb = 0 * rxn%kb
		
		reac = isomerList(rxn%isomerList(1))
		prod = isomerList(rxn%isomerList(2))
		
		if (reac%numSpecies == 1 .and. prod%numSpecies == 1) then
    
			! Calculate forward rate coefficient via inverse Laplace transform
			call rateILT(rxn%E0 - reac%E0, reac%densStates, &
				rxn%kinetics, T, E, rxn%kf)
				
			! Calculate backward rate coefficient via detailed balance
			do r = 1, nGrains
				if (prod%densStates(r) /= 0) then
					rxn%kb(r) = rxn%kf(r) * reac%densStates(r) / prod%densStates(r)
				end if
			end do
			
		elseif (reac%numSpecies == 1 .and. prod%numSpecies > 1) then
			
			! Calculate forward rate coefficient via inverse Laplace transform
			call rateILT(rxn%E0 - reac%E0, reac%densStates, &
				rxn%kinetics, T, E, rxn%kf)
				
			! Calculate equilibrium constant
			Keq = equilCoeff(prod, reac, T, P)
			
			! Calculate backward rate coefficient via detailed balance
			! assuming multimolecular isomer is fully thermalized
			do r = 1, nGrains
				rxn%kb(r) = Keq * rxn%kf(r) * &
					reac%densStates(r) * exp(-E(r) * 1000 / 8.314472 / T) / reac%Q * dE
			end do

		elseif (reac%numSpecies > 1 .and. prod%numSpecies == 1) then
			
			! Convert to dissocation so that a similar algorithm to above can be used
			reac = isomerList(rxn%isomerList(2))
			prod = isomerList(rxn%isomerList(1))
		
			! Calculate forward rate coefficient via inverse Laplace transform
			call rateILT(rxn%E0 - reac%E0, reac%densStates, &
				rxn%kinetics, T, E, rxn%kb)
				
			! Calculate equilibrium constant
			Keq = equilCoeff(prod, reac, T, P)
			
			! Calculate backward rate coefficient via detailed balance
			! assuming multimolecular isomer is fully thermalized
			do r = 1, nGrains
				rxn%kf(r) = Keq * rxn%kb(r) * &
					reac%densStates(r) * exp(-E(r) * 1000 / 8.314472 / T) / reac%Q * dE
			end do

		end if
		
		!write (*,*) rxn%isomerList(1), rxn%isomerList(2), rxn%kf(nGrains), rxn%kb(nGrains)

	end subroutine	
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!
	! Function: rateILT()
	!
	!   Calculates the microcanonical rate coefficient as a function of energy
	!	using inverse Laplace transform of modified Arrhenius parameters for
	!	k(T) at high pressure. The modified Arrhenius equation is 
	!
	!		k(T) = A * T^n * exp(-Ea / R / T)
	!
	!	and the result from ILT is
	!
	!		k(E) = A * T^n * rho(E - Ea) / rho(E)
	!
	!	where we have neglected the effect of T^n in the inverse transform
	!	because it's too hard.
	!
	! Parameters:
	!   rho = density of states for the reactant in (kJ/mol)^-1
	!   E0 = electronic + zero-point energies of ground state of transition 
	!       state in kJ/mol
	!	kinetics = Arrhenius kinetics for reaction in high P limit
	!	T = absolute temperature in K
	!   E = vector of energies at which k(E) will be evaluated in kJ/mol
	!   k = RRKM microcanonical rate coefficents evaluated at each E in s^-1
	subroutine rateILT(E0, rho, kinetics, T, E, k)

		! Provide parameter type-checking
		real(8), intent(in)					::	E0
		real(8), dimension(:), intent(in)	:: 	rho
		type(ArrheniusKinetics), intent(in)	::	kinetics
		real(8), intent(in)					::	T
		real(8), dimension(:), intent(in)	:: 	E
		real(8), dimension(:), intent(out)	:: 	k
		
		real(8) dE
		integer		:: r, s
		
		dE = E(2) - E(1)
		s = floor(kinetics%Ea / dE)
		if (s < 0) s = 0
		
		! Determine rate coefficients using inverse Laplace transform
		do r = 1, size(E)
			if (s >= r .or. rho(r) == 0) then
				k(r) = 0
			else
				k(r) = kinetics%A * (T ** kinetics%n) * rho(r - s) / rho(r)
			end if
		end do
		
	end subroutine

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!
	! Function: equilCoeff()
	!
	!   Calculates the macroscopic equilibrium coefficient for a reaction.
	!
	! Parameters:
	!   reac = reactant isomer
	!   prod = product isomer
	!	T = absolute temperature in K
	!   P = absolute pressure in bar
	function equilCoeff(reac, prod, T, P)

		type(Isomer), intent(in) 	::	reac
		type(Isomer), intent(in) 	::	prod
		real(8), intent(in)			::	T
		real(8), intent(in)			::	P
		real(8) 					:: 	equilCoeff
		
		equilCoeff = eqRatio(prod%dGf, prod%dHf, T) / eqRatio(reac%dGf, reac%dHf, T)
		!equilCoeff = equilCoeff * (P / 1.0)

	end function
	
end module
