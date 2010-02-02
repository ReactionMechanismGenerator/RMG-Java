!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	reaction.f90
!
!	Written by Josh Allen (jwallen@mit.edu)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module ReactionModule

	use SpeciesModule
	use IsomerModule
	
	implicit none
	
	type ArrheniusKinetics
		real(8)			::	A			! Arrhenius preexponential factor in s^-1
		real(8)			::	n			! Arrhenius temperature exponent
		real(8)			::	Ea			! Arrhenius activation energy in kJ/mol
	end type

	type Reaction
		character(len=256)					::	equation		! Chemical reaction equation
		integer								::	reac			! Isomer wells connected to transition state
		integer								::	prod			! Isomer wells connected to transition state
		real(8)								:: 	E0				! Ground-state electronic + zero-point energy of transition state in kJ/mol
		type(ArrheniusKinetics)				::	arrhenius		! Arrhenius kinetics for high pressure limit
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
	!
	subroutine microRates(rxn, T, E, isomerList, speciesList)
	
		type(Reaction), intent(inout) :: rxn
		real(8), intent(in) :: T
		real(8), dimension(:), intent(in) :: E
		type(Isomer), dimension(:), intent(in) :: isomerList
		type(Species), dimension(:), intent(in) :: speciesList
		
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
		
		reac = isomerList(rxn%reac)
		prod = isomerList(rxn%prod)
		
		if (size(reac%species) == 1 .and. size(prod%species) == 1) then
    
			! Calculate forward rate coefficient via inverse Laplace transform
			call rateILT(rxn%E0 - reac%E0, reac%densStates, &
				rxn%arrhenius, T, E, rxn%kf)
				
			! Calculate backward rate coefficient via detailed balance
			do r = 1, nGrains
				if (prod%densStates(r) /= 0) then
					rxn%kb(r) = rxn%kf(r) * reac%densStates(r) / prod%densStates(r)
				end if
			end do
			
		elseif (size(reac%species) == 1 .and. size(prod%species) > 1) then
			
			! Calculate forward rate coefficient via inverse Laplace transform
			call rateILT(rxn%E0 - reac%E0, reac%densStates, &
				rxn%arrhenius, T, E, rxn%kf)
				
			! Calculate equilibrium constant
			Keq = equilCoeff(reac, prod, T, speciesList)
			
			! Calculate backward rate coefficient via detailed balance
			! assuming multimolecular isomer is fully thermalized
			do r = 1, nGrains
				rxn%kb(r) = rxn%kf(r) / Keq * &
					reac%densStates(r) * exp(-E(r) / 8.314472 / T) / reac%Q * dE
			end do

		elseif (size(reac%species) > 1 .and. size(prod%species) == 1) then
			
			! Convert to dissocation so that a similar algorithm to above can be used
			reac = isomerList(rxn%prod)
			prod = isomerList(rxn%reac)
			
			! Calculate forward rate coefficient via inverse Laplace transform
			call rateILT(rxn%E0 - reac%E0, reac%densStates, &
				rxn%arrhenius, T, E, rxn%kb)
				
			! Calculate equilibrium constant
			Keq = equilCoeff(reac, prod, T, speciesList)
			
			! Calculate backward rate coefficient via detailed balance
			! assuming multimolecular isomer is fully thermalized
			do r = 1, nGrains
				rxn%kf(r) = rxn%kb(r) / Keq * &
					reac%densStates(r) * exp(-E(r) / 8.314472 / T) / reac%Q * dE
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
	function equilCoeff(reac, prod, T, speciesList)

		type(Isomer), intent(in) 	::	reac
		type(Isomer), intent(in) 	::	prod
		real(8), intent(in)			::	T
		type(Species), dimension(:), intent(in) :: speciesList
		real(8) 					:: 	equilCoeff
		
		real(8) dGrxn
		integer i
		
		! Determine free energy of reaction
		dGrxn = 0.0
		do i = 1, size(reac%species)
			dGrxn = dGrxn - freeEnergy(speciesList(reac%species(i))%thermo, T)
		end do
		do i = 1, size(prod%species)
			dGrxn = dGrxn + freeEnergy(speciesList(prod%species(i))%thermo, T)
		end do
		
		! Determine Ka
		equilCoeff = exp(-dGrxn / 8.314472 / T)
		! Determine Kc
		equilCoeff = equilCoeff * (100000.0 / 8.314472 / T)**(size(prod%species) - size(reac%species))
		
	end function

end module
