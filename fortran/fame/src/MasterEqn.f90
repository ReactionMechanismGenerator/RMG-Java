! ==============================================================================
!
! 	FullEGME.f90
!
! 	Written by Josh Allen (jwallen@mit.edu)
!
! ==============================================================================

module MasterEqnModule

	use SimulationModule
	use IsomerModule
	use ReactionModule

	implicit none
	
contains

	! ==========================================================================
	!
	! Function: masterEqn()
	!
	!   Calculates the full master equation matrix.
	!
	! Parameters:
	!   simData - The simulation parameters.
	!   uniData - The chemical data about each unimolecular isomer.
	!   multiData - The chemical data about each multimolecular source/sink.
	!   rxnData - The chemical data about each reaction.
	! 	Mi - Collisional transition matrix for each unimolecular well
	! 	Hn - Collisional transition matrix for each multimolecular well
	! 	Kij - Reactive transition from unimolecular well to unimolecular well
	! 	Fim - Reactive transition from multimolecular well to unimolecular well
	! 	Gnj - Reactive transition from unimolecular well to multimolecular well
	! 	Jnm - Reactive transition from multimolecular well to multimolecular well
	!	bi - Equilibrium distributions for unimolecular wells
	subroutine masterEqn(simData, uniData, multiData, rxnData, Mi, Hn, Kij, Fim, Gnj, Jnm, bi)

		type(Simulation), intent(in)				:: 	simData
		type(Isomer), dimension(:), intent(in)		:: 	uniData
		type(Isomer), dimension(:), intent(in)		:: 	multiData
		type(Reaction), dimension(:), intent(in)	:: 	rxnData
		real(8), dimension(:,:,:), intent(inout)	:: 	Mi
		real(8), dimension(:,:,:), intent(inout)	:: 	Hn
		real(8), dimension(:,:,:), intent(inout)	:: 	Kij
		real(8), dimension(:,:,:), intent(inout)	:: 	Fim
		real(8), dimension(:,:,:), intent(inout)	:: 	Gnj
		real(8), dimension(:,:,:), intent(inout)	:: 	Jnm
		real(8), dimension(:,:), intent(in)			:: 	bi
		
		Kij = 0 * Kij
		Fim = 0 * Fim
		Gnj = 0 * Gnj
		Mi = 0 * Mi
		Jnm = 0 * Jnm
		Hn = 0 * Hn
		
		! Reactive terms in master equation
		call ME_reaction(simData, uniData, multiData, rxnData, Mi, Hn, Kij, Fim, Gnj, Jnm, bi)
			
		! Collisional terms in master equation
		call ME_collision(simData, uniData, Mi, bi)
        
	end subroutine

	! ==========================================================================
	!
	! Function: ME_reaction()
	!
	!   Calculates the reaction components of the master equation matrix.
	!
	! Parameters:
	!   simData - The simulation parameters.
	!   uniData - The chemical data about each unimolecular isomer.
	!   multiData - The chemical data about each multimolecular source/sink.
	!   rxnData - The chemical data about each reaction.
	! 	Mi - Collisional transition matrix for each unimolecular well
	! 	Hn - Collisional transition matrix for each multimolecular well
	! 	Kij - Reactive transition from unimolecular well to unimolecular well
	! 	Fim - Reactive transition from multimolecular well to unimolecular well
	! 	Gnj - Reactive transition from unimolecular well to multimolecular well
	! 	Jnm - Reactive transition from multimolecular well to multimolecular well
	!	bi - Equilibrium distributions for unimolecular wells
	! 	bn - Equilibrium distributions for multimolecular wells
	subroutine ME_reaction(simData, uniData, multiData, rxnData, Mi, Hn, Kij, Fim, Gnj, Jnm, bi)

		type(Simulation), intent(in)				:: 	simData
		type(Isomer), dimension(:), intent(in)		:: 	uniData
		type(Isomer), dimension(:), intent(in)		:: 	multiData
		type(Reaction), dimension(:), intent(in)	:: 	rxnData
		real(8), dimension(:,:,:), intent(inout)	:: 	Mi
		real(8), dimension(:,:,:), intent(inout)	:: 	Hn
		real(8), dimension(:,:,:), intent(inout)	:: 	Kij
		real(8), dimension(:,:,:), intent(inout)	:: 	Fim
		real(8), dimension(:,:,:), intent(inout)	:: 	Gnj
		real(8), dimension(:,:,:), intent(inout)	:: 	Jnm
		real(8), dimension(:,:), intent(in)			:: 	bi
		
		! Indices
		integer		:: i, j, m, n, r, s, t
		
		real(8) Keq
		
		! Determine rate coefficients for each reaction
		! This is done by iterating over transition states
		do t = 1, simData%nRxn
			if (rxnData(t)%isomer(1) > simData%nUni .and. rxnData(t)%isomer(2) > simData%nUni) then			! bi <---> bi
				
				! Determine the wells involved
				n = rxnData(t)%isomer(1) - simData%nUni
				m = rxnData(t)%isomer(2) - simData%nUni
				
				! Not yet implemented
				
			elseif (rxnData(t)%isomer(1) > simData%nUni .or. rxnData(t)%isomer(2) > simData%nUni) then		! uni <---> bi
				
				! Determine the wells involved
				if (rxnData(t)%isomer(1) > simData%nUni) then
					n = rxnData(t)%isomer(1) - simData%nUni
					i = rxnData(t)%isomer(2)
				else
					i = rxnData(t)%isomer(1)
					n = rxnData(t)%isomer(2) - simData%nUni
				end if
				
				! Calculate rate coefficient for uni --> multi using ILT method
				call rateILT(rxnData(t)%E - uniData(i)%E(1), uniData(i)%densStates, &
					rxnData(t)%arrh_A, rxnData(t)%arrh_n, rxnData(t)%arrh_Ea, &
					simData%T, simData%dE, simData%E, Gnj(:,n,i))
				
 				! Calculate rate coefficient for multi --> uni via detailed balance
				! assuming multimolecular isomer is fully thermalized
				Keq = eqRatio(sum(uniData(i)%G), sum(uniData(i)%H), simData%T) / &
					eqRatio(sum(multiData(n)%G), sum(multiData(n)%H), simData%T)    
				Keq = Keq * (simData%P / 1.0)
				do r = 1, simData%nGrains
					Fim(r,i,n) = Keq * Gnj(r,n,i) * bi(r,i)
				end do
				
			else															! uni <---> uni
				
				! Determine the wells involved
				i = rxnData(t)%isomer(1)
				j = rxnData(t)%isomer(2)
				
				! Calculate rate coefficient for i --> j using ILT method
				call rateILT(rxnData(t)%E - uniData(i)%E(1), uniData(i)%densStates, &
					rxnData(t)%arrh_A, rxnData(t)%arrh_n, rxnData(t)%arrh_Ea, &
					simData%T, simData%dE, simData%E, Kij(:,j,i))
				
 				! Calculate rate coefficient for j --> i via detailed balance
				do r = 1, simData%nGrains
					if (uniData(j)%densStates(r) /= 0) then
						Kij(:,i,j) = Kij(:,j,i) * &
							uniData(i)%densStates(r) / uniData(j)%densStates(r)
					end if
				end do
				
			end if
		end do
			
		! Diagonals of unimolecular wells
		do i = 1, simData%nUni
			do j = 1, simData%nUni
				if (j /= i) then
					do r = 1, simData%nGrains
						Mi(r,r,i) = Mi(r,r,i) - Kij(r,j,i)
					end do
				end if
			end do
			do m = 1, simData%nMulti
				do r = 1, simData%nGrains
					Mi(r,r,i) = Mi(r,r,i) - Gnj(r,m,i)
				end do
			end do
		end do
		
		! Diagonals of multimolecular wells
		do n = 1, simData%nMulti
			do j = 1, simData%nUni
				do r = 1, simData%nGrains
					Hn(r,r,n) = Hn(r,r,n) - Fim(r,j,n)
				end do
			end do
			do m = 1, simData%nMulti
				if (m /= n) then
					do r = 1, simData%nGrains
						Hn(r,r,n) = Hn(r,r,n) - Jnm(r,m,n)
					end do
				end if
			end do
		end do
        
	end subroutine

	! --------------------------------------------------------------------------

	! Function: ME_collision()
	!
	!   Calculates the collision components of the master equation matrix.
	!
	! Parameters:
	!   simData - The simulation parameters.
	!   uniData - The chemical data about each unimolecular isomer.
	! 	Mi - Collisional transition matrix for each unimolecular well
	!	bi - Equilibrium distributions for unimolecular wells
	subroutine ME_collision(simData, uniData, Mi, bi)

		type(Simulation), intent(in)				:: 	simData
		type(Isomer), dimension(:), intent(in)		:: 	uniData
		real(8), dimension(:,:,:), intent(inout)	:: 	Mi
		real(8), dimension(:,:), intent(in)			:: 	bi
			
		real(8)			:: gasConc			! Gas concentration in molecules/m^3
		real(8)			:: mu				! Reduced mass in g/mol
		real(8)			:: kB = 1.381e-23	! Boltzmann's constant in J/K
		real(8)			:: w				! Collision frequency in s^-1
		
		integer			:: probRange		! Range of significant transfer probability
		integer			:: lb				! Lower bound of significant transfer probability
		integer			:: ub				! Upper bound of significant transfer probability
		real(8)			:: colSum			! Used for column normalization
		
		real(8), dimension(:,:), allocatable	::	P	! Transfer probability matrix
		real(8), dimension(:), allocatable		::	C	! Vector of normalization coefficients
		
		integer			:: i, r, s
		integer 		:: start
		
		gasConc = simData%P * 1e5 / kB / simData%T
		
		allocate( P(1:simData%nGrains, 1:simData%nGrains) )
		allocate( C(1:simData%nGrains) )
		
		! Collisional energy transfer contributions
		do i = 1, simData%nUni
		
			start = ceiling((uniData(i)%E(1) - simData%Emin) / simData%dE)
			
			! Determine collision frequency for the current isomer
			mu = 1/(1/uniData(i)%MW(1) + 1/simData%bathGas%MW) / 6.022e26
			call collisionFrequency(simData%T, 0.5 * (uniData(i)%sigma(1) + simData%bathGas%sigma), &
				0.5 * (uniData(i)%eps(1) + simData%bathGas%eps), mu, gasConc, w)
			
			P = 0 * P
			C = 0 * C
			
			! Determine unnormalized entries in collisional tranfer probability matrix for the current isomer
			probRange = ceiling(16 * simData%alpha / simData%dE)
			do r = start, simData%nGrains
				lb = max(r - probRange, start)
				ub = min(r + probRange, simData%nGrains)
				do s = lb, ub
					P(s,r) = transferRate(s, r, simData%E(s), simData%E(r), simData%alpha, &
						uniData(i)%E(1), uniData(i)%densStates, simData%T)
				end do
			end do
			
			! Normalize using detailed balance
			do r = simData%nGrains, start, -1
				C(r) = (1 - sum(C(r+1:simData%nGrains) * P(r+1:simData%nGrains,r))) / sum(P(1:r,r))
			end do
			do r = start, simData%nGrains
				P(r,1:r-1) = P(r,1:r-1) * C(r)
				P(1:r-1,r) = P(1:r-1,r) * C(r)
				P(r,r) = P(r,r) * C(r) - 1
			end do
			
			! Add to ME matrix
			Mi(:,:,i) = Mi(:,:,i) + w * P
			
		end do

		deallocate( P, C )
		
	end subroutine
	
	! --------------------------------------------------------------------------

	! Subroutine: transferRate() 
	!
	!   Computes the unnormalized probability of a collisional energy transfer
	!   from state j (with energy Ej) to state i (with energy Ei) using a
	!   single-exponental-down model.
	!
	! Parameters:
	!	i = index of destination state
	!	j = index of source state
	!   Ei = energy of destination state i in cm^-1
	!   Ej = energy of source state j in cm^-1
	!   alpha = parameter in exponential-down model in (cm^-1)^-1
	!   E0 = electronic + zero-point energy of ground state in cm^-1
	!   f = equilibrium distribution
	!   rate = collisional transfer probability in s^-1
	function transferRate(i, j, Ei, Ej, alpha, E0, rho, T)
	
		! Provide parameter type-checking
		integer, intent(in)					:: 	i
		integer, intent(in)					:: 	j
		real(8), intent(in)					:: 	Ei
		real(8), intent(in)					:: 	Ej
		real(8), intent(in)					:: 	alpha
		real(8), intent(in)					:: 	E0
		real(8), dimension(:), intent(in)	:: 	rho
		real(8), intent(in)					:: 	T
		real(8)								:: 	transferRate
		
		real(8) R

		R = 8.314472 / 1000
		
		! Evaluate collisional transfer probability
		if (Ej < E0 .or. Ei < E0 .or. rho(j) .eq. 0) then
			transferRate = 0.
		elseif (Ej >= Ei) then
			transferRate = exp(-(Ej - Ei) / alpha)
		else
			transferRate = exp(-(Ei - Ej) / alpha) * rho(i) / rho(j) * exp( -(Ei - Ej) / (R * T))
		end if

	end function

	! --------------------------------------------------------------------------
		
	! collisionFrequency() 
	!
	!   Computes the Lennard-Jones (12-6) collision frequency.
	!
	! Input:
	!   T = absolute temperature in K
	!   sigma = effective Lennard-Jones well-minimum parameter for collision in m
	!   eps = effective Lennard-Jones well-depth parameter for collision in J
	!   mu = reduced mass of molecules involved in collision in kg
	!   Na = molecules of A per m^3
	!
	! Output:
	!   omega = collision frequency in molecules m^-3 s^-1
	!
	subroutine collisionFrequency(T, sigma, eps, mu, Na, omega)
	
		! Provide parameter type-checking
		real(8), intent(in)					:: 	T
		real(8), intent(in)					:: 	sigma
		real(8), intent(in)					:: 	eps
		real(8), intent(in)					:: 	mu
		real(8), intent(in)					:: 	Na
		real(8), intent(out)				:: 	omega
		
		real(8)		:: kB = 1.3806504e-23
		real(8)		:: collisionIntegral
		
		collisionIntegral = 1.16145 / T**0.14874 + 0.52487 / exp(0.77320 * T) + 2.16178 / exp(2.43787 * T) &
			-6.435/10000 * T**0.14874 * sin(18.0323 * T**(-0.76830) - 7.27371)
    
		! Evaluate collision frequency
		omega = collisionIntegral *	sqrt(8 * kB * T / 3.141592654 / mu) * 3.141592654 * sigma**2 * Na
	
	end subroutine
	
	! --------------------------------------------------------------------------
		
	! Function: rate_rrkm()
	!
	!   Calculates the microcanonical rate coefficient as a function of energy
	!	using RRKM theory.
	!
	! Parameters:
	!   N = density of states for the reactant in (kJ/mol)^-1
	!   G = sum of states for the transition state
	!   E0 = electronic + zero-point energies of ground state of transition 
	!       state in kJ/mol
	!   E = vector of energies at which k(E) will be evaluated in kJ/mol
	!   k = RRKM microcanonical rate coefficents evaluated at each E in s^-1
	subroutine rate_rrkm(N, G, E0, E, k)

		! Provide parameter type-checking
		real(8), dimension(:), intent(in)	:: 	N
		real(8), dimension(:), intent(in)	:: 	G
		real(8), intent(in)					::	E0
		real(8), dimension(:), intent(in)	:: 	E
		real(8), dimension(:), intent(out)	:: 	k
		
		integer		:: r
		
		! Convert G(E) and N(E) into RRKM rate coefficients
		do r = 1, size(E)
			if (E(r) < E0 .or. N(r) == 0) then
				k(r) = 0
			else
				k(r) = G(r) / N(r) * 1000 / 6.626e-34 / 6.022e23
			end if
		end do
		
	end subroutine
	
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
	!		k(E) = A * T^n * N(E - Ea) / N(E)
	!
	!	where we have neglected the effect of T^n in the inverse transform
	!	because it's too hard.
	!
	! Parameters:
	!   N = density of states for the reactant in (kJ/mol)^-1
	!   E0 = electronic + zero-point energies of ground state of transition 
	!       state in kJ/mol
	!	arrh_A = Arrhenius preexponential factor in s^-1
	!	arrh_n = Arrhenius temperature exponent
	!	arrh_Ea = Arrhenius activation energy in kJ/mol
	!	T = absolute temperature in K
	!	dE = Grain spacing in kJ/mol
	!   E = vector of energies at which k(E) will be evaluated in kJ/mol
	!   k = RRKM microcanonical rate coefficents evaluated at each E in s^-1
	subroutine rateILT(E0, N, arrh_A, arrh_n, arrh_Ea, T, dE, E, k)

		! Provide parameter type-checking
		real(8), intent(in)					::	E0
		real(8), dimension(:), intent(in)	:: 	N
		real(8), intent(in)					::	arrh_A
		real(8), intent(in)					::	arrh_n
		real(8), intent(in)					::	arrh_Ea
		real(8), intent(in)					::	T
		real(8), intent(in)					::	dE
		real(8), dimension(:), intent(in)	:: 	E
		real(8), dimension(:), intent(out)	:: 	k
		
		integer		:: r, s
		 
		s = floor(arrh_Ea / dE)
		if (s < 0) s = 0
		
		! Determine rate coefficients using inverse Laplace transform
		do r = 1, size(E)
			if (s >= r .or. N(r) == 0) then
				k(r) = 0
			else
				k(r) = arrh_A * (T ** arrh_n) * N(r - s) / N(r)
			end if
		end do
		
	end subroutine

	! ==========================================================================
	!
	! Function: fitReverseKinetics()
	!
	!   Estimates the Arrhenius expression for the reverse reaction.
	!
	!
	! Parameters:
	!   reac = The reactant isomer for the forward reaction
	!   prod = The product isomer for the forward reaction
	!   rxn = The forward reaction
	!
	! Returns:
	!	arrh2_A = Arrhenius preexponential factor in s^-1 for reverse reaction
	!	arrh2_n = Arrhenius temperature exponent for reverse reaction
	!	arrh2_Ea = Arrhenius activation energy in kJ/mol for reverse reaction
	subroutine fitReverseKinetics(reac, prod, rxn, arrh2_A, arrh2_Ea, arrh2_n)
	
		type(Isomer), intent(in)		:: reac
		type(Isomer), intent(in)		:: prod
		type(Reaction), intent(in)		:: rxn
		real(8), intent(out)			:: arrh2_A
		real(8), intent(out)			:: arrh2_Ea
		real(8), intent(out)			:: arrh2_n
				
		real(8) Hrxn, Grxn, Srxn, R, T0
		
		T0 = 298.15
		Hrxn = (sum(prod%H) - sum(reac%H)) * 1000.
		Grxn = (sum(prod%G) - sum(reac%G)) * 1000.
		Srxn = (Hrxn - Grxn) / T0
		R = 8.314472
		
		arrh2_Ea = rxn%arrh_Ea - Hrxn / 1000.
		arrh2_n = 0
		arrh2_A = (rxn%arrh_A * T0 ** rxn%arrh_n) * exp(-Srxn / R)   
	
	end subroutine

end module
