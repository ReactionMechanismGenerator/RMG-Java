! ==============================================================================
!
! 	StrongCollision.f90
!
! 	Written by Josh Allen (jwallen@mit.edu)
!
! ==============================================================================

module StrongCollisionModule

	use SimulationModule
	use SpeciesModule
	use IsomerModule
	use ReactionModule

	implicit none
	
contains

	! --------------------------------------------------------------------------
	
	! Subroutine: ssmscRates()
	! 
	! Uses the steady state/modified strong collision approach of Chang,
	! Bozzelli, and Dean to estimate the phenomenological rate coefficients 
	! from a full ME matrix for the case of a chemically activated system.
	!
	! Parameters:
	!
	subroutine ssmscRates(T, P, simData, speciesList, isomerList, rxnList, K)
	
		! Provide parameter type checking of inputs and outputs
		real(8), intent(in)							:: 	T
		real(8), intent(in)							:: 	P
		type(Simulation), intent(in)				:: 	simData
		type(Species), dimension(:), intent(in)		:: 	speciesList
		type(Isomer), dimension(:), intent(in)		:: 	isomerList
		type(Reaction), dimension(:), intent(in)	:: 	rxnList
		real(8), dimension(:,:), intent(out) 		:: 	K
		
		! Steady-state populations
		real(8), dimension(:,:), allocatable		:: 	pa
		
		! Steady-state matrix and vector
		real(8), dimension(:,:), allocatable		:: 	A
		real(8), dimension(:), allocatable			:: 	b
		! Collision efficiencies (the "modified" in modified strong collision)
		real(8), dimension(:), allocatable			:: 	beta
		! Number of active-state energy grains for each unimolecular isomer
		integer, dimension(:), allocatable				:: 	nAct
		! Indices i and j represent sums over unimolecular wells
		integer											::	i, j
		! Indices m and n represent sums over bimolecular sources/sinks
		integer											::	m, n
		! Indices r and s represent sums over energy grains
		integer											::	r, s
		integer											::	u
		! Variables for BLAS and LAPACK
		integer, dimension(:), allocatable				::	iPiv
		integer											::	info
		integer	src, start
		real(8)	val
		
		integer nUni, reac, prod
		
		! Zero rate coefficient matrix
		do i = 1, simData%nIsom
			do j = 1, simData%nIsom
				K(i,j) = 0.0
			end do
		end do
		
		! Determine number of unimolecular wells
		nUni = 0
		do i = 1, simData%nIsom
			if (isomerList(i)%numSpecies == 1) nUni = nUni + 1
		end do
		
		! Find steady-state populations at each grain
		allocate( A(1:nUni, 1:nUni) )
		allocate( b(1:nUni) )
		allocate( iPiv(1:nUni) )
		allocate( pa(1:simData%nGrains, 1:nUni) )
			
		allocate( beta(1:simData%nIsom) )
			
		! Determine collision efficiency
		call collisionEfficiencies(T, simData, isomerList, rxnList, beta)

		do src = 1, simData%nIsom

			! Determine starting grain
			start = activeSpaceStart(simData, isomerList, rxnList, src)
			
			do r = 1, start-1
				do i = 1, nUni
					pa(r,i) = 0.0
				end do
			end do
			
			do r = start, simData%nGrains

				! Zero A matrix and b vector
				do i = 1, nUni
					do j = 1, nUni
						A(i,j) = 0.0
					end do
					b(i) = 0.0
				end do
				
				! Collisional deactivation
				do i = 1, nUni
					A(i,i) = - beta(i) * isomerList(i)%omega
				end do
				
				! Reactions
				do u = 1, simData%nRxn
					reac = rxnList(u)%isomerList(1)
					prod = rxnList(u)%isomerList(2)
					if (isomerList(reac)%numSpecies == 1 .and. isomerList(prod)%numSpecies == 1) then
						A(prod,reac) = rxnList(u)%kf(r)
						A(reac,reac) = A(reac,reac) - rxnList(u)%kf(r)
						A(reac,prod) = rxnList(u)%kb(r)
						A(prod,prod) = A(prod,prod) - rxnList(u)%kb(r)
					elseif (isomerList(reac)%numSpecies == 1 .and. isomerList(prod)%numSpecies > 1) then
						A(reac,reac) = A(reac,reac) - rxnList(u)%kf(r)
					elseif (isomerList(reac)%numSpecies > 1 .and. isomerList(prod)%numSpecies == 1) then
						A(prod,prod) = A(prod,prod) - rxnList(u)%kb(r)
					end if
				end do
					
				! Activation
				if (isomerList(src)%numSpecies == 1) then
					b(src) = beta(src) * isomerList(src)%omega * isomerList(src)%eqDist(r)
				else
					do u = 1, simData%nRxn
						reac = rxnList(u)%isomerList(1)
						prod = rxnList(u)%isomerList(2)
						if (reac == src) then
							b(prod) = rxnList(u)%kf(r)
						elseif (prod == src) then
							b(reac) = rxnList(u)%kb(r)
						end if
					end do
				end if
				
				! Solve for steady-state population
				call DGESV( nUni, 1, A, nUni, iPiv, b, nUni, info )
				if (info > 0) then
					write (*,*), "A singular matrix was encountered! Aborting."
					stop
				end if
				pa(r,:) = -b
				
			end do
			
			! Calculate rates
			
			! Stabilization rates (i.e.) R + R' --> Ai or M --> Ai
			do i = 1, nUni
				if (i /= src) then
					val = beta(i) * isomerList(i)%omega * sum(pa(:,i))
					K(i,src) = K(i,src) + val
					K(src,src) = K(src,src) - val
				end if
			end do
			
			! Dissociation rates (i.e.) R + R' --> Bn + Cn or M --> Bn + Cn
			do u = 1, simData%nRxn
				reac = rxnList(u)%isomerList(1)
				prod = rxnList(u)%isomerList(2)
				if (reac <= nUni .and. prod > nUni .and. prod /= src) then
					val = sum(rxnList(u)%kf * pa(:,reac))
					K(prod,src) = K(prod,src) + val
					K(src,src) = K(src,src) - val
				elseif (reac > nUni .and. reac /= src .and. prod <= nUni) then
					val = sum(rxnList(u)%kb * pa(:,prod))
					K(reac,src) = K(reac,src) + val
					K(src,src) = K(src,src) - val
				end if
			end do
			
		end do

		! Clean up
		deallocate( A, b, iPiv, pa, beta )

	end subroutine

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!
	! collisionEfficiencies() 
	!
	!   Computes the fraction of collisions that result in deactivation.
	!
	subroutine collisionEfficiencies(T, simData, isomerList, rxnList, beta)
	
		! Provide parameter type checking of inputs and outputs
		real(8), intent(in)							::	T
		type(Simulation), intent(in)				:: 	simData
		type(Isomer), dimension(:), intent(in)		:: 	isomerList
		type(Reaction), dimension(:), intent(in)	:: 	rxnList
		real(8), dimension(:), intent(out)			::	beta
		
		real(8) E0
		integer i, r
		
		do i = 1, simData%nIsom
			if (isomerList(i)%numSpecies == 1) then
			
				E0 = simData%Emax
	 			do r = 1, size(rxnList)
	 				if (rxnList(r)%isomerList(1) == i .or. rxnList(r)%isomerList(2) == i) then
	 					if (rxnList(r)%E0 < E0) E0 = rxnList(r)%E0
	 				end if
	 			end do
	 
	 			beta(i) = efficiency(T, &
	 				simData%alpha * (1000 / 6.022e23), &
	 				isomerList(i)%densStates / (1000 / 6.022e23), &
	 				simData%E * (1000 / 6.022e23), &
	 				E0 * (1000 / 6.022e23))
			else
				beta(i) = 0.
			end if
		end do
		
		!write (*,*) beta
		
	
	end subroutine
	
	! --------------------------------------------------------------------------
	! 
	! efficiency() 
	!
	!   Computes the fraction of collisions that result in deactivation. All
	!	parameters are assumed to be in SI units.
	!
	function efficiency(T, alpha, rho, E, E0)

		! Provide parameter type checking of inputs and outputs
		real(8), intent(in)					:: 	T
		real(8), intent(in)					:: 	alpha
		real(8), dimension(:), intent(in)	:: 	rho
		real(8), dimension(:), intent(in)	:: 	E
		real(8), intent(in)					:: 	E0
		real(8)								::	efficiency
		
		real(8) kB, Delta, dE
		real(8) Fe, FeNum, FeDen
		real(8) Delta1, Delta2, DeltaN
		real(8) temp
		integer r
		
		kB = 1.381e-23						! [=] J/K
		Delta = 1
		dE = E(2) - E(1)					! [=] J

		FeNum = 0
		FeDen = 0
		do r = 1, size(E)
			temp = rho(r) * exp(-E(r) / kB / T) * dE
			if (E(r) > E0) then
				FeNum = FeNum + temp * dE
				if (FeDen == 0) FeDen = temp * kB * T
			end if
		end do

		Fe = FeNum / FeDen

		Delta1 = 0
		Delta2 = 0
		DeltaN = 0
		do r = 1, size(E)
			temp = rho(r) * exp(-E(r) / kB / T) * dE
			if (E(r) < E0) then
				Delta1 = Delta1 + temp * dE
				Delta2 = Delta2 + temp * exp(-(E0 - E(r)) / (Fe * kB * T)) * dE
			end if
			DeltaN = DeltaN + temp * dE
		end do
		Delta1 = Delta1 / DeltaN
		Delta2 = Delta2 / DeltaN

		Delta = Delta1 - (Fe * kB * T) / (alpha + Fe * kB * T) * Delta2
		
		efficiency = (alpha / (alpha + Fe * kB * T))**2 / Delta
		
		! Place bounds on efficiency range (should probably also warn user!)
		if (efficiency < 1e-8) then
			efficiency = 1e-8
		elseif (efficiency > 1.0) then
			efficiency = 1.0
		end if
		
		!write (*,*), FeNum, FeDen, Fe, Delta, efficiency
		
	end function
	
	! --------------------------------------------------------------------------

	! Subroutine: activeSpaceStart()
	! 
	! Determines the grain below which the reservoir approximation will be used
	! and above which the pseudo-steady state approximation will be used by
	! examining the energies of the transition states connected to each 
	! unimolecular isomer.
	!
	! Parameters:
	!   simData - The simulation parameters.
	!   uniData - The chemical data about each unimolecular isomer.
	!   rxnData - The chemical data about each transition state.
	!
	! Returns:
	!	nRes - The reservoir cutoff grains for each unimolecular isomer.
	function activeSpaceStart(simData, isomerList, rxnList, well)
	
		type(Simulation), intent(in)				:: 	simData
		type(Isomer), dimension(:), intent(in)		:: 	isomerList
		type(Reaction), dimension(:), intent(in)	:: 	rxnList
		integer, intent(in)							::	well
		integer										::	activeSpaceStart
		
		real(8) Eres
		
		integer t, start
		
		start = simData%nGrains
		
		Eres = simData%Emax
		do t = 1, simData%nRxn
			if (rxnList(t)%isomerList(1) == well .or. rxnList(t)%isomerList(2) == well) then
				if (rxnList(t)%E0 < Eres) then
					Eres = rxnList(t)%E0
				end if
			end if
		end do
		activeSpaceStart = ceiling((Eres - simData%Emin) / simData%dE) + 1

	end function
	
	! --------------------------------------------------------------------------

end module
