! ==============================================================================
!
!	fame.f90
!
!	Written primarily by Josh Allen (jwallen@mit.edu) for use by RMG (Reaction 
!	Mechanism Generator).
!
!	Copyright (c) 2008 by the William H. Green research group.
!
! ==============================================================================

!	Program: fame
!
!	Version: 0.3.1							Date: 13 Mar 2009
!
! 	Calculates phenomenological rate coefficients k(T, P) for a given potential
! 	energy surface using provided thermochemical and density-of-state data for
! 	the unimolecular wells and bimolecular sources/sinks. The method used is
! 	the steady state/reservoir state approximation to the energy-grained master
! 	equation, a method introduced in 
!
! 	N. J. B. Green and Z. A. Bhatti. ``Steady-state Master Equation Methods.''
!		Phys. Chem. Chem. Phys. 9, p. 4275-4290 (2007).
!
!	Input to this program is provided by a set of text files that should be
!	placed in a folder named input/ in the directory containing the FAME
!	executable. Output from this program will be placed in a text file in a
!	folder named output/ in the directory containing the FAME executable.
!
program fame

	use SimulationModule
	use SpeciesModule
	use IsomerModule
	use ReactionModule
	use InputModule
	use DensityOfStatesModule
	use StrongCollisionModule
	use ReservoirStateModule
	use RateModelModule
	
	implicit none

	! Verbosity of console output flag
	integer verbose

	! The simulation parameters
	type(Simulation) 								::	simData
	! The species data
	type(Species), dimension(:), allocatable		:: 	speciesList	
	! The bimolecular source/sink data
	type(Isomer), dimension(:), allocatable			:: 	isomerList
	! The reaction data
	type(Reaction), dimension(:), allocatable		:: 	rxnList
	! K(i,j) = Phenomenological rate coefficients for transitions from species i to species j
	real(8), dimension(:,:,:,:), allocatable 		:: 	K
	! chebCoeff(i,j,t,p) = Coefficient matrix for product of Chebyshev polynomials phi_t(Tred) * phi_p(Pred) for reaction j --> i 
	real(8), dimension(:,:,:,:), allocatable 		:: 	chebCoeff
	! (T, P) = Current temperature and pressure
	real(8)											::	T
	real(8)											::	P
	! Indices
	integer u, v, w, i, j, n, m
	
	! Other variables
	real(8) conc
	integer found
	integer badRate
	
	verbose = 0

	! Load data from files on disk
	if (verbose >= 1) write (*,*), 'Reading input...'
	call loadNetwork('fame_input.txt', simData, speciesList, isomerList, rxnList, verbose)
    
	! Calculate density of states for each (unimolecular) well
	if (verbose >= 1) write (*,*), 'Calculating density of states...'
	do i = 1, simData%nIsom
		if (isomerList(i)%numSpecies == 1) then
			call densityOfStates(simData, isomerList(i), speciesList)
		end if
	end do
	
	allocate( 	K( 1:size(simData%Tlist), 1:size(simData%Plist), &
		1:(simData%nIsom), 1:(simData%nIsom) )	)

	if (verbose >= 1) write (*,*), 'Calculating k(T, P)...'
	
	do u = 1, size(simData%Tlist)
		
		T = simData%Tlist(u)
		
		! Calculate the equilibrium (Boltzmann) distributions
		if (verbose >= 2) write (*,*), '\tDetermining equilibrium distributions at T =', T, 'K...'
		do i = 1, simData%nIsom
			if (isomerList(i)%numSpecies == 1) then
				call eqDist(isomerList(i), simData%E, T)
			end if
		end do
		
		do v = 1, size(simData%Plist)

			P = simData%Plist(v)

			if (verbose >= 2) write (*,*), '\tCalculating k(T, P) at T =', T, 'K, P =', P, 'bar...'

			! Calculate collision frequencies
			do i = 1, simData%nIsom
				if (isomerList(i)%numSpecies == 1) then
					isomerList(i)%omega = collisionFrequency(T, P, simData, &
						speciesList(isomerList(i)%speciesList(1)))
				end if
			end do
			
			! Calculate microcanonical rate coefficients at the current conditions
			do w = 1, simData%nRxn
				call microRates(T, P, rxnList(w), simData%E, isomerList)
			end do
			
			if (simData%mode == 1) then
				! Apply steady state/reservoir state approximations
				if (verbose >= 3) write (*,*), '\t\tApplying steady-state/modified strong collision approximation...'
				call ssmscRates(T, P, simData, speciesList, isomerList, rxnList, K(u,v,:,:))
			elseif (simData%mode == 2) then
				if (verbose >= 3) write (*,*), '\t\tApplying steady-state/reservoir-state approximation...'
				call ssrsRates(T, P, simData, speciesList, isomerList, rxnList, K(u,v,:,:))
			else
				write (*,*), 'ERROR: An invalid solution mode was provided!'
				stop
			end if
			
			! Convert multimolecular k(T, P) to appropriate units (i.e. cm^3 mol^-1 s^-1)
			conc = P * 100000. / 8.314472 / T / 1e6
			do i = 1, simData%nIsom
				do j = 1, simData%nIsom
					if (isomerList(j)%numSpecies > 1) K(u,v,i,j) = K(u,v,i,j) / conc
				end do
			end do
			
		end do
	end do

   	!do i = 1, simData%nIsom
	!	write (*,*) K(4,3,i,:)
	!end do

	!do i = 1, size(simData%Tlist)
	!	write (*,*) K(i,:,2,1)
	!end do

	! Fit k(T, P) to approximate formula
	! Also test for validity of fitted rate coefficients
	if (verbose >= 1) write (*,*), 'Fitting k(T,P) to model...'
	allocate( chebCoeff(simData%nChebT, simData%nChebP, &
		simData%nIsom, simData%nIsom) )
	badRate = 0
	do i = 1, simData%nIsom
		do j = 1, simData%nIsom
			if (i /= j) then
				found = 0
				do u = 1, size(simData%Tlist)
					do v = 1, size(simData%Plist)
						if (K(u,v,i,j) <= 0) then
							found = 1
						end if
					end do
				end do
				if (found == 1) then
					badRate = 1
					chebCoeff(:,:,i,j) = 0 * chebCoeff(:,:,i,j)
				else
					if (verbose >= 2) then
						write (*,*), '\tFitting k(T,P) for isomers', i, 'and', j
					end if
					call fitRateModel(K(:,:,i,j), simData%Tlist, simData%Plist, simData%nChebT, simData%nChebP, chebCoeff(:,:,i,j))
				end if
			end if
		end do
	end do
	
	if (badRate == 1) then
		write(*,*), 'Warning: One ore more rate coefficients not properly estimated!'
	end if			
	
	! Write output file
	if (verbose >= 1) write (*,*), 'Saving results...'
	call saveResults('fame_output.txt', simData, K, chebCoeff)
	
	! Free memory
! 	do i = 1, simData%nSpecies
! 		deallocate( speciesList(i)%vibFreq, speciesList(i)%rotFreq, &
! 			speciesList(i)%hindFreq, speciesList(i)%hindBarrier )
! 	end do
	do i = 1, simData%nIsom
		deallocate( isomerList(i)%speciesList )
		if (isomerList(i)%numSpecies == 1) then
			deallocate( isomerList(i)%densStates, isomerList(i)%eqDist )
			if (allocated(isomerList(i)%Mcoll)) then
				deallocate(isomerList(i)%Mcoll)
			end if
		end if
	end do
	do i = 1, simData%nRxn
		deallocate( rxnList(i)%kf, rxnList(i)%kb )
	end do
		
	if (verbose >= 1) write (*,*), 'DONE!'
	

end program

! ==============================================================================
