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
!	Version: 0.2.0							Date: 14 Jan 2009
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
	use ReactionModule
	use DensityOfStatesModule
	use FullEGMEModule
	use FastEGMEModule
	use StrongCollisionModule
	use RateModelModule
	
	! Verbosity of console output flag
	integer verbose

	! The simulation parameters
	type(Simulation) 								::	simData
	! The unimolecular isomer data
	type(UniWell), dimension(:), allocatable		:: 	uniData	
	! The bimolecular source/sink data
	type(MultiWell), dimension(:), allocatable		:: 	multiData
	! The reaction data
	type(Reaction), dimension(:), allocatable		:: 	rxnData
	! Mi(r,s,i) = Collisional transition from energy s to energy r for unimolecular well i
	real(8), dimension(:,:,:), allocatable			:: 	Mi
	! Hn(r,s,n) = Collisional transition from energy s to energy r for bimolecular well n
	real(8), dimension(:,:,:), allocatable			:: 	Hn
	! Kij(r,i,j) = Reactive transition from unimolecular well j to unimolecular well i at energy grain r
	real(8), dimension(:,:,:), allocatable			:: 	Kij
	! Fim(r,i,m) = Reactive transition from bimolecular well n to unimolecular well i at energy grain r
	real(8), dimension(:,:,:), allocatable			:: 	Fim
	! Gnj(r,n,j) = Reactive transition from unimolecular well j to bimolecular well n at energy grain r
	real(8), dimension(:,:,:), allocatable			:: 	Gnj
	! Jnm(r,n,m) = Reactive transition from bimolecular well m to bimolecular well n at energy grain r
	real(8), dimension(:,:,:), allocatable			:: 	Jnm
	! bi(r,i) = Equilibrium distribution at energy r for unimolecular well i
	real(8), dimension(:,:), allocatable			:: 	bi
	! bn(r,n) = Equilibrium distribution at energy r for bimolecular well n
	real(8), dimension(:,:), allocatable			:: 	bn
	! K(i,j) = Phenomenological rate coefficients for transitions from species i to species j
	real(8), dimension(:,:,:,:), allocatable 		:: 	K
	! The reservoir cutoff grains for each unimolecular isomer
	integer, dimension(:), allocatable 				:: 	nRes
	! chebCoeff(i,j,t,p) = Coefficient matrix for product of Chebyshev polynomials phi_t(Tred) * phi_p(Pred) for reaction j --> i 
	real(8), dimension(:,:,:,:), allocatable 		:: 	chebCoeff
	
	! Indices
	integer t, p, i, j
	integer found
	
	verbose = 1

	! Part I: Load data from files on disk
	! ------------------------------------

	! Open file for reading; fail if unable to open
	open(1, iostat=ios, file='fame_input.txt', action='read', status='old', access='sequential', recl=128)
	if (ios /= 0) then
		print *, 'Error: Unable to open file "', path, '".'
		return
	end if

	if (verbose >= 1) write (*,*), 'Reading from file "input.txt"...'

	! Load simulation data from input file
	if (verbose >= 2) write (*,*), '\tReading simulation parameters...'
	call loadSimulationData(simData)
			
	! Load unimolecular isomer well data from input file
	if (verbose >= 2) write (*,*), '\tReading unimolecular isomer data...'
	allocate (uniData(1:simData%nUni))
	do i = 1, simData%nUni
		call loadUniWellData(uniData(i), simData%nGrains)
	end do

	! Load multimolecular isomer well data from input file
	if (verbose >= 2) write (*,*), '\tReading multimolecular isomer data...'
	allocate (multiData(1:simData%nMulti))
	do i = 1, simData%nMulti
		call loadMultiWellData(multiData(i), simData%nGrains)
	end do

	! Load reaction data from input file
	if (verbose >= 2) write (*,*), '\tReading reaction data...'
	allocate (rxnData(1:simData%nRxn))
	do i = 1, simData%nRxn
		call loadReactionData(rxnData(i), simData%nGrains)
	end do
	
	! DEBUG: output parameters to console
	if (verbose >= 3) then
		call writeSimulationData(simData)
		call writeIsomerData(uniData, multiData)
		call writeReactionData(rxnData)
	end if

	! Close file
	close(1)
 
	! Calculate density of states for each well
	if (verbose >= 1) write (*,*), 'Calculating density of states...'
	call calcDensityOfStates(simData%E, uniData, multiData, verbose)

	! Allocate memory for master equation matrices
	allocate( 	Mi( 1:simData%nGrains, 1:simData%nGrains, 1:simData%nUni)	)
	allocate( 	Hn( 1:simData%nGrains, 1:simData%nGrains, 1:simData%nMulti)	)
	allocate( 	Kij(1:simData%nGrains, 1:simData%nUni,    1:simData%nUni) 	)
	allocate( 	Fim(1:simData%nGrains, 1:simData%nUni,    1:simData%nMulti) )
	allocate( 	Gnj(1:simData%nGrains, 1:simData%nMulti,  1:simData%nUni) 	)
	allocate( 	Jnm(1:simData%nGrains, 1:simData%nMulti,  1:simData%nMulti)	)
	! Allocate memory for Boltzmann vectors
	allocate( 	bi( 1:simData%nGrains, 1:simData%nUni)		)
	allocate( 	bn( 1:simData%nGrains, 1:simData%nMulti)	)
	! Allocate memory for reservoir cutoffs	
	allocate( 	nRes(1:simData%nUni) 	)
	! Allocate memory for phenomenological rate coefficient matrix	
	allocate( 	K( size(simData%Tlist), size(simData%Plist), &
		1:(simData%nUni+simData%nMulti), 1:(simData%nUni+simData%nMulti) )	)

	if (verbose >= 1) write (*,*), 'Calculating k(T, P)...'
	
	do t = 1, size(simData%Tlist)
		
		simData%T = simData%Tlist(t)
		
		! Calculate the equilibrium (Boltzmann) distributions
		if (verbose >= 2) write (*,*), '\tDetermining equilibrium distributions at T =', simData%T, 'K...'
		call calcEqDists(uniData, multiData, simData%E, simData%T, bi, bn)
				
		do p = 1, size(simData%Plist)

			! Part II: Construct full ME matrix
			! ---------------------------------
	
			simData%P = simData%Plist(p)

			if (verbose >= 2) write (*,*), '\tCalculating k(T, P) at T =', simData%T, 'K, P =', simData%P, 'bar...'

			! Reactive terms in master equation
			if (verbose >= 3) write (*,*), '\t\tDetermining reactive terms in full master equation...'
			call ME_reaction(simData, uniData, multiData, rxnData, Mi, Hn, Kij, Fim, Gnj, Jnm, bi, bn)
			
			! Collisional terms in master equation
			if (verbose >= 3) write (*,*), '\t\tDetermining collisional terms in full master equation...'
			call ME_collision(simData, uniData, Mi, bi)
			
			! Part III: Determine approximate phenomenological rate coefficients
			! ------------------------------------------------------------------
		
			! Determine reservoir cutoff grains for each unimolecular isomer
			if (verbose >= 3) write (*,*), '\t\tDetermining reservoir cutoff grains...'
			call getReservoirCutoffs(simData, uniData, rxnData, nRes)
			do i = 1, size(nRes)
				j = ceiling(uniData(i)%E / simData%dE) + 1
				if (nRes(i) < j .or. nRes(i) >= simData%nGrains) then
					write (*,*), 'ERROR: Invalid reservoir grain.'
					stop
				end if
			end do
			
			if (simData%mode == 1) then
				! Apply steady state/reservoir state approximations
				if (verbose >= 3) write (*,*), '\t\tApplying steady-state/modified strong collision approximation...'
				call ssmscRates(simData, uniData, multiData, rxnData, Kij, Fim, Gnj, Jnm, bi, bn, K(t,p,:,:))
			elseif (simData%mode == 2) then
				! Apply steady state/reservoir state approximations
				if (verbose >= 3) write (*,*), '\t\tApplying steady-state/reservoir-state approximation...'
				call ssrsRates(simData, uniData, rxnData, nRes, Mi, Hn, Kij, Fim, Gnj, Jnm, bi, bn, K(t,p,:,:))
			else
				write (*,*), 'ERROR: An invalid solution mode was provided!'
				stop
			end if

		end do
	end do

 	!do i = 1, simData%nUni + simData%nMulti
	!	write (*,*), K(4,3,i,:)
	!end do
					
	! Fit k(T, P) to approximate formula
	! Also test for validity of fitted rate coefficients
	if (verbose >= 1) write (*,*), 'Fitting k(T,P) to model...'
	allocate( chebCoeff(simData%nChebT, simData%nChebP, &
		simData%nUni+simData%nMulti, simData%nUni+simData%nMulti) )
	do i = 1, simData%nUni + simData%nMulti
		do j = 1, simData%nUni + simData%nMulti
			if (i /= j) then
				found = 0
				do t = 1, size(simData%Tlist)
					do p = 1, size(simData%Tlist)
						if (K(t,p,i,j) == 0.) then
							found = 1
						end if
					end do
				end do
				if (found == 1) then
					!write(*,*), 'ERROR in FAME execution: Rate coefficient(s) not properly estimated!'
					!stop
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

	! Write output file
	if (verbose >= 1) write (*,*), 'Saving results...'
	call saveResults('fame_output.txt', simData, chebCoeff)
	
	! Free memory
	deallocate( Fim, Gnj, Kij, Jnm, Mi, Hn, bi, bn, Nres, K, chebCoeff )

	if (verbose >= 1) write (*,*), 'DONE!'
	

end program

! ==============================================================================
