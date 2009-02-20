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
	use IsomerModule
	use ReactionModule
	use InputModule
	use DensityOfStatesModule
	use MasterEqnModule
	use StrongCollisionModule
	use ReservoirStateModule
	use RateModelModule
	
	implicit none

	! Verbosity of console output flag
	integer verbose

	! The simulation parameters
	type(Simulation) 								::	simData
	! The unimolecular isomer data
	type(Isomer), dimension(:), allocatable			:: 	uniData	
	! The bimolecular source/sink data
	type(Isomer), dimension(:), allocatable			:: 	multiData
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
	! K(i,j) = Phenomenological rate coefficients for transitions from species i to species j
	real(8), dimension(:,:,:,:), allocatable 		:: 	K
	! The reservoir cutoff grains for each unimolecular isomer
	integer, dimension(:), allocatable 				:: 	nRes
	! chebCoeff(i,j,t,p) = Coefficient matrix for product of Chebyshev polynomials phi_t(Tred) * phi_p(Pred) for reaction j --> i 
	real(8), dimension(:,:,:,:), allocatable 		:: 	chebCoeff
	
	real(8) Eref
    real(8) conc
	
    ! Indices
	integer t, p, i, j, n, m
	integer found
	integer badRate
	
	verbose = 0

	! Load data from files on disk
	if (verbose >= 1) write (*,*), 'Reading input...'
	call loadNetwork('fame_input.txt', simData, uniData, multiData, rxnData, verbose)
    
	! Calculate density of states for each well
	if (verbose >= 1) write (*,*), 'Calculating density of states...'
	do i = 1, simData%nUni
		if (verbose >= 2) write (*,*), '\tUnimolecular isomer', i
        call densityOfStates(simData, uniData(i))
	end do
	
	! Allocate memory for master equation matrices
	allocate( 	Mi( 1:simData%nGrains, 1:simData%nGrains, 1:simData%nUni)	)
	allocate( 	Hn( 1:simData%nGrains, 1:simData%nGrains, 1:simData%nMulti)	)
	allocate( 	Kij(1:simData%nGrains, 1:simData%nUni,    1:simData%nUni) 	)
	allocate( 	Fim(1:simData%nGrains, 1:simData%nUni,    1:simData%nMulti) )
	allocate( 	Gnj(1:simData%nGrains, 1:simData%nMulti,  1:simData%nUni) 	)
	allocate( 	Jnm(1:simData%nGrains, 1:simData%nMulti,  1:simData%nMulti)	)
	! Allocate memory for Boltzmann vectors
	allocate( 	bi( 1:simData%nGrains, 1:simData%nUni)		)
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
		do i = 1, simData%nUni
			call eqDist(uniData(i), simData%E, simData%T, bi(:,i), uniData(i)%Q)
		end do
			
		do p = 1, size(simData%Plist)

			! Part II: Construct full ME matrix
			! ---------------------------------
	
			simData%P = simData%Plist(p)

			if (verbose >= 2) write (*,*), '\tCalculating k(T, P) at T =', simData%T, 'K, P =', simData%P, 'bar...'

			Kij = 0 * Kij
			Fim = 0 * Fim
			Gnj = 0 * Gnj
			Mi = 0 * Mi
			Jnm = 0 * Jnm
			Hn = 0 * Hn
			
			! Reactive terms in master equation
			call ME_reaction(simData, uniData, multiData, rxnData, Mi, Hn, Kij, Fim, Gnj, Jnm, bi)
			
			if (simData%mode == 1) then
				! Apply steady state/reservoir state approximations
				if (verbose >= 3) write (*,*), '\t\tApplying steady-state/modified strong collision approximation...'
				call ssmscRates(simData, uniData, multiData, rxnData, Kij, Fim, Gnj, Jnm, bi, K(t,p,:,:))
			elseif (simData%mode == 2) then
				! Collisional terms in master equation
				call ME_collision(simData, uniData, Mi, bi)
				! Determine reservoir cutoff grains for each unimolecular isomer
				if (verbose >= 3) write (*,*), '\t\tDetermining reservoir cutoff grains...'
				nRes = reservoirCutoffs(simData, uniData, rxnData)
				! Apply steady state/reservoir state approximations
				if (verbose >= 3) write (*,*), '\t\tApplying steady-state/reservoir-state approximation...'
				call ssrsRates(simData, uniData, multiData, rxnData, nRes, Mi, Hn, Kij, Fim, Gnj, Jnm, bi, K(t,p,:,:))
			else
				write (*,*), 'ERROR: An invalid solution mode was provided!'
				stop
			end if
			
			! Convert multimolecular k(T, P) to appropriate units (i.e. cm^3 mol^-1 s^-1)
			conc = simData%P * 100000. / 8.314472 / simData%T / 1e6
			do n = simData%nUni+1, simData%nUni+simData%nMulti
				do i = 1, simData%nUni+simData%nMulti
					K(t,p,i,n) = K(t,p,i,n) / conc
				end do
			end do
			
		end do
	end do

 	!do i = 1, simData%nUni + simData%nMulti
 	!	write (*,*) K(4,3,i,:)
 	!end do

	! Fit k(T, P) to approximate formula
	! Also test for validity of fitted rate coefficients
	if (verbose >= 1) write (*,*), 'Fitting k(T,P) to model...'
	allocate( chebCoeff(simData%nChebT, simData%nChebP, &
		simData%nUni+simData%nMulti, simData%nUni+simData%nMulti) )
	badRate = 0
	do i = 1, simData%nUni + simData%nMulti
		do j = 1, simData%nUni + simData%nMulti
			if (i /= j) then
				found = 0
				do t = 1, size(simData%Tlist)
					do p = 1, size(simData%Plist)
						!if (K(t,p,i,j) <= 0. .or. isnan(K(t,p,i,j)) .or. &
						!	abs(K(t,p,i,j)) <= huge(K(t,p,i,j))) then
						if (K(t,p,i,j) <= 0) then
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
	call saveResults('fame_output.txt', simData, chebCoeff)
	
	! Free memory
	deallocate( Fim, Gnj, Kij, Jnm, Mi, Hn, bi, Nres, K, chebCoeff )
	!deallocate( bn )
	
	if (verbose >= 1) write (*,*), 'DONE!'
	

end program

! ==============================================================================
