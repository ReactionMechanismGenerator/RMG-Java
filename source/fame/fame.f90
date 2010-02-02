!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	fame.f90
!
!	Copyright (c) 2008-2009 by Josh Allen.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program fame

	use SpeciesModule
	use IsomerModule
	use ReactionModule
	use NetworkModule
	use InputModule
	use StrongCollisionModule
	use ReservoirStateModule
	use ModelModule
	use OutputModule
	
	implicit none

	! The activated reaction network
	type(Network) net
	! K(t,p,i,j) = Phenomenological rate coefficients for transitions from species i to species j
	real(8), dimension(:,:,:,:), allocatable 		:: 	K
	! chebCoeff(t,p,i,j) = Coefficient matrix for product of Chebyshev polynomials phi_t(Tred) * phi_p(Pred) for reaction j --> i 
	real(8), dimension(:,:,:,:), allocatable 		:: 	chebyshevCoefficients
	! arrhenius(p,i,j) = Arrhenius kinetics matrix for logPInterpolate fitting model
	type(ArrheniusKinetics), dimension(:,:,:), allocatable 		:: 	pDepArrhenius
	! (T, P) = Current temperature and pressure
	real(8)											::	T, P
	! Indices
	integer u, v, w, i, j, n, m

	! Read network information from stdin
	call loadNetwork(net)
	
	! Determine energy range of calculation
	call selectEnergyGrains(net)
		
	! Calculate density of states for each (unimolecular) well
	do i = 1, size(net%isomers)
		if (size(net%isomers(i)%species) == 1) then
			call densityOfStates(net%isomers(i), net%species, net%E)
		end if
	end do
	! Allocate memory for k(T, P) data
	allocate( K( 1:size(net%Tlist), 1:size(net%Plist), 1:size(net%isomers), 1:size(net%isomers) ) )
	
	do u = 1, size(net%Tlist)
		
		T = net%Tlist(u)
		
		! Calculate the equilibrium (Boltzmann) distributions for each (unimolecular) well
		do i = 1, size(net%isomers)
			if (size(net%isomers(i)%species) == 1) then
				call eqDist(net%isomers(i), net%E, T)
			end if
		end do
		
		! Calculate microcanonical rate coefficients at the current conditions
		do w = 1, size(net%reactions)
			call microRates(net%reactions(w), T, net%E, net%isomers, net%species)
		end do
		
		do v = 1, size(net%Plist)

			P = net%Plist(v)

			! Calculate collision frequencies
			do i = 1, size(net%isomers)
				if (size(net%isomers(i)%species) == 1) then
					net%isomers(i)%collFreq = collisionFrequency(T, P, net%bathGas, &
						net%species(net%isomers(i)%species(1))%general)
				end if
			end do
			
			! Determine phenomenological rate coefficients
			if (net%method == 1) then
				call ssmscRates(net, T, P, K(u,v,:,:))
			elseif (net%method == 2) then
				call ssrsRates(net, T, P, K(u,v,:,:))
			end if
			
			! Check for validity of phenomenological rate coefficients
			do i = 1, size(net%isomers)
				do j = 1, size(net%isomers)
					if (i /= j) then
						if (K(u,v,i,j) <= 0) then
							write (0, fmt='(A)') 'One ore more rate coefficients not properly estimated.'
							stop
						end if
					end if
				end do
			end do
			
			! Convert units of bimolecular k(T, P) to cm^3/mol*s
			do j = 1, size(net%isomers)
				if (size(net%isomers(j)%species) > 1) then
					do i = 1, size(net%isomers)
						K(u,v,i,j) = K(u,v,i,j) * 1.0e6
					end do
				end if
			end do
			
			
		end do
		
	end do
	
	! Fit phenomenological rate coefficients to selected model
	if (net%model == 1) then
		allocate( chebyshevCoefficients(net%numChebT, net%numChebP, size(net%isomers), size(net%isomers)) )
	elseif (net%model == 2) then
		allocate( pDepArrhenius(size(net%Plist), size(net%isomers), size(net%isomers)) )
	end if
	do i = 1, size(net%isomers)
		do j = 1, size(net%isomers)
			if (i /= j) then
				if (net%model == 1) then
					call fitChebyshevModel(K(:,:,i,j), net%Tlist, net%Plist, net%numChebT, net%numChebP, chebyshevCoefficients(:,:,i,j))
				elseif (net%model == 2) then
					call fitPDepArrheniusModel(K(:,:,i,j), net%Tlist, net%Plist, pDepArrhenius(:,i,j))
				end if
			end if
		end do
	end do
	
	! Write results to stdout
	call saveNetwork(net, K, chebyshevCoefficients, pDepArrhenius)
	
end program

! ==============================================================================
