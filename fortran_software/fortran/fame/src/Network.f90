!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	network.f90
!
!	Written by Josh Allen (jwallen@mit.edu)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module NetworkModule

	use SpeciesModule
	use IsomerModule
	use ReactionModule
	use StatesModule
	
	implicit none
	
	type Network
		integer										:: method
		real(8), dimension(:), allocatable			:: Tlist
		real(8), dimension(:), allocatable			:: Plist
		real(8), dimension(:), allocatable			:: E
		integer										:: numGrains
		real(8)										:: grainSize
		integer										:: model
		integer										:: numChebT
		integer										:: numChebP
		integer										:: collisionModel
		real(8), dimension(:), allocatable			:: collisionParameters
		type(GeneralData)							:: bathGas
		type(Species), dimension(:), allocatable	:: species
		type(Isomer), dimension(:), allocatable		:: isomers
		type(Reaction), dimension(:), allocatable	:: reactions
		
	end type

contains

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!
	! Subroutine: selectEnergyGrains()
	!
	! Determines a suitable set of energy grains to use for the calculations.
	! The minimum is chosen to be just below the lowest ground-state energy in
	! the system. The maximum is chosen such that the equilibrium distribution
	! of the highest-energy isomer has tailed off to negligible levels.
	!
	subroutine selectEnergyGrains(net)
	
		type(Network), intent(inout)		:: 	net
		
		real(8) Tmax, Emax0, mult, value, tol, Emin, Emax, dE
		integer isom, i, r, done, maxindex, nE
		
		! For the purposes of finding the maximum energy we will use 201 grains
		nE = 201
		dE = 0.0
		
		! Use the maximum temperature because that is where equilibrium distribution is widest
		Tmax = maxval(net%Tlist)
		
		! Determine minimum energy and isomer with maximum ground-state energy
		Emin = net%isomers(1)%E0
		isom = 1
		do i = 1, size(net%isomers)
			if (net%isomers(i)%E0 < Emin) then
				Emin = net%isomers(i)%E0
			elseif (net%isomers(i)%E0 > Emax0) then
				Emax0 = net%isomers(i)%E0
				isom = i
			end if
		end do
		
		! Round Emin down to nearest integer
		Emin = floor(Emin)
	
		! (Try to) purposely overestimate Emax using arbitrary multiplier
		! This is to (hopefully) avoid multiple density of states calculations
		mult = 50
		done = 0
		do while (done == 0 .and. mult <= 250)
		
			Emax = ceiling(Emax0 + mult * 8.314472 * Tmax)
			
			call setEnergyGrains(net, Emin, Emax, dE, nE)
			call densityOfStates(net%isomers(isom), net%species, net%E)
			call eqDist(net%isomers(isom), net%E, Tmax)
			
			! Find maximum of distribution
			maxindex = 0
			value = 0.0
			do r = 1, nE
				if (net%isomers(isom)%eqDist(r) > value) then
					value = net%isomers(isom)%eqDist(r)
					maxindex = r
				end if
			end do
			
			tol = 1e-4
			if (maxindex == 0) then
				! Couldn't find a maximum; this suggests an error in the
				! density of states
				write (*,*) 'ERROR: Unable to determine maximum of equilibrium distribution.'
				write (*,*) 'Density of states is:'
				write (*,*) net%isomers(isom)%densStates
				write (*,*) 'Equilibrium distribution is:'
				write (*,*) net%isomers(isom)%eqDist
				deallocate( net%E )
				deallocate( net%isomers(isom)%densStates )
				deallocate( net%isomers(isom)%eqDist )
				stop
			elseif (net%isomers(isom)%eqDist(nE) / value < tol) then
				! If tail of distribution is much lower than the maximum, then we've found bounds for Emax
				r = nE - 1
				do while (r > 0 .and. done == 0)
					if (net%isomers(isom)%eqDist(r) / value > tol) then
						done = 1
					else
						r = r - 1
					end if
				end do
				Emax = net%E(r)
			else
				mult = mult + 50
			end if
			
			deallocate( net%E )
			deallocate( net%isomers(isom)%densStates )
			deallocate( net%isomers(isom)%eqDist )
			
		end do
		
		if (mult > 250) then
			Emax = Emax0 + 100 * 8.314472 * Tmax
		end if
		
		! Round Emax up to nearest integer
		Emax = ceiling(Emax)		
		
		! Set energy range for good
		call setEnergyGrains(net, Emin, Emax, net%grainSize, net%numGrains)
		
	end subroutine
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!
	! Subroutine: setEnergyGrains()
	!
	! Set the energy grains to use for the calculations.
	!
	subroutine setEnergyGrains(net, Emin, Emax, dE0, nGrains0)
	
	    type(Network), intent(inout)				:: 	net
		real(8), intent(in)							::	Emin
		real(8), intent(in)							::	Emax
		real(8), intent(in)							::	dE0
		integer, intent(in)							::	nGrains0
		
		real(8) dE
		integer i, nGrains
		
		if (nGrains0 == 0) then                      ! nGrains is undefined - assume nGrains0 is desired
			nGrains = aint((Emax - Emin) / dE0) + 1  ! dE is undefined - assume dE0 is desired
			dE = dE0
		else
			nGrains = nGrains0
			dE = (Emax - Emin) / (nGrains - 1)
		end if
		write (*,*) '#DEBUG: In setEnergyGrains, Emax=',Emax,'Emin=',Emin
		write (*,*) '#DEBUG: In setEnergyGrains, nGrains0=',nGrains0,'nGrains=',nGrains
		write (*,*) '#DEBUG: In setEnergyGrains, dE0=',dE0,'dE=',dE
		
		if (allocated(net%E)) then
			deallocate(net%E)
		end if
		
		allocate (net%E(1:nGrains))
        do i = 1, nGrains
            net%E(i) = dE * (i - 1) + Emin
        end do
		
		write (*,*) '#DEBUG: at end of setEnergyGrains, net%E is',size(net%E),'long and starts:',net%E(1:10)
	
	end subroutine
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!
	! Function: collision()
	!
	! Calculates the collisional transfer rate matrix for a given isomer.
	!
	function collision(net, T, P, isom)

		type(Network), intent(in)					:: 	net
		real(8), intent(in)							::	T
		real(8), intent(in)							::	P
		type(Isomer), intent(in)					::	isom
		real(8), dimension(1:size(net%E),1:size(net%E)) :: collision
		
		real(8)			:: gasConc			! Gas concentration in molecules/m^3
		
		type(Species)	:: spec
		
		real(8), dimension(:,:), allocatable		::	prob
		real(8), dimension(:), allocatable		::	C		! Vector of normalization coefficients
		
		real(8) Emin, Emax, dE
		
		integer			:: bandwidth, halfbandwidth
		integer			:: lb, ub, start
		integer			:: r, s
		
		! Determine (ideal) gas concentration
		gasConc = P / 1.381e-23 / T
		
		! Allocate collision matrix
		collision = 0 * collision
		
		! Allocate/zero vector of normalization coefficients
		allocate( C(1:size(net%E)) )
		C = 0 * C
		
		Emin = minval(net%E)
		Emax = maxval(net%E)
		dE = net%E(2) - net%E(1)
		
		! Get species of isomer (isom should be unimolecular!)
		spec = net%species(isom%species(1))
		
		! Determine bandwidth (at which transfer probabilities are so low that they can be truncated
		! with negligible error)
		halfbandwidth = getHalfBandwidth(net%collisionParameters(1), dE)
		bandwidth = 2 * halfbandwidth + 1
		
		! Determine start grain (corresponding to isomer ground-state energy)
		start = 0
		do r = 1, size(net%E)
			if (start == 0 .and. isom%densStates(r) > 0) then
				start = r
			end if
		end do
		if (start == 0) then
			write (0, fmt='(A)') 'Unable to determine starting energy grain for isomer collision matrix calculation.'
			stop
		end if
		
		! Determine unnormalized entries in collisional tranfer probability matrix for the current isomer
		do r = 1, size(net%E)
			do s = 1, size(net%E)
				lb = max(r - halfbandwidth, start)
				ub = min(r + halfbandwidth, size(net%E))
				if (r >= start .and. s >= lb .and. s <= ub) then
					collision(s,r) = transferRate(s, r, net%E(s), net%E(r), net%collisionParameters(1), &
						isom%E0, isom%densStates, T)
				else
					collision(s,r) = 0
				end if
			end do
		end do
		
		! Normalize using detailed balance
		do r = size(net%E), start, -1
			C(r) = (1 - sum(C(r+1:size(net%E)) * collision(r+1:size(net%E),r))) / sum(collision(1:r,r))
		end do
		!do r = start, size(net%E)
		!	C(r) = (1 - sum(C(start:r) * collision(start:r,r))) / sum(collision(r:size(net%E),r))
		!end do
		
		! Check for normalization consistency (i.e. all numbers are positive)
		!do r = start, size(net%E)
		!	if (C(r) <= 0) then
		!		write (*,*) 'Error normalizing collisional transfer probabilities matrix!'
		!		stop
		!	end if
		!end do
			
		do r = start, size(net%E)
			collision(r,1:r-1) = collision(r,1:r-1) * C(r)
			collision(1:r-1,r) = collision(1:r-1,r) * C(r)
			collision(r,r) = collision(r,r) * C(r) - 1
		end do
		
		! Multiply by collision frequency to determine collision rates
		collision = isom%collFreq * collision
		
		deallocate( C )
		
	end function
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!
	! Function: getHalfBandwidth
	!
	! Calculate the number of collisional transfer probabilities to include
	! before truncation.
	!
	function getHalfBandwidth(alpha, dE)
	
		real(8), intent(in) 		::	alpha
		real(8), intent(in) 		::	dE
		integer						::	getHalfBandwidth
	
		getHalfBandwidth = ceiling(16 * alpha / dE)
		
	end function
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!
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

		R = 8.314472
		
		! Evaluate collisional transfer probability - theoretically correct way
		!if (Ej < E0 .or. Ei < E0 .or. rho(j) .eq. 0) then
		!	transferRate = 0.
		!elseif (Ej >= Ei) then
		!	transferRate = exp(-(Ej - Ei) / alpha)
		!else
		!	transferRate = exp(-(Ei - Ej) / alpha) * rho(i) / rho(j) * exp( -(Ei - Ej) / (R * T))
		!end if
		
		! Evaluate collisional transfer probability - move density of states ratio to other side
		! so that it is always less than 1
		if (Ej < E0 .or. Ei < E0 .or. rho(j) .eq. 0) then
			transferRate = 0.
		elseif (Ej >= Ei) then
			transferRate = exp(-(Ej - Ei) / alpha) * rho(i) / rho(j)
		else
			transferRate = exp(-(Ei - Ej) / alpha) * exp( -(Ei - Ej) / (R * T))
		end if

	end function

end module
