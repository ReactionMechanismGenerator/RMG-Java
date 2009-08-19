!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	reservoirState.f90
!
!	Written by Josh Allen (jwallen@mit.edu)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module ReservoirStateModule

	use SpeciesModule
	use IsomerModule
	use ReactionModule
	use NetworkModule
	
	implicit none
	
contains

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!
	! Subroutine: ssrsRates()
	! 
	! Uses the steady state/reservoir state method of Green and Bhatti to
	! estimate the phenomenological rate coefficients from a full ME matrix.
	!
	subroutine ssrsRates(net, T, P, K)

		! Provide parameter type checking of inputs and outputs
		type(Network), intent(inout)				:: 	net
		real(8), intent(in)							::	T
		real(8), intent(in)							::	P
		real(8), dimension(:,:), intent(out) 		:: 	K
		
		! Number of reservoir and active-state energy grains for each isomer
		integer, dimension(:), allocatable			:: 	nRes, nAct

		! Number of unimolecular wells
		integer nUni, row
		
		! Accounting matrix
		integer, dimension(:,:), allocatable		::	indices
		
		! Active state grain matrix and RHS vectors
		real(8), dimension(:,:), allocatable	:: 	L
		real(8), dimension(:,:), allocatable	:: 	Z
		
		! Pseudo-steady state grain populations
		real(8), dimension(:,:,:), allocatable	:: 	pa
		
		! Variables for LAPACK
		integer, dimension(:), allocatable				::	iPiv
		integer											::	info
		
		! Indices
		integer i, n, r, s, reac, prod
		
		! Determine number of unimolecular wells
		nUni = 0
		do i = 1, size(net%isomers)
			if (size(net%isomers(i)%species) == 1) then
				nUni = nUni + 1
			end if
		end do
		
		! Collisional terms in master equation
		do i = 1, size(net%isomers)
			if (size(net%isomers(i)%species) == 1) then
				if (.not. allocated(net%isomers(i)%Mcoll)) then
					allocate( net%isomers(i)%Mcoll(1:size(net%E), 1:size(net%E)) )
				end if
				net%isomers(i)%Mcoll(:,:) = collision(net, T, P, net%isomers(i))
			end if
		end do
		
		! Determine reservoir cutoff grains for each unimolecular isomer
		allocate( nRes(1:size(net%isomers)), nAct(1:size(net%isomers)) )
		nRes = reservoirCutoffs(net)
		do i = 1, size(net%isomers)
			nAct(i) = size(net%E) - nRes(i)
		end do

		! Determine pseudo-steady state populations of active state
		allocate( pa(1:size(net%E), 1:size(net%isomers), 1:nUni) )
		pa = 0 * pa
		if (nUni == 1) then
			! Not worth it to call banded solver when only one well
			call activeStateFull(net, nUni, nRes, nAct, pa)
		else
			! Very worth it to call banded solver when more than one well
			call activeStateBanded(net, nUni, nRes, nAct, pa)
		end if
		
		! Check that PSSA populations are all nonnegative; fail if not
		do i = 1, nUni
			do n = 1, size(net%isomers)
				do r = 1, size(net%E)
					if (pa(r,n,i) < 0.0) then
						write (*,*) 'ERROR: Negative steady-state populations encountered during ReservoirState method.'
						stop
					end if
				end do
			end do
		end do
		
		! Initialize phenomenological rate coefficient matrix
		do i = 1, size(net%isomers)
			do n = 1, size(net%isomers)
				K(i,n) = 0
			end do
		end do
		
		! Determine phenomenological rate coefficients
		do i = 1, nUni
			do n = 1, size(net%isomers)
				do r = 1, nRes(i)
					K(i,n) = K(i,n) + sum(net%isomers(i)%Mcoll(r,:) * pa(:,n,i))
				end do
			end do
		end do
		do r = 1, size(net%reactions)
			reac = net%reactions(r)%reac
			prod = net%reactions(r)%prod
			if  (reac <= nUni .and. prod > nUni) then	! Dissociation
				do n = 1, size(net%isomers)
					K(prod,n) = sum(net%reactions(r)%kf * pa(:,n,reac))
				end do
			elseif  (reac > nUni .and. prod <= nUni) then	! Combination
				do n = 1, size(net%isomers)
					K(reac,n) = sum(net%reactions(r)%kb * pa(:,n,prod))
				end do
			end if
		end do
		do n = 1, size(net%isomers)
			K(n,n) = 0.
		end do
        
        ! Clean up
		deallocate( nRes, nAct, pa )

	end subroutine

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!
	! Subroutine: activeStateFull()
	! 
	! Determine the pseudo-steady state populations for the active state
	! grains using a full matrix linear solve.
	!
	! Parameters:
	!
	subroutine activeStateFull(net, nUni, nRes, nAct, pa)
	
		! Provide parameter type checking of inputs and outputs
		type(Network), intent(in)					:: 	net
		integer, intent(in)							:: 	nUni
		integer, dimension(:), intent(in)			:: 	nRes
		integer, dimension(:), intent(in)			:: 	nAct
		real(8), dimension(:,:,:), intent(out) 		:: 	pa
		
		!! Accounting matrix
		integer, dimension(:,:), allocatable		::	indices
		
		! Active state grain matrix and RHS vectors
		real(8), dimension(:,:), allocatable	:: 	L
		real(8), dimension(:,:), allocatable	:: 	Z
		
		! Variables for LAPACK
		integer, dimension(:), allocatable				::	iPiv
		integer											::	info
		
		! Indices
		integer i, n, r, s, reac, prod
		
		! Construct accounting matrix
		! Row is grain number, column is well number, value is index into active-state matrix
		allocate( indices(1:size(net%E), 1:nUni ) )
		call accountingMatrix(size(net%E), nUni, nRes, indices)
		
		! Create and zero active-state matrix and RHS vectors
		allocate( L(1:sum(nAct), 1:sum(nAct)), Z(1:sum(nAct), size(net%isomers)) )
		L = 0 * L
		Z = 0 * Z
		
		! Collisional terms in active-state matrix and RHS vectors
		do i = 1, nUni
			do r = nRes(i)+1, size(net%E)
				do s = nRes(i)+1, size(net%E)
					L(indices(r,i), indices(s,i)) = net%isomers(i)%Mcoll(r,s)
				end do
				Z(indices(r,i), i) = sum(net%isomers(i)%Mcoll(r,1:nRes(i)) * &
					net%isomers(i)%eqDist(1:nRes(i)))
			end do
		end do
		
		! Reactive terms in active-state matrix and RHS vectors
		do n = 1, size(net%reactions)
			reac = net%reactions(n)%reac
			prod = net%reactions(n)%prod
			if (reac <= nUni .and. prod <= nUni) then		! Isomerization
				do r = max(nRes(reac), nRes(prod))+1, size(net%E)
					L(indices(r,reac), indices(r,prod)) = net%reactions(n)%kb(r)
					L(indices(r,reac), indices(r,reac)) = L(indices(r,reac), indices(r,reac)) - net%reactions(n)%kf(r)
					L(indices(r,prod), indices(r,reac)) = net%reactions(n)%kf(r)
					L(indices(r,prod), indices(r,prod)) = L(indices(r,prod), indices(r,prod)) - net%reactions(n)%kb(r)
				end do
			elseif  (reac <= nUni .and. prod > nUni) then	! Dissociation
				do r = nRes(reac)+1, size(net%E)
					L(indices(r,reac), indices(r,reac)) = L(indices(r,reac), indices(r,reac)) - net%reactions(n)%kf(r)
					Z(indices(r,reac), prod) = net%reactions(n)%kb(r)
				end do
			elseif  (reac > nUni .and. prod <= nUni) then	! Combination
				do r = nRes(prod)+1, size(net%E)
					L(indices(r,prod), indices(r,prod)) = L(indices(r,prod), indices(r,prod)) - net%reactions(n)%kb(r)
					Z(indices(r,prod), reac) = net%reactions(n)%kf(r)
				end do
			end if
		end do
		
		Z = -Z
		
		! Solve for pseudo-steady state populations of active state
		allocate( iPiv(1:sum(nAct)) )
		call DGESV(sum(nAct), size(net%isomers), L, sum(nAct), iPiv, Z, sum(nAct), info)
		if (info /= 0) then
			write (*,*) 'ERROR: Active-state matrix is singular!'
			stop
		end if
		deallocate( iPiv )
		
		! Convert solution to pseudo-steady state populations
		do r = minval(nRes)+1, size(net%E)
			do n = 1, size(net%isomers)
				do i = 1, nUni
					if (indices(r,i) > 0) pa(r,n,i) = Z(indices(r,i), n)
				end do
			end do
		end do
		
		! Clean up
		deallocate( indices, L, Z )
		
		
	end subroutine
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!
	! Subroutine: activeStateBanded()
	! 
	! Determine the pseudo-steady state populations for the active state
	! grains using a full matrix linear solve.
	!
	! Parameters:
	!
	subroutine activeStateBanded(net, nUni, nRes, nAct, pa)
	
		! Provide parameter type checking of inputs and outputs
		type(Network), intent(in)					:: 	net
		integer, intent(in)							:: 	nUni
		integer, dimension(:), intent(in)			:: 	nRes
		integer, dimension(:), intent(in)			:: 	nAct
		real(8), dimension(:,:,:), intent(out) 		:: 	pa
		
		!! Accounting matrix
		integer, dimension(:,:), allocatable		::	indices
		
		! Active state grain matrix and RHS vectors
		real(8), dimension(:,:), allocatable		:: 	L
		real(8), dimension(:,:), allocatable		:: 	Z
		
		integer										::	bandwidth, halfbandwidth
		real(8)										::	dE
		
		! Variables for LAPACK
		integer, dimension(:), allocatable			::	iPiv
		integer										::	info
		
		! Indices
		integer i, n, r, s, reac, prod
		
		! Construct accounting matrix
		! Row is grain number, column is well number, value is index into active-state matrix
		allocate( indices(1:size(net%E), 1:nUni ) )
		indices = 0 * indices
		call accountingMatrix(size(net%E), nUni, nRes, indices)
		
		! Determine bandwidth (at which transfer probabilities are so low that they can be truncated
		! with negligible error)
		dE = net%E(2) - net%E(1)
		halfbandwidth = getHalfBandwidth(net%collisionParameters(1), dE) * nUni
		bandwidth = 2 * halfbandwidth + 1
		
		! Create and zero active-state matrix and RHS vectors
		allocate( L(3 * halfbandwidth + 1, 1:sum(nAct)), Z(1:sum(nAct), 1:size(net%isomers)) )
		L = 0 * L
		Z = 0 * Z
		
		! Collisional terms in active-state matrix and RHS vectors
		do i = 1, nUni
			do s = nRes(i)+1, size(net%E)
				do r = max(nRes(i)+1, s - halfbandwidth), min(size(net%E), s + halfbandwidth)
					L(bandwidth + indices(r,i) - indices(s,i), indices(s,i)) = net%isomers(i)%Mcoll(r,s)
				end do
				Z(indices(s,i), i) = sum(net%isomers(i)%Mcoll(s,1:nRes(i)) * &
					net%isomers(i)%eqDist(1:nRes(i)))
			end do
		end do
		
		! Reactive terms in active-state matrix and RHS vectors
		do n = 1, size(net%reactions)
			reac = net%reactions(n)%reac
			prod = net%reactions(n)%prod
			if (reac <= nUni .and. prod <= nUni) then		! Isomerization
				do r = max(nRes(reac), nRes(prod))+1, size(net%E)
					L(bandwidth + indices(r,reac) - indices(r,prod), indices(r,prod)) = net%reactions(n)%kb(r)
					L(bandwidth, indices(r,prod)) = L(bandwidth, indices(r,prod)) - net%reactions(n)%kb(r)
					L(bandwidth + indices(r,prod) - indices(r,reac), indices(r,reac)) = net%reactions(n)%kf(r)
					L(bandwidth, indices(r,reac)) = L(bandwidth, indices(r,reac)) - net%reactions(n)%kf(r)
				end do
			elseif  (reac <= nUni .and. prod > nUni) then	! Dissociation
				do r = nRes(reac)+1, size(net%E)
					L(bandwidth, indices(r,reac)) = L(bandwidth, indices(r,reac)) - net%reactions(n)%kf(r)
					Z(indices(r,reac), prod) = net%reactions(n)%kb(r)
				end do
			elseif  (reac > nUni .and. prod <= nUni) then	! Combination
				do r = nRes(prod)+1, size(net%E)
					L(bandwidth, indices(r,prod)) = L(bandwidth, indices(r,prod)) - net%reactions(n)%kb(r)
					Z(indices(r,prod), reac) = net%reactions(n)%kf(r)
				end do
			end if
		end do
		
		Z = -Z
		
		! Solve for pseudo-steady state populations of active state
		allocate( iPiv(1:sum(nAct)) )
		call DGBSV(sum(nAct), halfbandwidth, halfbandwidth, size(net%isomers), &
			L, 3 * halfbandwidth + 1, iPiv, Z, sum(nAct), info)
		if (info /= 0) then
			write (*,*) 'ERROR: Active-state matrix is singular!'
			stop
		end if
		deallocate( iPiv )
		
		! Convert solution to pseudo-steady state populations
		do r = minval(nRes)+1, size(net%E)
			do n = 1, size(net%isomers)
				do i = 1, nUni
					if (indices(r,i) > 0) pa(r,n,i) = Z(indices(r,i), n)
				end do
			end do
		end do
		
		! Clean up
		deallocate( indices, L, Z )
		
		
	end subroutine

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!
	! Subroutine: accountingMatrix()
	! 
	! Determines the grain below which the reservoir approximation will be used
	! and above which the pseudo-steady state approximation will be used by
	! examining the energies of the transition states connected to each 
	! unimolecular isomer.
	!
	! Parameters:
	!   simData - The simulation parameters.
	!   isomerList - The chemical data about each isomer.
	!   rxnList - The chemical data about each transition state.
	!
	! Returns:
	!	nRes - The reservoir cutoff grains for each unimolecular isomer.
	subroutine accountingMatrix(nGrains, nUni, nRes, indices)
	
		integer, intent(in)						::	nGrains
		integer, intent(in)						::	nUni
		integer, dimension(:), intent(in)		::	nRes
		integer, dimension(:,:), intent(out)	::	indices
		
		integer i, r, row
		
		! Construct accounting matrix
		! Row is grain number, column is well number, value is index into active-state matrix
		row = 1
		do r = minval(nRes)+1, nGrains
			do i = 1, nUni
				if (r > nRes(i)) then
					indices(r,i) = row
					row = row + 1
				end if
			end do
		end do
		
	
	end subroutine
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!
	! Subroutine: reservoirCutoffs()
	! 
	! Determines the grain below which the reservoir approximation will be used
	! and above which the pseudo-steady state approximation will be used by
	! examining the energies of the transition states connected to each 
	! unimolecular isomer.
	!
	! Parameters:
	!   simData - The simulation parameters.
	!   isomerList - The chemical data about each isomer.
	!   rxnList - The chemical data about each transition state.
	!
	! Returns:
	!	nRes - The reservoir cutoff grains for each unimolecular isomer.
	function reservoirCutoffs(net)
	
		type(Network), intent(in)					:: 	net
		integer, dimension(1:size(net%isomers))		::	reservoirCutoffs
		
		real(8) Eres, Emin, Emax, dE
		integer i, j, t, start
		
		Emin = minval(net%E)
		Emax = maxval(net%E)
		dE = net%E(2) - net%E(1)
		
		! Determine reservoir cutoffs by looking at transition state energies
		do i = 1, size(net%isomers)
			if (size(net%isomers(i)%species) == 1) then
			
				start = ceiling((net%isomers(i)%E0 - Emin) / dE) + 1
				Eres = Emax
				
				do t = 1, size(net%reactions)
					if (net%reactions(t)%reac == i .or. net%reactions(t)%prod == i) then
						if (net%reactions(t)%E0 < Eres) Eres = net%reactions(t)%E0
					end if
				end do
				
				reservoirCutoffs(i) = floor((Eres - 10 * net%collisionParameters(1) - Emin) / dE)
				
				! Make sure reservoirCutoffs(i) is valid (i.e. points to grain at or above ground state)
				if (reservoirCutoffs(i) < start) then
					reservoirCutoffs(i) = start
				end if
				
			else
				reservoirCutoffs(i) = size(net%E)
			end if
		end do

	end function
	
	

end module
