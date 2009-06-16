! ==============================================================================
!
! 	ReservoirState.f90
!
! 	Written by Josh Allen (jwallen@mit.edu)
!
! ==============================================================================

module ReservoirStateModule

	use SimulationModule
	use IsomerModule
	use ReactionModule

	implicit none
	
contains

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!
	! Subroutine: ssrsRates()
	! 
	! Uses the steady state/reservoir state method of Green and Bhatti to
	! estimate the phenomenological rate coefficients from a full ME matrix.
	!
	subroutine ssrsRates(T, P, simData, speciesList, isomerList, rxnList, K)

		! Provide parameter type checking of inputs and outputs
		real(8), intent(in)							::	T
		real(8), intent(in)							::	P
		type(Simulation), intent(in)				:: 	simData
		type(Species), dimension(:), intent(in)		:: 	speciesList
		type(Isomer), dimension(:), intent(inout)   ::  isomerList
        type(Reaction), dimension(:), intent(in)	:: 	rxnList
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
		do i = 1, simData%nIsom
			if (isomerList(i)%numSpecies == 1) then
				nUni = nUni + 1
			end if
		end do
		
		! Collisional terms in master equation
		do i = 1, simData%nIsom
			if (isomerList(i)%numSpecies == 1) then
				call collision(T, P, simData, isomerList(i), speciesList)
			end if
		end do
		
		! Determine reservoir cutoff grains for each unimolecular isomer
		allocate( nRes(1:simData%nIsom), nAct(1:simData%nIsom) )
		nRes = reservoirCutoffs(simData, isomerList, rxnList)
		do i = 1, simData%nIsom
			nAct(i) = simData%nGrains - nRes(i)
		end do

		! Determine pseudo-steady state populations of active state
		allocate( pa(1:simData%nGrains, 1:simData%nIsom, 1:nUni) )
		pa = 0 * pa
		!call activeStateFull(simData, isomerList, rxnList, nUni, nRes, nAct, pa)
		call activeStateBanded(simData, isomerList, rxnList, nUni, nRes, nAct, pa)
		
		! Check that PSSA populations are all nonnegative; fail if not
		do i = 1, nUni
			do n = 1, simData%nIsom
				do r = 1, simData%nGrains
					if (pa(r,n,i) < 0.0) then
						write (*,*) 'ERROR: Negative steady-state populations encountered during ReservoirState method.'
						stop
					end if
				end do
			end do
		end do
		
		
		! Initialize phenomenological rate coefficient matrix
		do i = 1, simData%nIsom
			do n = 1, simData%nIsom
				K(i,n) = 0
			end do
		end do
		
		! Determine phenomenological rate coefficients
		do i = 1, nUni
			do n = 1, simData%nIsom
				do r = 1, nRes(i)
					K(i,n) = K(i,n) + sum(isomerList(i)%Mcoll(r,:) * pa(:,n,i))
				end do
			end do
		end do
		do r = 1, simData%nRxn
			reac = rxnList(r)%isomerList(1)
			prod = rxnList(r)%isomerList(2)
			if  (reac <= nUni .and. prod > nUni) then	! Dissociation
				do n = 1, simData%nIsom
					K(prod,n) = sum(rxnList(r)%kf * pa(:,n,reac))
				end do
			elseif  (reac > nUni .and. prod <= nUni) then	! Combination
				do n = 1, simData%nIsom
					K(reac,n) = sum(rxnList(r)%kb * pa(:,n,prod))
				end do
			end if
		end do
		do n = 1, simData%nIsom
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
	subroutine activeStateFull(simData, isomerList, rxnList, nUni, nRes, nAct, pa)
	
		! Provide parameter type checking of inputs and outputs
		type(Simulation), intent(in)				:: 	simData
		type(Isomer), dimension(:), intent(in)   	::  isomerList
        type(Reaction), dimension(:), intent(in)	:: 	rxnList
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
		allocate( indices(1:simData%nGrains, 1:nUni ) )
		call accountingMatrix(simData%nGrains, nUni, nRes, indices)
		
		! Create and zero active-state matrix and RHS vectors
		allocate( L(1:sum(nAct), 1:sum(nAct)), Z(1:sum(nAct), 1:simData%nIsom) )
		L = 0 * L
		Z = 0 * Z
		
		! Collisional terms in active-state matrix and RHS vectors
		do i = 1, nUni
			do r = nRes(i)+1, simData%nGrains
				do s = nRes(i)+1, simData%nGrains
					L(indices(r,i), indices(s,i)) = isomerList(i)%Mcoll(r,s)
				end do
				Z(indices(r,i), i) = sum(isomerList(i)%Mcoll(r,1:nRes(i)) * &
					isomerList(i)%eqDist(1:nRes(i)))
			end do
		end do
		
		! Reactive terms in active-state matrix and RHS vectors
		do n = 1, simData%nRxn
			reac = rxnList(n)%isomerList(1)
			prod = rxnList(n)%isomerList(2)
			if (reac <= nUni .and. prod <= nUni) then		! Isomerization
				do r = max(nRes(reac), nRes(prod))+1, simData%nGrains
					L(indices(r,reac), indices(r,prod)) = rxnList(n)%kb(r)
					L(indices(r,reac), indices(r,reac)) = L(indices(r,reac), indices(r,reac)) - rxnList(n)%kf(r)
					L(indices(r,prod), indices(r,reac)) = rxnList(n)%kf(r)
					L(indices(r,prod), indices(r,prod)) = L(indices(r,prod), indices(r,prod)) - rxnList(n)%kb(r)
				end do
			elseif  (reac <= nUni .and. prod > nUni) then	! Dissociation
				do r = nRes(reac)+1, simData%nGrains
					L(indices(r,reac), indices(r,reac)) = L(indices(r,reac), indices(r,reac)) - rxnList(n)%kf(r)
					Z(indices(r,reac), prod) = rxnList(n)%kb(r)
				end do
			elseif  (reac > nUni .and. prod <= nUni) then	! Combination
				do r = nRes(prod)+1, simData%nGrains
					L(indices(r,prod), indices(r,prod)) = L(indices(r,prod), indices(r,prod)) - rxnList(n)%kb(r)
					Z(indices(r,prod), reac) = rxnList(n)%kf(r)
				end do
			end if
		end do
		
		Z = -Z
		
		! Solve for pseudo-steady state populations of active state
		allocate( iPiv(1:sum(nAct)) )
		call DGESV(sum(nAct), simData%nIsom, L, sum(nAct), iPiv, Z, sum(nAct), info)
		if (info /= 0) then
			write (*,*) 'ERROR: Active-state matrix is singular!'
			stop
		end if
		deallocate( iPiv )
		
		! Convert solution to pseudo-steady state populations
		do r = minval(nRes)+1, simData%nGrains
			do n = 1, simData%nIsom
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
	subroutine activeStateBanded(simData, isomerList, rxnList, nUni, nRes, nAct, pa)
	
		! Provide parameter type checking of inputs and outputs
		type(Simulation), intent(in)				:: 	simData
		type(Isomer), dimension(:), intent(in)   	::  isomerList
        type(Reaction), dimension(:), intent(in)	:: 	rxnList
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
		
		! Variables for LAPACK
		integer, dimension(:), allocatable			::	iPiv
		integer										::	info
		
		! Indices
		integer i, n, r, s, reac, prod
		
		! Construct accounting matrix
		! Row is grain number, column is well number, value is index into active-state matrix
		allocate( indices(1:simData%nGrains, 1:nUni ) )
		indices = 0 * indices
		call accountingMatrix(simData%nGrains, nUni, nRes, indices)
		
		! Determine bandwidth (at which transfer probabilities are so low that they can be truncated
		! with negligible error)
		halfbandwidth = getHalfBandwidth(simData%alpha, simData%dE) * nUni
		bandwidth = 2 * halfbandwidth + 1
		
		! Create and zero active-state matrix and RHS vectors
		allocate( L(3 * halfbandwidth + 1, 1:sum(nAct)), Z(1:sum(nAct), 1:simData%nIsom) )
		L = 0 * L
		Z = 0 * Z
		
		! Collisional terms in active-state matrix and RHS vectors
		do i = 1, nUni
			do s = nRes(i)+1, simData%nGrains
				do r = max(nRes(i)+1, s - halfbandwidth), min(simData%nGrains, s + halfbandwidth)
					L(bandwidth + indices(r,i) - indices(s,i), indices(s,i)) = isomerList(i)%Mcoll(r,s)
				end do
				Z(indices(s,i), i) = sum(isomerList(i)%Mcoll(s,1:nRes(i)) * &
					isomerList(i)%eqDist(1:nRes(i)))
			end do
		end do
		
		! Reactive terms in active-state matrix and RHS vectors
		do n = 1, simData%nRxn
			reac = rxnList(n)%isomerList(1)
			prod = rxnList(n)%isomerList(2)
			if (reac <= nUni .and. prod <= nUni) then		! Isomerization
				do r = max(nRes(reac), nRes(prod))+1, simData%nGrains
					L(bandwidth + indices(r,reac) - indices(r,prod), indices(r,prod)) = rxnList(n)%kb(r)
					L(bandwidth, indices(r,prod)) = L(bandwidth, indices(r,prod)) - rxnList(n)%kb(r)
					L(bandwidth + indices(r,prod) - indices(r,reac), indices(r,reac)) = rxnList(n)%kf(r)
					L(bandwidth, indices(r,reac)) = L(bandwidth, indices(r,reac)) - rxnList(n)%kf(r)
				end do
			elseif  (reac <= nUni .and. prod > nUni) then	! Dissociation
				do r = nRes(reac)+1, simData%nGrains
					L(bandwidth, indices(r,reac)) = L(bandwidth, indices(r,reac)) - rxnList(n)%kf(r)
					Z(indices(r,reac), prod) = rxnList(n)%kb(r)
				end do
			elseif  (reac > nUni .and. prod <= nUni) then	! Combination
				do r = nRes(prod)+1, simData%nGrains
					L(bandwidth, indices(r,prod)) = L(bandwidth, indices(r,prod)) - rxnList(n)%kb(r)
					Z(indices(r,prod), reac) = rxnList(n)%kf(r)
				end do
			end if
		end do
		
		Z = -Z
		
		! Solve for pseudo-steady state populations of active state
		allocate( iPiv(1:sum(nAct)) )
		call DGBSV(sum(nAct), halfbandwidth, halfbandwidth, simData%nIsom, &
			L, 3 * halfbandwidth + 1, iPiv, Z, sum(nAct), info)
		if (info /= 0) then
			write (*,*) 'ERROR: Active-state matrix is singular!'
			stop
		end if
		deallocate( iPiv )
		
		! Convert solution to pseudo-steady state populations
		do r = minval(nRes)+1, simData%nGrains
			do n = 1, simData%nIsom
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
	function reservoirCutoffs(simData, isomerList, rxnList)
	
		type(Simulation), intent(in)				:: 	simData
		type(Isomer), dimension(:), intent(in)		:: 	isomerList
		type(Reaction), dimension(:), intent(in)	:: 	rxnList
		integer, dimension(1:size(isomerList))		::	reservoirCutoffs
		
		real(8) Eres
		integer i, j, t, start
		
		! Determine reservoir cutoffs by looking at transition state energies
		do i = 1, simData%nIsom
			if (isomerList(i)%numSpecies == 1) then
			
				start = ceiling((isomerList(i)%E0 - simData%Emin) / simData%dE) + 1
				Eres = simData%Emax
				
				do t = 1, simData%nRxn
					if (rxnList(t)%isomerList(1) == i .or. rxnList(t)%isomerList(2) == i) then
						if (rxnList(t)%E0 < Eres) Eres = rxnList(t)%E0
					end if
				end do
				
				reservoirCutoffs(i) = floor((Eres - 10 * simData%alpha - simData%Emin) / simData%dE)
				
				! Make sure reservoirCutoffs(i) is valid (i.e. points to grain at or above ground state)
				if (reservoirCutoffs(i) < start) then
					reservoirCutoffs(i) = start
				end if
				
			else
				reservoirCutoffs(i) = simData%nGrains
			end if
		end do

	end function
	
	

end module
