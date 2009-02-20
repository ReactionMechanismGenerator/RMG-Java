! ==============================================================================
!
! 	FastEGME.f90
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

	! Subroutine: getReservoirCutoffs()
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
	function reservoirCutoffs(simData, uniData, rxnData)
	
		type(Simulation), intent(in)				:: 	simData
		type(Isomer), dimension(:), intent(in)		:: 	uniData
		type(Reaction), dimension(:), intent(in)	:: 	rxnData
		integer, dimension(1:size(uniData))			::	reservoirCutoffs
		
		real(8) Eres
		integer i, j, t, start
		
		! Determine reservoir cutoffs by looking at transition state energies
		do i = 1, simData%nUni
			start = ceiling((uniData(i)%E(1) - simData%Emin) / simData%dE) + 1
			Eres = simData%Emax
			do t = 1, simData%nRxn
				if (rxnData(t)%isomer(1) == i .or. rxnData(t)%isomer(2) == i) then
					if (rxnData(t)%E < Eres) Eres = rxnData(t)%E
				end if
			end do
			reservoirCutoffs(i) = floor((Eres - 10 * simData%alpha - simData%Emin) / simData%dE)
			! Make sure reservoirCutoffs(i) is valid (i.e. points to grain at or above ground state)
			if (reservoirCutoffs(i) < start) then
				reservoirCutoffs(i) = start
			end if
		end do
		
		do i = 1, size(reservoirCutoffs)
			j = ceiling((uniData(i)%E(1) - simData%Emin) / simData%dE) + 1
			if (reservoirCutoffs(i) < j .or. reservoirCutoffs(i) >= simData%nGrains) then
				write (*,*), 'ERROR: Invalid reservoir grain.'
				stop
			end if
		end do

	end function
	
	! --------------------------------------------------------------------------
	
	! Subroutine: ssrsRates()
	! 
	! Uses the steady state/reservoir state method of Green and Bhatti to
	! estimate the phenomenological rate coefficients from a full ME matrix.
	!
	! Parameters:
	!	nRes - The reservoir cutoff grains for each unimolecular isomer.
	! 	Mi - Collisional transition matrix for each unimolecular well
	! 	Hn - Collisional transition matrix for each bimolecular well
	! 	Kij - Reactive transition from unimolecular well to unimolecular well
	! 	Fim - Reactive transition from bimolecular well to unimolecular well
	! 	Gnj - Reactive transition from unimolecular well to bimolecular well
	! 	Jnm - Reactive transition from bimolecular well to bimolecular well
	!	bi - Equilibrium distributions for unimolecular wells
	! 	bn - Equilibrium distributions for bimolecular wells
	!	K - The calculated phenomenologicate rate coefficient matrix in s^-1.
	subroutine ssrsRates(simData, uniData, multiData, rxnData, nRes, Mi, Hn, Kij, Fim, Gnj, Jnm, bi, K)

		! Provide parameter type checking of inputs and outputs
		type(Simulation), intent(in)				:: 	simData
		type(Isomer), dimension(:), intent(in)		:: 	uniData
		type(Isomer), dimension(:), intent(in)      ::  multiData
        type(Reaction), dimension(:), intent(in)	:: 	rxnData
		integer, dimension(:), intent(in) 			:: 	nRes
		real(8), dimension(:,:,:), intent(in)		:: 	Mi
		real(8), dimension(:,:,:), intent(in)		:: 	Hn
		real(8), dimension(:,:,:), intent(in)		:: 	Kij
		real(8), dimension(:,:,:), intent(in)		:: 	Fim
		real(8), dimension(:,:,:), intent(in)		:: 	Gnj
		real(8), dimension(:,:,:), intent(in)		:: 	Jnm
		real(8), dimension(:,:), intent(in)			:: 	bi
		real(8), dimension(:,:), intent(out) 		:: 	K
		
		real(8), dimension(:,:), allocatable		:: 	ai
		real(8), dimension(:,:), allocatable		:: 	an
		
		! Active state grain matrix
		real(8), dimension(:,:), allocatable	:: 	L
		! Number of active-state energy grains for each unimolecular isomer
		integer, dimension(:), allocatable				:: 	nAct
		! Indices i and j represent sums over unimolecular wells
		integer											::	i, j
		! Indices m and n represent sums over bimolecular sources/sinks
		integer											::	m, n
		! Indices r and s represent sums over energy grains
		integer											::	r, s
		! Dummy vectors
		real(8), dimension(:), allocatable				:: 	tempV1, tempV2, tempV3	
		! Variables for BLAS and LAPACK
		real(8)											:: 	one, zero
		integer, dimension(:), allocatable				::	iPiv
		integer											::	info
		real(8), dimension(:), allocatable				::	work
		! Variables for shorthand
		integer	nUni, nMulti, nGrains
		integer found
        
		nUni = simData%nUni
		nMulti = simData%nMulti
		nGrains = simData%nGrains
		
		! Load nAct
		allocate( nAct(1:nUni) )
		do i = 1, nUni
			nAct(i) = nGrains - nRes(i)
		end do

		! Load L
		allocate( L(1:sum(nAct), 1:sum(nAct)) )
		do i = 1, nUni
			! Active-state collisional interconversion + reactive loss to all other wells
			L( sum(nAct(1:i-1))+1:sum(nAct(1:i)), sum(nAct(1:i-1))+1:sum(nAct(1:i)) ) = &
				Mi( nRes(i)+1:nGrains, nRes(i)+1:nGrains, i )
			! Reactive gain from all other unimolecular wells
			do j = 1, nUni
				if (i /= j) then
					s = min(nAct(i), nAct(j));
            		do r = 0, s-1
						L( sum(nAct(1:i)) - r, sum(nAct(1:j)) - r ) = &
                    		Kij( nGrains - r, i, j )
					end do
				end if
			end do
		end do
		
		! Invert L
		call symmetricInverse(L, nAct, simData, uniData, rxnData, bi)
!		call fullInverse(L, nAct)

		! Initialize phenomenological rate coefficient matrix
		do i = 1, nUni + nMulti
			do j = 1, nUni + nMulti
				K(i,j) = 0
			end do
		end do
		
		! Initialize temporary vectors used in determining phenomenological rate coefficient matrix
		allocate( tempV1(1:nGrains), tempV2(1:nGrains), tempV3(1:nGrains) )
		
		! Set variables used for BLAS
		one = 1.0
		zero = 0.0

		! Unimolecular rows of phenomenological rate coefficient matrix
		! -------------------------------------------------------------
		do i = 1, nUni
		
			! Reservoir rearrangement terms (on diagonal)
			! || M_irr * b_ir ||
 			!call DGEMV('N', nRes(i), nRes(i), &
 			!	one, Mi(1:nRes(i), 1:nRes(i), i), nRes(i), &
 			!	bi(1:nRes(i), i), 1, &
 			!	zero, tempV1(1:nRes(i)), 1)
 			!K(i,i) = K(i,i) + sum( tempV1(1:nRes(i)) )	
			
			! Unimolecular reaction/collision terms
			do j = 1, nUni
            	if (i /= j) then
					! M_jar * b_jr
					call DGEMV('N', nAct(j), nRes(j), &
						one, Mi(nRes(j)+1:nGrains, 1:nRes(j), j), nAct(j), &
						bi(1:nRes(j), j), 1, &
						zero, tempV3(1:nAct(j)), 1)
					! L_ij * M_jar * b_jr
					call DGEMV('N', nAct(i), nAct(j), &
						one,  L( sum(nAct(1:i-1))+1:sum(nAct(1:i)), sum(nAct(1:j-1))+1:sum(nAct(1:j)) ), nAct(i), &
						tempV3(1:nAct(j)), 1, &
						zero, tempV2(1:nAct(i)), 1)
					! M_ira * L_ij * M_jar * b_jr
					call DGEMV('N', nRes(i), nAct(i), &
						one,  Mi(1:nRes(i), nRes(i)+1:nGrains, i), nRes(i), &
						tempV2(1:nAct(i)), 1, &
						zero, tempV1(1:nRes(i)), 1)
					! || M_ira * L_ij * M_jar * b_jr || 
					K(i,j) = K(i,j) - sum( tempV1(1:nRes(i)) )
				end if
            end do
		
			! Bimolecular reaction terms
			do m = 1, nMulti
				do j = 1, nUni
					if (Fim(nGrains,j,m) /= 0) then
						! F_jma * b_m
						tempV3(1:nAct(j)) =  Fim(nRes(j)+1:nGrains, j, m)
						! L_ij * F_jma * b_m
						call DGEMV('N', nAct(i), nAct(j), &
							one,  L( sum(nAct(1:i-1))+1:sum(nAct(1:i)), sum(nAct(1:j-1))+1:sum(nAct(1:j)) ), nAct(i), &
							tempV3(1:nAct(j)), 1, &
							zero, tempV2(1:nAct(i)), 1)
						! M_ira * L_ij * F_jma * b_m
						call DGEMV('N', nRes(i), nAct(i), &
							one,  Mi(1:nRes(i), nRes(i)+1:nGrains, i), nRes(i), &
							tempV2(1:nAct(i)), 1, &
							zero, tempV1(1:nRes(i)), 1)
						! || M_ira * L_ij * F_jma * b_m || 
						K(i,nUni+m) = K(i,nUni+m) - sum( tempV1(1:nRes(i)) )
					end if
				end do
			end do

		end do
		
        ! Bimolecular rows of phenomenological rate coefficient matrix
		! ------------------------------------------------------------
		do n = 1, nMulti
		
			! Loss term (on diagonal)
			! || H_n * b_n ||
 			!call DGEMV('N', nGrains, nGrains, &
 			!	one, Hn(:, :, n), nGrains, &
 			!	bn(:, n), 1, &
 			!	zero, tempV1(1:nGrains), 1)
 			!K(nUni+n,nUni+n) = K(nUni+n,nUni+n) + sum( tempV1(1:nGrains) )	
			
			! Bimolecular-bimolecular interconversion terms - commented because it was assumed that there weren't any of these
			! || J_nm * b_m ||
! 			do m = 1, nMulti
! 				if (m /= n) then
! 					call DGEMV('N', nGrains, nGrains, &
! 						one, Jnm(:, n, m), nGrains, &
! 						bn(:, m), 1, &
! 						zero, tempV1(1:nGrains), 1)
! 					K(nUni+n, nUni+m) = K(nUni+n, nUni+m) + sum( tempV1(1:nGrains) )
! 				end if
! 			end do
			
            ! Unimolecular reaction terms
			do j = 1, nUni
				do i = 1, nUni
					if (Gnj(nGrains,n,i) /= 0) then
						! M_jar * b_jr
						call DGEMV('N', nAct(j), nRes(j), &
							one, Mi(nRes(j)+1:nGrains, 1:nRes(j), j), nAct(j), &
							bi(1:nRes(j), j), 1, &
							zero, tempV3(1:nAct(j)), 1)
						! L_ij * M_jar * b_jr
						call DGEMV('N', nAct(i), nAct(j), &
							one,  L( sum(nAct(1:i-1))+1:sum(nAct(1:i)), sum(nAct(1:j-1))+1:sum(nAct(1:j)) ), nAct(i), &
							tempV3(1:nAct(j)), 1, &
							zero, tempV2(1:nAct(i)), 1)
						! G_nia * L_ij * M_jar * b_jr
						tempV1(1:nAct(i)) = Gnj( nRes(i)+1:nGrains, n, i ) * tempV2(1:nAct(i))
						! || G_nia * L_ij * M_jar * b_jr || 
						K(nUni+n,j) = K(nUni+n,j) - sum( tempV1(1:nAct(i)) )
					end if
				end do
			end do
			
			! Bimolecular reaction terms
			do m = 1, nMulti
				do i = 1, nUni
					do j = 1, nUni
						if (Fim(nGrains,j,m) /= 0 .and. Gnj(nGrains,n,i) /= 0 .and. n /= m) then
							! F_jma * b_m
							tempV3(1:nAct(j)) =  Fim(nRes(j)+1:nGrains, j, m)
							! L_ij * F_jma * b_m
							call DGEMV('N', nAct(i), nAct(j), &
								one,  L( sum(nAct(1:i-1))+1:sum(nAct(1:i)), sum(nAct(1:j-1))+1:sum(nAct(1:j)) ), nAct(i), &
								tempV3(1:nAct(j)), 1, &
								zero, tempV2(1:nAct(i)), 1)
							! G_nia * L_ij * F_jma * b_m
							tempV1(1:nAct(i)) = Gnj( nRes(i)+1:nGrains, n, i ) * tempV2(1:nAct(i))
							! || G_nia * L_ij * F_jma * b_m || 
							K(nUni+n,nUni+m) = K(nUni+n,nUni+m) - sum( tempV1(1:nAct(i)) )
						end if
					end do
				end do
			end do
		
		end do
        
        ! Clean up
		deallocate( tempV1, tempV2, tempV3 )	
		deallocate( nAct, L )

	end subroutine

	! Subroutine: fullInverse()
	! 
	! Inverts the active-state matrix.
	!
	! Parameters:
	!
	subroutine fullInverse(L, nAct)
	
		real(8), dimension(:,:), intent(inout)		:: 	L
		integer, dimension(:), intent(in)			:: 	nAct
		! Variables for BLAS and LAPACK
		integer, dimension(:), allocatable			::	iPiv
		integer										::	info
		real(8), dimension(:), allocatable			::	work
		
		! Invert L (requires LAPACK)
		allocate( iPiv(1:sum(nAct)), work(1:sum(nAct)) )
		call DGETRF( sum(nAct), sum(nAct), L, sum(nAct), iPiv, info )
		if (info > 0) then
			write (*,*), "Active-state grain matrix is singular!"
			stop
		end if
		call DGETRI( sum(nAct), L, sum(nAct), iPiv, work, sum(nAct), info )
		if (info > 0) then
			write (*,*), "Active-state grain matrix is singular!"
			stop
		end if
		deallocate( iPiv, work )
		
	end subroutine

	! --------------------------------------------------------------------------
	
	! Subroutine: symmetricInverse()
	! 
	! Inverts the active-state matrix efficiently via a prior symmetrization step.
	!
	! Parameters:
	!
	subroutine symmetricInverse(L, nAct, simData, uniData, rxnData, bi)
	
		real(8), dimension(:,:), intent(inout)		:: 	L
		integer, dimension(:), intent(in)			:: 	nAct
		type(Simulation), intent(in)				:: 	simData
		type(Isomer), dimension(:), intent(in)		:: 	uniData
		type(Reaction), dimension(:), intent(in)	:: 	rxnData
		real(8), dimension(:,:), intent(in)			:: 	bi
		
		integer i, j, nUni, nGrains, offI, offJ, r, s
		! Variables for BLAS and LAPACK
		integer, dimension(:), allocatable			::	iPiv
		integer										::	info
		real(8), dimension(:), allocatable			::	work
		integer										::	block
		real(8)	Keq
		real(8) Rgas
		real(8) temp
		
		Rgas = 8.314472 / 1000.
		
		nUni = simData%nUni
		nGrains = simData%nGrains
		
		! Symmetrize upper triangular half of matrix L
		do i = 1, nUni
			offI = sum(nAct(1:i))
			do j = 1, nUni
				offJ = sum(nAct(1:j))
				if (i == j) then
					! Blocks on diagonal are full blocks, so symmetrize all entries above diagonal
					do r = 0, nAct(i)-1
						do s = 0, r-1
							if (L(offI-r,offJ-s) /= 0) then
								L(offI-r,offJ-s) = L(offI-r,offJ-s) * sqrt( &
									uniData(j)%densStates(nGrains-s) / uniData(i)%densStates(nGrains-r) * &
									uniData(i)%Q / uniData(j)%Q * &
									exp(-(simData%E(nGrains-s) - simData%E(nGrains-r)) / (Rgas * simData%T)))
							end if
						end do
					end do
				elseif (i < j) then
					! Blocks off diagonal are themselves diagonal, so take advantage of this to speed symmetrization
					do r = 0, min(nAct(i),nAct(j))-1
						if (L(offI-r,offJ-r) /= 0) then
							L(offI-r,offJ-r) = L(offI-r,offJ-r) * sqrt( &
								uniData(i)%densStates(nGrains-r) / uniData(j)%densStates(nGrains-r) * &
								uniData(j)%Q / uniData(i)%Q)
						end if
					end do
				end if
			end do
		end do
		
		! Invert L (requires LAPACK)
! 		block = ILAENV(1, "DSYTRF", 'U', sum(nAct), sum(nAct), sum(nAct), sum(nAct))
! 		write (*,*), block, sum(nAct), block*sum(nAct)
		block = 8
		allocate( iPiv(1:sum(nAct)) )
		allocate( work(1:block*sum(nAct)) )
		call DSYTRF( 'U', sum(nAct), L, sum(nAct), iPiv, work, block*sum(nAct), info )
		if (info > 0) then
			write (*,*), "Active-state grain matrix is singular!"
			return
		end if
		call DSYTRI( 'U', sum(nAct), L, sum(nAct), iPiv, work, info )
		if (info > 0) then
			write (*,*), "Active-state grain matrix is singular!"
			return
		end if
		deallocate( iPiv, work )
		
		! Fill in lower triangle of matrix L
		do i = 1, sum(nAct)
			L(i,1:i-1) = L(1:i-1,i)
		end do
		
		! Unsymmetrize
		do i = 1, nUni
			offI = sum(nAct(1:i))
			do j = 1, nUni
				offJ = sum(nAct(1:j))
				do r = 0, nAct(i)-1
					do s = 0, nAct(j)-1
						if (uniData(j)%densStates(nGrains-s) /= 0 .and. &
							abs(simData%E(nGrains-r) - simData%E(nGrains-s)) / (Rgas * simData%T) < 700) then
							L(offI-r,offJ-s) = L(offI-r,offJ-s) * sqrt( &
								uniData(i)%densStates(nGrains-r) / uniData(j)%densStates(nGrains-s) * &
								uniData(j)%Q / uniData(i)%Q * &
								exp(-(simData%E(nGrains-r) - simData%E(nGrains-s)) / (Rgas * simData%T)))
						end if
					end do
				end do
			end do
		end do		
		
	end subroutine
	
	! --------------------------------------------------------------------------

end module
