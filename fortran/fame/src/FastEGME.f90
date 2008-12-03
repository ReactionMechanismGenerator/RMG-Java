! ==============================================================================
!
! 	FastEGME.f90
!
! 	Written by Josh Allen (jwallen@mit.edu)
!
! ==============================================================================

module FastEGMEModule

	use SimulationModule
	use SpeciesModule
	use ReactionModule

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
	!   rxnData - The chemical data about each transition state.
	!	nRes - The reservoir cutoff grains for each unimolecular isomer.
	
	subroutine getReservoirCutoffs(simData, rxnData, nRes)
	
		type(Simulation), intent(in)				:: 	simData
		type(Reaction), dimension(:), intent(in)	:: 	rxnData
		integer, dimension(:), intent(out)			::	nRes
		
		integer i, t
		
		! Determine reservoir cutoffs by looking at transition state energies
		do i = 1, simData%nUni
			Eres = simData%Emax
			do t = 1, simData%nRxn
				if (rxnData(t)%isomer(1) == i .or. rxnData(t)%isomer(2) == i) then
					if (rxnData(t)%E < Eres) Eres = rxnData(t)%E
				end if
			end do
			nRes(i) = floor((Eres - 10 * simData%alpha) / simData%dE)
		end do

	end subroutine
	
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
	subroutine ssrsRates(simData, nRes, Mi, Hn, Kij, Fim, Gnj, Jnm, bi, bn, K)

		! Provide parameter type checking of inputs and outputs
		type(Simulation), intent(in)				:: 	simData
		integer, dimension(:), intent(in) 			:: 	nRes
		real(8), dimension(:,:,:), intent(in)		:: 	Mi
		real(8), dimension(:,:,:), intent(in)		:: 	Hn
		real(8), dimension(:,:,:), intent(in)		:: 	Kij
		real(8), dimension(:,:,:), intent(in)		:: 	Fim
		real(8), dimension(:,:,:), intent(in)		:: 	Gnj
		real(8), dimension(:,:,:), intent(in)		:: 	Jnm
		real(8), dimension(:,:), intent(inout)		:: 	bi
		real(8), dimension(:,:), intent(in)			:: 	bn
		real(8), dimension(:,:), intent(out) 		:: 	K
		
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
				if (i /= j .and. Kij( nGrains, i, j ) > 0) then
					s = min(nAct(i), nAct(j));
            		do r = 0, s-1
						L( sum(nAct(1:i)) - r, sum(nAct(1:j)) - r ) = &
                    		Kij( nGrains - r, i, j )
					end do
				end if
			end do
		end do
		
		! Invert L (requires LAPACK)
		allocate( iPiv(1:sum(nAct)), work(1:sum(nAct)) )
		call DGETRF( sum(nAct), sum(nAct), L, sum(nAct), iPiv, info )
		if (info > 0) then
			write (*,*), "Active-state grain matrix is singular!"
			return
		end if
		call DGETRI( sum(nAct), L, sum(nAct), iPiv, work, sum(nAct), info )
		if (info > 0) then
			write (*,*), "Active-state grain matrix is singular!"
			return
		end if
		deallocate( iPiv, work )

		! Renormalize unimolecular distributions such that the reservoir grains sum to unity
		do i = 1, nUni
    		bi(:,i) = bi(:,i) / sum(bi(1:nRes(i),i));
		end do

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
			call DGEMV('N', nRes(i), nRes(i), &
				one, Mi(1:nRes(i), 1:nRes(i), i), nRes(i), &
				bi(1:nRes(i), i), 1, &
				zero, tempV1(1:nRes(i)), 1)
			K(i,i) = K(i,i) + sum( tempV1(1:nRes(i)) )	
			
			! Unimolecular reaction/collision terms
			do j = 1, nUni
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
			end do
		
			! Bimolecular reaction terms
			do m = 1, nMulti
				do j = 1, nUni
					if (Fim(nGrains, j, m) > 0) then
				
						! F_jma * b_m
						tempV3(1:nAct(j)) =  Fim(nRes(j)+1:nGrains, j, m) * bn(nRes(j)+1:nGrains, m)
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
			call DGEMV('N', nGrains, nGrains, &
				one, Hn(:, :, n), nGrains, &
				bn(:, n), 1, &
				zero, tempV1(1:nGrains), 1)
			K(nUni+n,nUni+n) = K(nUni+n,nUni+n) + sum( tempV1(1:nGrains) )	
			
			! Bimolecular-bimolecular interconversion terms - commented because it was assumed that there weren't any of these
			! || J_nm * b_m ||
! 			do m = 1, nMulti
! 				if (m /= n .and. Jnm(nGrains, n, m) > 0) then
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
					if (Gnj(nGrains, n, i) > 0) then
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
						if (Fim(nGrains, j, m) > 0 .and. Gnj(nGrains, n, i) > 0) then
							! F_jma * b_m
							tempV3(1:nAct(j)) =  Fim(nRes(j)+1:nGrains, j, m) * bn(nRes(j)+1:nGrains, m)
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


end module
