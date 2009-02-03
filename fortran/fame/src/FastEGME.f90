! ==============================================================================
!
! 	FastEGME.f90
!
! 	Written by Josh Allen (jwallen@mit.edu)
!
! ==============================================================================

module FastEGMEModule

	use SimulationModule
	use IsomerModule
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
	!   uniData - The chemical data about each unimolecular isomer.
	!   rxnData - The chemical data about each transition state.
	!
	! Returns:
	!	nRes - The reservoir cutoff grains for each unimolecular isomer.
	subroutine getReservoirCutoffs(simData, uniData, rxnData, nRes)
	
		type(Simulation), intent(in)				:: 	simData
		type(Isomer), dimension(:), intent(in)		:: 	uniData
		type(Reaction), dimension(:), intent(in)	:: 	rxnData
		integer, dimension(:), intent(out)			::	nRes
		
		integer i, t, start
		
		! Determine reservoir cutoffs by looking at transition state energies
		do i = 1, simData%nUni
			start = ceiling((uniData(i)%E(1) - simData%Emin) / simData%dE) + 1
			Eres = simData%Emax
			do t = 1, simData%nRxn
				if (rxnData(t)%isomer(1) == i .or. rxnData(t)%isomer(2) == i) then
					if (rxnData(t)%E < Eres) Eres = rxnData(t)%E
				end if
			end do
			nRes(i) = floor((Eres - simData%Emin - 10 * simData%alpha) / simData%dE)
			! Make sure nRes(i) is valid (i.e. points to grain at or above ground state)
			if (nRes(i) < start) then
				nRes(i) = start
			end if
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
	subroutine ssrsRates(simData, uniData, rxnData, nRes, Mi, Hn, Kij, Fim, Gnj, Jnm, bi, bn, K)

		! Provide parameter type checking of inputs and outputs
		type(Simulation), intent(in)				:: 	simData
		type(Isomer), dimension(:), intent(in)		:: 	uniData
		type(Reaction), dimension(:), intent(in)	:: 	rxnData
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
! 		call fullInverse(L, nAct)

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
					if (Fim(nGrains,j,m) /= 0) then
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
						if (Fim(nGrains,j,m) /= 0 .and. Gnj(nGrains,n,i) /= 0) then
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

! 		! Sign check on K matrix (not sure why this is necessary...)
!         do i = 1, nUni+nMulti
!             do j = 1, nUni+nMulti
!                 if (i == j) then
!                     if (K(i,j) > 0) K(i,j) = -K(i,j)
!                 else
!                     if (K(i,j) < 0) K(i,j) = abs(K(i,j))
!                 end if
!             end do  
!         end do
        
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
			return
		end if
		call DGETRI( sum(nAct), L, sum(nAct), iPiv, work, sum(nAct), info )
		if (info > 0) then
			write (*,*), "Active-state grain matrix is singular!"
			return
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
		
		real(8), dimension(:,:), allocatable		:: 	ai
		real(8), dimension(:), allocatable			:: 	eqRatio
		integer i, j, nUni, nGrains, offI, offJ, r, s
		! Variables for BLAS and LAPACK
		integer, dimension(:), allocatable			::	iPiv
		integer										::	info
		real(8), dimension(:), allocatable			::	work
		integer										::	block
		
		nUni = simData%nUni
		nGrains = simData%nGrains
		
		! Determine equilibrium ratio and apply to equilibrium distributions
		allocate( ai(1:nGrains,1:nUni) )
		allocate( eqRatio(1:nUni) )
		call calcUniEqRatios(simData, uniData, rxnData, bi, eqRatio)
		do i = 1, nUni
			ai(:,i) = (eqRatio(i) * bi(:,i))**0.5
		end do
		
		! Symmetrize upper triangular half of matrix L
		do i = 1, nUni
			offI = sum(nAct(1:i))
			do j = 1, nUni
				offJ = sum(nAct(1:j))
				if (i == j) then
					! Blocks on diagonal are full blocks, so symmetrize all entries above diagonal
					do r = 0, nAct(i)-1
						do s = 0, r-1
							L(offI-r,offJ-s) = L(offI-r,offJ-s) * ai(nGrains-s,j) / ai(nGrains-r,i)
						end do
					end do
				elseif (i < j) then
					! Blocks off diagonal are themselves diagonal, so take advantage of this to speed symmetrization
					do r = 0, min(nAct(i),nAct(j))-1
						L(offI-r,offJ-r) = L(offI-r,offJ-r) * ai(nGrains-r,j) / ai(nGrains-r,i)
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
		
		! Unsymmetrize and restore full matrix L
		do i = 1, nUni
			offI = sum(nAct(1:i))
			do j = 1, nUni
				offJ = sum(nAct(1:j))
				if (i == j) then
					do r = 0, nAct(i)-1
						do s = 0, r-1
							L(offI-r,offJ-s) = L(offI-r,offJ-s) * ai(nGrains-r,i) / ai(nGrains-s,j)
						end do
						do s = r, nAct(j)-1
							L(offI-r,offJ-s) = L(offJ-s,offI-r) * ai(nGrains-r,i) / ai(nGrains-s,j)
						end do
					end do
				elseif (i < j) then
					do r = 0, nAct(i)-1
						do s = 0, nAct(j)-1
							L(offI-r,offJ-s) = L(offI-r,offJ-s) * ai(nGrains-r,i) / ai(nGrains-s,j)
							L(offJ-s,offI-r) = L(offI-r,offJ-s) * ai(nGrains-s,j) / ai(nGrains-r,i)
						end do
					end do
				end if
			end do
		end do
		
		! Clean up
		deallocate( ai, eqRatio )
		
		
	end subroutine
	
	! --------------------------------------------------------------------------
	
	! Subroutine: calcUniEqRatios()
	! 
	! Calculates the ratios of the concentrations of each unimolecular well
	! at equilibrium.
	!
	! Parameters:
	!
	subroutine calcUniEqRatios(simData, uniData, rxnData, bi, C)
	
		type(Simulation), intent(in)				:: 	simData
		type(Isomer), dimension(:), intent(in)		:: 	uniData
		type(Reaction), dimension(:), intent(in)	:: 	rxnData
		real(8), dimension(:,:), intent(in)			:: 	bi
		real(8), dimension(:), intent(out)			:: 	C
		
		real(8), dimension(:), allocatable			:: 	nAct
		real(8), dimension(:,:), allocatable		:: 	ai
		real(8), dimension(:,:), allocatable		:: 	A
		real(8), dimension(:), allocatable			:: 	b
		integer i, j, t, ind, nUni, nGrains
		real(8) Hrxn, Grxn, Keq_298, Keq
		integer, dimension(:), allocatable			:: iPiv
		integer info
		
		nUni = simData%nUni
		nGrains = simData%nGrains
		
		! Nothing to do if only one unimolecular well
		if (nUni == 1) then
			C(1) = 1
			return
		end if
		
		allocate( A(1:nUni,1:nUni) )
		allocate( b(1:nUni) )
		A = 0.0 * A
		b = 0.0 * b
		
		! Determine rate coefficients for each reaction
		! This is done by iterating over transition states
		ind = 1
		do t = 1, simData%nRxn
			if (rxnData(t)%isomer(1) <= simData%nUni .and. rxnData(t)%isomer(2) <= simData%nUni .and. ind < nUni) then			! bi <---> bi
				
				! Determine the wells involved
				i = rxnData(t)%isomer(1)
				j = rxnData(t)%isomer(2)
				
				! Determine equilibrium constant at rxn conditions from thermochemical data
				Hrxn = sum(uniData(i)%H) - sum(uniData(j)%H)
				Grxn = sum(uniData(i)%G) - sum(uniData(j)%G)
				Keq_298 = exp(-Grxn * 1000. / 8.314472 / 298.) 
				Keq = Keq_298 * exp(-Hrxn * 1000. / 8.314472 * (1./simData%T - 1./298.))
								
				A(ind,i) = -1
				A(ind,j) = Keq
				do k = 1, simData%nUni
					if (k /= i .and. k /= j) A(ind,k) = 0
				end do
				b(ind) = 0
				ind = ind + 1
				
			end if
		end do
		
        if (ind /= nUni) then
			write (*,*), 'ERROR: Incorrect number of isomerization reactions!'
			stop
		end if
		
		do i = 1, nUni-1
			A(nUni,i) = 0
		end do
		A(nUni,nUni) = 1
		b(nUni) = 1

		allocate( iPiv(1:nUni) )
		call DGESV( nUni, 1, A, nUni, iPiv, b, nUni, info )
		if (info > 0) then
			write (*,*), "Equilibrium ratio matrix is singular!"
			stop
		end if
		
		C = b / sum(b)
		
		deallocate( A, b, iPiv )
		
	end subroutine

end module
