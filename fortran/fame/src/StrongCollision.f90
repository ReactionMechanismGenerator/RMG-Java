! ==============================================================================
!
! 	StrongCollision.f90
!
! 	Written by Josh Allen (jwallen@mit.edu)
!
! ==============================================================================

module StrongCollisionModule

	use SimulationModule
	use SpeciesModule
	use ReactionModule

contains

	! --------------------------------------------------------------------------
	
	! Subroutine: ssmscRates()
	! 
	! Uses the steady state/modified strong collision approach of Chang,
	! Bozzelli, and Dean to estimate the phenomenological rate coefficients 
	! from a full ME matrix for the case of a chemically activated system.
	!
	! Parameters:
	!
	subroutine ssmscRates(simData, uniData, multiData, rxnData, Kij, Fim, Gnj, Jnm, bi, bn, K)

		! Provide parameter type checking of inputs and outputs
		type(Simulation), intent(in)				:: 	simData
		type(UniWell), dimension(:), intent(in)		:: 	uniData
		type(MultiWell), dimension(:), intent(in)	:: 	multiData
		type(Reaction), dimension(:), intent(in)	:: 	rxnData
		real(8), dimension(:,:,:), intent(in)		:: 	Kij
		real(8), dimension(:,:,:), intent(in)		:: 	Fim
		real(8), dimension(:,:,:), intent(in)		:: 	Gnj
		real(8), dimension(:,:,:), intent(in)		:: 	Jnm
		real(8), dimension(:,:), intent(inout)		:: 	bi
		real(8), dimension(:,:), intent(in)			:: 	bn
		real(8), dimension(:,:), intent(out) 		:: 	K
		
		! Steady-state populations
		real(8), dimension(:,:), allocatable	:: 	p
		
		! Steady-state matrix and vector
		real(8), dimension(:,:), allocatable	:: 	A
		real(8), dimension(:), allocatable		:: 	b
		! Collision frequency and efficiency
		real(8), dimension(:), allocatable				::	w
		real(8) eps, mu
		! Gas concentration
		real(8) gasConc
		! Number of active-state energy grains for each unimolecular isomer
		integer, dimension(:), allocatable				:: 	nAct
		! Indices i and j represent sums over unimolecular wells
		integer											::	i, j
		! Indices m and n represent sums over bimolecular sources/sinks
		integer											::	m, n
		! Indices r and s represent sums over energy grains
		integer											::	r, s
		! Variables for BLAS and LAPACK
		integer, dimension(:), allocatable				::	iPiv
		integer											::	info
		integer	src, start
		real(8)	val
		
		! Fraction of collisions that result in deactivation
		eps = 0.05

		! Gas concentration in molecules/m^3
		gasConc = simData%P * 1e5 / 1.381e-23 / simData%T

		! Determine collision frequency for each isomer
		allocate( w(1:simData%nUni) )
		do i = 1, simData%nUni
			mu = 1.0 / ( 1.0 / uniData(i)%MW + 1.0 / simData%bathGas%MW ) / 6.022e26
			call collisionFrequency(simData%T, 0.5 * (uniData(i)%sigma + simData%bathGas%sigma), &
				0.5 * (uniData(i)%eps + simData%bathGas%eps), mu, gasConc, w(i))
		end do
		
		! Zero rate coefficient matrix
		do i = 1, simData%nUni + simData%nMulti
			do j = 1, simData%nUni + simData%nMulti
				K(i,j) = 0
			end do
		end do
		
		! Find steady-state populations at each grain
		allocate( A(1:simData%nUni, 1:simData%nUni) )
		allocate( b(1:simData%nUni) )
		allocate( iPiv(1:simData%nUni) )
		allocate( p(1:simData%nGrains, 1:simData%nUni) )
			
		do src = 1, simData%nUni+simData%nMulti

			! Determine starting grain
			call getActiveSpaceStart(simData, uniData, rxnData, src, start)
			
			do r = 1, start-1
				do i = 1, simData%nUni
					p(r,i) = 0
				end do
			end do
			
			do r = start, simData%nGrains

				! Collisional deactivation
				do i = 1, simData%nUni
					A(i,i) = - eps * w(i);
				end do
				
				! Isomerization
				do i = 1, simData%nUni
					do j = 1, simData%nUni
						if (i /= j .and. Kij(r,i,j) > 1.0e-6) then
							A(i,j) = Kij(r,i,j)
							A(i,i) = A(i,i) - Kij(r,j,i)
						end if
					end do
				end do
				
				! Dissociation
				do i = 1, simData%nUni
					do n = 1, simData%nMulti
						A(i,i) = A(i,i) - Gnj(r,n,i)
						b(i) = 0
					end do
				end do
				
				! Activation
				if (src > simData%nUni) then
					n = src - simData%nUni
					do i = 1, simData%nUni
						b(i) = Fim(r,i,n) * bn(r,n)
					end do
				else
					b(src) = eps * w(src) * bi(r,src)
				end if
				
! 				do i = 1, simData%nUni
! 					write (*,*), A(i,:)
! 				end do
				
				! Solve for steady-state population
				call DGESV( simData%nUni, 1, A, simData%nUni, iPiv, b, simData%nUni, info )
				if (info > 0) then
					write (*,*), "A singular matrix was encountered! Aborting."
					stop
				end if
				p(r,:) = -b
				
			end do

			! Calculate rates
			
			! Stabilization rates (i.e.) R + R' --> Ai or M --> Ai
			do i = 1, simData%nUni
				val = sum(eps * w(i) * p(:,i))
				val = abs(val)
				K(i,src) = K(i,src) + val
				K(src,src) = K(src,src) - val
			end do
			
			! Dissociation rates (i.e.) R + R' --> Bn + Cn or M --> Bn + Cn
			if (src <= simData%nUni) then
				! Thermal activation
				do i = 1, simData%nUni
					do n = 1, simData%nMulti
						if (Gnj(simData%nGrains,n,i) > 0) then
							val = sum(Gnj(:,n,i) * p(:,i))
							val = abs(val)
							K(n+simData%nUni,src) = K(n+simData%nUni,src) + val
							K(src,src) = K(src,src) - val
						end if
					end do
				end do
			else
				! Chemical activation
				do i = 1, simData%nUni
					do n = 1, simData%nMulti
						if (n /= src - simData%nUni .and. Gnj(simData%nGrains,n,i) > 0) then
							val = sum(Gnj(:,n,i) * p(:,i))
							val = abs(val)
							K(n+simData%nUni,src) = K(n+simData%nUni,src) + val
							K(src,src) = K(src,src) - val
						end if
					end do
				end do
			end if
			
		end do

		! DEBUG: Sign check
! 		do i = 1, simData%nUni + simData%nMulti
! 			do j = 1, simData%nUni + simData%nMulti
! 				if (i /= j) then
! 					if (K(i,j) < 0) write (*,*), K(i,j), 'WARNING: Sign at', i, ',', j, 'is negative.'
! 				end if
! 			end do
! 		end do

		! Clean up
		deallocate( w, A, b, iPiv, p )

	end subroutine

	! --------------------------------------------------------------------------

	! collisionFrequency() 
	!
	!   Computes the Lennard-Jones (12-6) collision frequency.
	!
	! Input:
	!   T = absolute temperature in K
	!   sigma = effective Lennard-Jones well-minimum parameter for collision in m
	!   eps = effective Lennard-Jones well-depth parameter for collision in J
	!   mu = reduced mass of molecules involved in collision in kg
	!   Na = molecules of A per m^3
	!
	! Output:
	!   omega = collision frequency in molecules m^-3 s^-1
	!
	subroutine collisionFrequency(T, sigma, eps, mu, Na, omega)
	
		! Provide parameter type-checking
		real(8), intent(in)					:: 	T
		real(8), intent(in)					:: 	sigma
		real(8), intent(in)					:: 	eps
		real(8), intent(in)					:: 	mu
		real(8), intent(in)					:: 	Na
		real(8), intent(out)				:: 	omega
		
		real(8)		:: kB = 1.3806504e-23
		real(8)		:: collisionIntegral
		
		collisionIntegral = 1.16145 / T**0.14874 + 0.52487 / exp(0.77320 * T) + 2.16178 / exp(2.43787 * T) &
			-6.435/10000 * T**0.14874 * sin(18.0323 * T**(-0.76830) - 7.27371)
    
		! Evaluate collision frequency
		omega = collisionIntegral *	sqrt(8 * kB * T / 3.141592654 / mu) * 3.141592654 * sigma**2 * Na
	
	end subroutine
	
	! --------------------------------------------------------------------------

	! Subroutine: getActiveSpaceStart()
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
	subroutine getActiveSpaceStart(simData, uniData, rxnData, well, start)
	
		type(Simulation), intent(in)				:: 	simData
		type(UniWell), dimension(:), intent(in)		:: 	uniData
		type(Reaction), dimension(:), intent(in)	:: 	rxnData
		integer, intent(in)							::	well
		integer, intent(out)						::	start
		
		integer i, t, start0
		
		start = simData%nGrains
			
! 		! Determine reservoir cutoffs by looking at transition state energies
! 		do i = 1, simData%nUni
! 			start = ceiling(uniData(i)%E / simData%dE) + 1
! 			Eres = simData%Emax
! 			do t = 1, simData%nRxn
! 				if (rxnData(t)%isomer(1) == i .or. rxnData(t)%isomer(2) == i) then
! 					if (rxnData(t)%E < Eres) Eres = rxnData(t)%E
! 				end if
! 			end do
! 			start0 = floor(Eres / simData%dE)
! 			write (*,*), start0
! 			if (start0 < start) then
! 				start = start0
! 			end if
! 		end do
		
		i = well
		Eres = simData%Emax
		do t = 1, simData%nRxn
			if (rxnData(t)%isomer(1) == i .or. rxnData(t)%isomer(2) == i) then
				if (rxnData(t)%E < Eres) Eres = rxnData(t)%E
			end if
		end do
		start = floor(Eres / simData%dE)

	end subroutine
	
	! --------------------------------------------------------------------------

end module
