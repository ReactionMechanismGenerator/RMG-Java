! ==============================================================================
!
!	Isomer.f90
!
!	Written by Josh Allen (jwallen@mit.edu)
!
! ==============================================================================

module IsomerModule

	use SimulationModule
	use SpeciesModule

	implicit none

	! Struct: Isomer
	! 
	! Contains the information for a (unimolecular or multimolecular) isomer well.
	type Isomer	
		integer									::	numSpecies
		integer, dimension(:), allocatable		::	speciesList		! List of indices to species in the well
		real(8)									:: 	E0				! Ground-state electronic + zero-point energy in kJ/mol
		real(8)									:: 	dHf				! Standard enthalpy of formation in kJ/mol
		real(8)									:: 	dGf				! Standard Gibbs free energy of formation in kJ/mol
		real(8)									:: 	MW				! Molecular weight in g/mol
		real(8), dimension(:), allocatable		:: 	densStates		! Convolved density of states for all species in well
		real(8), dimension(:), allocatable		:: 	eqDist			! (Boltzmann) equilibrium distribution for well
		real(8)									::	Q				! Partition function at current temperature
		real(8)									::	omega			! Collisional frequency in s^-1
		real(8), dimension(:,:), allocatable	::	Mcoll			! Collisional transfer probability matrix
	end type

contains

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!
	! Function: eqDist()
	! 
	! Calculates the equilibrium distribution and partition function from
	! density-of-states data for an isomer.
	!
	! Parameters:
	!		isom - The isom to calculate the distribution for.
	!		E - The grain energies in kJ/mol.
	!		T - The absolute temperature in K.
	!		eqDist - The vector in which to place the equilibrium distribution.
	subroutine eqDist(isom, E, T)
		
		type(Isomer), intent(inout)				:: isom
		real(8), dimension(:), intent(in)		:: E
		real(8), intent(in)						:: T
		
		integer nGrains
		real(8)	R 					! Gas constant in J mol^-1 K^-1
		integer	s					! Dummy index
		real(8) dE
		
		nGrains = size(E)
		dE = E(2) - E(1)
		
		R = 8.314472 / 1000.
	
		if (.not. allocated(isom%eqDist)) then
			allocate( isom%eqDist(1:size(E)) )
		end if

		! Calculate unnormalized eqDist
		do s = 1, size(E)
			isom%eqDist(s) = isom%densStates(s) * exp(-E(s) / (R * T))
		end do
		
		! Normalize eqDist
		isom%Q = sum(isom%eqDist) * dE			
		isom%eqDist = isom%eqDist / sum(isom%eqDist)
		
	end subroutine
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!
	! Function: eqRatio()
	!
	!   Calculates the equilibrium ratio value for a single isomer.
	!
	! Parameters:
	!   dGf = standard free energy of formation in kJ/mol
	!   dHf = standard enthalpy of formation in kJ/mol
	!	T = absolute temperature in K
	function eqRatio(dGf, dHf, T)

		real(8), intent(in)			:: dGf
		real(8), intent(in)			:: dHf
		real(8), intent(in)			:: T
		real(8)						:: eqRatio
		
		real(8) R, T0
		
		R = 8.314472 / 1000
		T0 = 298

		eqRatio = exp(- dGf / (R * T0) - dHf / R * (1. / T - 1. / T0))

	end function
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!
	! Function: collision()
	!
	! Calculates the collisional transfer rate matrix for a given isomer.
	!
	subroutine collision(T, P, simData, isom, speciesList)

		real(8), intent(in)							::	T
		real(8), intent(in)							::	P
		type(Simulation), intent(in)				:: 	simData
		type(Isomer), intent(inout)					:: 	isom
		type(Species), dimension(:), intent(in)		:: 	speciesList
		
		real(8)			:: gasConc			! Gas concentration in molecules/m^3
		
		type(Species)	:: spec
		
		real(8), dimension(:,:), allocatable		::	prob
		real(8), dimension(:), allocatable		::	C		! Vector of normalization coefficients
		
		integer			:: bandwidth, halfbandwidth
		integer			:: lb, ub, start
		integer			:: r, s
		
		! Determine (ideal) gas concentration
		gasConc = P * 1e5 / 1.381e-23 / T
		
		! Allocate collision matrix
		if (.not. allocated(isom%Mcoll)) then
			allocate( isom%Mcoll(1:simData%nGrains, 1:simData%nGrains) )
		end if
		isom%Mcoll = 0 * isom%Mcoll
		
		! Allocate/zero vector of normalization coefficients
		allocate( C(1:simData%nGrains) )
		C = 0 * C
		
		! Get species of isomer (isom should be unimolecular!)
		spec = speciesList(isom%speciesList(1))
		
		! Determine bandwidth (at which transfer probabilities are so low that they can be truncated
		! with negligible error)
		halfbandwidth = getHalfBandwidth(simData%alpha, simData%dE)
		bandwidth = 2 * halfbandwidth + 1
		
		! Determine start grain (corresponding to isomer ground-state energy)
		start = ceiling((isom%E0 - simData%Emin) / simData%dE) + 1
		
		! Determine unnormalized entries in collisional tranfer probability matrix for the current isomer
		do r = start, simData%nGrains
			lb = max(r - halfbandwidth, start)
			ub = min(r + halfbandwidth, simData%nGrains)
			do s = lb, ub
				isom%Mcoll(s,r) = transferRate(s, r, simData%E(s), simData%E(r), simData%alpha, &
					isom%E0, isom%densStates, T)
			end do
		end do
		
		! Normalize using detailed balance
		do r = simData%nGrains, start, -1
			C(r) = (1 - sum(C(r+1:simData%nGrains) * isom%Mcoll(r+1:simData%nGrains,r))) / sum(isom%Mcoll(1:r,r))
		end do
		!do r = start, simData%nGrains
		!	C(r) = (1 - sum(C(start:r) * isom%Mcoll(start:r,r))) / sum(isom%Mcoll(r:simData%nGrains,r))
		!end do
		
		! Check for normalization consistency (i.e. all numbers are positive)
		!do r = start, simData%nGrains
		!	if (C(r) <= 0) then
		!		write (*,*) 'Error normalizing collisional transfer probabilities matrix!'
		!		stop
		!	end if
		!end do
			
		do r = start, simData%nGrains
			isom%Mcoll(r,1:r-1) = isom%Mcoll(r,1:r-1) * C(r)
			isom%Mcoll(1:r-1,r) = isom%Mcoll(1:r-1,r) * C(r)
			isom%Mcoll(r,r) = isom%Mcoll(r,r) * C(r) - 1
		end do
		
		! Multiply by collision frequency to determine collision rates
		isom%Mcoll = isom%omega * isom%Mcoll
		
		deallocate( C )
		
	end subroutine
	
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
	
	! --------------------------------------------------------------------------

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

		R = 8.314472 / 1000
		
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

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!
	! collisionFrequency() 
	!
	! Computes the Lennard-Jones (12-6) collision frequency.
	!
	! Parameters:
	!   T = absolute temperature in K
	!	P = absolute pressure in bar
	!	simData = simulation data
	!	spec = species of interest
	function collisionFrequency(T, P, simData, spec)
	
		! Provide parameter type-checking
		real(8), intent(in)					:: 	T
		real(8), intent(in)					:: 	P
		type(Simulation), intent(in)		:: 	simData
		type(Species), intent(in)			:: 	spec
		real(8)								::	collisionFrequency
		
		real(8)		::	kB = 1.3806504e-23
		real(8)		::	collisionIntegral
		real(8)		:: 	sigma
		real(8)		:: 	eps
		real(8)		:: 	mu
		real(8)		:: 	gasConc
		
		collisionIntegral = 1.16145 / T**0.14874 + 0.52487 / exp(0.77320 * T) + 2.16178 / exp(2.43787 * T) &
			-6.435/10000 * T**0.14874 * sin(18.0323 * T**(-0.76830) - 7.27371)
    
		gasConc = P * 1e5 / kB / T
		mu = 1 / (1/spec%MW + 1/simData%bathGas%MW) / 6.022e26
		sigma = 0.5 * (spec%sigma + simData%bathGas%sigma)
		eps = 0.5 * (spec%eps + simData%bathGas%eps)
		
		! Evaluate collision frequency
		collisionFrequency = collisionIntegral * &
			sqrt(8 * kB * T / 3.141592654 / mu) * 3.141592654 * sigma**2 * gasConc
	
	end function	
	
end module
