! ==============================================================================
!
!	DensityOfStates.f90
!
!	Written by Josh Allen (jwallen@mit.edu)
!
! ==============================================================================

module DensityOfStatesModule

	use SimulationModule
    use IsomerModule
	
	implicit none

contains

	! --------------------------------------------------------------------------
	!
	! Subroutine: densityOfStates()
	! 
	! Estimates the density of states for each isomer in the network.
	!
	! Parameters:
	!   simData - The simulation data.
	!   isomer - The chemical data for a single isomer.
	subroutine densityOfStates(simData, isom)
		
		! Provide parameter type checking of inputs and outputs
		type(Simulation), intent(in)	:: 	simData
		type(Isomer), intent(inout)		:: 	isom
		
		real(8) conv		! Conversion factor: cm^-1 per kJ/mol
		real(8) Emin		! Minimum energy of range in cm^-1
		real(8) Emax		! Maximum energy of range in cm^-1
		real(8) dE			! Grain size for density of states calculation in cm^-1
		real(8), dimension(:), allocatable :: E		! Energies for density of states calculation in cm^-1
		integer nE			! Number of energy grains
		integer dn			! Number of grains of E per grain of energies
		real(8), dimension(:), allocatable :: rho1		! Temporary array for density of states calculation
		real(8), dimension(:), allocatable :: rho2		! Temporary array for density of states calculation
		real(8) E1
		integer i, j, r
        
		real(8) Eisom
		integer start, length
		integer vibCount, rotCount, hindCount, speCount
		
		! Conversion factor: cm^-1 per kJ/mol
		conv = 1000 / 6.022e23 / 6.626e-34 / 2.9979e10
		
		speCount = 0
		do i = 1, isom%numSpecies
			
			vibCount = 0
			rotCount = 0
			hindCount = 0
			do j = 1, size(isom%vibFreq(i,:))
				if (isom%vibFreq(i,j) > 0) vibCount = vibCount + 1
			end do
			do j = 1, size(isom%rotFreq(i,:))
				if (isom%rotFreq(i,j) > 0) rotCount = rotCount + 1
			end do
			do j = 1, size(isom%hindFreq(i,:))
				if (isom%hindFreq(i,j) > 0) hindCount = hindCount + 1
			end do
    
			if (vibCount > 0 .or. rotCount > 0 .or. hindCount > 0) then
		
				speCount = speCount + 1
				
				! Get energy range of density of states
				Emin = 0
				Emax = simData%Emax - simData%Emin
				dE = simData%dE
				dn = 10
				
				! Convert Emin, Emax, and dE to units of cm^-1
				Emax = Emax * conv
				Emin = Emin * conv
                dE = dE * conv / dn
				
				! Create energy grains for density of states calculation
				nE = ceiling((Emax - Emin) / dE) + 1
				allocate( E(1:nE) )
				do j = 1, nE
					E(j) = (j - 1) * dE + Emin
				end do
        
				! Create temporary density of states arrays
				if (.not. allocated(rho1)) allocate( rho1(1:nE) )
				allocate( rho2(1:nE) )
		
				! Calculate the density of states for each multimolecular well
				rho2 = states(E + E(1), isom%vibFreq(1,:), isom%rotFreq(1,:), &
					isom%hindFreq(1,:), isom%hindBarrier(1,:), isom%symmNum(1))
				if (speCount == 1) then
					rho1 = rho2
				else
					rho1 = convolve(rho1, rho2, dE)
				end if
				
				deallocate( rho2, E )
				
			end if
		end do
		
		! Remove intermediate grains used only for density of states estimation
		rho1(1:simData%nGrains) = rho1(1:size(rho1):dn)
		
		! Shift density of states to appropriate energy range
		Eisom = sum(isom%E)
		start = 0
		do r = 1, simData%nGrains
			if (start == 0 .and. Eisom < simData%E(r)) then
				Eisom = simData%E(r)
				start = r
			end if
		end do
		length = simData%nGrains - start + 1

		if (.not. allocated(rho1)) then
			write (*,*) 'ERROR: Unable to calculate density of states. Stopping.'
			stop
		end if
		
		isom%densStates = 0 * isom%densStates
		isom%densStates(start:simData%nGrains) = rho1(1:length)
		
		deallocate( rho1 )
		
	end subroutine

	! --------------------------------------------------------------------------
	!
	! Function: states()
	! 
	! Estimates the density of states from the provided information.
	!
	! Parameters:
	!   E - The energy grains in cm^-1.
	!   vibFreq - The harmonic oscillator vibrational frequencies in cm^-1.
	!   rotFreq - The free rotor rotational frequencies in cm^-1.
	!   hindFreq - The hindered rotor rotational frequencies in cm^-1.
	!   hindBarrier - The hindered rotor barrier heights in cm^-1.
	!   symmNum - The vibrational frequencies in cm^-1.
	!
	! Returns:
	!	rho - The density of states in (cm^-1)^-1.
	function states(E, vibFreq, rotFreq, hindFreq, hindBarrier, symmNum)
		
		! Provide parameter type checking of inputs and outputs
		real(8), dimension(:), intent(in)		:: 	E
		real(8), dimension(:), intent(in)		:: 	vibFreq
		real(8), dimension(:), intent(in)		:: 	rotFreq
		real(8), dimension(:), intent(in)		:: 	hindFreq
		real(8), dimension(:), intent(in)		:: 	hindBarrier
		integer, intent(in)						:: 	symmNum
		real(8), dimension(1:size(E))			:: 	states
		
		real(8) dE
		real(8), dimension(:), allocatable		:: rho_RR
		real(8), dimension(:), allocatable		:: rho_HR
		real(8), dimension(:), allocatable		:: rho1
		real(8), dimension(:), allocatable		:: rho
		integer N
		integer i, j, k, found
		
		N = size(E)
		dE = E(2) - E(1)
		
		! Free rotors
		if (rotFreq(1) > 0) then
			allocate( rho_RR(1:N) )
			rho_RR = rho_RR * 0
			rho_RR = calcFreeRotorStates(E, rotFreq, symmNum)
		end if
		
        ! Hindered rotors
		if (hindFreq(1) > 0) then
			allocate( rho1(1:N) )
			do i = 1, size(hindFreq)
				if (hindFreq(i) > 0 .and. hindBarrier(i) > 0) then
					rho1 = calcHinderedRotorStates(E, hindFreq(i), hindBarrier(i))
					if (allocated(rho_HR)) then
						rho_HR = convolve(rho_HR, rho1, dE)
					else
						allocate( rho_HR(1:N) )
						rho_HR = rho1
					end if
				end if
			end do
		end if
		
		! Total rotational component of density of states
		if (allocated(rho_RR) .and. allocated(rho_HR)) then
			states = convolve(rho_RR, rho_HR, dE)
			deallocate(rho_RR, rho_HR)
		elseif (allocated(rho_RR)) then
			states = rho_RR
			deallocate(rho_RR)
		elseif (allocated(rho_HR)) then
			states = rho_HR
			deallocate(rho_HR)
		else
			states = 0 * states
		end if
        
		! Add vibrational states via Beyer-Swinehart algorithm
        call beyerSwinehart(E, vibFreq, states)
        
		! For this approach the density of states must be smooth
		i = 1
        do while (i <= size(E) - 1)
			if (states(i+1) == 0) then
				j = i + 2
				found = 1
				do while (j <= size(E) .and. found == 1)
					if (states(j) /= 0) then
						found = 0
					else
						j = j + 1
					end if
				end do
				if (found == 0) then
					do k = i + 1, j - 1	
						states(k) = (states(j) - states(i)) / (j - i) * (k - i) + states(i)
					end do
				else
					do k = i + 1, size(E)
						states(k) = states(k-1) * states(k-2) / states(k-3)
					end do
				end if
				i = j
			else
				i = i + 1
			end if
		end do
				
	end function
	
	! --------------------------------------------------------------------------
	!
	! Subroutine: convolve()
	! 
	! Convolves two density of states vectors.
	!
	! Parameters:
	!   rho1 - The first density of states in (cm^-1)^-1.
	!   rho2 - The second density of states in (cm^-1)^-1.
	!   dE - The grain size in cm^-1.
	!
	! Returns:
	!	rho1 - The convolved density of states in (cm^-1)^-1.
	function convolve(rho1, rho2, dE)
		
		! Provide parameter type checking of inputs and outputs
		real(8), dimension(:), intent(in)		:: 	rho1
		real(8), dimension(:), intent(in)		:: 	rho2
		real(8), intent(in)						:: 	dE
		real(8), dimension(1:min(size(rho1),size(rho2)))	:: 	convolve
		
		real(8), dimension(:), allocatable		:: 	rho
		real(8), dimension(:), allocatable		:: 	f
		integer i, j
        integer length
		
		
		length = min(size(rho1), size(rho2))
		
		allocate( f(1:length) )
		
		convolve = 0 * convolve
		do i = 2, length
			do j = 1, i
				f(j) = rho2(i-j+1) * rho1(j) * dE
			end do
			convolve(i) = sum(f(1:i))
		end do
		
		deallocate(f)
		
	end function

	! --------------------------------------------------------------------------
	!
	! Function: calcFreeRotorStates()
	! 
	! Convolves two density of states vectors.
	!
	! Parameters:
	!   E - The energy grains in cm^-1.
	!   E0 - The ground-state energy in cm^-1.
	!   rotFreq - The free rotor rotational constants in cm^-1.
	!   symmNum - The symmetry number.
	!
	! Returns:
	!	rho - The free rotor density of states in (cm^-1)^-1.
	function calcFreeRotorStates(E, rotFreq, symmNum)
		
		! Provide parameter type checking of inputs and outputs
		real(8), dimension(:), intent(in)		:: 	E
		real(8), dimension(:), intent(in)		:: 	rotFreq
		integer, intent(in)						:: 	symmNum
		real(8), dimension(1:size(E))			:: 	calcFreeRotorStates
		
		real(8) dE
		integer nRot, i
		
		dE = E(2) - E(1)
		
		nRot = 0
		do i = 1, size(rotFreq)
			if (rotFreq(i) > 0) nRot = nRot + 1
		end do
		
		calcFreeRotorStates = 0 * calcFreeRotorStates
		do i = 1, size(E)
			if (nRot == 1) then
				calcFreeRotorStates(i) = 1.0 / symmNum / rotFreq(1)
			elseif (nRot == 2) then
				calcFreeRotorStates(i) = 2.0 / symmNum * sqrt(E(i)) / sqrt(product(rotFreq(1:nRot) * rotFreq(1)))
			elseif (nRot == 3) then
				calcFreeRotorStates(i) = 2.0 / symmNum * sqrt(E(i)) / sqrt(product(rotFreq(1:nRot)))
			end if
		end do
		
	end function
	
	! --------------------------------------------------------------------------
	!
	! Function: calcHinderedRotorStates()
	! 
	! Convolves two density of states vectors.
	!
	! Parameters:
	!   E - The energy grains in cm^-1.
	!   E0 - The ground-state energy in cm^-1.
	!   freq - The hindered rotor rotational constant in cm^-1.
	!   barrier - The hindered rotor barrier height in cm^-1.
	!
	! Returns:
	!	rho - The hindered rotor density of states in (cm^-1)^-1.
	function calcHinderedRotorStates(E, freq, barrier)
		
		! Provide parameter type checking of inputs and outputs
		real(8), dimension(:), intent(in)		:: 	E
		real(8), intent(in)						:: 	freq
		real(8), intent(in)						:: 	barrier
		real(8), dimension(1:size(E))			:: 	calcHinderedRotorStates
		
		integer i, n
		real(8) m, K, tol
		
		tol = 1.0e-7
		
		do i = 1, size(E)
			if (E(i) < 0) then
				calcHinderedRotorStates(i) = 0
			elseif (E(i) < freq/2.d0) then
				m = 0.d0
				call celliptic(m, tol, K, n)
				calcHinderedRotorStates(i) = K / sqrt(barrier)
			elseif (E(i) - freq/2.d0 < barrier) then
				m = sqrt((E(i) - freq/2.d0) / barrier)
				call celliptic(m, tol, K, n)
				calcHinderedRotorStates(i) = K / sqrt(barrier)
			elseif (E(i) - freq/2.d0 == barrier) then
				m = 1.d0 - tol
				call celliptic(m, tol, K, n)
				calcHinderedRotorStates(i) = (K / sqrt(barrier) +  K/sqrt(E(i) - freq/2)) / 2.d0
			else
				m = sqrt(barrier / (E(i) - freq/2.d0))
				call celliptic(m, tol, K, n)
				calcHinderedRotorStates(i) = K / sqrt(E(i) - freq/2.d0)
			end if
		end do
		
	end function
	
	! --------------------------------------------------------------------------
	!
	! Subroutine: celliptic()
	! 
	! Evaluates the complete elliptic integral of the first kind. This is needed
	! by the hindered rotor states calculation algorithm.
	!
	! Parameters:
	!   m - The value at which to evaluate the integral; 0 <= m < 1
	!   tol - The desired uncertainty tolerance
	!
	! Returns:
	!   K - The value of the integral evaluated at m with tolerance tol.
	!	n - The number of iterations needed to evaulate the integral.
	
	subroutine celliptic(m, tol, K, n)

		! Provide parameter type checking of inputs and outputs
		real(8), intent(in)		:: 	m
		real(8), intent(in)		:: 	tol
		real(8), intent(out)	:: 	K
		integer, intent(out)	:: 	n
		
		real(8) A0, B0, A, B
		real(8) pi
		
		pi = 4.d0 * datan(1.d0)

		A = 1 + m 
		B = 1 - m
		n = 0;

		! Skip if m is outside of range (0, 1) or tolerance is negative
		if (m < 0 .or. m > 1 .or. tol < 0) then
			return
		end if
		
		do while (abs(A - B) > tol)
			n = n + 1
			! Generate improved values
			A0 = A
			B0 = B
			A = (A0 + B0) / 2.d0;
			B = sqrt(A0 * B0);
		end do

		! Result of integral
		K = pi / 2.d0 / A;
		
	end subroutine
	
	! --------------------------------------------------------------------------
	!
	! Subroutine: beyerSwinehart()
	! 
	! Convolves vibrational modes into a density of states vector using the
	! Beyer-Swinehart algorithm.
	!
	! Parameters:
	!   E - The energy grains in cm^-1.
	!   E0 - The ground-state energy in cm^-1.
	!   vibFreq - The harmonic oscillator vibrational constants in cm^-1.
	!
	! Returns:
	!	rho - The convolved density of states in (cm^-1)^-1.
	subroutine beyerSwinehart(E, vibFreq, rho)
		
		! Provide parameter type checking of inputs and outputs
		real(8), dimension(:), intent(in)		:: 	E
		real(8), dimension(:), intent(in)		:: 	vibFreq
		real(8), dimension(:), intent(inout)	:: 	rho
		
		real(8) dE
		integer nVib, dGrain, i, j
		
		dE = E(2) - E(1)
		
        ! Count nonzero vibrational modes
		nVib = 0
		do i = 1, size(vibFreq)
			if (vibFreq(i) > 0) nVib = nVib + 1
		end do
        
		! The Beyer-Swinehart algorithm
		rho(1) = 1.0 / dE
		do i = 1, nVib
			dGrain = nint(vibFreq(i) / dE)
			do j = 1, size(E)-dGrain
				rho(j + dGrain) = rho(j + dGrain) + rho(j)
			end do
		end do
        
	end subroutine

	! --------------------------------------------------------------------------

end module
