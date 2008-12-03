! ==============================================================================
!
!	DensityOfStates.f90
!
!	Written by Josh Allen (jwallen@mit.edu)
!
! ==============================================================================

module DensityOfStatesModule

	use SpeciesModule

contains

	! --------------------------------------------------------------------------
	!
	! Subroutine: calcDensityOfStates()
	! 
	! Estimates the density of states for each isomer in the network.
	!
	! Parameters:
	!   energies - The energy grains in kJ/mol (from simData%E).
	!   uniData - The chemical data about each unimolecular isomer.
	!   multiData - The chemical data about each multimolecular source/sink.
	subroutine calcDensityOfStates(energies, uniData, multiData)
		
		! Provide parameter type checking of inputs and outputs
		real(8), dimension(:), intent(in)				:: 	energies
		type(UniWell), dimension(:), intent(inout)		:: 	uniData
		type(MultiWell), dimension(:), intent(inout)	:: 	multiData
		
		real(8) conv		! Conversion factor: cm^-1 per kJ/mol
		real(8) Emax		! Maximum energy of range in cm^-1
		real(8) dE			! Grain size for density of states calculation in cm^-1
		real(8), dimension(:), allocatable :: E		! Energies for density of states calculation in cm^-1
		integer nE			! Number of energy grains
		integer dn			! Number of grains of E per grain of energies
		real(8), dimension(:), allocatable :: rho1		! Temporary array for density of states calculation
		real(8), dimension(:), allocatable :: rho2		! Temporary array for density of states calculation
		real(8) E1
		
		conv = 1000 / 6.022e23 / 6.626e-34 / 2.9979e10
		
		Emax = ceiling(maxval(energies) * conv)
		dE = 0.1 * conv
		nE = Emax / dE + 1
		allocate( E(1:nE) )
		do i = 1, nE
			E(i) = (i - 1) * dE
		end do
		allocate( rho1(1:nE) )
		allocate( rho2(1:nE) )
		
		dn = anint((energies(2) - energies(1)) / dE * conv)
		
		! Calculate the density of states for each unimolecular well
		do i = 1, size(uniData)
			call densityOfStates(E, &
				uniData(i)%E * conv, &
				uniData(i)%vibFreq, &
				uniData(i)%rotFreq, &
				uniData(i)%hindFreq, &
				uniData(i)%hindBarrier, &
				uniData(i)%symmNum, &
				rho1)
			uniData(i)%densStates = rho1(1:nE:dn) / conv
		end do
		
		! Calculate the density of states for each multimolecular well
		do n = 1, size(multiData)
			call densityOfStates(E, &
				multiData(n)%E(1) * conv, &
				multiData(n)%vibFreq(1,:), &
				multiData(n)%rotFreq(1,:), &
				multiData(n)%hindFreq(1,:), &
				multiData(n)%hindBarrier(1,:), &
				multiData(n)%symmNum(1), &
				rho1)
			E1 = multiData(n)%E(1)
			do i = 2, multiData(n)%numSpecies
				call densityOfStates(E, &
					multiData(n)%E(i) * conv, &
					multiData(n)%vibFreq(i,:), &
					multiData(n)%rotFreq(i,:), &
					multiData(n)%hindFreq(i,:), &
					multiData(n)%hindBarrier(i,:), &
					multiData(n)%symmNum(i), &
					rho2)
				call convolveStates(rho1, E1 * conv, rho2, multiData(n)%E(i) * conv, dE)
				E1 = E1 + multiData(n)%E(i)
			end do
			multiData(n)%densStates = rho1(1:nE:dn) / conv
		end do
		
		deallocate (rho1, rho2, E)
		
	end subroutine

	! --------------------------------------------------------------------------
	!
	! Subroutine: densityOfStates()
	! 
	! Estimates the density of states from the provided information.
	!
	! Parameters:
	!   E - The energy grains in cm^-1.
	!   E0 - The ground-state energy in cm^-1.
	!   vibFreq - The harmonic oscillator vibrational frequencies in cm^-1.
	!   rotFreq - The free rotor rotational frequencies in cm^-1.
	!   hindFreq - The hindered rotor rotational frequencies in cm^-1.
	!   hindBarrier - The hindered rotor barrier heights in cm^-1.
	!   symmNum - The vibrational frequencies in cm^-1.
	!
	! Returns:
	!	rho - The density of states in (cm^-1)^-1.
	subroutine densityOfStates(E, E0, vibFreq, rotFreq, hindFreq, hindBarrier, symmNum, rho)
		
		! Provide parameter type checking of inputs and outputs
		real(8), dimension(:), intent(in)		:: 	E
		real(8), intent(in)						:: 	E0
		real(8), dimension(:), intent(in)		:: 	vibFreq
		real(8), dimension(:), intent(in)		:: 	rotFreq
		real(8), dimension(:), intent(in)		:: 	hindFreq
		real(8), dimension(:), intent(in)		:: 	hindBarrier
		integer, intent(in)						:: 	symmNum
		real(8), dimension(:), intent(out)		:: 	rho
		
		real(8) dE
		real(8), dimension(:), allocatable		:: rho_RR
		real(8), dimension(:), allocatable		:: rho_HR
		integer N
		
		N = size(E)
		dE = E(2) - E(1)
		
		! Free rotors
		if (size(rotFreq) > 0) then
			allocate( rho_RR(1:N) )
			rho_RR = 0 * rho_RR
			call calcFreeRotorStates(E, E0, rotFreq, symmNum, rho_RR)
		end if
		
		! Hindered rotors
		if (size(hindFreq) > 0) then
			allocate( rho_HR(1:N) )
			rho_HR = 0 * rho_HR
			call calcHinderedRotorStates(E, E0, hindFreq(1), hindBarrier(1), symmNum, rho_HR)
		end if
		
		! Total rotational component of density of states
		if (allocated(rho_RR) .and. allocated(rho_HR)) then
			call convolveStates(rho_RR, E0, rho_HR, E0, dE)
			rho = rho_RR
			deallocate(rho_RR, rho_HR)
		elseif (allocated(rho_RR)) then
			rho = rho_RR
			deallocate(rho_RR)
		elseif (allocated(rho_HR)) then
			rho = rho_HR
			deallocate(rho_HR)
		else
			rho = 0 * rho
		end if
		
		! Add vibrational states via Beyer-Swinehart algorithm
		call beyerSwinehart(E, E0, vibFreq, rho)		
		
	end subroutine
	
	! --------------------------------------------------------------------------
	!
	! Subroutine: convolveStates()
	! 
	! Convolves two density of states vectors.
	!
	! Parameters:
	!   rho1 - The first density of states in (cm^-1)^-1.
	!   E1 - The first ground-state energy in cm^-1.
	!   rho2 - The second density of states in (cm^-1)^-1.
	!   E2 - The second ground-state energy in cm^-1.
	!   dE - The grain size in cm^-1.
	!
	! Returns:
	!	rho1 - The convolved density of states in (cm^-1)^-1.
	subroutine convolveStates(rho1, E1, rho2, E2, dE)
		
		! Provide parameter type checking of inputs and outputs
		real(8), dimension(:), intent(inout)	:: 	rho1
		real(8), intent(in)						:: 	E1
		real(8), dimension(:), intent(in)		:: 	rho2
		real(8), intent(in)						:: 	E2
		real(8), intent(in)						:: 	dE
		
		integer n1, n2
		real(8), dimension(:), allocatable		:: 	rho
		real(8), dimension(:), allocatable		:: 	f
		
		n1 = ceiling(E1 / dE) + 1
		n2 = ceiling(E2 / dE) + 1
		
		allocate( f(1:size(rho1)) )
		allocate( rho(1:size(rho1)) )
		
		do i = n2, size(rho1)
			f = 0 * f
			do j = n1, i
				f(j) = rho2(i-j+1) * rho1(j) * dE
			end do
			rho(i) = sum(f(n1:i))
		end do
		rho1 = rho
		
		deallocate(f, rho)
		
	end subroutine

	! --------------------------------------------------------------------------
	!
	! Subroutine: calcFreeRotorStates()
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
	subroutine calcFreeRotorStates(E, E0, rotFreq, symmNum, rho)
		
		! Provide parameter type checking of inputs and outputs
		real(8), dimension(:), intent(in)		:: 	E
		real(8), intent(in)						:: 	E0
		real(8), dimension(:), intent(in)		:: 	rotFreq
		integer, intent(in)						:: 	symmNum
		real(8), dimension(:), intent(out)		:: 	rho
		
		real(8) dE
		integer start
		integer nRot
		
		dE = E(2) - E(1)
		start = ceiling(E0 / dE) + 1
		
		nRot = 0
		do i = 1, size(rotFreq)
			if (rotFreq(i) > 0) nRot = nRot + 1
		end do
		
		rho = 0 * rho
		do i = start, size(E)
			if (nRot == 1) then
				rho(i) = 1.0 / symmNum / rotFreq(1)
			elseif (nRot == 2) then
				rho(i) = 2.0 / symmNum * sqrt(E(i) - E0) / sqrt(product(rotFreq(1:nRot) * rotFreq(1)))
			elseif (nRot == 3) then
				rho(i) = 2.0 / symmNum * sqrt(E(i) - E0) / sqrt(product(rotFreq(1:nRot)))
			end if
		end do
		
	end subroutine
	
	! --------------------------------------------------------------------------
	!
	! Subroutine: calcHinderedRotorStates()
	! 
	! Convolves two density of states vectors.
	!
	! Parameters:
	!   E - The energy grains in cm^-1.
	!   E0 - The ground-state energy in cm^-1.
	!   freq - The hindered rotor rotational constant in cm^-1.
	!   barrier - The hindered rotor barrier height in cm^-1.
	!   symmNum - The symmetry number.
	!
	! Returns:
	!	rho - The hindered rotor density of states in (cm^-1)^-1.
	subroutine calcHinderedRotorStates(E, E0, freq, barrier, symmNum, rho)
		
		! Provide parameter type checking of inputs and outputs
		real(8), dimension(:), intent(in)		:: 	E
		real(8), intent(in)						:: 	E0
		real(8), intent(in)						:: 	freq
		real(8), intent(in)						:: 	barrier
		integer, intent(in)						:: 	symmNum
		real(8), dimension(:), intent(out)		:: 	rho
		
		integer i, n
		real(8) m, K, tol
		
		tol = 1.0e-7
		
		do i = 1, size(E)
			if (E(i) - E0 < freq/2.d0) then
				m = 0.d0
				call celliptic(m, tol, K, n)
				rho(i) = K / sqrt(barrier)
			elseif (E(i) - E0 - freq/2.d0 < barrier) then
				m = sqrt((E(i) - E0 - freq/2.d0) / barrier)
				call celliptic(m, tol, K, n)
				rho(i) = K / sqrt(barrier)
			elseif (E(i) - E0 - freq/2.d0 == barrier) then
				m = 1.d0 - tol
				call celliptic(m, tol, K, n)
				rho(i) = (K / sqrt(barrier) +  K/sqrt(E(i) - E0 - freq/2)) / 2.d0
			else
				m = sqrt(barrier / (E(i) - E0 - freq/2.d0))
				call celliptic(m, tol, K, n)
				rho(i) = K / sqrt(E(i) - E0 - freq/2.d0)
			end if
		end do
		
	end subroutine
	
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
	subroutine beyerSwinehart(E, E0, vibFreq, rho)
		
		! Provide parameter type checking of inputs and outputs
		real(8), dimension(:), intent(in)		:: 	E
		real(8), intent(in)						:: 	E0
		real(8), dimension(:), intent(in)		:: 	vibFreq
		real(8), dimension(:), intent(inout)	:: 	rho
		
		real(8) dE
		integer start
		integer nVib
		integer grain
		
		dE = E(2) - E(1)
		start = ceiling(E0 / dE) + 1
		
		! Count nonzero vibrational modes
		nVib = 0
		do i = 1, size(vibFreq)
			if (vibFreq(i) > 0) nVib = nVib + 1
		end do
		
		! The Beyer-Swinehart algorithm
		rho(start) = 1 / dE
		do i = 1, nVib
			do j = start, size(E)
				grain = j + vibFreq(i) / dE			! Destination bin
				if (grain <= size(E)) then
					rho(grain) = rho(grain) + rho(j)
				end if
			end do
		end do
		
	end subroutine

	! --------------------------------------------------------------------------

end module
