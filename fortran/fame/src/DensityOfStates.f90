!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	states.f90
!
!	Written by Josh Allen (jwallen@mit.edu)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module StatesModule

	use SpeciesModule
	use IsomerModule
	
	implicit none
	
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
		
		R = 8.314472
	
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
	
	! --------------------------------------------------------------------------
	!
	! Subroutine: densityOfStates()
	! 
	! Estimates the density of states for a given isomer in the network.
	!
	subroutine densityOfStates(isom, speciesList, E0)
		
		! Provide parameter type checking of inputs and outputs
		type(Isomer), intent(inout)					:: 	isom
		type(Species), dimension(:), intent(in)		:: 	speciesList
		real(8), dimension(:), intent(in)			::	E0
		
		type(Species) :: spec
		
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
		
		! Conversion factor: cm^-1 per J/mol
		conv = 1.0 / 6.022e23 / 6.626e-34 / 2.9979e10
		
		speCount = 0
		do i = 1, size(isom%species)
			
			spec = speciesList(isom%species(i))
			
			vibCount = size(spec%spectral%vibFreq)
			rotCount = size(spec%spectral%rotFreq)
			hindCount = size(spec%spectral%hindFreq)
			if (.not. allocated(spec%spectral%vibFreq)) vibCount = 0
			if (.not. allocated(spec%spectral%rotFreq)) rotCount = 0
			if (.not. allocated(spec%spectral%hindFreq)) hindCount = 0
			
			if (vibCount > 0 .or. rotCount > 0 .or. hindCount > 0) then
		
				speCount = speCount + 1
				
				! Get energy range of density of states
				Emin = 0
				Emax = maxval(E0) - minval(E0)
				dE = E0(2) - E0(1)
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
				rho2 = states(E + E(1), spec%spectral%vibFreq, spec%spectral%rotFreq, &
					spec%spectral%hindFreq, spec%spectral%hindBarrier, spec%spectral%symmNum)
				
				if (speCount == 1) then
					rho1 = rho2
				else
					rho1 = convolve(rho1, rho2, dE)
				end if
				
				deallocate( rho2, E )
				
			end if
		end do
		
		! Remove intermediate grains used only for density of states estimation
		rho1(1:size(E0)) = rho1(1:size(rho1):dn)
		
		! Shift density of states to appropriate energy range
		Eisom = isom%E0
		start = 0
		do r = 1, size(E0)
			if (start == 0 .and. Eisom < E0(r)) then
				Eisom = E0(r)
				start = r
			end if
		end do
		length = size(E0) - start + 1

		if (.not. allocated(rho1)) then
			write (*,*) 'ERROR: Unable to calculate density of states. Stopping.'
			stop
		end if
		
		allocate(isom%densStates(1:size(E0)))
		isom%densStates = 0 * isom%densStates
		isom%densStates(start:size(E0)) = rho1(1:length) * conv
		
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
		real(8), dimension(:), allocatable, intent(in)		:: 	vibFreq
		real(8), dimension(:), allocatable, intent(in)		:: 	rotFreq
		real(8), dimension(:), allocatable, intent(in)		:: 	hindFreq
		real(8), dimension(:), allocatable, intent(in)		:: 	hindBarrier
		integer, intent(in)						:: 	symmNum
		real(8), dimension(1:size(E))			:: 	states
		
		real(8) dE
		real(8), dimension(:), allocatable		:: rho0
		real(8), dimension(:), allocatable		:: rho
		integer N
		integer i, j, k, found
		
		N = size(E)
		dE = E(2) - E(1)
		
		allocate( rho0(1:N) )
		
		! Free rotors
		if (size(rotFreq) > 0) then
			allocate( rho(1:N) )
			rho = calcFreeRotorStates(E, rotFreq, symmNum)
		end if
		
        ! Hindered rotors
		if (size(hindFreq) > 0) then
			do i = 1, size(hindFreq)
				if (hindFreq(i) > 0 .and. hindBarrier(i) > 0) then
					rho0 = calcHinderedRotorStates(E, hindFreq(i), hindBarrier(i))
					if (allocated(rho)) then
						rho = convolve(rho, rho0, dE)
					else
						allocate( rho(1:N) )
						rho = rho0
					end if
				end if
			end do
		end if
		
		! Add vibrational states via Beyer-Swinehart algorithm
		if (.not. allocated(rho)) then
			allocate( rho(1:N) )
			rho = 0 * rho
		end if
		if (size(vibFreq) > 0) then
			call beyerSwinehart(E, vibFreq, rho)
		end if
        states = rho
		
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
		
		deallocate( rho0, rho )
		
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
		real(8) m, K, tol, factor, pi
		
		tol = 1.0e-7
		pi = 3.1415926535897932
		
		do i = 1, size(E)
			factor = 2.0 / sqrt(pi**3) * sqrt(pi / freq)
			
			if (E(i) < 0) then
				calcHinderedRotorStates(i) = 0
			elseif (E(i) <= barrier) then
				m = sqrt(E(i) / barrier)
				call celliptic(m, tol, K, n)
				calcHinderedRotorStates(i) = factor * K / sqrt(barrier)
			else
				m = sqrt(barrier / E(i))
				call celliptic(m, tol, K, n)
				calcHinderedRotorStates(i) = factor * K / sqrt(E(i))
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

		!A = 1 + m 
		!B = 1 - m
		A = 1
		B = sqrt(1 - m)
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
