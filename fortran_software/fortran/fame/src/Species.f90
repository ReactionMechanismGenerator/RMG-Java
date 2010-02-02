!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	species.f90
!
!	Written by Josh Allen (jwallen@mit.edu)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module SpeciesModule

	implicit none
	
	type GeneralData
		real(8)	molWt			! The molecular weight of the bath gas in kg/mol
		real(8)	sigma			! The Lennard-Jones sigma parameter of the bath gas in m
		real(8)	eps				! The Lennard-Jones epsilon parameter of the bath gas in J
	end type
	
	type SpectralData
		real(8), dimension(:), allocatable		:: 	vibFreq		! 1D harmonic oscillator frequencies in cm^-1
		real(8), dimension(:), allocatable		:: 	rotFreq		! 1D rigid rotor frequencies in cm^-1
		real(8), dimension(:), allocatable		:: 	hindFreq	! 1D hindered rotor frequencies in cm^-1
		real(8), dimension(:), allocatable		:: 	hindBarrier	! 1D hindered rotor barrier heights in cm^-1
		integer									:: 	symmNum		! Symmetry number
	end type
	
	type ThermoData
		real(8)								:: 	H298			! Standard enthalpy of formation in J/mol
		real(8)								:: 	S298			! Standard entropy of formation in J/mol*K
		real(8), dimension(:), allocatable	:: 	Cp				! Heat capacity table in J/mol*K
	end type
	
	type Species
		character(len=128)					::	name			! Species name
		real(8)								:: 	E0				! Ground-state electronic + zero-point energy in J/mol
		type(GeneralData)					::	general			! General gas data
		type(ThermoData)					::	thermo			! Thermodynamics data
		type(SpectralData)					::	spectral		! Spectroscopic data
	end type	

contains

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
	function collisionFrequency(T, P, bathGas, spec)
	
		! Provide parameter type-checking
		real(8), intent(in)					:: 	T
		real(8), intent(in)					:: 	P
		type(GeneralData), intent(in)		:: 	bathGas
		type(GeneralData), intent(in)		:: 	spec
		real(8)								::	collisionFrequency
		
		real(8)		::	kB
		real(8)		::	collisionIntegral
		real(8)		:: 	sigma
		real(8)		:: 	eps
		real(8)		:: 	mu
		real(8)		:: 	gasConc
		
		collisionIntegral = 1.16145 / T**0.14874 + 0.52487 / exp(0.77320 * T) + 2.16178 / exp(2.43787 * T) &
			-6.435/10000 * T**0.14874 * sin(18.0323 * T**(-0.76830) - 7.27371)
    
		kB = 1.3806504e-23
	
		gasConc = P / kB / T
		mu = 1 / (1/spec%molWt + 1/bathGas%molWt) / 6.022e23
		sigma = 0.5 * (spec%sigma + bathGas%sigma)
		eps = 0.5 * (spec%eps + bathGas%eps)
		
		! Evaluate collision frequency
		collisionFrequency = collisionIntegral * &
			sqrt(8 * kB * T / 3.141592654 / mu) * 3.141592654 * sigma**2 * gasConc
	
	end function	


	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!
	! Function: heatCapacity()
	!
	! Return the heat capacity corresponding to a given set of thermo data at a 
	! given temperature.
	!
	function heatCapacity(thermo, T)
	
		type(ThermoData), intent(in)	:: thermo
		real(8), intent(in)				:: T
		real(8)	heatCapacity
	
		real(8) slope, intercept
	
		if (T < 298.0) then
			write (0, fmt='(A)') 'Invalid temperature for heat capacity calculation.'
			stop
		elseif (T < 300.0) then
			heatCapacity = thermo%Cp(1)
		elseif (T < 400.0) then
			slope = (thermo%Cp(2) - thermo%Cp(1)) / (400.0 - 300.0)
			intercept = (thermo%Cp(1) * 400.0 - thermo%Cp(2) * 300.0) / (400.0 - 300.0)
			heatCapacity = intercept + slope * (T - 300.0)
		elseif (T < 500.0) then
			slope = (thermo%Cp(3) - thermo%Cp(2)) / (500.0 - 400.0)
			intercept = (thermo%Cp(2) * 500.0 - thermo%Cp(3) * 400.0) / (500.0 - 400.0)
			heatCapacity = intercept + slope * (T - 400.0)
		elseif (T < 600.0) then
			slope = (thermo%Cp(4) - thermo%Cp(3)) / (600.0 - 500.0)
			intercept = (thermo%Cp(3) * 600.0 - thermo%Cp(4) * 500.0) / (600.0 - 500.0)
			heatCapacity = intercept + slope * (T - 500.0)
		elseif (T < 800.0) then
			slope = (thermo%Cp(5) - thermo%Cp(4)) / (800.0 - 600.0)
			intercept = (thermo%Cp(4) * 800.0 - thermo%Cp(5) * 600.0) / (800.0 - 600.0)
			heatCapacity = intercept + slope * (T - 600.0)
		elseif (T < 1000.0) then
			slope = (thermo%Cp(6) - thermo%Cp(5)) / (1000.0 - 800.0)
			intercept = (thermo%Cp(5) * 1000.0 - thermo%Cp(6) * 800.0) / (1000.0 - 800.0)
			heatCapacity = intercept + slope * (T - 800.0)
		elseif (T < 1500.0) then
			slope = (thermo%Cp(7) - thermo%Cp(6)) / (1500.0 - 1000.0)
			intercept = (thermo%Cp(6) * 1500.0 - thermo%Cp(7) * 1000.0) / (1500.0 - 1000.0)
			heatCapacity = intercept + slope * (T - 1000.0)
		else
			heatCapacity = thermo%Cp(7)
		end if
	
	end function
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!
	! Function: enthalpy()
	!
	! Return the enthalpy corresponding to a given set of thermo data at a given
	! temperature.
	!
	function enthalpy(thermo, T)
	
		type(ThermoData), intent(in)	:: thermo
		real(8), intent(in)				:: T
		real(8)	enthalpy
	
		real(8) slope, intercept
	
		enthalpy = thermo%H298
		
		if (T < 298.0) then
			write (0, fmt='(A)') 'Invalid temperature for enthalpy calculation.'
			stop
		end if
		
		if (T > 300.0) then
			slope = (thermo%Cp(2) - thermo%Cp(1)) / (400.0 - 300.0)
			intercept = (thermo%Cp(1) * 400.0 - thermo%Cp(2) * 300.0) / (400.0 - 300.0)
			if (T < 400.0) then
				enthalpy = enthalpy + 0.5 * slope * (T**2 - 300.0**2) + intercept * (T - 300.0)
			else
				enthalpy = enthalpy + 0.5 * slope * (400.0**2 - 300.0**2) + intercept * (400.0 - 300.0)
			end if
		end if
		if (T > 400.0) then
			slope = (thermo%Cp(3) - thermo%Cp(2)) / (500.0 - 400.0)
			intercept = (thermo%Cp(2) * 500.0 - thermo%Cp(3) * 400.0) / (500.0 - 400.0)
			if (T < 500.0) then
				enthalpy = enthalpy + 0.5 * slope * (T**2 - 400.0**2) + intercept * (T - 400.0)
			else
				enthalpy = enthalpy + 0.5 * slope * (500.0**2 - 400.0**2) + intercept * (500.0 - 400.0)
			end if
		end if
		if (T > 500.0) then
			slope = (thermo%Cp(4) - thermo%Cp(3)) / (600.0 - 500.0)
			intercept = (thermo%Cp(3) * 600.0 - thermo%Cp(4) * 500.0) / (600.0 - 500.0)
			if (T < 600.0) then
				enthalpy = enthalpy + 0.5 * slope * (T**2 - 500.0**2) + intercept * (T - 500.0)
			else
				enthalpy = enthalpy + 0.5 * slope * (600.0**2 - 500.0**2) + intercept * (600.0 - 500.0)
			end if
		end if
		if (T > 600.0) then
			slope = (thermo%Cp(5) - thermo%Cp(4)) / (800.0 - 600.0)
			intercept = (thermo%Cp(4) * 800.0 - thermo%Cp(5) * 600.0) / (800.0 - 600.0)
			if (T < 800.0) then
				enthalpy = enthalpy + 0.5 * slope * (T**2 - 600.0**2) + intercept * (T - 600.0)
			else
				enthalpy = enthalpy + 0.5 * slope * (800.0**2 - 600.0**2) + intercept * (800.0 - 600.0)
			end if
		end if
		if (T > 800.0) then
			slope = (thermo%Cp(6) - thermo%Cp(5)) / (1000.0 - 800.0)
			intercept = (thermo%Cp(5) * 1000.0 - thermo%Cp(6) * 800.0) / (1000.0 - 800.0)
			if (T < 1000.0) then
				enthalpy = enthalpy + 0.5 * slope * (T**2 - 800.0**2) + intercept * (T - 800.0)
			else
				enthalpy = enthalpy + 0.5 * slope * (1000.0**2 - 800.0**2) + intercept * (1000.0 - 800.0)
			end if
		end if
		if (T > 1000.0) then
			slope = (thermo%Cp(7) - thermo%Cp(6)) / (1500.0 - 1000.0)
			intercept = (thermo%Cp(6) * 1500.0 - thermo%Cp(7) * 1000.0) / (1500.0 - 1000.0)
			if (T < 1500.0) then
				enthalpy = enthalpy + 0.5 * slope * (T**2 - 1000.0**2) + intercept * (T - 1000.0)
			else
				enthalpy = enthalpy + 0.5 * slope * (1500.0**2 - 1000.0**2) + intercept * (1500.0 - 1000.0)
			end if
		end if
		if (T > 1500.0) then
			slope = (thermo%Cp(7) - thermo%Cp(6)) / (1500.0 - 1000.0)
			intercept = (thermo%Cp(6) * 1500.0 - thermo%Cp(7) * 1000.0) / (1500.0 - 1000.0)
			enthalpy = enthalpy + 0.5 * slope * (T**2 - 1500.0**2) + intercept * (T - 1500.0)
		end if
	
	end function
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!
	! Function: entropy()
	!
	! Return the entropy corresponding to a given set of thermo data at a given
	! temperature.
	!
	function entropy(thermo, T)
	
		type(ThermoData), intent(in)	:: thermo
		real(8), intent(in)				:: T
		real(8)	entropy
	
		real(8) slope, intercept
	
		entropy = thermo%S298
		
		if (T < 298.0) then
			write (0, fmt='(A)') 'Invalid temperature for entropy calculation.'
			stop
		end if
		
		if (T > 300.0) then
			slope = (thermo%Cp(2) - thermo%Cp(1)) / (400.0 - 300.0)
			intercept = (thermo%Cp(1) * 400.0 - thermo%Cp(2) * 300.0) / (400.0 - 300.0)
			if (T < 400.0) then
				entropy = entropy + slope * (T - 300.0) + intercept * log(T / 300.0)
			else
				entropy = entropy + slope * (400.0 - 300.0) + intercept * log(400.0 / 300.0)
			end if
		end if
		if (T > 400.0) then
			slope = (thermo%Cp(3) - thermo%Cp(2)) / (500.0 - 400.0)
			intercept = (thermo%Cp(2) * 500.0 - thermo%Cp(3) * 400.0) / (500.0 - 400.0)
			if (T < 500.0) then
				entropy = entropy + slope * (T - 400.0) + intercept * log(T / 400.0)
			else
				entropy = entropy + slope * (500.0 - 400.0) + intercept * log(500.0 / 400.0)
			end if
		end if
		if (T > 500.0) then
			slope = (thermo%Cp(4) - thermo%Cp(3)) / (600.0 - 500.0)
			intercept = (thermo%Cp(3) * 600.0 - thermo%Cp(4) * 500.0) / (600.0 - 500.0)
			if (T < 600.0) then
				entropy = entropy + slope * (T - 500.0) + intercept * log(T / 500.0)
			else
				entropy = entropy + slope * (600.0 - 500.0) + intercept * log(600.0 / 500.0)
			end if
		end if
		if (T > 600.0) then
			slope = (thermo%Cp(5) - thermo%Cp(4)) / (800.0 - 600.0)
			intercept = (thermo%Cp(4) * 800.0 - thermo%Cp(5) * 600.0) / (800.0 - 600.0)
			if (T < 800.0) then
				entropy = entropy + slope * (T - 600.0) + intercept * log(T / 600.0)
			else
				entropy = entropy + slope * (800.0 - 600.0) + intercept * log(800.0 / 600.0)
			end if
		end if
		if (T > 800.0) then
			slope = (thermo%Cp(6) - thermo%Cp(5)) / (1000.0 - 800.0)
			intercept = (thermo%Cp(5) * 1000.0 - thermo%Cp(6) * 800.0) / (1000.0 - 800.0)
			if (T < 1000.0) then
				entropy = entropy + slope * (T - 800.0) + intercept * log(T / 800.0)
			else
				entropy = entropy + slope * (1000.0 - 800.0) + intercept * log(1000.0 / 800.0)
			end if
		end if
		if (T > 1000.0) then
			slope = (thermo%Cp(7) - thermo%Cp(6)) / (1500.0 - 1000.0)
			intercept = (thermo%Cp(6) * 1500.0 - thermo%Cp(7) * 1000.0) / (1500.0 - 1000.0)
			if (T < 1500.0) then
				entropy = entropy + slope * (T - 1000.0) + intercept * log(T / 1000.0)
			else
				entropy = entropy + slope * (1500.0 - 1000.0) + intercept * log(1500.0 / 1000.0)
			end if
		end if
		if (T > 1500.0) then
			slope = (thermo%Cp(7) - thermo%Cp(6)) / (1500.0 - 1000.0)
			intercept = (thermo%Cp(6) * 1500.0 - thermo%Cp(7) * 1000.0) / (1500.0 - 1000.0)
			entropy = entropy + slope * (T - 1500.0) + intercept * log(T / 1500.0)
		end if
	
	end function
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!
	! Function: freeEnergy()
	!
	! Return the Gibbs free energy corresponding to a given set of thermo data 
	! at a given temperature.
	!
	function freeEnergy(thermo, T)
	
		type(ThermoData), intent(in)	:: thermo
		real(8), intent(in)				:: T
		real(8)	freeEnergy
	
		freeEnergy = enthalpy(thermo, T) - T * entropy(thermo, T)
	
	end function

end module
