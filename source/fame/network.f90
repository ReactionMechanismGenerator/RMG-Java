!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!    network.f90
!
!    Written by Josh Allen (jwallen@mit.edu)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

implicit none

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module NetworkModule

    type GeneralData
        real(8) :: molWt    ! The molecular weight of the bath gas in kg/mol
        real(8) :: sigma    ! The Lennard-Jones sigma parameter of the bath gas in m
        real(8) :: eps      ! The Lennard-Jones epsilon parameter of the bath gas in J
    end type

    type SpectralData
        real(8), dimension(:), allocatable :: vib        ! 1D harmonic oscillator frequencies in cm^-1
        real(8), dimension(:), allocatable :: rot        ! 1D rigid rotor frequencies in cm^-1
        real(8), dimension(:,:), allocatable :: hind     ! 1D hindered rotor frequencies and barrier heights in cm^-1
        integer ::     symm                              ! Symmetry number
    end type

    type ThermoData
        real(8) :: H298                             ! Standard enthalpy of formation in J/mol
        real(8) :: S298                             ! Standard entropy of formation in J/mol*K
        real(8), dimension(:), allocatable :: Cp    ! Heat capacity table in J/mol*K
    end type

    type Species
        character(len=128) :: name          ! Species name
        real(8) :: E0                       ! Ground-state electronic + zero-point energy in J/mol
        type(GeneralData) :: general        ! General gas data
        type(ThermoData) :: thermo          ! Thermodynamics data
        type(SpectralData) :: spectral      ! Spectroscopic data
    end type

    type Isomer
        integer, dimension(:), allocatable :: species       ! List of indices to species in the well
        real(8), dimension(:), allocatable :: densStates    ! Density of states in mol/J
        real(8) :: E0                                       ! Ground state energy in J/mol
        real(8) :: Q                                        ! Partition function
        real(8), dimension(:), allocatable :: eqDist        ! Equilibrium distribution
        real(8) :: collFreq                                 ! Collision frequency in Hz
        real(8), dimension(:,:), allocatable :: Mcoll       ! Collisional transfer probability matrix
    end type

    type ArrheniusKinetics
        real(8) :: A        ! Arrhenius preexponential factor in s^-1
        real(8) :: n        ! Arrhenius temperature exponent
        real(8) :: Ea       ! Arrhenius activation energy in kJ/mol
    end type

    type Reaction
        character(len=256) :: equation		        ! Chemical reaction equation
        integer	:: reac                             ! Isomer wells connected to transition state
        integer	:: prod                             ! Isomer wells connected to transition state
        real(8) :: E0                               ! Ground-state electronic + zero-point energy of transition state in kJ/mol
        type(ArrheniusKinetics) :: arrhenius        ! Arrhenius kinetics for high pressure limit
        real(8), dimension(:), allocatable :: kf    ! Microcanonical forward reaction rate in s^-1
        real(8), dimension(:), allocatable :: kb    ! Microcanonical backward reaction rate in s^-1
    end type

    type Network
        type(Species), dimension(:), allocatable :: species
        type(Isomer), dimension(:), allocatable :: isomers
        type(Reaction), dimension(:), allocatable :: reactions
        type(GeneralData) :: bathGas
    end type

contains

    function species_getHeatCapacity(thermo, T) result(Cp)

        type(ThermoData), intent(in)	:: thermo
        real(8), intent(in)				:: T
        real(8)	Cp

        real(8) slope, intercept

        if (T < 298.0) then
            write (0, fmt='(A)') 'Invalid temperature for heat capacity calculation.'
            stop
        elseif (T < 300.0) then
            Cp = thermo%Cp(1)
        elseif (T < 400.0) then
            slope = (thermo%Cp(2) - thermo%Cp(1)) / (400.0 - 300.0)
            intercept = (thermo%Cp(1) * 400.0 - thermo%Cp(2) * 300.0) / (400.0 - 300.0)
            Cp = intercept + slope * (T - 300.0)
        elseif (T < 500.0) then
            slope = (thermo%Cp(3) - thermo%Cp(2)) / (500.0 - 400.0)
            intercept = (thermo%Cp(2) * 500.0 - thermo%Cp(3) * 400.0) / (500.0 - 400.0)
            Cp = intercept + slope * (T - 400.0)
        elseif (T < 600.0) then
            slope = (thermo%Cp(4) - thermo%Cp(3)) / (600.0 - 500.0)
            intercept = (thermo%Cp(3) * 600.0 - thermo%Cp(4) * 500.0) / (600.0 - 500.0)
            Cp = intercept + slope * (T - 500.0)
        elseif (T < 800.0) then
            slope = (thermo%Cp(5) - thermo%Cp(4)) / (800.0 - 600.0)
            intercept = (thermo%Cp(4) * 800.0 - thermo%Cp(5) * 600.0) / (800.0 - 600.0)
            Cp = intercept + slope * (T - 600.0)
        elseif (T < 1000.0) then
            slope = (thermo%Cp(6) - thermo%Cp(5)) / (1000.0 - 800.0)
            intercept = (thermo%Cp(5) * 1000.0 - thermo%Cp(6) * 800.0) / (1000.0 - 800.0)
            Cp = intercept + slope * (T - 800.0)
        elseif (T < 1500.0) then
            slope = (thermo%Cp(7) - thermo%Cp(6)) / (1500.0 - 1000.0)
            intercept = (thermo%Cp(6) * 1500.0 - thermo%Cp(7) * 1000.0) / (1500.0 - 1000.0)
            Cp = intercept + slope * (T - 1000.0)
        else
            Cp = thermo%Cp(7)
        end if

    end function

    function species_getEnthalpy(thermo, T) result(H)

        type(ThermoData), intent(in)	:: thermo
        real(8), intent(in)				:: T
        real(8)	H

        real(8) slope, intercept

        H = thermo%H298

        if (T < 298.0) then
            write (0, fmt='(A)') 'Invalid temperature for enthalpy calculation.'
            stop
        end if

        if (T > 300.0) then
            slope = (thermo%Cp(2) - thermo%Cp(1)) / (400.0 - 300.0)
            intercept = (thermo%Cp(1) * 400.0 - thermo%Cp(2) * 300.0) / (400.0 - 300.0)
            if (T < 400.0) then
                H = H + 0.5 * slope * (T**2 - 300.0**2) + intercept * (T - 300.0)
            else
                H = H + 0.5 * slope * (400.0**2 - 300.0**2) + intercept * (400.0 - 300.0)
            end if
        end if
        if (T > 400.0) then
            slope = (thermo%Cp(3) - thermo%Cp(2)) / (500.0 - 400.0)
            intercept = (thermo%Cp(2) * 500.0 - thermo%Cp(3) * 400.0) / (500.0 - 400.0)
            if (T < 500.0) then
                H = H + 0.5 * slope * (T**2 - 400.0**2) + intercept * (T - 400.0)
            else
                H = H + 0.5 * slope * (500.0**2 - 400.0**2) + intercept * (500.0 - 400.0)
            end if
        end if
        if (T > 500.0) then
            slope = (thermo%Cp(4) - thermo%Cp(3)) / (600.0 - 500.0)
            intercept = (thermo%Cp(3) * 600.0 - thermo%Cp(4) * 500.0) / (600.0 - 500.0)
            if (T < 600.0) then
                H = H + 0.5 * slope * (T**2 - 500.0**2) + intercept * (T - 500.0)
            else
                H = H + 0.5 * slope * (600.0**2 - 500.0**2) + intercept * (600.0 - 500.0)
            end if
        end if
        if (T > 600.0) then
            slope = (thermo%Cp(5) - thermo%Cp(4)) / (800.0 - 600.0)
            intercept = (thermo%Cp(4) * 800.0 - thermo%Cp(5) * 600.0) / (800.0 - 600.0)
            if (T < 800.0) then
                H = H + 0.5 * slope * (T**2 - 600.0**2) + intercept * (T - 600.0)
            else
                H = H + 0.5 * slope * (800.0**2 - 600.0**2) + intercept * (800.0 - 600.0)
            end if
        end if
        if (T > 800.0) then
            slope = (thermo%Cp(6) - thermo%Cp(5)) / (1000.0 - 800.0)
            intercept = (thermo%Cp(5) * 1000.0 - thermo%Cp(6) * 800.0) / (1000.0 - 800.0)
            if (T < 1000.0) then
                H = H + 0.5 * slope * (T**2 - 800.0**2) + intercept * (T - 800.0)
            else
                H = H + 0.5 * slope * (1000.0**2 - 800.0**2) + intercept * (1000.0 - 800.0)
            end if
        end if
        if (T > 1000.0) then
            slope = (thermo%Cp(7) - thermo%Cp(6)) / (1500.0 - 1000.0)
            intercept = (thermo%Cp(6) * 1500.0 - thermo%Cp(7) * 1000.0) / (1500.0 - 1000.0)
            if (T < 1500.0) then
                H = H + 0.5 * slope * (T**2 - 1000.0**2) + intercept * (T - 1000.0)
            else
                H = H + 0.5 * slope * (1500.0**2 - 1000.0**2) + intercept * (1500.0 - 1000.0)
            end if
        end if
        if (T > 1500.0) then
            slope = (thermo%Cp(7) - thermo%Cp(6)) / (1500.0 - 1000.0)
            intercept = (thermo%Cp(6) * 1500.0 - thermo%Cp(7) * 1000.0) / (1500.0 - 1000.0)
            H = H + 0.5 * slope * (T**2 - 1500.0**2) + intercept * (T - 1500.0)
        end if

    end function

    function species_getEntropy(thermo, T) result(S)

        type(ThermoData), intent(in)	:: thermo
        real(8), intent(in)				:: T
        real(8)	S

        real(8) slope, intercept

        S = thermo%S298

        if (T < 298.0) then
            write (0, fmt='(A)') 'Invalid temperature for entropy calculation.'
            stop
        end if

        if (T > 300.0) then
            slope = (thermo%Cp(2) - thermo%Cp(1)) / (400.0 - 300.0)
            intercept = (thermo%Cp(1) * 400.0 - thermo%Cp(2) * 300.0) / (400.0 - 300.0)
            if (T < 400.0) then
                S = S + slope * (T - 300.0) + intercept * log(T / 300.0)
            else
                S = S + slope * (400.0 - 300.0) + intercept * log(400.0 / 300.0)
            end if
        end if
        if (T > 400.0) then
            slope = (thermo%Cp(3) - thermo%Cp(2)) / (500.0 - 400.0)
            intercept = (thermo%Cp(2) * 500.0 - thermo%Cp(3) * 400.0) / (500.0 - 400.0)
            if (T < 500.0) then
                S = S + slope * (T - 400.0) + intercept * log(T / 400.0)
            else
                S = S + slope * (500.0 - 400.0) + intercept * log(500.0 / 400.0)
            end if
        end if
        if (T > 500.0) then
            slope = (thermo%Cp(4) - thermo%Cp(3)) / (600.0 - 500.0)
            intercept = (thermo%Cp(3) * 600.0 - thermo%Cp(4) * 500.0) / (600.0 - 500.0)
            if (T < 600.0) then
                S = S + slope * (T - 500.0) + intercept * log(T / 500.0)
            else
                S = S + slope * (600.0 - 500.0) + intercept * log(600.0 / 500.0)
            end if
        end if
        if (T > 600.0) then
            slope = (thermo%Cp(5) - thermo%Cp(4)) / (800.0 - 600.0)
            intercept = (thermo%Cp(4) * 800.0 - thermo%Cp(5) * 600.0) / (800.0 - 600.0)
            if (T < 800.0) then
                S = S + slope * (T - 600.0) + intercept * log(T / 600.0)
            else
                S = S + slope * (800.0 - 600.0) + intercept * log(800.0 / 600.0)
            end if
        end if
        if (T > 800.0) then
            slope = (thermo%Cp(6) - thermo%Cp(5)) / (1000.0 - 800.0)
            intercept = (thermo%Cp(5) * 1000.0 - thermo%Cp(6) * 800.0) / (1000.0 - 800.0)
            if (T < 1000.0) then
                S = S + slope * (T - 800.0) + intercept * log(T / 800.0)
            else
                S = S + slope * (1000.0 - 800.0) + intercept * log(1000.0 / 800.0)
            end if
        end if
        if (T > 1000.0) then
            slope = (thermo%Cp(7) - thermo%Cp(6)) / (1500.0 - 1000.0)
            intercept = (thermo%Cp(6) * 1500.0 - thermo%Cp(7) * 1000.0) / (1500.0 - 1000.0)
            if (T < 1500.0) then
                S = S + slope * (T - 1000.0) + intercept * log(T / 1000.0)
            else
                S = S + slope * (1500.0 - 1000.0) + intercept * log(1500.0 / 1000.0)
            end if
        end if
        if (T > 1500.0) then
            slope = (thermo%Cp(7) - thermo%Cp(6)) / (1500.0 - 1000.0)
            intercept = (thermo%Cp(6) * 1500.0 - thermo%Cp(7) * 1000.0) / (1500.0 - 1000.0)
            S = S + slope * (T - 1500.0) + intercept * log(T / 1500.0)
        end if

    end function

    function species_getFreeEnergy(thermo, T) result(G)

        type(ThermoData), intent(in)	:: thermo
        real(8), intent(in)				:: T
        real(8)	G

        G = species_getEnthalpy(thermo, T) - T * species_getEntropy(thermo, T)

    end function

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    function isomer_collisionFrequency(T, P, bathGas, spec)

        ! Provide parameter type-checking
        real(8), intent(in)	:: T
        real(8), intent(in)	:: P
        type(GeneralData), intent(in) :: bathGas
        type(GeneralData), intent(in) :: spec
        real(8) :: collisionFrequency

        real(8)		::	kB, collisionIntegral, sigma, eps, mu, gasConc

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

    subroutine isomer_eqDist(isom, E, T)

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

        write (1,*) '#DEBUG: Sum of un-normalised equilibrium distribution is:', sum(isom%eqDist)
        write (1,*) '#DEBUG: Un-normalised Equilibrium distribution is',size(isom%eqDist),'long and starts:',isom%eqDist(1:10)
        write (1,*) '#DEBUG: Denisty of states is',size(isom%densStates ),'long and starts:',isom%densStates (1:10)
        write (1,*) '#DEBUG: E grains are',size(E),'long and starts:',E(1:10)

        isom%eqDist = isom%eqDist / sum(isom%eqDist)

    end subroutine

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine reaction_calculateMicrocanonicalRates(rxn, T, E, isomerList, speciesList)

        type(Reaction), intent(inout) :: rxn
        real(8), intent(in) :: T
        real(8), dimension(:), intent(in) :: E
        type(Isomer), dimension(:), intent(in) :: isomerList
        type(Species), dimension(:), intent(in) :: speciesList

        type(Isomer) :: reac, prod
        integer nGrains, r
        real(8) dE, Keq

        nGrains = size(E)
        dE = E(2) - E(1)

        ! Allocate and zero rate coefficients
        if (.not. allocated(rxn%kf)) then
            allocate( rxn%kf(1:nGrains), rxn%kb(1:nGrains) )
        end if
        rxn%kf = 0 * rxn%kf
        rxn%kb = 0 * rxn%kb

        reac = isomerList(rxn%reac)
        prod = isomerList(rxn%prod)

        if (size(reac%species) == 1 .and. size(prod%species) == 1) then

            ! Calculate forward rate coefficient via inverse Laplace transform
            call rateILT(rxn%E0 - reac%E0, reac%densStates, &
                rxn%arrhenius, T, E, rxn%kf)

            ! Calculate backward rate coefficient via detailed balance
            do r = 1, nGrains
                if (prod%densStates(r) /= 0) then
                    rxn%kb(r) = rxn%kf(r) * reac%densStates(r) / prod%densStates(r)
                end if
            end do

        elseif (size(reac%species) == 1 .and. size(prod%species) > 1) then

            ! Calculate forward rate coefficient via inverse Laplace transform
            call rateILT(rxn%E0 - reac%E0, reac%densStates, &
                rxn%arrhenius, T, E, rxn%kf)

            ! Calculate equilibrium constant
            Keq = equilCoeff(reac, prod, T, speciesList)

            ! Calculate backward rate coefficient via detailed balance
            ! assuming multimolecular isomer is fully thermalized
            do r = 1, nGrains
                rxn%kb(r) = rxn%kf(r) / Keq * &
                    reac%densStates(r) * exp(-E(r) / 8.314472 / T) / reac%Q * dE
            end do

        elseif (size(reac%species) > 1 .and. size(prod%species) == 1) then

            ! Convert to dissocation so that a similar algorithm to above can be used
            reac = isomerList(rxn%prod)
            prod = isomerList(rxn%reac)

            ! Calculate forward rate coefficient via inverse Laplace transform
            call rateILT(rxn%E0 - reac%E0, reac%densStates, &
                rxn%arrhenius, T, E, rxn%kb)

            ! Calculate equilibrium constant
            Keq = equilCoeff(reac, prod, T, speciesList)

            ! Calculate backward rate coefficient via detailed balance
            ! assuming multimolecular isomer is fully thermalized
            do r = 1, nGrains
                rxn%kf(r) = rxn%kb(r) / Keq * &
                    reac%densStates(r) * exp(-E(r) / 8.314472 / T) / reac%Q * dE
            end do

        end if

        !write (*,*) rxn%isomerList(1), rxn%isomerList(2), rxn%kf(nGrains), rxn%kb(nGrains)

    end subroutine

    subroutine reaction_kineticsILT(E0, rho, kinetics, T, E, k)

        ! Provide parameter type-checking
        real(8), intent(in)					::	E0
        real(8), dimension(:), intent(in)	:: 	rho
        type(ArrheniusKinetics), intent(in)	::	kinetics
        real(8), intent(in)					::	T
        real(8), dimension(:), intent(in)	:: 	E
        real(8), dimension(:), intent(out)	:: 	k

        real(8) dE
        integer		:: r, s

        dE = E(2) - E(1)
        s = floor(kinetics%Ea / dE)
        if (s < 0) s = 0

        ! Determine rate coefficients using inverse Laplace transform
        do r = 1, size(E)
            if (s >= r .or. rho(r) == 0) then
                k(r) = 0
            else
                k(r) = kinetics%A * (T ** kinetics%n) * rho(r - s) / rho(r)
            end if
        end do

    end subroutine

    function reaction_equilCoeff(reac, prod, T, speciesList)

        type(Isomer), intent(in) 	::	reac
        type(Isomer), intent(in) 	::	prod
        real(8), intent(in)			::	T
        type(Species), dimension(:), intent(in) :: speciesList
        real(8) 					:: 	reaction_equilCoeff

        real(8) dGrxn
        integer i

        ! Determine free energy of reaction
        dGrxn = 0.0
        do i = 1, size(reac%species)
            dGrxn = dGrxn - freeEnergy(speciesList(reac%species(i))%thermo, T)
        end do
        do i = 1, size(prod%species)
            dGrxn = dGrxn + freeEnergy(speciesList(prod%species(i))%thermo, T)
        end do

        ! Determine Ka
        reaction_equilCoeff = exp(-dGrxn / 8.314472 / T)
        ! Determine Kc
        reaction_equilCoeff = reaction_equilCoeff * (100000.0 / 8.314472 / T)**(size(prod%species) - size(reac%species))

    end function

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine network_calculateRateCoefficients(net, nIsom, nReac, nProd, &
        Elist, nGrains, Tlist, nT, Plist, nP, method, K)

        integer, intent(in) :: nIsom, nReac, nProd, nGrains, nT, nP
        type(Network), intent(in) :: net
        real(8), dimension(1:nGrains), intent(in) :: Elist
        real(8), dimension(1:nT), intent(in) :: Tlist
        real(8), dimension(1:nP), intent(in) :: Plist
        character(len=*), intent(in) :: method
        real(8), dimension(1:nT,1:nP,1:nIsom+nReac+nProd,1:nIsom+nReac+nProd), intent(out) :: K

        real(8) T, P
        integer i, j, u, v, w

        do u = 1, size(net%Tlist)

            T = Tlist(u)

            ! Calculate the equilibrium (Boltzmann) distributions for each (unimolecular) well
            do i = 1, nIsom
                eqDist(i,:) = isomer_eqDist(net%isomers(i), Elist, nGrains, T)
            end do

            ! Calculate microcanonical rate coefficients at the current conditions
            do w = 1, size(net%reactions)
                call reaction_calculateMicrocanonicalRates(net%reactions(w), Elist, nGrains, T, net%species)
            end do

            do v = 1, size(net%Plist)

                P = Plist(v)

                ! Calculate collision frequencies
                do i = 1, nIsom
                    net%isomers(i)%collFreq = collisionFrequency(T, P, net%bathGas, &
                        net%species(net%isomers(i)%species(1))%general)
                end do

                ! Determine phenomenological rate coefficients
                call network_applyApproximateMethod(net, nIsom, nReac, nProd, Elist, nGrains, T, P, method, K(u,v,:,:))


            end do

        end do

    end subroutine

    subroutine network_applyApproximateMethod(net, nIsom, nReac, nProd, &
        Elist, nGrains, T, P, method, K)

        integer, intent(in) :: nIsom, nReac, nProd, nGrains, nT, nP
        type(Network), intent(in) :: net
        real(8), dimension(1:nGrains), intent(in) :: Elist
        real(8), intent(in) :: T, P
        character(len=*), intent(in) :: method
        real(8), dimension(1:nIsom+nReac+nProd,1:nIsom+nReac+nProd), intent(out) :: K

        real(8), dimension(1:nIsom,1:nGrains) densStates
        real(8), dimension(1:nIsom+nReac+nProd) Eres
        real(8), dimension(1:nIsom,1:nIsom,1:nGrains) Kij
        real(8), dimension(1:nReac+nProd,1:nIsom,1:nGrains) Gnj
        real(8), dimension(1:nIsom,1:nReac,1:nGrains) Fim

        real(8), dimension(1:nIsom) collFreq


        write(1,*) 'Applying method at', T, 'K,', P/1e5, 'bar...'

        dE = Elist[1] - Elist[0]

        ! Density of states per partition function (i.e. normalized density of
        ! states with respect to Boltzmann weighting factor) for each isomer
        do i = 1, nIsom
            densStates(i,:) = net%isomers(i)%densStates * dE / net%isomers(i)%Q
        end do


        ! Active-state energy of each isomer
        do i = 1, nIsom+nReac+nProd
            Eres(i) = isomer_getActiveSpaceEnergy(net%isomers(i), net%reactions)
        end do

        ! Isomerization, dissociation, and association microcanonical rate
        ! coefficients, respectively
        do r = 1, size(net%reactions)
            i = net%reactions(r)%reactant
            j = net%reactions(r)%product
            if (i <= nIsom .and j <= nIsom) then ! Isomerization
                Kij(j,i,:) = net%reactions(r)%kf
                Kij(i,j,:) = net%reactions(r)%kb
            elseif (i <= nIsom .and j > nIsom) then ! Dissociation
                Gnj(j-nIsom,i,:) = net%reactions(r)%kf
                if (j - nIsom <= nReac) Fim(i,j-nIsom,:) = net%reactions(r)%kb
            elseif (i > nIsom .and j <= nIsom) then ! Association
                if (i - nIsom <= nReac) Fim(j,i-nIsom,:) = net%reactions(r)%kf
                Gnj(i-nIsom,j,:) = net%reactions(r)%kb
            end if
        end do


        if (method == 1) then ! Modified strong collision

            ! Modified collision frequency of each isomer
            do i = 1, nIsom
                collFreq(i) = net%isomers(i)%collFreq * &
                    isomer_calculateCollisionEfficiency(net%isomers(i), T, net%reactions, net%dEdown, Elist)
            end do

            ! Apply modified strong collision method
            call estimateratecoefficients_msc(T, P, Elist, collFreq, densStates, Eres, &
                Kij, Fim, Gnj, nIsom, nReac, nProd, nGrains, K, msg)

        elseif (method == 2) then ! Reservoir state

            ! Average energy transferred in a deactivating collision
            dEdown = net%dEdown

            ! Ground-state energy for each isomer
            do i = 1, nIsom
                E0(i) = net%isomers(i)%E0
                collFreq(i) = net%isomers(i)%collFreq
                densStates0 = net%isomers(i)%densStates
                call collisionmatrix(T, P, Elist, collFreq, densStates0, E0[i], dEdown, Mcoll(i,:,:), msg)
            end do

            call estimateratecoefficients_rs(T, P, Elist, Mcoll, densStates, E0, Eres, &
                Kij, Fim, Gnj, dEdown, nIsom, nReac, nProd, nGrains, K, msg)

        end if

        ! Check for validity of phenomenological rate coefficients
        invalidRate = 0
        do i = 1, nIsom+nReac+nProd
            do j = 1, nIsom+nReac
                if (i /= j) then
                    if (K(i,j) <= 0) invalidRate = 1
                end if
            end do
        end do
        if (invalidRate == 1) then
            write (1,fmt='(A)') 'ERROR: One ore more rate coefficients not properly estimated.'
            write (1,*) 'Temperature =', T, 'K, Pressure =', P, 'Pa, Rates ='
            do i = 1, nIsom+nReac+nProd
                write (1,*) K(i,:)
            end do
            stop
        end if

    end subroutine

end module
