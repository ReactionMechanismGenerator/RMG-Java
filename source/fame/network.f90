!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!    network.f90
!
!    Written by Josh Allen (jwallen@mit.edu)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module NetworkModule

    implicit none

    ! dEdown = alpha * (T / T0)^n
    type DeltaEDown
        real(8) :: alpha
        real(8) :: T0
        real(8) :: n
    end type
    
    type GeneralData
        real(8) :: molWt    ! The molecular weight of the bath gas in kg/mol
        real(8) :: sigma    ! The Lennard-Jones sigma parameter of the bath gas in m
        real(8) :: eps      ! The Lennard-Jones epsilon parameter of the bath gas in J
        type(DeltaEDown) :: dEdown   ! Average energy transferred in a deactivating collision in J/mol
    end type

    type SpectralData
        real(8), dimension(:), allocatable :: vibFreq        ! 1D harmonic oscillator frequencies in cm^-1
        real(8), dimension(:), allocatable :: rotFreq        ! 1D rigid rotor frequencies in cm^-1
        real(8), dimension(:), allocatable :: hindFreq       ! 1D hindered rotor frequencies and barrier heights in cm^-1
        real(8), dimension(:), allocatable :: hindBarrier    ! 1D hindered rotor frequencies and barrier heights in cm^-1
        integer :: symmNum                                   ! Symmetry number
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
        character(len=256) :: name                          ! Isomer name
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
        character(len=256) :: equation              ! Chemical reaction equation
        integer :: reac                             ! Isomer wells connected to transition state
        integer :: prod                             ! Isomer wells connected to transition state
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

        type(ThermoData), intent(in) :: thermo
        real(8), intent(in) :: T
        real(8) Cp

        real(8) slope, intercept

        Cp = 0.0

        if (T < 280.0) then
            write (0, fmt='(A)') 'Invalid temperature for heat capacity calculation. (Tmin = 280K)'
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

        type(ThermoData), intent(in) :: thermo
        real(8), intent(in) :: T
        real(8) H

        real(8) slope, intercept

        H = thermo%H298

        if (T < 280.0) then
            write (0, fmt='(A)') 'Invalid temperature for enthalpy calculation. (Tmin = 280K)'
            stop
        end if
        
        if (T < 300.0) then
            H = H + thermo%Cp(1) * (T - 298.0)
        else
            H = H + thermo%Cp(1) * (300.0 - 298.0)
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

        type(ThermoData), intent(in) :: thermo
        real(8), intent(in) :: T
        real(8) S

        real(8) slope, intercept

        S = thermo%S298

        if (T < 280.0) then
            write (0, fmt='(A)') 'Invalid temperature for entropy calculation. (Tmin = 280K)'
            stop
        end if

        if (T < 300.0) then
            S = S + thermo%Cp(1) * log(T / 298.0)
        else
            S = S + thermo%Cp(1) * log(300.0 / 298.0)
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

        type(ThermoData), intent(in) :: thermo
        real(8), intent(in) :: T
        real(8) G

        G = species_getEnthalpy(thermo, T) - T * species_getEntropy(thermo, T)

    end function

    subroutine species_getDensityOfStates(spec, Elist, nGrains, densStates)

        integer, intent(in) :: nGrains
        type(Species), intent(in) :: spec
        real(8), dimension(1:nGrains), intent(in) :: Elist
        real(8), dimension(1:nGrains), intent(out) :: densStates

        real(8), dimension(1:size(spec%spectral%vibFreq)) :: vib
        real(8), dimension(1:size(spec%spectral%rotFreq)) :: rot
        real(8), dimension(1:size(spec%spectral%hindFreq), 1:2) :: hind
        integer linear, symm
        character(len=128) msg

        integer Nvib, Nrot, Nhind, i

        Nvib = size(spec%spectral%vibFreq)
        Nrot = size(spec%spectral%rotFreq)
        Nhind = size(spec%spectral%hindFreq)

        ! Prepare inputs for density of states function
        do i = 1, Nvib
            vib(i) = spec%spectral%vibFreq(i)
        end do
        do i = 1, Nrot
            rot(i) = spec%spectral%rotFreq(i)
        end do
        do i = 1, Nhind
            hind(i,1) = spec%spectral%hindFreq(i)
            hind(i,2) = spec%spectral%hindBarrier(i)
        end do

        linear = 0
        if (size(rot) == 1) linear = 1

        symm = spec%spectral%symmNum

        ! Calculate the density of states
        call densityOfStates(Elist, nGrains, vib, Nvib, rot, Nrot, &
            hind, Nhind, symm, linear, densStates, msg)

    end subroutine

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    function isomer_getCollisionFrequency(T, P, bathGas, spec) result(collFreq)

        ! Provide parameter type-checking
        real(8), intent(in) :: T
        real(8), intent(in) :: P
        type(GeneralData), intent(in) :: bathGas
        type(GeneralData), intent(in) :: spec
        real(8) :: collFreq

        real(8) :: kB, collisionIntegral, sigma, mu, gasConc

        collisionIntegral = 1.16145 / T**0.14874 + 0.52487 / exp(0.77320 * T) + 2.16178 / exp(2.43787 * T) &
            -6.435/10000 * T**0.14874 * sin(18.0323 * T**(-0.76830) - 7.27371)

        kB = 1.3806504e-23

        gasConc = P / kB / T
        mu = 1 / (1/spec%molWt + 1/bathGas%molWt) / 6.022e23
        sigma = 0.5 * (spec%sigma + bathGas%sigma)
        !eps = 0.5 * (spec%eps + bathGas%eps)

        ! Evaluate collision frequency
        collFreq = collisionIntegral * &
            sqrt(8 * kB * T / 3.141592654 / mu) * 3.141592654 * sigma**2 * gasConc

    end function


    function isomer_getCollisionEfficiency(Elist, nGrains, T, isom, isomE0, densStates, reactions, dEdown) result(beta)

        integer, intent(in) :: nGrains
        real(8), dimension(1:nGrains), intent(in) :: Elist
        integer, intent(in) :: isom
        real(8), intent(in) :: isomE0
        real(8), dimension(1:nGrains), intent(in) :: densStates
        real(8), intent(in) :: T
        type(Reaction), dimension(:), intent(in) :: reactions
        real(8), intent(in) :: dEdown

        real(8) beta

        real(8) E0

        ! Determine the "barrier height" as the minimum transition state energy connected to that well
        E0 = isomer_getActiveSpaceEnergy(isom, reactions)

        ! Ensure that the barrier height is sufficiently above the ground state
        ! Otherwise invalid efficiencies are observed
        if (E0 - isomE0 < 100000) E0 = E0 + 100000

        ! Calculate efficiency
        beta = collisionEfficiency(Elist, nGrains, T, E0, dEdown, densStates)

    end function

    function collisionEfficiency(Elist, nGrains, T, E0, dEdown, densStates) result(beta)

        integer, intent(in) :: nGrains
        real(8), dimension(1:nGrains), intent(in) :: Elist
        real(8), intent(in) :: T
        real(8), intent(in) :: E0
        real(8), intent(in) :: dEdown
        real(8), dimension(1:nGrains), intent(in) :: densStates

        real(8) beta

        real(8) dE, Fe, FeNum, FeDen, Delta1, Delta2, DeltaN, Delta, val

        integer r

        dE = Elist(2) - Elist(1)

        FeNum = 0
        FeDen = 0
        Delta1 = 0
        Delta2 = 0
        DeltaN = 0
        Delta = 1
        val = 0.0

        do r = 1, nGrains
            val = densStates(r) * exp(-Elist(r) / 8.314472 / T)
            if (Elist(r) > E0) then
                FeNum = FeNum + val * dE
                if (FeDen == 0) FeDen = val * 8.314472 * T
            end if
        end do

        if (FeDen == 0) then
            beta = 1.0
        else

            Fe = FeNum / FeDen
            if (Fe > 1000000) Fe = 1000000
            do r = 1, nGrains
                val = densStates(r) * exp(-Elist(r) / 8.314472 / T)
                if (Elist(r) < E0) then
                    Delta1 = Delta1 + val * dE
                    Delta2 = Delta2 + val * dE * exp(-(E0 - Elist(r)) / (Fe * 8.314472 * T))
                end if
                DeltaN = DeltaN + val * dE
            end do

            Delta1 = Delta1 / DeltaN
            Delta2 = Delta2 / DeltaN

            Delta = Delta1 - (Fe * 8.314472 * T) / (dEdown + Fe * 8.314472 * T) * Delta2

            beta = (dEdown / (dEdown + Fe * 8.314472 * T))**2 / Delta

        end if

        if (beta < 0 .or. beta > 1) then
            write (1,*), 'Warning: Invalid collision efficiency', beta, 'calculated at', T, 'K.'
            if (beta < 0) beta = 0
            if (beta > 1) beta = 1
        end if

    end function

    subroutine isomer_getEqDist(isom, E, nGrains, T)

        integer, intent(in) :: nGrains
        type(Isomer), intent(inout) :: isom
        real(8), dimension(1:nGrains), intent(in) :: E
        real(8), intent(in) :: T

        real(8) R                   ! Gas constant in J mol^-1 K^-1
        integer s                   ! Dummy index
        real(8) dE

        dE = E(2) - E(1)

        R = 8.314472

        if (allocated(isom%eqDist)) deallocate(isom%eqDist)
        allocate( isom%eqDist(1:size(E)) )
        
        ! Calculate unnormalized eqDist
        do s = 1, size(E)
            isom%eqDist(s) = isom%densStates(s) * exp(-E(s) / (R * T))
        end do

        ! Normalize eqDist
        isom%Q = sum(isom%eqDist) * dE
        
        ! Check that we are going to get a sensible result for eqDist
        ! (i.e. to avoid divide-by-zero error)
        if (isom%Q == 0.) then
            write (1,fmt='(A)') 'ERROR: Partition function is zero, which would give NaN for eqDist. Check density of states.'
            write (*,fmt='(A)') 'ERROR: Partition function is zero, which would give NaN for eqDist. Check density of states.'
            stop
        end if

        isom%eqDist = isom%eqDist / sum(isom%eqDist)

    end subroutine

    function isomer_getActiveSpaceEnergy(isom, reactions) result(Eres)

        integer, intent(in) :: isom
        type(Reaction), dimension(:), intent(in) :: reactions
        real(8) :: Eres

        integer r

        Eres = 1.0e20

        do r = 1, size(reactions)
            if (reactions(r)%reac == isom .or. reactions(r)%prod == isom) then
                if (reactions(r)%E0 < Eres) Eres = reactions(r)%E0
            end if
        end do

    end function

    subroutine isomer_getDensityOfStates(net, isom, Elist, nGrains)

        integer, intent(in) :: nGrains
        type(Network), intent(in) :: net
        type(Isomer), intent(inout) :: isom
        real(8), dimension(1:nGrains), intent(in) :: Elist

        real(8), dimension(1:nGrains) :: densStates0, densStates

        integer index, i, r

        ! Initialize density of states arrays
        if (allocated(isom%densStates)) deallocate(isom%densStates)
        allocate(isom%densStates(1:nGrains))
        do r = 1, nGrains
            isom%densStates(r) = 0.0
            densStates(r) = 0.0
            densStates0(r) = 0.0
        end do

        ! Calculate the density of states for each species, convolving when
        ! multiple species are present
        do i = 1, size(isom%species)
            call species_getDensityOfStates(net%species(isom%species(i)), Elist, nGrains, densStates0)
            call convolve(densStates, densStates0, Elist, nGrains)
        end do

        isom%densStates = densStates

    end subroutine

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine reaction_calculateMicrocanonicalRates(rxn, T, E, nE, isomerList, speciesList, nIsom, nReac, nProd)

        integer, intent(in) :: nE
        integer, intent(in) :: nIsom
        integer, intent(in) :: nReac
        integer, intent(in) :: nProd
        type(Reaction), intent(inout) :: rxn
        real(8), intent(in) :: T
        real(8), dimension(1:nE), intent(in) :: E
        type(Isomer), dimension(:), intent(in) :: isomerList
        type(Species), dimension(:), intent(in) :: speciesList

        type(Isomer) :: reac, prod
        integer nGrains, r
        real(8) dE, Keq, kf0, kf, Keq0, kb
        type(ArrheniusKinetics) :: arrhenius
        real(8), dimension(1:23) :: Tlist0
        
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

        kf0 = 0
        kf = 0
        Keq0 = 0
        Keq = 0
        kb = 0
        
        ! Calculate forward rate coefficient via inverse Laplace transform
        call reaction_kineticsILT(rxn%E0, reac%densStates, rxn%arrhenius, T, E, rxn%kf)

        if (rxn%reac <= nIsom .and. rxn%prod <= nIsom) then

            ! Calculate equilibrium constant
            Keq = reaction_getEquilibriumConstant(reac, prod, T, speciesList)

            ! Calculate backward rate coefficient via detailed balance
            do r = 1, nGrains
                if (prod%densStates(r) /= 0) then
                    rxn%kb(r) = rxn%kf(r) / Keq * &
                        (reac%densStates(r) / reac%Q) / (prod%densStates(r) / prod%Q)
                end if
            end do

            ! Check that forward and reverse rates integrate to give the proper k(E) values
            kf0 = rxn%arrhenius%A * T ** rxn%arrhenius%n * exp(-rxn%arrhenius%Ea / 8.314472 / T)
            kf = sum(rxn%kf * reac%densStates * exp(-E / 8.314472 / T) / reac%Q * dE)
            kb = sum(rxn%kb * prod%densStates * exp(-E / 8.314472 / T) / prod%Q * dE)           
            Keq0 = reaction_getEquilibriumConstant(reac, prod, T, speciesList)
            Keq = kf / kb        

        elseif (rxn%reac <= nIsom .and. rxn%prod > nIsom) then

            ! Calculate equilibrium constant
            Keq = reaction_getEquilibriumConstant(reac, prod, T, speciesList)

            ! Calculate backward rate coefficient via detailed balance
            ! assuming multimolecular isomer is fully thermalized
            do r = 1, nGrains
                rxn%kb(r) = rxn%kf(r) / Keq * &
                    reac%densStates(r) * exp(-E(r) / 8.314472 / T) / reac%Q * dE
            end do

            ! Check that forward and reverse rates integrate to give the proper k(E) values
            kf0 = rxn%arrhenius%A * T ** rxn%arrhenius%n * exp(-rxn%arrhenius%Ea / 8.314472 / T)
            kf = sum(rxn%kf * reac%densStates * exp(-E / 8.314472 / T) / reac%Q * dE)
            kb = sum(rxn%kb)
            Keq0 = reaction_getEquilibriumConstant(reac, prod, T, speciesList)
            Keq = kf / kb   

        elseif (rxn%reac > nIsom .and. rxn%prod <= nIsom) then

            ! Calculate equilibrium constant
            Keq = reaction_getEquilibriumConstant(reac, prod, T, speciesList)

            ! Calculate backward rate coefficient via detailed balance
            do r = 1, nGrains
                if (prod%densStates(r) /= 0) then
                    rxn%kb(r) = rxn%kf(r) / Keq * &
                        (reac%densStates(r) / reac%Q) / (prod%densStates(r) / prod%Q)
                end if
            end do

            ! Include equilibrium distribution in forward direction
            do r = 1, nGrains
                rxn%kf(r) = rxn%kf(r) * reac%densStates(r) * exp(-E(r) / 8.314472 / T) / reac%Q * dE
            end do

            ! Check that forward and reverse rates integrate to give the proper k(E) values
            kf0 = rxn%arrhenius%A * T ** rxn%arrhenius%n * exp(-rxn%arrhenius%Ea / 8.314472 / T)
            kf = sum(rxn%kf)
            kb = sum(rxn%kb * prod%densStates * exp(-E / 8.314472 / T) / prod%Q * dE)           
            Keq0 = reaction_getEquilibriumConstant(reac, prod, T, speciesList)
            Keq = kf / kb  

        end if
        
        ! If the reaction is endothermic and barrierless, it is possible that the
        ! forward k(E) will have a nonzero value at an energy where the product
        ! density of states is zero (but the reactant density of states is not),
        ! which violates detailed balance
        ! To fix, we set the forward k(E) to zero wherever this is true
        ! (This is correct within the accuracy of discretizing the energy grains)
        do r = 1, Ngrains
            if (rxn%kf(r) == 0 .or. rxn%kb(r) == 0) then
                rxn%kf(r) = 0
                rxn%kb(r) = 0
            end if
        end do
        
        if (Keq / Keq0 > 2.0 .or. Keq / Keq0 < 0.5) then
            write (1,*) 'Warning: k(E) values do not satisfy detailed balance!'
            write (1,*) 'Reaction:', rxn%equation
            write (1,*) '    Expected Keq at', T, 'K =', Keq0  
            write (1,*) '      Actual Keq at', T, 'K =', Keq
        end if
        if (kf / kf0 > 10.0 .or. kf / kf0 < 0.1) then
            write (1,*) 'Warning: k(E) values do not match high-pressure-limit!'
            write (1,*) 'Reaction:', rxn%equation
            write (1,*) '    Expected kf at', T, 'K =', kf0 
            write (1,*) '      Actual kf at', T, 'K =', kf
        end if

        !write (*,*) rxn%reac, rxn%prod, rxn%kf(nGrains), rxn%kb(nGrains)

    end subroutine

    subroutine reaction_kineticsILT(E0, rho, kinetics, T, E, k)

        ! Provide parameter type-checking
        real(8), intent(in) :: E0
        real(8), dimension(:), intent(in) :: rho
        type(ArrheniusKinetics), intent(in) :: kinetics
        real(8), intent(in) :: T
        real(8), dimension(:), intent(in) :: E
        real(8), dimension(:), intent(out) :: k

        real(8), dimension(:), allocatable :: phi
        real(8) dE, A, n, Ea
        integer :: r, s

        real(8) r8_gamma

        ! We can't use a negative activation energy for this method, so we
        ! put it in the preexponential if it is encountered.
        A = kinetics%A
        n = kinetics%n
        Ea = kinetics%Ea
        if (Ea < 0) then
            A = A * exp(-Ea / 8.314472 / T)
            Ea = 0.0
        end if
        ! We can't use a negative temperature exponent for this method, so we
        ! put it in the preexponential if it is encountered.
        if (n < 0) then
            A = A * T**n
            n = 0.0
        end if

        dE = E(2) - E(1)
        s = floor(Ea / dE)
        if (s < 0) s = 0

        ! Determine rate coefficients using inverse Laplace transform
        if (n < 0.001) then

            do r = 1, size(E)
                if (s >= r .or. rho(r) == 0) then
                    k(r) = 0
                else
                    k(r) = A * (T ** n) * rho(r - s) / rho(r)
                end if
            end do

        else

            allocate( phi(1:size(E)) )
            ! Evaluate the inverse Laplace transform of the T**n piece, which only
            ! exists for n >= 0
            phi(1) = 0
            do r = 2, size(E)
                phi(r) = (E(r) - E(1))**(n-1) / (8.314472**n * r8_gamma(n))
            end do
            ! Evaluate the convolution
            call convolve(phi, rho, E, size(E))

            ! Apply to determine the microcanonical rate
            do r = s+1, size(E)
                if (E(r) > E0 .and. rho(r) /= 0 .and. phi(r-s) > 0) &
                    k(r) = A * phi(r-s) / rho(r)
            end do
            
            deallocate(phi)

        end if

    end subroutine

    function reaction_getEquilibriumConstant(reac, prod, T, speciesList) result(Keq)

        type(Isomer), intent(in) :: reac
        type(Isomer), intent(in) :: prod
        real(8), intent(in) :: T
        type(Species), dimension(:), intent(in) :: speciesList
        real(8) :: Keq

        real(8) dGrxn
        integer i

        ! Determine free energy of reaction
        dGrxn = 0.0
        do i = 1, size(reac%species)
            dGrxn = dGrxn - species_getFreeEnergy(speciesList(reac%species(i))%thermo, T)
        end do
        do i = 1, size(prod%species)
            dGrxn = dGrxn + species_getFreeEnergy(speciesList(prod%species(i))%thermo, T)
        end do
        
        ! Determine Ka
        Keq = exp(-dGrxn / 8.314472 / T)
        ! Determine Kc
        Keq = Keq * (100000.0 / 8.314472 / T)**(size(prod%species) - size(reac%species))

    end function

    function reaction_fitReverseKinetics(rxn, Tlist, nT, speciesList, isomerList) result(arrhenius)

        integer, intent(in) :: nT
        type(Reaction), intent(in) :: rxn
        real(8), dimension(1:nT), intent(in) :: Tlist
        type(Species), dimension(:), intent(in) :: speciesList
        type(Isomer), dimension(:), intent(in) :: isomerList
        type(ArrheniusKinetics) :: arrhenius

        real(8), dimension(1:nT,1:3) :: A
        real(8), dimension(1:nT) :: b
        
        real(8) Keq
        integer i
        
        real(8), dimension(1:64) :: work
        integer info

        do i = 1, nT
            A(i,1) = 1.0
            A(i,2) = log(Tlist(i))
            A(i,3) = -1.0 / 8.314472 / Tlist(i)

            Keq = reaction_getEquilibriumConstant( &
                isomerList(rxn%reac), isomerList(rxn%prod), Tlist(i), speciesList)
            b(i) = rxn%arrhenius%A * Tlist(i) ** rxn%arrhenius%n * exp(-rxn%arrhenius%Ea / 8.314472 / Tlist(i))
            b(i) = log(b(i) / Keq)
        end do

        call DGELS('N', nT, 3, 1, A, nT, b, nT, work, 64, info )
        if (info /= 0) then
            write (*,*) "Error fitting reverse kinetics!"
            stop
        end if

        arrhenius%A = exp(b(1))
        arrhenius%n = b(2)
        arrhenius%Ea = b(3)

    end function

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine network_shiftToZeroEnergy(net)

        type(Network), intent(inout) :: net

        real(8) Emin
        integer i, j

        Emin = 1.0e20

        do i = 1, size(net%isomers)
            net%isomers(i)%E0 = 0.0
            do j = 1, size(net%isomers(i)%species)
                net%isomers(i)%E0 = net%isomers(i)%E0 + net%species(net%isomers(i)%species(j))%E0
            end do
            if (net%isomers(i)%E0 < Emin) Emin = net%isomers(i)%E0
        end do
        do i = 1, size(net%isomers)
            net%isomers(i)%E0 = net%isomers(i)%E0 - Emin
        end do
        do i = 1, size(net%reactions)
            net%reactions(i)%E0 = net%reactions(i)%E0 - Emin
        end do

    end subroutine

    subroutine network_getEnergyGrains(Emin, Emax, dE0, nGrains0, Elist)

        real(8), intent(in) :: Emin, Emax, dE0
        integer, intent(in) :: nGrains0
        real(8), dimension(:), allocatable, intent(inout) :: Elist

        real(8) dE
        integer nGrains, useGrainSize, r

        if (nGrains0 <= 0 .and. dE0 > 0.0) then
            ! Use the grain size, since the number of grains was not given
            useGrainSize = 1
        elseif (nGrains0 > 0 .and. dE0 <= 0.0) then
            ! Use the number of grains, since the grain size was not given
            useGrainSize = 0
        else
            ! Choose the tighter constraint
            useGrainSize = 0
            dE = (Emax - Emin) / (nGrains0 - 1)
            if (dE > dE0) useGrainSize = 1
        end if

        if (useGrainSize == 0) then
            ! Use the number of grains
            nGrains = nGrains0
            dE = (Emax - Emin) / (nGrains0 - 1)
        else
            ! Use the grain size
            nGrains = int((Emax - Emin) / dE0) + 1
            dE = dE0
        end if

        if (allocated(Elist)) deallocate(Elist)
        allocate(Elist(1:nGrains))
        do r = 1, nGrains
            Elist(r) = Emin + (r - 1) * dE
        end do

    end subroutine

    subroutine network_determineEnergyGrains(net, nIsom, grainSize, nGrains, Tmax, Elist)

        type(Network), intent(in) :: net
        integer, intent(in) :: nIsom
        real(8), intent(in) :: grainSize
        integer, intent(in) :: nGrains
        real(8), intent(in) :: Tmax
        real(8), dimension(:), allocatable, intent(inout) :: Elist

        integer maxIndex
        real(8) Emin, Emax0, Emax
        
        integer i, r

        ! Determine minimum energy and isomer with minimum ground-state energy
        ! Also check reactant and product channels, in case they represent
        ! the minimum energy
        ! if network_shiftToZeroEnergy() has been called, then Emin should be 0
        Emin = net%isomers(1)%E0
        do i = 2, size(net%isomers)
            if (Emin > net%isomers(i)%E0) Emin = net%isomers(i)%E0
        end do
        do i = 2, size(net%reactions)
            if (Emin > net%reactions(i)%E0) Emin = net%reactions(i)%E0
        end do
        
        ! Determine maximum energy and isomer with maximum ground-state energy
        Emax0 = net%isomers(1)%E0
        do i = 2, size(net%isomers)
            if (Emax0 < net%isomers(i)%E0) Emax0 = net%isomers(i)%E0
        end do
        do i = 2, size(net%reactions)
            if (Emax0 < net%reactions(i)%E0) Emax0 = net%reactions(i)%E0
        end do
        
        ! Choose the actual Emax as many kB * T above the maximum energy on the PES
        ! You should check that this is high enough so that the Boltzmann distributions have trailed off to negligible values
        Emax = ceiling(Emax0 + 40 * 8.314472 * Tmax)

        ! Return the chosen energy grains
        call network_getEnergyGrains(Emin, Emax, grainSize, nGrains, Elist)

    end subroutine

    subroutine network_mapDensityOfStates(E0, densStates0, Elist0, nGrains0, densStates, Elist, nGrains)

        real(8), intent(in) :: E0        
        integer, intent(in) :: nGrains0
        real(8), dimension(1:nGrains0), intent(in) :: densStates0
        real(8), dimension(1:nGrains0), intent(in) :: Elist0
        integer, intent(in) :: nGrains
        real(8), dimension(1:nGrains), intent(out) :: densStates
        real(8), dimension(1:nGrains), intent(in) :: Elist

        integer :: r, s

        do r = 1, nGrains
            densStates(r) = 0.
        end do

        do r = 1, nGrains
            if (Elist(r) >= E0) then
                do s = 2, nGrains0
                    if (E0 + Elist0(s) >= Elist(r)) then
                        if (densStates0(s-1) > 0 .and. densStates0(s) > 0) then
                            densStates(r) = densStates0(s) * (densStates0(s-1) / densStates0(s)) ** &
                                ((Elist(r) - E0 - Elist0(s)) / (Elist0(s-1) - Elist0(s)))
                        else
                            densStates(r) = densStates0(s) + (densStates0(s-1) - densStates0(s)) * &
                                (Elist(r) - E0 - Elist0(s)) / (Elist0(s-1) - Elist0(s))
                        end if
                        exit
                    end if
                end do
            end if
        end do

    end subroutine

    subroutine network_calculateRateCoefficients(net, nIsom, nReac, nProd, &
        Elist0, nGrains0, Tlist, nT, Plist, nP, grainSize, numGrains, method, K)

        integer, intent(in) :: nIsom, nReac, nProd, nGrains0, nT, nP
        type(Network), intent(inout) :: net
        real(8), dimension(1:nGrains0), intent(in) :: Elist0
        real(8), dimension(1:nT), intent(in) :: Tlist
        real(8), dimension(1:nP), intent(in) :: Plist
        real(8), intent(in) :: grainSize
        integer, intent(in) :: numGrains
        integer, intent(in) :: method
        real(8), dimension(1:nT,1:nP,1:nIsom+nReac+nProd,1:nIsom+nReac+nProd), intent(out) :: K

        real(8), dimension(:), allocatable :: Elist
        real(8), dimension(1:nIsom+nReac+nProd,1:nGrains0) :: densStates0
        real(8) T, P
        integer i, j, u, v, w, nGrains

        ! Save the original densities of states for each isomer
        ! (We'll be overwriting the ones on each isomer at each temperature)
        do i = 1, nIsom+nReac+nProd
            densStates0(i,:) = net%isomers(i)%densStates
            deallocate(net%isomers(i)%densStates)
        end do

        do u = 1, nT

            T = Tlist(u)

            ! Choose energy grains
            call network_determineEnergyGrains(net, nIsom, grainSize, numGrains, T, Elist)
            nGrains = size(Elist)
            write (1,*) '    Using', nGrains, 'grains of size', Elist(2) - Elist(1), 'J/mol in range', minval(Elist), &
                ' to', maxval(Elist), 'J/mol'

            do i = 1, nIsom+nReac+nProd
                allocate( net%isomers(i)%densStates(1:nGrains) )
                call network_mapDensityOfStates(net%isomers(i)%E0, densStates0(i,:), Elist0, nGrains0, &
                    net%isomers(i)%densStates, Elist, nGrains)
            end do

            ! Calculate the equilibrium (Boltzmann) distributions for each (unimolecular) well
            do i = 1, nIsom+nReac+nProd
                call isomer_getEqDist(net%isomers(i), Elist, nGrains, T)
            end do

            ! Calculate microcanonical rate coefficients at the current conditions
            do w = 1, size(net%reactions)
                allocate(net%reactions(w)%kf(1:nGrains))
                allocate(net%reactions(w)%kb(1:nGrains))
                call reaction_calculateMicrocanonicalRates(net%reactions(w), T, Elist, nGrains, net%isomers, net%species, &
                    nIsom, nReac, nProd)
            end do

            do v = 1, nP

                P = Plist(v)

                ! Calculate collision frequencies
                do i = 1, nIsom
                    net%isomers(i)%collFreq = isomer_getCollisionFrequency(T, P, net%bathGas, &
                        net%species(net%isomers(i)%species(1))%general)
                end do

                ! Determine phenomenological rate coefficients
                call network_applyApproximateMethod(net, nIsom, nReac, nProd, Elist, nGrains, T, P, method, K(u,v,:,:))

            end do

            do i = 1, nIsom+nReac+nProd
                deallocate(net%isomers(i)%densStates)
                deallocate(net%isomers(i)%eqDist)
            end do
            do w = 1, size(net%reactions)
                deallocate(net%reactions(w)%kf)
                deallocate(net%reactions(w)%kb)
            end do

        end do

    end subroutine

    subroutine network_applyApproximateMethod(net, nIsom, nReac, nProd, &
        Elist, nGrains, T, P, method, K)

        integer, intent(in) :: nIsom, nReac, nProd, nGrains
        type(Network), intent(inout) :: net
        real(8), dimension(1:nGrains), intent(in) :: Elist
        real(8), intent(in) :: T, P
        integer, intent(in) :: method
        real(8), dimension(1:nIsom+nReac+nProd,1:nIsom+nReac+nProd), intent(out) :: K

        real(8), dimension(1:nIsom,1:nGrains) :: densStates
        real(8), dimension(1:nIsom+nReac+nProd) :: Eres
        real(8), dimension(1:nIsom) :: E0
        real(8), dimension(1:nIsom,1:nGrains,1:nGrains) :: Mcoll
        real(8), dimension(1:nIsom,1:nIsom,1:nGrains) :: Kij
        real(8), dimension(1:nReac+nProd,1:nIsom,1:nGrains) :: Gnj
        real(8), dimension(1:nIsom,1:nReac,1:nGrains) :: Fim

        real(8), dimension(1:nIsom) :: collFreq

        real(8) dE, dEdown
        character(len=128) :: msg

        integer i, j, r, s

        write(1,*) 'Applying method at', T, 'K,', P/1e5, 'bar...'

        dE = Elist(2) - Elist(1)

        ! Density of states per partition function (i.e. normalized density of
        ! states with respect to Boltzmann weighting factor) for each isomer
        do i = 1, nIsom
            densStates(i,:) = net%isomers(i)%densStates * dE / net%isomers(i)%Q
        end do

        ! Active-state energy of each isomer
        do i = 1, nIsom+nReac+nProd
            Eres(i) = isomer_getActiveSpaceEnergy(i, net%reactions)
        end do

        ! Zero collision matrix    
        do i = 1, nIsom
            do r = 1, Ngrains
                do s = 1, Ngrains
                    Mcoll(i,r,s) = 0.0
                end do
            end do
        end do
        ! Zero rate coefficient matrices
        do r = 1, nGrains
            do i = 1, nIsom
                do j = 1, nIsom
                    Kij(i,j,r) = 0.0
                end do
            end do
            do i = 1, nIsom
                do j = 1, nReac+nProd
                    Gnj(j,i,r) = 0.0
                end do
            end do
            do i = 1, nIsom
                do j = 1, nReac
                    Fim(i,j,r) = 0.0
                end do
            end do
        end do
        
        ! Isomerization, dissociation, and association microcanonical rate
        ! coefficients, respectively
        do r = 1, size(net%reactions)
            i = net%reactions(r)%reac
            j = net%reactions(r)%prod
            if (i <= nIsom .and. j <= nIsom) then ! Isomerization
                Kij(j,i,:) = net%reactions(r)%kf
                Kij(i,j,:) = net%reactions(r)%kb
            elseif (i <= nIsom .and. j > nIsom) then ! Dissociation
                Gnj(j-nIsom,i,:) = net%reactions(r)%kf
                if (j - nIsom <= nReac) Fim(i,j-nIsom,:) = net%reactions(r)%kb
            elseif (i > nIsom .and. j <= nIsom) then ! Association
                if (i - nIsom <= nReac) Fim(j,i-nIsom,:) = net%reactions(r)%kf
                Gnj(i-nIsom,j,:) = net%reactions(r)%kb
            end if
        end do

        ! Average energy transferred in a deactivating collision
        dEdown = net%bathGas%dEdown%alpha * (T / net%bathGas%dEdown%T0) ** net%bathGas%dEdown%n

        if (method == 1) then ! Modified strong collision

            ! Modified collision frequency of each isomer
            do i = 1, nIsom
                collFreq(i) = net%isomers(i)%collFreq * &
                    isomer_getCollisionEfficiency(Elist, nGrains, T, i, &
                    net%isomers(i)%E0, net%isomers(i)%densStates, net%reactions, dEdown)
            end do

            ! Apply modified strong collision method
            call estimateratecoefficients_msc(T, P, Elist, collFreq, densStates, Eres, &
                Kij, Fim, Gnj, nIsom, nReac, nProd, nGrains, K, msg)

        elseif (method == 2) then ! Reservoir state

            ! Ground-state energy for each isomer
            do i = 1, nIsom
                E0(i) = net%isomers(i)%E0
                collFreq(i) = net%isomers(i)%collFreq
                call collisionMatrix(T, P, Elist, collFreq, net%isomers(i)%densStates, E0(i), dEdown, nGrains, Mcoll(i,:,:), msg)
            end do

            call estimateratecoefficients_rs(T, P, Elist, Mcoll, densStates, E0, Eres, &
                Kij, Fim, Gnj, dEdown, nIsom, nReac, nProd, nGrains, K, msg)

        end if

        if (msg /= '') then
            write (*,fmt='(A)') 'ERROR: One or more rate coefficients not properly estimated. See fame.log for details.'
            write (*,fmt='(A)') 'The message returned was:', msg
            write (1,fmt='(A)') 'ERROR: One or more rate coefficients not properly estimated.'
            write (1,fmt='(A)') 'The message returned was:', msg
            write (1,*) 'Temperature =', T, 'K, Pressure =', P, 'Pa, Rates ='
            do i = 1, nIsom+nReac+nProd
                write (1,*) K(i,:)
            end do
            stop
        end if

    end subroutine

end module
