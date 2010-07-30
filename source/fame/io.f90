!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	io.f90
!
!	Written by Josh Allen (jwallen@mit.edu)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module IOModule

    use NetworkModule

contains

    subroutine readInput(net, Tlist, Plist, Tmin, Tmax, Pmin, Pmax, &
      grainSize, numGrains, method, model, modelOptions)


        type(Network), intent(inout) :: net
        real(8), dimension(:), allocatable, intent(inout) :: Tlist, Plist
        real(8), intent(out) :: Tmin, Tmax, Pmin, Pmax
        real(8), intent(out) :: grainSize
        integer, intent(out) :: numGrains
        integer, intent(out) :: method, model
        integer, dimension(:) :: modelOptions

        ! The current line of text from stdin
        character(len=1024) line
        character(len=64) units
        ! String for tokens
        character(len=256) token
        ! Numbers of species, isomers, and reactions
        integer numSpecies, numIsomers, numReactants, numProducts, numReactions
        ! Counters for loops
        integer i

        ! Read method
        line = readMeaningfulLine()
        if (index(line(1:23), 'modifiedstrongcollision') /= 0) then
            method = 1
            write (1,*) "Method set to modified strong collision"
        elseif (index(line(1:14), 'reservoirstate') /= 0) then
            method = 2
            write (1,*) "Method set to reservoir state"
        else
            write (0, fmt='(a)') 'Unable to determine method to use. Should be "ModifiedStrongCollision" or "ReservoirState".'
            stop
        end if

        ! Read temperatures; convert temperature units to K
        line = readMeaningfulLine()
        call processNumberListWithRange(line, units, Tlist, Tmin, Tmax)
        if (index(units(1:1), 'c') /= 0) then
            do i = 1, size(Tlist)
                Tlist(i) = Tlist(i) + 273.15
            end do
            Tmin = Tmin + 273.15
            Tmax = Tmax + 273.15
        elseif (index(units(1:1), 'f') /= 0) then
            do i = 1, size(Tlist)
                Tlist(i) = (Tlist(i) + 459.67) * 5.0 / 9.0
            end do
            Tmin = (Tmin + 459.67) * 5.0 / 9.0
            Tmax = (Tmax + 459.67) * 5.0 / 9.0
        elseif (index(units(1:1), 'r') /= 0) then
            do i = 1, size(Tlist)
                Tlist(i) = Tlist(i) * 5.0 / 9.0
            end do
            Tmin = Tmin * 5.0 / 9.0
            Tmax = Tmax * 5.0 / 9.0
        elseif (index(units(1:1), 'k') == 0) then
            write (0, fmt='(a)') 'Invalid units for temperature. Should be K, C, F, or R.'
            stop
        end if
        write (1,*) "Temperatures =", Tlist

        ! Read pressures; convert pressure units to Pa
        line = readMeaningfulLine()
        call processNumberListWithRange(line, units, Plist, Pmin, Pmax)
        if (index(units(1:3), 'bar') /= 0) then
            do i = 1, size(Plist)
                Plist(i) = Plist(i) * 100000.
            end do
            Pmin = Pmin * 100000.
            Pmax = Pmax * 100000.
        elseif (index(units(1:3), 'atm') /= 0) then
            do i = 1, size(Plist)
                Plist(i) = Plist(i) * 101325.
            end do
            Pmin = Pmin * 101325.
            Pmax = Pmax * 101325.
        elseif (index(units(1:4), 'torr') /= 0) then
            do i = 1, size(Plist)
                Plist(i) = Plist(i) * 101325. / 760.
            end do
            Pmin = Pmin * 101325. / 760.
            Pmax = Pmax * 101325. / 760.
        elseif (index(units(1:2), 'pa') == 0) then
            write (0, fmt='(a)') 'Invalid units for pressure. Should be bar, atm, torr, or Pa.'
            stop
        end if
        write (1,fmt=*) "Pressures =", Plist

        ! Read interpolation model
        line = readMeaningfulLine()
        call getFirstToken(line, token)
        if (index(token, 'none') /= 0) then
            write (1,*) "Model set to none"
            model = 0
            modelOptions(1) = 0
            modelOptions(2) = 0
        elseif (index(token, 'chebyshev') /= 0) then
            write (1,*) "Model set to Chebyshev"
            model = 1
            call getFirstToken(line, token)
            read(token, *), modelOptions(1)
            call getFirstToken(line, token)
            read(token, *), modelOptions(2)
        elseif (index(token, 'pdeparrhenius') /= 0) then
            write (1,*) "Model set to PDepArrhenius"
            model = 2
        else
            write (0, fmt='(a)') 'Invalid interpolation model specification. Should be "None", "Chebyshev", or "PDepArrhenius".'
            stop
        end if

        ! Read number of grains/grain size
        line = readMeaningfulLine()
        call getFirstToken(line, token)
        if (index(token, 'numgrains') /= 0) then
            call getFirstToken(line, token)
            read(token, *), numGrains
            grainSize = 0
            write (1,fmt=*) "Number of grains =", numGrains
        elseif (index(token, 'grainsize') /= 0) then
            call getFirstToken(line, token)
            read(token, *), grainSize
            numGrains = 0
            write (1,fmt=*) "Grain size =", grainSize, "J/mol"
        else
            write (0, fmt='(a)') 'Invalid grain size specification. Should be "NumGrains" or "GrainSize".'
            stop
        end if

        ! Read collisional transfer probability model
        line = readMeaningfulLine()
        call getFirstToken(line, token)
        if (index(token, 'singleexpdown') /= 0) then
            call processQuantity(line, units, net%bathGas%dEdown)
            if (index(units(1:6), 'kj/mol') /= 0) then
                net%bathGas%dEdown = net%bathGas%dEdown * 1000
            elseif (index(units(1:7), 'cal/mol') /= 0) then
                net%bathGas%dEdown = net%bathGas%dEdown * 4.184
            elseif (index(units(1:8), 'kcal/mol') /= 0) then
                net%bathGas%dEdown = net%bathGas%dEdown * 4184
            elseif (index(units(1:5), 'cm^-1') /= 0) then
                net%bathGas%dEdown = net%bathGas%dEdown * 2.9979e10 * 6.626e-34 * 6.022e23
            elseif (index(units(1:5), 'j/mol') == 0) then
                write (0, fmt='(a)') 'Invalid units for single exponential down parameter.'
                stop
            end if
            write (1,fmt=*) "dEdown =", net%bathGas%dEdown, "J/mol"
        else
            write (0, fmt='(a)') 'Invalid collisional transfer probability model specification. Should be "SingleExpDown".'
            stop
        end if

        ! Read bath gas parameters
        call readGasParameters(net%bathGas%molWt, net%bathGas%sigma, net%bathGas%eps)

        ! Read species
        line = readMeaningfulLine()
        call getFirstToken(line, token)
        read(token, *), numSpecies
        write (1,fmt=*) "Found", numSpecies, "species"
        allocate( net%species(1:numSpecies) )
        do i = 1, numSpecies
            call readSpecies(net%species(i))
            write (1,fmt=*) "    Species", i, "is ", net%species(i)%name
        end do

        ! Read isomers
        line = readMeaningfulLine()
        call getFirstToken(line, token)
        read(token, *), numIsomers
        write (1,fmt=*) "Found", numIsomers, "isomers"
        line = readMeaningfulLine()
        call getFirstToken(line, token)
        read(token, *), numReactants
        write (1,fmt=*) "Found", numReactants, "reactants"
        line = readMeaningfulLine()
        call getFirstToken(line, token)
        read(token, *), numProducts
        write (1,fmt=*) "Found", numProducts, "products"
        
        
        allocate( net%isomers(1:numIsomers+numReactants+numProducts) )
        do i = 1, numIsomers+numReactants+numProducts
            call readIsomer(net%isomers(i), net%species)
            if (size(net%isomers(i)%species) == 1) then
                write (1,fmt=*) "    Isomer", i, "is", net%isomers(i)%species(1)
            else
                write (1,fmt=*) "    Isomer", i, "is", net%isomers(i)%species(1), "and", net%isomers(i)%species(2)
            end if
        end do

        ! Read reactions
        line = readMeaningfulLine()
        call getFirstToken(line, token)
        read(token, *), numReactions
        write (1,fmt=*) "Found", numReactions, "path reactions"
        allocate( net%reactions(1:numReactions) )
        do i = 1, numReactions
            call readReaction(net%reactions(i), net%isomers)
            write (1,fmt=*) "    Reaction", i, "is ", net%reactions(i)%equation
        end do

    end subroutine

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! Subroutine: readReaction()
    !
    ! Read data for one reaction from stdin.
    !
    subroutine readReaction(rxn, isomerList)

        type(Reaction), intent(inout) :: rxn
        type(Isomer), dimension(:), intent(in) :: isomerList

        character(len=64) units
        character(len=128) token
        character(len=1024) line
        
        ! Read reaction equation
        line = readMeaningfulLine()
        call getFirstToken(line, token)
        rxn%equation = token

        ! Read reactant and product isomer indices
        line = readMeaningfulLine()
        call getFirstToken(line, token)
        read(token, *) rxn%reac
        call getFirstToken(line, token)
        read(token, *) rxn%prod

        ! Read ground-state energy of transition state
        line = readMeaningfulLine()
        call processQuantity(line, units, rxn%E0)
        if (index(units(1:6), 'kj/mol') /= 0) then
            rxn%E0 = rxn%E0 * 1000
        elseif (index(units(1:7), 'cal/mol') /= 0) then
            rxn%E0 = rxn%E0 * 4.184
        elseif (index(units(1:8), 'kcal/mol') /= 0) then
            rxn%E0 = rxn%E0 * 4184
        elseif (index(units(1:5), 'cm^-1') /= 0) then
            rxn%E0 = rxn%E0 * 2.9979e10 * 6.626e-34 * 6.022e23
        elseif (index(units(1:5), 'j/mol') == 0) then
            write (0, fmt='(a)') 'Invalid units for ground state energy of transition state.'
            stop
        end if

        ! Read kinetics model for high pressure limit
        line = readMeaningfulLine()
        call getFirstToken(line, token)
        if (index(token, 'arrhenius') /= 0) then

            line = readMeaningfulLine()
            call processQuantity(line, units, rxn%arrhenius%A)
            if (index(units(1:4), 's^-1') == 0) then
                write (0, fmt='(a)') 'Invalid units for Arrhenius preexponential.'
                stop
            end if

            line = readMeaningfulLine()
            call processQuantity(line, units, rxn%arrhenius%Ea)
            if (index(units(1:6), 'kj/mol') /= 0) then
                rxn%arrhenius%Ea = rxn%arrhenius%Ea * 1000
            elseif (index(units(1:7), 'cal/mol') /= 0) then
                rxn%arrhenius%Ea = rxn%arrhenius%Ea * 4.184
            elseif (index(units(1:8), 'kcal/mol') /= 0) then
                rxn%arrhenius%Ea = rxn%arrhenius%Ea * 4184
            elseif (index(units(1:5), 'cm^-1') /= 0) then
                rxn%arrhenius%Ea = rxn%arrhenius%Ea * 2.9979e10 * 6.626e-34 * 6.022e23
            elseif (index(units(1:5), 'j/mol') == 0) then
                write (0, fmt='(a)') 'Invalid units for Arrhenius activation energy.'
                stop
            end if

            line = readMeaningfulLine()
            call getFirstToken(line, token)
            read(token, *) rxn%arrhenius%n

        else
            write (0, fmt='(a)') 'Invalid high-pressure-limit kinetics model.'
            stop
        end if

    end subroutine

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! Subroutine: readIsomer()
    !
    ! Read data for one isomer from stdin.
    !
    subroutine readIsomer(isom, speciesList)

        type(Isomer), intent(inout) :: isom
        type(Species), dimension(:), intent(in) :: speciesList

        character(len=128) token
        character(len=1024) line
        integer numSpecies
        integer i, j

        ! Read species list
        line = readMeaningfulLine()
        call getFirstToken(line, token)
        read(token, *), numSpecies

        allocate( isom%species(1:numSpecies) )
        isom%E0 = 0

        do i = 1, numSpecies
            call getFirstToken(line, token)
            do j = 1, size(speciesList)
                if (index(speciesList(j)%name, token) /= 0) then
                    isom%species(i) = j
                    isom%E0 = isom%E0 + speciesList(j)%E0
                end if
            end do
        end do



    end subroutine

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! Subroutine: readSpecies()
    !
    ! Read data for one species from stdin.
    !
    subroutine readSpecies(spec)

        type(Species), intent(inout) :: spec

        character(len=64) units
        character(len=128) token
        character(len=4096) line
        integer i

        ! Read species identifier string
        line = readMeaningfulLine()
        call getFirstToken(line, token)
        spec%name = token

        ! Read ground-state energy
        line = readMeaningfulLine()
        call processQuantity(line, units, spec%E0)
        if (index(units(1:6), 'kj/mol') /= 0) then
            spec%E0 = spec%E0 * 1000
        elseif (index(units(1:7), 'cal/mol') /= 0) then
            spec%E0 = spec%E0 * 4.184
        elseif (index(units(1:8), 'kcal/mol') /= 0) then
            spec%E0 = spec%E0 * 4184
        elseif (index(units(1:5), 'cm^-1') /= 0) then
            spec%E0 = spec%E0 * 2.9979e10 * 6.626e-34 * 6.022e23
        elseif (index(units(1:5), 'j/mol') == 0) then
            write (0, fmt='(a)') 'Invalid units for ground state energy.'
            stop
        end if

        ! Read enthalpy and entropy of formation
        line = readMeaningfulLine()
        call processQuantity(line, units, spec%thermo%H298)
        if (index(units(1:6), 'kj/mol') /= 0) then
            spec%thermo%H298 = spec%thermo%H298 * 1000
        elseif (index(units(1:7), 'cal/mol') /= 0) then
            spec%thermo%H298 = spec%thermo%H298 * 4.184
        elseif (index(units(1:8), 'kcal/mol') /= 0) then
            spec%thermo%H298 = spec%thermo%H298 * 4184
        elseif (index(units(1:5), 'j/mol') == 0) then
            write (0, fmt='(a)') 'Invalid units for enthalpy of formation.'
            stop
        end if

        line = readMeaningfulLine()
        call processQuantity(line, units, spec%thermo%S298)
        if (index(units(1:8), 'kj/mol*f') /= 0 .or. index(units(1:8), 'kj/mol*r') /= 0) then
            spec%thermo%S298 = spec%thermo%S298 * 1000 * 5.0 / 9.0
        elseif (index(units(1:9), 'cal/mol*f') /= 0 .or. index(units(1:9), 'cal/mol*r') /= 0) then
            spec%thermo%S298 = spec%thermo%S298 * 4.184 * 5.0 / 9.0
        elseif (index(units(1:10), 'kcal/mol*f') /= 0 .or. index(units(1:10), 'kcal/mol*r') /= 0) then
            spec%thermo%S298 = spec%thermo%S298 * 4184 * 5.0 / 9.0
        elseif (index(units(1:7), 'j/mol*f') /= 0 .or. index(units(1:7), 'j/mol*r') /= 0) then
            spec%thermo%S298 = spec%thermo%S298 * 1000 * 5.0 / 9.0
        elseif (index(units(1:8), 'kj/mol*k') /= 0 .or. index(units(1:8), 'kj/mol*c') /= 0) then
            spec%thermo%S298 = spec%thermo%S298 * 1000
        elseif (index(units(1:9), 'cal/mol*k') /= 0 .or. index(units(1:9), 'cal/mol*c') /= 0) then
            spec%thermo%S298 = spec%thermo%S298 * 4.184
        elseif (index(units(1:10), 'kcal/mol*k') /= 0 .or. index(units(1:10), 'kcal/mol*c') /= 0) then
            spec%thermo%S298 = spec%thermo%S298 * 4184
        elseif (index(units(1:7), 'j/mol*k') == 0 .or. index(units(1:7), 'j/mol*c') /= 0) then
            write (0, fmt='(a)') 'Invalid units for entropy of formation.'
            stop
        end if

        ! Read list of heat capacities
        line = readMeaningfulLine()
        call processNumberList(line, units, spec%thermo%Cp)
        if (index(units(1:8), 'kj/mol*f') /= 0 .or. index(units(1:8), 'kj/mol*r') /= 0) then
            do i = 1, size(spec%thermo%Cp)
                spec%thermo%Cp(i) = spec%thermo%Cp(i) * 1000 * 5.0 / 9.0
            end do
        elseif (index(units(1:9), 'cal/mol*f') /= 0 .or. index(units(1:9), 'cal/mol*r') /= 0) then
            do i = 1, size(spec%thermo%Cp)
                spec%thermo%Cp(i) = spec%thermo%Cp(i) * 4.184 * 5.0 / 9.0
            end do
        elseif (index(units(1:10), 'kcal/mol*f') /= 0 .or. index(units(1:10), 'kcal/mol*r') /= 0) then
            do i = 1, size(spec%thermo%Cp)
                spec%thermo%Cp(i) = spec%thermo%Cp(i) * 4184 * 5.0 / 9.0
            end do
        elseif (index(units(1:7), 'j/mol*f') /= 0 .or. index(units(1:7), 'j/mol*r') /= 0) then
            do i = 1, size(spec%thermo%Cp)
                spec%thermo%Cp(i) = spec%thermo%Cp(i) * 1000 * 5.0 / 9.0
            end do
        elseif (index(units(1:8), 'kj/mol*k') /= 0 .or. index(units(1:8), 'kj/mol*c') /= 0) then
            do i = 1, size(spec%thermo%Cp)
                spec%thermo%Cp(i) = spec%thermo%Cp(i) * 1000
            end do
        elseif (index(units(1:9), 'cal/mol*k') /= 0 .or. index(units(1:9), 'cal/mol*c') /= 0) then
            do i = 1, size(spec%thermo%Cp)
                spec%thermo%Cp(i) = spec%thermo%Cp(i) * 4.184
            end do
        elseif (index(units(1:10), 'kcal/mol*k') /= 0 .or. index(units(1:10), 'kcal/mol*c') /= 0) then
            do i = 1, size(spec%thermo%Cp)
                spec%thermo%Cp(i) = spec%thermo%Cp(i) * 4184
            end do
        elseif (index(units(1:7), 'j/mol*k') == 0 .or. index(units(1:7), 'j/mol*c') /= 0) then
            write (0, fmt='(a)') 'Invalid units for heat capacity.'
            stop
        end if

        ! Read species gas parameters
        call readGasParameters(spec%general%molWt, spec%general%sigma, spec%general%eps)

        ! Read list of harmonic oscillator frequencies
        line = readMeaningfulLine()
        call processNumberList(line, units, spec%spectral%vibFreq)
        if (index(units(1:2), 'Hz') /= 0) then
            do i = 1, size(spec%spectral%vibFreq)
                spec%spectral%vibFreq(i) = spec%spectral%vibFreq(i) / 2.9979e10
            end do
        elseif (index(units(1:5), 'cm^-1') == 0) then
            write (0, fmt='(a)') 'Invalid units for harmonic oscillator frequencies.'
            stop
        end if

        ! Read list of rigid rotor frequencies
        line = readMeaningfulLine()
        call processNumberList(line, units, spec%spectral%rotFreq)
        if (index(units(1:2), 'Hz') /= 0) then
            do i = 1, size(spec%spectral%rotFreq)
                spec%spectral%rotFreq(i) = spec%spectral%rotFreq(i) / 2.9979e10
            end do
        elseif (index(units(1:5), 'cm^-1') == 0) then
            write (0, fmt='(a)') 'Invalid units for rigid rotor frequencies.'
            stop
        end if

        ! Read list of hindered rotor frequencies
        line = readMeaningfulLine()


        call processNumberList(line, units, spec%spectral%hindFreq)
        if (index(units(1:2), 'Hz') /= 0) then
            do i = 1, size(spec%spectral%hindFreq)
                spec%spectral%hindFreq(i) = spec%spectral%hindFreq(i) / 2.9979e10
            end do
        elseif (index(units(1:5), 'cm^-1') == 0) then
            write (0, fmt='(a)') 'Invalid units for hindered rotor frequencies.'
            stop
        end if

        ! Read list of hindered rotor barriers
        line = readMeaningfulLine()
        call processNumberList(line, units, spec%spectral%hindBarrier)
        if (index(units(1:6), 'kj/mol') /= 0) then
            do i = 1, size(spec%spectral%hindBarrier)
                spec%spectral%hindBarrier(i) = spec%spectral%hindBarrier(i) * 1000
            end do
        elseif (index(units(1:7), 'cal/mol') /= 0) then
            do i = 1, size(spec%spectral%hindBarrier)
                spec%spectral%hindBarrier(i) = spec%spectral%hindBarrier(i) * 4.184
            end do
        elseif (index(units(1:8), 'kcal/mol') /= 0) then
            do i = 1, size(spec%spectral%hindBarrier)
                spec%spectral%hindBarrier(i) = spec%spectral%hindBarrier(i) * 4184
            end do
        elseif (index(units(1:5), 'cm^-1') /= 0) then
            do i = 1, size(spec%spectral%hindBarrier)
                spec%spectral%hindBarrier(i) = spec%spectral%hindBarrier(i) * 2.9979e10 * 6.626e-34 * 6.022e23
            end do
        elseif (index(units(1:5), 'j/mol') /= 0) then
            write (0, fmt='(a)') 'Invalid units for hindered rotor barriers.'
            stop
        end if
        do i = 1, size(spec%spectral%hindBarrier)
            spec%spectral%hindBarrier(i) = spec%spectral%hindBarrier(i) / (2.9979e10 * 6.626e-34 * 6.022e23)
        end do

        ! Read symmetry number
        line = readMeaningfulLine()
        call getFirstToken(line, token)
        read(token, *) spec%spectral%symmNum


    end subroutine

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! Subroutine: readGasParameters()
    !
    ! Read a line from stdin containing the gas parameters: molecular weight,
    ! LJ sigma, and LJ epsilon.
    !
    subroutine readGasParameters(molWt, sigma, eps)

        real(8), intent(out) :: molWt
        real(8), intent(out) :: sigma
        real(8), intent(out) :: eps

        character(len=64) units
        character(len=1024) line

        line = readMeaningfulLine()
        call processQuantity(line, units, molWt)
        if (index(units(1:1), 'u') /= 0 .or. index(units(1:5), 'g/mol') /= 0) then
            molWt = molWt / 1000.
        else
            write (0, fmt='(a)') 'Invalid units for molecular weight. Should be u or g/mol.'
            stop
        end if

        line = readMeaningfulLine()
        call processQuantity(line, units, sigma)
        if (index(units(1:1), 'a') /= 0) then
            sigma = sigma * 1.0e-10
        elseif (index(units(1:1), 'm') == 0) then
            write (0, fmt='(a)') 'Invalid units for Lennard-Jones sigma parameter. Should be m or A.'
            stop
        end if

        line = readMeaningfulLine()
        call processQuantity(line, units, eps)
        if (index(units(1:1), 'k') /= 0) then
            eps = eps * 1.381e-23
        elseif (index(units(1:1), 'j') == 0) then
            write (0, fmt='(a)') 'Invalid units for Lennard-Jones epsilon parameter. Should be J or K.'
            stop
        end if

    end subroutine

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! Function: readMeaningfulLine()
    !
    ! Reads a meaningful line from stdin.
    !
    function readMeaningfulLine() result(line)

        character(len=1024) line

        ! The i/o status flag
        integer ios, found

        ios = 0
        found = 0
        do while (ios == 0 .and. found == 0)

            ! Read one line from the file
            read (*, fmt='(a1024)', iostat=ios), line

            ! Print the input line (as a comment) to the output, for debugging
            ! WRITE(*,fmt='(A,A)') '#IN: ', trim(line)

            ! Skip if comment line
            if (index(line(1:1), '#') /= 0) cycle

            ! Skip if blank line
            if (len(trim(line)) == 0) cycle

            ! Otherwise return it
            found = 1

        end do

        line = toLowercase(adjustl(line))

    end function

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! Function: toLowercase()
    !
    ! Converts a string to all lowercase.
    !
    function toLowercase(string)

        character(len=*), intent(in) :: string
        character(len=len(string)) :: toLowercase

        character(len=26) lower, upper
        integer i, n

        lower = 'abcdefghijklmnopqrstuvwxyz'
        upper = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

        toLowercase = string
        do i = 1, len(string)
            n = index(upper, string(i:i))
            if (n /= 0) toLowercase(i:i) = lower(n:n)
        end do

    end function

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! Subroutine: processNumberList()
    !
    ! Reads a list of numbers from a string with the format 'count (units) value1 value2 ...'
    !
    subroutine processNumberList(string, units, values)

        character(len=*), intent(inout) :: string
        character(len=64), intent(inout) :: units
        real(8), dimension(:), allocatable, intent(inout) :: values

        character(len=256) token
        integer numValues, n
        
        call getFirstToken(string, token)
        read(token, *), numValues

        call getFirstToken(string, token)
        units = token

        allocate(values(1:numValues))
        do n = 1, numValues
            string = readMeaningfulLine()
            call getFirstToken(string, token)
            read(token, *), values(n)
        end do

    end subroutine

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! Subroutine: processNumberListWithRange()
    !
    ! Reads a list of numbers from a string with the format 'count (units) value1 value2 ...'
    !
    subroutine processNumberListWithRange(string, units, values, Nmin, Nmax)

        character(len=*), intent(inout) :: string
        character(len=64), intent(inout) :: units
        real(8), dimension(:), allocatable, intent(inout) :: values
        real(8), intent(out) :: Nmin
        real(8), intent(out) :: Nmax
        
        character(len=256) token
        integer numValues, n
        
        call getFirstToken(string, token)
        read(token, *), numValues

        call getFirstToken(string, token)
        units = token

        call getFirstToken(string, token)
        read(token, *), Nmin

        call getFirstToken(string, token)
        read(token, *), Nmax

        allocate(values(1:numValues))
        do n = 1, numValues
            string = readMeaningfulLine()
            call getFirstToken(string, token)
            read(token, *), values(n)
        end do

    end subroutine

    subroutine processQuantity(string, units, value)

        character(len=*), intent(inout) :: string
        character(len=64), intent(inout) :: units
        real(8), intent(inout) :: value

        character(len=256) token

        call getFirstToken(string, token)
        units = token

        call getFirstToken(string, token)
        read(token, *), value

    end subroutine

    subroutine getFirstToken(string, token)

        character(len=*), intent(inout) :: string
        character(len=*), intent(inout) :: token

        integer i

        string = adjustl(string)
        i = getFirstWhitespaceIndex(string)
        token = string(1:i-1)
        string = adjustl(string(i:))
        !write (*,fmt='(a32,a64)') token, string

    end subroutine

    function getFirstWhitespaceIndex(string)

        character(len=*), intent(in) :: string
        integer getFirstWhitespaceIndex

        getFirstWhitespaceIndex = index(string, ' ')
        if (index(string, '\t') > 0 .and. (getFirstWhitespaceIndex == 0 .or. getFirstWhitespaceIndex > index(string, '\t'))) then
            getFirstWhitespaceIndex = index(string, '\t')
        end if
        if (index(string, '\r') > 0 .and. (getFirstWhitespaceIndex == 0 .or. getFirstWhitespaceIndex > index(string, '\r'))) then
            getFirstWhitespaceIndex = index(string, '\r')
        end if
        if (index(string, '\n') > 0 .and. (getFirstWhitespaceIndex == 0 .or. getFirstWhitespaceIndex > index(string, '\n'))) then
            getFirstWhitespaceIndex = index(string, '\n')
        end if

    end function

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine writeOutputIntro(net, nIsom, nReac, nProd, Tlist, nT, Plist, nP, Elist, nGrains, &
        method, K, model, modelOptions, chebyshev, pDepArrhenius)

        integer :: nIsom, nReac, nProd, nGrains, nT, nP
        type(Network), intent(in) :: net
        real(8), dimension(1:nT), intent(in) :: Tlist
        real(8), dimension(1:nP), intent(in) :: Plist
        real(8), dimension(1:nGrains), intent(in) :: Elist
        integer, intent(in) :: method, model
        real(8), dimension(1:nT,1:nP,1:nIsom+nReac+nProd,1:nIsom+nReac+nProd), intent(in) :: K
        integer, dimension(1:10), intent(in) :: modelOptions
        real(8), dimension(:,:,:,:), intent(in) :: chebyshev
        type(ArrheniusKinetics), dimension(1:nP,1:nIsom+nReac+nProd,1:nIsom+nReac+nProd), intent(in) :: pDepArrhenius

        ! Header
        write (*, fmt='(A)'), '################################################################################'
        write (*, fmt='(A)'), '#'
        write (*, fmt='(A)'), '#   FAME output'
        write (*, fmt='(A)'), '#'
        write (*, fmt='(A)'), '################################################################################'
        write (*, fmt='(A)'), ''

        ! Method
        write (*, fmt='(A)'), '# The method used to calculate the phenomenological rate coefficients k(T, P)'
        if (method == 1) then
            write (*, *), 'ModifiedStrongCollision'
        elseif (method == 2) then
            write (*, *), 'ReservoirState'
        end if
        write (*, fmt='(A)'), ''

        ! Temperatures
        write (*, fmt='(A)'), '# The temperatures at which the k(T, P) were estimated'
        write (*, *), nT, 'K', Tlist
        write (*, fmt='(A)'), ''
        ! Pressures
        write (*, fmt='(A)'), '# The pressures at which the k(T, P) were estimated'
        write (*, *), nP, 'Pa', Plist
        write (*, fmt='(A)'), ''

        ! Interpolation model
        write (*, fmt='(A)'), '# The interpolation model used to fit the k(T, P)'
        if (model == 1) then
            write (*, *), 'Chebyshev', modelOptions(1), modelOptions(2)
        elseif (model == 2) then
            write (*, *), 'PDepArrhenius'
        else
            write (*, *), 'None'
        end if
        write (*, fmt='(A)'), ''

        ! Number of species
        write (*, fmt='(A)'), '# The number of species in the network'
        write (*, *), size(net%species)
        write (*, fmt='(A)'), ''
        ! Number of isomers
        write (*, fmt='(A)'), '# The number of isomers in the network'
        write (*, *), nIsom
        write (*, fmt='(A)'), '# The number of reactant channels in the network'
        write (*, *), nReac
        write (*, fmt='(A)'), '# The number of product channels in the network'
        write (*, *), nProd
        write (*, fmt='(A)'), ''
        ! Number of reactions
        write (*, fmt='(A)'), '# The number of reactions in the network'
        write (*, *), size(net%reactions)
        write (*, fmt='(A)'), ''
        ! Number of k(T, P) values
        write (*, fmt='(A)'), '# The number of phenomenological rate coefficients in the network'
        write (*, *), ((nIsom+nReac)*(nIsom+nReac+nProd-1))
        write (*, fmt='(A)'), ''

    end subroutine

end module