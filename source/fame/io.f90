!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	io.f90
!
!	Written by Josh Allen (jwallen@mit.edu)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Subroutine: loadNetwork()
!
! Loads information about the activated reaction network from stdin.
!
subroutine loadNetwork(net)

    type(Network), intent(inout)		:: net

    ! The i/o status flag
    integer ios
    ! The current line of text from stdin
    character(len=1024) line
    character(len=64) units
    ! String for tokens
    character(len=256) token
    ! Numbers of species, isomers, and reactions
    integer numSpecies, numIsomers, numReactions
    ! Counters for loops
    integer i

    ! Read method
    line = readMeaningfulLine()
    if (index(line(1:23), 'modifiedstrongcollision') /= 0) then
        net%method = 1
    elseif (index(line(1:14), 'reservoirstate') /= 0) then
        net%method = 2
    else
        write (0, fmt='(a)') 'Unable to determine method to use. Should be "ModifiedStrongCollision" or "ReservoirState".'
        stop
    end if

    ! Read temperatures; convert temperature units to K
    line = readMeaningfulLine()
    call processNumberList(line, units, net%Tlist)
    if (index(units(1:1), 'c') /= 0) then
        do i = 1, size(net%Tlist)
            net%Tlist(i) = net%Tlist(i) + 273.15
        end do
    elseif (index(units(1:1), 'f') /= 0) then
        do i = 1, size(net%Tlist)
            net%Tlist(i) = (net%Tlist(i) + 459.67) * 5.0 / 9.0
        end do
    elseif (index(units(1:1), 'r') /= 0) then
        do i = 1, size(net%Tlist)
            net%Tlist(i) = net%Tlist(i) * 5.0 / 9.0
        end do
    elseif (index(units(1:1), 'k') == 0) then
        write (0, fmt='(a)') 'Invalid units for temperature. Should be K, C, F, or R.'
        stop
    end if

    ! Read pressures; convert pressure units to Pa
    line = readMeaningfulLine()
    call processNumberList(line, units, net%Plist)
    if (index(units(1:3), 'bar') /= 0) then
        do i = 1, size(net%Plist)
            net%Plist(i) = net%Plist(i) * 100000.
        end do
    elseif (index(units(1:3), 'atm') /= 0) then
        do i = 1, size(net%Plist)
            net%Plist(i) = net%Plist(i) * 101325.
        end do
    elseif (index(units(1:4), 'torr') /= 0) then
        do i = 1, size(net%Plist)
            net%Plist(i) = net%Plist(i) * 101325. / 760.
        end do
    elseif (index(units(1:2), 'pa') == 0) then
        write (0, fmt='(a)') 'Invalid units for pressure. Should be bar, atm, torr, or Pa.'
        stop
    end if

    ! Read interpolation model
    line = readMeaningfulLine()
    call getFirstToken(line, token)
    if (index(token, 'none') /= 0) then
        net%model = 0
        net%numChebT = 0
        net%numChebP = 0
    elseif (index(token, 'chebyshev') /= 0) then
        net%model = 1
        call getFirstToken(line, token)
        read(token, *), net%numChebT
        call getFirstToken(line, token)
        read(token, *), net%numChebP
    elseif (index(token, 'pdeparrhenius') /= 0) then
        net%model = 2
        net%numChebT = 0
        net%numChebP = 0
    else
        write (0, fmt='(a)') 'Invalid interpolation model specification. Should be "None", "Chebyshev", or "PDepArrhenius".'
        stop
    end if

    ! Read number of grains/grain size
    line = readMeaningfulLine()
    call getFirstToken(line, token)
    if (index(token, 'numgrains') /= 0) then
        call getFirstToken(line, token)
        read(token, *), net%numGrains
        net%grainSize = 0
    elseif (index(token, 'grainsize') /= 0) then
        call getFirstToken(line, token)
        read(token, *), net%grainSize
        net%numGrains = 0
    else
        write (0, fmt='(a)') 'Invalid grain size specification. Should be "NumGrains" or "GrainSize".'
        stop
    end if

    ! Read collisional transfer probability model
    line = readMeaningfulLine()
    call getFirstToken(line, token)
    if (index(token, 'singleexpdown') /= 0) then
        net%collisionModel = 1
        allocate(net%collisionParameters(1:1))
        call processQuantity(line, units, net%collisionParameters(1))
        if (index(units(1:6), 'kj/mol') /= 0) then
            net%collisionParameters(1) = net%collisionParameters(1) * 1000
        elseif (index(units(1:7), 'cal/mol') /= 0) then
            net%collisionParameters(1) = net%collisionParameters(1) * 4.184
        elseif (index(units(1:8), 'kcal/mol') /= 0) then
            net%collisionParameters(1) = net%collisionParameters(1) * 4184
        elseif (index(units(1:5), 'cm^-1') /= 0) then
            net%collisionParameters(1) = net%collisionParameters(1) * 2.9979e10 * 6.626e-34 * 6.022e23
        elseif (index(units(1:5), 'j/mol') == 0) then
            write (0, fmt='(a)') 'Invalid units for single exponential down parameter.'
            stop
        end if
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
    allocate( net%species(1:numSpecies) )
    do i = 1, numSpecies
        call readSpecies(net%species(i))
    end do

    ! Read isomers
    line = readMeaningfulLine()
    call getFirstToken(line, token)
    read(token, *), numIsomers
    allocate( net%isomers(1:numIsomers) )
    do i = 1, numIsomers
        call readIsomer(net%isomers(i), net%species)
    end do

    ! Read reactions
    line = readMeaningfulLine()
    call getFirstToken(line, token)
    read(token, *), numReactions
    allocate( net%reactions(1:numReactions) )
    do i = 1, numReactions
        call readReaction(net%reactions(i), net%isomers)
    end do

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Subroutine: readReaction()
!
! Read data for one reaction from stdin.
!
subroutine readReaction(rxn, isomerList)

    type(Reaction), intent(inout)				:: rxn
    type(Isomer), dimension(:), intent(in)		:: isomerList

    character(len=64) units
    character(len=128) token
    character(len=1024) line
    integer numIsomers
    integer i, j

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

    type(Isomer), intent(inout)				:: isom
    type(Species), dimension(:), intent(in)	:: speciesList

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

    type(Species), intent(inout)	:: spec

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

    real(8), intent(out)	:: molWt
    real(8), intent(out)	:: sigma
    real(8), intent(out)	:: eps

    character(len=64) units
    character(len=128) token
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
function readMeaningfulLine()

    character(len=1024) readMeaningfulLine

    ! The i/o status flag
    integer ios, found

    ios = 0
    found = 0
    do while (ios == 0 .and. found == 0)

        ! Read one line from the file
        read (*, fmt='(a1024)', iostat=ios), readMeaningfulLine

        ! Print the input line (as a comment) to the output, for debugging
        ! WRITE(*,fmt='(A,A)') '#IN: ', trim(readMeaningfulLine)

        ! Skip if comment line
        if (index(readMeaningfulLine(1:1), '#') /= 0) cycle

        ! Skip if blank line
        if (len(trim(readMeaningfulLine)) == 0) cycle

        ! Otherwise return it
        found = 1

    end do

    readMeaningfulLine = toLowercase(adjustl(readMeaningfulLine))

end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Function: toLowercase()
!
! Converts a string to all lowercase.
!
function toLowercase(string)

    character(len=*), intent(in) 	:: string
    character(len=len(string)) 		:: toLowercase

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

    character(len=*), intent(inout)			:: string
    character(len=64), intent(inout)		:: units
    real(8), dimension(:), allocatable, intent(inout)	:: values

    character(len=256) token
    integer numValues, i, n
    real(8) value

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

subroutine processQuantity(string, units, value)

    character(len=*), intent(inout)			:: string
    character(len=64), intent(inout)		:: units
    real(8), intent(inout)					:: value

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
		
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Subroutine: saveNetwork()
!
! Saves information about the activated reaction network to stdout.
!
subroutine saveNetwork(net, K, chebyshev, pDepArrhenius)

    type(Network), intent(in)		:: net
    real(8), dimension(:,:,:,:), intent(in)					:: K
    real(8), dimension(:,:,:,:), intent(in)					:: chebyshev
    type(ArrheniusKinetics), dimension(:,:,:), intent(in) 	:: 	pDepArrhenius

    character(len=64) fmtStr
    integer i, j, t, p

    ! Header
    write (*, fmt='(A)'), '################################################################################'
    write (*, fmt='(A)'), '#'
    write (*, fmt='(A)'), '#   FAME output'
    write (*, fmt='(A)'), '#'
    write (*, fmt='(A)'), '################################################################################'
    write (*, fmt='(A)'), ''

    ! Method
    write (*, fmt='(A)'), '# The method used to calculate the phenomenological rate coefficients k(T, P)'
    if (net%method == 1) then
        write (*, *), 'ModifiedStrongCollision'
    elseif (net%method == 2) then
        write (*, *), 'ReservoirState'
    end if
    write (*, fmt='(A)'), ''

    ! Temperatures
    write (*, fmt='(A)'), '# The temperatures at which the k(T, P) were estimated'
    write (*, *), size(net%Tlist), 'K', net%Tlist
    write (*, fmt='(A)'), ''
    ! Pressures
    write (*, fmt='(A)'), '# The pressures at which the k(T, P) were estimated'
    write (*, *), size(net%Plist), 'Pa', net%Plist
    write (*, fmt='(A)'), ''

    ! Interpolation model
    write (*, fmt='(A)'), '# The interpolation model used to fit the k(T, P)'
    if (net%model == 1) then
        write (*, *), 'Chebyshev', net%numChebT, net%numChebP
    elseif (net%model == 2) then
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
    write (*, *), size(net%isomers)
    write (*, fmt='(A)'), ''
    ! Number of reactions
    write (*, fmt='(A)'), '# The number of reactions in the network'
    write (*, *), size(net%reactions)
    write (*, fmt='(A)'), ''
    ! Number of k(T, P) values
    write (*, fmt='(A)'), '# The number of phenomenological rate coefficients in the network'
    write (*, *), (size(net%isomers) * (size(net%isomers) - 1))
    write (*, fmt='(A)'), ''

    ! The phenomenological rate coefficients
    do j = 1, size(net%isomers)
        do i = 1, size(net%isomers)
            if (i /= j) then
                write (*, fmt='(A)'), '# The reactant and product isomers'
                write (*,*), j, i
                write (*, fmt='(A)'), '# Table of phenomenological rate coefficients'
                write (fmtStr,*), '(A8,', size(net%Plist), 'ES14.2E2)'
                write (*, fmtStr), 'T \ P', net%Plist
                write (fmtStr,*), '(F8.1,', size(net%Plist), 'ES14.4E3)'
                do t = 1, size(net%Tlist)
                    write (*, fmt=fmtStr), net%Tlist(t), K(t,:,i,j)
                end do
                if (net%model == 1) then
                    write (*, fmt='(A)'), '# The fitted Chebyshev polynomial model'
                    write (fmtStr,*), '(', net%numChebT, 'ES14.4E3)'
                    do t = 1, net%numChebT
                        write (*,fmt=fmtStr), chebyshev(t,:,i,j)
                    end do
                elseif (net%model == 2) then
                    write (*, fmt='(A)'), '# The fitted pressure-dependent Arrhenius model'
                    do p = 1, size(net%Plist)
                        write (*, fmt='(ES8.2E2,ES14.4E3,F14.4,F10.4)'), net%Plist(p), &
                            pDepArrhenius(p,i,j)%A, pDepArrhenius(p,i,j)%Ea, pDepArrhenius(p,i,j)%n
                    end do
                end if
                write (*, fmt='(A)'), ''
            end if
        end do
    end do
    !write (fmtStr,*), '(I4,A,', size(net%Tlist), 'F8.1)'
    !write (*, fmt=fmtStr), size(net%Tlist), ' K', net%Tlist



end subroutine
