!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   fame.f90
!
!   Copyright (c) 2008-2009 by Josh Allen.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program fame

    use NetworkModule
    use IOModule
    use ModelModule

    implicit none

    type(Network) net
    real(8), dimension(:), allocatable :: Tlist, Plist, Elist
    real(8) :: grainSize, Tmin, Tmax, Pmin, Pmax
    integer :: numGrains
    integer :: method, model
    integer, dimension(1:10) :: modelOptions
    real(8), dimension(:,:,:,:), allocatable :: K
    real(8), dimension(:,:,:,:), allocatable :: chebyshevCoeffs
    type(ArrheniusKinetics), dimension(:,:,:), allocatable :: pDepArrhenius
    integer :: nIsom, nReac, nProd, nGrains, nT, nP, done

    integer i, j, t, p, reac, prod
    character(len=64) fmtStr

    integer invalidRate

    ! Use unit 1 for logging; file will be called fame.log or fort.1
    open(1, file='fame.log')

    ! Log header

    ! Read network information (from stdin)
    write (unit=1,fmt='(A)') 'Reading network information...'
    call readInput(net, Tlist, Plist, Tmin, Tmax, Pmin, Pmax, &
        grainSize, numGrains, method, model, modelOptions)
    nIsom = 0
    nReac = 0
    nProd = 0
    done = 0
    do i = 1, size(net%isomers)
        if (size(net%isomers(i)%species) > 1) then
            done = 1
            nReac = nReac + 1
        elseif (done == 0) then
            nIsom = nIsom + 1
        else
            nProd = nProd + 1
        end if
    end do
    nT = size(Tlist)
    nP = size(Plist)

    ! Shift network such that lowest-energy isomer has a ground state of 0.0
    write (1,fmt='(A)') 'Zeroing lowest energy isomer...'
    call network_shiftToZeroEnergy(net)

    ! Determine energy grains
    write (1,fmt='(A)') 'Determining energy grains...'
    call network_determineEnergyGrains(net, nIsom, grainSize, numGrains, maxval(Tlist), Elist)
    nGrains = size(Elist)

    ! Calculate density of states for all isomers in network
    write (1,fmt='(A)') 'Calculating densities of states...'
    do i = 1, nIsom
        call isomer_getDensityOfStates(net, net%isomers(i), Elist, nGrains)
    end do

    ! Determine phenomenological rate coefficients
    write (1,fmt='(A)') 'Calculating phenomenological rate coefficients...'
    allocate( K( 1:nT, 1:nP, 1:nIsom+nReac+nProd, 1:nIsom+nReac+nProd) )
    call network_calculateRateCoefficients(net, nIsom, nReac, nProd, &
        Elist, nGrains, Tlist, nT, Plist, nP, method, K)

    ! Adjust units of k(T,P) values to combinations of cm^3, mol, and s
    ! Prior to this point FAME has been thinking in terms of m^3, mol, and s,
    ! but RMG thinks in terms of cm^3, mol, and s
    ! This particularly affects net reactions of the form A + B -> products,
    ! which would be off by 10^6 without this
    write (1,fmt='(A)') 'Converting phenomenological rate coefficients to cm^3, mol, and s...'
    do i = 1, nIsom+nReac+nProd
        K(:,:,:,i) = K(:,:,:,i) * 1.0e6 ** (size(net%isomers(i)%species) - 1)
    end do

    ! Write results (to stdout)
    write (1,fmt='(A)') 'Writing results header...'
    call writeOutputIntro(net, nIsom, nReac, nProd, Tlist, nT, Plist, nP, Elist, nGrains, &
        method, K, model, modelOptions, chebyshevCoeffs, pDepArrhenius)

    ! Fit interpolation model
    allocate( chebyshevCoeffs(1:modelOptions(1), 1:modelOptions(2), 1:nIsom+nReac+nProd, 1:nIsom+nReac+nProd) )
    allocate( pDepArrhenius(1:nP, 1:nIsom+nReac+nProd, 1:nIsom+nReac+nProd) )
    if (model == 1) then
        write (1,fmt='(A)') 'Fitting Chebyshev interpolation models...'
    elseif (model == 2) then
        write (1,fmt='(A)') 'Fitting PDepArrhenius interpolation models...'
    end if

    ! The phenomenological rate coefficients
    ! The phenomenological rate coefficients
    ! Note that, although we have generated the k(T,P) values in both
    ! directions, we only need them in one direction (the other comes from
    ! equilibrium)
    ! Thus, we will only save the faster direction as output
    do j = 1, nIsom+nReac+nProd
        do i = 1, j-1
            ! Decide which direction to save
            if (i > nIsom+nReac .and. j > nIsom+nReac) then
                ! Reaction is product channel -> product channel, which by
                ! definition has zero rate, so we skip it
                exit
            elseif (sum(K(:,:,i,j) - K(:,:,j,i)) > 0) then
                reac = j; prod = i
            else
                reac = i; prod = j
            end if
            
            ! Check for validity of phenomenological rate coefficients
            invalidRate = 0
            do t = 1, nT
                do p = 1, nP
                    if (K(t,p,prod,reac) <= 0) then
                        !K(t,p,prod,reac) = 10**(-300)
                        invalidRate = 1
                    end if
                end do
            end do
            if (invalidRate /= 0) then
                write (1,fmt='(A)') 'Warning: One or more k(T,P) values for a net reaction was zero.'
                write (1,fmt='(A)') 'These have been set to 1e-300 to allow for k(T,P) interpolation model fitting.'
                write (1,fmt='(A)') trim(net%isomers(reac)%name)//' -> '//trim(net%isomers(prod)%name)
                do t = 1, nT
                    write (1,fmt=*) K(t,:,prod,reac)
                end do
                stop
            end if
            
            write (*, fmt='(A)'), '# The reactant and product isomers'
            write (*,*), reac, prod
            write (*, fmt='(A)'), '# Table of phenomenological rate coefficients (cm3, mol, s)'
            write (fmtStr,*), '(A8,', nP, 'ES14.2E2)'
            write (*, fmtStr), 'T \ P', Plist
            write (fmtStr,*), '(F8.1,', nP, 'ES14.4E3)'
            do t = 1, nT
                write (*, fmt=fmtStr), Tlist(t), K(t,:,prod,reac)
            end do
            if (model == 1) then
                call fitChebyshevModel(K(:,:,prod,reac), Tlist, Plist, Tmin, Tmax, &
                    Pmin, Pmax, modelOptions(1), modelOptions(2), chebyshevCoeffs(:,:,prod,reac))
                write (*, fmt='(A)'), '# The fitted Chebyshev polynomial model (cm3, mol, s)'
                write (fmtStr,*), '(', modelOptions(1), 'ES14.4E3)'
                do t = 1, modelOptions(1)
                    write (*,fmt=fmtStr), chebyshevCoeffs(t,:,prod,reac)
                end do
            elseif (model == 2) then
                call fitPDepArrheniusModel(K(:,:,prod,reac), Tlist, Plist, pDepArrhenius(:,prod,reac))
                write (*, fmt='(A)'), '# The fitted pressure-dependent Arrhenius model (cm3, mol, s)'
                do p = 1, nP
                    write (*, fmt='(ES8.2E2,ES14.4E3,F14.4,F10.4)'), Plist(p), &
                        pDepArrhenius(p,prod,reac)%A, pDepArrhenius(p,prod,reac)%Ea, pDepArrhenius(p,prod,reac)%n
                end do
            end if
            write (*, fmt='(A)'), ''

        end do
    end do

    ! Close log file
    close(1)

    deallocate(K, chebyshevCoeffs, pDepArrhenius)

end program

