!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	fame.f90
!
!	Copyright (c) 2008-2009 by Josh Allen.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program fame

    use NetworkModule
    use IOModule
    
    implicit none

    type(Network) net
    real(8), dimension(:), allocatable :: Tlist, Plist, Elist
    real(8) :: grainSize
    integer :: numGrains
    integer :: method, model
    integer, dimension(1:10) :: modelOptions
    real(8), dimension(:,:,:,:), allocatable :: K
    real(8), dimension(:,:,:,:), allocatable :: chebyshev
    type(ArrheniusKinetics), dimension(:,:,:), allocatable :: pDepArrhenius
    integer :: nIsom, nReac, nProd, nGrains, nT, nP

    integer i

    ! Use unit 1 for logging; file will be called fame.log or fort.1
    open(unit=1, file='fame.log')

    ! Log header

    ! Read network information (from stdin)
    write (1,fmt='(A)') 'Reading network information...'
    call readInput(net, Tlist, Plist, grainSize, numGrains, method, model, modelOptions)
    nIsom = 0
    nReac = 0
    nProd = 0
    do i = 1, size(net%isomers)
        if (size(net%isomers(i)%species) > 1) then
            nReac = nReac + 1
        else
            nIsom = nIsom + 1
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

    ! Fit interpolation model
    !write (1,fmt='(A)') 'Fitting interpolation model...'
    !call network_fitInterpolationModel(net, Tlist, nT, Plist, nP, K, nIsom+nReac+nProd, model, modelOptions)

    ! Write results (to stdout)
    write (1,fmt='(A)') 'Writing results...'
    call writeOutput(net, nIsom+nReac+nProd, Tlist, nT, Plist, nP, Elist, nGrains, &
        method, K, model, modelOptions, chebyshev, pDepArrhenius)

    ! Close log file
    close(1)

    deallocate(K)

end program

