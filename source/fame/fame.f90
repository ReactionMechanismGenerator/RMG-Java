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

    integer i, j, t, p
    character(len=64) fmtStr
    
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
    do j = 1, nIsom+nReac
        do i = 1, nIsom+nReac+nProd
            if (i /= j) then
                write (*, fmt='(A)'), '# The reactant and product isomers'
                write (*,*), j, i
                write (*, fmt='(A)'), '# Table of phenomenological rate coefficients'
                write (fmtStr,*), '(A8,', nP, 'ES14.2E2)'
                write (*, fmtStr), 'T \ P', Plist
                write (fmtStr,*), '(F8.1,', nP, 'ES14.4E3)'
                do t = 1, nT
                    write (*, fmt=fmtStr), Tlist(t), K(t,:,i,j)
                end do
                if (model == 1) then
                    call fitChebyshevModel(K(:,:,i,j), Tlist, Plist, Tmin, Tmax, &
                        Pmin, Pmax, modelOptions(1), modelOptions(2), chebyshevCoeffs(:,:,i,j))
                    write (*, fmt='(A)'), '# The fitted Chebyshev polynomial model'
                    write (fmtStr,*), '(', modelOptions(1), 'ES14.4E3)'
                    do t = 1, modelOptions(1)
                        write (*,fmt=fmtStr), chebyshevCoeffs(t,:,i,j)
                    end do
                elseif (model == 2) then
                    call fitPDepArrheniusModel(K(:,:,i,j), Tlist, Plist, pDepArrhenius(:,i,j))
                    write (*, fmt='(A)'), '# The fitted pressure-dependent Arrhenius model'
                    do p = 1, nP
                        write (*, fmt='(ES8.2E2,ES14.4E3,F14.4,F10.4)'), Plist(p), &
                            pDepArrhenius(p,i,j)%A, pDepArrhenius(p,i,j)%Ea, pDepArrhenius(p,i,j)%n
                    end do
                end if
                write (*, fmt='(A)'), ''
            end if
        end do
    end do
    
    ! Close log file
    close(1)

    deallocate(K, chebyshevCoeffs, pDepArrhenius)

end program

