!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	model.f90
!
!	Written by Josh Allen (jwallen@mit.edu)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module ModelModule

    use NetworkModule

    implicit none

contains

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! Subroutine: fitChebyshevModel()
    !
    ! Fits a set of rate coefficients k(T, P) fitted at multiple temperatures
    ! and pressures to a model of the form
    !
    !		log k(T, P) = sum_i sum_j alpha_ij * phi_i(Tred) * phi_j(Pred)
    !
    ! where
    !
    !		Tred = (2/T - 1/Tmin - 1/Tmax) / (1/Tmax - 1/Tmin)
    !		Pred = (2 log P - log Pmin - log Pmax) / (log Pmax - log Pmin)
    !
    ! and phi_n(x) is the nth Chebyshev polynomial evaluated at x.
    !
    ! Parameters:
    !		K - Matrix of phenomenological rate coefficients.
    !		Tlist - Vector of absolute temperatures in K.
    !		Plist - Vector of pressures in bar.
    !		nChebT - Number of temperatures to use to fit Chebyshev polynomials
    !		nChebP - Number of pressures to use to fit Chebyshev polynomials
    !		alpha - Matrix of coefficients.
    subroutine fitChebyshevModel(K, Tlist, Plist, Tmin, Tmax, Pmin, Pmax, nChebT, nChebP, alpha)

        real(8), dimension(:,:), intent(in) :: K
        real(8), dimension(:), intent(in) :: Tlist
        real(8), dimension(:), intent(in) :: Plist
        real(8), intent(in) :: Tmin, Tmax, Pmin, Pmax
        integer, intent(in) :: nChebT
        integer, intent(in) :: nChebP
        real(8), dimension(:,:), intent(inout) :: alpha

        real(8), dimension(:), allocatable :: Tred
        real(8), dimension(:), allocatable :: Pred

        real(8), dimension(:,:), allocatable :: A
        real(8), dimension(:,:), allocatable :: b

        integer nT, nP, t, p, t1, p1, t2, p2
        real(8) phiT, phiP

        ! Variables for BLAS and LAPACK
        integer :: info
        integer, dimension(:), allocatable :: work
        integer :: one
        character :: N

        ! Determine number of temperatures and pressures used
        nT = size(Tlist)
        nP = size(Plist)

        allocate( Tred(1:nT), Pred(1:nP) )
        allocate( A(1:nT*nP, 1:nChebT*nChebP), b(1:nT*nP,1) )

        ! Calculate reduced temperatures and pressures
        do t = 1, nT
            Tred(t) = (2.0 / Tlist(t) - 1.0 / Tmin - 1.0 / Tmax) / &
                (1.0 / Tmax - 1.0 / Tmin)
        end do
        do p = 1, nP
            Pred(p) = (2.0 * log(Plist(p)) - log(Pmin) - log(Pmax)) / &
                (log(Pmax) - log(Pmin))
        end do

        ! Populate A matrix and b vector
        ! A_ij = phi_j1(T_i) * phi_j2(P_i)
        ! b_i = log10(k(T_i, P_i))
        do t1 = 1, nT
            do p1 = 1, nP
                do t2 = 1, nChebT
                    do p2 = 1, nChebP
                        call chebyshev(t2-1, Tred(t1), phiT)
                        call chebyshev(p2-1, Pred(p1), phiP)
                        A((p1-1)*nT+t1, (p2-1)*nChebT+t2) = phiT * phiP
                    end do
                end do
            b((p1-1)*nT+t1,1) = log10(K(t1,p1))
            end do
        end do

        ! Solve Ax = b to determine coefficients alpha_ij
        allocate( work(1:8*nChebT*nChebP) )
        one = 1
        N = 'N'
        call DGELS(N, nT*nP, nChebT*nChebP, one, A, nT*nP, b, nT*nP, work, 8*nChebT*nChebP, info)
        if (info > 0) then
            write (*,*), "Chebyshev fit matrix is singular!"
            stop
        end if
        do t = 1, nChebT
            do p = 1, nChebP
                alpha(t,p) = b((p-1)*nChebT+t,1)
            end do
        end do

        ! Clean up
        deallocate( Tred, Pred, A, b, work )


    end subroutine

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! Subroutine: fitPDepArrheniusModel()
    !
    ! Fits a set of rate coefficients k(T, P) fitted at multiple temperatures
    ! and pressures to a model of the form
    !
    !		log k(T, P) = A(P) * T**n(P) * exp( - Ea(P) / (R * T) )
    !
    ! Parameters:
    !		K - Matrix of phenomenological rate coefficients.
    !		Tlist - Vector of absolute temperatures in K.
    !		Plist - Vector of pressures in bar.
    !		arrhenius - Matrix of coefficients.
    subroutine fitPDepArrheniusModel(K, Tlist, Plist, arrhenius)

        real(8), dimension(:,:), intent(in) :: K
        real(8), dimension(:), intent(in) :: Tlist
        real(8), dimension(:), intent(in) :: Plist
        type(ArrheniusKinetics), dimension(:), intent(inout) :: arrhenius

        real(8), dimension(:,:), allocatable :: A
        real(8), dimension(:), allocatable :: b

        integer t, p, nT, nP
        real(8) R

        ! Variables for BLAS and LAPACK
        integer :: info
        integer, dimension(:), allocatable :: work
        integer :: one
        character :: N

        R = 8.314472    ! [=] J/mol*K

        nT = size(Tlist)
        nP = size(Plist)

        allocate( A(1:nT, 1:3), b(1:nT) )
        allocate( work(1:8*nT) )

        ! At each pressure, fit a modified Arrhenius kinetics model
        do p = 1, nP

            ! Populate A matrix and b vector for linear least-squares fit
            do t = 1, nT
                A(t,1) = 1.0
                A(t,2) = log(Tlist(t))
                A(t,3) = -1.0 / (R * Tlist(t))
                b(t) = log(K(t,p))
            end do

            ! Solve Ax = b to determine Arrhenius coefficients
            one = 1
            N = 'N'
            call DGELS(N, nT, 3, one, A, nT, b, nT, work, 8*nT, info)
            if (info > 0) then
                write (*,*), "Log P interpolate fit matrix is singular!"
                stop
            end if

            ! Extract Arrhenius coefficients
            arrhenius(p)%A = exp(b(1))
            arrhenius(p)%n = b(2)
            arrhenius(p)%Ea = b(3)

        end do

        ! Clean up
        deallocate( work, A, b )

    end subroutine


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! Subroutine: chebyshev()
    !
    ! Evaluates the nth Chebyshev polynomial of the first kind at the value x.
    !
    ! Parameters:
    !		n - Index of Chebyshev polynomial
    !		x - Value to evaluate the polynomial at.
    !		phi - Index of nth Chebyshev polynomial
    subroutine chebyshev(n, x, phi)

        integer, intent(in) :: n
        real(8), intent(in) :: x
        real(8), intent(out) :: phi

        phi = cos(n * acos(x));

    end subroutine



end module
