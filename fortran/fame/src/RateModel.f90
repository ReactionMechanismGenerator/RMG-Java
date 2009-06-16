! ==============================================================================
!
!	RateModel.f90
!
!	Written by Josh Allen (jwallen@mit.edu)
!
! ==============================================================================

module RateModelModule

	use SimulationModule

	implicit none
	
contains

	! --------------------------------------------------------------------------
	!
	! Subroutine: fitRateModel()
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
	subroutine fitRateModel(K, Tlist, Plist, nChebT, nChebP, alpha)
		
		real(8), dimension(:,:), intent(in)		:: K
		real(8), dimension(:), intent(in)		:: Tlist
		real(8), dimension(:), intent(in)		:: Plist
		integer, intent(in)						:: nChebT
		integer, intent(in)						:: nChebP
		real(8), dimension(:,:), intent(inout)	:: alpha
		
		real(8), dimension(:), allocatable		:: Tred
		real(8), dimension(:), allocatable		:: Pred
		
		real(8), dimension(:,:), allocatable	:: A
		real(8), dimension(:,:), allocatable	:: b
		
		integer nT, nP, t, p, t1, p1, t2, p2
		real(8) Tmax, Tmin, Pmax, Pmin
		real(8) phiT, phiP
		
		! Variables for BLAS and LAPACK
		integer									::	info
		integer, dimension(:), allocatable		::	work
		integer									::	one
		character								::	N
		
		! Determine number of temperatures and pressures used
		nT = size(Tlist)
		nP = size(Plist)
		
		allocate( Tred(1:nT), Pred(1:nP) )
		allocate( A(1:nT*nP, 1:nChebT*nChebP), b(1:nT*nP,1) )
		
		! Determine temperature and pressure ranges
		Tmax = maxval(Tlist)
		Tmin = minval(Tlist)
		Pmax = maxval(Plist)
		Pmin = minval(Plist)
		
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
		
	! --------------------------------------------------------------------------
	
	! Subroutine: chebyshev()
	! 
	! Evaluates the nth Chebyshev polynomial of the first kind at the value x.
	!
	! Parameters:
	!		n - Index of Chebyshev polynomial
	!		x - Value to evaluate the polynomial at.
	!		phi - Index of nth Chebyshev polynomial
	subroutine chebyshev(n, x, phi)
		
		integer, intent(in)		:: n
		real(8), intent(in)		:: x
		real(8), intent(out)	:: phi
		
		phi = cos(n * acos(x));

	end subroutine
		
	! --------------------------------------------------------------------------
	
	! Subroutine: saveResults()
	! 
	! Saves the results to the output file
	!
	! Parameters:
	!		path - Path to output file.
	!		simData - The simulation parameters.
	!		alpha - Matrix of Chebyshev coefficients for each reaction.
	subroutine saveResults(path, simData, K, alpha)
		
		character(len=*), intent(in)				:: path
		type(Simulation), intent(in)				:: simData
		real(8), dimension(:,:,:,:), intent(in)		:: K
		real(8), dimension(:,:,:,:), intent(in)		:: alpha
		
		real(8), dimension(:), allocatable			:: rate
		
		integer i, j, t, p, ind
		integer ios
		
		character(len=64) fmtStr
		
		! Create format string
		write (fmtStr,*), '(', size(simData%Plist), 'ES16.4E3)'
		
		write (*,*), '# FAME output'
		write (*,*), 'Number of wells                     ', simData%nIsom
		write (*,*), 'Number of Chebyshev temperatures    ', simData%nChebT
		write (*,*), 'Number of Chebyshev pressures       ', simData%nChebP
		write (*,*), 'Temperature range of fit            ', minval(simData%Tlist), 'K    ', maxval(simData%Tlist), 'K'
		write (*,*), 'Pressure range of fit               ', minval(simData%Plist), 'bar  ', maxval(simData%Plist), 'bar'
		write (*,*), ''
		
		! The order of each row is (1,2), (1,3), ..., (2,1), (2,3), ...
		! (i,j) represents the reaction j --> i
		do i = 1, simData%nIsom
			do j = 1, simData%nIsom
				if (i /= j) then
					write (*,*), '# Reaction', j, '-->', i
					! Actual fitted rate coefficients
					do t = 1, size(simData%Tlist)
						write (*, fmtStr), K(t,:,i,j)
					end do
					! Chebyshev fits
					do t = 1, simData%nChebT
						write (*,*), alpha(t,:,i,j)
					end do
					write (*,*), ''
				end if
			end do
		end do

	end subroutine
		
	! --------------------------------------------------------------------------
	
end module
