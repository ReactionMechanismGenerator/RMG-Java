c	Calculate coefficients for a Chebyshev polynomial that approximates
c	k(T,P). Actually the fitting is for log10(k(T,P)) as a fucntion of
c	log10(P) and 1/T. The fitting is done by direct calculation of the 
c	coefficients rather than least-squares fitting. This requires k(T,P)
c	be known at the Gauss-Chebyshev values of T,P (calculated in cdgets.f).
c
c	pey 4/6/04
c
	subroutine chebyshev(iwell,iprod)
c	calculate the chebyshev coefficents and store them in the common array, chebcoeffs	

	implicit none
        include 'cdparams.fh'
	include 'cdcheb.fh'
	include 'cdrange.fh'
	include 'cdrates.fh'

	integer iwell, iprod, i, j, m, n
	real*8 amn, factor, phi, logk(mxTPts,mxPPts)
c	end of declarations

c	calculate log10(k) and store
	do i=1,ntemps
		do j=1,npres
			logk(i,j) = log10(RK(iwell,iprod,i,j))
		enddo
	enddo

c	calculate the chebyshev coefficients
	factor = 4.0/(ntemps*npres)

	do i=1,nchebT
		do j=1,nchebP
			amn = 0.0
			do m=1,ntemps
				do n=1,npres
				 amn = amn + logk(m,n)*phi(i,tcheb(m))*phi(j,pcheb(n))
				enddo
			enddo
			chebcoeffs(i,j) = factor*amn
		enddo
	enddo

c	correct factor
	do i=1,nchebT
		chebcoeffs(i,1) = chebcoeffs(i,1)/2.0
	enddo

	do j=1,nchebP
		chebcoeffs(1,j) = chebcoeffs(1,j)/2.0
	enddo

	return
	end
c ***
c ***
c ***
c ***
	real*8 function phi(i,x)
c	calculate the ith order Chebyshev basis function
	implicit none
	real*8 x
	integer i

	phi = cos(float(i-1)*acos(x))

	return
	end
