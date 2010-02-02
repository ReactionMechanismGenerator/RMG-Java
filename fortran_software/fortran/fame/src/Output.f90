!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	output.f90
!
!	Written by Josh Allen (jwallen@mit.edu)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module OutputModule

	use SpeciesModule
	use IsomerModule
	use ReactionModule
	use NetworkModule
	
	implicit none
	
contains

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
		
end module
