
      program callfit3pbnd

      implicit double precision (a-h,o-z), integer(i-n)
      parameter (mamax=100, mg=1, mdw=mamax+mg, n=3, k=mamax+mg) 
      parameter (lin=10, lout=11)
      dimension rate(mamax), temp(mamax)
      dimension w(mdw,n+1), ws(2*n+k+(mg+2)*(n+7)), x(n), ip(mg+2*n+2)

c     Purpose: calling program to test fit3pbnd subroutine

c     The input.dat_orig file has rate data for the rxn h + o2 = oh + o
c     from Miller (1997). T has units of Kelvin and k has units of
c     L/mol-s. The rate parameters should be A = 6.71e8 L/mol-s, 
c     n = 0.55, Ea = 11.86 kcal/mol. Source: NIST kinetic database.

c     The input.dat file has data which give a negative activation energy
c     when fit with the original unbounded fit3p.f subroutine.
c     fit3p.f results (unbounded):  A=6.71e8,  n=2.000, Ea=-3.97
c     fit3pbnd.f results (bounded): A=1.25e12, n=1.153, Ea=0.0

c     open files for i/o
      open(lin, status='unknown', file='input.dat')
      open(lout, status='unknown', file='output.dat')

c     read input data (column 1 is T, column 2 is k)
      i = 1
 5    if (i.le.mamax) then ! only read first mamax points
         read(lin, *, end=15) temp(i), rate(i) ! break when eof reached
         i = i+1
         goto 5
      endif
 15   ma = i-1

c     calculate A, n, Tarr
      call fit3pbnd(ma, temp, rate, afactor, an, tarr, 
     >              w, ws, x, ip, mg, mdw, n, mode)

c     check if solution was found
      if (mode.gt.0) then
         write(lout,*) 'No solution found in fit3pbnd: mode =', mode
         stop
      endif

c     convert Arrhenius temperautre to activation energy
      r  = 1.987d-3 ! kcal/mol-K
      ea = tarr*r   ! kcal/mol

c     print output and check goodness of fit
      write(lout,'(A, e13.5e3)') ' A [L/mol-s] = ', afactor
      write(lout,'(A, f12.4)')  ' n = ', an
      write(lout,'(A, f12.3)') ' Ea [kcal/mol] = ', ea
      write(lout,*)

      write(lout,'(A)') ' Temp     k           kfit'
      do i=1,ma
         t     = temp(i)
         rk    = rate(i)
         rkfit = afactor * t**an * exp(-tarr/t)

         write(lout,'(f7.1, e12.3, e12.3)') t, rk, rkfit
      end do

c     close files for i/o
      close(lin)
      close(lout)
         
      stop
      end
c
c
c
c
      subroutine fit3pbnd(ma, temp, ak,  afactor, an, tarr, 
     >                    w, ws, x, ip, mg, mdw, n, mode)

      implicit double precision (a-h,o-z), integer(i-n)
      dimension ak(*), temp(*), w(mdw,*), ws(*), x(*), ip(*), prgopt(1)

c     Purpose: Fit Arrhenius parameters (A, n, Tarr) given data
c              for the rate constant at various temperatures. 
c
c              Calls the DLSEI subroutine from the SLATEC library to solve
c              the constrained linear least-squares problem. The Arrhenius
c              temperature (Ea/R) is constrained to be non-negative. DLSEI.f
c              was downloaded from www.netlib.org on 8 March 04.
c
c     Inputs:  ma   = the number of data points for k and T
c              temp = a vector of temperatures
c              ak   = a vector of rate constants
c              w    = least-squares matrix and constraints
c              ws   = real workspace for LSEI
c              x    = solution vector x=[ln(A),n,Tarr]'
c              ip   = integer workspace for LSEI
c              mg   = number of inequality constraints
c              mdw  = first dimension of w (should be >= ma+me+mg)
c              n    = number of parameters in fitting
c
c     Outputs: afactor = the pre-exponential factor, A [same units as k]
c              an      = the temperature exponent, n
c              tarr    = Arrhenius temperature [same units as temp]
c              mode    = value returned by LSEI indicating if sol'n was found
c
c     Paul Yelvington, 8 March 04

c     construct the w matrix

c     equations to be solved in a least-squares sense
      do i=1,ma
         w(i,1) = 1.0d0
         w(i,2) = log(temp(i))
         w(i,3) = -1.0d0/temp(i)
         w(i,4) = log(ak(i))
      enddo

c     inequaility constraint (Tarr>=0)
      w(ma+1,1) = 0
      w(ma+1,2) = 0
      w(ma+1,3) = 1
      w(ma+1,4) = 0
      
c     parameters for LSEI
      me        = 0   ! number of equality constraints (equals 0)
      prgopt(1) = 1   ! use default parameters for LSEI

c     call LSEI to solve the constrained lls problem
      call dlsei(w, mdw, me, ma, mg, n, prgopt, x, rnorme, 
     >           rnorml, mode, ws, ip)

c     calculate fitted Arrhenius parameters
      afactor = exp(x(1))
      an      = x(2)
      tarr    = x(3)

      return
      end
