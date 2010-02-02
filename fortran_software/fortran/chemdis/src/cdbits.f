c this is all matrix solving stuff from press et al., numerical recipies

c*****ludcmp*************************************************************
      subroutine ludcmp(A,n,np,indx,d,vv,ierr,lout)
c    this is adaped from numerical recipes by press et al.  <ayc 3/93>
c    this is an lu decomposer and the output is passed to susequent routines

c    Note: change of inputs - esp array vv so it can be fixed in size
c    from above

c    A is input array of n equations and is physically dimensioned np x np
c    note input matrix A is changed by this routine to give LU decomp'ed A

c    indx is output 1-D array of pivot info - this along with variable D is
c    passed to companion subroutines

c    note: code assumes no-trip do's (fortran 77 standard)

      implicit none
      real*8 tiny
      parameter(tiny = 1.0D-70)
	
      integer n,np,indx(np),imax,ierr,lout
      integer i,j,k
      real*8 A(np,np),vv(np),d,aamax,sum,dum
	
      ierr = 0
        
      d = 1.d0
      do 120 i = 1,n
         aamax = 1.d0
           
         do 110 j = 1,n
            if (dabs(A(i,j)).gt.aamax) aamax = dabs(A(i,j))
110      continue

         if (aamax.le.tiny) then
            write(lout,*) 'ERROR: singular matrix.'
            ierr = 1
            return
         endif
	   
         vv(i) = 1.d0/aamax
120   continue

      do 190 j = 1,n
         do 140 i = 1,j-1
            sum = A(i,j)
              
            do 130 k= 1,i-1
               sum = sum - A(i,k)*A(k,j)
130         continue

            A(i,j) = sum
140      continue

         aamax = 0.d0
         do 160 i = j,n
            sum = A(i,j)
               
            do 150 k = 1,j-1
               sum = sum - A(i,k)*A(k,j)
150         continue

            A(i,j) = sum
            dum = vv(i)*dabs(sum)
               
            if (dum.ge.aamax) then
               imax = i
               aamax = dum
            endif
160      continue

         if (j.ne.imax) then
            do 170 k = 1,n
               dum = A(imax,k)
               A(imax,k) = A(j,k)
               A(j,k) = dum
170         continue

            d = -d
            vv(imax) = vv(j)
         endif
            
         indx(j) = imax
         if (dabs(A(j,j)).le.tiny) A(j,j) = dsign(tiny,A(j,j))
            
         if (j.ne.n) then
            dum = 1.d0/A(j,j)
               
            do 180 i = j+1,n
               A(i,j) = A(i,j)*dum
180         continue
         endif
          
190   continue

      return
      end
         
c*****lubksb*************************************************************
      subroutine lubksb(A,n,np,indx,B)
c    this is taken from numerical recipes by press et al.   <ayc 3/93>
         
c    A is LU decomp'ed array of n equations, physically dimensioned np x np
c    and was obtained by ludcmp

c    indx is 1-D array of pivot info from lubksb

c    B is input column vector (A X = B) but is changed on output to give X

      implicit none
	
      integer n,np,indx(np),ii,ll
      integer i,j
      real*8 A(np,np),B(np),sum
        
      ii = 0
      do 120 i = 1,n
         ll = indx(i)
         sum = B(ll)
         B(ll) = B(i)
                 
         if (ii.ne.0) then
            do 110 j = ii,i-1
                sum = sum - A(i,j)*B(j)
110         continue

         elseif (sum.ne.0.d0) then
            ii = i
         endif
           
         B(i) = sum
120   continue

      do 140 i = n,1,-1
         sum = B(i)
           
         do 130 j = i+1,n
            sum = sum - A(i,j)*B(j)
130      continue

         B(i) = sum/A(i,i)
140   continue

      return
      end

c*****lfit*************************************************************
      subroutine lfit(x,y,sig,ndata,a,ma,lista,mfit,covar,ncvm,
     2   chisq,funcs,lerr)
c    this is taken from numerical recipes by press et al.   <ayc 4/93>
c    its pretty complicated so you'd better see description in book
c    comments in code are mine <ayc>

c    this program allows you to only adjust mfit of ma coefficients, 
c    while keeping the rest fixed - lista is an input vector specifying
c    what parameters you want adjusted

      implicit none
      
c    mmax is max # of coefficients
      integer mmax
      parameter (mmax = 10)

      integer ndata,ma,lista(ma),mfit,ncvm
      integer kk,ihit,i,j,k   
      real*8 x(ndata),y(ndata),sig(ndata),a(ma),covar(ncvm,ncvm),
     2   beta(mmax),afunc(mmax),chisq
      real*8 ym,sig2i,wt,sum
      external funcs
      
c    the following are for ldu routines which we now use
c    instead of gaussj
      integer indx(mmax),ierr,lerr
      real*8 d,vv(mmax)
c    done declarations
     
      kk = mfit + 1
      
      do 12 j = 1,ma
         ihit = 0
         
         do 11 k = 1,mfit
            if (lista(k).eq.j) ihit = ihit + 1
11       continue

         if (ihit.eq.0) then
            lista(kk) = j
            kk = kk + 1
         else if (ihit.gt.1) then
            stop 'improper set in lista'
         endif
12    continue

      if (kk.ne.(ma+1)) stop 'improper set in lista'
      
      do 14 j = 1,mfit
         do 13 k = 1,mfit
            covar(j,k) = 0.d0
13       continue
         beta(j) = 0.d0
14    continue

      do 18 i = 1,ndata
         call funcs(x(i),afunc,ma)
         ym = y(i)
         
         if (mfit.lt.ma) then
            do 15 j = mfit+1,ma
               ym = ym - a(lista(j))*afunc(lista(j))
15          continue
         endif
         
         sig2i = 1.d0/sig(i)**2
         
         do 17 j = 1,mfit
            wt = afunc(lista(j))*sig2i
            do 16 k = 1,j
               covar(j,k) = covar(j,k) + wt*afunc(lista(k))
16          continue
            beta(j) = beta(j) + ym*wt
17       continue
18    continue

      if (mfit.gt.1) then
         do 21 j = 2,mfit
            do 19 k = 1,j-1
               covar(k,j) = covar(j,k)
19          continue
21       continue
      endif

c    we have replaced this with call to ldu solver      
c     call gaussj(covar,mfit,ncvm,beta,1,1)
      
c    we use lud routines instead of gaussj
c    we aren't going to do anything with inverse
c    so we just solve to get solution vector

c    NOTE: we have changed some of the inputs of these
c    routines as compared to book!!
      call ludcmp(covar,mfit,ncvm,indx,d,vv,ierr,lerr)
      call lubksb(covar,mfit,ncvm,indx,beta)
      
c    this has totally messed up covar but we don't
c    need it anymore - see book to get the 3 lines
c    of code to compute inverse of covar from here

      do 22 j = 1,mfit
         a(lista(j)) = beta(j)
22    continue

      chisq = 0.d0
      do 24 i = 1,ndata
         call funcs(x(i),afunc,ma)
         
         sum = 0.d0      
         do 23 j = 1,ma
            sum = sum + a(j)*afunc(j)
23       continue
         chisq = chisq + ((y(i)-sum)/sig(i))**2
24    continue

c     call covsrt(covar,ncvm,ma,lista,mfit)
c    this sorts covariant matrix to original ordering
c    of fitting coefficients; here we neither changed
c    ordering nor do we care about this matrix so
c    we don't worry about it - note with gaussj,
c    covar would have been the inverse of the original;
c    here its become befuddled - but if we cared we
c    could fix it

      return
      end   
      
c*****gaussj*************************************************************

c    again, same source - we were having problems with this before
c    but it works now - still, ldu routines are faster, since
c    they don't compute inverse which we don't use anyway
 
c    we typed this in so we'll leave it here as an alternative

      subroutine gaussj(a,n,np,b,m,mp)
      implicit none   
      include 'cdparams.fh'
c     this is HORRIBLE!  dmm ... fixed 06-17-2003.  THIS is FINALLY 
c     the HIDDEN 12 ISOMER PROBLEM !!! ARRGH! 
      
c    note - set nmax to >= mxwells
c      integer nmax
c      parameter (nmax = 50)

      integer n,np,m,mp,i,j,k,irow,icol,l,ll
c      integer ipiv(nmax),indxr(nmax),indxc(nmax)      
c      real*8 a(np,np),b(np,mp)
c      real*8 big,dum,pivinv

      integer ipiv(mxWells),indxr(mxWells),indxc(mxWells)      
      real*8 a(np,np),b(np,mp)
      real*8 big,dum,pivinv
      
      do 11 j = 1,n
         ipiv(j) = 0
11    continue

      do 22 i = 1,n
         big = 0.d0
         do 13 j = 1,n
         
            if (ipiv(j).ne.1) then
               do 12 k = 1,n
               
                  if (ipiv(k).eq.0) then
                     if (abs(a(j,k)).ge.big) then
                        big = abs(a(j,k))
                        irow = j
                        icol = k
                     endif
                     
                  else if (ipiv(k).gt.1) then
                     return
                  endif
12             continue
            endif   
13       continue

         ipiv(icol) = ipiv(icol) + 1
         if (irow.ne.icol) then
            do 14 l = 1,n 
               dum = a(irow,l)
               a(irow,l) = a(icol,l)
               a(icol,l) = dum
14          continue   
           
            do 15 l = 1,m
               dum = b(irow,l)
               b(irow,l) = b(icol,l)
               b(icol,l) = dum
15          continue
         endif
         
         indxr(i) = irow
         indxc(i) = icol
         if (a(icol,icol).eq.0.d0) return
         
         pivinv = 1.d0/a(icol,icol) 
         a(icol,icol) = 1.d0        
         do 16 l = 1,n
            a(icol,l) = a(icol,l)*pivinv
16       continue  
         
         do 17 l = 1,m
            b(icol,l) = b(icol,l)*pivinv
17       continue

         do 21 ll = 1,n
            if (ll.ne.icol) then
               dum = a(ll,icol)
               a(ll,icol) = 0.d0
               
               do 18 l = 1,n
                  a(ll,l) = a(ll,l) - a(icol,l)*dum
18             continue

               do 19 l = 1,m
                  b(ll,l) = b(ll,l) - b(icol,l)*dum
19             continue
            endif
21       continue
22    continue

      do 24 l = n,1,-1
         if (indxr(l).ne.indxc(l)) then
            do 23 k = 1,n
               dum = a(k,indxr(l))
	       a(k,indxr(l)) = a(k,indxc(l))
               a(k,indxc(l)) = dum
23          continue
         endif
24    continue

      return
      end
            
     
c $Id$
c $Author$
c $Date$
c $Source$
c $Revision$
c $Log$
c Revision 1.1  2007-02-20 23:10:23  sandeeps
c Initial revision
c
c Revision 1.3  2003/06/17 20:10:27  dmmatheu
c fixed UNBELIEVABLE HIDDEN '12' error in gaussj!  I have been looking
c for this !@#$I#" error for over 3 years ...
c
c Revision 1.2  2000/04/04 17:34:00  dmmatheu
c changed 'tiny' to 1.d-70 instead of 1.d-175 to avoid compile time
c warning
c
c Revision 1.1  2000/04/04 17:33:14  dmmatheu
c Initial revision
c
