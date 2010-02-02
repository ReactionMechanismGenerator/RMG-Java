      SUBROUTINE mrqmin(x,y,sig,ndata,a,ia,ma,covar,alpha,nca,chisq,
     *funcs,alamda, amin, amax)

c     modified by dmmatheu to take double precision variables and the
c     amin, amax bounds for the parameters

      INTEGER ma,nca,ndata,ia(ma),MMAX
      REAL*8 alamda,chisq,funcs,a(ma),alpha(nca,nca),covar(nca,nca),
     *sig(ndata),x(ndata),y(ndata), amin(ma), amax(ma)
      EXTERNAL funcs
      PARAMETER (MMAX=20)
CU    USES covsrt,gaussj,mrqcof
      INTEGER j,k,l,mfit
      REAL*8 ochisq,atry(MMAX),beta(MMAX),da(MMAX)
      integer ibound
      SAVE ochisq,atry,beta,da,mfit
      if(alamda.lt.0.)then
        mfit=0
        do 11 j=1,ma
          if (ia(j).ne.0) mfit=mfit+1
11      continue
        alamda=0.001
        call mrqcof(x,y,sig,ndata,a,ia,ma,alpha,beta,nca,chisq,funcs)
        ochisq=chisq
        do 12 j=1,ma
          atry(j)=a(j)
12      continue
      endif
      do 14 j=1,mfit
        do 13 k=1,mfit
          covar(j,k)=alpha(j,k)
13      continue
        covar(j,j)=alpha(j,j)*(1.+alamda)
        da(j)=beta(j)
14    continue
      call gaussj(covar,mfit,nca,da,1,1)
      if(alamda.eq.0.)then
        call covsrt(covar,nca,ma,ia,mfit)
        call covsrt(alpha,nca,ma,ia,mfit)
        return
      endif
      j=0
      do 15 l=1,ma
        if(ia(l).ne.0) then
          j=j+1
          atry(l)=a(l)+da(j)
        endif
15    continue
      call mrqcof(x,y,sig,ndata,atry,ia,ma,covar,da,nca,chisq,funcs)
      call checkbounds(atry, amin, amax, ma, ibound)

c      if( (chisq.lt.ochisq).or.(ibound.gt.0))then
      if( (chisq.lt.ochisq).and.(ibound.eq.0))then

c        alamda=0.1*alamda
        alamda=0.5*alamda

c     the following logic copied from the original therfit.  I have no
c     idea why this should work or if it does at all.  dmm
c        if( ibound.gt.0) then
c           alamda = ((0.5)**(dfloat(ibound-1)))*alamda
c        endif
        ochisq=chisq
        do 17 j=1,mfit
          do 16 k=1,mfit
            alpha(j,k)=covar(j,k)
16        continue
          beta(j)=da(j)
17      continue
        do 18 l=1,ma
          a(l)=atry(l)
18      continue
      else
c        alamda=10.*alamda
        alamda=2.*alamda
        if( ibound.gt.0) then
           alamda = ((2.0)**(dfloat(ibound-1)))*alamda
        endif

        chisq=ochisq
      endif
      return
      END


c $Id$
c $Author$
c $Date$
c $Source$
c $Revision$
c $Log$
c Revision 1.1  2007-02-20 23:10:22  sandeeps
c Initial revision
c
c Revision 1.1  2000/06/07 15:09:40  dmmatheu
c Initial revision
c
