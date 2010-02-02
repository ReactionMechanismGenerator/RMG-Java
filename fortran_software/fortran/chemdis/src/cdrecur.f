c****RECUR1**************************************
      function recur1(Eavail,iwell)
c    to conform with MORONIC FORTRAN RULES: must have mxfreqs-2 versions
c    of this subroutine
c    each version of recur[iver] handles case of freq(iver+1)
c    version mxfreqs-2 calls subroutine RECURN to finishup

c    TO CLONE: copy; then replace every occurence of RECURn with
c      RECURn+1 and increment iver in parameter statement

      implicit none
      include 'cdparams.fh'
      include 'cdisprop.fh'
      
      integer iver
      parameter (iver = 1)
      
c    local variables
      real*8 recur1,recur2,termV,termVR
      real*8 Eavail,sum,qRes1,vu,xx
      integer iwell,i
c    end declarations

      vu = 2.859d-3*freq(iwell,iver+1)
      sum = 0.d0
      
c    last recursion; IMPORTANT!:  special treatment handles 1st freq
c    since this is the one the counting is based on

      if (iver+1.eq.nfreqs(iwell)) then
            
         do 10 i = 0,int(Eavail/vu)         
            qRes1 =  (Eavail - vu*dfloat(i)) / (2.859d-3*freq(iwell,1))

            sum = sum + termV(i, degen(iwell,iver+1))
     2         * termVR(int(qRes1), degen(iwell,1), nRot)
10       continue

c    general recursion 
 
      else            
         do 20 i = 0,int(Eavail/vu)
            xx = degen(iwell, iver+1)
            sum = sum + termV(i, degen(iwell,iver+1))
     2          * recur2(Eavail - vu*dfloat(i), iwell)   
20       continue
      endif

      recur1 = sum
      
      return
      end
      
c****RECUR2**************************************
      function recur2(Eavail,iwell)
c    to conform with MORONIC FORTRAN RULES: must have mxfreqs-2 versions
c    of this subroutine
c    each version of recur[iver] handles case of freq(iver+1)
c    version mxfreqs-2 calls subroutine RECURN to finishup

c    TO CLONE: copy; then replace every occurence of RECURn with
c      RECURn+1 and increment iver in parameter statement


      implicit none
      include 'cdparams.fh'
      include 'cdisprop.fh'
      
      integer iver
      parameter (iver = 2)

c    local variables
      real*8 recur2,recur3,termV,termVR
      real*8 Eavail,sum,qRes1,vu
      integer iwell,i
c    end declarations

      vu = 2.859d-3*freq(iwell,iver+1)
      sum = 0.d0
      
c    last recursion; IMPORTANT!:  special treatment handles 1st freq
c    since this is the one the counting is based on

      if (iver+1.eq.nfreqs(iwell)) then
            
         do 10 i = 0,int(Eavail/vu)         
            qRes1 =  (Eavail - vu*dfloat(i)) / (2.859d-3*freq(iwell,1))
                 
            sum = sum + termV(i, degen(iwell,iver+1))
     2         * termVR(int(qRes1), degen(iwell,1), nRot)
10       continue

c    general recursion 
 
      else            
         do 20 i = 0,int(Eavail/vu)       
            sum = sum + termV(i, degen(iwell,iver+1))
     2          * recur3(Eavail - vu*dfloat(i), iwell)   
20       continue
      endif

      recur2 = sum
      
      return
      end
      
c****RECUR3**************************************
      function recur3(Eavail,iwell)
c    to conform with MORONIC FORTRAN RULES: must have mxfreqs-2 versions
c    of this subroutine
c    each version of recur[iver] handles case of freq(iver+1)
c    version mxfreqs-2 calls subroutine RECURN to finishup

c    TO CLONE: copy; then replace every occurence of RECURn with
c      RECURn+1 and increment iver in parameter statement
   
      implicit none
      include 'cdparams.fh'
      include 'cdisprop.fh'

      integer iver
      parameter (iver = 3)
      
c    local variables
      real*8 recur3,recur4,termV,termVR
      real*8 Eavail,sum,qRes1,vu
      integer iwell,i
c    end declarations

      vu = 2.859d-3*freq(iwell,iver+1)
      sum = 0.d0
      
c    last recursion; IMPORTANT!:  special treatment handles 1st freq
c    since this is the one the counting is based on

      if (iver+1.eq.nfreqs(iwell)) then
            
         do 10 i = 0,int(Eavail/vu)         
            qRes1 =  (Eavail - vu*dfloat(i)) / (2.859d-3*freq(iwell,1))
                 
            sum = sum + termV(i, degen(iwell,iver+1))
     2         * termVR(int(qRes1), degen(iwell,1), nRot)
10       continue

c    general recursion 
 
      else            
         do 20 i = 0,int(Eavail/vu)       
            sum = sum + termV(i, degen(iwell,iver+1))
     2          * recur4(Eavail - vu*dfloat(i), iwell)   
20       continue
      endif

      recur3 = sum
      
      return
      end
      
c****RECUR4**************************************
      function recur4(Eavail,iwell)
c    to conform with MORONIC FORTRAN RULES: must have mxfreqs-2 versions
c    of this subroutine
c    each version of recur[iver] handles case of freq(iver+1)
c    version mxfreqs-2 calls subroutine RECURN to finishup

c    TO CLONE: copy; then replace every occurence of RECURn with
c      RECURn+1 and increment iver in parameter statement

      implicit none
      include 'cdparams.fh'
      include 'cdisprop.fh'
      
      integer iver
      parameter (iver = 4)
      
c    local variables
      real*8 recur4,recurN,termV,termVR
      real*8 Eavail,sum,qRes1,vu
      integer iwell,i
c    end declarations

      vu = 2.859d-3*freq(iwell,iver+1)
      sum = 0.d0
      
c    last recursion; IMPORTANT!:  special treatment handles 1st freq
c    since this is the one the counting is based on

      if (iver+1.eq.nfreqs(iwell)) then
            
         do 10 i = 0,int(Eavail/vu)         
            qRes1 =  (Eavail - vu*dfloat(i)) / (2.859d-3*freq(iwell,1))
                 
            sum = sum + termV(i, degen(iwell,iver+1))
     2         * termVR(int(qRes1),degen(iwell,1), nRot)
10       continue

c    general recursion 
 
      else            
         do 20 i = 0,int(Eavail/vu)       
            sum = sum + termV(i, degen(iwell,iver+1))
     2          * recurN(Eavail - vu*dfloat(i), iwell)   
20       continue
      endif

      recur4 = sum
      
      return
      end
      
c****RECURN**************************************
      function recurN(Eavail,iwell)

c    This is a special case: the last possible recursive subroutine
c     - called from the last recur[iver] in the series (corresponding
c     to iver = mxfreqs-2) 
c    this routine is effectively iver = mxfreqs-1 
      
      implicit none
      include 'cdparams.fh'
      include 'cdisprop.fh'
      
      integer iver
      parameter (iver = mxfreqs-1)
      
c    local variables
      real*8 recurN,termV,termVR
      real*8 Eavail,sum,qRes1,vu
      integer iwell,i
c    end declarations

      vu = 2.859d-3*freq(iwell,iver+1)
      sum = 0.d0
      
c    LAST ITERATION!

      if (iver+1.eq.nfreqs(iwell)) then
            
         do 10 i = 0,int(Eavail/vu)         
            qRes1 =  (Eavail - vu*dfloat(i)) / (2.859d-3*freq(iwell,1))
                 
            sum = sum + termV(i, degen(iwell,iver+1))
     2         * termVR(int(qRes1), degen(iwell,1), nRot)
10       continue

      else
         stop ' NO RECURx'
      endif

      recurN = sum
      
      return
      end
      
