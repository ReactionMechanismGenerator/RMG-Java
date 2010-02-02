c genernal string subroutines for Fortran    AYC 11/10/88
c fixed bug in concat <ayc 2/94>
c
c
c PARSTR
      subroutine parstr(string,nget,substr,ntotal)
c
c    parses input string
c
c    input:   string = list separated by commas or spaces
c             nget   = number of element in list to be returned
c    output:  ntotal = total number of elements
c             substr = substring specified by nget
c
      integer nget,ntotal,iword,istr
      character preschar,prevchar
      character*(*) string,substr
      character*256 tempstr
c
      substr= ' '
      prevchar= ' '
      iword= 1
c
      do 50 istr=1,len(string)
          preschar= string(istr:istr)
          if ((preschar.ne.',').and.(preschar.ne.' ')) then
             if(iword.eq.nget) then
                tempstr= substr
                call concat(tempstr,preschar,substr)
             endif
          else
              if ((prevchar.ne.',').and.(prevchar.ne.' ')) iword=iword+1
          endif
          prevchar=preschar
50    continue
c
      ntotal= iword-1
      call ljust(substr)
c
      return
      end
c
c
c UPCASE
       subroutine upcase(string)
c
c     converts elements of string to upper case
c
       integer i,int
       character*(*) string
c
       do 10 i= 1, len(string)
          int= ichar(string(i:i))
          if ((int.ge.97).and.(int.le.122)) int= int-32
          string(i:i)= char(int)
10     continue
c
       return
       end
c
c
c LJUST
       subroutine ljust(s1)
c
c     left justifies by stripping out leading blanks of string
c
       integer i,istart
       character*(*) s1
c
       do 10 i=1,len(s1)
          istart= i
          if(s1(i:i).ne.' ') goto 15
10     continue
15     continue
c
       s1= s1(istart:len(s1))
c
       return
       end
c
c
c RJUST
       subroutine rjust(s1)
c
c     right justifies by stripping out trailing blanks of string
c
       integer i,iend
       character*(*) s1
c
       do 10 i=len(s1),1,-1
          iend= i
          if(s1(i:i).ne.' ') goto 15
10     continue
15     continue
c
       do 20 i=len(s1),1,-1
          if (i.gt.(len(s1)-iend)) then
              s1(i:i)= s1(i-(len(s1)-iend):i-(len(s1)-iend))
          else
              s1(i:i) = ' '
          endif
20     continue
c
       return
       end
c
c
c CONCAT
       subroutine concat(s1,s2,s3)
       implicit none
c
c     returns s3 = s1 + s2
c     (concatenates s1 and s2, stripping out trailing blanks of s1)
c
       integer i,len2,len1,idum
       character*(*) s1,s2,s3
c
       do 10 i=len(s1),1,-1
          len1= i
          if(s1(i:i).ne.' ') goto 15
10     continue
15     continue
c
       len2 = len(s3) - len1
       if (len2.gt.0) then
          idum = min(len2,len(s2))
          s3= s1(1:len1)//s2(1:idum)
       else
          idum = min(len(s3),len(s1))
          s3 = s1(1:idum)
       endif
c
       return
       end
c
