!
! Fortran 90 program to read XML-formatted input file - pey 3/4/04
!
subroutine readxml(kk,ksym,lin,lout,reac,kerr,icase)

   use xmlparse

   integer, parameter :: lendata = 30, ndata = 20

   integer                          :: kk,lin,lout,icase
   logical                          :: ierr, kerr
   character(len=*), dimension(*)   :: ksym
   double precision, dimension(*)   :: reac
   double precision                 :: ctot1, ctot2, diffc
   double precision, dimension(1)   :: value
   character(len=80)                :: line

   logical           :: mustread
   type(XML_PARSE)   :: info

   character(len=30)                           :: tag
   logical                                     :: endtag
   character(len=30), dimension(1:2,1:20)      :: attribs
   integer                                     :: no_attribs
   character(len=lendata), dimension(1:ndata)  :: data
   integer                                     :: no_data
   integer                                     :: i

   double precision :: RTOL,ATOL,TEMP,PRES,TZERO,TSTOP,SPD,CADZERO,RC,RLA,VC

   COMMON /INPUTS/ RTOL,ATOL,TEMP,PRES,TZERO,TSTOP,SPD,CADZERO,RC,RLA,VC

!  initialize info data structure
   mustread = .true.
   call xml_init( info, lin, mustread )
   call xml_options( info, ignore_whitespace = .true. , no_data_truncation = .true.)

   if ( info%error ) then
      kerr = .true.
      write(lout,'(a)') ' Error: xml_init in readxml.f90 failed.'
      return
   endif

   write(lout,'(/,a,/)') 'Echo of xml input file:'

!  loop thru xml input file
   do
	  ierr = .false.

      call xml_get( info, tag, endtag, attribs, no_attribs, data, no_data )
      if ( info%error ) then
	     kerr = .true.
         write(lout,'(a)') ' Error: xml_get in readxml.f90 failed'
         exit
      endif

      if (.not. endtag) then

         ! echo xml read
         write(lout,*) tag
         do i = 1,no_attribs
            write(lout,*) i, '>', attribs(1,i), '<=', trim(attribs(2,i))
         enddo
         write(lout,*) (i, '>',trim(data(i)), '<', i=1,no_data)

         ! match tag
		 if (trim(tag).eq. 'reactortype') then
			if(trim(data(1)).eq.'ic_engine') then
			   icase=1
			elseif (trim(data(1)).eq.'constant_tp') then
			   icase=2
			endif
         elseif (trim(tag).eq.'starttime') then
            call ckxnum(data(1), 1, lout, nval, TZERO, ierr)
         elseif (trim(tag).eq.'endtime') then
		    call ckxnum(data(1), 1, lout, nval, TSTOP, ierr)
	     elseif (trim(tag).eq.'rtol') then
            call ckxnum(data(1), 1, lout, nval, RTOL, ierr)
         elseif (trim(tag).eq.'atol') then
		    call ckxnum(data(1), 1, lout, nval, ATOL, ierr)
	     elseif (trim(tag).eq.'temperature') then
            call ckxnum(data(1), 1, lout, nval, TEMP, ierr)
         elseif (trim(tag).eq.'pressure') then
            call ckxnum(data(1), 1, lout, nval, PRES, ierr)
         elseif (trim(tag).eq.'amount') then
            line = attribs(2,2) // ' ' // data(1)
            call cksnum(line, 1, lout, ksym, kk, kspec, nval, value, ierr)
            if (.not. ierr) reac(kspec) = value(1)
       !  else    
       !     write (lout, '(3a)') ' Warning: Tag ' ,trim(tag), ' was not found.'
         endif

	     if (ierr) then
	        kerr = .true.
		    write(lout, '(3a,/)') ' Error: Reading data for tag ', trim(tag),' failed.'
	     endif

	  endif ! check if not endtag

      if (info%eof) exit ! break loop if eof reached

	  if (.not. xml_ok(info)) then
	     kerr = .true.
	     write(lout,'(a)') ' Error: xml_ok(info) in readxml.f90 returned false.'
		 exit
	  endif

   enddo
!  end of input

   if (kerr) return

!  Check that concentrations are consistent with T and P
   ctot1 = PRES/(8.314d6*TEMP)

   ctot2 = 0.0d0
   do k=1,kk
      ctot2 = ctot2 + reac(k)
   enddo

   diffc = abs(ctot1-ctot2)/(0.5d0*(ctot1+ctot2))
   if(diffc .gt. 0.1) then
      kerr =.true.
      write (lout, '(a)') ' Warning: Concentrations, T, P inconsistent.'
   endif

   return
end subroutine readxml
