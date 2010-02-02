! xmlparse.for - Simple, limited XML parser in Fortran
!
! Arjen Markus
!
! General information:
! The module reads XML files by:
! - Identifying the tag and all attributes and data belonging
!   to the tag.
! - Returning to the calling subprogram to let it take care of
!   the tag, attributes and data.
! - If the tag is actually an ending tag, then this is flagged
!   too.
! - Handling all the data is left to the calling subprogram,
!   the module merely facilitates in the parsing.
!
! Note:
! The module in its current version has a number of limitations:
! - It does not handle escape sequences (like &gt. to signify
!   a ">" sign)
! - It does not handle tags with attributes that are spread
!   over more than one line
! - The maximum length of a line is 1000 characters
! - It may report too many lines of data (empty lines)
! - No DOM support nor support for an object tree
! - It is probably not very robust in detecting malformed XML files
!
! Some questions:
! - What to do with leading blanks?
!
! Update - several ideas:
! - Introduce at least two options (via xml_options):
!   - ignore_whitespace  - remove leading blanks and leading and trailing
!                          empty lines from the PCDATA
!   - no_data_truncation - consider truncation of data (more
!                          attributes or lines of character data than
!                          can be stored) a read error
! - Introduce convenience functions and subroutines:
!   - xml_ok()           - all is well, reading can continue
!   - xml_data_trunc()   - was there truncation of the data?
!   - xml_find_attrib()  - find an attribute by name
!
! Further ideas:
!   - simple checking via a table: parent, tag, id, min, max
!
module xmlparse

   implicit none

   integer, parameter :: XML_BUFFER_LENGTH = 1000
   !
   ! Define the data type that holds the parser information
   !
   type XML_PARSE
      integer          :: lun                ! LU-number of the XML-file
      integer          :: level              ! Indentation level (output)
      logical          :: ignore_whitespace  ! Ignore leading blanks etc.
      logical          :: no_data_truncation ! Do not allow data truncation
      logical          :: too_many_attribs   ! More attributes than could be stored?
      logical          :: too_many_data      ! More lines of data than could be stored?
      logical          :: eof                ! End of file?
      logical          :: error              ! Invalid XML file or other error?
      character(len=XML_BUFFER_LENGTH) :: line  ! Buffer
   end type XML_PARSE

contains

subroutine xml_open( info, fname, mustread )

   character(len=*), intent(in)     :: fname
   logical,          intent(in)     :: mustread
   type(XML_PARSE),  intent(out)    :: info

   integer                          :: i
   integer                          :: k
   integer                          :: ierr
   logical                          :: opend
   logical                          :: exists

   info%lun = 10
   info%ignore_whitespace  = .false.
   info%no_data_truncation = .false.
   info%too_many_attribs   = .false.
   info%too_many_data      = .false.
   info%eof                = .false.
   info%error              = .false.

   do i = 10,99
      inquire( unit = i, opened = opend )
      if ( .not. opend ) then
         info%lun = i
         inquire( file = fname, exist = exists )
         if ( .not. exists .and. mustread ) then
            info%lun = -1
            info%error = .true.
         else
            open( unit = info%lun, file = fname )
         endif
         exit
      endif
   enddo

   if ( .not. info%error .and. mustread ) then
      k = 1
      do while ( k .ge. 1 )
         read( info%lun, '(a)', iostat = ierr ) info%line
         if ( ierr .eq. 0 ) then
            info%line = adjustl(  info%line )
            k         = index( info%line, '<?' )
         else
            call xml_close( info )
            info%error = .true.
            exit
         endif
      enddo
   endif
end subroutine xml_open

! initialize the info data structure for a file with unit number = lunit (pey 3/4/04)
subroutine xml_init( info, lunit, mustread )

   logical,          intent(in)     :: mustread
   type(XML_PARSE),  intent(out)    :: info

   integer                          :: lunit
   integer                          :: i
   integer                          :: k
   integer                          :: ierr
   logical                          :: opend
   logical                          :: exists

   info%lun = lunit
   info%ignore_whitespace  = .false.
   info%no_data_truncation = .false.
   info%too_many_attribs   = .false.
   info%too_many_data      = .false.
   info%eof                = .false.
   info%error              = .false.

   if ( .not. info%error .and. mustread ) then
      k = 1
      do while ( k .ge. 1 )
         read( info%lun, '(a)', iostat = ierr ) info%line
         if ( ierr .eq. 0 ) then
            info%line = adjustl(  info%line )
            k         = index( info%line, '<?' )
         else
            call xml_close( info )
            info%error = .true.
            exit
         endif
      enddo
   endif
end subroutine xml_init

subroutine xml_close( info )
   type(XML_PARSE),  intent(inout)    :: info

   close( info%lun )

   info%lun              = -1
   info%too_many_attribs = .false.
   info%too_many_data    = .false.
   info%eof              = .false.
   info%error            = .false.
end subroutine xml_close

subroutine xml_get( info, tag, endtag, attribs, no_attribs, &
                   data, no_data )
   type(XML_PARSE),  intent(inout)               :: info
   character(len=*), intent(out)                 :: tag
   logical,          intent(out)                 :: endtag
   character(len=*), intent(out), dimension(:,:) :: attribs
   integer,          intent(out)                 :: no_attribs
   character(len=*), intent(out), dimension(:)   :: data
   integer,          intent(out)                 :: no_data

   integer         :: kspace
   integer         :: kend
   integer         :: keq
   integer         :: kfirst
   integer         :: ksecond
   integer         :: idxat
   integer         :: idxdat
   integer         :: ierr
   logical         :: close_bracket
   character(len=XML_BUFFER_LENGTH) :: nextline

   !
   ! Initialise the output
   !
   endtag     = .false.
   no_attribs = 0
   no_data    = 0

   info%too_many_attribs = .false.
   info%too_many_data    = .false.

   if ( info%lun .lt. 0 ) then
      return
   endif

   !
   ! From the previous call or the call to xmlopen we have
   ! the line that we need to parse already in memory:
   ! <tag attrib1="..." attrib2="..." />
   !
   close_bracket = .false.
   kspace        = index( info%line, ' ' )
   kend          = index( info%line, '>' )
   do while ( kend .le. 0 )
      read( info%lun, '(a)', iostat = ierr ) nextline
      if ( ierr .eq. 0 ) then
         info%line = trim(info%line) // ' ' // adjustl(nextline)
      else
         info%error = .true.
         call xml_close( info )
         return
      endif
      kend = index( info%line, '>' )
   enddo
   if ( kend .gt. kspace ) then
      kend = kspace
   else
      close_bracket = .true.
   endif

   if ( info%line(1:2) .eq. '</' ) then
      endtag = .true.
      tag    = info%line(3:kend-1)
   else
      if ( info%line(1:1) .eq. '<' ) then
         tag    = info%line(2:kend-1)
      else
         kend   = 0 ! Beginning of data!
      endif
   endif

   info%line = adjustl( info%line(kend+1:) )

   idxat     = 0
   idxdat    = 0

   do while ( info%line .ne. ' ' .and. .not. close_bracket )

      keq  = index( info%line, '=' )
      kend = index( info%line, '>' )
      if ( keq .gt. kend ) keq = 0 ! Guard against multiple tags
                                   ! with attributes on one line

      !
      ! No attributes any more?
      !
      if ( keq .lt. 1 ) then
         kend = index( info%line, '/>' )
         if ( kend .ge. 1 ) then
            kend   = kend + 1 ! To go beyond the ">" character
            endtag = .true.
         else
            kend = index( info%line, '>' )
            if ( kend .lt. 1 ) then
               info%error = .true. ! Wrong ending of line!
               call xml_close( info )
               return
            else
               close_bracket = .true.
            endif
         endif
         if ( kend .ge. 1 ) then
            info%line = adjustl( info%line(kend+1:) )
         endif
         exit
      endif

      idxat = idxat + 1
      if ( idxat .le. size(attribs,2) ) then
         no_attribs = idxat
         attribs(1,idxat) = adjustl(info%line(1:keq-1)) ! Use adjustl() to avoid
                                                     ! multiple spaces, etc
         info%line = adjustl( info%line(keq+1:) )

         !
         ! We have almost found the start of the attribute's value
         !
         kfirst  = index( info%line, '"' )
         if ( kfirst .lt. 1 ) then
            info%error = .true. ! Wrong form of attribute-value pair
            call xml_close( info )
            return
         endif

         ksecond = index( info%line(kfirst+1:), '"' ) + kfirst
         if ( ksecond .lt. 1 ) then
            info%error = .true. ! Wrong form of attribute-value pair
            call xml_close( info )
            return
         endif

         attribs(2,idxat) = info%line(kfirst+1:ksecond-1)
         info%line = adjustl( info%line(ksecond+1:) )
      endif

     !!!!idxat = idxat + 1
      if ( idxat .gt. size(attribs,2) ) then
         info%too_many_attribs = .true.
         info%line             = ' '
         exit
      endif
   enddo

   !
   ! Now read the data associated with the current tag
   ! - all the way to the next "<" character
   !
   ! To do: reduce the number of data lines - empty ones
   ! at the end should not count.
   !
   do
      kend   = index( info%line, '<' )
      idxdat = idxdat + 1
      if ( idxdat .le. size(data) ) then
         no_data = idxdat
         if ( kend .ge. 1 ) then
            data(idxdat) = info%line(1:kend-1)
            info%line    = info%line(kend:)
         else
            data(idxdat) = info%line
         endif
      else
         info%too_many_data = .true.
         exit
      endif

      !
      ! No more data? Otherwise, read on
      !
      if ( kend .ge. 1 ) then
         exit
      else
         read( info%lun, '(a)', iostat = ierr ) info%line
         if ( ierr .lt. 0 ) then
            info%eof = .true.
         elseif ( ierr .gt. 0 ) then
            info%error = .true.
         endif
         if ( ierr .ne. 0 ) then
            exit
         endif
      endif
   enddo

   !
   ! Compress the data?
   !
   if ( info%ignore_whitespace ) then
      call xml_compress_( data, no_data )
   endif

end subroutine xml_get

subroutine xml_put( info, tag, endtag, attribs, no_attribs, &
                   data, no_data )
   type(XML_PARSE),  intent(inout)               :: info
   character(len=*), intent(in)                  :: tag
   logical,          intent(in)                  :: endtag
   character(len=*), intent(in), dimension(:,:)  :: attribs
   integer,          intent(in)                  :: no_attribs
   character(len=*), intent(in), dimension(:)    :: data
   integer,          intent(in)                  :: no_data

   integer         :: i

   character(len=300), parameter :: indent = ' '

   !
   ! Start by writing the tag
   !
   write( info%lun, '(3a)', advance = 'no' ) &
      indent(1:3*info%level), '<', adjustl(tag)
   write( info%lun, '(5a)', advance = 'no' ) &
      (' ',trim(attribs(1,i)),'="',trim(attribs(2,i)),'"' ,i=1,no_attribs )

   !
   ! The logic here should be more complex!
   !
   if ( .not. endtag ) then
      write( info%lun, '(a)' ) '>'
      info%level = info%level + 1
   else
      if ( no_data .eq. 0 ) then
         write( info%lun, '(a)' ) '/>'
      else
         write( info%lun, '(a)' ) '>'
      endif
   endif

   if ( no_data .gt. 0 ) then
      write( info%lun, '(2a)' ) &
         ( indent(1:3*info%level), trim(data(i)), i=1,no_data )
      if ( endtag ) then
         write( info%lun, '(4a)' ) &
            indent(1:3*info%level),'</', tag, '>'
      endif
   endif
   if ( endtag ) then
      info%level = info%level - 1
   endif
end subroutine xml_put

subroutine xml_compress_( data, no_data )
   character(len=*), intent(inout), dimension(:)    :: data
   integer,          intent(inout)                  :: no_data

   integer :: i
   integer :: j

   !
   ! TODO: keep empty lines in the middle!
   !
   j = 0
   do i = 1,no_data
      if ( len_trim(data(i)) .ne. 0 ) then
         j       = j + 1
         data(j) = adjustl(data(i))
      endif
   enddo

   no_data = j

end subroutine xml_compress_

subroutine xml_options( info, ignore_whitespace, no_data_truncation )
   type(XML_PARSE),  intent(inout)               :: info
   logical, intent(in), optional                 :: ignore_whitespace
   logical, intent(in), optional                 :: no_data_truncation

   if ( present(ignore_whitespace) ) then
      info%ignore_whitespace = ignore_whitespace
   endif
   if ( present(no_data_truncation) ) then
      info%no_data_truncation = no_data_truncation
   endif
end subroutine xml_options

logical function xml_ok( info )
   type(XML_PARSE),  intent(in)               :: info

   xml_ok = info%eof .or. info%error .or. &
            ( info%no_data_truncation .and.    &
                 ( info%too_many_attribs .or. info%too_many_data ) )
   xml_ok = .not. xml_ok
end function xml_ok

logical function xml_data_trunc( info )
   type(XML_PARSE),  intent(in)               :: info

   xml_data_trunc = info%too_many_attribs .or. info%too_many_data
end function xml_data_trunc

integer function xml_find_attrib( attribs, no_attribs, name, value )
   character(len=*), dimension(:,:)  :: attribs
   integer                           :: no_attribs
   character(len=*)                  :: name
   character(len=*)                  :: value

   integer :: i

   xml_find_attrib = -1
   do i = 1,no_attribs
      if ( name .eq. attribs(1,i) ) then
         value           = attribs(2,i)
         xml_find_attrib = i
         exit
      endif
   enddo

end function xml_find_attrib

end module xmlparse
