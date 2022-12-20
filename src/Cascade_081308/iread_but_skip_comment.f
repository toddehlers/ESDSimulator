c iread_but_skip_comment

      subroutine iread_but_skip_comment (iunit,n,x)

c to read in an array of integer*4 numbers of length n
c while skipping potential comment lines starting with the character '#'

c INPUT: n             = length of array to be read
c        iunit         = unit number where data is read from

c OUTPUT: x(n)         = array of length n read from unit iunit
 
c subroutines called:
c NONE

      common /vocal/ ivocal

      character*1 comment
      integer x(n)
 
    1 read (iunit,'(a1)',end=999,err=998) comment
      if (comment.eq.'#') goto 1
      backspace (iunit)
      read (iunit,*,end=999,err=998) x
 
      return

999   print*,'End of file reached on unit ',iunit
      stop '(iread_but_skip_comment)'

998   print*,'Wrong data type on unit ',iunit
      stop '(iread_but_skip_comment)'

      end
