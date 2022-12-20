c debug

      subroutine debug (message,iflag)

c subroutine to help debug cascade
c it keeps track of the subroutine cascade is in when the program
c is stopped (by the user or by a run-time error)
c to use the debugger you must set ivocal to 1 in cascade
c the information is stored in a file called debug.out

c INPUT: message       = character string containing the message to be 
c                        saved in debug.out; the last character must
c                        be a dollar sign ($)
c        iflag         = if = 1 the message is added to debug.out
c                        if = 0 debug.out is first rewound

c subroutines called:
c NONE

      common /vocal/ ivocal

      character message*256

      nmessage=0
        do i=1,256
          if (message(i:i).eq.'$') then
          nmessage=i-1
          goto 1
          endif
        enddo
       
1      continue
 
      if (iflag.eq.0) then
      rewind (89)
      write (89,'(a)') 'Program stopped after entering routine '
     &                 //message(1:nmessage)
      else
      write (89,'(a)') 'and after returning in '//message(1:nmessage)
      endif
      call flush(89)

      return
      end
