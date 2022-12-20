      subroutine random (x,n)
      common /seed/idum

c subroutine used to generate random numbers (derived from
c numerical recipes; see book for further information)

c INPUT: n    = number of random numbers to be generated

c OUTPUT: x   = array of random numbers of length n

c subroutines called:
c NONE

      real    x(n)
c  idum is the seed 
c
          do 100 i=1,n
          x(i)=ran0(idum)
100       continue
      return
      end
      function ran0(idum)
      integer idum,ia,im,iq,ir,mask
      real ran0,am
      parameter (ia=16807,im=2147483647,am=1./im,
     *iq=127773,ir=2836,mask=123459876)
      integer k
      idum=ieor(idum,mask)
      k=idum/iq
      idum=ia*(idum-k*iq)-ir*k
      if(idum.lt.0) idum=idum+im
      ran0=am*idum
      idum=ieor(idum,mask)
      return
      end

