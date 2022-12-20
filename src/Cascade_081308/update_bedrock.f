c update_bedrock

      subroutine update_bedrock (h,h0,nnode)

c subroutine to update bedrock position in case of bedrock incision

c INPUT: h       = present topography
c        h0      = bedrock-alluvial interface
c        nnode   = number of nodes

c OUTPUT: h0      = updated bedrock-alluvial interface

c subroutines called:
c NONE

      common /vocal/ ivocal

      real     h(nnode),h0(nnode)

        do i=1,nnode
        if (h(i).lt.h0(i)) h0(i)=h(i)+1.e-6
        enddo

      return
      end
