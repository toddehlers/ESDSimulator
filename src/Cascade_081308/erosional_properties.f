c erosional_properties

      subroutine erosional_properties(x,y,h,h0,hi,nnode,
     &                                param,nparam,nnodemax,
     &                                xkf,xlf_BR,xkdiff)

c subroutine to determine erosional properties (or lithology)
c written by Peter
c (comments by Jean)

c the main ingredients are the present topography (h)
c the present position of the bedrock-alluvion interface (h0)
c and the initial topography (hi)
c for every node

c note that erosional properties are nodal values
c and include param(*,1)=xkf (fluvial constant)
c             param(*,2)=xlf_BR (fluvial length scale)
c             param(*,3)=xkdiff (diffusional constant)

c INPUT: x,y      = x-,y- coordinates of nodes
c        h        = present topography
c        h0       = bedrock-alluvion interface
c        hi       = initial topography (at time 0)
c        nnode    = number of nodes
c        param    = erosional parameters
c        nparam   = numper ob erosional parameters
c        nnodemax = maximum number of nodes

c OUTPUT: param is updated

c subroutines called:
c NONE

      common /vocal/ ivocal

      real     x(nnode),y(nnode),h(nnode)
      real     h0(nnode),hi(nnode)
      real     param(nnodemax,nparam)
      real     xkf,xlf_BR,xkdiff

        do inode=1,nnode
c changed xkf from 1e-5 
c        param(inode,1)=5.e-4
        param(inode,1)=xkf
c changed xlf from 10 
c        param(inode,2)=1000.
        param(inode,2)=xlf_BR
c changed xkdiff from 1e-6
c        param(inode,3)=1.e-5
        param(inode,3)=xkdiff
        enddo

      return
      end
