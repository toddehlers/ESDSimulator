c find_order

      subroutine find_order (ibucket,ndon,fix,iorder,
     &                       nnode,norder)

c subroutine to calculate the proper order in which the nodes must be
c stored to perform the river erosion calculations

c INPUT: ibucket     = working array
c        ndon        = donor array
c        fix         = boundary condition array
c        nnode       = total number of nodes

c OUTPUT: iorder      = node ordering
c         norder      = number of nodes in the ordering

c subroutines called:
c NONE

      common /vocal/ ivocal

      integer ibucket(*),ndon(*),iorder(*)
      real    fix(*)
      logical more

      norder=0

c find self-donors
      do i=1,nnode
       ibucket(i)=1
       if (ndon(i).eq.i) ibucket(i)=0
      enddo

c beginning of pseudo time stepping

c neat
    1 continue

      do i=1,nnode
       if (ibucket(i).ne.0) ibucket(ndon(i))=-1
      enddo

      more=.false.
      do i=1,nnode
       ib=ibucket(i)
       if (ib.gt.0) then
        norder=norder+1
        iorder(norder)=i
        ib=0
        more=.true.
       endif
       ibucket(i)=ib**2
      enddo

      if (more) goto 1

c puts self-donors at the end of the list
      do i=1,nnode
       if (ndon(i).eq.i) then
        norder=norder+1
        iorder(norder)=i
       endif
      enddo

      return
      end
