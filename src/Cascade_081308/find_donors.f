c find_donors

      subroutine find_donors (x,y,h,slope,length,nnode,
     &                        ndon,nb,nn,nbmax,delta)

c this subroutine simply finds the lowest neighbour of each point among
c its natural neighbours.
c It also calculates the length of the connection and the slope.
c Note that in the case one has made a periodic
c boundary condition, one needs to check the length of rivers across
c the whole grid.

c INPUT: x,y     = x- and y-nodal locations
c        h       = present topography
c        nnode   = number of nodes
c        nb      = number of neighbours per node
c        nn      = list of neighbours per node
c        nbmax   = maxiumum number of neighbours per node
c        delta   = mean nodal spacing (in km)

c OUTPUT: slope   = slope of stream segment originating at node i
c                    (downstream slope)
c         length  = length of stream segment originating at node i
c         ndon    = donor of node i (ie where the segment ends)

c subroutines called:
c NONE

      common /vocal/ ivocal

      integer ndon(nnode),nb(nnode),nn(nbmax,nnode)
      real    x(nnode),y(nnode),h(nnode)
      real    slope(nnode),length(nnode)

c loop to find donors

c loop over the nodes
        do inode=1,nnode
        dhmin=0.
        ndon(inode)=inode
        nbb=nb(inode)
c loop over the neighbours
          do jnode=1,nbb
          nnn=nn(jnode,inode)
          if (nnn.le.nnode .and. nnn.ne.inode) then
          dh=h(nnn)-h(inode)
          dx=x(nnn)-x(inode)
          dy=y(nnn)-y(inode)
          dd=sqrt(dx**2+dy**2)
            if (dh/dd.lt.dhmin) then
            ndon(inode)=nnn
            dhmin=dh/dd
            endif
          endif
          enddo
        enddo

c loop to calculate lengths and slopes

        do i=1,nnode
        if (ndon(i).eq.i) then
        slope(i)=0.
        length(i)=0.
        else
        dx=x(i)-x(ndon(i))
        dy=y(i)-y(ndon(i))
        length(i)=sqrt(dx**2+dy**2)
        slope(i)=(h(ndon(i))-h(i))/length(i)
        endif
        enddo

      return
      end
