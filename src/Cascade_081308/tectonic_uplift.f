c tectonic_uplift

      subroutine tectonic_uplift(x,y,h,h0,hi,nnode,
     &                           fix,dt,time,
     &                           influx,surface,
     &                           uplift_rate)

c this routine defines the tectonic uplift function
c it may vary in space and time

c INPUT: x,y       = x- and y-nodal coordinates
c        h         = present topography
c        h0        = bedrock-alluvials interface
c        hi        = original topograpy
c        nnode     = number of nodes
c        fix       = boundary conditions
c        dt        = time step length
c        time      = current time
c        surface   = surface associated with each node
c        uplift_rate = rate in m/yr - defined w/ general params

c OUTPUT:  h        = update topography
c          hi       = updated original topograpy
c          h0       = updated bedrock-alluvials interface
c          influx   = updated influx of material by tectonic uplift

c subroutines called:
c NONE

      common /vocal/ ivocal

      real     x(nnode),y(nnode),h(nnode),h0(nnode),hi(nnode)
      real     fix(nnode)
      real     surface(nnode),influx
      real     uplift_rate

      real     dij

c the example is a square topography that stops growing after a while
c note that you can use the array fix to prevent the base level to move
c with tectonic uplift

c commented out 6/15/1, now in initialize_general_parameters.f
c      uplift_rate=2.E-3

      do inode=1,nnode

c uplift the topography keeping the boundary pinned
       if (fix(inode).gt.0.5) then

        dh=uplift_rate*dt

c uplift interior 
        h(inode) = h(inode) + dh

c added 6/15/1
c  rectangular uplift pattern
c        if (x(inode).gt.3.0.and.x(inode).lt.30.-3.0.and.
c     &      y(inode).gt.3.0.and.y(inode).lt.30.-3.0) then
c         h(inode)=h(inode)+dh
c        endif
c

c  rectangular uplift pattern
c      if (x(inode).gt.2.0.and.x(inode).lt.8.0) then
c       h(inode) = h(inode) + dh
c      endif

c changes added 6/14/1
c     circular uplift pattern
c	  dij = sqrt((x(inode)-25.)**2+(y(inode)-25.)**2)
c	   if (dij.lt.10.) then
c       h(inode)=h(inode)+dh
c      endif
      
       endif

c do not touch the following lines
c they update h0, hi and calculate the influx of material into the
c landscape
        h0(inode)=h0(inode)+dh
        hi(inode)=hi(inode)+dh
        influx=influx+dh*surface(inode)

       enddo

      return
      end
