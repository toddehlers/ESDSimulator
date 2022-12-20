      subroutine check_for_removal (surface,dh,
     &                              nb,nn,fluxmax_erosion,surfmin,dt,
     &                              iremove)

c subroutine where the condition for node removal is defined

c INPUT: surface     = surface attached to current node
c        dh          = incremental erosion over the current time step
c        nb          = number of neighbours per node
c        nn          = list of neighbours
c        fluxmax_erosion = maximum erosive flux permitted
c        surfmin     = minimum surface permitted
c        dt          = current time step length

c OUTPUT:
c        iremove     = 0 keep it
c                    = 1 remove it

c subroutines called:
c NONE

      integer nn(*)

      common /vocal/ ivocal

      iremove=0
      flux=min(0.,surface*dh/dt)
c        do j=1,nb
c        jj=nn(j)
c        flux=flux+min(0.,surface*dh/dt)
c        enddo

c original code:
c      if(flux.gt.-fluxmax_erosion/10. .or.
c     &   surface.lt.surfmin/2.) iremove=1

      if (surface.lt.surfmin/4.) iremove=1
      return
      end
