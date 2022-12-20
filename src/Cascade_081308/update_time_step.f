c update_time_step

      subroutine update_time_step (dt,dt0,dtold,xkf,nnode,iadjust,
     &                             delta,iaction,
     &                             dhminfluv,dhmaxfluv,
     &                             dhmindiff,dhmaxdiff,
     &                             time,shorttime,writetime,endtime,
     &                             itime,dtmin,dtmax)

c subroutine to update time step based on Peter's modifications

c INPUT: dt           = current time step length
c        dt0          = original time step length
c        dtold        = used to avoid missing a writing time step
c        xkf          = fluvial constant []
c        nnode        = number of nodes
c        iadjust      = if iadjust=0, dt=dt0 always
c        delta        = mean nodal spacing
c        iaction      = 0 initialization phase
c                       1 update pahse
c        dhminfluv    = minimum erosion by fluvial action
c        dhmaxfluv    = maximum erosion by fluvial action
c        dhmindiff    = minimum erosion by diffusion
c        dhmaxdiff    = maximum erosion by diffusion
c        time         = current time
c        shorttime    = used to avoid missing a writing time step
c        writetime    = time at which the next saving si required
c        endtime      = time at which this run terminates
c        itime        = current time step

c OUTPUT:  dt          = new time step length
 
c subroutines called:
c NONE

      common /vocal/ ivocal

      real xkf(*)

c initialization phase
      if (iaction.eq.0) then

        xkfmax=xkf(1)
        do i=1,nnode
          xkfmax=max(xkf(i),xkfmax)
        enddo
        dtmax=100.*(0.5*delta)/xkfmax
        dtmin=dtmax/1.e4
        dt=dt0
        dtold=0.
        time=0.
        shorttime=0.
        itime=1

c update phase
      else

        if ((dhminfluv.lt.-1.e-3).or.(dhmindiff.lt.-1e-3)) then
        dt=dt/2.
        else
          if ((dhmaxfluv.lt.1.e-5).and.(dhmaxdiff.lt.1e-5)) dt=dt*2.
        endif

        if ((dhmaxfluv.gt.1.e-3).or.(dhmaxdiff.gt.1e-3)) then
        dt=dt/2.
        else
          if ((dhminfluv.gt.-1.e-5).and.(dhmindiff.gt.-1e-5)) dt=dt*2.
        endif

c dt must be bound to some finite values otherwise you may run into trouble

        if (dt.gt.dtmax) then
        dt=dtmax
        endif

        if (dt.lt.dtmin) then
        dt=dtmin
        endif

c shortcircuit the dynamic time step adjustment if iadjust=0
        if (iadjust.eq.0) dt=dt0

c updates time step to make sure that we dont miss an output time
c or the end of the run

        if ((shorttime+dt).gt.writetime) then
        dtold=dt
        dt=writetime-shorttime
        endif

        if ((time+dt).gt.endtime) dt=endtime-time

        itime=itime+1

      endif

      return
      end
