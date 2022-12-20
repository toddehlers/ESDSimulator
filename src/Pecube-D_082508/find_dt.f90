      subroutine find_dt (zl,diffusivity,nsurf,zsurf,zsurfp, &
                 nz,Pecletz,timesurf,timesurfp,istep,eps,ilog, &
                 dt,ntime,istatic)

      implicit real*8 (a-h,o-z)

      real*8 zsurf(nsurf),zsurfp(nsurf)

! first constraint on time step from conduction

      dt1=zl**2/diffusivity/100.

! second constrain on time step from advection

!      dt2=dt1
!      dzmin=zl
!        do i=1,nsurf
!        dzmin=min(dzmin,(1.-sqrt(1.-1./(nz-1)))*(zl+zsurfp(i)))
!        dzmin=min(dzmin,(1.-sqrt(1.-1./(nz-1)))*(zl+zsurf(i)))
!        enddo
!      if (abs(Pecletz).gt.eps) dt2=dzmin/abs(Pecletz)/2.

      dt2=dt1
      if (abs(Pecletz).gt.eps) dt2=zl/abs(Pecletz)/100.

! third constrain on time step from surface lowering

      dt3=dt1
        if (istep.ne.0) then
        dzmax=0.
          do i=1,nsurf
          dzmax=max(dzmax,zsurfp(i)-zsurf(i))
          enddo
        Pesurf=dzmax/(timesurf-timesurfp)
!        if (abs(Pesurf).gt.eps) dt3=dzmin/Pesurf/5.
        if (abs(Pesurf).gt.eps) dt3=zl/Pesurf/5.
        endif

! find optimum time step and number of steps

      dt=min(dt1,dt2)
      dt=min(dt,dt3)
      if (istep.ne.0) dt=min(dt,timesurf-timesurfp)

      ntime=int((timesurf-timesurfp)/dt)
      ntime=ntime+1
      dt=(timesurf-timesurfp)/ntime
      istatic=0

        if (istep.eq.0) then
        ntime=1
        dt=0.
        istatic=1
        endif

      if (ilog.eq.1) write (9,*) 'ntime/dt/dt1/dt2/dt3= ',ntime,dt,dt1,dt2,dt3

      return
      end
