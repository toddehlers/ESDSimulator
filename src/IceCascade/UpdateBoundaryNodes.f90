! tectonic_uplift
        module m_UpdateBoundaryNodes
          contains
          subroutine UpdateBoundaryNodes(configData)
                use rt_param
                use cascade_globals
                implicit none

!      subroutine tectonic_uplift(x,y,h,h0,hi,nnode,fix,dt,time,influx,surface,uplift_rate,shelf)

! this routine defines the tectonic uplift function
! it may vary in space and time

! INPUT: x,y       = x- and y-nodal coordinates (in km's)
!        h         = present topography (in m's)
!        h0        = bedrock-alluvials interface (in m's)
!        hi        = original topograpy (in m's)
!        nnode     = number of nodes
!        fix       = boundary conditions
!        dt        = time step length
!        time      = current time
!        surface   = surface associated with each node
!        uplift_rate = rate in m/yr - defined w/ general params
!       shelf       =array controlling where uplift occurs
! OUTPUT:  h        = update topography
!          hi       = updated original topograpy
!          h0       = updated bedrock-alluvials interface
!          influx   = updated influx of material by tectonic uplift

! subroutines called:
! NONE

!      real     x(nnode),y(nnode),h(nnode),h0(nnode),hi(nnode)
!      real     fix(nnode),shelf(nnode)
!      real     surface(nnode),influx
!      real     uplift_rate

!      real     dij

  type(config) :: configData
  integer(4) :: inode,count
  real(8) :: dh




! the example is a square topography that stops growing after a while
! note that you can use the array fix to prevent the base level to move
! with tectonic uplift

! commented out 6/15/1, now in initialize_general_parameters.f
!      uplift_rate=2.E-3
  !print*,"Rear Boundary Uplifting"
  !print*,"Max Elevation [m]:",configData%imposeRearBoundaryUpliftMaxElevation*1e3_8
  !print*,"Uplift Rate [m/yr]:",configData%imposeRearBoundaryUpliftRate / 1e3_8
  !print*,"dt:",global_cascade_dt
  count=0

  do inode=1,configData%nnode

! uplift the topography keeping the boundary pinned
!     if (memory(inode, 5) > 0.5_8) then

!       dh=configData%uplift_rate*global_cascade_dt*shelf(inode)

! ! uplift interior 
!       h(inode) = h(inode) + (dh*shelf(inode))
!        if (shelf(inode) > 0.0 ) then
!            print *, "h(inode), dh, uplift_rate, dt, shelf: ", &
!            h(inode), dh, configData%uplift_rate, global_cascade_dt, shelf(inode)
!        endif

!       print*,dh,shelf(inode),h(inode),dt,uplift_rate
! added 6/15/1
!  rectangular uplift pattern
!        if (x(inode).gt.3.0.and.x(inode).lt.30.-3.0.and.
!     &      y(inode).gt.3.0.and.y(inode).lt.30.-3.0) then
!         h(inode)=h(inode)+dh
!        endif
!

!  rectangular uplift pattern
!      if (x(inode).gt.2.0.and.x(inode).lt.8.0) then
!       h(inode) = h(inode) + dh
!      endif

! changes added 6/14/1
!     circular uplift pattern
!     dij = sqrt((x(inode)-25.)**2+(y(inode)-25.)**2)
!      if (dij.lt.10.) then
!       h(inode)=h(inode)+dh
!      endif

    ! x and y are the x- and y-coordinates of the nodes in km
! h is the current topography
! hi is the initial topography
! h0 is the location of the bedrock interface
! all h's are in m, iceth(nnodemax)

    !configData%imposeRearBoundaryUpliftMaxElevation is in km (CONVERT TO m's)
    !configData%imposeRearBoundaryUpliftRate is in mm/yr (or km/My) (CONVERT TO m's)

    ! Elseif condition Added by Victoria M Buford Parks Jan 2020
    ! For nodes that are at the rear boundary (eg y=200)
        if (y(inode).eq.configData%sidey) then
        !    print*, "inode:",inode,"y:",int(y(inode)),"h:",h(inode)
            if (h(inode).lt.configData%imposeRearBoundaryUpliftMaxElevation*1e3_8) then
                dh = configData%imposeRearBoundaryUpliftRate / 1e3_8 * global_cascade_dt
                h(inode)=h(inode)+dh
                
                if (h(inode).gt.configData%imposeRearBoundaryUpliftMaxElevation*1e3_8) then
                    h(inode)=configData%imposeRearBoundaryUpliftMaxElevation*1e3_8
                endif
                
              !  print*, "inode:",inode,"y:",y(inode),"new h:",h(inode),"dh",dh
                
            endif
        count=count+1
        ! Return to original code VMBP Jan 2020  
        ! else
        !     dh = 0.0_8
        
    ! do not touch the following lines
    ! they update h0, hi and calculate the influx of material into the
    ! landscape
        h0(inode)=h0(inode)+(dh)
        hi(inode)=hi(inode)+(dh)
        influx=influx+dh*memory(inode,7)
        endif

    enddo
    !print*,"# nodes at rear boundary", count
  return
  end subroutine UpdateBoundaryNodes
end module m_UpdateBoundaryNodes
