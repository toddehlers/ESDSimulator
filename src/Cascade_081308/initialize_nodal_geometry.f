c initialize_nodal_geometry

      subroutine initialize_nodal_geometry
     & (nnodemax,nnode,x,y,h,fix,delta,surfscale,
     &  run_name,nrun_name,surfmin,nx,ny,sidex,sidey,bdry)

c subroutine to read in initial nodal information

c INPUT: nnodemax      = maximum number of nodes that can be read
c        run_name      = name of the subdirectory where the run output
c                        will be stored
c        nrun_name     = length of the character string run_name

c OUTPUT: nnode        = number of nodes read
c         x(nnode)     = x-location of nodes (in km)
c         y(nnode)     = y-location of nodes (in km)
c         h(nnode)     = height of nodes (in m)
c         fix(nnode)   = boundary condition (=0 means node height is fixed;
c                        =1 means node height is free)
c         delta        = mean nodal spacing (in km)
c         surfscale    = mean nodal surface (a nodal surface is the surface
c                        of the part of the landscape associated with a node)
c                        (in km**2)
c         surfmin
c         bdry         = flag to identify moving bdry nodes (added 11/28/1 DS)
c                        0: not on bdry, 1: on bdry


c subroutines called:
c - debug
c - random
c - iread_but_skip_comment
c - read_but_skip_comment

      common /vocal/ ivocal

      real          x(*),y(*),h(*)
      real          fix(*)
      integer       bdry(*)
      character     run_name*256

      integer       nx,ny,nnode
      real          sidex,sidey
      real          pi

c # of valleys in initial topography, amplitude of initial valleys
      real          nv,A

      iread=0
      pi = 3.1415926
c Number of valleys requested
      nv = 1.
c Ridge-valley relief [m]
      A = 5000.

      if (iread.eq.0) then

c Nodal geometry now defined in initialize_general_parameters.f
c      nx=41
c      ny=41
c      nnode=nx*ny
c      sidex=40.
c      sidey=40.

      if (nnode.gt.nnodemax) then
      stop 'Too many nodes...'
      endif

      delta=sidex/float(nx-1)
      surfscale=sidex*sidey/nnode
      surfmin=surfscale/4

c this is a bit of random noise put on the initial grid
c note that the noise is not used for the nodes along the boundary

      if (ivocal.eq.1) call debug ('random$',0)
      call random (x,nnode)
      call random (y,nnode)
      call random (h,nnode)
      if (ivocal.eq.1) call debug ('initialize_nodal_geometry$',1)
        do j=1,ny
        do i=1,nx

        inode=(j-1)*nx+i

c set random topography perturbation flag
        ishake=1
        if (i.eq.1 .or.i.eq.nx .or.j.eq.1 .or.j.eq.ny) ishake=0

c for flat topography, multiply float() by 5 to decrease x,y noise
        x(inode)=((x(inode)-.5)/(float(nx-1))*ishake
     &            +1.*float(i-1)/float(nx-1))*sidex
        y(inode)=((y(inode)-.5)/(float(ny-1))*ishake
     &            +1.*float(j-1)/float(ny-1))*sidey
        fix(inode)=1.

c generate nv valleys (with zero elevation around entire model domain)
c modified by dwhipp 01/08
c parallel to x-axis
c      h(inode) =  h(inode) + (A/2)*cos(2.*pi*nv*y(inode)/sidey-pi)
c     & + (A/2)
c parallel to y-axis
c      h(inode) =  h(inode) + (A/2)*cos(2.*pi*nv*x(inode)/sidex-pi)
c     & + (A/2)

c pin the y=0, y=sidey sides of the model
c      if (j.eq.1.or.j.eq.ny)then
c        h(inode)=0.
c        fix(inode)=0.
c      endif

c set bdry flag
c      bdry(inode) = 0
c      if (i.eq.1.or.i.eq.nx) bdry(inode) = 1

c pin the x=0,x=sidex sides of the model
c      if (i.eq.1.or.i.eq.nx)then
c       h(inode)=0.0
c       fix(inode)=0.
c      endif

c pin all four sides of the model and also set the value of h(inode) to zero
c if it lies on the boundary
      if (i.eq.1 .or.i.eq.nx .or.j.eq.1.or.j.eq.ny) then
        fix(inode)=0.
        h(inode)=0.
      endif

c original code: initial "hat" topography
c        dh=y(inode)/sidey*2.*10.
c        if (y(inode).gt.sidey/2.) dh=200.-dh
c        h(inode)=h(inode)/1.e3+dh

c generate an inclined plane, approximate flank of shield volcano
c        theta=(5*3.141592654)/180
c        dh=y(inode)*tan(theta)
c        h(inode)=h(inode)/1.e3+dh
c convert h to meters from kilometers in keeping with Jean's format
c         h(inode)=h(inode)*1000

c generate a flat initial topography with noise
c         h(inode)=h(inode)/1.e3

c generate nv valleys (with sinusoidal boundaries)
c modified by dwhipp 01/08
c parallel to x-axis
c      h(inode) =  h(inode) + (A/2)*cos(2.*pi*nv*y(inode)/sidey-pi)
c     & + (A/2)
c parallel to y-axis
c      h(inode) =  h(inode) + (A/2)*cos(2.*pi*nv*x(inode)/sidex-pi)
c     & + (A/2)

c original code: just pin the x sides of the model
c      if (j.eq.1.or.j.eq.ny)then
c      h(inode)=0.0
c      fix(inode)=0.
c      endif

        enddo
        enddo

      else

      open (7,file=run_name(1:nrun_name)//
     &      '/cascade.node.in',status='old')

      call iread_but_skip_comment (7,1,nnode)
        if (nnode.gt.nnodemax) then
        stop 'Too many nodes...'
        endif

      call read_but_skip_comment (7,nnode,x)
      call read_but_skip_comment (7,nnode,y)
      call read_but_skip_comment (7,nnode,h)
      call read_but_skip_comment (7,nnode,fix)
      xmin=x(1)
      xmax=x(1)
      ymin=y(1)
      ymax=y(1)
        do i=1,nnode
        xmin=min(xmin,x(i))
        xmax=max(xmax,x(i))
        ymin=min(ymin,y(i))
        ymax=max(ymax,y(i))
        enddo
      sidex=xmax-xmin
      sidey=ymax-ymin
      snode=sqrt(float(nnode))
      delta=(sidex/snode+sidey/snode)/2.
      surfscale=sidex*sidey/nnode

      close (7)

      endif

        if (nnode.lt.3) then
        stop 'nnode too small...'
        endif

        if (delta.le.0.) then
        stop 'delta must be greater than 0...'
        endif

        if (surfscale.le.0.) then
        stop 'surfscale must be greater than 0...'
        endif

        if (sidex.le.0.) then
        stop 'sidex must be greater than 0...'
        endif

        if (sidey.le.0.) then
        stop 'sidey must be greater than 0...'
        endif

      return
      end
