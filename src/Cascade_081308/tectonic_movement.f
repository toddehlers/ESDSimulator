c tectonic_movement

      subroutine tectonic_movement (x,y,nnode,
     &                              fix,dt,itime,
     &                              points,vertices,neighbour,
     &                              surface,newsurface,nt,
     &                              nn,nb,nbmax,nswaps,
     &                              mask,mask_e,xy,pp,aa,bb,surfscale,
     &                              delta,tcheck,sidex,sidey,
     &                              advec_vel,
     &                              nn2,nb2,cell)

c this routine allows the points of the grid to be moved around
c due to horizontal tectonic movements

c INPUT: x, y      = x- and y-nodal coordinates
c        nnode     = number of nodes
c        fix       = boundary conditions
c        dt        = time step length
c        itime     = time step
c        points    = working array
c        vertices  = connectivity list
c        neighbour = triangle neighbour list
c        surface   = surface attached to each node
c        newsurface= to determine whether nodal surface has to be reestimated
c        nt        = number of triangles
c        nn        = neighbour list
c        nb        = number of neighbour per node
c        nbmax     = maximum number of neighbours
c        mask      = working array
c        mask_e    = working array
c        xy        = working array
c        pp        = working array
c        aa        = working array
c        bb        = working array
c
c DS 7/27/1
c        delta     = reference grid spacing [km]
c        tcheck    = time since last check [yrs]
c        sidex     = length in x direction [km]
c        sidey     = length in y direction [km]
c        advec_vel = maximum advection velocity [km/yr]

c OUTPUT: x,y       = updated nodal coordinates
c         nt        =number of triangles
c        vertices  = connectivity list
c        neighbour = triangle neighbour list
c        surface   = surface attached to each node
c        nn        = neighbour list
c        nb        = number of neighbour per node
c        nswaps    = number of triangle side swaps operated during the
c                    grid reorganisation

c subroutines called:
c - debug
c - del_flip
c - find_surface

      common /vocal/ ivocal

      real*4    x(nnode),y(nnode)
      real*8    points(2,*)
      real      fix(nnode)
      integer   vertices(3,*),neighbour(3,*)
      real*4    surface(nnode)
      real      newsurface(nnode)
      integer   nb(nnode),nn(nbmax,nnode)
      logical   mask(*),mask_e(3,*)
      real      xy(2),pp(2,nbmax),aa(nbmax,2),bb(nbmax)

      integer   nb2(nnode),nn2(nbmax,nnode),cell(nnode,nbmax,2)

      integer   nn3(nbmax,nnode)

c This is pure shear shortening in Y with velocity (advec_vel) D.S 7/27/1
      do i=1,nnode
       if(fix(i).gt.0.5)then
        y(i)=y(i)-advec_vel*(y(i)/sidey)*dt
       endif
      enddo

c this is the spinning landscape
c      omega=2.*3.1415/1.e6
c        do i=1,nnode
c        x0=x(i)-50.
c        y0=y(i)-50.
c          if (x0**2+y0**2.lt.25.**2) then
c          x(i)=x(i)-omega*dt*y0
c          y(i)=y(i)+omega*dt*x0
c          endif
c        enddo

c perform a check to see if mesh needs to be updated DS 11/27/1 
c  remember endif statement at end if commenting out
      if (advec_vel*tcheck.gt.delta/10.) then
 
c reset tcheck
      tcheck = 0.

c from here on you should not change anything...
      do i=1,nnode
       points(1,i)=x(i)
       points(2,i)=y(i)
      enddo

      do i=1,nnode
       newsurface(i)=0.
      enddo

      if (ivocal.eq.1) call debug ('del_flip$',0)
      call del_flip (points,neighbour,vertices,nt,mask,mask_e,nswaps)
        do it=1,nt
        if (mask(it)) then
        newsurface(vertices(1,it))=1.
        newsurface(vertices(2,it))=1.
        newsurface(vertices(3,it))=1.
        endif
        enddo
      if (ivocal.eq.1) call debug ('tectonic_movement$',1)

c if there were any swapping the natural neighbours and surfaces need to be 
c recalculated...

      if (nswaps.ne.0) then

        do i=1,nnode
        nb(i)  = 0
        nb2(i) = 0
        enddo
 
c
c find nn,nb, nb2,nn2,cell
c
        do it=1,nt
        i1=vertices(1,it)
        i2=vertices(2,it)
        i3=vertices(3,it)
 
c i1
        nb(i1)=nb(i1)+1
        nn(nb(i1),i1)=i2
        if (nb(i1).gt.nbmax) stop 'nbmax too small...'
 
        nb2(i1)=nb2(i1)+1
        nn2(nb2(i1),i1)=i2
        cell(i1,int((1.+real(nb2(i1)))/2.),1)=i2
 
        nb2(i1)=nb2(i1)+1
        nn2(nb2(i1),i1)=i3
        cell(i1,int((1.+real(nb2(i1)))/2.),2)=i3
        if (nb2(i1).gt.nbmax) stop 'nbmax too small...1'
 
c i2
        nb(i2)=nb(i2)+1
        nn(nb(i2),i2)=i3
        if (nb(i2).gt.nbmax) stop 'nbmax too small...'
 
        nb2(i2)=nb2(i2)+1
        nn2(nb2(i2),i2)=i3
        cell(i2,int((1.+real(nb2(i2)))/2.),1)=i1
 
        nb2(i2)=nb2(i2)+1
        nn2(nb2(i2),i2)=i1
        cell(i2,int((1.+real(nb2(i2)))/2.),2)=i3
        if (nb2(i2).gt.nbmax) stop 'nbmax too small...2'
 
c i3
        nb(i3)=nb(i3)+1
        nn(nb(i3),i3)=i1
        if (nb(i3).gt.nbmax) stop 'nbmax too small...'
 
        nb2(i3)=nb2(i3)+1
        nn2(nb2(i3),i3)=i1
        cell(i3,int((1.+real(nb2(i3)))/2.),1)=i1
 
        nb2(i3)=nb2(i3)+1
        nn2(nb2(i3),i3)=i2
        cell(i3,int((1.+real(nb2(i3)))/2.),2)=i2
        if (nb2(i3).gt.nbmax) stop 'nbmax too small...3'

        enddo
 
c
c added directly from find_neighbours.f
c
        do i=1,nnode
 
c classement de nn2 de i le plus gd au plus pt
         do ie=1,nb2(i)-1
          ja=nn2(ie,i)
          jb=ie
          do jess=ie+1,nb2(i)
           jess2=nn2(jess,i)
           if(ja.lt.jess2) then
            ja=jess2
            jb=jess
           endif
          enddo
          nn2(jb,i)=nn2(ie,i)
          nn2(ie,i)=ja
         enddo
 
c elimination des pts doubles de enl
         i0=1
         nn3(1,i)=nn2(1,i)
         do iess=2,nb2(i)-1
          if(nn2(iess,i).ne.nn3(i0,i)) then
           nn3(i0+1,i)=nn2(iess,i)
           i0=i0+1
          endif
         enddo
         do j=1,i0
          nn2(j,i)=nn3(j,i)
         enddo
         nb2(i)=i0
        enddo
c
c back to regular code
c

      if (ivocal.eq.1) call debug ('find_surface$',0)
      call find_surface (nn,nb,surface,nbmax,nnode,
     &                   x,y,xy,pp,aa,bb,newsurface,surfscale)
      if (ivocal.eq.1) call debug ('tectonic_movement$',1)

        do i=1,nnode
        nb(i)=nb(i)+1
          if (nb(i).gt.nbmax) then
          stop 'nbmax too small...'
          endif
        nn(nb(i),i)=i
        enddo

      endif

      endif

      return
      end
