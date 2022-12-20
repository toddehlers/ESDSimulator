c diffusion_erosion

      subroutine diffusion_erosion (xkdiff,
     &                              x,y,h,hp,dh,nnode,ntri,ncmax,
     &                              icon,kcon,jcon,ncon,
     &                              a1,a2,b,fix,
     &                              diag,dt,dhmin,dhmax,
     &                              sea_level,outflux,surface)

c subroutine to calculate diffusion erosion

c INPUT: xkdiff     = nodal diffusivity
c        x, y       = x- and y- nodal coordinates
c        h          = topography
c        hp         = working array (initial topography)
c        dh         = total height removed/added
c        nnode      = number of nodes
c        ntri       = number of triangles (here finite elements)
c        ncmax      = maximum number of elements attached to a node
c        icon       = connectivity matrix
c        kcon       = list of elements attached to a node
c        jcon       = node number in the equivalent element
c        ncon       = number of elements attached to a node
c        a1         = elemental conductivity matrix
c        a2         = elemental mass matrix
c        b          = right hand side vector
c        fix        = boundary condition array
c        diag       = working array (diagonal of matrix to be solved
c                     by Gauss-siedel method)
c        dt         = time step length (in yrs)
c        dhmin      = minimum amount eroded/deposited
c        dhmax      = maximum amount eroded/deposited
c        sea_level  = sea-level in meters
c        outflux    = flux out of the system through base level
c        surface    = array containing the surface attached to each node

c OUTPUT: a range of arrays/variables are updated:
c        h, dh, dhmin, dhmax, outflux
c         working arrays are modified

c subroutines calles:
c - debug
c - build_a
c - solve_diffusion

      common /vocal/ ivocal

      real    x(nnode),y(nnode),h(nnode),hp(nnode),dh(nnode)
      real    xkdiff(nnode)
      integer icon(3,ntri)
      real    fix(nnode)
      integer kcon(ncmax,nnode),jcon(ncmax,nnode),ncon(nnode)
      real    a1(6,ntri),a2(6,ntri),b(nnode),diag(nnode)
      real    surface(nnode)

c stores initial height

        do i=1,nnode
        hp(i)=h(i)
        dh(i)=h(i)
        enddo

c finds inverse connectivities

        do i=1,nnode
        ncon(i)=0
        enddo

        do it=1,ntri
        if (icon(1,it).le.nnode .and. icon(2,it).le.nnode
     &.and. icon(3,it).le.nnode) then
          do i=1,3
          ic=icon(i,it)
          ncon(ic)=ncon(ic)+1
            if (ncon(ic).gt.ncmax) then
            print*,ic,ncon(ic)
            print*,(kcon(ii,ic),ii=1,ncon(ic)-1)
            stop 'too tight connectivity...'
            endif
          kcon(ncon(ic),ic)=it
          jcon(ncon(ic),ic)=i
          enddo
        endif
        enddo

c builds matrices and right-hand sides vectors

c        alpha      = time integration parameter (0 = explicit)
c                                                (1 = implicit)
      alpha=1.
      if (ivocal.eq.1) call debug ('build_a$',0)
        do it=1,ntri
        i1=icon(1,it)
        i2=icon(2,it)
        i3=icon(3,it)
        xk=(xkdiff(i1)+xkdiff(i2)+xkdiff(i3))*(dt/3.)
          if (i1.le.nnode.and.i2.le.nnode.and.i3.le.nnode) then
          call  build_a (x(i1),x(i2),x(i3),y(i1),y(i2),y(i3),
     &                   xk,alpha,a1(1,it),a2(1,it))
          endif
        enddo
      if (ivocal.eq.1) call debug ('diffusion_erosion$',1)

        if (ivocal.eq.1) call debug ('solve_diffusion$',0)
        call solve_diffusion
     &                    (a1,a2,ntri,nnode,b,h,hp,icon,fix,
     &                     kcon,jcon,ncon,ncmax,diag,
     &                     sea_level,outflux,surface)
        if (ivocal.eq.1) call debug ('diffusion_erosion$',1)

c dh contains original heights
        dhmin=h(1)-dh(1)
        dhmax=h(1)-dh(1)
        do i=1,nnode

c this lines don't enforce any bc's
         dh(i)=fix(i)*(h(i)-dh(i))

         dhmin=amin1(dhmin,dh(i))
         dhmax=amax1(dhmax,dh(i))
        enddo

      return
      end
