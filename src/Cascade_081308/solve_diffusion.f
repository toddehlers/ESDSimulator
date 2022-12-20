c solve_diffusion

      subroutine solve_diffusion (a1,a2,ntri,nnode,b,h,hp,icon,fix,
     &                            jc,kc,nc,ncmax,diag,
     &                            sea_level,outflux,surface)

c this subroutine solves iteratively (using Gauss-Siedel + overrelaxation)
c an algebraic system. The matrix is stored EBE in condensed form (a1 and a2)
c the right hand-side vector (b) is formed at the begining of iterations
c as a2*h0
c the solution is returned in h
c fix stores the b.c. (fix=0 means fixed height; fix=1 means variable height)
c 
c jc(ncmax,i) is an array giving the list of triangles to which node i
c belongs; there is a maximum of ncmax triangles to which any node may belong;
c node i belongs to nc(i) triangles
c for each node, kc(ncmax,i) gives the number (1-3) that node i has in triangle
c jc(ncmax,i).

c INPUT: a1        = elemental conductivity matrix
c        a2        = elemental mass matrix
c        ntri      = number of triangles (elements)
c        nnode     = number of nodes
c        b         = right hand-side vector
c        hp        = working array
c        icon      = connectivity matric (nodes to elements)
c        fix       = boundary conditions
c        jc        = working array
c        kc        = working array
c        nc        = working array
c        ncmax     = maximum number of elements to which a node may be attached
c        diag      = working array
c        sea_level = sea level
c        surface   = surface associated to each node

c OUPUT: h         = updated topography
c        outflux   = updated flux through the boundaries

c subroutines called:
c NONE

      common /vocal/ ivocal

      real     a1(6,ntri),a2(6,ntri),b(nnode)
      integer  icon(3,ntri)
      real     fix(nnode)
      real     h(nnode),hp(nnode),diag(nnode)
      integer  jc(ncmax,nnode),kc(ncmax,nnode),nc(nnode)
      real     surface(nnode)

      beta=1.2
      tol=1.e-4
      itmax=1000

      if (ivocal.eq.1) open (19,file='interation.out',status='unknown')

      hmax=0.
        do i=1,nnode
        b(i)=0.
        hp(i)=h(i)
        diag(i)=1.
        hmax=amax1(hmax,abs(h(i)))
        enddo

        do it=1,ntri
        i1=icon(1,it)
        i2=icon(2,it)
        i3=icon(3,it)
        if (i1.le.nnode .and. i2.le.nnode .and. i3.le.nnode) then
        b(i1)=b(i1)+a2(1,it)*hp(i1)+a2(2,it)*hp(i2)+a2(4,it)*hp(i3)
        b(i2)=b(i2)+a2(2,it)*hp(i1)+a2(3,it)*hp(i2)+a2(5,it)*hp(i3)
        b(i3)=b(i3)+a2(4,it)*hp(i1)+a2(5,it)*hp(i2)+a2(6,it)*hp(i3)
        diag(i1)=diag(i1)+a1(1,it)
        diag(i2)=diag(i2)+a1(3,it)
        diag(i3)=diag(i3)+a1(6,it)
        endif
        enddo

      htol=max(hmax,1.)*tol

      ncalc=0
110   mcalc=0
      ncalc=ncalc+1
      if (ncalc.gt.itmax) goto 999
      dhmin=hmax
      dhmax=0.
        do i=1,nnode
        sup=b(i)
c loop thru all triangles connected to node i
          do it=1,nc(i)

c get id of triangle
          jt=jc(it,i)

c get position of node i in triangle
          k=kc(it,i)
          if (k.eq.1) then
          sup=sup-(a1(1,jt)*h(icon(1,jt))
     &            +a1(2,jt)*h(icon(2,jt))
     &            +a1(4,jt)*h(icon(3,jt)))
          elseif (k.eq.2) then
          sup=sup-(a1(2,jt)*h(icon(1,jt))
     &            +a1(3,jt)*h(icon(2,jt))
     &            +a1(5,jt)*h(icon(3,jt)))
          else
          sup=sup-(a1(4,jt)*h(icon(1,jt))
     &            +a1(5,jt)*h(icon(2,jt))
     &            +a1(6,jt)*h(icon(3,jt)))
          endif

        enddo

c this dh is the dh during the solution iteration
        dhcont=beta*sup/diag(i)
        dh=fix(i)*dhcont

        if (h(i).lt.sea_level) dh=0.

c enforce constant elevation boundary conditions
        h(i)=hp(i)+dh

        if (ivocal.eq.1) then
         dhmin=min(dhmin,abs(dh))
         dhmax=max(dhmax,abs(dh))
        endif

        outflux=outflux+(fix(i)-1.)*dhcont*surface(i)

        if (abs(dh).gt.abs(htol)) mcalc=mcalc+1

        enddo

c update height working array
        do i=1,nnode
        hp(i)=h(i)
        enddo

      if (ivocal.eq.1) write (19,*) ncalc,mcalc,dhmin,dhmax,htol

      if (mcalc.ne.0) goto 110

      if (ivocal.eq.1) close (19)

      return

999   if (ivocal.eq.1) print*,'Have a look in iteration.out'
      close (19)
      stop 'too many iterations...'

      end
