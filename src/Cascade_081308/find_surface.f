c find_surface

      subroutine find_surface (nn,nb,surface,nbmax,nnode,
     &                         x,y,xx,pp,a,b,newsurface,surfscale)

c subroutine to calculate the surface area (part of the landscape)
c that is attached to each node.
c note that some voronoi cell are unbounded (because their node is
c on the convex hull) and others are abnomalously large (because)
c their node is within one "row" of being on the convex hull
c I have imbeded a simple test that checks for infinite voronoi cells
c and for voronoi cells that are bigger than 10 times the average.
c Both types of nodes are given a "mean" surface area.

c INPUT:  nn         = natural neighbour list
c         nb         = number of natural neighbours per node
c         nbmax      = maximum number of natural neighbours
c         nnode      = number of nodes
c         x,y        = x- and y- coordinates of nodes
c         cc         = working array
c         pp         = working array
c         a          = working array
c         b          = working array
c         newsurface = which nodal surface area has to be updated
c         surfscale  = mean nodal surface

c OUTPUT: surface   = surface area attahced to each node
c         newsurface is reset to 0.

c subroutines called:
c - debug
c - first_voronoi

      common /vocal/ ivocal

      real    x(nnode),y(nnode),surface(nnode)
      integer nb(nnode),nn(nbmax,nnode)
      real    xx(2),pp(2,nbmax),a(nbmax,2),b(nbmax)
      real    newsurface(nnode)

      n=2
      if (ivocal.eq.1) call debug ('first_voronoi$',0)
        do i=1,nnode
        if (newsurface(i).gt.0.5) then
        xx(1)=x(i)
        xx(2)=y(i)
          do j=1,nb(i)
          pp(1,j)=x(nn(j,i))
          pp(2,j)=y(nn(j,i))
          enddo
        m=nb(i)
        call first_voronoi (xx,pp,n,m,nbmax,2,a,b,vol)

c I have decided to remove the extra row of nodes around the perimeter
c of the real landscape; instead I check here for voronoi cells that
c have been given a volume of 0 (because they were unbounded) or
c for voronoi cells that are larger than 10 times the average value
c In both cases, the surface area attached to the node is set to some
c average value (surfscale=sidex*sidey/nnode)
c this modification may cause some unwanted results and needs to be
c fully tested
c alternatively I can come up with a better idea to remove the spurious
c voronoi cells that are generated for the nodes on the convex hull
c (infinite voronoi cells) or in its vicinity (very large voronoi cells)
          if (vol.gt.surfscale*10.) then
          surface(i)=surfscale
          elseif (vol.gt.0.) then
          surface(i)=vol
          else
          surface(i)=surfscale
          endif
        endif
        enddo
      if (ivocal.eq.1) call debug ('find_surface$',1)

      return
      end


