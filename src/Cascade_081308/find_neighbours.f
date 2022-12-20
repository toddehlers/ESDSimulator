c find_neighbours
c from peter van der beek (summmer,2001)

      subroutine find_neighbours (x,y,nn,nb,nnode,nbmax,nn2,nb2,
     &                            points,vertices,neighbour,nodes,
     &                            vis_tlist,vis_elist,add_tlist,nt,
     &                            surface,newsurface,eps,
     &                            xx,pp,aa,bb,surfscale,cell)

c subroutine to find the list of natural neighbours attached to each
c nodes. It computes the nn and nb lists but also the surfaces
c attached to the nodes.

c INPUT: x,y       = x,y nodal coordinates
c        nnode     = number of nodes
c        nbmax     = maximum number of neighbour per node
c        points    = working array
c        nodes     = working array
c        vis_tlist = working array
c        vis_elist = working array
c        add_tlist = working array
c        newsurface= here it is just used as a working array
c        eps       = precision
c        xx        = working array
c        pp        = working array
c        aa        = working array
c        bb        = working array
c        surfscale = average nodal surface

c OUTPUT: nn        = neighbour array
c         nb        = number of neighbour for each node
c         vertices  = triangle list
c         neighbour = neighbour list
c         nt        = number of triangle
c         surface   = nodal surface (voronoi cell surface area)
c		  nn2		= added by Dffait, unknown purpose (TAE 5/01)
c		  nb2		= added by Duffait, unknown purpose (TAE 5/01)
c		  cell		= added by Duffait, unknown purpose (TAE 5/01)

c subroutines called:
c NONE

      common /vocal/ ivocal

      real          x(*),y(*),surface(*)
      integer       nn(nbmax,*),nn2(nbmax,*)
      integer       nb(*),nb2(*)
      integer       cell(nnode,nbmax,2)

      real*8		points(2,*),eps
      integer		vertices(3,*)
      integer		neighbour(3,*)
      integer       nodes(*)
      integer       vis_tlist(*)
      integer       vis_elist(*)
      integer       add_tlist(*)
c      integer		ccw
c      integer           nn3(nbmax,*)

      real          xx(2),pp(2,nbmax),aa(nbmax,2),bb(nbmax)
      real          newsurface(*)

      integer       nn3(nbmax,nnode)
 
      do i=1,nnode
       points(1,i)=dble(x(i))
       points(2,i)=dble(y(i))
       nodes(i)=i
      enddo

      np=nnode

      nv_max=nnode

c sorting by x

      if (ivocal.eq.1) call debug ('indexx$',0)
      call indexx(np,points,nodes)
      if (ivocal.eq.1) call debug ('find_neighbours$',1)

c ensure initial triangle is in ccw order
 
      if(ccw(points(1,nodes(1)),
     &       points(1,nodes(2)),
     &       points(1,nodes(3)),k).eq.-1)then
             itemp = nodes(1)
             nodes(1) = nodes(2)
             nodes(2) = itemp
      end if

c						Call the routine that 
c						does the work

      mode=0
      eps=0.d0
      if (ivocal.eq.1) call debug ('delaun$',0)
      call delaun (points,np,neighbour,vertices,nt,2*np,
     &             vis_tlist,vis_elist,add_tlist,eps,nv_max,
     &             mode,inactive,bfirst,itstart,subset)
      if (ivocal.eq.1) call debug ('find_neighbours$',1)

      nt_max=nnode*3
        if(nt.gt.nt_max)then
        write(6,*)' Error number of triangles calculated is larger'
        write(6,*)' than array sizes.'
        write(6,*)' Number of triangles    = ',nt
        write(6,*)' Maximum value (nt_max) = ',nt_max
        stop
        endif

        do i=1,nnode
        nb(i)=0
        nb2(i)=0
        enddo

        do it=1,nt
        
c vertices are ordered ccw
         i1=vertices(1,it)
         i2=vertices(2,it)
         i3=vertices(3,it)
         
          nb(i1)=nb(i1)+1
          nn(nb(i1),i1)=i2

          nb2(i1)=nb2(i1)+1
          nn2(nb2(i1),i1)=i2

          if (nb2(i1).gt.nbmax) stop 'nbmax too small...'

          cell(i1,int((1.+real(nb2(i1)))/2.),1)=i2
          nb2(i1)=nb2(i1)+1
          nn2(nb2(i1),i1)=i3
          cell(i1,int((1.+real(nb2(i1)))/2.),2)=i3

          nb(i2)=nb(i2)+1
          nn(nb(i2),i2)=i3

          nb2(i2)=nb2(i2)+1
            if (nb2(i2).gt.nbmax) then
            stop 'nbmax too small...'
            endif

c testing
c          nn2(nb2(i2),i2)=i3
c          cell(i2,int((1.+real(nb2(i2)))/2.),1)=i3
c          nb2(i2)=nb2(i2)+1
c          nn2(nb2(i2),i2)=i1
c          cell(i2,int((1.+real(nb2(i2)))/2.),2)=i1

c original
          nn2(nb2(i2),i2)=i3
          cell(i2,int((1.+real(nb2(i2)))/2.),1)=i1
          nb2(i2)=nb2(i2)+1
          nn2(nb2(i2),i2)=i1
          cell(i2,int((1.+real(nb2(i2)))/2.),2)=i3

          nb(i3)=nb(i3)+1
          nn(nb(i3),i3)=i1

          nb2(i3)=nb2(i3)+1
            if (nb2(i3).gt.nbmax) then
            stop 'nbmax too small...'
            endif
          nn2(nb2(i3),i3)=i1
          cell(i3,int((1.+real(nb2(i3)))/2.),1)=i1
          nb2(i3)=nb2(i3)+1
          nn2(nb2(i3),i3)=i2
          cell(i3,int((1.+real(nb2(i3)))/2.),2)=i2

        enddo
        
c sort neighbor id's
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

        do i=1,nnode
        newsurface(i)=1.
        enddo

      if (ivocal.eq.1) call debug ('find_surface$',0)
      call find_surface (nn,nb,surface,nbmax,nnode,
     &                   x,y,xx,pp,aa,bb,newsurface,surfscale)
      if (ivocal.eq.1) call debug ('find_neighbours$',1)

c add id of node to end of neighbors list
        do i=1,nnode
        nb(i)=nb(i)+1
          if (nb(i).gt.nbmax) then
          stop 'nbmax too small...'
          endif
        nn(nb(i),i)=i
        enddo

      return
      end
