c
c------------------------------------------------------------------------
c
c	nn2Dfd - calculates natural neighbour interpolation and 
c		 1st derivatives at point x when point is inside convex hull 
c		 (using method based on a closed 2-D formula)
c
c	Input:
c		x(2)			co-ordinates of input point
c		points(2,np)		array of node points	
c	        vertices(3,nt)		array of triangle vertices	
c	        neighbour(3,nt)		array of neighbouring triangles.	
c					Neighbour(i,j) is the triangle
c					opposite node i in triangle j,
c					stored counterclockwise about j.
c               centres(3,nt)           centres(j,i) (j=1,2) contains the
c                                       co-ordinates of the centre of 
c                                       circumcircle about Delaunay 
c                                       triangle i, (i=1,...,nt),
c                                       j=3 contains squared radius of circle.
c	        loc			index of triangle containing 
c					input point.
c		data			data values at each node
c		nnpn_max		maximum number of neighbours per node
c					(depends on the point distribution,
c                                        set to ~50 in calling program)
c		tstore			integer work array of size nnpn_max
c		index			integer work array of size nnpn_max
c		vpairs			integer work array of size 2*nnpn_max
c		angles			real work array of size nnpn_max
c		pverts			real*8 work array of size 2*nnpn_max
c		vverts			real*8 work array of size 2*nnpn_max
c		dx_verts		real*8 work array of size 2*nnpn_max
c		dy_verts		real*8 work array of size 2*nnpn_max
c		ltri			logical work array of size nt_max
c		lnode			logical work array of size nt_max
c
c	Output:
c		f			interpolated data value
c		df(1)			interpolated derivative df/dx
c		df(2)			interpolated derivative df/dy
c		v			area of voronoi cell about x
c
c	Comments:
c		 On input point x must be in triangle loc. 
c		 Note: The input triangle contains x and so it's 
c		 circumcircle must also contain x.
c
c		 In order to use the direct formula to
c		 calculate the nn co-ordinate of the input point with
c		 respect to its ith neighbour we must find the set of
c		 nodes that are both neighbours of i and x, where 
c		 the neighbours of i are determined before addition of x.
c		 This is done using a LIFO stack to avoid an 
c		 expensive global search over all triangles.
c
c		 It is assumed that the logical arrays ltri(nt_max) 
c		 and lnode are initialized to false on input.
c
c		 This version does not calculate derivatives.
c
c		 Calls are made to: plot_tc, stackpairinit,
c				    poppair,pushpair & stackpairempty 
c				    indexx & second_v_area_d.
c
c					M. Sambridge, RSES, May 1996.
c
c------------------------------------------------------------------------
c
	Subroutine nn2Dfd
     &  (x,points,vertices,neighbour,centres,loc,data,nnpn_max,
     &   pverts,vverts,tstore,vpairs,angles,index,dx_verts,dy_verts,
     &   ltri,lnode,f,df,v)

	real*8		points(2,*)
	real*8		centres(3,*)
	real*8		data(*)
	real*8		x(2),dist,dx,dy
	real*8		f,v,df(2),deriv(2),dvx,dvy
        real*8          pverts(2,nnpn_max)
        real*8          vverts(2,nnpn_max)
        real*8          dx_verts(2,nnpn_max)
        real*8          dy_verts(2,nnpn_max)
        real*8          dpverts(4,2)
        real            angles(nnpn_max)
	integer		vertices(3,*)
	integer		neighbour(3,*)
	integer		tri,pos,tcount
	integer		tstore(nnpn_max)
	integer		index(nnpn_max)
	integer		vpairs(2,nnpn_max)
	logical		nnwrite
	logical*1	ltri(*)
	logical*1	lnode(*)
c
	integer		cyc(3)
        data		cyc/2,3,1/

        

	common/nnswitches/nnwrite,lud

c					Initialize variables
	f = 0.d0
	df(1) = 0.d0
	df(2) = 0.d0
	v = 0.d0
	dvx = 0.d0
	dvy = 0.d0
        tcount = 1

	if(nnwrite)write(lud,*)loc
c					plot input circum-triangle 
c					and circum-circle

 	call plot_tc(loc,points,vertices,centres)

c					record circum-triangle 
        tstore(tcount) = loc
        ltri(loc) = .true.
c
c					Find all other circumcircles 
c					containing input point using a 
c					directed walk
c
c					initialize stack
	call stackpairinit


c					put input triangle's neighbouring 
c					triangles on LIFO stack
c					together with position of input
c					triangle in their neighbour list

	i = neighbour(1,loc)
	j = neighbour(2,loc)
	k = neighbour(3,loc)
	if(i.ne.0)then
	   if(neighbour(1,i).eq.loc)then
              call pushpair(i,1)
	   else if(neighbour(2,i).eq.loc)then
              call pushpair(i,2)
	   else if(neighbour(3,i).eq.loc)then
              call pushpair(i,3)
           end if
	end if
	if(j.ne.0)then
	   if(neighbour(1,j).eq.loc)then
              call pushpair(j,1)
	   else if(neighbour(2,j).eq.loc)then
              call pushpair(j,2)
	   else if(neighbour(3,j).eq.loc)then
              call pushpair(j,3)
           end if
	end if
	if(k.ne.0)then
	   if(neighbour(1,k).eq.loc)then
              call pushpair(k,1)
	   else if(neighbour(2,k).eq.loc)then
              call pushpair(k,2)
	   else if(neighbour(3,k).eq.loc)then
              call pushpair(k,3)
           end if
	end if

 10	call stackpairempty(k)
c					if stack empty then finish
	if(k.eq.1)go to 100
c					take triangle from stack
	call poppair(tri,pos)

c					test if x is in circumcircle
        dx = x(1)-centres(1,tri)
        dy = x(2)-centres(2,tri)
	dist = dx*dx + dy*dy 

 	if(dist.lt.centres(3,tri))then
           if(nnwrite)write(lud,*)tri

c					plot circum-triangle and circum-circle

 	   call plot_tc(tri,points,vertices,centres)

	   if(pos.eq.1)then
	      new1 = neighbour(2,tri)
	      new2 = neighbour(3,tri)
	   else if(pos.eq.2)then
	      new1 = neighbour(1,tri)
	      new2 = neighbour(3,tri)
	   else 
	      new1 = neighbour(1,tri)
	      new2 = neighbour(2,tri)
	   end if
c					record current circum-triangle
c					and pseudo angle to x
           tcount = tcount + 1
           tstore(tcount) = tri
           ltri(tri) = .true.

	   if(new1.ne.0)then
	      if(neighbour(1,new1).eq.tri)then
                 call pushpair(new1,1)
	      else if(neighbour(2,new1).eq.tri)then
                 call pushpair(new1,2)
	      else if(neighbour(3,new1).eq.tri)then
                 call pushpair(new1,3)
              end if
	   end if
	   if(new2.ne.0)then
	      if(neighbour(1,new2).eq.tri)then
                 call pushpair(new2,1)
	      else if(neighbour(2,new2).eq.tri)then
                 call pushpair(new2,2)
	      else if(neighbour(3,new2).eq.tri)then
                 call pushpair(new2,3)
              end if
	   end if
	end if
	go to 10
 
 100    continue
        call stackpairflush()

c						check size of tcount
        if(tcount+2.gt.nnpn_max)then
           write(*,*)' '
           write(*,*)' Error: work arrays in subroutine nn2Dfd',
     &               ' not big enough.'
           write(*,*)' Too many circum-triangles for this node'
           write(*,*)' Remedy: Increase size of parameter'
           write(*,*)'         nnpn_max in calling program'
           write(*,*)'         current value = ',nnpn_max
           write(*,*)'         required value >= ',tcount+2
	   stop
        end if
c						find list of vertices
c						of new voronoi cell about
c						x and store the pair of
c						nodes associated with each
c						vertex
c       write(*,*)' Circum-triangles of x:'
        nv = 0
        do 110 i=1,tcount
           it = tstore(i)
c          write(*,*)it
           do j=1,3
              node = vertices(j,it)
              j1 = cyc(j)
              j2 = cyc(j1)
              na = vertices(j1,it)
              nb = vertices(j2,it)
              iopp = neighbour(j,it)
              if(iopp.eq.0.or..not.ltri(iopp))then
                 nv = nv + 1
c
c						find circum-centre of
c						x, node, and nodeb 
c
                 call circum_d (x,points(1,na),points(1,nb),
     &                          vverts(1,nv),
     &                          dx_verts(1,nv),dy_verts(1,nv))

                 vpairs(1,nv) = na
                 vpairs(2,nv) = nb
                 angles(nv) = pangle(x,vverts(1,nv))
c
c                write(*,*)' na =',na,' nb=',nb
c                write(*,*)' dx_verts =',dx_verts(1,nv)
c                write(*,*)' dx_verts =',dx_verts(2,nv)
c                write(*,*)' dy_verts =',dy_verts(1,nv)
c                write(*,*)' dy_verts =',dy_verts(2,nv)
c
              end if
           end do
           
 110    continue
        
c						sort vertices 
c
        call indexx(nv,angles,index)
c						debug I/O
c       write(*,*)' Number of vertices = ',nv
c       write(*,*)' index '
c	do i=1,nv
c          write(*,*)' i',i,' index',index(i),' a',angles(i)
c       end do
c
c						find vertices of each 
c						second-order Voronoi cell
c
	ncount = 0
c						loop over neighbours of x
        do 120 it=1,tcount
	   tri = tstore(it)
           do 130 in = 1,3
              node = vertices(in,tri)
              if(lnode(node))go to 130
              lnode(node) = .true.
c             write(*,*)' Neighbour node =',node
c
c						find pair of `external' vertices
c						of second-order voronoi cell  
c
              do i=1,nv
                 k = index(i)
                 if(vpairs(1,k).eq.node.or.
     &              vpairs(2,k).eq.node)then
                    pverts(1,1) = vverts(1,k) - x(1)
                    pverts(2,1) = vverts(2,k) - x(2)
                    dpverts(1,1) = dx_verts(1,k)
                    dpverts(2,1) = dx_verts(2,k)
                    dpverts(3,1) = dy_verts(1,k)
                    dpverts(4,1) = dy_verts(2,k)
                    ip = i+1
                    im = i-1
                    if(i.eq.nv)ip=1
                    if(i.eq.1)im=nv
                    kp = index(ip)
                    km = index(im)
                    if(vpairs(1,kp).eq.node.or.
     &                 vpairs(2,kp).eq.node)then
                       pverts(1,2) = vverts(1,kp) - x(1)
                       pverts(2,2) = vverts(2,kp) - x(2)
                       dpverts(1,2) = dx_verts(1,kp)
                       dpverts(2,2) = dx_verts(2,kp)
                       dpverts(3,2) = dy_verts(1,kp)
                       dpverts(4,2) = dy_verts(2,kp)
                    else if(vpairs(1,km).eq.node.or.
     &                      vpairs(2,km).eq.node)then
                       pverts(1,1) = vverts(1,km) - x(1)
                       pverts(2,1) = vverts(2,km) - x(2)
                       pverts(1,2) = vverts(1,k) - x(1)
                       pverts(2,2) = vverts(2,k) - x(2)
                       dpverts(1,1) = dx_verts(1,km)
                       dpverts(2,1) = dx_verts(2,km)
                       dpverts(3,1) = dy_verts(1,km)
                       dpverts(4,1) = dy_verts(2,km)
                       dpverts(1,2) = dx_verts(1,k)
                       dpverts(2,2) = dx_verts(2,k)
                       dpverts(3,2) = dy_verts(1,k)
                       dpverts(4,2) = dy_verts(2,k)
                    else
                       write(*,*)' Error in nn2Dfd'
                       write(*,*)
     &                 ' Can not find neighbouring vertices' 
                       write(*,*)' node=',node,' i =',i
                    end if
c
                    go to 150
                 end if
              end do
 150          continue
              nverts = 2
c
c						find `internal' vertices of 
c						second-order voronoi cell  
	      do 140 jt = 1,tcount
                 jtri = tstore(jt)
                 do jn = 1,3
                    jnode = vertices(jn,jtri)
                    if(jnode.eq.node)then
                       nverts = nverts + 1
c
c						record circum-centre
c
                       pverts(1,nverts) = centres(1,jtri) - x(1)
                       pverts(2,nverts) = centres(2,jtri) - x(2)

c                      write(*,*)' found common triangle',jtri

                       go to 140
                    end if
                 end do
 140          continue
c						Calculate volume of 
c						second-order voronoi cell 
c						for current node using
c						direct 2-D formula 
c
              call second_v_area_d
     &        (x,nverts,pverts,dpverts,angles,area,deriv)
c
              v = v + area
              f = f + area*data(node)
              df(1) = df(1) + data(node)*deriv(1)
              df(2) = df(2) + data(node)*deriv(2)
              dvx = dvx + deriv(1)
              dvy = dvy + deriv(2)

 130       continue
 120    continue
 	if(v.ne.0.0)then
           f = f/v
           df(1) = (df(1) - f*dvx)/v
           df(2) = (df(2) - f*dvy)/v
        end if
        v = v*0.5
c
c						reset logical arrays
        do i=1,tcount
           tri = tstore(i)
           ltri(tri) = .false.
           do j = 1,3
              lnode(vertices(j,tri)) = .false.
           end do
        end do 

	return
	end
c
c------------------------------------------------------------------------
c
c	nn2Dfdd - calculates natural neighbour interpolation 
c		  and 1st & 2nd  derivatives at point x when 
c		  point is inside convex hull 
c		  (using method based on a closed 2-D formula)
c
c	Input:
c		x(2)			co-ordinates of input point
c		points(2,np)		array of node points	
c	        vertices(3,nt)		array of triangle vertices	
c	        neighbour(3,nt)		array of neighbouring triangles.	
c					Neighbour(i,j) is the triangle
c					opposite node i in triangle j,
c					stored counterclockwise about j.
c               centres(3,nt)           centres(j,i) (j=1,2) contains the
c                                       co-ordinates of the centre of 
c                                       circumcircle about Delaunay 
c                                       triangle i, (i=1,...,nt),
c                                       j=3 contains squared radius of circle.
c	        loc			index of triangle containing 
c					input point.
c		data			data values at each node
c		nnpn_max		maximum number of neighbours per node
c					(depends on the point distribution,
c                                        set to ~50 in calling program)
c		tstore			integer work array of size nnpn_max
c		index			integer work array of size nnpn_max
c		vpairs			integer work array of size 2*nnpn_max
c		angles			real work array of size nnpn_max
c		pverts			real*8 work array of size 2*nnpn_max
c		vverts			real*8 work array of size 2*nnpn_max
c		dx_verts		real*8 work array of size 2*nnpn_max
c		dy_verts		real*8 work array of size 2*nnpn_max
c		dxx_verts		real*8 work array of size 2*nnpn_max
c		dyy_verts		real*8 work array of size 2*nnpn_max
c		dxy_verts		real*8 work array of size 2*nnpn_max
c		ltri			logical work array of size nt_max
c		lnode			logical work array of size nt_max
c
c	Output:
c		f			interpolated data value
c		df(1)			interpolated derivative df/dx
c		df(2)			interpolated derivative df/dy
c		ddf(1)			interpolated derivative d2f/dxx
c		ddf(2)			interpolated derivative d2f/dyy
c		ddf(3)			interpolated derivative d2f/dxy
c		v			area of voronoi cell about x
c
c	Comments:
c		 On input point x must be in triangle loc. 
c		 Note: The input triangle contains x and so it's 
c		 circumcircle must also contain x.
c
c		 In order to use the direct formula to
c		 calculate the nn co-ordinate of the input point with
c		 respect to its ith neighbour we must find the set of
c		 nodes that are both neighbours of i and x, where 
c		 the neighbours of i are determined before addition of x.
c		 This is done using a LIFO stack to avoid an 
c		 expensive global search over all triangles.
c
c		 It is assumed that the logical arrays ltri(nt_max) 
c		 and lnode are initialized to false on input.
c
c		 This version does not calculate derivatives.
c
c		 Calls are made to: plot_tc, stackpairinit,
c				    poppair,pushpair & stackpairempty 
c				    indexx & second_v_area_dd.
c
c					M. Sambridge, RSES, May 1996.
c
c------------------------------------------------------------------------
c
	Subroutine nn2Dfdd
     &  (x,points,vertices,neighbour,centres,loc,data,nnpn_max,
     &   pverts,vverts,tstore,vpairs,angles,index,dx_verts,dy_verts,
     &   dxx_verts,dyy_verts,dxy_verts,ltri,lnode,f,df,ddf,v)

	real*8		points(2,*)
	real*8		centres(3,*)
	real*8		data(*)
	real*8		x(2),dist,dx,dy
	real*8		f,v,df(2),ddf(3)
        real*8		deriv(2),dderiv(3),dvx,dvy,dvxx,dvyy,dvxy
        real*8          pverts(2,nnpn_max)
        real*8          dpverts(4,2)
        real*8          ddpverts(6,2)
        real*8          vverts(2,nnpn_max)
        real*8          dx_verts(2,nnpn_max)
        real*8          dy_verts(2,nnpn_max)
        real*8          dxx_verts(2,nnpn_max)
        real*8          dyy_verts(2,nnpn_max)
        real*8          dxy_verts(2,nnpn_max)
	integer		vertices(3,*)
	integer		neighbour(3,*)
	integer		tri,pos,tcount
	integer		tstore(nnpn_max)
	integer		index(nnpn_max)
	integer		vpairs(2,nnpn_max)
	logical		nnwrite
	logical*1	ltri(*)
	logical*1	lnode(*)
        real            angles(nnpn_max)
	integer		cyc(3)
        data		cyc/2,3,1/

	common/nnswitches/nnwrite,lud
c					Initialize variables
	f = 0.d0
	df(1) = 0.d0
	df(2) = 0.d0
	ddf(1) = 0.d0
	ddf(2) = 0.d0
	ddf(3) = 0.d0
	v = 0.d0
	dvx = 0.d0
	dvy = 0.d0
	dvxx = 0.d0
	dvyy = 0.d0
	dvxy = 0.d0
        tcount = 1

	if(nnwrite)write(lud,*)loc
c					plot input circum-triangle 
c					and circum-circle

 	call plot_tc(loc,points,vertices,centres)

c					record circum-triangle 
        tstore(tcount) = loc
        ltri(loc) = .true.
c
c					Find all other circumcircles 
c					containing input point using a 
c					directed walk
c
c					initialize stack
	call stackpairinit


c					put input triangle's neighbouring 
c					triangles on LIFO stack
c					together with position of input
c					triangle in their neighbour list

	i = neighbour(1,loc)
	j = neighbour(2,loc)
	k = neighbour(3,loc)
	if(i.ne.0)then
	   if(neighbour(1,i).eq.loc)then
              call pushpair(i,1)
	   else if(neighbour(2,i).eq.loc)then
              call pushpair(i,2)
	   else if(neighbour(3,i).eq.loc)then
              call pushpair(i,3)
           end if
	end if
	if(j.ne.0)then
	   if(neighbour(1,j).eq.loc)then
              call pushpair(j,1)
	   else if(neighbour(2,j).eq.loc)then
              call pushpair(j,2)
	   else if(neighbour(3,j).eq.loc)then
              call pushpair(j,3)
           end if
	end if
	if(k.ne.0)then
	   if(neighbour(1,k).eq.loc)then
              call pushpair(k,1)
	   else if(neighbour(2,k).eq.loc)then
              call pushpair(k,2)
	   else if(neighbour(3,k).eq.loc)then
              call pushpair(k,3)
           end if
	end if

 10	call stackpairempty(k)
c					if stack empty then finish
	if(k.eq.1)go to 100
c					take triangle from stack
	call poppair(tri,pos)

c					test if x is in circumcircle
        dx = x(1)-centres(1,tri)
        dy = x(2)-centres(2,tri)
	dist = dx*dx + dy*dy 

 	if(dist.lt.centres(3,tri))then
           if(nnwrite)write(lud,*)tri

c					plot circum-triangle and circum-circle

 	   call plot_tc(tri,points,vertices,centres)

	   if(pos.eq.1)then
	      new1 = neighbour(2,tri)
	      new2 = neighbour(3,tri)
	   else if(pos.eq.2)then
	      new1 = neighbour(1,tri)
	      new2 = neighbour(3,tri)
	   else 
	      new1 = neighbour(1,tri)
	      new2 = neighbour(2,tri)
	   end if
c					record current circum-triangle
c					and pseudo angle to x
           tcount = tcount + 1
           tstore(tcount) = tri
           ltri(tri) = .true.

	   if(new1.ne.0)then
	      if(neighbour(1,new1).eq.tri)then
                 call pushpair(new1,1)
	      else if(neighbour(2,new1).eq.tri)then
                 call pushpair(new1,2)
	      else if(neighbour(3,new1).eq.tri)then
                 call pushpair(new1,3)
              end if
	   end if
	   if(new2.ne.0)then
	      if(neighbour(1,new2).eq.tri)then
                 call pushpair(new2,1)
	      else if(neighbour(2,new2).eq.tri)then
                 call pushpair(new2,2)
	      else if(neighbour(3,new2).eq.tri)then
                 call pushpair(new2,3)
              end if
	   end if
	end if
	go to 10
 
 100    continue
        call stackpairflush()

c						check size of tcount
        if(tcount+2.gt.nnpn_max)then
           write(*,*)' '
           write(*,*)' Error: work arrays in subroutine nn2Dfdd',
     &               ' not big enough.'
           write(*,*)' Too many circum-triangles for this node'
           write(*,*)' Remedy: Increase size of parameter'
           write(*,*)'         nnpn_max in calling program'
           write(*,*)'         current value = ',nnpn_max
           write(*,*)'         required value >= ',tcount+2
	   stop
        end if
c
c						find list of vertices
c						of new voronoi cell about
c						x and store the pair of
c						nodes associated with each
c						vertex
c       write(*,*)' Circum-triangles of x:'
        nv = 0
        do 110 i=1,tcount
           it = tstore(i)
c          write(*,*)it
           do j=1,3
              node = vertices(j,it)
              j1 = cyc(j)
              j2 = cyc(j1)
              na = vertices(j1,it)
              nb = vertices(j2,it)
              iopp = neighbour(j,it)
              if(iopp.eq.0.or..not.ltri(iopp))then
                 nv = nv + 1
c
c						find circum-centre of
c						x, node, and nodeb 
c
                 call circum_dd (x,points(1,na),points(1,nb),
     &                          vverts(1,nv),
     &                          dx_verts(1,nv),dy_verts(1,nv),
     &                          dxx_verts(1,nv),dyy_verts(1,nv),
     &                          dxy_verts(1,nv))

                 vpairs(1,nv) = na
                 vpairs(2,nv) = nb
                 angles(nv) = pangle(x,vverts(1,nv))
c                write(*,*)' na =',na,' nb=',nb
 
              end if
           end do
           
 110    continue
        
c						sort vertices 
c
        call indexx(nv,angles,index)
c						debug I/O
c       write(*,*)' Number of vertices = ',nv
c	do i=1,nv
c          write(*,*)' i',i,' index',index(i),' a',angles(i)
c       end do
c
c						find vertices of each 
c						second-order Voronoi cell
c
	ncount = 0
c						loop over neighbours of x
        do 120 it=1,tcount
	   tri = tstore(it)
           do 130 in = 1,3
              node = vertices(in,tri)
              if(lnode(node))go to 130
              lnode(node) = .true.
c             write(*,*)' Neighbour node =',node
c
c						find pair of `external' vertices
c						of second-order voronoi cell  
c
              do i=1,nv
                 k = index(i)
                 if(vpairs(1,k).eq.node.or.
     &              vpairs(2,k).eq.node)then
                    pverts(1,1) = vverts(1,k) - x(1)
                    pverts(2,1) = vverts(2,k) - x(2)
                    dpverts(1,1) = dx_verts(1,k)
                    dpverts(2,1) = dx_verts(2,k)
                    dpverts(3,1) = dy_verts(1,k)
                    dpverts(4,1) = dy_verts(2,k)
                    ddpverts(1,1) = dxx_verts(1,k)
                    ddpverts(2,1) = dxx_verts(2,k)
                    ddpverts(3,1) = dyy_verts(1,k)
                    ddpverts(4,1) = dyy_verts(2,k)
                    ddpverts(5,1) = dxy_verts(1,k)
                    ddpverts(6,1) = dxy_verts(2,k)
                    ip = i+1
                    im = i-1
                    if(i.eq.nv)ip=1
                    if(i.eq.1)im=nv
                    kp = index(ip)
                    km = index(im)
                    if(vpairs(1,kp).eq.node.or.
     &                 vpairs(2,kp).eq.node)then
                       pverts(1,2) = vverts(1,kp) - x(1)
                       pverts(2,2) = vverts(2,kp) - x(2)
                       dpverts(1,2) = dx_verts(1,kp)
                       dpverts(2,2) = dx_verts(2,kp)
                       dpverts(3,2) = dy_verts(1,kp)
                       dpverts(4,2) = dy_verts(2,kp)
                       ddpverts(1,2) = dxx_verts(1,kp)
                       ddpverts(2,2) = dxx_verts(2,kp)
                       ddpverts(3,2) = dyy_verts(1,kp)
                       ddpverts(4,2) = dyy_verts(2,kp)
                       ddpverts(5,2) = dxy_verts(1,kp)
                       ddpverts(6,2) = dxy_verts(2,kp)
                    else if(vpairs(1,km).eq.node.or.
     &                      vpairs(2,km).eq.node)then
                       pverts(1,1) = vverts(1,km) - x(1)
                       pverts(2,1) = vverts(2,km) - x(2)
                       pverts(1,2) = vverts(1,k) - x(1)
                       pverts(2,2) = vverts(2,k) - x(2)
                       dpverts(1,1) = dx_verts(1,km)
                       dpverts(2,1) = dx_verts(2,km)
                       dpverts(3,1) = dy_verts(1,km)
                       dpverts(4,1) = dy_verts(2,km)
                       dpverts(1,2) = dx_verts(1,k)
                       dpverts(2,2) = dx_verts(2,k)
                       dpverts(3,2) = dy_verts(1,k)
                       dpverts(4,2) = dy_verts(2,k)
                       ddpverts(1,1) = dxx_verts(1,km)
                       ddpverts(2,1) = dxx_verts(2,km)
                       ddpverts(3,1) = dyy_verts(1,km)
                       ddpverts(4,1) = dyy_verts(2,km)
                       ddpverts(5,1) = dxy_verts(1,km)
                       ddpverts(6,1) = dxy_verts(2,km)
                       ddpverts(1,2) = dxx_verts(1,k)
                       ddpverts(2,2) = dxx_verts(2,k)
                       ddpverts(3,2) = dyy_verts(1,k)
                       ddpverts(4,2) = dyy_verts(2,k)
                       ddpverts(5,2) = dxy_verts(1,k)
                       ddpverts(6,2) = dxy_verts(2,k)
                    else
                       write(*,*)' Error in nn2Dfdd'
                       write(*,*)
     &                 ' Can not find neighbouring vertices' 
                       write(*,*)' node=',node,' i =',i
                    end if
c
                    go to 150
                 end if
              end do
 150          continue
              nverts = 2
c
c						find `internal' vertices of 
c						second-order voronoi cell  
	      do 140 jt = 1,tcount
                 jtri = tstore(jt)
                 do jn = 1,3
                    jnode = vertices(jn,jtri)
                    if(jnode.eq.node)then
                       nverts = nverts + 1
c
c						record circum-centre
c
                       pverts(1,nverts) = centres(1,jtri) - x(1)
                       pverts(2,nverts) = centres(2,jtri) - x(2)

c                      write(*,*)' found common triangle',jtri

                       go to 140
                    end if
                 end do
 140          continue
c						Calculate volume of 
c						second-order voronoi cell 
c						for current node using
c						direct 2-D formula 
c
              call second_v_area_dd
     &             (x,nverts,pverts,dpverts,ddpverts,
     &              angles,area,deriv,dderiv)
c
              v = v + area
              f = f + area*data(node)
              df(1) = df(1) + data(node)*deriv(1)
              df(2) = df(2) + data(node)*deriv(2)
              dvx = dvx + deriv(1)
              dvy = dvy + deriv(2)
              dvxx = dvxx + dderiv(1)
              dvyy = dvyy + dderiv(2)
              dvxy = dvxy + dderiv(3)
              ddf(1) = ddf(1) + data(node)*dderiv(1)
              ddf(2) = ddf(2) + data(node)*dderiv(2)
              ddf(3) = ddf(3) + data(node)*dderiv(3)

c             write(*,*)' Second order Voronoi-cell between x and ',node,
c    &                  ' has ',nverts,' vertices'
 130       continue
 120    continue
 	if(v.ne.0.0)then
           f = f/v
           df(1) = (df(1) - f*dvx)/v
           df(2) = (df(2) - f*dvy)/v
           ddf(1) = (ddf(1) - 2.*dvx*df(1) - dvxx*f)/v
           ddf(2) = (ddf(2) - 2.*dvy*df(2) - dvyy*f)/v
           ddf(3) = (ddf(3) - dvy*df(1) - dvx*df(2) - dvxy*f)/v
        end if
        v = v*0.5
c
c						reset logical arrays
        do i=1,tcount
           tri = tstore(i)
           ltri(tri) = .false.
           do j = 1,3
              lnode(vertices(j,tri)) = .false.
           end do
        end do 

	return
	end
c
c------------------------------------------------------------------------
c
c	first_voronoi - calculates the area or volume of the voronoi	
c			polygon (polyhedron) about point x.
c
c	Input:
c		x(n)			input point	
c		p(n,m)			m neighbouring points about x
c		n			number of dimensions
c		m			number of neighbours
c		nmax			maximum number of dimensions
c		mmax			maximum number of neighbours
c		a(m,n)			work array 
c		b(m)			work array 
c
c	Output:
c		vol			volume (n=3) or area (n=2) of 
c					voronoi cell about x.
c
c	Comments:
c		 The first-order voronoi cell is defined by the perpendicular
c		 bisectors of x and each of its neighbours p(n,i).
c		 Each bisector forms a linear inequality constraint and 
c		 the voronoi cell is the region which satisfies all 
c		 constraints.
c
c		 This routine determines the matrix a(i,j) and vector b(j)   
c		 of the linear system and calls the recursive C routine
c		 `volume' to calculate the volume (area) of the voronoi cell.
c
c
c					J. Braun,
c					& M. Sambridge, RSES, 1995.
c
c------------------------------------------------------------------------
c
      Subroutine first_voronoi (x,p,n,m,mmax,nmax,a,b,vol)

      real  x(nmax),p(nmax,mmax)
      real  a(mmax,nmax),b(mmax)

c					Set up coefficients 
c					of linear system

      do 5 j=1,m
         do 6 i=1,n
            a(j,i)=p(i,j)-x(i)
 6       continue
         b(j)=0.
         do 7 i=1,n
            b(j)=b(j)+p(i,j)*p(i,j)-x(i)*x(i)
 7       continue
         b(j)=b(j)/2.
 5    continue

      call volume (a,b,m,n,mmax,nmax,vol)

      if(vol.eq.-1)then
c         write(*,*)' Warning in routine first_voronoi'
c         write(*,*)' Voronoi cell is unbounded'
c         write(*,*)' volume is set to zero'
         vol = 0.
      else if(vol.eq.0.)then
         write(*,*)' Warning in routine first_voronoi'
         write(*,*)' An inconsistent set of constraints has been'
         write(*,*)' detected in calculating the volume of '
         write(*,*)' a Voronoi cell.'
         write(*,*)' volume is set to zero'
      end if

      return
      end
c
c------------------------------------------------------------------------
c
c	second_voronoi - calculates the area or volume of the second-order
c			 voronoi polygon (polyhedron) between x and its
c			 neighbour p(n,1).
c
c	Input:
c		x(n)			input point	
c		p(n,m)			neighbours shared by x and p(n,1)
c		n			number of dimensions
c		m			number of neighbours
c		nmax			maximum number of dimensions
c		mmax			maximum number of neighbours
c		xo(n)			a local origin point.
c		a(m,n)			work array 
c		b(m)			work array 
c
c	Output:
c		vol			volume (n=3) or area (n=2) of 
c					voronoi cell about x.
c
c	Comments:
c		 The second-order voronoi cell between x and its neighbour p 
c		 is defined by the perpendicular bisector between x and p, 
c		 and the bisectors between p and the nodes which are neighbours 
c		 of both p and x. Each bisector forms a linear inequality 
c		 constraint and the second-order voronoi cell is the region 
c		 which satisfies all constraints.
c
c		 This routine determines the matrix a(i,j) and vector b(j)   
c		 of the linear system and calls the recursive C routine
c		 `volumef' to calculate the volume (area) of the voronoi cell.
c
c		 The point p is stored in the first row of input array p 
c		 (p(i,1),i=1,n). xo is any `nearby' point used as the origin
c		 for the calculation. It's use is only to help reduce
c		 roundoff error. For perfect arithmetic it will not affect  
c		 the calculated voronoi volume.
c
c					J. Braun, 
c					& M. Sambridge, RSES, 1995.
c
c------------------------------------------------------------------------
c
      Subroutine second_voronoi (x,p,n,m,mmax,nmax,xo,a,b,vol)

      real*8	x(nmax)
      real*8	p(nmax,mmax)
      real*8    xo(nmax)
      real  	a(mmax,nmax)
      real	b(mmax)
      real*8	dx1,dx2

c					Set up coefficients 
c					of linear system

      do 5 i=1,n
         a(1,i)=sngl(p(i,1)-x(i))
 5    continue
      b(1)=0.
      do 6 i=1,n
         dx1 = p(i,1)-xo(i)
         dx2 = x(i)-xo(i)
         b(1)=b(1) + sngl(dx1*dx1 - dx2*dx2)
 6    continue
      b(1)=b(1)/2.

      do 7 j=2,m
         do 8 i=1,n
            a(j,i)=sngl(p(i,j)-p(i,1))
 8       continue 
         b(j)=0.
         do 9 i=1,n
            dx1 = p(i,j)-xo(i)
            dx2 = p(i,1)-xo(i)
            b(j)=b(j) + sngl(dx1*dx1 - dx2*dx2)
 9       continue
         b(j)=b(j)/2.
 7    continue


      call volumef (a,b,m,n,mmax,nmax,vol)

      if(vol.eq.-1)then
         write(*,*)' Warning in routine second_voronoi'
         write(*,*)' Second order Voronoi cell is unbounded'
         write(*,*)' volume is set to zero'
         vol = 0.
      else if(vol.eq.0.)then
         write(*,*)' Warning in routine second_voronoi'
         write(*,*)' An inconsistent set of constraints has been'
         write(*,*)' detected in calculating the volume of '
         write(*,*)' the second-order Voronoi cell.'
         write(*,*)' volume is set to zero'
      end if

      return
      end
c
c------------------------------------------------------------------------
c
c	second_voronoi_d - calculates the area or volume of the second-order
c		 	   voronoi polygon (polyhedron) between x and its
c			   neighbour p(n,1). It also calculates the
c			   derivatives of the area (volume) with respect to
c			   the co-ordinates of x. 
c
c	Input:
c		x(n)			input point	
c		p(n,m)			neighbours shared by x and p(n,1)
c		n			number of dimensions
c		m			number of neighbours
c		nmax			maximum number of dimensions
c		mmax			maximum number of neighbours
c		xo(n)			a local origin point.
c		a(m,n)			work array 
c		b(m)			work array 
c
c	Output:
c		vol			volume (n=3) or area (n=2) of 
c					voronoi cell about x.
c		deriv(n)		first derivatives dv/dx_i (i=1,n)
c
c	Comments:
c		 See comments for routine second_voronoi.
c
c		 Calls C routines volumef, volumeb and dvdaf.
c
c					J. Braun, 
c					& M. Sambridge, RSES, 1995.
c
c------------------------------------------------------------------------
c
      subroutine second_voronoi_d (x,p,n,m,mmax,nmax,xo,a,b,vol,deriv)

      real*8	x(nmax)
      real*8	p(nmax,mmax)
      real*8    xo(nmax)
      real  	a(mmax,nmax)
      real	b(mmax)
      real	deriv(nmax)
      real*8    dx1,dx2

c					Set up coefficients 
c					of linear system

      do 5 i=1,n
         a(1,i)=sngl(p(i,1)-x(i))
 5    continue
      b(1)=0.
      do 6 i=1,n
         dx1 = p(i,1)-xo(i)
         dx2 = x(i)-xo(i)
         b(1)=b(1) + sngl(dx1*dx1 - dx2*dx2)
 6    continue
      b(1)=b(1)/2.

      do 7 j=2,m
         do 8 i=1,n
            a(j,i)=sngl(p(i,j)-p(i,1))
 8       continue 
         b(j)=0.
         do 9 i=1,n
            dx1 = p(i,j)-xo(i)
            dx2 = p(i,1)-xo(i)
            b(j)=b(j) + sngl(dx1*dx1 - dx2*dx2)
 9       continue
         b(j)=b(j)/2.
 7    continue

c     write(6,*)' A and b'
c     do 100 j=1,m
c        write(6,*)(a(j,i),i=1,n),' : ',b(j)
c100  continue
c     write(6,*)' '
c						calculate volume and
c						derivative d(vol)/db(1)

      call volumeb (a,b,m,n,mmax,nmax,1,vol,dvdb)

      if(vol.eq.-1)then
         write(*,*)' Warning in routine second_voronoi_d'
         write(*,*)' Second order Voronoi cell is unbounded'
         write(*,*)' volume and derivatives set to zero'
         vol = 0.
         do 10 k=1,n
            deriv(k) = 0.
 10      continue
      else if(vol.eq.0.)then
         write(*,*)' Warning in routine second_voronoi_d'
         write(*,*)' An inconsistent set of constraints has been'
         write(*,*)' detected in calculating the volume of '
         write(*,*)' the second-order Voronoi cell.'
         write(*,*)' volume and derivatives set to zero'
         do 11 k=1,n
            deriv(k) = 0.
 11      continue
      else
c						calculate derivatives
c						d(vol)/da(1,j), (j=1,n),
c						and use formulae to find
c						d(vol)/dx(k).
         do 20 k=1,n
   
           call dvdaf (a,b,m,n,mmax,nmax,k,deriva,icode)

           if(icode.eq.-1)then
              write(*,*)' Warning in routine second_voronoi_d'
              write(*,*)' Second order Voronoi cell is unbounded'
              write(*,*)' derivatives set to zero'
           else if(icode.eq.1)then
              write(*,*)' Warning in routine second_voronoi_d'
              write(*,*)' An inconsistent set of constraints has been'
              write(*,*)' detected in calculating the volume of '
              write(*,*)' the second-order Voronoi cell.'
              write(*,*)' All derivatives set to zero'
              do 30 i=1,n
                 deriv(i) = 0.
 30           continue
           end if

           deriv(k) = -1.*(dvdb*(x(k)-xo(k)) + deriva)
 20      continue 
      end if

      return
      end
c
c------------------------------------------------------------------------
c
c	Heapsort - Sorts an array into ascending order.
c		   Modified from numerical recipes 2, to use a 
c		   two dimensional double precision array.
c
c		   Sorting is done on index ka (ka=1 or 2).
c
c------------------------------------------------------------------------
c
c
      SUBROUTINE hpsort_d(n,ka,ra)
      INTEGER n
c     REAL ra(n)
      REAL*8 ra(2,n)
      INTEGER i,ir,j,l,ka
c     REAL rra
      REAL*8 rra(2)
      kb = 1
      if(ka.eq.1)kb = 2
      if (n.lt.2) return
      l=n/2+1
      ir=n
10    continue
        if(l.gt.1)then
          l=l-1
          rra(ka)=ra(ka,l)
          rra(kb)=ra(kb,l)
        else
          rra(ka)=ra(ka,ir)
          rra(kb)=ra(kb,ir)
          ra(ka,ir)=ra(ka,1)
          ra(kb,ir)=ra(kb,1)
          ir=ir-1
          if(ir.eq.1)then
            ra(ka,1)=rra(ka)
            ra(kb,1)=rra(kb)
            return
          endif
        endif
        i=l
        j=l+l
20      if(j.le.ir)then
          if(j.lt.ir)then
            if(ra(ka,j).lt.ra(ka,j+1))j=j+1
          endif
          if(rra(ka).lt.ra(ka,j))then
            ra(ka,i)=ra(ka,j)
            ra(kb,i)=ra(kb,j)
            i=j
            j=j+j
          else
            j=ir+1
          endif
        goto 20
        endif
        ra(ka,i)=rra(ka)
        ra(kb,i)=rra(kb)
      goto 10
      END
c
c------------------------------------------------------------------------
c
c	pangle - pseudo angle routine 
c
c		 returns a number between 0 and 360 which is NOT the
c		 angle made by the line from p1 to p2 with the horizontal
c		 but which has the same order properties as that angle,
c		 i.e. has the same order of angles as arctan dy/dx.
c	 	 This function involves only simple products and quotients.
c
c		 From Sedgewick (1990) `Algorithms in C' (Addison Wesley)
c
c						M. Sambridge 1996.
c
c------------------------------------------------------------------------
c
	Function pangle(p1,p2)

	real*8		p1(2)
	real*8		p2(2)

	dx = p2(1) - p1(1)
        ax = abs(dx)
	dy = p2(2) - p1(2)
        ay = abs(dy)
        t = 0.
        a = ax+ay
	if(a.ne.0.)then
           t = dy/a
        end if
        if(dx.lt.0.)then
          t = 2-t
        else if(dy.lt.0.)then 
          t = 4+t
        end if
        pangle = t*90

	return
	end
c
c------------------------------------------------------------------------
c
c	Circum - calculates circum-centre of three points
c
c	Input:
c		pa,pb,pc		array of input points	
c
c	Output:
c               centre(3)               centre(j) (j=1,2) contains the
c                                       co-ordinates of the centre of 
c                                       circle passing through the
c					three input points. 
c	Comments:
c
c		 Solves 3x3 linear system of equations.
c
c		 No calls to other routines.
c
c					M. Sambridge, RSES, April 1996.
c
c------------------------------------------------------------------------
c
	Subroutine circum(pa,pb,pc,centre)
c
	real*8		pa(2),pb(2),pc(2)
	real*8		centre(2)
	real*8		x1,x2,x3,y1,y2,y3
	real*8		dx2m1,dx2p1,dy2m1,dy2p1
	real*8		dx3m1,dx3p1,dy3m1,dy3p1
	real*8		denom
c						Find centre of circum-circle
	   x1 = pa(1)
	   x2 = pb(1)
	   x3 = pc(1)
	   y1 = pa(2)
	   y2 = pb(2)
	   y3 = pc(2)

           dx2m1 = x2-x1
           dx2p1 = x2+x1
           dy2m1 = y2-y1
           dy2p1 = y2+y1
           dx3m1 = x3-x1
           dx3p1 = x3+x1
           dy3m1 = y3-y1
           dy3p1 = y3+y1
           denom = dx2m1*dy3m1-dx3m1*dy2m1

	   centre(1) = ((dx2m1*dx2p1 + dy2m1*dy2p1)*dy3m1
     &         -(dx3m1*dx3p1 + dy3m1*dy3p1)*dy2m1)/
     &         (denom)

	   centre(2) = (dx2m1*(dx3m1*dx3p1 + dy3m1*dy3p1) 
     &         -dx3m1*(dx2m1*dx2p1 + dy2m1*dy2p1))/ 
     &         (denom)
c
           centre(1) = centre(1)*0.5d0
           centre(2) = centre(2)*0.5d0

	return
	end
c
c------------------------------------------------------------------------
c
c	Circum_d - calculates circum-centre of three points and 
c	           derivatives of circum-centre with respect to 
c		   co-ordinates of first point.
c
c	Input:
c		pa,pb,pc		array of input points	
c
c	Output:
c               centre(3)               centre(j) (j=1,2) contains the
c                                       co-ordinates of the centre of 
c                                       circle passing through the
c					three input points. 
c		vx(2)			derivative of centre with respect
c					to x-component of pa
c		vy(2)			derivative of centre with respect
c					to y-component of pa
c	Comments:
c
c		 Solves 3x3 linear system of equations and calculates
c		 derivatives of circum-centre with respect to co-ordinates
c		 of input vector pa(2).
c
c		 No calls to other routines.
c
c					M. Sambridge, RSES, May 1996.
c
c------------------------------------------------------------------------
c
	Subroutine circum_d(pa,pb,pc,centre,vx,vy)
c
	real*8		pa(2),pb(2),pc(2)
	real*8		centre(2)
	real*8		x1,x2,x3,y1,y2,y3
	real*8		dx2m1,dx2p1,dy2m1,dy2p1
	real*8		dx3m1,dx3p1,dy3m1,dy3p1
	real*8		denom
	real*8		vx(2),vy(2)
c						Find centre of circum-circle
	   x1 = pa(1)
	   x2 = pb(1)
	   x3 = pc(1)
	   y1 = pa(2)
	   y2 = pb(2)
	   y3 = pc(2)

           dx2m1 = x2-x1
           dx2p1 = x2+x1
           dy2m1 = y2-y1
           dy2p1 = y2+y1
           dx3m1 = x3-x1
           dx3p1 = x3+x1
           dy3m1 = y3-y1
           dy3p1 = y3+y1
           denom = dx2m1*dy3m1-dx3m1*dy2m1

	   centre(1) = ((dx2m1*dx2p1 + dy2m1*dy2p1)*dy3m1
     &         -(dx3m1*dx3p1 + dy3m1*dy3p1)*dy2m1)/
     &         (denom)

	   centre(2) = (dx2m1*(dx3m1*dx3p1 + dy3m1*dy3p1) 
     &         -dx3m1*(dx2m1*dx2p1 + dy2m1*dy2p1))/ 
     &         (denom)
c
           centre(1) = centre(1)*0.5d0
           centre(2) = centre(2)*0.5d0
c						X-derivative
c
           dcxm1 = centre(1) - x1
           dcym1 = centre(2) - y1

           denum1 = (dy3m1 - dy2m1)/denom
           denum2 = (dx2m1 - dx3m1)/denom

	   vx(1) = dcxm1*denum1
	   vx(2) = dcxm1*denum2

c						Y-derivative
c
	   vy(1) = dcym1*denum1
	   vy(2) = dcym1*denum2

	return
	end
c
c------------------------------------------------------------------------
c
c	Circum_dd - calculates circum-centre of three points and 
c		    1st and 2nd derivatives of circum-centre with 
c		    respect to co-ordinates of first point.
c
c	Input:
c		pa,pb,pc		array of input points	
c
c	Output:
c               centre(3)               centre(j) (j=1,2) contains the
c                                       co-ordinates of the centre of 
c                                       circle passing through the
c					three input points. 
c		vx(2)			derivative of centre with respect
c					to x-component of pa
c		vy(2)			derivative of centre with respect
c					to y-component of pa
c	Comments:
c
c		 Solves 3x3 linear system of equations and calculates
c		 derivatives of circum-centre with respect to co-ordinates
c		 of input vector pa(2).
c
c		 No calls to other routines.
c
c					M. Sambridge, RSES, May 1996.
c
c------------------------------------------------------------------------
c
	Subroutine circum_dd(pa,pb,pc,centre,vx,vy,vxx,vyy,vxy)
c
	real*8		pa(2),pb(2),pc(2)
	real*8		centre(2)
	real*8		x1,x2,x3,y1,y2,y3
	real*8		dx2m1,dx2p1,dy2m1,dy2p1
	real*8		dx3m1,dx3p1,dy3m1,dy3p1
	real*8		denom
	real*8		vx(2),vy(2),vxx(2),vyy(2),vxy(2)
c
c						Find centre of circum-circle
	   x1 = pa(1)
	   x2 = pb(1)
	   x3 = pc(1)
	   y1 = pa(2)
	   y2 = pb(2)
	   y3 = pc(2)

           dx2m1 = x2-x1
           dx2p1 = x2+x1
           dy2m1 = y2-y1
           dy2p1 = y2+y1
           dx3m1 = x3-x1
           dx3p1 = x3+x1
           dy3m1 = y3-y1
           dy3p1 = y3+y1
           denom = dx2m1*dy3m1-dx3m1*dy2m1

	   centre(1) = ((dx2m1*dx2p1 + dy2m1*dy2p1)*dy3m1
     &         -(dx3m1*dx3p1 + dy3m1*dy3p1)*dy2m1)/
     &         (denom)

	   centre(2) = (dx2m1*(dx3m1*dx3p1 + dy3m1*dy3p1) 
     &         -dx3m1*(dx2m1*dx2p1 + dy2m1*dy2p1))/ 
     &         (denom)
c
           centre(1) = centre(1)*0.5d0
           centre(2) = centre(2)*0.5d0
c						X-derivative
c
           dcxm1 = centre(1) - x1
           dcym1 = centre(2) - y1

           denum1 = (dy3m1 - dy2m1)/denom
           denum2 = (dx2m1 - dx3m1)/denom

	   vx(1) = dcxm1*denum1
	   vx(2) = dcxm1*denum2
c						Y-derivative
c
	   vy(1) = dcym1*denum1
	   vy(2) = dcym1*denum2
c						Second derivatives
	   f11 = 2*vx(1) - 1.
	   f22 = 2*vy(2) - 1.
	   f12 = vy(1) + vx(2)

           vxx(1) = f11*denum1
           vxx(2) = f11*denum2
           vyy(1) = f22*denum1
           vyy(2) = f22*denum2
           vxy(1) = f12*denum1
           vxy(2) = f12*denum2

	return
	end
c
c------------------------------------------------------------------------
c
c						heapsort modified to 
c						sort two arrays  
c
c------------------------------------------------------------------------
c
      SUBROUTINE hpsort_two(n,ra,rb)
      INTEGER n
      REAL ra(n)
      REAL*8 rb(2,n)
      INTEGER i,ir,j,l
      REAL rra
      REAL*8 rrb1,rrb2
      if (n.lt.2) return
      l=n/2+1
      ir=n
10    continue
        if(l.gt.1)then
          l=l-1
          rra=ra(l)
          rrb1=rb(1,l)
          rrb2=rb(2,l)
        else
          rra=ra(ir)
          rrb1=rb(1,ir)
          rrb2=rb(2,ir)
          ra(ir)=ra(1)
          rb(1,ir)=rb(1,1)
          rb(2,ir)=rb(2,1)
          ir=ir-1
          if(ir.eq.1)then
            ra(1)=rra
            rb(1,1)=rrb1
            rb(2,1)=rrb2
            return
          endif
        endif
        i=l
        j=l+l
20      if(j.le.ir)then
          if(j.lt.ir)then
            if(ra(j).lt.ra(j+1))j=j+1
          endif
          if(rra.lt.ra(j))then
            ra(i)=ra(j)
            rb(1,i)=rb(1,j)
            rb(2,i)=rb(2,j)
            i=j
            j=j+j
          else
            j=ir+1
          endif
        goto 20
        endif
        ra(i)=rra
        rb(1,i)=rrb1
        rb(2,i)=rrb2
      goto 10
      END
c
c------------------------------------------------------------------------
c
c					Numerical Recipes routine index
c
c------------------------------------------------------------------------
c
      SUBROUTINE indexx(n,arr,indx)
      INTEGER n,indx(n),M,NSTACK
      REAL arr(n)
      PARAMETER (M=7,NSTACK=50)
      INTEGER i,indxt,ir,itemp,j,jstack,k,l,istack(NSTACK)
      REAL a
      do 11 j=1,n
        indx(j)=j
11    continue
      jstack=0
      l=1
      ir=n
1     if(ir-l.lt.M)then
        do 13 j=l+1,ir
          indxt=indx(j)
          a=arr(indxt)
          do 12 i=j-1,1,-1
            if(arr(indx(i)).le.a)goto 2
            indx(i+1)=indx(i)
12        continue
          i=0
2         indx(i+1)=indxt
13      continue
        if(jstack.eq.0)return
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        k=(l+ir)/2
        itemp=indx(k)
        indx(k)=indx(l+1)
        indx(l+1)=itemp
        if(arr(indx(l+1)).gt.arr(indx(ir)))then
          itemp=indx(l+1)
          indx(l+1)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l)).gt.arr(indx(ir)))then
          itemp=indx(l)
          indx(l)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l+1)).gt.arr(indx(l)))then
          itemp=indx(l+1)
          indx(l+1)=indx(l)
          indx(l)=itemp
        endif
        i=l+1
        j=ir
        indxt=indx(l)
        a=arr(indxt)
3       continue
          i=i+1
        if(arr(indx(i)).lt.a)goto 3
4       continue
          j=j-1
        if(arr(indx(j)).gt.a)goto 4
        if(j.lt.i)goto 5
        itemp=indx(i)
        indx(i)=indx(j)
        indx(j)=itemp
        goto 3
5       indx(l)=indx(j)
        indx(j)=indxt
        jstack=jstack+2
        if(jstack.gt.NSTACK)pause 'NSTACK too small in indexx'
        if(ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        endif
      endif
      goto 1
      END
C  (C) Copr. 1986-92 Numerical Recipes Software '%1&9p#!.
c
c------------------------------------------------------------------------
c
c	second_v_area - calculates the area of a second-order Voronoi 
c			cell using an un-ordered list of vertices and a
c		        a closed formula. 
c
c	Input:
c		x(2)			input point	
c		p(2,n)			vertices of polygon.
c		n			number of vertices
c		a(n)			pseudo angles of nodes with 
c					respect to input point
c
c	Output:
c		area			area of polygon X 2
c
c	Comments:
c		 The vertices may be input in any order but they are 
c		 re-ordered anti-clockwise upon output.
c		
c		 Calls are made to pangle, hpsort_two.
c
c					M. Sambridge, RSES, May 1996.
c
c------------------------------------------------------------------------
c
      Subroutine second_v_area(x,n,p,a,area)

      real*8	x(2)
      real*8	p(2,*)
      real*4    theta
      real*4    a(*)
c	      					calculate pseudo angles
c						from first node
c
      do i=2,n
        a(i) = pangle(p(1,1),p(1,i))
      end do
      a(1) = -1
c
      theta = a(2)
      do i=2,n
         a(i) = a(i) - theta
        if(a(i).lt.0)a(i) = a(i) + 360.
      end do
      theta = a(2)
c
c     write(*,*)' unsorted angles'
c     do i=1,n
c        write(*,*)i,' :',a(i),p(1,i),p(2,i)
c     end do
c						sort these nodes by 
c						pseudo angle from first edge
      call hpsort_two(n,a,p)

      if(a(2).ne.theta)write(*,*)' ERROR: first angle moved ?'

c     write(*,*)' sorted angles'
c     do i=1,n
c        write(*,*)i,' :',a(i),p(1,i),p(2,i)
c     end do
c						calculate area of 
c						second-order Voronoi cell.
      area = 0.
      do i=1,n-1
         j = i+1
         area = area + p(1,i)*p(2,j) - p(2,i)*p(1,j)
      end do
      area = area + p(1,n)*p(2,1) - p(2,n)*p(1,1)
      area = abs(area)
      
      return
      end
c
c------------------------------------------------------------------------
c
c	second_v_area_d - calculates the area of a second-order Voronoi 
c			  cell using an un-ordered list of vertices and a
c		          a closed formula. 
c
c	Input:
c		x(2)			input point	
c		p(2,n)			vertices of polygon.
c		dp(4,2)			derivatives of first two vertices 
c					of polygon.
c		n			number of vertices
c		a(n)			pseudo angles of nodes with 
c					respect to input point
c
c	Output:
c		area			area of polygon X 2
c		df(2)			first derivatives of area X 2
c					      df(1) = df/dx,
c					      df(2) = df/dy
c
c	Comments:
c		 The first two vertices are assumed to be the vertices
c		 dependent on x, i.e. the vertices on the voronoi cell about x.
c		 These two vertices must be in clockwise order.
c		 The remaining vertices may be input in any order but they are 
c		 re-ordered anti-clockwise upon output.
c		
c		 Calls are made to pangle, hpsort_two.
c
c					M. Sambridge, RSES, April 1996.
c
c------------------------------------------------------------------------
c
      Subroutine second_v_area_d(x,n,p,dp,a,area,df)

      real*8	x(2)
      real*8	p(2,*)
      real*8	dp(4,2)
      real*8	df(2)
      real*4    theta
      real*4    a(*)
c
c	      					calculate pseudo angles
c						from first node
c
      do i=2,n
        a(i) = pangle(p(1,1),p(1,i))
      end do
      a(1) = -1
c
c     write(*,*)' first pseudo angle =',theta 
c
      theta = a(2)
      do i=2,n
         a(i) = a(i) - theta
        if(a(i).lt.0)a(i) = a(i) + 360.
      end do
      theta = a(2)
c
c     write(*,*)' unsorted angles'
c     do i=1,n
c        write(*,*)i,' :',a(i),p(1,i),p(2,i)
c     end do
c						sort these nodes by 
c						pseudo angle from first edge
      call hpsort_two(n,a,p)

      if(a(2).ne.theta)write(*,*)' ERROR: first angle moved ?'

c						plot nodes and v-cell
c     write(*,*)' sorted angles'
c     do i=1,n
c        write(*,*)i,' :',a(i),p(1,i),p(2,i)
c     end do
c						calculate area of 
c						second-order Voronoi cell.
      area = 0.
      do i=1,n-1
         j = i+1
         area = area + p(1,i)*p(2,j) - p(2,i)*p(1,j)
      end do
      area = area + p(1,n)*p(2,1) - p(2,n)*p(1,1)
c
c						calculate 1st derivatives
      d222n = p(2,2) - p(2,n) 
      d1n12 = p(1,n) - p(1,2) 
      d2321 = p(2,3) - p(2,1)
      d1113 = p(1,1) - p(1,3)
      df(1) =   dp(1,1)*d222n + dp(2,1)*d1n12 
     &        + dp(1,2)*d2321 + dp(2,2)*d1113
      df(2) =   dp(3,1)*d222n + dp(4,1)*d1n12 
     &        + dp(3,2)*d2321 + dp(4,2)*d1113

      if(area.lt.0.)then
        df(1) = -df(1)
        df(2) = -df(2)
      end if

      area = abs(area)
       
      return
      end
c
c------------------------------------------------------------------------
c
c	second_v_area_dd - calculates the area of a second-order Voronoi 
c			   cell using an un-ordered list of vertices and a
c		           a closed formula. 
c
c	Input:
c		x(2)			input point	
c		p(2,n)			vertices of polygon.
c		dp(4,2)			derivatives of first 
c					two vertices of polygon.
c		n			number of vertices
c		a(n)			pseudo angles of nodes with 
c					respect to input point
c
c	Output:
c		area			area of polygon X 2
c		df(2)			first derivatives of area X 2
c					      df(1) = df/dx
c					      df(2) = df/dy
c		ddf(3)			second derivatives of area X 2
c					      ddf(1) = d2f/dxx
c					      ddf(2) = d2f/dyy
c					      ddf(3) = d2f/dxy
c
c	Comments:
c		 The first two vertices are assumed to be the vertices
c		 dependent on x, i.e. the vertices on the voronoi cell about x.
c		 These two vertices must be in clockwise order.
c		 The remaining vertices may be input in any order but they are 
c		 re-ordered anti-clockwise upon output.
c		
c		 Calls are made to pangle, hpsort_two and xplot routines.
c
c					M. Sambridge, RSES, April 1996.
c
c------------------------------------------------------------------------
c
      Subroutine second_v_area_dd(x,n,p,dp,ddp,a,area,df,ddf)

      real*8	x(2)
      real*8	p(2,*)
      real*8	dp(4,2)
      real*8	ddp(6,2)
      real*8	df(2),ddf(3)
      real*4    theta
      real*4    a(*)
c	      					calculate pseudo angles
c						from first node
c
      do i=2,n
        a(i) = pangle(p(1,1),p(1,i))
      end do
      a(1) = -1
c
c     write(*,*)' first pseudo angle =',theta 
c
      theta = a(2)
      do i=2,n
        a(i) = a(i) - theta
        if(a(i).lt.0)a(i) = a(i) + 360.
      end do
      theta = a(2)
c
c     write(*,*)' unsorted angles'
c     do i=1,n
c        write(*,*)i,' :',a(i),p(1,i),p(2,i)
c     end do
c						sort these nodes by 
c						pseudo angle from first edge
      call hpsort_two(n,a,p)

      if(a(2).ne.theta)write(*,*)' ERROR: first angle moved ?'

c						plot nodes and v-cell
c     write(*,*)' sorted angles'
c     do i=1,n
c        write(*,*)i,' :',a(i),p(1,i),p(2,i)
c     end do
c						calculate area of 
c						second-order Voronoi cell.
      area = 0.
      do i=1,n-1
         j = i+1
         area = area + p(1,i)*p(2,j) - p(2,i)*p(1,j)
      end do
      area = area + p(1,n)*p(2,1) - p(2,n)*p(1,1)
c
c						calculate 1st derivatives
      d222n = p(2,2) - p(2,n) 
      d1n12 = p(1,n) - p(1,2) 
      d2321 = p(2,3) - p(2,1)
      d1113 = p(1,1) - p(1,3)
      df(1) =   dp(1,1)*d222n + dp(2,1)*d1n12 
     &        + dp(1,2)*d2321 + dp(2,2)*d1113
      df(2) =   dp(3,1)*d222n + dp(4,1)*d1n12 
     &        + dp(3,2)*d2321 + dp(4,2)*d1113

c						calculate 2nd derivatives
     
      ddf(1) =   ddp(1,1)*d222n + ddp(2,1)*d1n12 
     &         + ddp(1,2)*d2321 + ddp(2,2)*d1113
     &         + 2.*(dp(1,1)*dp(2,2) - dp(2,1)*dp(1,2))
      ddf(2) =   ddp(3,1)*d222n + ddp(4,1)*d1n12 
     &         + ddp(3,2)*d2321 + ddp(4,2)*d1113
     &         + 2.*(dp(3,1)*dp(4,2) - dp(4,1)*dp(3,2))
      ddf(3) =   ddp(5,1)*d222n + ddp(6,1)*d1n12 
     &         + ddp(5,2)*d2321 + ddp(6,2)*d1113
     &         + dp(1,1)*dp(4,2) - dp(2,1)*dp(3,2)
     &         + dp(3,1)*dp(2,2) - dp(4,1)*dp(1,2)

      if(area.lt.0.)then
        df(1) = -df(1)
        df(2) = -df(2)
        ddf(1) = -ddf(1)
        ddf(2) = -ddf(2)
        ddf(3) = -ddf(3)
      end if

      area = abs(area)
       
      return
      end
