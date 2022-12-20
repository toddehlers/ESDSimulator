c nn_remove

c
c	nn_remove - Removes a node from the triangulation and updates
c		    vertex and adjacency array using only local search.
c		    Optionally re-numbers the nodes from 1 to np-1.
c
c	Input:	
c               node                    node to be removed 
c               np                      number of nodes
c	        nt			number of triangles
c               points(2,np)            array of node co-ordinates
c	        vertices(3,nt)		array of triangle vertices	
c		neighbour(3,nt)		adjacency array
c		nnpn_max		maximum number of neighbours per node
c					(used to set array sizes)
c		vis_tlist		work array used by `delaun'
c		vis_elist		work array used by `delaun'
c		add_elist		work array used by `delaun'
c		nv_max			size of work arrays (<=nnpn_max)
c		c_list			work array used by `find_ctriangles'
c		nodelist		work array used by `find_ctriangles'
c		tlist			work array used by `update_vn_lists'
c		loc			a first guess at a triangle attached 
c					to the input node. Can be any value
c					between 1<= loc <= nt.
c		renumber_nodes		logical. If true then points are repacked
c
c	Output:
c		nt		        number of triangles 
c	        vertices(3,nt)		array of triangle vertices	
c		neighbour(3,nt)		adjacency array
c
c
c	Comments:
c		 If the logical renumber_nodes is set to true then
c		 the node numbers will also be renumbered from 1 to np-1.
c		 This is done by 1) renaming node np as the input node if
c		 node is not already equal to np, and 2) set np = np-1.
c		 In this mode it is assumed that node np exists in the
c		 current triangles. An error will occur if the node np
c		 is removed without renumbering node numbers and subsequently
c		 a call is made with node renumbering. i.e.  nn_remove is 
c		 called with node=np and renumber_nodes = .false., and then 
c		 a later call is made with renumber_nodes = .true.
c		 Be careful in mixing calls with renumber_nodes changing.
c
c		 Assumes input list of vertices in anticlockwise order
c
c		 The arrays must be dimensioned in the following way:
c
c		 real*8  points(np_max)
c		 integer vertices(3,nt_max)
c		 integer neighbour(3,nt_max)
c
c		 integer nodelist(nnpn_max)    : work array
c		 integer tlist(nnpn_max)       : work array
c		 logical c_list(2*nnpn_max)    : work array
c		 integer v_local(3,2*nnpn_max) : work array
c		 integer n_local(3,2*nnpn_max) : work array
c      		 integer vis_tlist(nv_max)     : work array
c		 integer vis_elist(nv_max)     : work array
c	         integer add_tlist(nv_max)     : work array
c
c		 where np_max = maximum number of points, nt_max = maximum
c		 number of triangles (<2*np_max) and nnpn_max is the maximum 
c		 number of neighbours about any node.
c
c		 The work arrays vis_tlist, vis_elist,and add_tlist must be 
c		 dimensioned size nv_max, (the maximum number of
c		 triangles visible from any point from outside of the hull). 
c		 In this case delaun is only being used for a maximum of 
c		 nnpn_max points and so nv_max will always be <= nnpn_max.
c
c		 Calls to find_ctriangles and update_vn_lists.
c
c					M. Sambridge, RSES, Nov. 1995.
c
c------------------------------------------------------------------------
c
	Subroutine nn_remove 
     &             (node,np,nt,points,vertices,neighbour,loc,
     &              nnpn_max,nv_max,vis_tlist,vis_elist,add_tlist,
     &              v_local,n_local,c_list,nodelist,tlist,
     &              renumber_nodes)

      common /vocal/ ivocal,iprofile

	integer		vertices(3,*)
	integer		neighbour(3,*)
	integer		nodelist(*)
	integer		tlist(*)
        integer         vis_tlist(*)
        integer         vis_elist(*)
        integer         add_tlist(*)
	integer		v_local(3,*)
	integer		n_local(3,*)
        logical         c_list(*)
        logical         on_hull
        logical         renumber_nodes
        real*8		points(2,*)
        logical		debug_mode
        common		/debug_nn_remove/debug_mode

c						Check error conditions
	if(np.eq.3)then
           write(*,100)
           stop
        end if
c						find list of triangles 
c						and nodes attached to input node
c
	call find_neighbours_r
     &       (node,points,nnpn_max,vertices,neighbour,
     &        loc,on_hull,n_neighbours,n_triangles,nodelist,tlist)

c						find the circum-triangles
c						of the current node 

c
        np_local = n_neighbours
	call find_ctriangles 
     &       (node,np,nt,points,vertices,neighbour,nodelist,
     &        np_local,on_hull,nnpn_max,nv_max,vis_tlist,
     &        vis_elist,add_tlist,nt_local,v_local,n_local,c_list)

c						replace current node and
c						its attached triangles with
c						circum-triangles in
c						vertices and neighbour lists
c
 	call update_vn_lists
     &      (node,np,nt,vertices,neighbour,tlist,
     &       on_hull,nnpn_max,np_local,
     &       nt_local,v_local,n_local,c_list)
c
	if(.not.renumber_nodes)return
c
c						renumber nodes from 1 to np-1
c
	if(node.eq.np)then
          np = np -1
          return
        end if
c						renumber node np to the
c						value of the input node in
c						all triangles attached to it 
        node_old = np
        loc2 = nt
	call renumber_node
     &       (node_old,node,points,vertices,neighbour,
     &        loc2,on_hull,n_neighbours,n_triangles)
        np = np - 1


 100    format(/'Error in subroutine nn_remove:'//,
     &          'Only three nodes in triangulation'/, 
     &          'Cannot remove a node'/) 

	return
	end
c
c------------------------------------------------------------------------
c
c	renumber_node  - renumbers node `node' by `node_new' in
c			 arrays vertices and neighbour and points.
c
c	Input:	
c               node                    index of node to be removed
c               node_new                index of node to be added 
c               points(2,np)            array of node co-ordinates
c	        vertices(3,nt)		array of triangle vertices	
c		neighbour(3,nt)		adjacency array
c		loc			a first guess at a triangle
c					attached to input node `node'
c
c	Output:
c	        vertices(3,nt)		updated 
c		n_neighbours		number of neighbours of `node' 
c		n_triangles		number of triangles about `node'
c		on_hull			logical = .true. if node is 
c					on the convex hull
c
c	Comments:
c		 This routine returns the vertex and adjacency arrays with
c		 all references to node `node' replaced by `node_new'.
c		 This is done using only local loops (i.e. no loop
c		 over all triangles). Extra information in the form
c		 of n_neighbours, n_triangles, on_hull is, returned and 
c		 may be useful. 
c
c		 Calls are made to node_loc.
c
c					M. Sambridge, RSES, Jan. 1996.
c
c------------------------------------------------------------------------
c
	Subroutine renumber_node
     &             (node,node_new,points,vertices,neighbour,
     &              loc,on_hull,n_neighbours,n_triangles)

	real*8		points(2,*)
	integer		vertices(3,*)
	integer		neighbour(3,*)
	integer		t,t0,tnew
	integer		p,p0,pnew
        integer		c1(3)
        integer		c2(3)
        integer		edg
        logical		on_hull
        logical		debug_mode
        common		/debug_nn_remove/debug_mode
        data		c1/2,3,1/
        data		c2/3,1,2/

c					find a single triangle 
c					attached to the input node

        t0 = loc

        if(debug_mode)then
           write(*,*)' '
           write(*,*)' start renumber_node node = ',node
           write(*,*)' input to node_loc t0 =',t0
        end if

	call node_loc(node,points,vertices,neighbour,t0,p0)

        if(debug_mode)
     &  write(*,*)' output of node loc t0 =',t0,' pos =',p0

c					find remaining nodes and
c					triangles attached to input node
        on_hull = .false.
        k = 2
        t = t0
        p = p0
        vertices(p,t) = node_new
  5     tnew = neighbour(c1(p),t)
        if(debug_mode)then
           write(*,*)'first k',k
           write(*,*)'first tnew',tnew
        end if
        if(tnew.eq.0)then
           if(debug_mode)write(*,*)' reached hull on first sweep'
           on_hull = .true.
           t = t0
           p = p0
  6        tnew = neighbour(c2(p),t)
           if(debug_mode)then
              write(*,*)'second k',k
              write(*,*)'second tnew',tnew
              if(tnew.eq.0)write(*,*)' reached hull on second sweep'
           end if
           if(tnew.eq.0)go to 10
           pnew = edg(tnew,node,vertices)
           k = k + 1
           t = tnew
           p = pnew
           vertices(p,t) = node_new
           go to 6
        else if(tnew.eq.t0)then
           k = k - 1
           go to 10
        end if
        pnew = edg(tnew,node,vertices)
        k = k + 1
        t = tnew
        p = pnew
        vertices(p,t) = node_new
        go to 5
 10     continue

        n_neighbours = k

c						renumber node in 
c						points array
        points(1,node_new) = points(1,node)
        points(2,node_new) = points(2,node)

        if(debug_mode)then
           write(*,*)' '
           write(*,*)' done renumber_node '
           write(*,*)' '
        end if

	return
	end
c
c------------------------------------------------------------------------
c
c	find_neighbours_r - finds the list of nodes that are neighbours
c			    of an input node and the triangles attached to
c			    that node.
c
c	Input:	
c               node                    input node 
c               points(2,np)            array of node co-ordinates
c	        vertices(3,nt)		array of triangle vertices	
c		neighbour(3,nt)		adjacency array
c		loc			a first guess at a triangle
c					attached to input node
c
c	Output:
c		n_neighbours		number of neighbours of node 
c		n_triangles		number of triangles attached to node
c	        nodelist		list of neighbours to input node 
c	        tlist		        list of triangles attached to input node
c		on_hull			logical = .true. if node is on the
c					convex hull
c
c	Comments:
c		 Arrays must be dimensioned in the following way
c
c		 integer nodelist(nnpn_max)    
c		 integer tlist(nnpn_max)      
c
c		 where nnpn_max is the maximum number of neighbours per node.
c
c		 An error is reported if the number of neighbours of
c		 the input node is greater than nnpn_max. This is 
c		 necessary because in this case work arrays may not be 
c		 large enough.
c
c		 Calls are made to node_loc.
c
c					M. Sambridge, RSES, Dec. 1995.
c
c------------------------------------------------------------------------
c
	Subroutine find_neighbours_r
     &             (node,points,nnpn_max,vertices,neighbour,
     &              loc,on_hull,n_neighbours,n_triangles,nodelist,tlist)

	real*8		points(2,*)
	integer		vertices(3,*)
	integer		neighbour(3,*)
	integer		tlist(*)
	integer		nodelist(*)
	integer		t,t0,tnew
	integer		p,p0,pnew
        integer		c1(3)
        integer		c2(3)
        integer		edg
        logical		on_hull
        data		c1/2,3,1/
        data		c2/3,1,2/
        logical		debug_mode
        common		/debug_nn_remove/debug_mode

c					find a single triangle 
c					attached to the input node

        t0 = loc
        if(debug_mode)then
           write(*,*)' '
           write(*,*)' start find_neighbours '
           write(*,*)' start find_neighbours node = ',node
           write(*,*)' input to node_loc t0 =',t0
        end if

	call node_loc(node,points,vertices,neighbour,t0,p0)

        if(debug_mode)
     &  write(*,*)' output of node loc t0 =',t0,' pos =',p0

c					find remaining nodes and
c					triangles attached to input node
        on_hull = .false.
	k = 1
        tlist(k) = t0
        nodelist(k) = vertices(c1(p0),t0)
        k = k + 1
        nodelist(k) = vertices(c2(p0),t0)
        t = t0
        p = p0
  5     tnew = neighbour(c1(p),t)
        if(debug_mode)then
           write(*,*)'first k',k,' tlist',(tlist(kk),kk=1,k-1)
           write(*,*)'first nodelist',(nodelist(kk),kk=1,k)
           write(*,*)'first tnew',tnew
        end if
        if(tnew.eq.0)then
           if(debug_mode)write(*,*)' reached hull on first sweep'
           on_hull = .true.
           t = t0
           p = p0
  6        tnew = neighbour(c2(p),t)
           if(debug_mode)then
              write(*,*)'second k',k,' tlist',(tlist(kk),kk=1,k-1)
              write(*,*)'second nodelist',(nodelist(kk),kk=1,k)
              write(*,*)'second tnew',tnew
              if(tnew.eq.0)write(*,*)' reached hull on second sweep'
           end if
           if(tnew.eq.0)go to 10
           tlist(k) = tnew
           pnew = edg(tnew,node,vertices)
           k = k + 1
           nodelist(k) = vertices(c1(pnew),tnew)
           t = tnew
           p = pnew
           go to 6
        else if(tnew.eq.t0)then
           nodelist(k) = 0
           k = k - 1
           go to 10
        end if
        tlist(k) = tnew
        pnew = edg(tnew,node,vertices)
        k = k + 1
        nodelist(k) = vertices(c2(pnew),tnew)
        t = tnew
        p = pnew
        go to 5
 10     continue

        n_neighbours = k
c						check size of nnpn_max

        call check_nnpn_max(n_neighbours,node,nnpn_max)

        if(on_hull)then
           n_triangles = n_neighbours-1
        else
           n_triangles = n_neighbours
        end if

        if(debug_mode)then
           write(*,*)' on_hull =',on_hull
           write(*,*)' n_neighbours =',n_neighbours
           write(*,*)' nodes:',(nodelist(i),i=1,n_neighbours)
           write(*,*)' triangles:',(tlist(i),i=1,n_triangles)
           write(*,*)' '
           write(*,*)' done find_neighbours '
           write(*,*)' '
        end if

	return
	end
c
c------------------------------------------------------------------------
c
c	Node_loc - locates a triangle attached to a particular node
c
c	Input:
c		node			input node
c		points(2,np)		array of node points	
c	        vertices(3,nt)		array of triangle vertices	
c	        neighbour(3,nt)		array of neighbouring triangles.	
c					Neighbour(i,j) is the triangle
c					opposite node i in triangle j,
c					stored counterclockwise about j.
c	        loc			first guess of triangle containing
c					input node
c
c	Output:
c	        loc			index of triangle containing 
c					input point.
c	        pos			position of node in triangle
c
c	Comments:
c
c		 No calls to other routines.
c
c					M. Sambridge, RSES, April 1994.
c
c------------------------------------------------------------------------
c
	Subroutine node_loc(node,points,vertices,neighbour,loc,pos)
c
	real*8		points(2,*)
	integer		vertices(3,*)
	integer		neighbour(3,*)
	integer		p1,p2
        integer		pos
	real*8		x,y,del1,del2
        integer		c1(3)
        data		c1/2,3,1/
c
c						check starting triangle
        do i=1,3
           k = vertices(i,loc)
           if(node.eq.k)then
             pos = i
             return
           end if
        end do

        x = points(1,node)
        y = points(2,node)

 10     continue

        do 20 i=1,3
	   j = c1(i)
	   k = c1(j)
           p1 = vertices(i,loc)
           p2 = vertices(j,loc)
	   del1 = (points(2,p1)-y)*(points(1,p2)-x)
	   del2 = (points(1,p1)-x)*(points(2,p2)-y)
	   if(del1.gt.del2)then
	      if(neighbour(k,loc).eq.0)then
                 write(*,*)' '
                 write(*,*)' Error in subroutine node_loc'
                 write(*,*)' '
                 write(*,*)
     &           ' Cannot find a triangle attached to node',node
                 write(*,*)' Node is possibly outside of convex hull'
                 write(*,*)' or already removed from the list ?'
                 stop
	      else
	         loc = neighbour(k,loc)
	      end if
	      go to 10
	   end if
 20     continue
	
c						find position in triangle
        do i=1,3
           k = vertices(i,loc)
           if(node.eq.k)then
             pos = i
             return
           end if
        end do

        write(*,*)' '
	write(*,*)' Error in subroutine node_loc'
        write(*,*)' '
        write(*,*)' Cannot find a triangle attached to node',node
        write(*,*)' Has node already removed from the list ?'
	write(*,*)' Last triangle encountered : ',loc
	write(*,*)' with vertex nodes         : ',
     &            (vertices(k,loc),k=1,3)

        stop
	end
c
c------------------------------------------------------------------------
c
c	find_ctriangles - finds the circum-triangles of input node among
c			  its set of neighbouring nodes (stored in nodelist).
c
c	Input:	
c               node                    input node 
c               np                      number of nodes
c               points(2,np)            array of node co-ordinates
c	        nt			number of triangles
c	        vertices(3,nt)		array of triangle vertices	
c		neighbour(3,nt)		adjacency array
c	        nodelist		list of neighbours to input node 
c		vis_tlist		work array used by delaun
c		vis_elist		work array used by delaun
c		add_elist		work array used by delaun
c		nv_max			size of work arrays
c		np_local		number of neighbours of node 
c		on_hull			logical: true if point is on convex hull
c
c	Output:
c		v_local			array of triangle vertices between
c					neighbours of node
c		n_local			local adjacency array 
c		c_list			list of triangles in v_local
c					which are circumcircles of node
c		nt_local		number of triangles attached to node
c
c	Comments:
c		 This routine calculates the Delaunay triangulation of
c		 the neighbours of the input node. The result is put into
c		 v_local and n_local. It also determines which of these
c		 local triangles are not circum-triangles of the input node
c		 and flags these in the logical array c_list.
c		 This information can be used to remove the input node
c		 from the vertex and neighbour lists of the original 
c		 triangulation.
c
c		 Assumes input list of vertices in anticlockwise order
c
c		 The arrays must be dimensioned in the following way:
c
c		 integer nodelist(nnpn_max)    
c		 integer tlist(nnpn_max)      
c		 logical c_list(2*nnpn_max)   
c		 integer v_local(3,2*nnpn_max) 
c		 integer n_local(3,2*nnpn_max) 
c      		 integer vis_tlist(nv_max)     : work array
c		 integer vis_elist(nv_max)     : work array
c	         integer add_tlist(nv_max)     : work array
c
c		 where nnpn_max is the maximum number of neighbours per node. 
c		 (Note a strict upper bound on the number of triangles 
c		 attached to a node is 2*nnpn_max-2.) 
c
c		 nv_max is the size of the work arrays vis_tlist,
c		 vis_elist, add_tlist, used by delaun. In this case
c		 delaun is only being used for a maximum of nnpn_max points
c		 and so nv_max <= nnpn_max.
c
c		 Calls are made to delaun.
c
c					M. Sambridge, RSES, Jan. 1996.
c
c------------------------------------------------------------------------
c
	Subroutine find_ctriangles 
     &        (node,np,nt,points,vertices,neighbour,nodelist,
     &         np_local,on_hull,nnpn_max,nv_max,vis_tlist,vis_elist,
     &         add_tlist,nt_local,v_local,n_local,c_list)

	integer		vertices(3,*)
	integer		neighbour(3,*)
	integer		nodelist(*)
	integer		v_local(3,*)
	integer		n_local(3,*)
        integer         vis_tlist(*)
        integer         vis_elist(*)
        integer         add_tlist(*)
        logical         c_list(*)
        logical         on_hull
        real*8		points(2,*)
        real*8		circle(3)
        real*8		x,y,dx,dy,dist
        real*8		eps
        logical		debug_mode
        logical		ldummy
        common		/debug_nn_remove/debug_mode
	common 	        /roundoff/eps
c						check error conditions 
        if(node.gt.np)then
           write(*,100)node,np
           stop
        end if
        if(node.lt.1)then
           write(*,200)node
           stop
        end if

        if(eps.eq.0.d0)eps = 1.E-9
        nt_local_max = 2*nnpn_max

c						calculate Delaunay 
c						triangulation of neighbours

        if(np_local.gt.2)then
	   call delaun 
     &              (points,np_local,n_local,v_local,nt_local,
     &               nt_local_max,vis_tlist,vis_elist,
     &               add_tlist,eps,nv_max,2,ldummy,0,0,nodelist)
c
c	   call delaun_subset 
c    &              (points,np_local,n_local,v_local,nt_local,
c    &               nt_local_max,nodelist,vis_tlist,vis_elist,
c    &               add_tlist,eps,nv_max)
        else
           nt_local = 0
        end if

c                                               write local Delaunay and nodes
        if(debug_mode)then
	   write(*,*)' node = ',node
	   write(*,*)' np_local',np_local,' nt_local',nt_local
	   write(*,*)' neighbours '
	   write(*,*)(nodelist(i),i=1,np_local)
           do j=1,nt_local
	      write(*,*)' v:',(v_local(k,j),k=1,3),
     &                  ' n:',(n_local(k,j),k=1,3)
           end do

c                                               draw local Delaunay and nodes
c          call xpen(2)
c          call draw_nodes_sub
c    &          (np_local,nodelist,points,0.0075)
c          call xpen(3)
c          call xthick(2.)
c          call draw_delaunay
c    &          (np_local,points,nt_local,v_local)
        end if

c						If the number of triangles
c						differs from np-2 then
c						we must have a convex
c						set of neighbours and we
c						need to remove some
c
	if(nt_local.eq.np_local-2.and..not.on_hull)then
           do i=1,nt_local
              c_list(i) = .true.
           end do
           if(debug_mode)then
              write(*,*)' No extra triangles to remove'
              write(*,*)' The set of neighbours is convex'
           end if
        else
c						find local triangles 
c						whose circumcircles do not
c						contain node node
           do j=1,nt_local
c						determine if node node is
c						inside circum-circle of
c						current triangle

              call circum_centre(j,v_local,points,circle)
              x = points(1,node)
              y = points(2,node)
              dx = circle(1)-x
              dy = circle(2)-y
              dist = dx*dx + dy*dy
      
              if(dist.gt.circle(3))then
                 if(debug_mode)write(*,*)j,' not circum-triangle'
                 c_list(j) = .false.
              else
                 c_list(j) = .true. 
                 if(debug_mode)write(*,*)j,' circum-triangle'
              end if
           end do
c                                               write local Delaunay and nodes
	   if(debug_mode)then
              write(*,*)' reduced triangles'
              do j=1,nt_local
                 if(c_list(j))then
	            write(*,*)' v:',(v_local(k,j),k=1,3),
     &                        ' n:',(n_local(k,j),k=1,3)
                 end if
              end do
	   end if

	end if

 100    format(/' Error in subroutine find_ctriangles:'//, 
     &          ' Node number is greater than the number of nodes'/, 
     &          ' Input node:',i8,' Number of nodes:',i8/) 
 200    format(/' Error in subroutine find_ctriangles:'//, 
     &          ' Node number is negative '/, 
     &          ' Input node:',i8/)
 300    format(/' Error in subroutine find_ctriangles:'/, 
     &          ' Neighbour list is empty'/,
     &          ' Has Neighbour list been built ?'/) 

	return
	end
c
c------------------------------------------------------------------------
c
c	update_vn_lists - Removes the triangles attached to node `node'
c			  from vertices and replaces them with those in 
c			  v_local. It also updates the adjacency array 
c			  neighbour using n_local.
c
c	Input:	
c               node                    input node 
c	        nt			number of triangles
c	        vertices(3,nt)		array of triangle vertices	
c		neighbour(3,nt)		adjacency array
c		v_local			array of triangle vertices between
c					neighbours of node
c		n_local			adjacency array for v_local 
c		c_list			list of triangles in v_local
c					which are circumcircles of node
c		nt_local		number of local triangles
c	        tlist		        list of triangles attached to input node
c		nnpn_max		maximum number of neighbours per node
c
c	Output:
c	        nt			updated number of triangles
c	        vertices(3,nt)		updated array of triangle vertices	
c		neighbour(3,nt)		updated adjacency array
c
c	Comments:
c		 This routine updates the vertex and adjacency lists
c		 by removing the triangles in v_local. It is used in
c		 conjunction with the routine find_ctriangles to remove
c		 a node from a 2-D Delaunay triangulation.
c
c		 The algorithm updates vertices and neighbour using only 
c		 local loops, i.e. without a global loop over all triangles. 
c		 Only triangles in c_list are considered, because these are 
c		 the only ones which are circum_circles of the input node. 
c
c		 Assumes input list of vertices in anticlockwise order
c
c		 The arrays must be dimensioned in the following way:
c
c		 integer tlist(nnpn_max) 
c		 logical c_list(2*nnpn_max)    
c		 integer v_local(3,2*nnpn_max)
c		 integer n_local(3,2*nnpn_max) 
c
c		 It is assumed that the maximum size of the local arrays
c		 v_local(3,*), n_local(3,*) and c_list(*) 
c                (which store local triangles) is 2*nnpn_max, 
c		 where nnpn_max is the maximum number of 
c		 neighbours per node. (Note: a strict upper bound on 
c		 the number of triangles attached to a node is 2*nnpn_max-2) 
c
c		 Only uses local loops
c
c                Calls no other routines.
c
c					M. Sambridge, RSES, Dec. 1995.
c
c------------------------------------------------------------------------
c
	Subroutine update_vn_lists 
     &      (node,np,nt,vertices,neighbour,tlist,
     &       on_hull,nnpn_max,np_local,
     &       nt_local,v_local,n_local,c_list)

	integer		vertices(3,*)
	integer		neighbour(3,*)
	integer		v_local(3,*)
	integer		n_local(3,*)
        logical         c_list(*)
        integer         tlist(*)
        integer		c1(3)
        integer		c2(3)
        integer		t
        integer		edg
        logical         on_hull
        data		c1/2,3,1/
        data		c2/3,1,2/
        logical         debug_mode
        common		/debug_nn_remove/debug_mode
c						write out v and n lists
        if(debug_mode)then
           write(*,*)' '
           write(*,*)' Inside update_vn_lists'
           write(*,*)' '
           write(*,*)' node = ',node
           write(*,*)' nt =',nt
           do j=1,nt
	      write(*,*)' v:',(vertices(k,j),k=1,3),
     &                  ' n:',(neighbour(k,j),k=1,3)
           end do
           write(*,*)' '
           write(*,*)' nt_local =',nt_local
           do j=1,nt_local
	      write(*,*)' vl:',(v_local(k,j),k=1,3),
     &                  ' nl:',(n_local(k,j),k=1,3)
           end do
           write(*,*)' '
        end if
c
c						renumber local triangles 
c						to negative values to avoid
c						confusion with global triangles
c						having the same number
	do i=1,nt_local
           n_local(1,i) = -n_local(1,i)
           n_local(2,i) = -n_local(2,i)
           n_local(3,i) = -n_local(3,i)
        end do

c
c						Procedure to update neighbour
c						and vertices lists in 3 stages
c						
c						stage 1:
c						update n_local and neighbour
c						for exterior edges of
c						local triangulation
c						(using negative local indices)
        j1 = 1
        j2 = np_local
        if(on_hull)j2 = np_local-1
        nt_old = j2-j1+1
        if(debug_mode)
     &  write(*,*)' nt_old=',nt_old,' nt_local =',nt_local

	do i = 1,nt_local
           if(c_list(i))then
              do j=1,3
                 t = -n_local(j,i) 
                 if(t.eq.0.or..not.c_list(t))then
                   n1 = v_local(c1(j),i)
                   n2 = v_local(c2(j),i)
                   do ip = j1,j2
                      it = tlist(ip)
                      do jt = 1,3
                         m1 = vertices(c1(jt),it)
                         m2 = vertices(c2(jt),it)
                         if(m1.eq.n1.and.m2.eq.n2)then
                            nn = neighbour(jt,it)
                            n_local(j,i) = nn
                            if(nn.ne.0)then
                               if(debug_mode)
     &                         write(*,*)' found a match',m1,m2,n1,n2,
     &                         ' local t',i,' old t',it,
     &                         ' neighbour',n_local(j,i)
                               neighbour(edg(nn,it,neighbour),nn) = -i
                            end if
                         else if(m1.eq.n2.and.m2.eq.n1)then
                            nn = neighbour(jt,it)
                            n_local(j,i) = nn
                            if(nn.ne.0)then
                               if(debug_mode)
     &                         write(*,*)' found a match',m1,m2,n1,n2,
     &                         ' local t',i,' old t',it,
     &                         ' neighbour',n_local(j,i)
                               neighbour(edg(nn,it,neighbour),nn) = -i
                            end if
                         end if
                      end do
                   end do
                 end if
              end do
           end if
        end do
c						for points on the hull
c						set to zero the adjacency
c						entries for all neighbours of
c						triangles attached to input 
c						node that are not already
c						attached to new triangles 
        if(on_hull)then
           do i=1,np_local-1 
              t = tlist(i)
              j = neighbour(edg(t,node,vertices),t)
              do k=1,3
                 if(neighbour(k,j).eq.t)then
                    if(debug_mode)
     &              write(*,*)' reseting triangle ',j,
     &              ' neighbour',t,' to zero' 
                    neighbour(k,j)=0
                    go to 6
                 end if
              end do
 6            continue
           end do
        end if

	if(debug_mode)then
           write(*,*)' after stage 1'
           do j=1,nt
	      write(*,*)' v:',(vertices(m,j),m=1,3),
     &                  ' n:',(neighbour(m,j),m=1,3)
           end do
           do j=1,nt_local
	      write(*,*)' vl:',(v_local(m,j),m=1,3),
     &                  ' nl:',(n_local(m,j),m=1,3)
           end do
        end if
c						stage 2:
c						replace old triangles with new			
        k = j1
	do i = 1,nt_local
           if(c_list(i))then
              t = tlist(k)
              k = k + 1
              if(debug_mode)
     &        write(*,*)' triangle',-i,' becomes triangle',t
              do j=1,3
                 vertices(j,t) = v_local(j,i)
                 neighbour(j,t) = n_local(j,i)
                 ii = n_local(j,i)
                 if(ii.lt.0)then
                    n_local(edg(-ii,-i,n_local),-ii) = t
                 else if(ii.gt.0)then
                    neighbour(edg(ii,-i,neighbour),ii) = t
                 end if
              end do
            end if
        end do
        nt_done = k-1
        nt_rem = nt_old-k+j1

        if(debug_mode)then
           write(*,*)' Number of triangles removed =',
     &     k-j1,' nt_old=',nt_old
           write(*,*)' Number remaining =',nt_rem
	   write(*,*)' after stage 2'
           do j=1,nt
	      write(*,*)' v:',(vertices(m,j),m=1,3),
     &                  ' n:',(neighbour(m,j),m=1,3)
           end do
        end if
 
c						stage 3:
c						swap last triangles (nt,nt-1,..)
c						into remaining vacant triangle 
c						positions 
c						(but do not do so if they are
c						 already in the vacant list)
        k1 = k
        k2 = nt_rem+(k1-1)
        k3 = nt-nt_rem
        do 10 i=1,nt_rem
           it = nt-(i-1)
c						determine if triangle to be
c						moved is already in vacant list 
c						If so then ignore
           do j=k1,k2
              t = tlist(j)
              if(it.eq.t)go to 10
           end do
c						If not then find 
c						next vacant slot
           do j = k,k2
              t = tlist(j)
              if(t.le.k3)then
                 if(debug_mode)
     &           write(*,*)' triangle',it,' becomes triangle',t
                 do m=1,3
                    vertices(m,t) = vertices(m,it)
                    neighbour(m,t) = neighbour(m,it)
                    ii = neighbour(m,it)
                    if(ii.ne.0)neighbour(edg(ii,it,neighbour),ii) = t
                 end do
                 k = j+1
                 go to 10
              end if
           end do
           write(*,*)' Error in subroutine update_vn_lists'
           write(*,*)' Error no vacant slots'
           write(*,*)' t =',t,' k=',k,' k2=',k2
           stop
 10     continue

c						update number of triangles
        nt = it-1
c						write out v and n lists

        if(debug_mode)then
           write(*,*)' after stage 3: New nt = ',nt
           do j=1,nt
	      write(*,*)' v:',(vertices(m,j),m=1,3),
     &                  ' n:',(neighbour(m,j),m=1,3)
           end do
           write(*,*)' '
        end if

	return
	end
c
c------------------------------------------------------------------------
c
c       circum_centre - calculates circum-centre of triangle t 
c
c
c       Input:
c               t                       triangle 
c               points(2,np)            array of node points
c               vertices(3,nt)          array of triangle vertices
c
c       Output:
c               centre(3)               centres(j,i) (j=1,2) contains the
c                                       co-ordinates of the centre of
c                                       circumcircle about Delaunay
c                                       triangle i, (i=1,...,nt),
c                                       j=3 contains squared radius of circle.
c       Comments:
c
c                No calls to other routines.
c
c                                       M. Sambridge, RSES, April 1994.
c
c------------------------------------------------------------------------
c
        Subroutine circum_centre(t,vertices,points,centre)
c
        real*8          points(2,*)
        real*8          centre(3)
        real*8          x1,x2,x3,y1,y2,y3,x,y
        real*8          dx2m1,dx2p1,dy2m1,dy2p1
        real*8          dx3m1,dx3p1,dy3m1,dy3p1
        integer         vertices(3,*)
        integer         t
c                                               Find centres of 
c                                               Delaunay circumcircle
        x1 = points(1,vertices(1,t))
        x2 = points(1,vertices(2,t))
        x3 = points(1,vertices(3,t))
        y1 = points(2,vertices(1,t))
        y2 = points(2,vertices(2,t))
        y3 = points(2,vertices(3,t))
 
        dx2m1 = x2-x1
        dx2p1 = x2+x1
        dy2m1 = y2-y1
        dy2p1 = y2+y1
        dx3m1 = x3-x1
        dx3p1 = x3+x1
        dy3m1 = y3-y1
        dy3p1 = y3+y1
        x = ((dx2m1*dx2p1 + dy2m1*dy2p1)*dy3m1*0.5d0
     &      -(dx3m1*dx3p1 + dy3m1*dy3p1)*0.5d0*dy2m1)/
     &      (dx2m1*dy3m1-dx3m1*dy2m1)
 
        y = (dx2m1*(dx3m1*dx3p1 + dy3m1*dy3p1)*0.5d0
     &      -dx3m1*(dx2m1*dx2p1 + dy2m1*dy2p1)*0.5d0)/
     &      (dx2m1*dy3m1-dx3m1*dy2m1)
 
        centre(1) = x
        centre(2) = y
        x1 = x - x1
        y1 = y - y1
        centre(3) = x1*x1 + y1*y1
 
        return
        end
c
c------------------------------------------------------------------------
c
c	check_nnpn_max - performs an error check on the size of nnpn_max
c
c	Input:	
c	        numn			number of neighbours	        
c	        node			node index 
c		nnpn_max		maximum number of neighbours for 
c					all nodes
c	Output:
c
c	Comments:
c
c		 No calls to other routines.
c
c					M. Sambridge, RSES, Jan. 1996.
c
c------------------------------------------------------------------------
c
	Subroutine check_nnpn_max(numn,node,nnpn_max)

	if(numn.gt.nnpn_max)then
           write(*,*)' '
           write(*,*)' Error - maximum number of neighbours'
           write(*,*)'         per node is not large enough'
           write(*,*)' '
           write(*,*)' Node',node,' has ',numn,' neighbours'
           write(*,*)' nnpn_max currently set to',nnpn_max
           write(*,*)' Work arrays are not big enough'
           write(*,*)' '
           write(*,*)' Remedy - increase size of parameter'
           write(*,*)'          nnpn_max in main program'
           write(*,*)' '
           stop
        end if

	return
	end
c
