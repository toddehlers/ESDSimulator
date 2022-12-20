c------------------------------------------------------------------------
c
c	nn2d_setup - Performs all setup procedures for natural neighbour 
c		     interpolation routine nn2D.
c
c	Input:
c	        np			number of nodes
c	        nt_max			maximum number of triangles
c	        nh_max			maximum number of triangles on convex 
c					hull (max size of array hulltriangles)
c		np_max			maximum number of nodes 
c		nnpn_max		maximum number of neighbours per node
c					(depends on the point distribution,
c                                        set to ~20 in calling program)
c		nmax			maximum sum of the number of neighbours 
c					per node (should be set to 
c					3*nt_max + np_max in calling program)
c		points(2,np)		array of node co-ordinates
c	        dmode			Delaunay calculation mode (integer)
c	        nmode			NN setup mode
c	        clockwise		logical for the vertices input order  
c		data(np)		data values at each node
c		nnn			Integer work array used by build_nv
c		nnlist			Integer work array used by build_nv
c		ntrilist		Integer work array used by build_nv
c		nohalt_hull		determines error response in routine
c					calculate_hulltriangles
c		eps 			tolerance parameter used by delaun
c					(see delaun for details)
c		vis_tlist		Integer work array used by delaun
c		vis_elist		Integer work array used by delaun
c		add_tlist		Integer work array used by delaun
c		nv_max			size of delaun work arrays vis_tlist
c					vis_elist, & add_tlist (passed to 
c					delaun for error checking)
c
c	Output:
c	        nt			number of triangles
c	        vertices(3,nt)		array of triangle vertices 
c               centres(3,nt)           centres(j,i) (j=1,2) contains the
c                                       co-ordinates of the centre of 
c                                       circumcircle about Delaunay 
c                                       triangle i, (i=1,...,nt),
c                                       j=3 contains squared radius of circle.
c	        neighbour(3,nt)		array of neighbouring triangles.	
c					Neighbour(i,j) is the triangle
c					opposite node i in triangle j,
c					stored counterclockwise about j.
c	        nh			number of triangles with an edge
c					on the convex hull
c		hulltriangles(nh)	array of triangles with an edge 
c					on the convex hull
c		loc			an initial guess triangle for point 
c					location routine `Triloc' used by nn2D
c					(set somewhere near the centre)
c
c	Operation modes:
c
c		 The setup routine will perform different tasks depending
c		 on the input parameters dmode and nmode (see table below).
c		 Depending on the modes used some work arrays may be set 
c		 to size 1 to save memory. The "Memory Savings" column in the
c		 table below shows the dimension statement that may
c		 be used in the calling program if the routine is ONLY EVER 
c		 CALLED IN THE CORRESPONDING MODE. 
c
c		 PARAMETERS	ACTION	 		MEMORY SAVINGS
c
c		 nmode = 1	Delaunay only 		real*8 centres(3,1)
c							integer hulltriangles(1)
c		 nmode = 0	Delaunay + nn setup	
c		 nmode = -1	nn setup only		
c
c		 dmode > 0	Delaunay read in from   integer vis_tlist(1)
c				logical unit dmode.	integer vis_elist(1) 
c							integer add_tlist(1) 
c		 dmode = 0      Qhull used		Same as dmode > 0.
c		 dmode = -1	Delaun + X-sort         integer nnn(1)
c							integer nnlist(1)
c							integer ntrilist(1)
c		 dmode = -2	Delaun + no sort        Same as dmode=-1
c
c		 dmode = 0 & nmode=1			integer neighbour(3,1)
c		 
c		 A call with nmode = -1 can only be made after a call 
c		 with nmode = 1.
c
c	Comments:
c
c		 If the arrays are used then they should be dimensioned
c		 in the calling program in the following way:
c
c		 real*8		points(2,np_max)
c		 real*8		centres(3,nt_max)
c		 integer	vertices(3,nt_max)
c		 integer	neighbour(3,nt_max)
c		 integer	hulltriangles(nh_max)
c		 
c		 Note: nh_max can usually be dimensioned much less than nt_max
c		 because hulltriangles, stores in a compact form, all 
c		 triangles with an edge on the convex hull. Except for
c		 very irregular point distributions nh << nt. If nh is 
c		 determined to be > nh_max then an error is reported and the
c		 program is halted (unless nohalt parameter is set to 1). 
c		 The array hulltriangles is only used by nn2Do see routine 
c		 calculate_hulltriangles. If nh_max = 1 then hulltriangles 
c		 is not calculated.
c
c		 The initial guess triangle 'loc' is set in nn_setup but
c		 at each call it will be set to the triangle found during
c		 the previous call to nn2D. The user may modify its value
c		 if the input point (x,y) is known to be in, or near, a
c		 particular triangle.
c
c		 If dmode > 0 the the deluanay tessellation is read in from
c		 logical unit `dmode' instead of being calculated internally 
c		 This can be useful if qhullf fails because
c		 of precision errors. The Deluanay may be determined
c		 externally to this program using a double precision version
c		 or another algorithm, e.g. Fortune's sweepline method.
c
c		 If 50 > dmode > 0 then:
c		 It is assumed that the read in format has one triangle per
c		 line represented as a triplet of nodes numbered from ZERO, 
c		 which is the standard output format of codes qhull 
c		 (quickhull method) and voronoi (sweepline method).
c		 If clockwise = .true. (.false.) then the vertices are assumed 
c		 to be in clockwise (anti-clockwise) order. Note program
c		 qhull outputs vertices in anti-clockwise order while 
c		 voronoi in clockwise order. The internal format is 
c		 anti-clockwise and nodes numbered from ONE.
c
c		 If dmode => 50 then:
c		 It is assumed that the read in format has one triangle per
c		 line represented as a triplet of nodes numbered from ONE,
c		 which is the output format of program del (using delaun). 
c
c		 Three other work arrays are produced as a `by product'
c		 of the routine build_nv which calculates the neighbour
c		 array. These must be dimensioned in the calling program in 
c		 the following way (unless delaun is used for calculating the
c		 Delaunay because it already determines the neighbour array)
c
c		 integer nnn(np_max+1)  : number of neighbours per node
c		 integer nnlist(nmax)   : natural neighbours per node
c		 integer ntrilist(nmax) : triangles attached to each node 
c
c		 The value of nmax should be set to (3*nt_max + np_max)
c		 in the calling program.
c
c		 Each of these are useful lists that describe features of
c		 the Voronoi diagram. Both nnlist and ntrilist are stored in
c		 a compact format to avoid zeros. They are only used 
c		 in the setup routine and the memory may be freed once
c		 initialization is completed.
c		 
c
c		 Calls are made to: qhullf, ccentres, build_nv and 
c				    calculate_hulltriangles, delaun.
c		 
c					M. Sambridge, RSES, April 1994.
c  					        (Last modified 10/4/96)
c
c------------------------------------------------------------------------
c
	Subroutine nn2d_setup
     &             (np,nt_max,nh_max,np_max,nnpn_max,nmax,
     &              points,dmode,nmode,clockwise,data,nt,vertices,
     &              centres,neighbour,nh,hulltriangles,nohalt_hull,
     &              loc,nnn,nnlist,ntrilist,
     &              eps,nv_max,vis_tlist,vis_elist,add_tlist,
     &		    lt_work,ln_work)

	real*8		points(2,*)
	real*8		centres(3,*)
	real*8		data(*)
	real*8		eps
	integer		vertices(3,*)
	integer		neighbour(3,*)
        integer		hulltriangles(*)
	integer		nnn(*)
	integer		nnlist(*)
	integer		ntrilist(*)
        integer         vis_tlist(*)
        integer         vis_elist(*)
        integer         add_tlist(*)
        integer         dmode,nmode
	logical*1	lt_work(*)
	logical*1	ln_work(*)
	logical		nnwrite
	logical		clockwise
	logical		ldummy

	common/nnswitches/nnwrite,lud

        if(nmode.eq.1.or.nmode.eq.0)then

           if(dmode.eq.0)then
c                                       calculate Delaunay using qhull 
 
              call qhullf(np,2,2,nt_max,0,points,nt,vertices)

           else if(dmode.eq.-1.or.dmode.eq.-2)then

c					sort the points in ascending x order
	      if(dmode.eq.-1)then

c		  write(*,*)' X sort in progress'
	          call hpsort_d(np,1,points)
c		  write(*,*)' X sort done'

	      end if

c                                       calculate Delaunay using delaun 
 
              call delaun (points,np,neighbour,vertices,nt,nt_max,
     &                     vis_tlist,vis_elist,add_tlist,eps,nv_max,
     &                     0,ldummy,0,0,0)

	   else
c                                       read in Delaunay vertices

              nt = 0
              i1 = 1
              i2 = 2
              if(clockwise)then
                 i1 = 2
                 i2 = 1
              end if
              read(dmode,*)
  1           read(dmode,*,end=3,err=2)
     &        vertices(i1,nt+1),vertices(i2,nt+1),vertices(3,nt+1)
              nt = nt + 1
              if(nt.ge.nt_max)then
                 write(*,*) 'Error in nn_setup: too many triangles'
                 write(*,*) 'Remedy: increase size of parameter nt_max'
                 write(*,*) '        in calling program.'
                 stop 
              end if
              go to 1
  2           write(*,*)
     &        'Error in nn_setup: read error in Delaunay input file'
              stop
  3           continue
     
           end if
 
c					adjust array vertices to
c					range from nodes 1 to np
c
	   if(dmode.ge.0.and.dmode.lt.50)then
	      do 5 i = 1,nt
	         vertices(1,i) = vertices(1,i) + 1
	         vertices(2,i) = vertices(2,i) + 1
	         vertices(3,i) = vertices(3,i) + 1
 5            continue
	   end if

	end if
c
c					Perform set up for nn interpolation
c
        if(nmode.eq.0.or.nmode.eq.-1)then

c					set initial guess for 
c					triangle location procedure
           loc = nt/2

c                                       Calculate Circumcentres

           call ccentres(points,vertices,nt,centres)

c                                       Build neighbour matrix
c					(if not already built)

           if(dmode.ge.0)then
              call build_nv
     &        (np,vertices,nt,np_max,nmax,
     &         neighbour,nnn,nnlist,ntrilist)
           end if

c					calculate hulltriangles

           if(nh_max.gt.1) call calculate_hulltriangles
     &     (neighbour,nt,nh_max,hulltriangles,nh,nohalt_hull)

c					initialize logical work arrays 
           do i=1,nt
              lt_work(i) = .false.
           end do
           do i=1,np
              ln_work(i) = .false.
           end do

	end if

	return
	end
c
c------------------------------------------------------------------------
c
c	nn2D - calculates natural neighbour interpolation at point x,y
c
c	Input:
c		x,y			co-ordinates of input point
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
c		hulltriangles(nh)	array of triangles with an edge 
c					on the convex hull
c	        nh			number of triangles with an edge
c					on the convex hull
c		loc			an initial guess triangle for point 
c					location routine `Triloc'
c		data			data values at each node
c		nnext			logical if true then extend nn outside
c					of convex hull
c	        int_method		determines the type of interpolation
c					method. 0= Watson's method
c						1= recursive method
c						2= recursive method + df(2)
c						3= Linear interpolation + df(2)
c						4= closed formula method 
c						5= same as 4 + df(2)
c						6= same as 4 + df(2) + ddf(4)
c		nnpn_max		maximum number of neighbours per node
c					(depends on the point distribution,
c                                        set to ~50 in calling program)
c		work...                 work arrays of the size of 
c					nnpn_max.
c		lt_work			logical work array of size nt_max 
c		ln_work			logical work array of size np_max 
c
c	Output:
c		f			interpolated data value
c		df(1)			interpolated derivative df/dx
c		df(2)			interpolated derivative df/dy
c		v			area of voronoi cell about (x,y)
c		out			logical determining if (x,y) is i
c					inside the convex hull of points
c
c	Comments:
c		 The area of the voronoi cell about (x,y) can be used 
c		 as an inverse measure of the quality of the interpolation.
c		 If (x,y) is outside of the hull the voronoi cell is
c		 unbounded and v is the sum of the natural neighbour
c		 co-ordinates calculated with Watson's method.
c
c		 Note: plot_tc and plot_c are user supplied routines
c		 that can be dummy routines, or used with a plot
c		 program for plotting circum-triangles and circum-circles 
c		 found during the calculation. The file nnplot.f contains
c		 a set of dummy routines to be compiled with nn.f. 
c		 Alternatively all references to plot_tc and plot_c may
c		 commented out.
c
c		 A recursive method is available, as an alternative to
c		 Watson's method, for calculation of natural neighbour
c		 co-ordinates. The action of nn_setup is the same for 
c		 either method although the parameter nnpn_max is also 
c		 used to size arrays used in the recursive mode.
c
c		 Calls are made to: Triloc, nn2Dr, nn2Drd, nn2Do or nn2Di
c				    nn2Df, nn2Dfd, nn2Dfdd.
c
c					M. Sambridge, RSES, April 1994.
c
c------------------------------------------------------------------------
c
	Subroutine nn2D
     &  (x,y,points,vertices,neighbour,nnext,int_method,centres,
     &   hulltriangles,nh,loc,data,nnpn_max,work_r1,work_r2,
     &   work_d1,work_i1,work_i2,work_i3,work_d2,work_d3,work_d4,
     &   work_d5,work_d6,work_d7,lt_work,ln_work,out,f,df,ddf,v)

	real*8		points(2,*)
	real*8		centres(3,*)
	real*8		xd(2),x,y
	real*8		data(*)
	real*8		f,df(2),ddf(3),v
	integer		vertices(3,*)
	integer		neighbour(3,*)
        integer		hulltriangles(*)
	integer		int_method
	logical		out
	logical		nnext
	logical		nnwrite

        real*8          work_d1(2,nnpn_max)
        real*8          work_d2(2,nnpn_max)
        real*8          work_d3(2,nnpn_max)
        real*8          work_d4(2,nnpn_max)
        real*8          work_d5(2,nnpn_max)
        real*8          work_d6(2,nnpn_max)
        real*8          work_d7(2,nnpn_max)
        real            work_r1(nnpn_max,2)
        real            work_r2(nnpn_max)
	integer		work_i1(nnpn_max)
        integer         work_i2(2,nnpn_max)
        integer         work_i3(nnpn_max)
	logical*1	lt_work(*)
	logical*1	ln_work(*)

	common/nnswitches/nnwrite,lud
 
        xd(1) = x
        xd(2) = y
c					Use triangle walking routine
c					to locate triangle of input
c					point

        call Triloc(x,y,points,vertices,neighbour,loc,out)
c
c					If point is outside convex hull
c					call routine nn2Do

	if(out)then


 	   if(nnext.and.int_method.eq.0)then
	        if(nnwrite)write(*,*)' calling nn2Do'
                call nn2Do
     &               (x,y,points,vertices,neighbour,
     &                centres,hulltriangles,nh,data,f,df,v)
           else
               f = 0.d0
               df(1) = 0.d0
               df(2) = 0.d0
               v = 0.d0
	   end if


	else

c					use modified Watson's method 
	   if(int_method.eq.0)then

 	      if(nnwrite)write(*,100)

              call nn2Di
     &             (x,y,points,vertices,neighbour,
     &              centres,loc,data,f,df,v)

           else if(int_method.eq.1)then

 	      if(nnwrite)write(*,200)
c					use recursive method
              call nn2Dr
     &             (x,y,points,vertices,neighbour,centres,
     &              loc,data,nnpn_max,work_r1,work_r2,work_d1,
     &              work_i1,work_i2,work_i3,f,df,v)


           else if(int_method.eq.2)then

 	      if(nnwrite)write(*,300)
c					use recursive method and 
c					calculate derivatives
              call nn2Drd
     &             (x,y,points,vertices,neighbour,centres,
     &              loc,data,nnpn_max,work_r1,work_r2,work_d1,
     &              work_i1,work_i2,work_i3,f,df,v)

           else if(int_method.eq.3)then

 	      if(nnwrite)write(*,400)
c					use linear interpolation in triangles
	      call nn2DL
     &             (x,y,points,vertices,loc,data,f,df,v)

           else if(int_method.eq.4)then

 	      if(nnwrite)write(*,500)
c					use closed formula method
c					no derivatives
c
              call nn2Df
     &             (xd,points,vertices,neighbour,centres,
     &              loc,data,nnpn_max,work_d1,work_d2,work_i1,
     &              work_i2,work_r2,work_i3,lt_work,ln_work,f,v)

              df(1) = 0.d0
              df(2) = 0.d0

           else if(int_method.eq.5)then

 	      if(nnwrite)write(*,600)
c					use closed formula method
c					with 1st derivatives
c
              call nn2Dfd
     &             (xd,points,vertices,neighbour,centres,
     &              loc,data,nnpn_max,work_d1,work_d2,work_i1,
     &              work_i2,work_r2,work_i3,work_d3,work_d4,
     &              lt_work,ln_work,f,df,v)

           else if(int_method.eq.6)then

 	      if(nnwrite)write(*,600)
c					use closed formula method
c					with 1st and 2nd derivatives
c
              call nn2Dfdd
     &             (xd,points,vertices,neighbour,centres,
     &              loc,data,nnpn_max,work_d1,work_d2,work_i1,
     &              work_i2,work_r2,work_i3,work_d3,work_d4,
     &	            work_d5,work_d6,work_d7,
     &              lt_work,ln_work,f,df,ddf,v)

           else

              write(*,*)' '
              write(*,*)
     &        ' Error in nn2D: interpolation method not defined'
              write(*,*)'                int_method set to ',int_method
              write(*,*)' '
              write(*,*)' Valid options are:'
              write(*,*)'   0=NN-Watson'
              write(*,*)'   1=NN-recursive'
              write(*,*)'   2=NN-recursive+derivatives'
              write(*,*)'   3=Linear interpolation'
              write(*,*)'   4=NN-closed formula'
              write(*,*)'   5=NN-closed formula + df(2)'
              write(*,*)'   6=NN-closed formula + df(2) + ddf(3)'
              stop


           end if

	end if

c	if(nnwrite)then
c	   write(*,*)' Overlap voronoi area = ',v
c	   write(*,*)' Normalized value of f   =',f
c	   write(*,*)' Normalized value of dfx =',df(1)
c	   write(*,*)' Normalized value of dfy =',df(2)
c	end if

 100 	format(' calling nn2Di:   Watson method')
 200 	format(' calling nn2Dr:   recursive method')
 300 	format(' calling nn2Drd:  recursive method with derivatives')
 400 	format(' calling nn2DL:   Linear interpolation')
 500 	format(' calling nn2Df:   Closed formulae method')
 600 	format(' calling nn2Dfd:   Closed formulae method',
     &         'with first derivatives')
 700 	format(' calling nn2Dfdd:   Closed formulae method',
     &         'with first and second derivatives')

	return
	end
c
c------------------------------------------------------------------------
c
c       nn2DL - calculates linear triangle interpolation at point x,y
c               (This routine is an alternative to nn2D)
c
c       Input:
c               x,y                     co-ordinates of input points
c               points(2,np)            array of node points
c               vertices(3,nt)          array of triangle vertices
c               loc                     an initial guess triangle for point
c                                       location routine `Triloc'
c               data                    data values at each node
c
c       Output:
c               f                       interpolated data value
c               df(1)                   interpolated derivative df/dx
c               df(2)                   interpolated derivative df/dy
c               v                       area of voronoi cell about (x,y)
c
c       Comments:
c
c                Performs simple linear interpolation in triangles
c
c                                       M. Sambridge, RSES, April 1994.
c
c------------------------------------------------------------------------
c
        Subroutine nn2DL
     &  (x,y,points,vertices,loc,data,f,df,v)

        real*8          points(2,*)
        real*8          x,y
        real*8          x1,y1,x2,y2,x3,y3,A,a1,a2,a3
        real*8          dx12,dx13,dy12,dy13,df12,df13
        real*8          data(*)
        real*8          f,df(2),v
        integer         vertices(3,*)


c                                       calculate area of triangle loc
        i = vertices(1,loc)
        j = vertices(2,loc)
        k = vertices(3,loc)
        x1 = points(1,i)
        y1 = points(2,i)
        x2 = points(1,j)
        y2 = points(2,j)
        x3 = points(1,k)
        y3 = points(2,k)
        A = ((x2-x1)*(y3-y1) - (x3-x1)*(y2-y1))
        if(A.lt.0.0)
     &  write(6,*)' Error in nn2DL: negative area for T=',loc

c                                       calculate f and derivatives
        df12 = data(i) - data(j)
        df13 = data(i) - data(k)
        dx12 = x1 - x2
        dx13 = x1 - x3
        dy12 = y1 - y2
        dy13 = y1 - y3
        a1 = (df12*dy13 - df13*dy12)/A
        a2 = (df13*dx12 - df12*dx13)/A
        a3 = data(i) - a1*x1 - a2*y1
        f = a1*x + a2*y + a3
        df(1) = a1
        df(2) = a2
        v = 0.d0
c       write(*,*)' New value of f =',f
c       write(*,*)' New value of dfx =',df(1)
c       write(*,*)' New value of dfy =',df(2)

        return
        end
c
c
c------------------------------------------------------------------------
c
c	build_nv - Builds neighbour array for Delaunay triangulation in 2-D.
c
c	Input:	
c	        vertices(3,nt)		array of triangle vertices	
c	        nt			number of triangles
c		np_max			maximum number of nodes
c		nmax			maximum total number of neighbours 
c					per node (should be set to 
c					3*nt_max + np_max in calling program)
c
c	Output:
c		neighbour(3,nt)		array of neighbouring triangles
c
c	Comments:
c		 Assumes input list of vertices in anticlockwise sequence
c		 and produces an anticlockwise list of neighbour triangles.
c		 The value of neighbour(i,j) is the index of the neighbouring
c		 triangle opposite node i in triangle j.
c
c		 Three temporary work arrays are used and must be dimensioned
c		 in the calling program in the following way:
c
c		 integer nnn(np_max+1)  : number of neighbours per node
c		 integer nnlist(nmax)   : natural neighbours per node
c		 integer ntrilist(nmax) : triangles attached to node 
c
c		 The value of nmax should be set to (3*nt_max + np_max)
c		 in the calling program.
c
c		 No calls to other routines.
c
c					M. Sambridge, RSES, April 1994.
c					(using ideas by J.Braun)
c
c------------------------------------------------------------------------
c
	Subroutine build_nv
     &             (np,vertices,nt,np_max,nmax,
     &              neighbour,nnn,nnlist,ntrilist)
c
	integer		vertices(3,*)
	integer		neighbour(3,*)
	integer		nnn(*)
	integer		nnlist(*)
	integer		ntrilist(*)
	logical		nnwrite

	common/nnswitches/nnwrite,lud

	if(nnwrite)write(*,*)' Building neighbour v ...'
c
c					initialize neighbour list
	do 5 i = 1,3
	   do 5 j = 1,nt
	      neighbour(i,j) = 0
 5      continue
c					initialize work arrays
        do 6 i = 1,nmax
	   nnlist(i) = 0
	   ntrilist(i) = 0
 6      continue

        do 7 i = 1,np
	   nnn(i) = 0
 7      continue

	do 10 it = 1,nt
	   i1 = vertices(1,it)
	   i2 = vertices(2,it)
	   i3 = vertices(3,it)
	   nnn(i1) = nnn(i1) + 1
	   nnn(i2) = nnn(i2) + 1
	   nnn(i3) = nnn(i3) + 1
 10     continue

c					turn nnn into a running sum
	itemp = nnn(1)+1
	nnn(1) = 1
	do 20 j = 2,np+1
	   itemp2  = itemp 
	   itemp   = itemp + nnn(j)+1
	   nnn(j) = itemp2 + 1
 20     continue
c       write(*,*)' size of array =',nnn(np+1)-1
c       write(*,*)' 3nt+np        =',3*nt+np

	if(nnn(np+1).ge.nmax)then
           write(*,*)'Error: array sizes too small in subroutine '
     &               ,'build_nv'
           write(*,*)'       maximum number of neighbours for all nodes'
           write(*,*)'       is too small: current value =',nmax
           write(*,*)'       Increase size of parameter nmax'
           write(*,*)'       to at least',nnn(np+1)
           write(*,*)'       This will be satisfied if nmax is set'
           write(*,*)'       to 3*nt_max+np_max in calling program' 
	   stop
	end if

	do 25 it = 1,nt
	   i1 = vertices(1,it) 
	   i2 = vertices(2,it) 
	   i3 = vertices(3,it) 
c						compare neighbours i1 i2
c						(remove go to ?)
	   j1 = nnn(i1)
	   j2 = nnn(i1+1) - 1 
	   jt = 0
	   do 30 j = j1,j2
	      if(nnlist(j).eq.0)then
	         nnlist(j) = i2
	         ntrilist(j) = it
c						if we have recorded connection
c						then jump out of loop
	         go to 31
	      else if(nnlist(j).eq.i2.and.ntrilist(j).ne.it)then
                 jt = ntrilist(j)
	         go to 31
	      end if
  30       continue
  31       continue
c						if neighbours are found then
c						skip second loop 
	   if(jt.eq.0)then
	      j1 = nnn(i2)
	      j2 = nnn(i2+1) - 1 
	      do 32 j = j1,j2
	         if(nnlist(j).eq.0)then
	            nnlist(j) = i1
	            ntrilist(j) = it
c						if we have inserted connection
c						then jump out of loop
	            go to 33
	         end if
  32          continue
	   end if
  33       continue

	   if(jt.ne.0)then
c						found neighbours it,jt with
c						common nodes i1 and i2
	      neighbour(3,it) = jt
	      k1 = vertices(1,jt)
	      k2 = vertices(2,jt)
	      k3 = vertices(3,jt)
	      if(k1.ne.i1.and.k1.ne.i2)then 
	         neighbour(1,jt) = it
	      else if(k2.ne.i1.and.k2.ne.i2)then 
	         neighbour(2,jt) = it
	      else
	         neighbour(3,jt) = it
	      end if
           end if
c						compare neighbours i1 i3
	   jt = 0
	   j1 = nnn(i1)
	   j2 = nnn(i1+1) - 1 
	   do 130 j = j1,j2
	      if(nnlist(j).eq.0)then
	         nnlist(j) = i3
	         ntrilist(j) = it
c						if we have recorded connection
c						then jump out of loop
	         go to 131
	      else if(nnlist(j).eq.i3.and.ntrilist(j).ne.it)then
                 jt = ntrilist(j)
	         go to 131
	      end if
  130      continue
  131      continue
c						if neighbours are found then
c						skip second loop 
	   if(jt.eq.0)then
	      j1 = nnn(i3)
	      j2 = nnn(i3+1) - 1 
	      do 132 j = j1,j2
	         if(nnlist(j).eq.0)then
	            nnlist(j) = i1
	            ntrilist(j) = it
c						if we have inserted connection
c						then jump out of loop
	            go to 133
	         end if
  132         continue
	      end if
  133      continue
	   if(jt.ne.0)then
c						found neighbours it,jt with
c						common nodes i1 and i3
	     neighbour(2,it) = jt
	     k1 = vertices(1,jt) 
	     k2 = vertices(2,jt)
	     k3 = vertices(3,jt)
	     if(k1.ne.i1.and.k1.ne.i3)then 
	        neighbour(1,jt) = it
	     else if(k2.ne.i1.and.k2.ne.i3)then 
	        neighbour(2,jt) = it
	     else
	        neighbour(3,jt) = it
	     end if
           end if
c						compare neighbours i2 i3
	   jt = 0
	   j1 = nnn(i2)
	   j2 = nnn(i2+1) - 1 
	   do 230 j = j1,j2
	      if(nnlist(j).eq.0)then
	         nnlist(j) = i3
	         ntrilist(j) = it
c						if we have recorded connection
c						then jump out of loop
	         go to 231
 	      else if(nnlist(j).eq.i3.and.ntrilist(j).ne.it)then
                 jt = ntrilist(j)
 	         go to 231
	      end if
  230      continue
  231      continue
c						if neighbours are found then
c						skip second loop 
	   if(jt.eq.0)then
	      j1 = nnn(i3)
	      j2 = nnn(i3+1) - 1 
	      do 232 j = j1,j2
	         if(nnlist(j).eq.0)then
	            nnlist(j) = i2
	            ntrilist(j) = it
c						if we have inserted connection
c						then jump out of loop
	            go to 233
	         end if
  232         continue
	      end if
  233      continue
	   if(jt.ne.0)then
c						found neighbours it,jt with
c						common nodes i2 and i3
	     neighbour(1,it) = jt
	     k1 = vertices(1,jt) 
	     k2 = vertices(2,jt)
	     k3 = vertices(3,jt)
	     if(k1.ne.i2.and.k1.ne.i3)then 
	        neighbour(1,jt) = it
	     else if(k2.ne.i2.and.k2.ne.i3)then 
	        neighbour(2,jt) = it
	     else
	        neighbour(3,jt) = it
	     end if
           end if

 25     continue

	if(nnwrite)write(*,*)' built neighbour v'

	return
	end
c
c------------------------------------------------------------------------
c
c       Calculate_hulltriangles - finds all triangles with a face on
c                                 the convex hull by searching through
c                                 the entries in the array neighbour.
c
c       Input:
c               neighbour(3,nt)         array of neighbouring tetrahedra
c               nt                      number of tetrahedra
c	        nh_max			maximum number of triangles on convex 
c					hull (max size of array hulltriangles)
c		nohalt			determines error response
c       Output:
c		hulltriangles(nh)	array of triangles with an edge
c					on the convex hull
c               nh                      number of tetrahedra with an edge
c                                       on the convex hull
c       Comments:
c
c                This routine fills up the array hulltriangles which
c                is only used by routine nn2Do, i.e the `pseudo-extension' 
c		 Watson's nn-interpolation method to points outside of the 
c		 convex hull. If nnext is set to false then hulltriangles
c		 is never used and the array can be set to size 1.
c
c		 If nohalt = 0 then the routine will stop with an error
c		 message if nh > nh_max. If nohalt .ne. 0 and nh > nh_max
c		 then it will return nh = -1. 
c
c                No calls to other routines.
c
c                                       M. Sambridge, RSES, May 1995.
c
c------------------------------------------------------------------------
c
	Subroutine calculate_hulltriangles
     &             (neighbour,nt,nh_max,hulltriangles,nh,nohalt)
c
	integer		neighbour(3,*)
	integer		hulltriangles(*)

c                                               store list of triangles
c                                               which have an edge on the
c                                               convex hull.
c                                               (used by routine nn2D)
        nh = 1
        do 100 j = 1,nt
           if(neighbour(1,j).eq.0.or.
     &        neighbour(2,j).eq.0.or.
     &        neighbour(3,j).eq.0)then
              hulltriangles(nh) = j
              nh = nh + 1
              if(nh.gt.nh_max.and.nohalt.eq.0)then
                  write(*,*)' Error array storing outward facing '
                  write(*,*)' triangles on convex hull is too small.'
                  write(*,*)' Increase size of parameter nh_max'
                  stop
              else if(nh.gt.nh_max.and.nohalt.ne.0)then
                  nh = -1
                  return
              end if
           end if
 100    continue
        nh = nh -1
 
	return
	end
c
c------------------------------------------------------------------------
c
c	Triloc - locates the triangle containing point x,y
c
c	Input:
c		x,y			co-ordinates of input points	
c		points(2,np)		array of node points	
c	        vertices(3,nt)		array of triangle vertices	
c	        neighbour(3,nt)		array of neighbouring triangles.	
c					Neighbour(i,j) is the triangle
c					opposite node i in triangle j,
c					stored counterclockwise about j.
c	        loc			first guess of triangle containing
c					(x, y).
c
c	Output:
c	        loc			index of triangle containing 
c					input point.
c	        out			=true if (x,y) is outside of
c					the convex hull, otherwise = false. 
c
c	Comments:
c		 If (x,y) is outside convex hull loc is a `nearby' triangle
c		 on the hull and out is set to true.
c
c		 No calls to other routines.
c
c					M. Sambridge, RSES, April 1994.
c
c------------------------------------------------------------------------
c
	Subroutine Triloc(x,y,points,vertices,neighbour,loc,out)
c
	real*8		points(2,*)
	integer		vertices(3,*)
	integer		neighbour(3,*)
	integer		p1,p2
        logical		out
	real*8		x,y,del1,del2
        integer		c1(3)
        data		c1/2,3,1/
c
	out = .false.

 10     continue
c					point is outside convex hull
        if(out)return

        do 20 i=1,3
	   j = c1(i)
	   k = c1(j)
           p1 = vertices(i,loc)
           p2 = vertices(j,loc)
	   del1 = (points(2,p1)-y)*(points(1,p2)-x)
	   del2 = (points(1,p1)-x)*(points(2,p2)-y)
	   if(del1.gt.del2)then
	      if(neighbour(k,loc).eq.0)then
                 out = .true.
	      else
	         loc = neighbour(k,loc)
	      end if
	      go to 10
	   end if
 20     continue
	
	return
	end
c
c------------------------------------------------------------------------
c
c	Ccentres - calculates centres of all Delaunay circumcircles
c
c
c	Input:
c		points(2,np)		array of node points	
c	        vertices(3,nt)		array of triangle vertices	
c	        nt			number of triangles
c
c	Output:
c               centres(3,nt)           centres(j,i) (j=1,2) contains the
c                                       co-ordinates of the centre of 
c                                       circumcircle about Delaunay 
c                                       triangle i, (i=1,...,nt),
c                                       j=3 contains squared radius of circle.
c
c	Comments:
c
c		 No calls to other routines.
c
c					M. Sambridge, RSES, April 1994.
c
c------------------------------------------------------------------------
c
	Subroutine ccentres(points,vertices,nt,centres)
c
	real*8		points(2,*)
	real*8		centres(3,*)
	real*8		x1,x2,x3,y1,y2,y3,x,y
	real*8		dx2m1,dx2p1,dy2m1,dy2p1
	real*8		dx3m1,dx3p1,dy3m1,dy3p1
	real*8		denom
	integer		vertices(3,*)
c						Find centres of all
c						Delaunay Circumcircles
	do 5 i= 1,nt

	   x1 = points(1,vertices(1,i))
	   x2 = points(1,vertices(2,i))
	   x3 = points(1,vertices(3,i))
	   y1 = points(2,vertices(1,i))
	   y2 = points(2,vertices(2,i))
	   y3 = points(2,vertices(3,i))

           dx2m1 = x2-x1
           dx2p1 = x2+x1
           dy2m1 = y2-y1
           dy2p1 = y2+y1
           dx3m1 = x3-x1
           dx3p1 = x3+x1
           dy3m1 = y3-y1
           dy3p1 = y3+y1
           denom = dx2m1*dy3m1-dx3m1*dy2m1
	   x = ((dx2m1*dx2p1 + dy2m1*dy2p1)*dy3m1*0.5d0
     &         -(dx3m1*dx3p1 + dy3m1*dy3p1)*0.5d0*dy2m1)/
     &         (denom)

	   y = (dx2m1*(dx3m1*dx3p1 + dy3m1*dy3p1)*0.5d0 
     &         -dx3m1*(dx2m1*dx2p1 + dy2m1*dy2p1)*0.5d0)/ 
     &         (denom)

	   centres(1,i) = x
	   centres(2,i) = y
           x1 = x - x1
           y1 = y - y1
	   centres(3,i) = x1*x1 + y1*y1

 5	continue

	return
	end
c
c------------------------------------------------------------------------
c
c	nn2Di - calculates natural neighbour interpolation at point x,y
c		when point is inside convex hull 
c
c	Input:
c		x,y			co-ordinates of input point
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
c
c	Output:
c		f			interpolated data value
c		df(1)			interpolated derivative df/dx
c		df(2)			interpolated derivative df/dy
c		v			area of voronoi cell about (x,y)
c
c	Comments:
c		 On input point (x,y) must be in triangle loc. 
c		 Note: The input triangle contains (x,y) and so it's 
c		 circumcircle must also contain (x,y).
c
c		 See comments on plot routine plot_tc above.
c
c		 This routine uses Sloan's LIFO stack technique to 
c		 avoid an expensive global search over all triangles.
c		 (see Sambridge, Braun and McQueen, 1995; 
c		  Geophys. J. Int., vol 122, 837-857).
c
c		 Calls are made to: nn_tri. plot_tc, stackpairinit,
c				    poppair,pushpair & stackpairempty. 
c
c					M. Sambridge, RSES, April 1994.
c
c------------------------------------------------------------------------
c
	Subroutine nn2Di
     &             (x,y,points,vertices,neighbour,
     &              centres,loc,data,f,df,v)

	real*8		points(2,*)
	real*8		centres(3,*)
	real*8		data(*)
	real*8		x,y,dist,dx,dy
	real*8		f,df(2),v,dvx,dvy
	integer		vertices(3,*)
	integer		neighbour(3,*)
	integer		tri,pos
	logical		nnwrite

	common/nnswitches/nnwrite,lud

 	if(nnwrite)write(lud,*)' Result of nn2Di search:'

c					Initialize variables
	f = 0.d0
	df(1) = 0.d0
	df(2) = 0.d0
	v = 0.d0
	dvx = 0.d0
	dvy = 0.d0
c
c					Calculate natural neighbour
c					contribution from input triangle.
c
 	call nn_tri(x,y,points,vertices,centres,loc,
     &              data,f,df,v,dvx,dvy)

	if(nnwrite)write(lud,*)loc
c
c					Find all other circumcircles 
c					containing input point and update
c					natural neighbour contribution
c					from these triangles
c
c					plot circum-triangle and circum-circle

 	call plot_tc(loc,points,vertices,centres)

	   
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
	   if(neighbour(1,i).eq.loc)call pushpair(i,1)
	   if(neighbour(2,i).eq.loc)call pushpair(i,2)
	   if(neighbour(3,i).eq.loc)call pushpair(i,3)
	end if
	if(j.ne.0)then
	   if(neighbour(1,j).eq.loc)call pushpair(j,1)
	   if(neighbour(2,j).eq.loc)call pushpair(j,2)
	   if(neighbour(3,j).eq.loc)call pushpair(j,3)
	end if
	if(k.ne.0)then
	   if(neighbour(1,k).eq.loc)call pushpair(k,1)
	   if(neighbour(2,k).eq.loc)call pushpair(k,2)
	   if(neighbour(3,k).eq.loc)call pushpair(k,3)
	end if

 10	call stackpairempty(k)
c					if stack empty then finish
	if(k.eq.1)go to 100
c					take triangle from stack
	call poppair(tri,pos)

c					test if (x,y) is in circumcircle

         
        dx = x-centres(1,tri)
        dy = y-centres(2,tri)
	dist = dx*dx + dy*dy 

 	if(dist.lt.centres(3,tri))then
           if(nnwrite)write(lud,*)tri
c					Calculate natural neighbour
c					contribution from current triangle.
c
 	   call nn_tri(x,y,points,vertices,centres,tri
     &                 ,data,f,df,v,dvx,dvy)

c					plot triangle and circumcircle

 	   call plot_tc(tri,points,vertices,centres)
	   if(pos.eq.1)then
	      new1 = neighbour(2,tri)
	      new2 = neighbour(3,tri)
	   end if
	   if(pos.eq.2)then
	      new1 = neighbour(1,tri)
	      new2 = neighbour(3,tri)
	   end if
	   if(pos.eq.3)then
	      new1 = neighbour(1,tri)
	      new2 = neighbour(2,tri)
	   end if
	   if(new1.ne.0)then
	      if(neighbour(1,new1).eq.tri)call pushpair(new1,1)
	      if(neighbour(2,new1).eq.tri)call pushpair(new1,2)
	      if(neighbour(3,new1).eq.tri)call pushpair(new1,3)
	   end if
	   if(new2.ne.0)then
	      if(neighbour(1,new2).eq.tri)call pushpair(new2,1)
	      if(neighbour(2,new2).eq.tri)call pushpair(new2,2)
	      if(neighbour(3,new2).eq.tri)call pushpair(new2,3)
	   end if
	end if
	go to 10
 
 100    continue
        call stackpairflush()
	if(v.ne.0.0)then
	   f = f/v
	   df(1) = (df(1) - f*dvx)/v
	   df(2) = (df(2) - f*dvy)/v
c	   if(nnwrite)then
c	   write(*,*)' Overlap voronoi area = ',v
c	   write(*,*)' Normalized value of f   =',f
c	   write(*,*)' Normalized value of dfx =',df(1)
c	   write(*,*)' Normalized value of dfy =',df(2)
c	   end if
c	else
c	   if(nnwrite)then
c	   write(*,*)' Subroutine nn2Di: overlap voronoi area = 0'
c	   write(*,*)' Unnormalized value of f   =',f
c	   write(*,*)' Current value of dfx      =',df(1)
c	   write(*,*)' Current value of dfy      =',df(2)
c	   end if
	end  if

	return
	end
c
c------------------------------------------------------------------------
c
c	nn2Do - calculates natural neighbour interpolation at point x,y
c		when point is outside convex hull 
c
c	Input:
c		x,y			co-ordinates of input point
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
c		hulltriangles(nh)	array of triangles with an edge
c					on the convex hull
c	        nh			number of triangles with an edge
c					on the convex hull
c		data			data values at each node
c
c	Output:
c		f			interpolated data value
c		df(1)			interpolated derivative df/dx
c		df(2)			interpolated derivative df/dy
c		v			area of voronoi cell about (x,y)
c
c	Comments:
c		The following procedure works on the principle that if (x,y) is
c		inside the union of circumcircles then it must be in one 
c		of the circumcircles whose Delaunay triangle has an edge 
c		on the convex hull. It then proceeds to find all circumcircles
c		containing the point by repeating the triangle walking
c		algorithm described by Sloan (1987, and attributed to 
c		D. Watson) starting from each triangle on the hull whose
c		circumcircles contain the point. Note the repeated use of
c		the procedure is necessary because when (x,y) is outside of 
c		the hull the required triangles may not form a contiguous set.
c
c		See comments on plot routine plot_tc above.
c
c		Calls are made to: nn_tri. plot_tc, stackpairinit,
c				    poppair,pushpair & stackpairempty. 
c
c					M. Sambridge, RSES, April 1994.
c
c------------------------------------------------------------------------
c
	Subroutine nn2Do
     &  (x,y,points,vertices,neighbour,centres,
     &   hulltriangles,nh,data,f,df,v)

	real*8		points(2,*)
	real*8		centres(3,*)
	real*8		x,y,xc,yc,dist
	real*8		data(*)
	real*8		f,df(2),v,dvx,dvy
	integer		vertices(3,*)
	integer		neighbour(3,*)
	integer		tri,pos
        integer		hulltriangles(*)
	logical		inside_union
	logical		nnwrite

	common/nnswitches/nnwrite,lud

 	if(nnwrite)write(lud,*)' Result of nn2Do search:'
        inside_union = .false.

c					Initialize variables
	f = 0.d0
	df(1) = 0.d0
	df(2) = 0.d0
	v = 0.d0
	dvx = 0.d0
	dvy = 0.d0
c					Find all circumcircles containing
c					input point

c					test all circumcircles whose
c					Delaunay triangles have an edge
c					on the convex hull
	do 5 list = 1,nh
           tri = hulltriangles(list)
	   xc = centres(1,tri)
	   yc = centres(2,tri)
	   dist = (xc-x)*(xc-x) + (yc-y)*(yc-y)
	   if(dist.lt.centres(3,tri))then
	      if(nnwrite)write(lud,*)' circle test found ',tri
	      call plot_tc(tri,points,vertices,centres)
	      inside_union = .True.
	      loc = tri


	      if(nnwrite)write(lud,*)loc
c					Calculate natural neighbour
c					contribution from current triangle.

 	      call nn_tri(x,y,points,vertices,centres,loc,
     &                    data,f,df,v,dvx,dvy)


c					plot triangle and circumcircle

	      call plot_tc(loc,points,vertices,centres)
	   
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
	         if(neighbour(1,i).eq.loc)call pushpair(i,1)
	         if(neighbour(2,i).eq.loc)call pushpair(i,2)
	         if(neighbour(3,i).eq.loc)call pushpair(i,3)
	      end if
	      if(j.ne.0)then
	         if(neighbour(1,j).eq.loc)call pushpair(j,1)
	         if(neighbour(2,j).eq.loc)call pushpair(j,2)
	         if(neighbour(3,j).eq.loc)call pushpair(j,3)
	      end if
	      if(k.ne.0)then
	         if(neighbour(1,k).eq.loc)call pushpair(k,1)
	         if(neighbour(2,k).eq.loc)call pushpair(k,2)
	         if(neighbour(3,k).eq.loc)call pushpair(k,3)
	      end if
         
 10	      call stackpairempty(k)
c					if stack empty then finish
	      if(k.eq.1)go to 100
c					take triangle from stack
	      call poppair(tri,pos)

c					test if (x,y) is in circumcircle

	      dist = (x-centres(1,tri))*(x-centres(1,tri)) + 
     &               (y-centres(2,tri))*(y-centres(2,tri))
	      if(dist.lt.centres(3,tri))then
                 if(nnwrite)write(lud,*)tri
c					Calculate natural neighbour
c					contribution from current triangle.

 	         call nn_tri(x,y,points,vertices,centres,tri,
     &                       data,f,df,v,dvx,dvy)

c					plot triangle and circumcircle
         
	         call plot_tc(tri,points,vertices,centres)
	         if(pos.eq.1)then
	            new1 = neighbour(2,tri)
	            new2 = neighbour(3,tri)
	         end if
	         if(pos.eq.2)then
	            new1 = neighbour(1,tri)
	            new2 = neighbour(3,tri)
	         end if
	         if(pos.eq.3)then
	            new1 = neighbour(1,tri)
	            new2 = neighbour(2,tri)
	         end if
	         if(new1.ne.0)then
	            if(neighbour(1,new1).eq.tri)call pushpair(new1,1)
	            if(neighbour(2,new1).eq.tri)call pushpair(new1,2)
	            if(neighbour(3,new1).eq.tri)call pushpair(new1,3)
	         end if
	         if(new2.ne.0)then
	            if(neighbour(1,new2).eq.tri)call pushpair(new2,1)
	            if(neighbour(2,new2).eq.tri)call pushpair(new2,2)
	            if(neighbour(3,new2).eq.tri)call pushpair(new2,3)
	         end if
	      end if
	      go to 10
 100          continue
              call stackpairflush()
	   end if
 5      continue
	if(.not.inside_union.and.nnwrite)
     &     write(lud,*)' point is outside all circumcircles'
	if(.not.inside_union.and.nnwrite)
     &     write(*,*)' point is outside all circumcircles'

	if(v.ne.0.0)then
	   f = f/v
	   df(1) = (df(1) - f*dvx)/v
	   df(2) = (df(2) - f*dvy)/v
	end  if

	return
	end
c
c------------------------------------------------------------------------
c
c	nn_tri - calculates natural neighbour contribution from 
c		 triangle tri at point x,y
c
c	Input:
c		x,y			co-ordinates of input point
c		points(2,np)		array of node points	
c	        vertices(3,nt)		array of triangle vertices	
c               centres(3,nt)           centres(j,i) (j=1,2) contains the
c                                       co-ordinates of the centre of 
c                                       circumcircle about Delaunay 
c                                       triangle i, (i=1,...,nt),
c                                       j=3 contains squared radius of circle.
c		tri			index of current triangle
c	Output:
c		f			partial f interpolation
c		df(2)			partial derivative interpolation
c		dvx			partial sum of derivatives
c		dvy			partial sum of derivatives
c		v			partial area of voronoi cell about (x,y)
c
c	Comments:
c		 Note this routine only adds on the contribution to
c		 f and derivatives to sum from last call. All summation
c		 variables are initialized in the calling routine. 
c		 Nodes are assumed to be in counterclockwise direction.
c
c		 The terms dvx, dvy are needed in the calling routine to 
c		 complete the calculation of f and df. 
c
c		 See comments on plot routine plot_c above.
c
c		 Calls are made to: plot_c.
c
c					M. Sambridge, RSES, April 1994.
c
c------------------------------------------------------------------------
c
	Subroutine nn_tri
     &             (x1,x2,points,vertices,centres,tri,
     &              data,f,df,v,dvx,dvy)

	real*8		points(2,*)
	real*8		centres(3,*)
	real*8		x1,x2,f,df(2),v
	real*8		dvx,dvy,a,b,dxw,dyw,det
	real*8		p(2,3),c(2,3),f1(2,3),f2(2,3),cx,cy
	real*8		dcx(2,3),dcy(2,3)
	real*8		data(*)
	integer		vertices(3,*)
	integer		tri
	logical		nnwrite
        integer		c1(3)
        data		c1/2,3,1/

	common/nnswitches/nnwrite,lud



c					set position vectors of nodes
	do 5 i=1,2
	   do 6 j=1,3
	      p(i,j) = points(i,vertices(j,tri))
 6         continue
 5      continue
c					set circumcentre nodes
	cx = centres(1,tri)
	cy = centres(2,tri)
	do 10 i = 1,3
           j = c1(i)
           k = c1(j)

c					calculate centres of three 
c					new triangles

           c(1,i) = (((p(1,k)-p(1,j))*(p(1,k)+p(1,j))/2. 
     &               +(p(2,k)-p(2,j))*(p(2,k)+p(2,j))/2.)
     &               *(x2-p(2,j))
     &             -(((x1-p(1,j))*(x1+p(1,j))/2. 
     &               +(x2-p(2,j))*(x2+p(2,j))/2.)
     &               *(p(2,k)-p(2,j))))/
     &               ((p(1,k)-p(1,j))*(x2-p(2,j))
     &               -(x1-p(1,j))*(p(2,k)-p(2,j)))
 
           c(2,i) = ((p(1,k)-p(1,j))*((x1-p(1,j))
     &              *(x1+p(1,j))/2. 
     &              +(x2-p(2,j))*(x2+p(2,j))/2.)
     &             -((x1-p(1,j))*((p(1,k)-p(1,j))
     &              *(p(1,k)+p(1,j))/2. 
     &              +(p(2,k)-p(2,j))*(p(2,k)+p(2,j))/2.)))/
     &              ((p(1,k)-p(1,j))*(x2-p(2,j))
     &              -(x1-p(1,j))*(p(2,k)-p(2,j)))
 
           xs = c(1,i)
           ys = c(2,i)
           call plot_c(xs,ys,x1,x2)

	   f1(1,i) = c(1,i) - x1
	   f1(2,i) = f1(1,i) 
	   f2(1,i) = c(2,i) - x2
	   f2(2,i) = f2(1,i) 

	   dcx(1,i) = ((p(2,k)-x2)*f1(1,i)-(p(2,j)-x2)*f1(2,i))/
     &               ((p(1,j)-x1)*(p(2,k)-x2)-(p(2,j)-x2)*(p(1,k)-x1))
	   dcx(2,i) = ((p(1,j)-x1)*f1(2,i)-(p(1,k)-x1)*f1(1,i))/
     &               ((p(1,j)-x1)*(p(2,k)-x2)-(p(2,j)-x2)*(p(1,k)-x1))
	   dcy(1,i) = ((p(2,k)-x2)*f2(1,i)-(p(2,j)-x2)*f2(2,i))/
     &               ((p(1,j)-x1)*(p(2,k)-x2)-(p(2,j)-x2)*(p(1,k)-x1))
	   dcy(2,i) = ((p(1,j)-x1)*f2(2,i)-(p(1,k)-x1)*f2(1,i))/
     &               ((p(1,j)-x1)*(p(2,k)-x2)-(p(2,j)-x2)*(p(1,k)-x1))
 10	continue

	do 20 i = 1,3
           j = c1(i)
           k = c1(j)
	   node = vertices(i,tri)
	   det  = ((c(1,j)-cx)*(c(2,k)-cy)
     &            -(c(1,k)-cx)*(c(2,j)-cy))/2.d0
	   v = v + det
	   f = f + det*data(node)
	   a = (dcx(1,j)*(c(2,k)-cy)-(c(1,k)-cx)*dcx(2,j))
	   b = ((c(1,j)-cx)*dcx(2,k)-dcx(1,k)*(c(2,j)-cy))
	   dxw = (a + b)/2.d0
	   a = (dcy(1,j)*(c(2,k)-cy)-(c(1,k)-cx)*dcy(2,j))
	   b = ((c(1,j)-cx)*dcy(2,k)-dcy(1,k)*(c(2,j)-cy))
	   dyw = (a + b)/2.d0
	   df(1) = df(1) + dxw*data(node)
	   df(2) = df(2) + dyw*data(node)
	   dvx =  dvx  + dxw
	   dvy =  dvy  + dyw
 20	continue

c	if(nnwrite)write(*,*)'x',x1,' y',x2,' f ',f,' t ',tri
c	if(nnwrite)write(*,*)' df1',df(1),' df2 ',df(2),' v',v

	return
	end
c
c------------------------------------------------------------------------
c
c	nn2Dr - calculates natural neighbour interpolation at point x,y
c		when point is inside convex hull 
c		(using method based on recursive formula)
c
c	Input:
c		x,y			co-ordinates of input point
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
c
c	Output:
c		f			interpolated data value
c		df(1)			interpolated derivative df/dx
c		df(2)			interpolated derivative df/dy
c		v			area of voronoi cell about (x,y)
c
c	Comments:
c		 On input point (x,y) must be in triangle loc. 
c		 Note: The input triangle contains (x,y) and so it's 
c		 circumcircle must also contain (x,y).
c
c		 In order to use the recursive Lasserre formula to
c		 calculate the nn co-ordinate of the input point with
c		 respect to its ith neighbour we must find the set of
c		 nodes that are both neighbours of i and (x,y), where 
c		 the neighbours of i are determined before addition of 
c		 (x,y). This is done using an approach similar to Watson's
c		 algorithm for updating a Delaunay triangulation, although
c		 we avoid the complication of having to remove `double
c		 entires in the list' and also use Sloan's LIFO stack
c		 approach to avoid an expensive global search over
c		 all triangles.
c
c		 This version does not calculate derivatives.
c
c		 Calls are made to: plot_tc, stackpairinit,
c				    poppair,pushpair & stackpairempty. 
c
c					M. Sambridge, RSES, April 1995.
c
c------------------------------------------------------------------------
c
	Subroutine nn2Dr
     &             (x,y,points,vertices,neighbour,centres,loc,
     &              data,nnpn_max,a,b,p,nodes,pairs,tstore,f,df,v)

	real*8		points(2,*)
	real*8		centres(3,*)
	real*8		data(*)
	real*8		x,y,dist,dx,dy
	real*8		f,df(2),v
	real*8		xo(2)
        real*8  	xs(2)
        real*8          p(2,nnpn_max)
	integer		vertices(3,*)
	integer		neighbour(3,*)
	integer		tri,pos,pcount
	logical		nnwrite
	integer		pairs(2,nnpn_max)
	integer		nodes(nnpn_max)
	integer		tstore(nnpn_max)
        real            a(nnpn_max,2)
        real            b(nnpn_max)

        

	common/nnswitches/nnwrite,lud

 	if(nnwrite)write(lud,*)' Result of nn2Dr search:'

c					Initialize variables
	f = 0.d0
	df(1) = 0.d0
	df(2) = 0.d0
	v = 0.d0

c					Determine list of 
c					neighbours to input point.
	
c
c					Initialize pair list 
c

        pairs(1,1) = vertices(1,loc)
        pairs(2,1) = vertices(2,loc)
        pairs(1,2) = vertices(2,loc)
        pairs(2,2) = vertices(3,loc)
        pairs(1,3) = vertices(3,loc)
        pairs(2,3) = vertices(1,loc)
        pcount = 3

	if(nnwrite)write(lud,*)loc

c					plot input circum-triangle 
c					and circum-circle

 	call plot_tc(loc,points,vertices,centres)

c					record circum-triangle 
        tstore(1) = loc
        tstore(2) = loc
        tstore(3) = loc
c
c					Find all other circumcircles 
c					containing input point and update
c					pair list from these triangles
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
	   if(neighbour(1,i).eq.loc)call pushpair(i,1)
	   if(neighbour(2,i).eq.loc)call pushpair(i,2)
	   if(neighbour(3,i).eq.loc)call pushpair(i,3)
	end if
	if(j.ne.0)then
	   if(neighbour(1,j).eq.loc)call pushpair(j,1)
	   if(neighbour(2,j).eq.loc)call pushpair(j,2)
	   if(neighbour(3,j).eq.loc)call pushpair(j,3)
	end if
	if(k.ne.0)then
	   if(neighbour(1,k).eq.loc)call pushpair(k,1)
	   if(neighbour(2,k).eq.loc)call pushpair(k,2)
	   if(neighbour(3,k).eq.loc)call pushpair(k,3)
	end if

 10	call stackpairempty(k)
c					if stack empty then finish
	if(k.eq.1)go to 100
c					take triangle from stack
	call poppair(tri,pos)

c					test if (x,y) is in circumcircle

        dx = x-centres(1,tri)
        dy = y-centres(2,tri)
	dist = dx*dx + dy*dy 

 	if(dist.lt.centres(3,tri))then
           if(nnwrite)write(lud,*)tri

c					plot circum-triangle and circum-circle

 	   call plot_tc(tri,points,vertices,centres)

c					record each edge of current
c					triangles and triangle
c
	   if(pos.eq.1)then
	      new1 = neighbour(2,tri)
	      new2 = neighbour(3,tri)
              n1 = 3
              n2 = 1
              n3 = 1
              n4 = 2
	   else if(pos.eq.2)then
	      new1 = neighbour(1,tri)
	      new2 = neighbour(3,tri)
              n1 = 2
              n2 = 3
              n3 = 1
              n4 = 2
	   else 
	      new1 = neighbour(1,tri)
	      new2 = neighbour(2,tri)
              n1 = 2
              n2 = 3
              n3 = 3
              n4 = 1
	   end if
           pcount = pcount + 1
           pairs(1,pcount) = vertices(n1,tri)
           pairs(2,pcount) = vertices(n2,tri)
           tstore(pcount) = tri
           pcount = pcount + 1
           pairs(1,pcount) = vertices(n3,tri)
           pairs(2,pcount) = vertices(n4,tri)
           tstore(pcount) = tri

	   if(new1.ne.0)then
	      if(neighbour(1,new1).eq.tri)call pushpair(new1,1)
	      if(neighbour(2,new1).eq.tri)call pushpair(new1,2)
	      if(neighbour(3,new1).eq.tri)call pushpair(new1,3)
	   end if
	   if(new2.ne.0)then
	      if(neighbour(1,new2).eq.tri)call pushpair(new2,1)
	      if(neighbour(2,new2).eq.tri)call pushpair(new2,2)
	      if(neighbour(3,new2).eq.tri)call pushpair(new2,3)
	   end if
	end if
	go to 10
 
 100    continue
        call stackpairflush()

c						check size of pcount
        if(pcount.gt.nnpn_max)then
           write(*,*)' '
           write(*,*)' Error: work arrays in subroutine nn2Dr',
     &               ' not big enough.'
           write(*,*)' Remedy: Increase size of parameter'
           write(*,*)'         nnpn_max in calling program'
           write(*,*)'         current value = ',nnpn_max
           write(*,*)'         required value >= ',pcount
	   stop
        end if
c						write out neighbouring
c						nodes of input points
c
c       do 110 i=1,pcount
c          write(*,*)' neighbouring pair of nodes',pairs(1,i),pairs(2,i)
c110    continue
c
c						find set of common neighbours
c						between input point, and each
c						of its natural neighbours
c						and use this information to
c						set up p matrix of co-ordinates.
c
        j=0
        xs(1) = x
        xs(2) = y
        do 120 i=1,pcount
	   tri = tstore(i)
           do 121 ii=1,2
	      node = pairs(ii,i)
              do 122 kk=1,j
                 if(node.eq.nodes(kk))go to 121
 122          continue
              j = j + 1
	      nodes(j) = node
              p(1,1) = points(1,node)
              p(2,1) = points(2,node)
              jj = 1
              do 130 k=1,pcount
                 if(pairs(1,k).eq.node)then
c                   write(*,*)' Node:',node,' neighbour',pairs(2,k)
                    jj = jj + 1
                    p(1,jj) = points(1,pairs(2,k))
                    p(2,jj) = points(2,pairs(2,k))
                 end if
                 if(pairs(2,k).eq.node)then
c                   write(*,*)' Node:',node,' neighbour',pairs(1,k)
                    jj = jj + 1
                    p(1,jj) = points(1,pairs(1,k))
                    p(2,jj) = points(2,pairs(1,k))
                 end if
 130          continue
              
c						calculate volume of 
c						second-order voronoi cell 
c						for current node using
c						recursive formula and reset
c						origin to circumcentre
c
              xo(1) = centres(1,tri)
              xo(2) = centres(2,tri)

              call second_voronoi (xs,p,2,jj,nnpn_max,2,xo,a,b,vol)
c
              v = v + vol
              f = f + vol*data(node)
 121       continue
 120    continue

        do 140 k=1,j
           nodes(k) = 0
 140    continue

 	if(v.ne.0.0) f = f/v

	return
	end
c
c------------------------------------------------------------------------
c
c	nn2Drd - calculates natural neighbour interpolation at point x,y
c	 	 when point is inside convex hull 
c		 (using method based on recursive formulae)
c
c	Input:
c		x,y			co-ordinates of input point
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
c
c	Output:
c		f			interpolated data value
c		df(1)			interpolated derivative df/dx
c		df(2)			interpolated derivative df/dy
c		v			area of voronoi cell about (x,y)
c
c	Comments:
c		 On input point (x,y) must be in triangle loc. 
c		 Note: The input triangle contains (x,y) and so it's 
c		 circumcircle must also contain (x,y).
c
c		 In order to use the recursive Lasserre formula to
c		 calculate the nn co-ordinate of the input point with
c		 respect to its ith neighbour we must find the set of
c		 nodes that are both neighbours of i and (x,y), where 
c		 the neighbours of i are determined before addition of 
c		 (x,y). This is done using an approach similar to Watson's
c		 algorithm for updating a Delaunay triangulation, although
c		 we avoid the complication of having to remove `double
c		 entires in the list' and also use Sloan's LIFO stack
c		 approach to avoid an expensive global search over
c		 all triangles.
c
c		 This version also calculates derivatives.
c
c		 Calls are made to: plot_tc, stackpairinit,
c				    poppair,pushpair & stackpairempty. 
c
c					M. Sambridge, RSES, April 1995.
c
c------------------------------------------------------------------------
c
	Subroutine nn2Drd
     &             (x,y,points,vertices,neighbour,centres,loc,
     &              data,nnpn_max,a,b,p,nodes,pairs,tstore,f,df,v)

	real*8		points(2,*)
	real*8		centres(3,*)
	real*8		data(*)
	real*8		x,y,dist,dx,dy
	real*8		f,df(2),v,dvx,dvy
	real*8		xo(2)
        real*8  	xs(2)
        real*8          p(2,nnpn_max)
	integer		vertices(3,*)
	integer		neighbour(3,*)
	integer		tri,pos,pcount
	logical		nnwrite
	integer		pairs(2,nnpn_max)
	integer		nodes(nnpn_max)
	integer		tstore(nnpn_max)
        real  		deriv(2)
        real            a(nnpn_max,2)
        real            b(nnpn_max)

        

	common/nnswitches/nnwrite,lud

 	if(nnwrite)write(lud,*)' Result of nn2Dr search:'

c					Initialize variables
	f = 0.d0
	df(1) = 0.d0
	df(2) = 0.d0
	v = 0.d0
	dvx = 0.d0
	dvy = 0.d0
c					Initialize pair list 
c
        pairs(1,1) = vertices(1,loc)
        pairs(2,1) = vertices(2,loc)
        pairs(1,2) = vertices(2,loc)
        pairs(2,2) = vertices(3,loc)
        pairs(1,3) = vertices(3,loc)
        pairs(2,3) = vertices(1,loc)
        pcount = 3

	if(nnwrite)write(lud,*)loc

c					plot input circum-triangle 
c					and circum-circle

 	call plot_tc(loc,points,vertices,centres)

c					record circum-triangle 
        tstore(1) = loc
        tstore(2) = loc
        tstore(3) = loc
c
c					Find all other circumcircles 
c					containing input point and update
c					pair list from these triangles
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
	   if(neighbour(1,i).eq.loc)call pushpair(i,1)
	   if(neighbour(2,i).eq.loc)call pushpair(i,2)
	   if(neighbour(3,i).eq.loc)call pushpair(i,3)
	end if
	if(j.ne.0)then
	   if(neighbour(1,j).eq.loc)call pushpair(j,1)
	   if(neighbour(2,j).eq.loc)call pushpair(j,2)
	   if(neighbour(3,j).eq.loc)call pushpair(j,3)
	end if
	if(k.ne.0)then
	   if(neighbour(1,k).eq.loc)call pushpair(k,1)
	   if(neighbour(2,k).eq.loc)call pushpair(k,2)
	   if(neighbour(3,k).eq.loc)call pushpair(k,3)
	end if

 10	call stackpairempty(k)
c					if stack empty then finish
	if(k.eq.1)go to 100
c					take triangle from stack
	call poppair(tri,pos)

c					test if (x,y) is in circumcircle

        dx = x-centres(1,tri)
        dy = y-centres(2,tri)
	dist = dx*dx + dy*dy 

 	if(dist.lt.centres(3,tri))then
           if(nnwrite)write(lud,*)tri

c					plot input circum-triangle 
c					and circum-circle

 	   call plot_tc(tri,points,vertices,centres)

c					record each edge 
c					of current triangle.
c
	   if(pos.eq.1)then
	      new1 = neighbour(2,tri)
	      new2 = neighbour(3,tri)
              n1 = 3
              n2 = 1
              n3 = 1
              n4 = 2
	   else if(pos.eq.2)then
	      new1 = neighbour(1,tri)
	      new2 = neighbour(3,tri)
              n1 = 2
              n2 = 3
              n3 = 1
              n4 = 2
	   else 
	      new1 = neighbour(1,tri)
	      new2 = neighbour(2,tri)
              n1 = 2
              n2 = 3
              n3 = 3
              n4 = 1
	   end if
           pcount = pcount + 1
           pairs(1,pcount) = vertices(n1,tri)
           pairs(2,pcount) = vertices(n2,tri)
           tstore(pcount) = tri
           pcount = pcount + 1
           pairs(1,pcount) = vertices(n3,tri)
           pairs(2,pcount) = vertices(n4,tri)
           tstore(pcount) = tri

	   if(new1.ne.0)then
	      if(neighbour(1,new1).eq.tri)call pushpair(new1,1)
	      if(neighbour(2,new1).eq.tri)call pushpair(new1,2)
	      if(neighbour(3,new1).eq.tri)call pushpair(new1,3)
	   end if
	   if(new2.ne.0)then
	      if(neighbour(1,new2).eq.tri)call pushpair(new2,1)
	      if(neighbour(2,new2).eq.tri)call pushpair(new2,2)
	      if(neighbour(3,new2).eq.tri)call pushpair(new2,3)
	   end if
	end if
	go to 10
 
 100    continue
        call stackpairflush()
c						check size of pcount
        if(pcount.gt.nnpn_max)then
           write(*,*)' '
           write(*,*)' Error: work arrays in subroutine nn2Drd',
     &               ' not big enough.'
           write(*,*)' Remedy: Increase size of parameter'
           write(*,*)'         nnpn_max in calling program'
           write(*,*)'         current value = ',nnpn_max
           write(*,*)'         required value >= ',pcount
	   stop
        end if

c						write out neighbouring
c						nodes of input points
c
c       do 110 i=1,pcount
c          write(*,*)' neighbouring pair of nodes',pairs(1,i),pairs(2,i)
c110    continue
c
c						find set of common neighbours
c						between input point, and each
c						of its natural neighbours.
c						and use this information to
c						set up p matrix of co-ordinates.
c
        j=0
        xs(1) = x
        xs(2) = y
        do 120 i=1,pcount
           tri = tstore(i)
           do 121 ii=1,2
	      node = pairs(ii,i)
              do 122 kk=1,j
                 if(node.eq.nodes(kk))go to 121
 122          continue
              j = j + 1
	      nodes(j) = node
              p(1,1) = points(1,node)
              p(2,1) = points(2,node)
              jj = 1
              do 130 k=1,pcount
                 if(pairs(1,k).eq.node)then
c                   write(*,*)' Node:',node,' neighbour',pairs(2,k)
                    jj = jj + 1
                    p(1,jj) = points(1,pairs(2,k))
                    p(2,jj) = points(2,pairs(2,k))
                 end if
                 if(pairs(2,k).eq.node)then
c                   write(*,*)' Node:',node,' neighbour',pairs(1,k)
                    jj = jj + 1
                    p(1,jj) = points(1,pairs(1,k))
                    p(2,jj) = points(2,pairs(1,k))
                 end if
 130          continue
              
c						calculate volume and derivatives
c						of second-order voronoi cell 
c						for current node using
c						recursive formula and reset
c                                               origin to circumcentre
c
              xo(1) = centres(1,tri)
              xo(2) = centres(2,tri)

              call second_voronoi_d 
     &             (xs,p,2,jj,nnpn_max,2,xo,a,b,vol,deriv)
c
              v = v + vol
              f = f + vol*data(node)
              df(1) = df(1) + data(node)*deriv(1)
              df(2) = df(2) + data(node)*deriv(2)
              dvx = dvx + deriv(1)
              dvy = dvy + deriv(2)
 121       continue
 120    continue

        do 140 k=1,j
           nodes(k) = 0
 140    continue

 	if(v.ne.0.0)then
 	   f = f/v
 	   df(1) = (df(1) - f*dvx)/v
 	   df(2) = (df(2) - f*dvy)/v
 	end  if

	return
	end
c
c------------------------------------------------------------------------
c
c	nn2Df - calculates natural neighbour interpolation at point x
c		when point is inside convex hull 
c		(using method based on a closed 2-D formula)
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
c		ltri			logical work array size of nt_max
c		lnode			logical work array size of nt_max
c
c	Output:
c		f			interpolated data value
c		v			area of voronoi cell about x
c
c	Comments:
c		 On input point x must be in triangle loc. 
c		 Note: The input triangle contains x and so it's 
c		 circumcircle must also contain x.
c
c		 In order to use the closed formula to
c		 calculate the nn co-ordinate of the input point with
c		 respect to its ith neighbour we must find the set of
c		 nodes that are both neighbours of i and x, where 
c		 the neighbours of i are determined before addition of 
c		 x. This is done using a LIFO stack to avoid an 
c		 expensive global search over all triangles.
c
c		 It is assumed that the logical arrays ltri(nt_max) 
c		 and lnode are initialized to false on input.
c
c		 This version does not calculate derivatives.
c
c		 Calls are made to: plot_tc, stackpairinit,
c				    poppair,pushpair & stackpairempty
c				    indexx & second_v_area.
c
c					M. Sambridge, RSES, April 1996.
c
c------------------------------------------------------------------------
c
	Subroutine nn2Df
     &  (x,points,vertices,neighbour,centres,loc,data,nnpn_max,
     &   pverts,vverts,tstore,vpairs,angles,index,ltri,lnode,f,v)

	real*8		points(2,*)
	real*8		centres(3,*)
	real*8		data(*)
	real*8		x(2),dist,dx,dy
	real*8		f,v
        real*8          pverts(2,nnpn_max)
        real*8          vverts(2,nnpn_max)
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
c
	integer		cyc(3)
        data		cyc/2,3,1/

        

	common/nnswitches/nnwrite,lud

c					Initialize variables
	f = 0.d0
	v = 0.d0
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
           write(*,*)' Error: work arrays in subroutine nn2Df',
     &               ' not big enough.'
           write(*,*)' Too many circum-triangles for this node'
           write(*,*)' Remedy: Increase size of parameter'
           write(*,*)'         nnpn_max in calling program'
           write(*,*)'         current value = ',nnpn_max
           write(*,*)'         required value >= ',tcount+2
	   stop
        end if
c
c       call xsave
c
c
c						find list of vertices
c						of new voronoi cell about
c						x and store the pair of
c						nodes associated with each
c						vertex
c       write(*,*)' Circum-triangles of x:'
c
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
                 call circum (x,points(1,na),
     &                points(1,nb),vverts(1,nv))

                 vpairs(1,nv) = na
                 vpairs(2,nv) = nb
                 angles(nv) = pangle(x,vverts(1,nv))
                  
              end if
           end do
           
 110    continue
        
c						sort vertices 
c
        call indexx(nv,angles,index)
c
c						debug I/O
c						
c
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
              nverts=0
c             write(*,*)' Neighbour node =',node
c
c						find pair of `external' vertices
c						of second-order voronoi cell  
c
              do i=1,nv
                 k = index(i)
                 if(vpairs(1,k).eq.node.or.
     &              vpairs(2,k).eq.node)then
                    nverts = nverts + 1
                    pverts(1,nverts) = vverts(1,k) - x(1)
                    pverts(2,nverts) = vverts(2,k) - x(2)
                    ip = i+1
                    im = i-1
                    if(i.eq.nv)ip=1
                    if(i.eq.1)im=nv
                    kp = index(ip)
                    km = index(im)
                    if(vpairs(1,kp).eq.node.or.
     &                 vpairs(2,kp).eq.node)then
                       nverts = nverts + 1
                       pverts(1,nverts) = vverts(1,kp) - x(1)
                       pverts(2,nverts) = vverts(2,kp) - x(2)
                    else if(vpairs(1,km).eq.node.or.
     &                      vpairs(2,km).eq.node)then
                       nverts = nverts + 1
                       pverts(1,1) = vverts(1,km) - x(1)
                       pverts(2,1) = vverts(2,km) - x(2)
                       pverts(1,2) = vverts(1,k) - x(1)
                       pverts(2,2) = vverts(2,k) - x(2)
                    else
                       write(*,*)' Error in nn2Df'
                       write(*,*)
     &                 ' Can not find neighbouring vertices' 
                       write(*,*)' node=',node,' i =',i
                    end if
                    go to 150
                 end if
              end do
 150          continue
c
c						find `internal' vertices of 
c						second-order voronoi cell  
	      do 140 jt = 1,tcount
                 jtri = tstore(jt)
                 do jn = 1,3
                    jnode = vertices(jn,jtri)
                    if(jnode.eq.node)then
                       nverts = nverts + 1
c						record circum-centre
c
                       pverts(1,nverts) = centres(1,jtri) - x(1)
                       pverts(2,nverts) = centres(2,jtri) - x(2)
                       go to 140
                    end if
                 end do
 140          continue
c						Calculate volume of 
c						second-order voronoi cell 
c						for current node using
c						direct 2-D formula 
c
              call second_v_area
     &             (x,nverts,pverts,angles,area)
c
           v = v + area
           f = f + area*data(node)

 130       continue
 120    continue
 	if(v.ne.0.0) f = f/v
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
