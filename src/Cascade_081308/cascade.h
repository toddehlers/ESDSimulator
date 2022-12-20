c cascade.h
c
c nx and ny are the linear dimensions of the initial rectangular grid
c note that there is no reason for the grid to be rectangular
c this is done here because it is still the most popular way of
c designing numerical meshes

c nnodemax is the maximum number of nodes that the grid is allowed to grow
c into. If you decide not to use the idynamic=1 option you should chose
c nnodemax as close as possible to nnode

c nbmax is the maximum number of "natural neighbours" that any point can have
c for a general set of points nbmax=20 seems appropriate. However it is
c possible that for arbitrary sets of points, nbmax could be > 20. If it is
c the case, a warning message will be produced by cascade

c ntmax is the maximum nuber of triangles the circumcircles of which contain
c a common point

c Note that it will be assumed that the maximum number of triangles is
c three times the maximum number of nodes (this is very safe as the number
c of triangles is usually of the order of twice the number of nodes).

c nparam is the number of geomorphic parameters to be carried
c by the nodes

c nmemory is the number of working arrays carried by the nodes

c nflex is the discretization used to solve the flexural problem
c it must be a power of 2

c      parameter (nnodemax=205*79,nparam=3,nmemory=8,
c     &           nbmax=50,ntmax=30,nflex=256)

c nodemax boosted up to run higher resolution models on petrarch cluster
c dwhipp 04/07

      parameter (nnodemax=1001*1001,nparam=3,nmemory=8,
     &           nbmax=50,ntmax=30,nflex=256)

c x and y are the x- and y-coordinates of the nodes in km
c h is the current topography
c hi is the initial topography
c h0 is the location of the bedrock interface
c all h's are in m

      real       x(nnodemax),y(nnodemax),hi(nnodemax)
      real       h(nnodemax),h0(nnodemax)
c
c param are geomorphic parameters attached to each nodes
c there are nparam of them
c
      real       param(nnodemax,nparam)
c
c param(*,1)=fluvial erosion constant
c param(*,2)=bedrock erosion length scale
c param(*,3)=diffusion erosion constant
c
c memory are variables that have to be stored from one step to the next
c for each node
c there are nmemory of them
c
      real       memory(nnodemax,nmemory)
c
c memory(*,1)=dhcrit
c memory(*,2)=dhfluvial
c memory(*,3)=dhdiff
c memory(*,4)=hiso
c memory(*,5)=fix
c memory(*,6)=newsurface
c memory(*,7)=surface
c memory(*,8)=dhlandslide (added by Ehlers 6/01)
c
c work is a working array
c
      real       work(nnodemax)

c water is the amount of water that drains down the landscape
c it is equivalent to the discharge
c sediment is the sediment load in the rivers

      real       water(nnodemax),sediment(nnodemax)

c slope is the slope between a node and its donor neighbour
c note that slopes are in meter per kilometer as our horizontal
c length unit is a kilometer while the horizontal unit is the meter

      real       slope(nnodemax),length(nnodemax)

c ndon is the name of the donor neighbour node
c nn is the list of neighbours
c  nb is the number of neighbours for each node

      integer    ndon(nnodemax),nn(nbmax,nnodemax)
      integer    nb(nnodemax)
      integer 	 nb2(nnodemax)
      integer 	 nn2(nbmax,nnodemax)

c ibucket is a working array that is used in the
c "pass the bucket" algorithm on which the cascade method is based
c to define the river network (the ndon array)

      integer    ibucket(nnodemax)

c the following arrays are also used in the cascade algorithm

      integer    iorder(nnodemax)
      integer    itype_node(nnodemax)

c nwork is a working array used in determing the catchment to which each
c node belongs; that is stored in the ncat array which has the name
c of the exiting node of the catchment (that the way catchments are named)

      integer    nwork(nnodemax),ncat(nnodemax)

c nsill and nempty are used in the algorithm that looks for
c sill nodes in case of local minima
c nlake is a flag that determined whether a node belong to a
c lake or not

      integer    nsill(nnodemax),nempty(nnodemax),nlake(nnodemax)

c influx is the flux of material into the landscape
c brought in the system by the tectonic uplift
c outflux is the flux out of the system, ie through the nodes
c where fix=0.

      real       influx,outflux

c rain_vel is the precipitation rate when considering a uniform rainfall
c (non orographic).

       real      rain_vel

c
c initialize precip. array - added by DW 10/06
c TAE 7/01
c NOTE: if you unccomment the following line the program will execute, but
c if start to pass prec to other subroutines you will seg fault.  I'm not
c sure why at this point? TAE
      real  prec(nnodemax)
      real  x_gr(nnodemax)
      real  y_gr(nnodemax)
      real  prec_gr(nnodemax,nnodemax)



c the following arrays are needed in the natural neighbour routines
c have a look inside the library routienes to figure out what
c their purpose is

      real*8     points(2,nnodemax)
      real*8     eps
      integer    vertices(3,nnodemax*3)
      integer    neighbour(3,nnodemax*3)
      integer    nodes(nnodemax)
      integer    vis_tlist(nnodemax)
      integer    vis_elist(nnodemax)
      integer    add_tlist(nnodemax)
      integer    nodelist(nbmax)
      integer    tlist(nbmax)
      logical    c_list(2*nbmax)
      integer    v_local(3,2*nbmax)
      integer    n_local(3,2*nbmax)
      logical    mask(nnodemax*3),mask_e(3,nnodemax*3)
      logical    inactive(nnodemax)

c the following arrays are used to calculate the surface of
c voronoi cells

      real       xy(2),pp(2,nbmax),aa(nbmax,2),bb(nbmax)

c the following arrays are used to solve the diffusion equation iteratively

      integer    kcon(ntmax,nnodemax),jcon(ntmax,nnodemax)
      real       hp(nnodemax)
      integer    nkcon(nnodemax)
      real       ael1(6,nnodemax*3),ael2(6,nnodemax*3)
      real       bel(nnodemax),diag(nnodemax)

c itadd and jtadd are used when dynamic remeshing is turned
c on they are used to determine where resolution has to be increased

      integer    itadd(nnodemax),jtadd(nnodemax*3)

c the following arrays are used in the flexural isostasy
c calculations
c nflex is the resolution at which the FFT are done to calculate the
c flexural response; nflex has to be a power of 2

      real       flex(nflex,nflex),work_flex(nflex,nflex)

      character  run_name*256
c Variables needed for landsliding routine landslide.f
c The following lines were added by Ehlers 6/01
	  real 		 smax(nnodemax),tt(nnodemax)
	  integer	 cell(nnodemax,nbmax,2)


c landslide time series variables DS 6/15/1
c  these store info between write_output calls
c  might be problems if there are more than nnodemax slides
c  during a period
c      integer    snum
c      real       stime(nnodemax),sx(nnodemax),sy(nnodemax)
c      real       svol(nnodemax),sarea(nnodemax)

c dslope:  downstream slope found after all erosion DS 11/18/1
       real       dslope(nnodemax)

       integer    bdry(nnodemax)

c timeint: needed for naming tecplot output files
       integer    timeint
