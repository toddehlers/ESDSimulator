c initialize_general_parameters.f

      subroutine initialize_general_parameters
     & (surfscale,dt,iadjust,endtime,ishow,writetime,nshortwrite,iflux,
     &  run_name,nrun_name,iadapt,ihorizontal,iflexure,hflex,
     &  ixflex,iyflex,thickflex,ym,pratio,rhocflex,rhoaflex,oro_length,
     &  oro_height,oro_scale,wind_direction,xlf_AL,sea_level,
     &  ideposition,idiffusion,rain_vel,width_c,thresh,calc_rain,
     &  iflag_oro,iflag_uni,pmax,cohes,rho,grav,xk0,xk1,dtc,distmax,
     &  ilandslide,lsmeth,nx,ny,nnode,sidex,sidey,uplift_rate,advec_vel,
     &  xkf,xlf_BR,xkdiff,tecflag)

c subroutine to initialize general parameters

c INPUT: surfscale       = scaling surface (in km**2)

c OUTPUT: iadjust        = allows for dynamic time stepping (=0
c                          means that time step is fixed; =1 means
c                          that time step is adjusted dynamically)
c         dt             = (initial) time step length (in yr)
c         endtime        = final time (in yr)
c         ishow          = frequency of graphic displays (in time steps)
c         writetime      = frequency of output (in yr)
c         nshortwrite    = frequency of short screen output (in time steps)
c         iflux          = frequency of flux saves (=0 no flux save;
c                          otherwise frequency in time steps); this output
c                          is of total flux through the fixed nodes
c         run_name       = name of this run; must also be the name
c                          of an existing folder where all output files
c                          will be stored
c         iadapt         = flag to allow dynamic remeshing (=0 no; =1 yes)
c         surfmin        = minimum surface area that can be reached by 
c                          remeshing/refining
c                          **(removed see initialize_nodal_geometry 6/14/1 DS)
c         ihorizontal    = flag to permit horizontal mesh movement (=0
c                          means no horizontal movement; =1 means horizontal
c                          movements permitted
c         iflexure       = flag to permit flexural isostasy (=0 no flexure;
c                          =1 flexure)
c         hflex          = size (in km) of the square mesh on which the
c                          thin elastic plate calculations are done
c         ixflex         = flag to permit flexure in the x-direction
c                          (=0 no elastic strength in x-direction;
c                          =1 means elastic strength in x-direction)
c         iyflex         = flag to permit flexure in the y-direction
c                          (=0 no elastic strength in y-direction;
c                          =1 means elastic strength in y-direction)
c         thickflex      = elastic thickness in km
c         ym             = Young Modulus (in Pa)
c         pratio         = Poissons's ratio
c         rhocflex       = crustal density (in kg/m**3)
c         rhoaflex       =asthenospheric density (in kg/m**3)
c         oro_length     = orographic length scale (in km)
c                          (=0 means no orographic control on precipitation)
c         oro_height     = orographic height scale (in m)
c         oro_scale      = background precipitation (in adequate units)
c         wind_direction = wind direction (0= along x-axis)
c         xlf_AL         = fluvial erosion length scale for alluvials
c                          (in km)
c         sea_level      = where sea-level is (in km above the 0 datum)
c                          below sea-level no sediment transport is
c                          allowed
c         ideposition    = flag to allow for fluvial deposition (=0 erosion
c                          only; =1 sedimentation allowed)
c         idiffusion     = flag to allow diffusion processes (=0 no
c                          diffusion; =1 diffusion allowed)i
c
c         FLUVIAL INCISION/TRANSPORT MODIFICATIONS (S. Willett ~ 6/00)
c         xkf            = fluvial erosion constant []
c         xlf_BR         = bedrock erosion length scale [m]
c         width_c        = coef controling channel width as fnct discharge
c                          [sqrt(yr/m)]
c         thresh         = discharge threshold for channel formation
c                          [m km2/yr]
c         rain_vel       = precipitation rate for a uniform rainfall model
c                          [m/yr]
c         xlf_AL         = fluvial erosion length scale for alluvial material
c                          [m]
c
c         LANDSLIDING PARAMS (via P. van der Beek and B. Champel)
c         pmax           = landslide critical slope angle of material
c         cohes          = hill slope effective cohesion
c         rho            = density of hiill slope material
c         grav           = accel. due to gravity
c         xk0            = landslide scaling parameter (?)
c         xk1            = landslide scaling parameter (?)
c         distmax        = maximum landslide size
c                          [m] (changed by DS 6/15/1)
c added by DS 11/13/1
c         lsmeth         = 1: probabilistic method from pvdb
c                          2: simple slope threshold
c                             pmax is the critical slope
c                             distmax is the maximum distance
c                          3: simpler slope threshold in landslide_simple.f
c                             pmax is the critical slope DS 11/17/1
c added by DS 7/4/1
c         sidex,sidey    = x,y lengths of domain [km]
c         nx,ny          = x,y number of nodes
c         nnode          = number of nodes

c subroutines called:
c - read_but_skip_comment
c -iread_but_skip_comment

      common /vocal/ ivocal

      character*256    run_name
      real Ne,D

      iread=0

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc Read in model input below if iread=0                                       cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if (iread.eq.0) then
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc Define run name                                                          ccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        do k=1,256
          run_name(k:k)=' '
        enddo
        run_name(1:14)="output/Cascade"
        do k=256,1,-1
          if (run_name(k:k).eq.' ') nrun_name=k-1
        enddo
        if (nrun_name.eq.0) then
          print*,'No run name available '
          stop
        endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc Define type of output wanted
ccc 0=All output files written
ccc 1=Only the tecplot formatted output files and the geometry and topography
ccc   files are written
      tecflag=0
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc Define model time step parameters                                        ccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Enable dynamic time step calculation
      iadjust=1
c Initial time step length [y]
      dt=50.0
c Total model run time [y]
      endtime=10.001e6
c Frequency of graphical display [time steps]
      ishow=25
c Frequency of writing model output [y]
      writetime=5.0e5
c Frequency of screen output [time steps]
      nshortwrite=1000
c Frequency of flux outputs [time steps]
      iflux=0
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc Define model mesh size, resolution and behavior                          ccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Number of nodes along x and y edges of model
      nx = 401
      ny = 155
      nnode = nx*ny
c Length of x and y sides of model [km]
      sidex = 400.
      sidey = 150.
c Enable adaptive remeshing
      iadapt=0
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc Define model tectonics                                                   ccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Uplift rate [mm/y]
      uplift_rate = 0.5
c Allow horizontal advection
      ihorizontal = 0
c Maximum advection velocity [mm/y ????]
      advec_vel = 2.5
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc Define Flexure Parameters  (added by TE 06/08, missing for some reason)  ccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      iflexure=0
      hflex=1000.
      ixflex=1
      iyflex=1
      thickflex=15.
      ym=1.e11
      pratio=0.25
      rhocflex=2750.
      rhoaflex=3300.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc Define erosional properties                                              ccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

ccc Fluvial erosion parameters ccccccccccccccccccccccccccccccccccccccccccccccccc

c Fluvial erosion constant []
      xkf = 2.e-4
c Bedrock erosion length scale [m]
      xlf_BR = 1000.
c Coefficient controling channel width as function of discharge [sqrt(yr/m)]
      width_c= 0.1
c Discharge threshold for channel formation [m km2/yr]
      thresh = 0.0
c Fluvial erosion length scale for alluvial material [m]
c NOTE (TAE) xlf_AL set to 100.0 in our old version
      xlf_AL=100.
c Enable fluvial depostion
      ideposition=0
c Define sea level elevation; no erosion below sea level [km]
      sea_level=0.

ccc Hillslope diffusion parameters ccccccccccccccccccccccccccccccccccccccccccccc
c Enable hillslope diffusion
      idiffusion=1
c Diffusion constant [km2/yr]
      xkdiff = 2.e-5

ccc Landslide erosion parameters ccccccccccccccccccccccccccccccccccccccccccccccc
c landsliding
c Enable landsliding
      ilandslide=0
c Landslide method to be used (see above)
      lsmeth = 3
c Threshold hillslope angle for landsliding (used for lsmeth=1-3)
      pmax=25.
      pmax=pmax*3.14159/180.
c Maximum transport distance (used for lsmeth=1,2) [m]
      distmax=1000.
c Parameters for lsmeth=1 below
c Hillslope effective cohesion
      cohes=6.e4
c Density of hillslope material [kg/m3]
      rho=2500.
c Acceleration due to gravity [m/s2]
      grav=9.8
c suggested by van der Beek: xk0=0.01, xk1=1., dtc = 100. yrs
c Landslide scaling parameter (?)
      xk0 =0.00
c Landslide scaling parameter (?)
      xk1 =1.
      dtc = 1.e2
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc Define precipitation model                                               ccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c Enable orographic precipitation
      iflag_oro = 1
c Enable uniform precipitation
      iflag_uni = 0
c Time between updates of precipitation field [yr]
      calc_rain=10000.

ccc Other precip model parameters are defined in rainmaker.f

ccc Nothing below here actually does anything to modify the precipitation field
ccc This section should probably be removed in future versions
ccc dwhipp - 01/08

ccc Orographic model ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Orographic length scale (zero means no orographic control on precip [km]
c NOTE: This does nothing.  Modify variables at the top of rainmaker.f to change
c   the effects of orographic length
      oro_length=1.
c Orographic height scale [m?]
c NOTE: This does nothing.  Modify variables at the top of rainmaker.f to change
c   the effects of orographic height
      oro_height=1.
c Background precipitation
c NOTE: This does nothing.  Modify variables at the top of rainmaker.f to change
c   the orographic scaling
      oro_scale=50.
c Wind direction (0 = along x-axis)
c NOTE: This does nothing.  Modify variables at the top of rainmaker.f to change
c   the prevailing wind direction
      wind_direction=180.

ccc Uniform preciptation model ccccccccccccccccccccccccccccccccccccccccccccccccc
c Precipitation rate for a uniform rainfall model [m/yr]
c NOTE: This does nothing.  Modify variables at the top of rainmaker.f to change
c   the precipitation rate
      rain_vel=4.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc NOT USED - IGNORE                                                        ccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c now in initialize_nodal_parameters
c      surfmin=surfscale/4.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc Conversions for units above                                              ccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c convert advec_vel from [mm/y] to [km/y]
      advec_vel=advec_vel*1e-6
c convert uplift_rate from [mm/y] to [m/y]
      uplift_rate=uplift_rate*1e-3


      else
        print*,'Enter run name >'
        read (*,'(a)') run_name
        do k=256,1,-1
          if (run_name(k:k).eq.' ') nrun_name=k-1
        enddo
        if (nrun_name.eq.0) then
          print*,'No run name available '
          stop
        endif
        open (7,file=run_name(1:nrun_name)//
     &       '/cascade.parameters.in',status='old',err=997)
        goto 998
997     print*,'Cant open input file'
        print*,'You must create a folder/directory named: '//
     &       run_name(1:nrun_name)
        print*,'where the input file cascade.parameters.in resides'
        stop
998     continue
c NOTE: Section below needs to be updated
c dwhipp - 01/08
        call read_but_skip_comment (7,1,dt)
        call iread_but_skip_comment (7,1,iadjust)
        call read_but_skip_comment (7,1,endtime)
        call iread_but_skip_comment (7,1,ishow)
        call read_but_skip_comment (7,1,writetime)
        call iread_but_skip_comment (7,1,nshortwrite)
        call iread_but_skip_comment (7,1,iflux)
        call iread_but_skip_comment (7,1,iadapt)
        call read_but_skip_comment (7,1,surfmin)
        call iread_but_skip_comment (7,1,ihorizontal)
        call iread_but_skip_comment (7,1,iflexure)
        call read_but_skip_comment (7,1,hflex)
        call iread_but_skip_comment (7,1,ixflex)
        call iread_but_skip_comment (7,1,iyflex)
        call read_but_skip_comment (7,1,thickflex)
        call read_but_skip_comment (7,1,ym)
        call read_but_skip_comment (7,1,pratio)
        call read_but_skip_comment (7,1,rhocflex)
        call read_but_skip_comment (7,1,rhoaflex)
        call read_but_skip_comment (7,1,oro_length)
        call read_but_skip_comment (7,1,oro_height)
        call read_but_skip_comment (7,1,oro_scale)
        call read_but_skip_comment (7,1,wind_direction)
        call read_but_skip_comment (7,1,xlf_AL)
        call read_but_skip_comment (7,1,sea_level)
        call iread_but_skip_comment (7,1,ideposition)
        call iread_but_skip_comment (7,1,idiffusion)
        close (7)
      endif

c if iecho = 1 the input parameters are echoed on the screen

      iecho=1

      print*,'nx = ',nx
      print*,'ny = ',ny
      print*,'nnode = ',nnode
      print*,'sidex = ',sidex
      print*,'sidey = ',sidey
      print*,'uplift_rate = ',uplift_rate
      print*,'xkf = ',xkf
      print*,'xlf_BR = ',xlf_BR
      print*,'xkdiff = ',xkdiff
      print*,'thresh = ',thresh
      print*,'width_c = ',width_c
      print*,'dt = ',dt
      print*,'iadjust = ',iadjust
      print*,'endtime = ',endtime
      print*,'ishow = ',ishow
      print*,'writetime = ',writetime
      print*,'nshortwrite = ',nshortwrite
      print*,'iflux = ',iflux
      print*,'run_name(1:6) = ',run_name(1:6)
      print*,'iadapt = ',iadapt
c      print*,'surfmin = ',surfmin
      print*,'ihorizontal = ',ihorizontal
      print*,'iflexure = ',iflexure
      if(iflexure.ne.0)then
        print*,'hflex = ',hflex
        print*,'ixflex = ',ixflex
        print*,'iyflex = ',iyflex
        print*,'thickflex = ',thickflex
        print*,'ym = ',ym
        print*,'pratio = ',pratio
        print*,'rhocflex = ',rhocflex
        print*,'rhoaflex = ',rhoaflex
      endif
      print*,'oro_length = ',oro_length
      print*,'oro_height = ',oro_height
      print*,'oro_scale = ',oro_scale
      print*,'wind_direction = ',wind_direction
      print*,'xlf_AL = ',xlf_AL
      print*,'sea_level = ',sea_level
      print*,'ideposition = ',ideposition
      print*,'idiffusion = ',idiffusion
      print*,'ilandslide = ',ilandslide
      if (ilandslide.ne.0) then
        print*,'lsmeth = ',lsmeth
        print*,'pmax = ',pmax
        print*,'cohes = ',cohes
        print*,'rho = ',rho
        print*,'k0 = ',xk0
        print*,'k1 = ',xk1
        print*,'dtc = ',dtc
        print*,'distmax = ',distmax
      endif
c print out Ne and D
      Ne = (xkf/xlf_BR)*(1000.*sidex)/(width_c*sqrt(uplift_rate/1000.))
      D  = xkdiff/(sidex*(uplift_rate/1000.))
      print *,'Ne = ',Ne, ', log10(Ne) = ',alog10(Ne)
      print *,'D  = ',D,  ', log10(D)  = ',alog10(D)

      return
      end
