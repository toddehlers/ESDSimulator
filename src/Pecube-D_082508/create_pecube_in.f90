! Subroutine that reads in input parameters from Pecube.in, modifies variables,
! and writes variables to scratch file for Pecube.f90 to read in

      subroutine create_pecube_in

! The user can modify this file at will but the output it produces (Pecube.in)
! must obey rules that are described in the user guide and in Pecube.f90

      real,dimension(:,:,:),allocatable::zNZ
      real,dimension(:,:),allocatable::z
      real,dimension(:),allocatable::timek,Pe,topomag,topooffset,x1f
      real,dimension(:),allocatable::y1f,x2f,y2f,def,dif,zoff,geoflag
      real,dimension(:),allocatable::Pe2
      real,dimension(:),allocatable::theta,phi,mft_ratein,mbt_ratein,mct_ratein,stf_ratein,thermflag
      integer,dimension(:),allocatable::nfnme,age_flags,nobsfile
      real,dimension(:),allocatable::zmin,zmax
      real*8,dimension(:),allocatable :: age_obs,dage_obs
      character,dimension(:),allocatable::fnme*300,obsfile*300
      real dummy1,dummy2,dummy3,dummy4,dummy5,dummy6,dummy7
      real*8 isdiff,shdiff,lhdiff,ghdiff,thdiff,ishp,shhp,lhhp,ghhp,thhp
      real*8 istc,shtc,lhtc,ghtc,thtc,isrho,shrho,lhrho,ghrho,thrho
      real*8 ishc,shhc,lhhc,ghhc,thhc
      character run*100,line*1024,temp_filename*300,topo_prefix*50,file_iter*4,det_calc*300
      integer num_topo_files,topo_file_flag,prefix_length,He_flag
      real xlonobs,xlatobs,heightobs
      integer shearint,nrun,ios,node_thresh
      real*8 fric
      character cascadedir*100

      open (54,file='input/Pecube-D/Pecube.in',status='unknown',iostat=ios)
      if (ios.ne.0) then
        open (54,file='Pecube.in',status='unknown')
      endif
      open (55,status='scratch')
    1 read (54,'(a1024)',end=2) line
      if (line(1:1).ne.'$'.and. line(1:1).ne.' ') write (55,'(a)') line
      goto 1
    2 close (54)
      rewind (55)

! run is the name of the run (assumes a directory of that name exists)
! should be 5 character long
      do i=1,100
        run(i:i)=' '
      enddo
      read (55,'(a100)') run
      do i=1,100
        if (run(i:i).ne.' ') nrun=i
      enddo
!         run='output/Pecube-D'

! Number of topography files to be loaded (added 2/7/2007 by cspath)
    read(55,*) num_topo_files
    
! Reads the topography flag value (1=User lists all files individually,
! 2=User sets prefix and Pecube iterates through all files with that prefix)
    read(55,*) topo_file_flag

! Checks that the program should be reading in topography file names
    if (num_topo_files .ne. 0) then
          allocate(fnme(num_topo_files),nfnme(num_topo_files))
          if (topo_file_flag .eq. 1) then
! Reads in the names of the topography files (added 2/7/2007 by cspath)
            do k=1,num_topo_files
              do i=1,300
                fnme(k)(i:i)=' '
              enddo
              read (55,'(a)') fnme(k)
              do i=1,300
                if (fnme(k)(i:i).ne.' ') nfnme(k)=i
              enddo
            enddo

          else
! Reads in topography prefix name
            read(55,*) topo_prefix
            do i=1,50
              if (topo_prefix(i:i).ne.' ') prefix_length=i
            enddo
          endif

        else
          allocate(fnme(1),nfnme(1))
          do i=1,300
            fnme(1)(i:i)=' '
          enddo
          read (55,'(a)') fnme(1)
          do i=1,300
            if (fnme(1)(i:i).ne.' ') nfnme(1)=i
          enddo
        endif

! nx0 and ny0 are the number of points in the input topographic file
! dx and dy are the data spacings in the x- and y-direction in the input topographic
! file (in degrees or meters)
! nskip is the number of points that are skipped when sampling the initial
! topo file
! xlat and xlon are the latitude/longitude of the bottom left corner of the
! data set (in deg. or meters)
! Reads in coordinate flag (coordflag)
! If coordflag=1, the input type is degrees
! If coordflag=2, the input type is utm
      read (55,*) coordflag
      read (55,*) nx0,ny0,dx,dy,nskip,xlon,xlat

! Convert necessary inputs from meters to kilometers if coordinate system is utm
      if(coordflag.eq.2) then
    dx=dx/1000.
    dy=dy/1000.
    xlon=xlon/1000.
    xlat=xlat/1000.
      endif

! nstep is the number of time step
! tau is the erosion time scale (assuming an exponential decvrease in topography) (in Myr)
      read (55,*) nstep,tau

! Compares the number of topography files to the number of time steps
! If there are less topography files than nstep+1, then it assumes the same topography
! for the last time steps
! If there are more topography files, then it quits
! Added 2/6/2007 by cspath
    if ((nstep+1) .gt. num_topo_files) then
      print *, 'Warning: Number of topography files less than number of time steps + 1.'
      print *, 'Assuming same topography for last ',((nstep+1)-num_topo_files),'files.'
    else if ((nstep+1) .lt. num_topo_files) then
      print *, 'Error: Number of topography files exceeds number of time steps + 1'
      print *, 'Reduce files to or below amount of time steps'
      stop
    endif 
    
! for each time step + 1
! timek is the time (in Myr)
! Pe is the exhumation rate (in km/Myr)
! topomag is the topographic amplification factor at this step
! topooffset is the topographic offset

      allocate (timek(nstep+1),Pe(nstep+1),topomag(nstep+1),topooffset(nstep+1))
      allocate (x1f(nstep+1),y1f(nstep+1),x2f(nstep+1),y2f(nstep+1),def(nstep+1),dif(nstep+1))
      allocate (geoflag(nstep+1),theta(nstep+1),phi(nstep+1),thermflag(nstep+1))
      allocate (mft_ratein(nstep+1),mbt_ratein(nstep+1),mct_ratein(nstep+1),stf_ratein(nstep+1))
      allocate (Pe2(nstep+1))

! Reads in time step information from Pecube.in
! Reads in 12 values no matter what geometry is set
! Checks the geometry and then determines which values are needed for later calculations
! Sets all unnecessary variables to 0
        do istep=1,nstep+1

        read (55,*) timek(istep),topomag(istep),topooffset(istep),&
                    thermflag(istep),geoflag(istep),Pe(istep),dummy1,&
                    dummy2,dummy3,dummy4,dummy5,dummy6,dummy7

    if (geoflag(istep) .eq. 1) then
      x1f(istep)=0.
      y1f(istep)=0.
      x2f(istep)=0.
      y2f(istep)=0.
      def(istep)=0.
      dif(istep)=0.
      theta(istep)=0.
      phi(istep)=0.
      mft_ratein(istep)=0.
          mbt_ratein(istep)=0.
      mct_ratein(istep)=0.
      stf_ratein(istep)=0.
	  Pe2(istep)=0.
    else if (geoflag(istep) .eq. 2) then
      theta(istep)=dummy1
      phi(istep)=dummy2
      x1f(istep)=0.
      y1f(istep)=0.
      x2f(istep)=0.
      y2f(istep)=0.
      def(istep)=0.
      dif(istep)=0.
      mft_ratein(istep)=0.
          mbt_ratein(istep)=0.
      mct_ratein(istep)=0.
      stf_ratein(istep)=0.
      Pe2(istep)=0.
    else if (geoflag(istep) .eq. 3) then
      x1f(istep)=dummy1
      y1f(istep)=dummy2
      x2f(istep)=dummy3
      y2f(istep)=dummy4
      def(istep)=dummy5
      dif(istep)=dummy6
      theta(istep)=0.
      phi(istep)=0.
      mft_ratein(istep)=0.
          mbt_ratein(istep)=0.
      mct_ratein(istep)=0.
      stf_ratein(istep)=0.
      Pe2(istep)=0.
    else if (geoflag(istep) .eq. 4) then
      mft_ratein(istep)=dummy1
      mbt_ratein(istep)=dummy2
      mct_ratein(istep)=dummy3
      stf_ratein(istep)=dummy4
      x1f(istep)=0.
      y1f(istep)=0.
      x2f(istep)=0.
      y2f(istep)=0.
      def(istep)=0.
      dif(istep)=0.
      theta(istep)=0.
      phi(istep)=0.
      Pe2(istep)=0.
	else if (geoflag(istep) .eq. 5) then
	  x1f(istep)=dummy1
	  y1f(istep)=dummy2
	  x2f(istep)=dummy3
	  y2f(istep)=dummy4
	  def(istep)=dummy5
	  dif(istep)=dummy6
	  theta(istep)=0.
	  phi(istep)=0.
	  mft_ratein(istep)=0.
          mbt_ratein(istep)=0.
	  mct_ratein(istep)=0.
	  stf_ratein(istep)=0.
	  Pe2(istep)=0.
	else if ((geoflag(istep) .eq. 6) .or. (geoflag(istep) .eq. 7)) then
	  x1f(istep)=dummy1
	  y1f(istep)=dummy2
	  x2f(istep)=dummy3
	  y2f(istep)=dummy4
	  def(istep)=dummy5
	  dif(istep)=dummy6
	  Pe2(istep)=dummy7
	  theta(istep)=0.
	  phi(istep)=0.
	  mft_ratein(istep)=0.
          mbt_ratein(istep)=0.
	  mct_ratein(istep)=0.
	  stf_ratein(istep)=0.
    endif

    enddo

! converts geological time into model time
        do istep=nstep+1,1,-1
          timek(istep)=timek(1)-timek(istep)
        enddo

! isostasy flag (0 no isostasy, 1 isostasy on)
! ** NOTE: rhoc and rhom are now defined below in the thermal properties section **
! young is the elastic plate young modulus (in Pa)
! poisson is poisson's ratio (dimensionless)
! thickness is the elastic thickness of the plate (in km)
! nxiso and nyiso are the resolutions in the x- and y-directions of the grid on
! which the isostatic (flexural) calculations are performed (including the FFT)
! note that these numbers must be powers of two.
      read (55,*) isoflag,young,poisson,thickness,nxiso,nyiso

! crustal thickness is the averaged crustal thickness (i.e. the depth at which the
! temperature is assumed to be constant) (in km)
! nzin is the number of points in the z-direction
! tcond is the thermal conductivity (in W/m K)
! heatcap is the heat capacity (in J/kg K)
! rhoc and rhom are the densities for the crust and mantle, respectively (in kg/m3)
! ** NOTE: these values are also used in the isostatic calculations! **
! -REMOVED- diffusivity is the heat diffusivity (in km2/Myr) -REMOVED-
! ** NOTE: Diffusivity is now calculated below, rather than input **
! tmax is the basal temperature (in C)
! tmsl is the temperature at the top of the model (at z=0)
! tlapse is the lapse rate (or change of temperature with height in the atmosphere)
! (in C/km)
! heatproduction is the rate of heat production (in uW/m^3)

      read (55,*) crustal_thickness,nzin,tcond,heatcap,rhoc,rhom
      read (55,*) tmax,tmsl,tlapse,hpc,efold,hpm,shearint,fric

! Thermal model parameters for the Nepal model geometry
! Each line lists the values below for the Indian shield, Sub-Himalaya
! Lesser Himalaya, Greater Himalaya and Tethyan Himalaya
! Line one is the volumetric heat production [uW/m^3]
! Line two is the thermal conductivity [W/m K]
! Line three is the rock density [kg/m^3]
! Line four is the specific heat capacity [J/kg K]
      read(55,*) ishp,shhp,lhhp,ghhp,thhp
      read(55,*) istc,shtc,lhtc,ghtc,thtc
      read(55,*) isrho,shrho,lhrho,ghrho,thrho
      read(55,*) ishc,shhc,lhhc,ghhc,thhc

! Convert input values to units used in Pecube
      ! crustal thickness (km --> m)
      crustal_thickness=crustal_thickness*1.d3
      ! number of seconds in 1 My
      yrsec=3.15569259747d7
      ! heat production (uW/m^3-->C/My)
      heatproduction=((hpc*1.d-6)/(rhoc*heatcap))*(yrsec*1.e6)
      ishp=((ishp*1.d-6)/(isrho*ishc))*(yrsec*1.e6)                             ! Indian Shield HP
      shhp=((shhp*1.d-6)/(shrho*shhc))*(yrsec*1.e6)                             ! Sub-Himalaya HP
      lhhp=((lhhp*1.d-6)/(lhrho*lhhc))*(yrsec*1.e6)                             ! Lesser Himalaya HP
      ghhp=((ghhp*1.d-6)/(ghrho*ghhc))*(yrsec*1.e6)                             ! Greater Himalaya HP
      thhp=((thhp*1.d-6)/(thrho*thhc))*(yrsec*1.e6)                             ! Tethyan Himalaya HP
      ! thermal diffusivity (km^2/My)
      diffusivity=(tcond/(rhoc*heatcap))*(1.d-3)**2*yrsec*1.e6
      isdiff=(istc/(isrho*ishc))*(1.d-3)**2*yrsec*1.e6                          ! Indian Shield diffusivity
      shdiff=(shtc/(shrho*shhc))*(1.d-3)**2*yrsec*1.e6                          ! Sub-Himalaya diffusivity
      lhdiff=(lhtc/(lhrho*lhhc))*(1.d-3)**2*yrsec*1.e6                          ! Lesser Himalaya diffusivity
      ghdiff=(ghtc/(ghrho*ghhc))*(1.d-3)**2*yrsec*1.e6                          ! Greater Himalaya diffusivity
      thdiff=(thtc/(thrho*thhc))*(1.d-3)**2*yrsec*1.e6                          ! Tethyan Himalaya diffusivity

! obsfile is the name of the observation file
        read (55,*) num_obsfiles
        allocate (obsfile(num_obsfiles),nobsfile(num_obsfiles))
        do j=1,num_obsfiles
          do i=1,300
            obsfile(j)(i:i)=' '
          enddo
          read (55,'(a)') obsfile(j)
          do i=1,300
            if (obsfile(j)(i:i).ne.' ') nobsfile(j)=i
          enddo
        enddo

! Read in flags for which ages to calculate and output
      allocate (age_flags(11))
      read (55,*) age_flags(1),age_flags(2),age_flags(3),&
                  age_flags(4),age_flags(5),age_flags(6),&
                  age_flags(7),age_flags(9),age_flags(8),&
                  age_flags(10),age_flags(11)

      ! Type of detrital calculation, if any
      ! Options are for no PDF calculation, Cascade catchments, and user specified catchments for non-cascade models
      read (55,*) det_calc

      ! Minimum number of nodes for a Cascade model catchment to be output
      read (55,*) node_thresh

      ! Directory where the Cascade output files are located for use in Pecube
      ! Only necessary if option '1' is selected for det_calc
      do i=1,100
        cascadedir(i:i)=' '
      enddo
      if (det_calc .eq. '1') read (55,'(a100)') cascadedir

      close (55)

! Opens all topography files and stores the values in 3D array (zNZ)
! Modified for multiple topography file inputs on 2/7/2007 by cspath
    allocate (zNZ(nstep+1,nx0,ny0))
    if (topo_file_flag.eq.1) then
      if (num_topo_files.eq.0) then
            zNZ=0.d0
      else
        do k=1, nstep+1
          if(k.le.num_topo_files) then
              open (45,file='input/Pecube-D/'//fnme(k)(1:nfnme(k)),status='old',iostat=ios)
              if (ios.ne.0) then
                open (45,file=fnme(k)(1:nfnme(k)),status='old')
              endif
              do m=1,ny0
                do n=1,nx0
                  read (45,*) zNZ(k,n,m)
                enddo
              enddo
          else
          open (45,file='input/Pecube-D/'//fnme(num_topo_files)(1:nfnme(num_topo_files)),status='old',iostat=ios)
              if (ios.ne.0) then
                open (45,file=fnme(num_topo_files)(1:nfnme(num_topo_files)),status='old')
              endif
              do m=1,ny0
                do n=1,nx0
                  read (45,*) zNZ(k,n,m)
                enddo
              enddo
          endif
              close (45)
        enddo
      endif
    elseif (topo_file_flag.eq.2) then
      if (num_topo_files.eq.0) then
            zNZ=0.d0
      else
        do k=0, nstep
          if (k.le.(num_topo_files-1)) then
          write (file_iter,'(i4)') k
          if (k.lt.1000) file_iter(1:1)='0'
              if (k.lt.100) file_iter(1:2)='00'
              if (k.lt.10) file_iter(1:3)='000'
              open (45,file='input/Pecube-D/'//topo_prefix(1:prefix_length)//file_iter//'.dat',status='old',iostat=ios)
              if (ios.ne.0) then
                open (45,file=topo_prefix(1:prefix_length)//file_iter//'.dat',status='old')
              endif
              fnme(k+1)=topo_prefix(1:prefix_length)//file_iter//'.dat'
              read (45,*) zNZ(k+1,:,:)
          else
          open (45,file='input/Pecube-D/'//topo_prefix(1:prefix_length)//file_iter//'.dat',status='old',iostat=ios)
              if (ios.ne.0) then
                open (45,file=topo_prefix(1:prefix_length)//file_iter//'.dat',status='old')
              endif
              read (45,*) zNZ(k+1,:,:)
          endif
              close (45)
        enddo
      endif
    endif

! Stores certain values of the topography files input based on nx, ny, nskip
! Modified for multiple topography file input on 2/7/2007 by cspath
      nx=(nx0-1)/nskip+1
      ny=(ny0-1)/nskip+1
      allocate (z(nstep+1,nx*ny),zoff(nx*ny))

      do k=1, nstep+1
      ij=0
        do j=1,ny0,nskip
          do i=1,nx0,nskip
            ij=ij+1
            z(k,ij)=zNZ(k,i,j)
      enddo
        enddo
      enddo
    
! Modified for multiple topography inputs on 2/7/2007 by cspath
      allocate(zmin(nstep+1),zmax(nstep+1))
      do k=1,nstep+1
        zmin(k)=minval(z(k,:))
    zmax(k)=maxval(z(k,:))
      enddo


! Calculates xl,yl based on coordinate flag (1=Degrees, 2=UTM) (added by cspath)
      if (coordflag .eq. 1) then      
      xl=dx*(nx-1)*nskip*111.11*cos((xlat+dy*ny0/2.)*3.141592654/180.)
      yl=dy*(ny-1)*nskip*111.11
      else if (coordflag .eq. 2) then
      xl=dx*(nx-1)*nskip
      yl=dy*(ny-1)*nskip
      else
      print *, 'Coordinate flag is set to an invalid value'
      endif

      zl=crustal_thickness/1.e3      

! Modified for multiple topography inputs on 2/7/2007 by cspath
    !print *, 'z1:',z(1,1)
      do k=1,nstep+1
        z(k,:)=(z(k,:)-zmin(k))/crustal_thickness*zl
      enddo
    !print *, 'z2:',z(1,1)

      ! Commented out because this was leading to incorrect heat
      !   production values.
      ! dwhipp - 10/07
      !heatproduction=heatproduction*diffusivity/zl**2*tmax

! Converts values to kilometers based on coordinate flag (added by cspath)
! Note: If coordinate system is utm, the values inputted here should be in km as specified in Pecube.in
      if (coordflag .eq. 1) then
        do istep=1,nstep+1
        x1f(istep)=(x1f(istep)-xlon)*111.11*cos((xlat+dy*ny0/2.)*3.141592654/180.)
        y1f(istep)=(y1f(istep)-xlat)*111.11
        x2f(istep)=(x2f(istep)-xlon)*111.11*cos((xlat+dy*ny0/2.)*3.141592654/180.)
        y2f(istep)=(y2f(istep)-xlat)*111.11
        enddo
      else if (coordflag .eq. 2) then
        do istep=1,nstep+1
    x1f(istep)=x1f(istep)-xlon
    y1f(istep)=y1f(istep)-xlat
    x2f(istep)=x2f(istep)-xlon
    y2f(istep)=y2f(istep)-xlat
    enddo
      endif

! Modified for multiple topography inputs on 2/7/2007 by cspath
      do k=1,nstep+1
    zmax(k)=maxval(z(k,:))
      enddo         

      open (7,file='Pecube.dat',status='unknown')
      !open (7,status='scratch')

      ilog=0
      iterative=1
      interpol=1

      nxs=nx/2
! 2011.07.25, WK: non standard format specifier
      write (7,'(a,i10)') run,nrun
      write (7,*) num_topo_files,det_calc,age_flags(1),age_flags(2),age_flags(3),age_flags(4),&
                  age_flags(5),age_flags(6),age_flags(7),age_flags(8),age_flags(9),&
                  age_flags(10),age_flags(11)
      write (7,'(A)') (fnme(k),k=1,num_topo_files)
      write (7,*) 4,nx*ny,nzin,(nx-1)*(ny-1),zl,diffusivity,heatproduction,efold,shearint,fric
      ! Write out Nepal model geometry info
      ! dwhipp 11/07
      write (7,*) isdiff,shdiff,lhdiff,ghdiff,thdiff,ishp,shhp,lhhp,ghhp,lhhp
      write (7,*) tmax,tmsl,tlapse,nstep,ilog,iterative,interpol
      write (7,*) isoflag,tau,rhoc,rhom
      write (7,*) nx,ny,nxiso,nyiso,nx0,ny0,dx,dy
      write (7,*) xl/(nx-1)*1.d3,yl/(ny-1)*1.d3,young,poisson,thickness*1.d3
      write (7,*) xlon,xlon+dx*(nx-1)*nskip,xlat,xlat+dy*(ny-1)*nskip


        do j=1,ny
          do i=1,nx
          x=xl*float(i-1)/float(nx-1)
          y=yl*float(j-1)/float(ny-1)
          write (7,*) x,y
          enddo
        enddo

        do j=1,ny-1
          do i=1,nx-1
          icon1=(j-1)*nx+i
          icon2=icon1+1
          icon3=icon1+nx+1
          icon4=icon1+nx
          write (7,*) icon1,icon2,icon3,icon4
          enddo
        enddo

! Modified for multiple topography inputs on 2/7/2007 by cspath
      do k=1,num_topo_files
    zmax(k)=maxval(z(k,:))
      enddo 

! Modified for multiple topography inputs on 2/7/2007 by cspath
! Writes all time step information to scratch file to be read by Pecube.f90
    do istep=0,nstep
        write (7,*) timek(istep+1),Pe(istep+1),0, &
                    x1f(istep+1),y1f(istep+1),x2f(istep+1),y2f(istep+1), &
                    def(istep+1),dif(istep+1),thermflag(istep+1),geoflag(istep+1),theta(istep+1),phi(istep+1),&
            mft_ratein(istep+1),mbt_ratein(istep+1),mct_ratein(istep+1),stf_ratein(istep+1),Pe2(istep+1)
    write (7,*) (z(istep+1,k)*topomag(istep+1)+topooffset(istep+1),k=1,nx*ny)
        enddo

    if (det_calc .eq. '1') write (7,*) node_thresh
    if (det_calc .eq. '1') write (7,*) cascadedir

! observations
      write (7,*) num_obsfiles
      allocate (age_obs(11),dage_obs(11))
      if (num_obsfiles.eq.0) then
        nobs=0
        write (7,*) nobs
      else if (obsfile(1)(1:nobsfile(1)).eq.'Nil') then
        nobs=0
        write (7,*) nobs
      else
        do j=1,num_obsfiles
          open (8,file='input/Pecube-D/'//obsfile(j)(1:nobsfile(j)),status='old',iostat=ios)
          if (ios.ne.0) then
            open (8,file=obsfile(j)(1:nobsfile(j)),status='old')
          endif

          read (8,*) nobs
          write (7,*) nobs
          do i=1,nobs

            read (8,*) xlonobs,xlatobs,heightobs,He_flag,age_obs(1),dage_obs(1),age_obs(9),dage_obs(9),&
                       age_obs(8),dage_obs(8),age_obs(10),dage_obs(10),age_obs(11),dage_obs(11)

            age_obs(He_flag) = age_obs(1)
            dage_obs(He_flag) = dage_obs(1)

            i1=int((xlonobs-xlon)/(dx*nskip))+1
            if (i1.eq.nx) i1=nx-1
            j1=int((xlatobs-xlat)/(dy*nskip))+1
            if (j1.eq.ny) j1=ny-1
            ieobs=i1+(j1-1)*(nx-1)
            r=(xlonobs-(i1-1)*dx*nskip-xlon)/(dx*nskip)
            r=-1.+2.*r
            s=(xlatobs-(j1-1)*dy*nskip-xlat)/(dy*nskip)
            s=-1.+2.*s
            wobs1=(1.-r)*(1.-s)/4.
            wobs2=(1.+r)*(1.-s)/4.
            wobs3=(1.+r)*(1.+s)/4.
            wobs4=(1.-r)*(1.+s)/4.

            write (7,*) He_flag
            write (7,*) xlonobs,xlatobs,heightobs,wobs1,wobs2,wobs3,wobs4,ieobs,&
                        age_obs(He_flag),dage_obs(He_flag),age_obs(8),dage_obs(8),&
                        age_obs(9),dage_obs(9),age_obs(10),dage_obs(10),&
                        age_obs(11),dage_obs(11)

          enddo
          close (8)
        enddo
      endif

      deallocate (zNZ,z,zoff,zmin,zmax)
      deallocate (timek,Pe,topomag,topooffset)
      deallocate (x1f,x2f,y1f,y2f,def,dif)
      deallocate (geoflag,theta,phi,mft_ratein,mbt_ratein,mct_ratein,stf_ratein,thermflag)
      deallocate (Pe2)
      deallocate (fnme,nfnme)
      deallocate (age_flags)
      deallocate (age_obs,dage_obs)
      deallocate (obsfile,nobsfile)

      return
      end
