! Please quote the following reference(s) when publishing results obtained
! with Pecube:

! Braun , J., 2003. Pecube: A new finite element code to solve the 3D heat
!  transport equation including the effects of a time-varying, finite
!  amplitude surface topography.  Computers and Geosciences, v.29, pp.787-794.

! Braun , J., 2002. Quantifying the effect of recent relief changes on age-elevation
!  relationships.  Earth and Planetary Science Letters, v.200, pp.331-343.

! Braun, J., 2002. Estimating exhumation rate and relief evolution by spectral
!  analysis of age-elevation datasets. Terra Nova, v.14, pp.210-214.

      program Pecube

      implicit real*8 (a-h,o-z)

      real*8,dimension(:),allocatable :: x,y,z,xp,yp,zp,dummy,dummyp,header_info
      real*8,dimension(:),allocatable :: t,tp,f,therm_his_val
      integer,dimension(:,:),allocatable :: icon
      real*8,dimension(:,:,:),allocatable :: ael
      real*8,dimension(:,:),allocatable :: bel
      integer, dimension(:),allocatable :: kfix,ielsurf
      real*8,dimension(:),allocatable :: xsurf,ysurf,zsurf,zsurfp
      real*8,dimension(:),allocatable :: topoa,topob,rsurf,rebound
      real*8,dimension(:),allocatable :: eheat,econd
      integer,dimension(:,:),allocatable :: iconsurf,neighbour
      real*8,dimension(:,:),allocatable :: xdepth,ydepth,zdepth,zsurf_all
      real*8,dimension(:,:),allocatable :: tprev
      real*8, dimension(:,:),allocatable :: xexhumation,yexhumation,zexhumation
      integer,dimension(:),allocatable :: proc,eproc,ice,age_flags
      real*4,dimension(:),allocatable :: t3,t4,t5,t6,edot,bg_edot
      real*4,dimension(:,:),allocatable :: t2,t2p,ztime,ztemp,edot_store
      real*4,dimension(:,:),allocatable :: xstore,ystore,zstore
      real*4,dimension(:,:,:),allocatable :: tshow
      real*4,dimension(:),allocatable :: topo_edot
      real*4,dimension(:,:),allocatable :: fta,ftag,ftm,ftmg,ftz,ftzg
      real*4,dimension(:,:,:),allocatable :: fte,fteg
      real*8 mft_ratein,mbt_ratein,mct_ratein,stf_ratein,tb_velo
      real*8 vx,vy,vz,iout,PecletN,PecletO
      real*4 times(1000),zjunk,t1
      real*8 isdiff,shdiff,lhdiff,ghdiff,thdiff,ishp,shhp,lhhp,ghhp,thhp
      real*8 istc,shtc,lhtc,ghtc,thtc,isrho,shrho,lhrho,ghrho,thrho
      real*8 ishc,shhc,lhhc,ghhc,thhc,cur_depth,xy_mean,fact,zmin,zmax
      real*8 scale_fact,zmin_orig
      integer ie,je,num_therm_his,niter,plane_store1,plane_store2
      integer,dimension(:),allocatable :: jrec
      integer::summ,irec_val,cur_rec,counter,nz,nzin,num_left,sum_ages,num_obsfiles
      character run*100,c3*3,header,det_calc*300,comp*3
      logical interpol_success,obflag,depth_flag
      real*8,dimension(:),allocatable :: age_obs,dage_obs
      real*8,dimension(:),allocatable :: xmisfit_tot,age_prd,sum_xmisfit
      real*8,dimension(:,:),allocatable :: xmisfit
      integer He_flag,nx0,ny0,node_thresh
      real*8 dx,dy,x_basin,y_basin
      real*8,dimension(:),allocatable :: pdf_ages,error
      real*8,dimension(:),allocatable :: x_basin1,y_basin1
      real*8,dimension(:),allocatable :: x_basin2,y_basin2
      character,dimension(:),allocatable :: data_file1*300
      character data_file*300,data_comp*3,cascadedir*100
      integer,dimension(:),allocatable :: age_type1,age_type2
      character,dimension(:),allocatable::fnme*300
      integer dataCount,nilCount,age_type,foundtwo,shearint,nrun,ios
      real*8,dimension(:,:),allocatable :: pages,perates
      character,dimension(:),allocatable :: data_comp1*3,data_comp2*3
      real*8 friction
      real*8 Peclet2

      call cpu_time (times(1))
      nproc=1
      iproc=0

      eps=tiny(eps)

! Pecube is a Finite Element solver of the 3D, transient heat transfer
! equation that allows for conduction, vertical advection and production of
! heat. Pecube also allows for the surface geometry to vary with time.

! This version is an improvement (we hope) of the original version
! available by download from the author's website in that it reads an
! input (topography) file as well as a set of cooling ages (He and FT in
! apatite of known location and calculates synthetic ages at the same
! locations for comparison.

! A flexural isostatic model has also been included to
! calculate the effect of isostasy in amplifying the exhumation caused
! by erosional unloading

! The tectonic uplift has been replaced by movement along a listric
! thrust fault

! To understand how this version works, read the information in the input
! file named Pecube.in

! Note that this input file is read in the subroutine create_pecube_in which
! generates another (scratch) file defined as unit 7 from which Pecube reads
! in the controlling data/information. Note that if you know what you are doing
! you can by pass this operation and create your own input file to be read directly
! by Pecube

      write (6,*) 'Begin Pecube Execution'

      if (iproc.eq.0) write (6,*) 'Reading input'

      call create_pecube_in

      call cpu_time (times(2))

! opens input files

      rewind (7)

! read in general information

! first line:
! run: 5 character string that will determine the name of the folder where the input file
!      will be copied to and where the output files (Pecube.out and Pecube.ptt) will be stored.
! num_topo_files: number of topography files to read
      allocate (age_flags(11))
! 2011.07.25, WK: format specifier not correct
      read (7,'(a,i10)') run,nrun
      read (7,*) num_topo_files,det_calc,age_flags(1),age_flags(2),age_flags(3),&		! modified 01/08
                 age_flags(4),age_flags(5),age_flags(6),age_flags(7),&
                 age_flags(8),age_flags(9),age_flags(10),age_flags(11)
      allocate (fnme(num_topo_files))
      read (7,'(A)') (fnme(k),k=1,num_topo_files)

! second line
! npe: number of nodes per surface (2D) elements (3 = triangular elements - 4 = rectangular elements)
! nsurf: number of nodes defining the surface geometry
! nzin: number of nodes in the vertical direction
! nelemsurf: number of surface 2D elements
! zl: thickness of crustal layer (= depth at which t is fixed) (in km)
! diffusivity: thermal diffusivity (in km^2/Myr)
! heatproduction: heat production (in degC/Myr)
! efold: e-folding depth for decrease in heat production
! Note: HP is constant above msl and decreases exponentially below
! 1/e decrease in HP at efold

      read (7,*) npe,nsurf,nzin,nelemsurf,zl,diffusivity,heatproduction,efold,&
                 shearint,friction

! Read in Nepal model geometry info
! This line lists the thermal diffusivity and heat production values for the
! Indian shield, Sub-Himalaya Lesser Himalaya, Greater Himalaya and Tethyan Himalaya
      read (7,*) isdiff,shdiff,lhdiff,ghdiff,thdiff,ishp,shhp,lhhp,ghhp,thhp
      mpe=2*npe

      allocate (iconsurf(npe,nelemsurf),ielsurf(nsurf))
      allocate (neighbour(npe,nelemsurf))
      allocate (xsurf(nsurf),ysurf(nsurf),zsurf(nsurf),zsurfp(nsurf))
      allocate (topoa(nsurf),topob(nsurf),rsurf(nsurf))

! second line:
! tmax: basal temperature (at z=-zl) (in degC)
! tmsl: temperature at mean sea level (z=0) (in degC)
! tlapse: lapse rate (indegC/km) 
! nstep: number of stages (the surface topo and exhumatoin rate can be specified at each time stage)
! ilog: redundant (dont use)
! iterative: iterative solver method (1 = Gauss Siedel with overrelaxation - 2 = Conjugate gradient)
! interpol: interpolation used (1 = linear - 2 = quadratic, not recommended)

! isoflag: flag for isostasy (0 no isostasy; 1 isostasy)
! tau: erosion time scale which determines the exponential rate qt which topography
! is transformed from one step to the other
! rhoc crustal density
! rhom mantle density
! nx, ny, nxiso,nyiso are the x and y discretization of the FE grid and the isostatic grid
! note that we have lost a bit of generality in Pecube when implementing the isostatic response
! in that we need (for computational efficiency) the FE grid to be rectangular
! note also that nxiso and nyiso need to be powers of 2
! xstep, ystep are the grid spacing
! young and poisson are young modulus and poisson's ratio
! thickness is effective elastic thickness

      read (7,*) tmax,tmsl,tlapse,nstep,ilog,iterative,interpol
      read (7,*) isoflag,tau,rhoc,rhom
      read (7,*) nx,ny,nxiso,nyiso,nx0,ny0,dx,dy
      read (7,*) xstep,ystep,young,poisson,thickness
      if (ilog.eq.1) open (9,file='Pecube.log',status='unknown')

! Note: the *_store arrays currently have a hard coded 1st column length of 500
! This is because this length is determined by how many sub-time intervals Pecube
! calculates for every time step
      allocate (xdepth(nstep,nsurf),ydepth(nstep,nsurf),zdepth(nstep,nsurf),tprev(nstep,nsurf))
      allocate (xexhumation(nstep,nsurf),yexhumation(nstep,nsurf))
      allocate (zexhumation(nstep,nsurf),zsurf_all(nstep,nsurf))
      allocate (therm_his_val(nstep+1))
      allocate (edot_store(nstep,nsurf))

! read in nodal geometry
      read (7,*) xlonmin,xlonmax,xlatmin,xlatmax

! xsurf, ysurf: x- and y-locations of the surface nodes
! iconsurf: connectivity matrix describing the 2D surface elements

      read (7,*) (xsurf(i),ysurf(i),i=1,nsurf)
      read (7,*) ((iconsurf(i,j),i=1,npe),j=1,nelemsurf)

      if (iproc.eq.0) write (6,*) 'Advecting rocks'

      xmin=minval(xsurf)
      xmax=maxval(xsurf)
      ymin=minval(ysurf)
      ymax=maxval(ysurf)

! go through the input file a first time to determine the depth of the
! points that will end up at the surface at the end of the model run

      niteradvec=10
      xexhumation=0.
      yexhumation=0.
      zexhumation=0.
        do istep=nstep,1,-1
        if (iproc.eq.0) write (6,*) 'Time step ',istep
        timesurfp=0.
        rewind (7)
        read (7,*)
        read (7,*) i1,det_calc,i3,i4,x1,x2,x3,x4,x5,x6,x7,x8,x9
        read (7,'(A)') (fnme(k),k=1,num_topo_files)
        read (7,*) i1,i2,i3,i4,x1,x2,x3,x4,x5,x6
        read (7,*) x1,x2,x3,x4,x5,x6,x7,x8,x9,x10
        read (7,*) x1,x2,x3,i1,i2,i3,i4
        read (7,*) x1,x2,x3,x4
        read (7,*) i1,i2,i3,i4,i5,i6,x1,x2
        read (7,*) x1,x2,x3,x4,x5
        read (7,*) x1,x2,x3,x4
        read (7,*) (xsurf(i),ysurf(i),i=1,nsurf)
        read (7,*) ((iconsurf(i,j),i=1,npe),j=1,nelemsurf)
          do kstep=0,istep-1
          read (7,*) timesurfp,Peclet,iout,x1f,y1f,x2f,y2f,def,dif,thermflag,geoflag,&
                     theta,phi,mft_ratein,mbt_ratein,mct_ratein,stf_ratein,Peclet2
          read (7,*) (zsurfp(i),i=1,nsurf)
          enddo
        read (7,*) timesurf,Peclet,iout,x1f,y1f,x2f,y2f,def,dif,thermflag,geoflag,&
                   theta,phi,mft_ratein,mbt_ratein,mct_ratein,stf_ratein,Peclet2
        read (7,*) (zsurf(i),i=1,nsurf)

! Calculate proper z-node spacing if nz is zero
! cspath and dwhipp 11/07
        if (nzin.eq.0) then                                                     ! If number of input z levels is zero, then code will calculate nz and spacing
          zmin=minval(zsurf)                                                    !   of z levels.  z levels will be spaced 1:1 with the input x and y spacing
          zmin_orig=zmin                                                        !   down to ~5 km below the model surface, 3:1 down to ~15 km below the surface
          zmax=maxval(zsurf)                                                    !   and ~9:1 for the rest of the model.
          xy_mean=(xstep+ystep)/2                                               ! This first portion of the variable z spacing code determines the number of
          cur_depth=zmin+zl                                                     !   z levels for the new geometry.
          nz=1
          depth_flag=.true.
          do while (cur_depth.gt.0.)                                            ! Work down from min elevation to base determining number of z levels needed
            nz=nz+1                                                             ! While still in this loop, increment nz
            if (cur_depth.gt.(zmin+zl)-5.) then                                 ! If within 5 km of model top surface, space node planes at xy_mean (1:1)
              cur_depth=cur_depth-(xy_mean/1000.)                               ! Subtract off new node plane spacing from remaining depth range
              plane_store1=nz-1                                                 ! Store number of planes used at 1:1 spacing
            else if (cur_depth.gt.zmin+zl-15) then                              ! If within 15 km of model surface, space node planes at 3*xy_mean (3:1)
              cur_depth=cur_depth-3*(xy_mean/1000.)                             ! Subtract off new node plane spacing from remaining depth range
              plane_store2=nz-1                                                 ! Store number of planes used at 3:1 spacing
            else                                                                ! If greater than 15 km from model surface, use ~9:1 node plane spacing
              if (depth_flag) then                                              ! Calculate node plane spacing on first click through this condition
                depth_flag=.false.                                              ! Set depth_flag to false to avoid repeating this calculation
                num_left=int((cur_depth)/(9*xy_mean/1000.))+1                   ! Number of remaining node planes is equal to the remaining model depth over
                fact=cur_depth/(num_left*(xy_mean/1000.))                       !   the 9:1 spacing increment, plus one.  This yields spacing of <=9:1.
              endif
              if (int(cur_depth-fact*(xy_mean/1000.)).eq.0.) then               ! If near the base of the model, set the cur_depth to zero
                cur_depth=0.
              else                                                              ! Subtract off new node spacing from remaining depth range
                cur_depth=cur_depth-fact*(xy_mean/1000.)
              endif
            endif
          enddo
        else if (nzin.gt.0) then                                                ! If number of input z levels is positive, then use that number for nz
          nz=nzin
        else                                                                    ! Stop if input nz value is negative
          print *,'Error in Pecube.in: nz must be zero or a positive integer'
          stop
        endif

! Stores the thermal flag values into array to be used to determine
! what thermal histories to write out at the end
! Moved so thermal history flag values are stored before time stepping (cspath 10/07)
        therm_his_val(istep)=thermflag

! isostatic rebound
        topoa=zsurfp
        topob=zsurf

          if (isoflag.eq.1) then
          call isostatic_rebound (topoa,topob,xsurf,ysurf,nsurf,rsurf, &
                                  rhoc,rhom,nx,ny,nxiso,nyiso, &
                                  xstep,ystep,young,poisson,thickness)
          else
          rsurf=0.
          endif

      do i=1,nsurf
        zexhumation(istep,i)=zl+zsurf(i)		! Adds the model thickness to the extra elevation
        yexhumation(istep,i)=ysurf(i)
        xexhumation(istep,i)=xsurf(i)
      enddo

        if (istep.eq.nstep) tfinal=timesurf

        ! Use correct Peclet number and diffusivity for Nepal model geometry
        ! Added by dwhipp 11/07
        if (geoflag.eq.4) then
          ! DAVE: Should include uplift in velo!
          tb_velo=(mft_ratein+mbt_ratein+mct_ratein-stf_ratein)
          diffusivityO=diffusivity
          diffusivity=max(isdiff,shdiff,lhdiff,ghdiff,thdiff)
          if (20.-tb_velo.ge.tb_velo) then                                      ! Max velocity is underthrusting in Nepal model geometry
            PecletN=20.-tb_velo
          else                                                                  ! Max velocity is overthrusting in Nepal model geometry
            PecletN=tb_velo
          endif
          PecletO=Peclet
          Peclet=PecletN
        endif

        call find_dt (zl,diffusivity,nsurf,zsurf,zsurfp, &
                      nz,Peclet,timesurf,timesurfp,istep,eps,ilog, &
                      dt,ntime,istatic)

        ! Reset Peclet for Nepal model geometry
        ! dwhipp 11/07
        if (geoflag.eq.4) then
          Peclet=PecletO
          diffusivity=diffusivityO
        endif

! Tracks surface nodes from all topography files back to bottom of model
! Stores these initial postions in xdepth,ydepth,zdepth
	do m=nstep,istep,-1
	  do ktime=1,ntime
          time=timesurf-dt*ktime
            do i=1,nsurf
            xxx=xexhumation(m,i)
            yyy=yexhumation(m,i)
            zzz=zexhumation(m,i)
              do k=1,niteradvec
              call find_velo (xxx,yyy,zzz,dt,Peclet,x1f,y1f,x2f,y2f,def,dif, &
                 time,xmin,xmax,ymin,ymax,zl,vx,vy,vz,geoflag,theta,phi,&
		 mft_ratein,mbt_ratein,mct_ratein,stf_ratein,Peclet2)
              dxxx=xexhumation(m,i)-dt*Peclet*vx/2.-xxx
              dyyy=yexhumation(m,i)-dt*Peclet*vy/2.-yyy
              dzzz=zexhumation(m,i)-dt*Peclet*vz/2.-zzz
              xxx=xxx+dxxx
              yyy=yyy+dyyy
              zzz=zzz+dzzz
              enddo
            xexhumation(m,i)=xexhumation(m,i)-dt*Peclet*vx
            yexhumation(m,i)=yexhumation(m,i)-dt*Peclet*vy
            zexhumation(m,i)=zexhumation(m,i)-dt*Peclet*vz-rsurf(i)/ntime
            enddo
          enddo
        xdepth(m,:)=xexhumation(m,:)
        ydepth(m,:)=yexhumation(m,:)
        zdepth(m,:)=zexhumation(m,:)
        zsurfp=zsurf
        timesurfp=timesurf
        enddo
        enddo

! Deallocations
! Moved 09/07 by cspath
      deallocate (xexhumation)
      deallocate (yexhumation)
      deallocate (zexhumation)

! reset the input file to its proper position

      if (iproc.eq.0) write (6,*) 'initializing'

      rewind (7)
! 2011.07.25, WK: format specifier not correct!
      read (7,'(a,i10)') run,nrun
      read (7,*) num_topo_files,det_calc,age_flags(1),age_flags(2),age_flags(3),&
                 age_flags(4),age_flags(5),age_flags(6),age_flags(7),&
                 age_flags(8),age_flags(9),age_flags(10),age_flags(11)
      read (7,'(A)') (fnme(k),k=1,num_topo_files)
      read (7,*) npe,nsurf,nzin,nelemsurf,zl,diffusivity,heatproduction,efold,&
                 shearint,friction

      ! Read in Nepal model geometry info
      ! This line lists the thermal diffusivity and heat production values for the
      ! Indian shield, Sub-Himalaya Lesser Himalaya, Greater Himalaya and Tethyan Himalaya
      read (7,*) isdiff,shdiff,lhdiff,ghdiff,thdiff,ishp,shhp,lhhp,ghhp,thhp

      read (7,*) tmax,tmsl,tlapse,nstep,ilog,iterative,interpol
      nstep=abs(nstep)
      read (7,*) isoflag,tau,rhoc,rhom
      read (7,*) nx,ny,nxiso,nyiso,nx0,ny0,dx,dy
      read (7,*) xstep,ystep,young,poisson,thickness
      read (7,*) xlonmin,xlonmax,xlatmin,xlatmax
      read (7,*) (xsurf(i),ysurf(i),i=1,nsurf)
      read (7,*) ((iconsurf(i,j),i=1,npe),j=1,nelemsurf)

      nnode=nsurf*nz
      nelem=nelemsurf*(nz-1)
      if (ilog.eq.1) write (9,*) 'nnode/nelem= ',nnode,nelem

! opens output files

! Pecube.out contains the temperature field at the end of each stage

      if (iproc.eq.0) open (8,file=run(1:nrun)//'/Pecube.out',status='unknown',access='direct', &
            recl=4*nnode)

! Pecube.ptt contains the depth-temperture-paths of all surface nodes
      open (10,status='scratch',access='direct', &
            recl=4*(1+nsurf*5))

      allocate (x(nnode),y(nnode),z(nnode),t(nnode))
      allocate (xp(nnode),yp(nnode),zp(nnode),tp(nnode))
      allocate (dummy(nnode),dummyp(nnode))
      allocate (icon(mpe,nelem))
      allocate (kfix(nnode))
      allocate (f(nnode),rebound(nnode))
      allocate (jrec(nstep+1))

! build 3D element connectivity

      ie=0
        do iesurf=1,nelemsurf
          do k=1,(nz-1)
          ie=ie+1
            do kk=1,npe
            icon(kk,ie)=(iconsurf(kk,iesurf)-1)*nz+k
            icon(kk+npe,ie)=icon(kk,ie)+1
            enddo
          enddo
        enddo

      if (ie.ne.nelem) then
      stop 'nelem mismatch'
      endif

! finds processor topology
! this is left over from an ancient parallelized version of Pecube

      allocate (proc(nnode),eproc(nelem))
      call define_proc (nproc,nsurf,nz,proc)
        do ie=1,nelem
        ido=0
          do k=1,mpe
          if (proc(icon(k,ie)).eq.iproc) ido=1
          enddo
        eproc(ie)=ido
        enddo

! allocate reduced number of elemental matrices

      nelemloc=0
        do ie=1,nelem
          if (eproc(ie).eq.1) then
          nelemloc=nelemloc+1
          endif
        enddo

      allocate (ael(mpe,mpe,nelemloc),bel(mpe,nelemloc),ice(nelemloc))

      je=0
        do ie=1,nelem
          if (eproc(ie).eq.1) then
          je=je+1
          ice(je)=ie
          endif
        enddo

      if (ilog.eq.1) then
      write (9,*) 'icon'
        do ie=1,nelem
        write (9,*) (icon(k,ie),k=1,mpe)
        enddo
      endif

! finds neighbour connectivity matrix

      call find_neighbours (iconsurf,neighbour,npe,nelemsurf,nsurf)
      ielsurf=1

! initialize global parameters
! alpha is the time integration parameter

      alpha=0.5
      time=0.
      timesurf=0.
      irec=0
      jrec=0
      niter=0

! allocate erosion rate arrays - dwhipp (08/07)
      allocate(edot(nsurf),bg_edot(nsurf),topo_edot(nsurf))

! allocate heat production and thermal conductivity arrays
! dwhipp (10/07)
      allocate(eheat(nelem),econd(nelem))

      if (iproc.eq.0) write (6,*) 'Start of time stepping'

!*********************************************
!*********************************************
!** begining of surface stepping (stageing) **
!*********************************************
!*********************************************

      do istep=0,nstep

      if (iproc.eq.0) write (6,*) 'Step ',istep

      if (ilog.eq.1) write (9,*) 'istep= ',istep

      if (istep.ne.0) zsurfp=zsurf
      timesurfp=timesurf

! read the step information

! timesurf is the time since the begining of the experiment
! Peclet is the exhumation velocity since the end of the last stage
! iout indicates if the T information is to be saved at the end of the stage

      read (7,*) timesurf,Peclet,iout,x1f,y1f,x2f,y2f,def,dif,thermflag,geoflag,theta,phi, &
                 mft_ratein,mbt_ratein,mct_ratein,stf_ratein,Peclet2
        if (istep.eq.0) then
          if (timesurf.gt.eps.and.iproc.eq.0) then
          write (6,*) 'timesurf= ',timesurf
          write (6,*) 'first topography record must be at time zero ...'
          stop
          endif
        endif
      read (7,*) (zsurf(i),i=1,nsurf)

! Stores the thermal flag values into array to be used to determine
! what thermal histories to write out at the end
        if(istep.ne.0) then
          therm_his_val(istep)=thermflag
        endif

! initial (or zeroth) step

        if (istep.eq.0) then
          zsurfp=zsurf
          in=0
          if (nzin.gt.0) then                                                   ! If the input nz value is positive, build geometry with constant z node plane
            do i=1,nsurf                                                        !   spacing equal to the total model thickness divided by the input nz
              do k=1,nz
                in=in+1
                fact=float(k-1)/float(nz-1)
                if (interpol.eq.2.and.fact.gt.eps) fact=sqrt(fact)
                xp(in)=xsurf(i)
                x(in)=xp(in)
                yp(in)=ysurf(i)
                y(in)=yp(in)
                zh=zsurf(i)+zl
                zp(in)=zh*fact
                kfix(in)=0
                ! Modified to identify top and bottom fixed T B/Cs (original
                ! version below) - dwhipp (09/07)
                !if (k.eq.1 .or. k.eq.nz) kfix(in)=1
                if (k.eq.1) kfix(in)=1
                if (k.eq.nz) kfix(in)=2
              enddo
            enddo
          else                                                                  ! If the input nz value is zero, build geometry with variable z node plane
            zmin=minval(zsurf)                                                  !   spacing.  The spacing is such that from the top surface of the model down
            zmax=maxval(zsurf)                                                  !   to 5 km below it, the spacing is 1:1 with the x and y node spacing (xy_mean)
            do i=1,nsurf                                                        !   Below that, down to 15 km below the surface, the spacing is ~3:1 (3*xy_mean)
              do k=1,nz                                                         !   Below that, the spacing is ~9:1 down to the base of the model
                in=in+1                                                         ! Set the topography scaling factor to zero if the top surface is flat
                if ((zmax-zmin).eq.0.) then
                  scale_fact=0.
                else                                                            ! Define scaling factor to shift nodes beneath the topography
                  scale_fact=(zsurf(i)-zmin)*(1/(zsurf(i)+zl))
                endif
                if (k.gt.nz-(plane_store1)) then
                  if (k.eq.nz) then
                    scale_fact=0.
                    cur_depth=zsurf(i)+zl
                  else
                    cur_depth=cur_depth+(xy_mean/1000.)
                  endif
                else if (k.gt.nz-(plane_store2)) then
                  cur_depth=cur_depth+3*(xy_mean/1000.)
                else
                  if (k.eq.1) then
                    cur_depth=0.
                  else
                    cur_depth=cur_depth+fact*(xy_mean/1000.)
                  endif
                endif
                xp(in)=xsurf(i)
                x(in)=xp(in)
                yp(in)=ysurf(i)
                y(in)=yp(in)
                zp(in)=cur_depth+cur_depth*scale_fact
                kfix(in)=0
                ! Modified to identify top and bottom fixed T B/Cs (original
                !   version below) - dwhipp (09/07)
                !if (k.eq.1 .or. k.eq.nz) kfix(in)=1
                if (k.eq.1) kfix(in)=1
                if (k.eq.nz) kfix(in)=2
              enddo
            enddo
          endif

! calculates initial temperature

          do i=1,nsurf
          tsurf=-zsurf(i)*tlapse+tmsl
            do k=1,nz
            in=(i-1)*nz+k
            zh=zsurf(i)+zl
            tp(in)=tsurf+(tmax-tsurf)*(zh-zp(in))/zh
            enddo
          enddo

          t=tp

          if (ilog.eq.1) then
            write (9,*) 'Nodal geometry'
            do i=1,nnode
              write (9,*) i,xp(i),yp(i),zp(i),tp(i),kfix(i)
            enddo
          endif

        endif

        if(istep.eq.0) then
          jrec(istep+1)=jrec(istep+1)
        else
          jrec(istep+1)=jrec(istep)
        endif

! calculates time step length to satisfy stability and accuracy criteria
! note that an additional criterion must be fullfilled for stability
! but it is independent of the time step: (dx u / kappa) < 1 (the local
! Peclet number must be smaller than 1) where dx is the size of an element
! in a direction where a velocity u is applied, with kappa being the thermal
! diffusivity

      ! Use correct Peclet number for Nepal model geometry
      ! Added by dwhipp 11/07
      if (geoflag.eq.4) then
        tb_velo=(mft_ratein+mbt_ratein+mct_ratein-stf_ratein)
        if (20.-tb_velo.ge.tb_velo) then                                      ! Max velocity is underthrusting in Nepal model geometry
          PecletN=20.-tb_velo
        else                                                                  ! Max velocity is overthrusting in Nepal model geometry
          PecletN=tb_velo
        endif
        PecletO=Peclet
        Peclet=PecletN
      endif

      call find_dt (zl,diffusivity,nsurf,zsurf,zsurfp, &
                    nz,Peclet,timesurf,timesurfp,istep,eps,ilog, &
                    dt,ntime,istatic)

!******************************
!* beginning of time substeps *
!******************************

        do itime=1,ntime
        call cpu_time (times(3))
        time=time+dt
        ftime=float(itime)/float(ntime)
        ftimep=float(itime-1)/float(ntime)
        ftime=(1.-exp(-ftime*tfinal/tau))/(1.-exp(-tfinal/tau))
        ftimep=(1.-exp(-ftimep*tfinal/tau))/(1.-exp(-tfinal/tau))
        if (ilog.eq.1) write (9,*) 'itime,time,Pe= ',itime,time,Peclet
        if (iproc.eq.0) write (6,'(a,i4,a,i4,a,i3,a,i3,a,i6,a)') &
           'Doing time step ',itime,' of ',ntime, &
                 ' in stage ',istep,' of ',nstep,' (',niter,')'

      ! Reset Peclet for Nepal model geometry
      ! dwhipp 11/07
      if (geoflag.eq.4) Peclet=PecletO

! build new node geometry

        z=zp
        do i=1,nsurf
          in=i*nz
          inp=i*nz-1
          zsurf1=zsurfp(i)+(zsurf(i)-zsurfp(i))*ftime+zl
          zsurf2=zsurfp(i)+(zsurf(i)-zsurfp(i))*ftimep+zl
          ! calculate rate of surface change
          ! dwhipp 08/07
          if (istep.eq.0) then
            topo_edot(i)=0.
          else
            topo_edot(i)=(zsurf1-zsurf2)/dt
          endif
          topoa(i)=zsurf2
          topob(i)=zsurf1
          z(in)=z(in)+zsurf1-zsurf2
          z(inp)=z(inp)+(zsurf1-zsurf2)/2.
        enddo

        if (iproc.eq.0) write (6,*) 'Min-Max Topo: ',minval(topoa)-zl,maxval(topoa)-zl

! isostatic rebound

        if (isoflag.eq.1) then
          call isostatic_rebound (topoa,topob,xsurf,ysurf,nsurf,rsurf, &
                                     rhoc,rhom,nx,ny,nxiso,nyiso, &
                                     xstep,ystep,young,poisson,thickness)
        else
          if (iproc.eq.0) write (6,*) 'Min-Max Rebound: ',minval(rsurf),maxval(rsurf)
          rsurf=0.
        endif

        do i=1,nsurf
          do k=1,nz
            rebound(k+(i-1)*nz)=rsurf(i)
          enddo
        enddo

        if (in.ne.nnode) then
          stop 'nnode mismatch'
        endif


! build local FE matrices

        if (iproc.eq.0) write (6,*) 'Building matrix '

        call cpu_time (times(7))
        je=0
        do je=1,nelemloc
          ie=ice(je)
          if (ie.eq.0) print*,'null element'

          call make_matrix (mpe,ael(1,1,je),bel(1,je),icon(1,ie),x,y,z,xp,&
                             yp,zp,kfix,diffusivity,heatproduction,Peclet, &
                             x1f,y1f,x2f,y2f,def,dif,alpha,dt,time,tp,nnode,&
                             istatic,zl,xmin,xmax,ymin,ymax,rebound,geoflag,&
                             theta,phi,mft_ratein,mbt_ratein,mct_ratein,&
                             stf_ratein,tlapse,tmsl,efold,rhoc,je,nelemloc,&
                             eheat,econd,isdiff,shdiff,lhdiff,ghdiff,thdiff,&
                             ishp,shhp,lhhp,ghhp,thhp,shearint,friction,Peclet2)

        enddo

        call cpu_time (times(8))

! build global RHS vector

        f=0.d0
        do je=1,nelemloc
          ie=ice(je)
          do k=1,mpe
            ic=icon(k,ie)
            f(ic)=f(ic)+bel(k,je)
          enddo
        enddo

! solve global FE equations

        if (iproc.eq.0) write (6,*) 'Solving matrix '

        call cpu_time (times(5))
        if (iterative.eq.1) then
          call solve_iterative (mpe,1,ael,f,t,kfix,icon, &
                                 nnode,nelemloc,niter,proc,ice, &
                                 z,tlapse,tmsl,zl)
          if (ilog.eq.1) write (9,*) niter,' iterations'
        elseif (iterative.eq.2) then
          stop 'solution strategy not implemented'
        else
          stop 'solution strategy not implemented'
        endif
        call cpu_time (times(6))

! sends short result to the screen

      if (iproc.eq.0) write (6,*) 'Min-Max Temperature: ',minval(t),maxval(t)

! stretch old grid

      if (nzin.gt.0) then                                                       ! Linearly shift nodes if the input nz value is positive
        do i=1,nsurf
          do k=1,nz
            in=(i-1)*nz+k
            fact=float(k-1)/float(nz-1)
            if (interpol.eq.2.and.fact.gt.eps) fact=sqrt(fact)
            zsurf1=zsurfp(i)+(zsurf(i)-zsurfp(i))*ftime
            zh=zsurf1+zl
            zp(in)=zh*fact
          enddo
        enddo
      else                                                                      ! Shift variably spaced z node planes for new topography
        zmin=minval(zsurfp+(zsurf-zsurfp)*ftime)                                ! Calculate new zmin and zmax for new topo
        zmax=maxval(zsurfp+(zsurf-zsurfp)*ftime)
        do i=1,nsurf
          do k=1,nz
            in=(i-1)*nz+k
            if ((zmax-zmin).eq.0.) then                                         ! Set scale_fact to zero if surface topo is flat
              if (zmin.eq.0.) then
                scale_fact=0.
              else
                scale_fact=(zsurfp(i)+(zsurf(i)-zsurfp(i))*ftime)*(1/zl)        ! Set scale factor to allow for topography shifts
              endif
            else                                                                ! Set scale factor to allow for topography shifts and relief scaling
                scale_fact=((zsurfp(i)+(zsurf(i)-zsurfp(i))*ftime)-zmin_orig)*(1/zl)
            endif
            if (k.gt.nz-(plane_store1)) then
              if (k.eq.nz) then
                cur_depth=(zsurfp(i)+(zsurf(i)-zsurfp(i))*ftime)+zl
                scale_fact=0.
              else
                cur_depth=cur_depth+(xy_mean/1000.)
              endif
            else if (k.gt.nz-(plane_store2)) then
              cur_depth=cur_depth+3.*(xy_mean/1000.)
            else
              if (k.eq.1) then
                cur_depth=0.
              else
                cur_depth=cur_depth+fact*(xy_mean/1000.)
              endif
            endif
            zp(in)=cur_depth+cur_depth*scale_fact
          enddo
        enddo
      endif

! interpolate result onto undeformed mesh

          do i=1,nsurf
            ij=(i-1)*nz+1
            call interpolate (t(ij),tp(ij),z(ij),zp(ij),nz)
          enddo


! The loop below runs through for the zeroth time step but there are only nstep
! amount of depths
! This statement simply allows all the x, y, and z depths to get updated at the zeroth
! time step
	if(istep.eq.0) then
	  last_val=1
	else
	  last_val=istep
	endif

! update ptt
! Updates every surface node for each topography that has not reached the surface yet
! Calculates new x,y,z positions and finds a new temperature at that location
        do j=nstep,last_val,-1
          do i=1,nsurf
            xxx=xdepth(j,i)
            yyy=ydepth(j,i)
            zzz=zdepth(j,i)
            do k=1,niteradvec
              call find_velo (xxx,yyy,zzz,dt,Peclet,x1f,y1f,x2f,y2f,def,dif, &
                             time,xmin,xmax,ymin,ymax,zl,vx,vy,vz,geoflag,theta,phi,&
                             mft_ratein,mbt_ratein,mct_ratein,stf_ratein,Peclet2)
              dxxx=xdepth(j,i)+dt*Peclet*vx/2.-xxx
              dyyy=ydepth(j,i)+dt*Peclet*vy/2.-yyy
              dzzz=zdepth(j,i)+dt*Peclet*vz/2.-zzz
              xxx=xxx+dxxx
              yyy=yyy+dyyy
              zzz=zzz+dzzz
            enddo
            xdepth(j,i)=xdepth(j,i)+Peclet*vx*dt
            ydepth(j,i)=ydepth(j,i)+Peclet*vy*dt
            zdepth(j,i)=zdepth(j,i)+Peclet*vz*dt+rsurf(i)
            ! calculate erosion rates at last itime step of each ntime
            ! added by dwhipp 08/07
            if (itime.eq.ntime) then
              bg_edot(i)=Peclet*vz
              edot(i)=bg_edot(i)-topo_edot(i)
            endif
            call find_element (xdepth(j,i),ydepth(j,i),zdepth(j,i),tnow,x,y,z,t, &
                                xsurf,ysurf,zsurf,ielsurf(i),neighbour,iconsurf, &
                                icon,nelemsurf,nelem,nsurf,nz,nnode,npe,mpe, &
                                interpol_success,obflag)

            ! Modified interpolation detection scheme to find whether points exit side/base of model versus not being found on model surface
            ! cspath/dwhipp 11/07
            if (interpol_success) then                                          ! If we found the tracked particle in an element, assign T
              tprev(j,i)=tnow
            else                                                                ! Unable to find tracked particle within an element
              if (obflag) then                                                  ! Point is out side or base of model
                tprev(j,i)=-tmax
              else                                                              ! Point is not found on upper model surface, assign surf T
                tprev(j,i)=-(zdepth(j,i)-zl)*tlapse+tmsl
              endif
            endif
          enddo
        enddo

      if (iproc.eq.0) write (6,*) 'PTt Updated'
! saves tt and pt
! Writes the x,y,z depths to the scratch file instead of storing to array (cspath 10/07)
! Writes out only thermal histories specified in Pecube.in (cspath 10/07)
        do i=1,nstep
          if (therm_his_val(i).eq.1) then
            jrec(istep+1)=jrec(istep+1)+1
            write (10,rec=jrec(istep+1)) sngl(tfinal-time),sngl(tprev(i,:)), &
                  sngl(xdepth(i,:)),sngl(ydepth(i,:)),sngl(zdepth(i,:)), &
                  sngl(zl+zsurfp+(zsurf-zsurfp)*ftime-zdepth(i,:))
          endif
        enddo

!************************
!* end of time stepping *
!************************

        call cpu_time (times(4))

        if (iproc.eq.0) write(6,*) 'This time step : ',times(4)-times(3)

        enddo

! write output

        if (iout.eq.1) then
        if (iproc.eq.0) write (6,*) 'Saving at step ',istep
        dummy=0.d0
          do i=1,nnode
          if (proc(i).eq.iproc) dummy(i)=zp(i)
          enddo
        dummyp=0.d0
        dummyp=dummy
        irec=irec+1
        if (iproc.eq.0) write (8,rec=irec) sngl(dummyp)
        dummy=0.d0
          do i=1,nnode
          if (proc(i).eq.iproc) dummy(i)=tp(i)
          enddo
        dummyp=0.d0
        dummyp=dummy
        irec=irec+1
        if (iproc.eq.0) write (8,rec=irec) sngl(dummyp)
        endif

! note that the following output creates a file containing the temperature
! in a format that can be read by an array visualizer

      if (iproc.eq.0) then
      allocate (tshow(nx,ny,nz))
        do k=1,nz
          do j=1,ny
            do i=1,nx
              tshow(i,j,k)=t(k+(i-1)*nz+(j-1)*(nx*nz))
            enddo
          enddo
        enddo
      continue
      write (c3,'(i3)') istep
      if (istep.lt.100) c3(1:1)='0'
      if (istep.lt.10) c3(1:2)='00'
      open (31,file=run(1:nrun)//'/Time'//c3//'.txt',status='unknown')
      write (31,'(a)') '#agl ascii data V1.0'
      write (31,'(a)') '#type AGL_FLOAT'
      write (31,'(a)') '#lang Fortran'
      write (31,'(a4,3i5)') '#dim',nx,ny,nz
      do i=1,nx
      do j=1,ny
      write (31,'(101f7.1)') (tshow(i,j,k),k=1,nz)
      enddo
      write (31,*)
      enddo
      close (31)
      ! Write out Tecplot formatted output files
      call tec_mat_output (x,y,z,t,nnode,icon,nelem,mpe,run,&
                     c3,dt,Peclet,x1f,y1f,x2f,y2f,def,dif,time,&
                     xmin,xmax,ymin,ymax,zl,vx,vy,vz,geoflag,theta,phi,&
                     mft_ratein,mbt_ratein,mct_ratein,stf_ratein,xlonmin,&
                     xlatmin,xlonmax,xlatmax,nrun,Peclet2)

      ! Call subroutine to write out surface erosion rates
      ! dwhipp - (08/07)
      call erates (nsurf,xxx,xlonmin,xlonmax,xsurf,xmin,xmax,yyy,&
                   xlatmin,xlatmax,ysurf,ymin,ymax,zsurf,run,edot,c3,&
                   bg_edot,topo_edot,nrun,nx,ny,nelemsurf)

      deallocate (tshow)
      endif

! Stores current elevation values for later output
      if (istep.ne.0) zsurf_all(istep,:)=zsurf

      if (istep.ne.0) edot_store(istep,:) = edot

!*****************************
!** end of surface stepping **
!*****************************

      enddo

! Deallocations
! Moved 09/07 by cspath
      deallocate (ielsurf,neighbour)
      deallocate (topoa,topob,rsurf,rebound)
      deallocate (x,y,z,t)
      deallocate (xp,yp,zp,tp,dummy,dummyp)
      deallocate (icon,ael,bel)
      deallocate (kfix,f)

! Finds the number of requested thermal history outputs (specified in Pecube.in)
! Added 10/07 by cspath
      num_therm_his=0
      do k=1,nstep
	if (therm_his_val(k).eq.1) num_therm_his=num_therm_his+1
      enddo

! clean up
! ie we write a last record to the unit 10 file
! Made t2 and t2p 2D to hold temperatures for surface nodes of each topography
      allocate(t2(nstep,nsurf),t2p(nstep,nsurf),t3(nsurf),t4(nsurf),t5(nsurf),&
               t6(nsurf))
      t2p=0.
        do krec=jrec(nstep+1),1,-1

! Added index counter to keep track of temperatures for surface nodes
! of the different surface topographies
	if(mod(krec,num_therm_his).eq.0) then
	  counter=num_therm_his
	else
	  counter=counter-1
        endif

        read (10,rec=krec) t1,t2(counter,:),t3,t4,t5,t6

          do i=1,nsurf
            !***
            ! Modified to fix thermal histories - dwhipp (10/07)
            ! Old version below
            !*****
            !if (t2(counter,i).lt.0.) t2(counter,i)=t2p(counter,i)
            if (t2(counter,i).eq.-tmax) t2(counter,i)=t2p(counter,i)
            t2p(counter,i)=t2(counter,i)
          enddo
        write (10,rec=krec) t1,t2(counter,:),t3,t4,t5,t6
        enddo
      deallocate (t2,t2p,t3,t4,t5,t6)

      summ=jrec(nstep+1)
      do i=1,nstep
      if (therm_his_val(i).eq.1) then
      summ=summ+1
      write (10,rec=summ) sngl(-1.d0),sngl(tprev(i,:)), &
           sngl(xdepth(i,:)),sngl(ydepth(i,:)),sngl(zdepth(i,:)), &
           sngl(zl+zsurfp+(zsurf-zsurfp)*ftime-zdepth(i,:))
      endif
      enddo

      if (iproc.eq.0) close (8)
      if (ilog.eq.1) close (9)

      write (6,*)
      write (6,*) 'Calculating ages and thermal history...'
      write (6,*) 'Please be patient. This may take awhile.'

! calculate FT ages
! Made age arrays 2,3-D to hold the age values of the surface nodes for
! different topography files as they are tracked
! fte is the Helium ages array, fta is fission track ages
! ftm is the muscovite ages
! The fte array is 3-D because of the extra Helium ages calculated within Mad_He
! in addition to finding ages for every time step
! xstore,ystore,zstore hold the x,y,z, depths read from the scratch file (cspath 10/07)
      allocate (fta(nstep,nsurf),fte(nstep,nsurf,8),ftag(nstep,nsurf),fteg(nstep,nsurf,8),&
               ftm(nstep,nsurf),ftmg(nstep,nsurf),header_info(6),ftz(nstep,nsurf),ftzg(nstep,nsurf))

      fta=0.
      fte=0.
      ftm=0.
      ftz=0.
      cnt=0
      do m=1,nstep
      if (therm_his_val(m).eq.1) then
          write (6,*) 'Timestep ',m
          cnt=cnt+1
	  irec_val=jrec(m+1)
          num_recs = irec_val/num_therm_his
          allocate (ztime(nsurf,num_recs),ztemp(nsurf,num_recs))
	  allocate (xstore(nsurf,num_recs),ystore(nsurf,num_recs),zstore(nsurf,num_recs))
        do i=1,nsurf
          ij=(i-1)*nz+1
 	  cur_rec=0
          if (proc(ij).eq.iproc) then
            do irec=cnt,irec_val,num_therm_his
	      cur_rec=cur_rec+1
              read (10,rec=irec) ztime(i,cur_rec),(ztemp(i,cur_rec),k=1,i),(zjunk,k=i+1,nsurf), &
                                  (xstore(i,cur_rec),k=1,i),(zjunk,k=i+1,nsurf), &
                                  (ystore(i,cur_rec),k=1,i),(zjunk,k=i+1,nsurf),(zstore(i,cur_rec),k=1,i), &
                                  (zjunk,k=i+1,nsurf),(zjunk,k=1,nsurf)
	    enddo

            sum_ages=0
            do j=1,8
              sum_ages = sum_ages + age_flags(j)
            enddo

            if (sum_ages.gt.0) then
              call Mad_He (ztime(i,:),ztemp(i,:),num_recs,fte(m,i,:),m,i,header_info,age_flags)
            endif

            if (age_flags(9).eq.1) then
              call Mad_Trax (ztime(i,:),ztemp(i,:),num_recs,0,2,&
                             fta(m,i),ftld,ftldmean,ftldsd)
            else
              fta(m,i)=0.
            endif

            if (age_flags(10).eq.1) then
	      call ZFT (ztemp(i,:),ztime(i,:),num_recs,ftz(m,i))
            else
              ftz(m,i)=0.
            endif

            if (age_flags(11).eq.1) then
              call Muscovite (ztemp(i,:),ztime(i,:),num_recs,ftm(m,i))
            else
              ftm(m,i)=0.
            endif

        endif
        enddo

       call therm_output(ztime,ztemp,num_recs,nsurf,xstore,&
                         ystore,zstore,run,nstep,m,nrun)

        deallocate (ztime,ztemp,xstore,ystore,zstore)
      endif
      enddo

      ftag=fta
      fteg=fte
      ftmg=ftm
      ftzg=ftz

      ! Write out ages for all timesteps

      call ages_header (nsurf,xlonmin,xlonmax,xsurf,xmin,xmax,xlatmin,xlatmax,ysurf,ymin,&
                          ymax,zsurf_all,fteg,ftag,run,nstep,ftmg,header_info,ftzg,&
                          therm_his_val,age_flags,nrun,nx,ny,nelemsurf)

      ! Checks if user specified to have detrital age calculations
      ! If det_calc is 1, then PDFs for all catchments that drain to edge of model are calculated
      ! Note: The catchments_output subroutine will only work for models from cascade because it uses
      ! the catchments found from a cascade output file (topo_tec*.dat)
      if (det_calc .eq. '1') then
        write (6,*)
        write (6,*) 'Calculating Cascade catchment ages and PDFs...'
        read (7,*) node_thresh
        read (7,'(a100)') cascadedir
        call catchments_output (ftag,fteg,ftmg,ftzg,nstep,nsurf,run,xdepth,ydepth,det_calc,&
                                edot_store,nrun,cascadedir,node_thresh,age_flags,header_info,&
                                fnme,num_topo_files)
      elseif (det_calc .ne. '0') then                                                                           ! If user specified a file with specific basin outlets
        write (6,*)
        write (6,*) 'Calculating PDFs for user defined basins...'
        open (104,file='input/Pecube-D/'//det_calc,status='old',iostat=ios)
        if (ios.ne.0) then
          open (104,file=det_calc,status='old')
        endif
        num_basins=0
        do
          read (104,*,end=20)
          num_basins=num_basins+1
        enddo
20      rewind(104)
        allocate (x_basin1(num_basins),y_basin1(num_basins),data_file1(num_basins),age_type1(num_basins))
        allocate (x_basin2(num_basins),y_basin2(num_basins),age_type2(num_basins))
        allocate (data_comp1(num_basins),data_comp2(num_basins))
        nilCount=0
        dataCount=0
        do j=1,num_basins                                                                                       ! Iterates through basin list file
          read (104,*) x_basin,y_basin,age_type,data_file,data_comp                                         ! x position of outlet, y position of outlet, age type (1-8), associated data file or Nil
          if (data_file.eq.'Nil') then
            nilCount=nilCount+1
            x_basin2(nilCount)=x_basin
            y_basin2(nilCount)=y_basin
            age_type2(nilCount)=age_type
            data_comp2(nilCount)=data_comp
          else
            dataCount=dataCount+1
            x_basin1(dataCount)=x_basin
            y_basin1(dataCount)=y_basin
            age_type1(dataCount)=age_type
            data_file1(dataCount)=data_file
            data_comp1(dataCount)=data_comp
          endif
        enddo
        close (104)
        if (nilCount .ne. 0) then
          allocate (pages(nilCount,nx0*ny0),perates(nilCount,nx0*ny0))
          call find_upstream_points (fnme,num_topo_files,x_basin2,y_basin2,nstep,nx0,ny0,dx,dy,&
               run,fteg,ftag,ftzg,ftmg,age_type2,nsurf,xdepth,ydepth,edot_store,&
               xlonmin,xlatmin,xlonmax,xlatmax,zsurf_all,xmin,xmax,ymin,ymax,nilCount,pages,perates,nrun)
        endif
        if (dataCount .ne. 0) then                                                                                                 ! If data file specified, read in data and error and output PDF
          do i=1,dataCount
            open (105,file='input/Pecube-D/'//data_file1(i),status='old',iostat=ios)
            if (ios.ne.0) then
              open (105,file=data_file1(i),status='old')
            endif
            number=0
            do
              read (105,*,end=21)
              number=number+1
            enddo
21          rewind(105)
            allocate (pdf_ages(number),error(number))
            do j=1,number
              read (105,*) pdf_ages(j),error(j)
            enddo
            call pdfmaker_for_data (pdf_ages,error,number,age_type1(i),run,x_basin1(i),y_basin1(i),nrun)
            if (data_comp1(i) .eq. 'yes' .or. data_comp1(i) .eq. 'Yes') then
              foundtwo=0
              do j=1,nilCount
                if (x_basin2(j).eq.x_basin1(i) .and. y_basin2(j).eq.y_basin1(i) .and. age_type2(j).eq.age_type1(i)) then
                  if (data_comp2(j) .eq. 'yes' .or. data_comp2(j) .eq. 'Yes') then
                    foundtwo=1
                    call detrital_mc (data_file1(i),number,pdf_ages,error,&
                         pages(j,:),perates(j,:),nx0*ny0,run,x_basin1(i),y_basin1(i),nrun)
                  endif
                endif
              enddo
              if (foundtwo .eq. 0) then
! 2011.07.25, WK: format specifier not correct!
                write (6,'(A,f8.4,a,f8.4,A)') 'Warning: Could not find second basin data for basin:',x_basin1(i),',',y_basin1(i),' to run Monte Carlo test'
                write (6,'(A)') 'Check that the two basins have the same x position, y position, age type, and that both have "Yes" for the Monte Carlo test to run'
                write (6,'(A)') 'Skipping Monte Carlo test for this basin'
              endif
            endif
            deallocate (pdf_ages,error)
            close (105)
          enddo
        endif
!         enddo
      endif

! writes ages in various format
        if (iproc.eq.0) then
        open (12,file=run(1:nrun)//'/Ages.txt',status='unknown')
        write (12,*) nx,ny
        open (13,file=run(1:nrun)//'/He.txt',status='unknown')
        write (13,'(a)') '#agl ascii data V1.0'
        write (13,'(a)') '#type AGL_FLOAT'
        write (13,'(a)') '#lang Fortran'
        write (13,'(a4,3i5)') '#dim',nx,ny
        open (14,file=run(1:nrun)//'/FT.txt',status='unknown')
        write (14,'(a)') '#agl ascii data V1.0'
        write (14,'(a)') '#type AGL_FLOAT'
        write (14,'(a)') '#lang Fortran'
        write (14,'(a4,3i5)') '#dim',nx,ny
        open (15,file=run(1:nrun)//'/Height.txt',status='unknown')
        write (15,'(a)') '#agl ascii data V1.0'
        write (15,'(a)') '#type AGL_FLOAT'
        write (15,'(a)') '#lang Fortran'
        write (15,'(a4,3i5)') '#dim',nx,ny
          do i=1,nsurf
          xxx=xlonmin+(xlonmax-xlonmin)*(xsurf(i)-xmin)/(xmax-xmin)
          yyy=xlatmin+(xlatmax-xlatmin)*(ysurf(i)-ymin)/(ymax-ymin)
          write (12,'(5f12.4)') xxx,yyy,zsurf(i),fteg(nstep,i,1),ftag(nstep,i)
	  enddo
          do i=1,nx
          write (13,'(501f7.1)') (fteg(nstep,(j-1)*nx+i,1),j=1,ny)
          write (14,'(501f7.1)') (ftag(nstep,(j-1)*nx+i),j=1,ny)
          write (15,'(501f7.1)') (zsurf((j-1)*nx+i),j=1,ny)
          enddo
        open (16,file=run(1:nrun)//'/HeHeight.txt',status='unknown')
        open (17,file=run(1:nrun)//'/FTHeight.txt',status='unknown')
        write (16,'(2f12.4)') (fteg(nstep,i,1),zsurf(i),i=1,nsurf)
        write (17,'(2f12.4)') (ftag(nstep,i),zsurf(i),i=1,nsurf)
        close (12)
        close (13)
        close (14)
        close (15)
        close (16)
        close (17)

        endif

! calculate misfit
! Modified to calculate misfit values for all different types of ages
! Also allows for multiple data files to be entered in Pecube.in
! Will give a Comparison file for each data file and a total
! of all data files in Comparison.txt
      allocate (age_obs(11),dage_obs(11),age_prd(11))
      if (iproc.eq.0) then
      read (7,*) num_obsfiles
      allocate (xmisfit(num_obsfiles,11),xmisfit_tot(num_obsfiles))
      if (num_obsfiles.ne.0) then
      open (14,file=run(1:nrun)//'/Comparison.txt',status='unknown')
      do j=1,num_obsfiles

        write (comp,'(i3)') j
        if (j.lt.100) comp(1:1)='0'
        if (j.lt.10) comp(1:2)='00'

        open (13,file=run(1:nrun)//'/Comparison_'//comp//'.txt',status='unknown')
        xmisfit(j,:)=0.

        read (7,*) nobs
        write (13,*) 'Number of observations: ',nobs
        write (14,*) 'Number of observations: ',nobs

        write (13,'(A55)',ADVANCE="no") ' Data Set Longitude  Latitude   Height Obs. Height Int.'
        write (13,'(A49)',ADVANCE="no") ' AHe Type AHe Age Obs. AHe Age Error AHe Age Prd.'
        write (13,'(A40)',ADVANCE="no") ' AFT Age Obs. AFT Age Error AFT Age Prd.'
        write (13,'(A40)',ADVANCE="no") ' ZHe Age Obs. ZHe Age Error ZHe Age Prd.'
        write (13,'(A40)',ADVANCE="no") ' ZFT Age Obs. ZFT Age Error ZFT Age Prd.'
        write (13,'(A40)',ADVANCE="no") ' MAr Age Obs. MAr Age Error MAr Age Prd.'
        write (13,*)

        write (14,'(A55)',ADVANCE="no") 'Data Set Longitude  Latitude   Height Obs. Height Int.'
        write (14,'(A49)',ADVANCE="no") ' AHe Type AHe Age Obs. AHe Age Error AHe Age Prd.'
        write (14,'(A40)',ADVANCE="no") ' AFT Age Obs. AFT Age Error AFT Age Prd.'
        write (14,'(A40)',ADVANCE="no") ' ZHe Age Obs. ZHe Age Error ZHe Age Prd.'
        write (14,'(A40)',ADVANCE="no") ' ZFT Age Obs. ZFT Age Error ZFT Age Prd.'
        write (14,'(A40)',ADVANCE="no") ' MAr Age Obs. MAr Age Error MAr Age Prd.'
        write (14,*)

        do iobs=1,nobs
          read (7,*) He_flag
          read (7,*) xlonobs,xlatobs,heightobs,wobs1,wobs2,wobs3,wobs4,ieobs,&
                     age_obs(He_flag),dage_obs(He_flag),age_obs(8),dage_obs(8),&
                     age_obs(9),dage_obs(9),age_obs(10),dage_obs(10),&
                     age_obs(11),dage_obs(11)

          if (age_flags(He_flag).eq.0 .and. age_obs(He_flag).gt.0) then
            print *, 'Warning: Apatite Helium age, type:',He_flag,'not being predicted; cannot compare to observed age'
            print *, 'For comparison enter 1 for prediction calculation in Pecube.in'
          endif

          if (iconsurf(1,ieobs).gt.nsurf .or. iconsurf(2,ieobs).gt.nsurf .or. iconsurf(3,ieobs).gt.nsurf .or. iconsurf(4,ieobs).gt.nsurf) then
            print *, 'Error: Surface node connectivity value out of range'
            print *, 'Check that thermochronological data file coordinates are same system (degrees/utm) as model'
          endif

          hei=wobs1*zsurf(iconsurf(1,ieobs)) &
             +wobs2*zsurf(iconsurf(2,ieobs)) &
             +wobs3*zsurf(iconsurf(3,ieobs)) &
             +wobs4*zsurf(iconsurf(4,ieobs))
          hei=hei*1000.

          write (13,'(i5,4f12.4)',ADVANCE="no") j,xlonobs,xlatobs,heightobs,hei
          write (14,'(i5,4f12.4)',ADVANCE="no") j,xlonobs,xlatobs,heightobs,hei

              age_prd(He_flag)=wobs1*fteg(nstep,iconsurf(1,ieobs),He_flag) &
                     +wobs2*fteg(nstep,iconsurf(2,ieobs),He_flag) &
                     +wobs3*fteg(nstep,iconsurf(3,ieobs),He_flag) &
                     +wobs4*fteg(nstep,iconsurf(4,ieobs),He_flag)

              age_prd(8)=wobs1*fteg(nstep,iconsurf(1,ieobs),8) &
                     +wobs2*fteg(nstep,iconsurf(2,ieobs),8) &
                     +wobs3*fteg(nstep,iconsurf(3,ieobs),8) &
                     +wobs4*fteg(nstep,iconsurf(4,ieobs),8)

              age_prd(9)=wobs1*ftag(nstep,iconsurf(1,ieobs)) &
                     +wobs2*ftag(nstep,iconsurf(2,ieobs)) &
                     +wobs3*ftag(nstep,iconsurf(3,ieobs)) &
                     +wobs4*ftag(nstep,iconsurf(4,ieobs))

              age_prd(10)=wobs1*ftzg(nstep,iconsurf(1,ieobs)) &
                     +wobs2*ftzg(nstep,iconsurf(2,ieobs)) &
                     +wobs3*ftzg(nstep,iconsurf(3,ieobs)) &
                     +wobs4*ftzg(nstep,iconsurf(4,ieobs))

              age_prd(11)=wobs1*ftmg(nstep,iconsurf(1,ieobs)) &
                     +wobs2*ftmg(nstep,iconsurf(2,ieobs)) &
                     +wobs3*ftmg(nstep,iconsurf(3,ieobs)) &
                     +wobs4*ftmg(nstep,iconsurf(4,ieobs))

              if (age_prd(He_flag).eq.0) age_prd(He_flag)=-999
              if (age_prd(8).eq.0) age_prd(8)=-999
              if (age_prd(9).eq.0) age_prd(9)=-999
              if (age_prd(10).eq.0) age_prd(10)=-999
              if (age_prd(11).eq.0) age_prd(11)=-999

              write (13,'(i7,3f14.4)',ADVANCE="no") He_flag,age_obs(He_flag),dage_obs(He_flag),age_prd(He_flag)
              write (13,'(3f13.4)',ADVANCE="no") age_obs(9),dage_obs(9),age_prd(9)
              write (13,'(3f13.4)',ADVANCE="no") age_obs(8),dage_obs(8),age_prd(8)
              write (13,'(3f13.4)',ADVANCE="no") age_obs(10),dage_obs(10),age_prd(10)
              write (13,'(3f14.4)',ADVANCE="no") age_obs(11),dage_obs(11),age_prd(11)

              write (14,'(i7,3f14.4)',ADVANCE="no") He_flag,age_obs(He_flag),dage_obs(He_flag),age_prd(He_flag)
              write (14,'(3f13.4)',ADVANCE="no") age_obs(9),dage_obs(9),age_prd(9)
              write (14,'(3f13.4)',ADVANCE="no") age_obs(8),dage_obs(8),age_prd(8)
              write (14,'(3f13.4)',ADVANCE="no") age_obs(10),dage_obs(10),age_prd(10)
              write (14,'(3f14.4)',ADVANCE="no") age_obs(11),dage_obs(11),age_prd(11)

            if (age_obs(He_flag).gt.0. .and. age_flags(He_flag).eq.1) then
              xmisfit(j,He_flag)=xmisfit(j,He_flag)+(age_obs(He_flag)-age_prd(He_flag))**2/dage_obs(He_flag)**2
            endif

            if (age_obs(9).gt.0. .and. age_flags(9).eq.1) then
              xmisfit(j,9)=xmisfit(j,9)+(age_obs(9)-age_prd(9))**2/dage_obs(9)**2
            endif

            if (age_obs(8).gt.0. .and. age_flags(8).eq.1) then
              xmisfit(j,8)=xmisfit(j,8)+(age_obs(8)-age_prd(8))**2/dage_obs(8)**2
            endif

            if (age_obs(10).gt.0. .and. age_flags(10).eq.1) then
              xmisfit(j,10)=xmisfit(j,10)+(age_obs(10)-age_prd(10))**2/dage_obs(10)**2
            endif

            if (age_obs(11).gt.0. .and. age_flags(11).eq.1) then
              xmisfit(j,11)=xmisfit(j,11)+(age_obs(11)-age_prd(11))**2/dage_obs(11)**2
            endif

            write (13,*)
            write (14,*)

          enddo

        if (age_flags(8).eq.0 .and. age_obs(8).gt.0) then
          print *, 'Warning: Zircon Helium age not being predicted; cannot compare to observed age'
          print *, 'For comparison enter 1 for prediction calculation in Pecube.in'
        endif
        if (age_flags(9).eq.0 .and. age_obs(9).gt.0) then
          print *, 'Warning: Apatite Fission Track age not being predicted; cannot compare to observed age'
          print *, 'For comparison enter 1 for prediction calculation in Pecube.in'
        endif
        if (age_flags(10).eq.0 .and. age_obs(10).gt.0) then
          print *, 'Warning: Zircon Fission Track age not being predicted; cannot compare to observed age'
          print *, 'For comparison enter 1 for prediction calculation in Pecube.in'
        endif
        if (age_flags(11).eq.0 .and. age_obs(11).gt.0) then
          print *, 'Warning: Muscovite age not being predicted; cannot compare to observed age'
          print *, 'For comparison enter 1 for prediction calculation in Pecube.in'
        endif

        xmisfit(j,1)=xmisfit(j,1)+xmisfit(j,2)+xmisfit(j,3)+xmisfit(j,4)+xmisfit(j,5)+xmisfit(j,6)+xmisfit(j,7)
        xmisfit_tot(j)=xmisfit(j,1)+xmisfit(j,8)+xmisfit(j,9)+xmisfit(j,10)+xmisfit(j,11)
        if (nobs.gt.0) xmisfit_tot(j)=sqrt(xmisfit_tot(j))

        if (xmisfit(j,1).gt.0.) then
            if (nobs.gt.0) xmisfit(j,1) = sqrt(xmisfit(j,1))
            write (13,*) "Misfit AHe:",xmisfit(j,1)
            write (14,*) "Misfit AHe:",xmisfit(j,1)
        endif

        if (xmisfit(j,9).gt.0.) then
            if (nobs.gt.0) xmisfit(j,9) = sqrt(xmisfit(j,9))
            write (13,*) "Misfit AFt:",xmisfit(j,9)
            write (14,*) "Misfit AFt:",xmisfit(j,9)
        endif

        if (xmisfit(j,8).gt.0.) then
            if (nobs.gt.0) xmisfit(j,8) = sqrt(xmisfit(j,8))
            write (13,*) "Misfit ZHe:",xmisfit(j,8)
            write (14,*) "Misfit ZHe:",xmisfit(j,8)
        endif

        if (xmisfit(j,10).gt.0.) then
            if (nobs.gt.0) xmisfit(j,10) = sqrt(xmisfit(j,10))
            write (13,*) "Misfit ZFt:",xmisfit(j,10)
            write (14,*) "Misfit ZFt:",xmisfit(j,10)
        endif

        if (xmisfit(j,11).gt.0.) then
            if (nobs.gt.0) xmisfit(j,11) = sqrt(xmisfit(j,11))
            write (13,*) "Misfit MAr:",xmisfit(j,11)
            write (14,*) "Misfit MAr:",xmisfit(j,11)
        endif

        write (13,*) 'Total Misfit: ',xmisfit_tot(j)
        write (14,*) 'Total Misfit: ',xmisfit_tot(j)
        write (14,*)

        close (13)
      enddo

      allocate (sum_xmisfit(6))
      write (14,*) 'Total Misfits for all data sets:'

      sum_xmisfit=0.0
      do n=1,num_obsfiles
        sum_xmisfit(1) = sum_xmisfit(1) + xmisfit(n,1)
        sum_xmisfit(2) = sum_xmisfit(2) + xmisfit(n,9)
        sum_xmisfit(3) = sum_xmisfit(3) + xmisfit(n,8)
        sum_xmisfit(4) = sum_xmisfit(4) + xmisfit(n,10)
        sum_xmisfit(5) = sum_xmisfit(5) + xmisfit(n,11)
        sum_xmisfit(6) = sum_xmisfit(6) + xmisfit_tot(n)
      enddo

      write (14,*) 'Misfit AHe:',sum_xmisfit(1)
      write (14,*) 'Misfit AFt:',sum_xmisfit(2)
      write (14,*) 'Misfit ZHe:',sum_xmisfit(3)
      write (14,*) 'Misfit ZFt:',sum_xmisfit(4)
      write (14,*) 'Misfit MAr:',sum_xmisfit(5)
      write (14,*) 'Total Misfit:',sum_xmisfit(6)

      endif
        close (14)
      endif

! terminates the job
      close (10)
      deallocate (fta,fte,ftag,fteg)
      deallocate (ftm,ftmg,ftz,ftzg)
      deallocate (age_obs,dage_obs,xmisfit,age_prd)
      deallocate (age_flags)
      deallocate (iconsurf)
      deallocate (xsurf,ysurf,zsurf)
      deallocate (proc,eproc,ice)
      deallocate (header_info,therm_his_val)
      deallocate (zsurf_all)
      deallocate (xdepth,ydepth,zdepth,tprev,zsurfp)

      !deallocate (eheat,econd,edot,bg_edot,topo_edot)
      !deallocate (edot_store)
      !deallocate (jrec,x_basin,y_basin,data_file)
      !deallocate (age_type,fnme)

      call cpu_time (times(9))

      if (iproc.eq.0) write (6,*) 'Total times : ',times(9)-times(1)

      close (6)

      end
