
      subroutine find_upstream_points (fnme,num_topo_files,xoutlet,youtlet,nstep,nx,ny,dx,dy,&
                 run,fteg,ftag,ftzg,ftmg,age_type,nsurf,xdepth,ydepth,edot,&
                 xlonmin,xlatmin,xlonmax,xlatmax,zsurf_all,xmin2,xmax,ymin2,&
                 ymax,num_basins,pages,perates,nrun)

      ! This subroutine is passed in x, y, and z coordinates and creates a flow path network
      ! The closest grid point to the specified basin outlet point is found and used for
      ! finding all points upstream of the basin outlet
      ! Once all points upstream are found, pdfmaker_for_pecube.f90 is called to calculate
      ! and output the PDF for the specified basin 
      !
      ! The UPSTREAMJ and UPSTREAMI arrays hold the x and y indices, respectively, of the points which
      ! flow into the associated point
      ! The surrounding points to a center point appear as follows:
      ! *******************************
      ! ******** D7   D8   D1 *********
      ! ******** D6   C    D2 *********
      ! ******** D5   D4   D3 *********
      ! *******************************
      ! So, the point D7 with center point (C) of (j,i) would be held in the UPSTREAM arrays as
      ! UPSTREAMJ(j,i,7) and UPSTREAMI(j,i,7). If it flows into the center point then the
      ! values of j-1 and i-1 are stored in the arrays, respectively
      ! The point j,i can have at most 8 values in the UPSTREAM arrays because at most 8 points
      ! can surround the center point (ie. (j-1,i-1); (j-1,i); (j+1,i+1); etc)
      ! This creates a flow network that can easily be tracked when given a point on the grid
      ! The subroutine getCatchment recursively finds the upstream points and appends a list
      ! of the x values, y values, ages, and erosion rates

      IMPLICIT NONE

      integer num_topo_files,i,nx,ny,nstep,j,k,m,nsurf,num_basins
      character fnme(num_topo_files)*300,run*100
      real*8 Z(nx,ny),xoutlet(num_basins),youtlet(num_basins),dx,dy,distance
      real DI,DIS,CAREA
      integer UPSTREAMI(nx,ny,8),UPSTREAMJ(nx,ny,8),age_type(num_basins)
      real*8 x(nx,ny),y(nx,ny),min_distance
      integer xmin,ymin
      real*4 fteg(nstep,nsurf,8),edot(nstep,nsurf),ages(nstep,nsurf)
      real*4 ftag(nstep,nsurf),ftzg(nstep,nsurf),ftmg(nstep,nsurf)
      real*8 xdepth(nstep,nsurf),ydepth(nstep,nsurf)
      integer md,ncp,num_xchars,num_ychars
      real*8 x2(nx*ny),y2(nx*ny),new_ages(nx*ny)
      real*8 wk(8*nsurf),ages2(nsurf),new_edot(nx*ny)
      integer iwk(max(31,27+4)*nsurf+nx*ny)
      real*8 pdf_ages(nx*ny),edot2(nsurf),pages(num_basins,nx*ny)
      real*8 new_ages2(nx,ny),new_edot2(nx,ny),pdf_edot(nx*ny),perates(num_basins,nx*ny)
      real*8 ZD,ZDD,CHAN
      integer counter,c,MFP,nrun,ios
      character x_basin_char*100,y_basin_char*100,t4*4
      logical contained(nx,ny),tapesg
      real*8 xlonmin,xlatmin,xlonmax,xlatmax,zsurf_all(nstep,nsurf)
      real*8 xmin2,xmax,ymin2,ymax,zmin
      real*8 xstore(nx*ny),ystore(nx*ny),zstore(nx*ny)

      md=1
      ncp=4
      ! Iterates through all time steps
      do i=1,nstep
        if (i.lt.num_topo_files) then
          open (203,file='input/Pecube-D/'//fnme(i+1),status="old",iostat=ios)              ! Opens the topography file for the time step
          if (ios.ne.0) then
            open (203,file=fnme(i+1),status='old')
          endif
        else
          open (203,file='input/Pecube-D/'//fnme(num_topo_files),status="old",iostat=ios)
          if (ios.ne.0) then
            open (203,file=fnme(num_topo_files),status='old')
          endif
        endif
        read (203,*) (Z(:,k),k=ny,1,-1)                                     ! Reads in all elevation values
        close (203)

        CAREA = 9e20                                            ! Maximum cross-grading area (m^2). Hard coded not sure if value is valid for all models
        DI = dx*1000.                                           ! Grid cell size (m)
        DIS = dx*1000.                                          ! Grid cell size (m)
        ZD = -999.                                              ! Value for missing data point in input file
        ZDD = -999.                                             ! Value to be output for missing data point
        MFP = 1                                                 ! Method of catchment area computations: 1=D8 or Rho8, 2=FD8/FRho8 (Multiple drainage
                                                                ! direction method using slope weighting algorithm), 3=Stream tube method (currently not implemented)
        CHAN = 1.0                                              ! Scaling factor to reduce elevations to absolute elevations in meters (typically 1.0 or 0.01)
        UPSTREAMI(:,:,:)=0                                      ! Initialize upstream y value indices to 0
        UPSTREAMJ(:,:,:)=0                                      ! Initialize upstream x value indices to 0

        ! Subroutine creates flow path network
        ! Stores x and y indices of surrounding points (max of 8) that flow into
        ! center point
        call DEPLSS(ny,nx,nx*ny,CAREA,DI,DIS,ZDD,ZD,&
                    CHAN,MFP,Z(:,:),UPSTREAMI,UPSTREAMJ)

        ! Interpolates the current calculated erosion rate distrbution (with skipping factor)
        ! onto the full resolution mesh (with no skipping factor)
        wk = 0.
        iwk = 0
        new_edot = 0.
        edot2 = edot(i,:)

        call idbvip (md,ncp,nsurf,xdepth(i,:),ydepth(i,:),edot2,nx*ny,&
                     x2,y2,new_edot,iwk,wk)

        x = 0.
        y = 0.
        x2 = 0.
        y2 = 0.
        counter=1
        do k=1,ny
          do j=1,nx
            x(j,ny-k+1) = xlonmin+dx*(nx-1)*float(j-1)/float(nx-1)
            y(j,ny-k+1) = xlatmin+dy*(ny-1)*float(k-1)/float(ny-1)
            !x(j,k)=xlonmin+(xlonmax-xlonmin)*(x(j,k)-xmin2)/(xmax-xmin2)
            !y(j,k)=xlatmin+(xlatmax-xlatmin)*(y(j,k)-ymin2)/(ymax-ymin2)
            x2(counter) = x(j,ny-k+1)
            y2(counter) = y(j,ny-k+1)
            counter = counter + 1
          enddo
        enddo

        do m=1,num_basins

          if (age_type(m).ge.1 .and. age_type(m).le.7) then
            ages = fteg(:,:,age_type(m))
          else if (age_type(m).eq.8) then
            ages = ftag
          else if (age_type(m).eq.9) then
            ages = fteg(:,:,8)
          else if (age_type(m).eq.10) then
            ages = ftzg
          else if (age_type(m).eq.11) then
            ages = ftmg
          endif

          ! Interpolates the current calculated age distribution set (with skipping factor)
          ! onto the full resolution mesh (with no skipping factor)
          wk = 0.
          iwk = 0
          new_ages = 0.
          ages2 = ages(i,:)

          call idbvip (md,ncp,nsurf,xdepth(i,:),ydepth(i,:),ages2,nx*ny,&
                       x2,y2,new_ages,iwk,wk)

          ! Converts the interpolated 1-D arrays of ages and erosion rates to 2-D
          ! to be used when finding the upstream points since the upstream indices
          ! are stored as (j,i) pairs in UPSTREAMJ and UPSTREAMI
          counter=1
          do j=1,ny
            do k=1,nx
              new_ages2(k,j) = new_ages(counter)
              new_edot2(k,j) = new_edot(counter)
              counter = counter + 1
            enddo
          enddo

          ! Calculates the x and y values for the full resolution Pecube model (skipping factor of 1)
          ! The 2-D arrays are used for ease in obtaining the x and y coordinates when finding
          ! upstream points of basin outlet
          ! The 1-D arrays are used for the interpolations (see below) of the ages and erosion rates
          ! onto the full resolution model
          ! Distance is always >= 0 so initialize to -1 to check if it is the first distance calculation
          ! Finds the closest point on the Pecube grid to the basin outlet specified
          ! This closest point is used for the upstream calculations
          min_distance = 100000.
          xmin = 0
          ymin = 0
          distance = 0.
          do k=1,ny
            do j=1,nx
              distance = sqrt ( ( xoutlet(m) - x(j,ny-k+1) )**2. + ( youtlet(m) - y(j,ny-k+1) )**2. )
              if ( ( distance .lt. min_distance ) ) then
                min_distance = distance
                ! The xmin and ymin variables store the indices of the x and y values of the closest point
                xmin = j
                ymin = ny-k+1
              endif
            enddo
          enddo

          ! If the distance between the specified basin outlet and the closest point is more than
          ! the distance between two nodes in the full resolution mesh, than skip that outlet
          if (min_distance .gt. sqrt(dx**2.+dy**2.)) then
            write (6,*)
            write (6,'(A25,f7.2,A1,f7.2,A48)') 'Error: The basin outlet (',xoutlet(m),',',youtlet(m),') is too far from its closest point in the model'
            write (6,'(A29,i1)') 'Skipping outlet for timestep ',i
          else
            ! xoutlet, youtlet, and t4 are used in the getCatchment.f90 subroutine
            ! for opening a new file for output of the x, y, and z positions
            ! of the basin points
            do j=1,100
              x_basin_char(j:j)=' '
              y_basin_char(j:j)=' '
            enddo
            write (x_basin_char,'(f100.2)') xoutlet(m)
            write (y_basin_char,'(f100.2)') youtlet(m)
            do j=1,100
              if (x_basin_char(j:j).ne.'0' .and. x_basin_char(j:j).ne.'1' .and. x_basin_char(j:j).ne.'2' &
                  .and. x_basin_char(j:j).ne.'3' .and. x_basin_char(j:j).ne.'4' .and. x_basin_char(j:j).ne.'5' &
                  .and. x_basin_char(j:j).ne.'6' .and. x_basin_char(j:j).ne.'7' .and. x_basin_char(j:j).ne.'8' &
                  .and. x_basin_char(j:j).ne.'9' .and. x_basin_char(j:j).ne.'.') num_xchars=j+1
              if (y_basin_char(j:j).ne.'0' .and. y_basin_char(j:j).ne.'1' .and. y_basin_char(j:j).ne.'2' &
                  .and. y_basin_char(j:j).ne.'3' .and. y_basin_char(j:j).ne.'4' .and. y_basin_char(j:j).ne.'5' &
                  .and. y_basin_char(j:j).ne.'6' .and. y_basin_char(j:j).ne.'7' .and. y_basin_char(j:j).ne.'8' &
                  .and. y_basin_char(j:j).ne.'9' .and. y_basin_char(j:j).ne.'.') num_ychars=j+1
            enddo

            write (t4,'(i4)') i

            if (i.lt.1000) t4(1:1)='0'
            if (i.lt.100) t4(1:2)='00'
            if (i.lt.10) t4(1:3)='000'

            ! The variable 'c' keeps track of the number of upstream points
            ! This subroutine compiles a list of all points upstream of the
            ! basin outlet and their cooresponding x values, y values, ages
            ! and erosion rates
            c=1
            call getCatchment(UPSTREAMI,UPSTREAMJ,nx,ny,new_ages2,x,y,Z,&
                              xmin,ymin,pdf_ages,new_edot2,pdf_edot,c,contained,&
                              xstore,ystore,zstore)

            open (105,file=run(1:nrun)//'/Timestep_'//t4//'_Basin_X_'//x_basin_char(num_xchars:100)//'_Y_'//y_basin_char(num_ychars:100)//'.dat',status='unknown')
            write (105,*) 'TITLE = "Pecube User Defined Basin"'
            write (105,*) 'VARIABLES = "x (km)" "y (km)" "z (km)"'
            write (105,*) 'ZONE T = "Basin"'
            write (105,*) 'I=',c,', J=1, K=1, ZONETYPE=Ordered'
            write (105,*) 'DATAPACKING=POINT'
            write (105,*) 'DT=(DOUBLE DOUBLE DOUBLE)'

            do j=1,c
              write (105,'(3f15.5)') xstore(j),ystore(j),zstore(j)
            enddo
            close (105)

            ! Subroutine to calculate and output PDF for the basin
            call pdfmaker_for_pecube(pdf_ages,c,age_type(m),i,-1,run,pdf_edot,&
                                     x_basin_char,y_basin_char,num_xchars,num_ychars,nrun)

            pages(m,:)=pdf_ages
            perates(m,:)=pdf_edot

          endif
        enddo
      enddo

      return
      end subroutine
