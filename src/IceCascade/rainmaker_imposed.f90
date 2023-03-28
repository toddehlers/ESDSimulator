!-----------------------------------------------------------------
module m_rainmaker_imposed
    use rt_param
    use cascade_globals

    implicit none
    contains
!##################################################################!
!-----------------------------------------------------------------
! Option for Imposed Precipitation      
! This module was added by Victoria M Buford Parks January 2020    
!-----------------------------------------------------------------
    subroutine rainmaker_imposed(configData)
    
! code to import external precipitation data apply to  CASCADE...
        ! only executes if using imposed precipitation (icecascade.in : A9=3)
        ! only called on initialization of cascade and when  
        !   imposed_oro_time counter is greater than configData%imposedPrecipUpdateTime  (value in cascade.f90)
        ! Defined in cascade.f90: oro_time = oro_time+global_cascade_dt

!     inputs: nnode             no. of nodes
!              x,y,h             x,y, and h of nodes
!              surface           surface area assoc. with each node
! outputs:
!              prec              precipitation rate at each node
!              water             precip *area at each node
!              x_gr,y_gr,prec_gr x,y,precip on reg. grid
        ! Notes on: Input precipitation grids
        ! The .txt referenced in icecascade.in / A10 must be 
        ! in this format:
            ! first line: depth in y dimension in km
            ! second line: number of nodes in y dimension
            ! third line: number of samples to take for interpolation
            ! fourth line: step size in x direction for interpolation
            ! fifth line to end of file: time in Myrs and filename
            !
            ! example:
            ! 10.0
            ! 5
            ! 10
            ! 1.0
            ! 8.10 file1.dat
            ! 3.2 file2.dat
            ! 0 file3.dat
            !
            ! The first four lines will be ignored herein
            ! The ages are the end age at which the precipitaiton grid will be applied:
            ! e.g. file1.dat will be applied from the model start time (input C3 in icecascade.in) until 8.1Ma

        ! each grid (e.g. filei.dat) must have a header and 
        ! three input columns: x y precip (units: km, km, and m/yr)
            ! example:
            ! x[km]     y[km]   P[m/yr]
            ! 0.000     0.000   1.200
            ! 0.076     0.001   0.250
            ! 15.000    22.250  4.600

        use rt_param
        use cascade_globals
        use m_check_var

        implicit none
! begin routine.
        type(config) :: configdata
        real(8), dimension(:), allocatable :: x_precip, y_precip, z_precip, &
            InputPrecipTimeStep, PTimeStep,x_grid,y_grid
        real(8) :: resolution
        character(255), dimension(:), allocatable :: PrecipGridFiles
        character(255) :: PrecipFileName_List, PFileName_Grid1,FileName_Precip, &
            FileName3,FileName4
        integer(4) ::  FileUnit1 = 92, FileUnit3 = 94, & 
            FileUnit4 = 95, Filelines, NumberOfPrecipFiles
        integer(4) :: i, j, k, jj, count, iflag_precip
        integer :: ios
        logical :: Test

        
        print*,'Module: rainmaker_imposed'

        PrecipFileName_List = trim('input/IceCascade/'//configData%precipname)
        ! determine the number of precip files 
        call LengthOfFile(PrecipFileName_List, FileLines)
        
        NumberofPrecipFiles = FileLines - 4
        FileLines = 0

        ! open list of precip files
        open(FileUnit1, file=PrecipFileName_List, status='old',action='read')
        ! Read first four lines & do nothing
            do i = 1,4
                read (FileUnit1, *)
            end do

            ! Read said input .txt file (list of precip files)
            allocate (InputPrecipTimeStep(NumberofPrecipFiles), &
              PrecipGridFiles(NumberofPrecipFiles),PTimeStep(NumberofPrecipFiles+1))
            ! Store the names of the precip files 
              ! & compute intended time (in cascade reference frame) of each file
              PTimeStep(i)=0.0_8
              print*,'Model Run Time [yrs]:',configData%endtime
              print*,"i=",1," PTime(i)= ",int(PTimeStep(1))
            do i = 1, NumberofPrecipFiles
                ! Read time steps and input file names; 
                read(FileUnit1, *, iostat=ios), InputPrecipTimeStep(i), PrecipGridFiles(i)
                if (ios/=0) exit
					
                ! reverse time line
                ! This only works if precip files start at the same time as the model start time
                !PTimeStep(i) = (InputPrecipTimeStep(1) - InputPrecipTimeStep(i)) * 1e6_8
                PTimeStep(i+1) = configData%endtime - (InputPrecipTimeStep(i) * 1e6_8) ! in years
                print*,"i=",i+1,' PTime(i) [y]: ',int(PTimeStep(i+1))
                print*,'input time (i) [Ma]: ',int(InputPrecipTimeStep(i)),'Grid:',PrecipGridFiles(i)
                call ShortenString(PrecipGridFiles(i), 4)
                !print *, "PrecipGridFiles(i): ", PrecipGridFiles(i)
            end do
        close(FileUnit1)

        resolution=4 ! points per km

        ! Create Precip Files in initialization
        if (cascade_istep.eq.1) then
            print*,'Initializing rainmaker_imposed'
            call Write_CascadePrecipFiles(PrecipGridFiles, NumberofPrecipFiles, &
                int(configData%sidex), int(configData%sidey), &
                configData%y_flip, PTimeStep,resolution)
            ! then read in first file and interpolate onto precip nodes
            FileName_Precip = 'input/IceCascade/'//trim(PrecipGridFiles(1))//'_interpolated.dat'
            ! go to 100 
        endif
        ! ! Determine which precip file to use, and interpolate it onto nodes.
         do i=2,NumberofPrecipFiles+1
            if     ( (time.ge.0) .and. (time.le.PTimeStep(2)) ) then
                FileName_Precip = 'input/IceCascade/'//trim(PrecipGridFiles(1))//'_interpolated.dat'
                jj=i-1
            elseif ( (time.gt.PTimeStep(i-1)) .and. (time.le.PTimeStep(i)) ) then
                FileName_Precip = 'input/IceCascade/'//trim(PrecipGridFiles(i-1))//'_interpolated.dat'
                jj=i-1
            endif
        enddo
        print *, "FileName_Precip: ", FileName_Precip

        call LengthOfFile(FileName_Precip,FileLines)

                    ! ! Trying to call and create precip files in case they get deleted mid-run;  
                    ! ! error rank mismatch somewhere.
                    ! inquire(file=FileName_Precip, exist=Test)
                    !     if (Test) then
                    !         !FileName(1:(FileName_StringLength - ReduceBy))
                    !         call ShortenString(FileName_Precip,17)
                    !         FileName_Precip=FileName_Precip//'.dat'
                    !         call Write_CascadePrecipFiles(FileName_Precip, int(1), &
                    !             int(configData%sidex), int(configData%sidey), &
                    !             configData%y_flip, PTimeStep,resolution)
                    !         call ShortenString(FileName_Precip,4)
                    !         FileName_Precip=FileName_Precip//'_interpolated.dat'
                    !     end if

        allocate(x_grid(FileLines-1), y_grid(FileLines-1), z_precip(FileLines-1))
        call Read_CascadePrecipFile(FileName_Precip,FileLines-1,x_grid,y_grid,z_precip)
                    ! Interpolate nice precip data onto CASCADE nodes
                    !-----------------------------------------   

                    ! !-----------------------------------------      
                    ! !-----------------------------------------      
                    ! !-----------------------------------------      
                    ! ! write a test precip file !! For testing only!
                    ! FileName3='input/IceCascade/test_input_interpolated.dat'
                    ! do j  = 1,configdata%sidey*resolution
                    !     y_precip(j) = y_grid(j)
                    !     !count=count+configdata%sidex*resolution
                    ! enddo! fix this above to read in y_precip data correctly

                    ! count=1
                    ! do j  = 1,configdata%sidex*resolution
                    !     x_precip(j) = x_grid(j+count)
                    !     count=count+configdata%sidey*resolution
                    ! enddo! fix this above to read in y_precip data correctlyx_precip=x_precip(1,configdata%sidex*resolution)
                    ! open(FileUnit3, file=FileName3, status='new')
                    ! do j=1,size(y_precip)
                    ! do kk=1,size(x_precip)
                    !     prec_gr(kk,j)=z_precip( (j-1)*(configdata%sidex*resolution) + kk)
                    !     write(FileUnit3, *) x_precip(kk), y_precip(j), prec_gr(kk,j)
                    ! end do
                    ! !count=count+configdata%sidex*resolution
                    ! end do
                    ! close(FileUnit3)
                    ! !! For testing only!
                    ! !-----------------------------------------      
                    ! !-----------------------------------------      
                    ! !-----------------------------------------       
        call interp2(configData%nnode,x_grid,y_grid,z_precip,resolution)
                    ! !-----------------------------------------      
                    ! !-----------------------------------------      
                    ! !-----------------------------------------    
                    ! ! write a test output precip file !! For testing only!
                    ! FileName4='input/IceCascade/test_output_interpolated.dat'
                    ! open(FileUnit4, file=FileName4, status='new')
                    ! ! count=0
                    !     do j=1,size(prec)
                    !         write(FileUnit4, *)  x(j), y(j), prec(j)
                    !     end do
                    ! close(FileUnit4)
                    ! !! For testing only!
                    ! !-----------------------------------------      
                    ! !-----------------------------------------      
                    ! !-----------------------------------------    

    deallocate(InputPrecipTimeStep, PrecipGridFiles, PTimeStep, z_precip, & 
    x_grid, y_grid) 

    ! multiply precip rate by nodal area to get rate of water input.
    !        surface_area     = surface associated with each node, [km2] = memory(1,7)

    do i =1,configData%nnode
        water(i) = memory(i, 7)*prec(i)
    enddo
    !z=zint

return 
end subroutine rainmaker_imposed

! memory(*,1)=dhcrit
! memory(*,2)=dhfluvial     (amount of material eroded/deposited over the time step)
! memory(*,3)=dhdiff        (total height removed/added)
! memory(*,4)=hiso
! memory(*,5)=fix           (Boundary condition array)
! memory(*,6)=newsurface    (which nodal surface area has to be updated)
! memory(*,7)=surface_area  [km2]
! memory(*,8)=dhlandslide   (added by Ehlers 6/01)
! memory(*,9)=dhglacier     (added by Yanites 10/09)

!        x,y         = x- and y-nodal positions [km]
!       h           = present topography [m]
!       dt          = time step length [yr]
!       water       = water discharge at each point
!       prec_gr [m/yr]
!       bckgr [m/yr]

!##################################################################!

!-----------------------------------------
subroutine Write_CascadePrecipFiles(PrecipGridFiles,NumberofPrecipFiles,size_x,size_y,flip,PTimeStep,resolution)
!-----------------------------------------
! added by vmbp Jan 13 2020     
    implicit none

    real(8), dimension(:), allocatable :: x_precip, y_precip, z_precip, &
        x_grid,y_grid,PTimeStep
    real(8), dimension(:,:), allocatable :: z
    real(8) :: resolution !, dimension(:)
    !real(8) :: prec_gr
    character(255), dimension(:) :: PrecipGridFiles
    character(255) ::  PFileName_Grid1, FileName_Precip
    integer(4) :: FileUnit5 = 74, Filelines, NumberOfPrecipFiles,size_x,size_y,i,j,k,kk
    integer :: ios
    logical :: Test,flip

    ! Create Precipitation Files
    do i = 1, NumberOfPrecipFiles
        ! create new Precip file
        FileName_Precip = trim('input/IceCascade/'//PrecipGridFiles(i))//'_interpolated.dat'

        ! Test whether file already exists
        inquire(file=FileName_Precip, exist=Test)
        if (Test) then
            cycle
        end if

        ! Initialize Precip
        ! Read in Precipitation Grid
        PFileName_Grid1 = 'input/IceCascade/'//trim(PrecipGridFiles(i))//'.dat'
        call LengthOfFile(PFileName_Grid1, FileLines)
        print*, "PFileName_Grid1:", PFileName_Grid1, 'FileLines= ', FileLines
        ! Each input precip data grid has A HEADER!
        allocate (x_precip(FileLines-1), y_precip(FileLines-1), z_precip(FileLines-1))
            
        call Read_InputPrecipGrid(PFileName_Grid1, FileLines, y_precip, x_precip, z_precip)
        
        ! Performing y-flip switch if enabled
        if (flip) then
            do j=1,FileLines-1
                x_precip(j) = size_y - x_precip(j)
            end do
        end if

        print*, 'Writing Precipitation file for time step: '
        print *, "i: ", i, "PtimeStep(i): ",nint(PTimeStep(i)),"years PrecipGridFiles(i): ", PrecipGridFiles(i)
                
        ! Create Regular Precipitation grid interpolated from import
        open(FileUnit5, file = FileName_Precip, status='new')
        print*, FileName_Precip
        write(FileUnit5, *), 'REGULAR GRID # x,y,precip,; units: km, mm/year'
        ! define regular grid.
        ! with 1/4km grid spacing 
        ! configdata%sidex = length of model side
        allocate (x_grid(size_x*resolution),y_grid(size_y*resolution))
        do k  = 1,size_x*resolution
                x_grid(k) = dble((k-1)/resolution)
        enddo
        do j  = 1,size_y*resolution
                y_grid(j) = dble((j-1))/resolution
        enddo
        print*, 'size(x_grid):', size(x_grid), 'size(y_grid):', size(y_grid)
        ! modify x_gr and y_gr to appropriate size
        print *, "x_precip size: ", size(x_precip), ", y_precip size: ", size(y_precip), &
        ", z_precip size: ", size(z_precip), ", size_x: ",size_x, ", size_y: ",size_y
        ! calculate/interpolate my precipitation on regular grid
            ! subroutine grid2(nnode,x,y,h,x_gr,y_gr,nxs,nys,z)
        allocate(z(size(x_grid), size(y_grid)))
        call grid2_notglobalvar(size(z_precip(:)), x_precip, y_precip,  z_precip, &
            x_grid,y_grid,size(x_grid), size(y_grid),z) 
        do k=1,size(x_grid)
            do kk=1, size(y_grid)
            write(FileUnit5,*), x_grid(k), y_grid(kk), z(k,kk)
            enddo
        enddo
        close(FileUnit5)

        deallocate (x_precip, y_precip, z_precip,x_grid,y_grid, z)
        
    end do
    print*, 'DONE WRITING IMPOSED PRECIP FILES!'
          
end subroutine Write_CascadePrecipFiles
!-----------------------------------------
! Subroutine that reads in Input Precip grids
!-----------------------------------------
! added by vmbp Jan 13 2020    
subroutine Read_InputPrecipGrid(FileName, Lines, x, y, z)
    implicit none

    character(255), intent(in) :: FileName
    integer(4) :: FileUnit11 = 12, ios, i
    integer(4), intent(in) :: Lines
    real(8), dimension(:), intent(out) :: x, y, z

    open(FileUnit11, file=FileName, status='old')
        read(FileUnit11,*) ! Skips 1st line: header
        do i = 1,Lines
                read(FileUnit11, *, iostat=ios), x(i), y(i), z(i)
        end do
    close(FileUnit11)
end subroutine Read_InputPrecipGrid

!-----------------------------------------------------------------------
! Subroutine to interpolate regularly gridded data onto Cascade mesh
!-----------------------------------------------------------------------
! added by vmbp Jan 13 2020    
subroutine interp2(nnode2,x_grid,y_grid,prec_grid,resolution)!,resolution)
! written by Victoria M. Buford Parks (vmbparks@gmail.com)
! for use with imposed orography, Jan 2020
    use cascade_globals
    use m_check_var

    implicit none
    ! subroutine to do linear interpolation from regular grid back onto nodes
    !     inputs: 
    ! (global)   nnode             number of nodes
    ! (global)   x,y               coords at nodes
    !            x_grid,y_grid         coordinates of regular grid
    ! (global)   prec_grid           precipitation on regular grid
    !            resolution        spacing of x_grid, y_grid values
    !
    !     outputs:
    !            prec              precipitation on nodes
    real(8), dimension(:), intent(in) ::  x_grid, y_grid !resolution
    real(8), dimension(:), intent(in) :: prec_grid(size(x_grid))
    !real(8), intent (out) :: prec(nnode)
    integer(4)  :: i,k, nnode2,&
        idx_xMax_yMax,idx_xMin_yMax,idx_xMax_yMin,idx_xMin_yMin
    real(8)     :: xqmax, xqmin, yqmax, yqmin, & 
                    xmx_ymx_wt, xmx_ymn_wt, xmn_ymx_wt, &
                    xmn_ymn_wt, resolution
    !resolution=abs(x_grid(2)-x_grid(1))
    ! begin loop over nodes
    do i=1, nnode2
        ! determine all nearest neighbors
        xqmax    = ceiling(x(i)*resolution)/resolution
        xqmin    = floor(x(i)*resolution)/resolution
        yqmax    = ceiling(y(i)*resolution)/resolution
        yqmin    = floor(y(i)*resolution)/resolution
        ! remove issues with nodes near the boundary of modle
        if (xqmax.lt.0) then 
            xqmax=0
        elseif (xqmin.lt.0) then
            xqmin=0
        elseif (yqmax.lt.0 )then
            yqmax=0
        elseif (yqmin.lt.0 )then
            yqmin=0    
        elseif (xqmax.gt.maxval(x_grid)) then
            xqmax=maxval(x_grid)
        elseif (xqmin.gt.maxval(x_grid)) then
            xqmin=maxval(x_grid)
        elseif (yqmin.gt.maxval(y_grid)) then
            yqmin=maxval(y_grid)
        elseif (yqmax.gt.maxval(y_grid)) then
            yqmax=maxval(y_grid)    
        end if
        ! Find nearest 4 grid locations to sample point
        do k=1,size(x_grid)
            if ((x_grid(k).eq.xqmax) .and. (y_grid(k).eq.yqmax)) then
                idx_xMax_yMax = k
                !cycle
            elseif ((x_grid(k).eq.xqmin) .and. (y_grid(k).eq.yqmax)) then
                idx_xMin_yMax = k
                !cycle
            elseif ((x_grid(k).eq.xqmax) .and. (y_grid(k).eq.yqmin)) then
                idx_xMax_yMin = k
                !cycle
            elseif ((x_grid(k).eq.xqmin) .and. (y_grid(k).eq.yqmin)) then
                idx_xMin_yMin = k
                !cycle
            endif
         end do
         ! intend to average the above points, but there is some error in weighting. thus, use only one.
        ! do k=1,size(x_grid)
        !     ! Find MOVE velocities within resolution range of regular grid point
        !     allocate(LocArray(size(x_grid)))
        !     do l=1,FileLines
        !         if ((Grid_y1(l) <= (RegGrid_X(k) + RegGrid_Resolution)) .and. &
        !             (Grid_y1(l) > (RegGrid_X(k) - RegGrid_Resolution)) .and. &
        !             (Grid_z1(l) <= (RegGrid_Z(k) + RegGrid_Resolution)) .and. &
        !             (Grid_z1(l) > (RegGrid_Z(k) - RegGrid_Resolution))) then
        !             LocArray(l) = l
        !         else
        !             LocArray(l) = 0.0_8
        !         end if
        !     end do
        !     LocArray = pack(LocArray, LocArray /= 0.0_8 .and. LocArray < FileLines)

        !     ! Calculate average vertical and horizontal velocities of found record.
        !     ! If record is empty, set velocities to zero.
        !     if (size(LocArray) == 0) then
        !         RegGrid_Xvel(k) = 0.0_8
        !         RegGrid_Zvel(k) = 0.0_8
        !         deallocate(LocArray)
        !         allocate(LocArray(1))
        !     else
        !         RegGrid_Xvel(k) = sum(Vel_y(LocArray)) / size(LocArray)
        !         RegGrid_Zvel(k) = sum(Vel_z(LocArray)) / size(LocArray)
        !     end if

        !     ! Write output to file
        !     write(FileUnit2,*), k,  0.0_8, RegGrid_X(k), RegGrid_Z(k), 0.0_8, RegGrid_Xvel(k), RegGrid_Zvel(k), 7
        !     deallocate(LocArray)
        ! end do
        
        ! compute linearly interpolated, weighted precipitaiton
        ! prec(i)  =  (prec_grid(idx_xMax_yMax) * xmx_ymx_wt + &
        !             prec_grid(idx_xMax_yMin)  * xmx_ymn_wt + &
        !             prec_grid(idx_xMin_yMax)  * xmn_ymx_wt + &
        !             prec_grid(idx_xMin_yMin)  * xmn_ymn_wt) /&
        !             (xmx_ymx_wt+xmx_ymn_wt+xmn_ymx_wt+xmn_ymn_wt)
         prec(i)  =  (prec_grid(idx_xMax_yMax))! + prec_grid(idx_xMax_yMin) + &
                    ! prec_grid(idx_xMin_yMax) + prec_grid(idx_xMin_yMin))/4

         if (i.eq.2) then
            print*,"prec(i=2)", prec(i)
            !exit
        end if
    enddo
end subroutine interp2   
    
!-----------------------------------------     
! Subroutine Read_CascadePrecipFile   
!-----------------------------------------       
! added by vmbp Jan 13 2020    
subroutine Read_CascadePrecipFile(FileName, Lines, loc_x, loc_y, loc_z)
! Reads an imposed precipitation file file
    implicit none

    character(255), intent(in) :: FileName
    integer(4) :: FileUnit6 = 13, ios, j
    integer, intent(in) :: Lines
    real(8), dimension(:), intent(out) :: loc_x, loc_y, loc_z
    
    open(FileUnit6, file=FileName, status='old')
    read(FileUnit6,*) ! read header and do nothing
        do j = 1,Lines-1
            read(FileUnit6, *, iostat=ios), loc_x(j),  loc_y(j), loc_z(j)
                
            if (ios/=0) then
                exit
            end if
        end do

    close(FileUnit6)
end subroutine Read_CascadePrecipFile
!-----------------------------------------     
! Subroutine grid2_notglobalvar   
!-----------------------------------------     
! borrowed from rainmaker.f90 
subroutine grid2_notglobalvar(nnode2,x_data,y_data,h_data,x_gr,y_gr,nxs2,nys2,z)
    implicit none
!     inputs: 
!             nnode2     number of nodes 
!             x_data,y_data        co-ords of nodes 
!             h_data          elevation of nodes 
!             x_gr, y_gr co-ords of reg grid 
!             nxs, nys   number of values in x_data and y_data direction
!             subsam    subsampling of input nodes.
!
!      outputs: 
!             z       elevation on regular grid 
    ! FIX: initialization /create precip files vs after precip files exist + time comparison

    integer(4), intent(in) :: nnode2, nxs2, nys2
    real(8), intent(in) :: x_data(nnode2), y_data(nnode2), h_data(nnode2), x_gr(nxs2), y_gr(nys2)
    real(8), intent(out) :: z(nxs2, nys2)

    integer(4) :: i, j, k, n1, n2, n3, n4
    real(8) :: radius, d, d1, d2, d3, d4
    real(8) :: factor1, factor2, factor3, factor4, factor_sum

    ! consider a radius of 10 km
    radius = 10.0_8
    n1 = 1
    n2 = 1
    n3 = 1
    n4 = 1

    print *, "x_data min: ", minval(x_data), &
    ", x_data max_data: ", maxval(x_data), &
    ", y_data min: ", minval(y_data), &
    ", y_data max_data: ", maxval(y_data)
    print *, "x_gr min: ", minval(x_gr), ", x_gr max_data: ", maxval(x_gr)
    print *, "y_gr min: ", minval(y_gr), ", y_gr max_data: ", maxval(y_gr)

    do i=1,nxs2
        do j=1,nys2

            ! now look for th_datae 4 nearest neigh_databors
            d1 = huge(d1)
            d2 = huge(d2)
            d3 = huge(d3)
            d4 = huge(d4)

            do k=1,nnode2
                d = sqrt((x_gr(i) - x_data(k))**2 + (y_gr(j) - y_data(k))**2)

                if (d < d1) then
                    d2 = d1
                    d3 = d2
                    d4 = d3
                    d1 = d

                    n2 = n1
                    n3 = n2
                    n4 = n3
                    n1 = k
                else if (d < d2) then
                    d3 = d2
                    d4 = d3
                    d2 = d

                    n3 = n2
                    n4 = n3
                    n2 = k
                else if (d < d3) then
                    d4 = d3
                    d3 = d

                    n4 = n3
                    n3 = k
                else if (d < d4) then
                    d4 = d

                    n4 = k
                endif
            enddo

            factor1 = exp(-d1/radius)
            factor2 = exp(-d2/radius)
            factor3 = exp(-d3/radius)
            factor4 = exp(-d4/radius)

            factor_sum = factor1 + factor2 + factor3 + factor4
            z(i, j) = ((h_data(n1)*factor1) + (h_data(n2)*factor2) + (h_data(n3)*factor3) + (h_data(n4)*factor4)) / factor_sum

        enddo
    enddo
end subroutine grid2_notglobalvar

!-----------------------------------------
! Subroutine that shortens a given string by an integer number of
!-----------------------------------------
! added by vmbp Jan 13 2020     
subroutine ShortenString(FileName, ReduceBy)
    implicit none

    character(255), intent(out) :: FileName
    integer(4), intent(in) :: ReduceBy
    integer(4) :: FileName_StringLength

    FileName_StringLength = len_trim(FileName)

    if (FileName_StringLength < ReduceBy) then
        print*, 'Subroutine ShortenString: String is too short to be shortened by ', ReduceBy, ' characters.'
    else
        FileName = FileName(1:(FileName_StringLength - ReduceBy))
    end if
end subroutine ShortenString
!-----------------------------------------
! Subroutine that determines the length of a file
!-----------------------------------------
! added by vmbp Jan 13 2020
subroutine LengthOfFile(FileName, LineCounter)
    implicit none

    character(255), intent(in) :: FileName
    integer(4) :: FileUnit2 = 11, ios
    integer(4), intent(out) :: LineCounter

    open(FileUnit2, file=FileName, status='old')
        LineCounter = 0
        do
            read(FileUnit2, *, iostat=ios)
            if (ios/=0) then
                exit
            end if
            LineCounter = LineCounter + 1
        end do
    close(FileUnit2)
end subroutine LengthOfFile  
end module m_rainmaker_imposed