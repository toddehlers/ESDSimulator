        module m_ice
            contains
            subroutine ice(configData)
            use rt_param
            use cascade_globals
            use m_total_topography
            use m_write_int
            use m_time_step
            use m_basal_temperature
            use m_calve
            use m_mb_ice
            use m_avalanche
            use m_update_height
            use m_update_height_comsol
            use m_interp_ice
            use m_write_output
            use m_export_to_comsol
            use m_import_from_comsol

            implicit none

            type(config) :: configData

            real(8), dimension(:,:), allocatable :: temps,tempb,tempm,junk
            real(8), dimension(:,:), allocatable :: dtemp
            real(8), dimension(:,:), allocatable :: Pe
            real(8), dimension(:,:), allocatable :: h_tmp, h_tmp_comsol ! [km]
            real(8), dimension(:,:), allocatable :: constriction
            real(8), dimension(:,:), allocatable :: t, tsmooth, ht, ice_mass_balance, d
            real(8), dimension(:,:), allocatable :: dh, sliding, sliding_comsol

            logical :: verbose_ice
            logical :: verbose_eros

            integer :: i,j,k,nx,ny, intime, ice_istep, com_run

            real(8) :: chmax, dhmax, dtice, dticep, dx, dy, lastdh
            real(8) :: ice_time_tmp
            real(8) :: timeice

!
            verbose_ice = .true.
            verbose_eros = .true.

            nx=nxice
            ny=nyice

            print *, "ice.f90, nx: ", nx
            print *, "ice.f90, ny: ", ny

            allocate(temps(nx,ny))
            allocate(tempb(nx,ny))
            allocate(tempm(nx,ny))
            allocate(junk(nx,ny))
            allocate(dtemp(nx,ny))
            allocate(Pe(nx,ny))
            allocate(constriction(nx,ny))
            allocate(t(nx,ny)) ! [km]
            allocate(ht(nx,ny)) ! [km]
            allocate(tsmooth(nx,ny))
            allocate(ice_mass_balance(nx,ny))
            allocate(d(nx,ny))
            allocate(dh(nx,ny))
            allocate(sliding(nx,ny))
            allocate(sliding_comsol(nx,ny))
            allocate(h_tmp(nx,ny))
            allocate(h_tmp_comsol(nx,ny))

! WK: init memory:

        do i=1,nx
            do j=1,ny
                    temps(i,j) = 0.0_8
                    tempb(i,j) = 0.0_8
                    tempm(i,j) = 0.0_8
                    junk(i,j) = 0.0_8
                    dtemp(i,j) = 0.0_8
                    Pe(i,j) = 0.0_8
                    constriction(i,j) = 0.0_8
                    t(i,j) = 0.0_8 ! is in km
                    ht(i,j) = 0.0_8 ! is in km
                    tsmooth(i,j) = 0.0_8
                    ice_mass_balance(i,j) = 0.0_8
                    d(i,j) = 0.0_8
                    dh(i,j) = 0.0_8
                    sliding(i,j) = 0.0_8
                    h_tmp(i,j) = 0.0_8 ! is in km
            enddo
        enddo

            verbose_ice = .FALSE.
            verbose_eros = .FALSE.

            i = 0
            j = 0
            k = 0
            intime = 0

            chmax = 0.0_8
            dhmax = 0.0_8
            dtice = 0.0_8
            dticep = 0.0_8
            dx = 0.0_8
            dy = 0.0_8
            lastdh = 0.0_8 ! does this need to be preserved between function calls ?
            ice_time_tmp = 0.0_8
            timeice = 0.0_8


! Define the spatial discretization:
! There are 2 spatial grids (see manual), nx and ny define the coarse grid, nx2,ny2 the high resolution grid.
! There is interpolation from one grid to another. 

! It also important the low resolution grid is a power of 2 (for the flexure calculation)

    
        
      ice_istep=0
!    tface=tzpt
! WK, 2011.07.13: unit conversion from km to m
! delx, dely in km
! dx, dy in m
      dx=delx*1000.0_8
      dy=dely*1000.0_8
! dx1 = delx
! dy1 = dely


      dh=0.0_8

! WK, 2011.07.13: tface == t == z -> regular grid [km].
! The grid is transformed from the iregular cascade grid (del_gr)
! to a regular grid (z) inside rainmaker, subroutine "grid".

    print *, "ICE.f90, z: ", minval(z), maxval(z)

    t=z ! t(nx=nxice=nxs, ny=nyice=nys), tface(nxi=200, nyi=120)
    dtice=configData%shallow_ice
       do i=1+1,nx-1
          do j=1+1,ny-1
          tsmooth(i,j)=(8.0_8*t(i,j)+2.0_8*t(i-1,j)+2.0_8*t(i+1,j)+2.0_8*t(i,j-1)+2.0_8*t(i,j+1)+t(i-1,j-1)+t(i-1,j+1)+&
         t(i+1,j-1)+t(i+1,j+1))/20.0_8
         enddo
      enddo

     do k=1,3
        do i=1+1,nx-1
              do j=1+1,ny-1
              tsmooth(i,j)=(8.0_8*tsmooth(i,j)+2.0_8*tsmooth(i-1,j)+2.0_8*tsmooth(i+1,j)+2.0_8*tsmooth(i,j-1)+&
              2.0_8*tsmooth(i,j+1)+tsmooth(i-1,j-1)+tsmooth(i-1,j+1)+&
             tsmooth(i+1,j-1)+tsmooth(i+1,j+1))/20.0_8
             enddo
          enddo
        enddo

! WK, 2011.07.14, output prec_gr

!    open(50, file="prec_gr.dat")
!    do i=1,ny
!        do j=1,nx
!            write (50, "(F16.8) ", advance="no") prec_gr(j,i)
!        enddo
!        write (50, "(A)") ""
!    enddo
!    close(50)


! OPEN output files
!       open (734,file=run//'/topography.txt',status='unknown')
!       write (734,*)dx2,dy2
!       open (735,file=run//'/stats.txt',status='unknown')

! beginning of time stepping

!avalanche_counter = 0

ice_time_tmp=0.0_8
ice_istep=0
intime=1
!tt_tmp=0.0_8
h_tmp=hforice ! in km
      
!ice_thickness_counter = 1

print *, "ICE.f90: temperature: ", temperature

11111 continue
!      tp=t
      ice_time_tmp=ice_time_tmp+configData%dt_icetime
      ice_istep=ice_istep+1
      timeice=0.0_8

!    h_tmp_max = maxval(h_tmp-h_tmp_old)

    print *, "------------------"
    tsys0 = system('date')
    print *, "t (z): ", minval(t), maxval(t)
    if (configData%ice_tfinal > 0.0_8) then
        print *, "ice_time_tmp: ", ice_time_tmp, " (", configData%ice_tfinal, "), ", 100.0_8 * ice_time_tmp / configData%ice_tfinal
    endif
    print *, "ice_istep: ", ice_istep
    print *, "dtice: ", dtice
    print *, "dticep: ", dticep
!    print *, "ice thickness diff: ", h_tmp_max
!    print *, "ice_thickness_counter: ", ice_thickness_counter
!    print *, "avalanche_counter: ", avalanche_counter, avalanche_counter / ice_thickness_counter
    print *, "------------------"

!    h_tmp_old = h_tmp
!    ice_thickness_counter = 1

!************* it first interpolates the high res topo grid on the low res topo grid**********************
 !     tcheck=t2
      if (intime == 1 .or.ice_time_tmp > 1000.0_8) then
        intime=0
!        tcheck=t2
!        call highres2lowres (nx,ny,nx2,ny2,t,t2,dx,dy,dx2,dy2)
!        t2=tcheck
           !here it filters the grid to make it smoother, needs to be done otherwise can be unstable
           do i=1+1,nx-1
              do j=1+1,ny-1
                t(i,j)=(8.0_8*t(i,j)+2.0_8*t(i-1,j)+2.0_8*t(i+1,j)+2.0_8*t(i,j-1)+2.0_8*t(i,j+1)+t(i-1,j-1)+t(i-1,j+1)+&
                t(i+1,j-1)+t(i+1,j+1))/20.0_8
             enddo
          enddo
        endif
        
        ! initialize mass balance

        call total_topography (nx,ny,h_tmp,t,ht)
        
!       print*,maxval(prec_gr),minval(prec_gr)
        call mb_ice(nx,ny,ht,temperature,ice_mass_balance,configData)
! The ice thickness is nil (to do so make sure that the erosion time step is large)
! This will need some thoughts to check whether it is completely valid...
!       print*,maxval(a),minval(a)
!        h=0.
!        h2=0.
        !sliding=0.0_8 ! sliding is initialized as 0, so is this necessary?
!        sliding2=0.
        tempm=0.0_8
        tempb=0.0_8


! HERE IT STARTS TO CALCULATE THE ICE THICKNESS*********************
22222 timeice=timeice+dtice
      dticep=dtice
!      ice_thickness_counter = ice_thickness_counter + 1
      if (timeice > configData%dt_icetime) dtice=configData%dt_icetime-(timeice-dtice)

      
! WK: inline the function to find 'double free' crash


      do i=1,nx
         do j=1,ny
            temps(i,j)=temperature-configData%xlapse_rate*ht(i,j) !-0.7576*xlatitude
            tempb(i,j)=temperature-configData%xlapse_rate*t(i,j) !-0.7576*xlatitude
         enddo
      enddo



! !   now calculating precip in cascade using rainmaker output
!       call precipitation_rate (nx,ny,olda,temps,tempb,tempm,friction_flux,qb,h,xmb,xpmin,xpmax)

      if (configData%itemp) then
       call basal_temperature (nx,ny,tempb,temps,tempm,dtemp, &
                              Pe,junk,ice_mass_balance, &
                              diffusivity,SEC_IN_YR,dtdzb,h_tmp)

! WK: adjust the melting tempertature:
! WK: removed call to external function and inlined it instead

!      call melting_temperature (nx,ny,tempm,h_tmp)

        do i=1,nx
            do j=1,ny
                if (h_tmp(i,j) > 0.0_8) then
                   tempm(i,j) = -8.7e-4_8 * h_tmp(i,j)
                else
                    tempm(i,j) = 0.0_8
                endif
            enddo
        enddo

      endif
    
! WK 2011.07.15: commented out calve to test code
      call calve(nx,ny,h_tmp,t,ht,configData%sea_level,dx,dy,dtice)
      call total_topography (nx,ny,h_tmp,t,ht)



!     call effective_diffusivity (nx,ny,di,ds,dd,d, &
!                                 ht,h,tempb,tempm, &
!                                 dx,dy,n,xns,c,cs,xalpha,t,&
!                                gamma,constriction)

!      call velocity (nx,ny,us,vs,ud,vd,u,v,h,ds,dd,tempb,tempm,ht, &
!                      t,d2tdx2,d2tdxdy,d2tdy2,normx,normy,constriction, &
!                      gamma,sliding,dx,dy,friction_flux,xalpha)
! 'a' is now defined in Cascade.  be careful BJY 101209 to use updated mass balance rules

! WK: added sanity check for dtice, as suggested by BJY

    if (dtice > 0.1_8) then
        dtice = 0.1_8
    endif

!      print *, "------------------"
!      print *, "timeice: ", timeice
!      print *, "ice.f90, maxval(h_tmp) 3: ", maxval(h_tmp)
!      print *, "ice.f90, ice_mass_balance 3: ", maxval(ice_mass_balance), minval(ice_mass_balance)



        call update_height (nx,ny,dh,ht,ice_mass_balance,h_tmp, &
                          dx,dy,dtice,configData%dh_allowed,tempb,tempm, &
                          configData%power_law_exp,configData%exp_sliding_law,&
                          GLOBAL_C,GLOBAL_CS,sliding,&
                          configData%constriction,constriction,t, &
                          tsmooth)
!                          

!      print *, "ice.f90, maxval(h_tmp) 4: ", maxval(h_tmp)
!      print *, "ice.f90, ice_mass_balance 4: ", maxval(ice_mass_balance), minval(ice_mass_balance)
!      print *, "------------------"

      call total_topography (nx,ny,h_tmp,t,ht)

! WK 2011.07.15: commented out avalanche to test code
      call avalanche(ht,dx,dy,h_tmp,configData%dh_allowed,nx,ny,dtice)

!      call apply_bc (nx,ny,h,dh)

      ! CONDITIONS ON H TO ENSURE THAT IT DOES NOT BECOME NEGATIVE
        do j=1,ny
          do i=1,nx
            if (h_tmp(i,j) < 0.0_8) h_tmp(i,j)=0.0_8
          enddo
        enddo

      dh=abs(dh)
      dhmax=maxval(dh)
!      dhmax=sum(dh)/nx/ny

!      print *, "ice.f90, dtice 2: ", dtice
! WK, 2011.07.12: disabled time step adaptation for now to test code
     call time_step (dhmax,configData%dh_allowed,dtice,&
           lastdh,chmax)
!      print *, "ice.f90, dtice 3: ", dtice
     call total_topography (nx,ny,h_tmp,t,ht)
        if((dtice.ge.5.0_8).and.(chmax.lt.0.05_8)) then
         print*,'SS Ice MaxDh %MaxDh,DT,TIMEICE',dhmax,chmax*100.0_8,dtice,timeice
          timeice=configData%dt_icetime
        endif
!*********************** END OF ICE THICKNESS CALCULATIONS******************************************
      if ((timeice < configData%dt_icetime)) goto 22222
      
!*********************** COMSOL/Higher-order Ice Physics ******************************************
! RH 2012.09.28: Added com_run so that if comsol doesn't run for the first timeice step,
! it won't run during any of the following steps because the thickness shouldn't change 
! enough to make comsol run properly if it doesn't initially work
! This takes the full ice domain from above and segments into the COMSOL loop

     com_run = 1  ! initially

        ! if com_run == 0, then comsol did not successfully run the last time, and it skips it for this time step\
        ! if com_run == 1, then comsol runs again
        ! if com_run == 2, then comsol successfully ran, and delh is minimized (?)

! this is called when comsol loops and hasn't minimized the change in ice thickness
     if (configData%use_comsol .AND. maxval(h_tmp) > 10) then
          ! ht is the surface topography & tsmooth is the smoothed total topo! originally using ht
          33333 continue
          h_tmp_comsol = h_tmp
            	
          ! find the basal temperature: eventually, this should go into comsol
          do i=1,nx
               do j=1,ny
                    temps(i,j)=temperature-configData%xlapse_rate*(t(i,j)+h_tmp_comsol(i,j)) !-0.7576*xlatitude
                    tempb(i,j)=temperature-configData%xlapse_rate*t(i,j) !-0.7576*xlatitude
               enddo
          enddo
            
          if (configData%itemp) then
               call basal_temperature (nx,ny,tempb,temps,tempm,dtemp, &
               Pe,junk,ice_mass_balance, &
               diffusivity,SEC_IN_YR,dtdzb,h_tmp)
          endif
            
          call export_to_comsol(nx, ny, configData%del_gr, tsmooth, tempb, ice_mass_balance, h_tmp_comsol, t)
          call update_height_comsol()
          call import_from_comsol(nx, ny, 1000*configData%del_gr, h_tmp_comsol, sliding, com_run)
          
          if (com_run == 1) then 
               h_tmp = h_tmp_comsol
               print*,'Comsol ran once, go on'
               goto 33333
          elseif (com_run == 2) then
               print *, "Nested domain in Comsol has converged."
               h_tmp = h_tmp_comsol
          else 
               print *, "Comsol fail."
          endif
     endif ! use_comsol
!*********************** End of COMSOL Section ******************************************

        
!       call velocity (nx,ny,us,vs,ud,vd,u,v,h,ds,dd,tempb,tempm,ht, &
!                      t,d2tdx2,d2tdxdy,d2tdy2,normx,normy,constriction, &
!                      gamma,sliding,dx,dy,friction_flux,xalpha,n,xns,c,cs)
!      sliding=ds
      dtice=dticep
!*******************      it goes back on the high resolution grid**********************************

    call total_topography (nx,ny,h_tmp,t,ht)
!    call lowres2highres  (nx,ny,nx2,ny2,ht,ht2,h,h2,sliding,&
!                          sliding2,dx,dy,dx2,dy2,us,vs,gamma,t2)

!******************** HERE IT CALCULATES EROSION************************************************

! 
!         icefree=1
!         where (h2 > 0.1) icefree=0
! 
        print *,"constriction: ", maxval(constriction), minval(constriction)
        do j=1,ny
          do i=1,nx
            if (ht(i,j).le.t(i,j)) then
              iceftoc(i,j)=1.0_8
            else
              iceftoc(i,j)=0.0_8
            endif
          enddo
        enddo

!      where(anthislope.eq.1) sliding=0.0_8
    
    
        call write_int (nx,ny,h_tmp, &
                        tempm, &
                        ice_mass_balance,ice_time_tmp,configData%dt_icetime, &
                        t, sliding,dx,dy)
            print*,'Temperature: ', temperature
            print *, "======================================================="

! WK, 2011.07.12: write ice data to tecplot ouput file:

!
! Changed 2011.08.29, Willi Kappler
! This condition is only valid when ICE is not called from Cascade but run alone
! So we set ice_tfinal to zero in the configuration file
      if (ice_time_tmp <= dble(configData%ice_tfinal)) goto 11111

!trimline data
!      call trimlines (ht2,icefree,nx2,ny2,dx2,dy2,xfit)
 
!      print*,xfit

      close(8)
      close(9)

! httoc is global [km], changes are going back to cascade
     htoc=h_tmp 
     httoc=ht ! ht = h_tmp + t
     constoc=constriction
!     antitoc=anthislope
     hforice=h_tmp ! in km
     ! slidetoc is a different size than other things
     slidetoc=sliding ! doesn't run past this part!!!!
            
            deallocate(temps)
            deallocate(tempb)
            deallocate(tempm)
            deallocate(junk)
            deallocate(dtemp)
            deallocate(Pe)
            deallocate(constriction)
            deallocate(t)
            deallocate(ht)
            deallocate(tsmooth)
            deallocate(ice_mass_balance)
            deallocate(d)
            deallocate(dh)
            deallocate(sliding)
            deallocate(h_tmp)

      return
!   deallocate (htoc,httoc, slidingtoc,icefreetoc)
          end subroutine ice
      end module m_ice

