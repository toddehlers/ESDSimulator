! initialize_nodal_geometry

module m_read_nodal_geometry2
   contains
   subroutine read_nodal_geometry2(configData)
      use rt_param
      use cascade_globals

      implicit none

! subroutine to read in initial nodal information
!
!

      type(config) :: configData

      integer(4), parameter :: FILE_ID = 60
      integer(4) :: i, j, num_of_nodes, pos1, pos2

      real(8) :: htemp, htemp0, landslide_erosion_rate, glacial_erosion_rate
      real(8) :: unused_data, delta_z, fluvial_erosion_rate, diffusion_erosion_rate
      real(8) :: dist
      real(8), dimension(:), allocatable :: delta_p

      character(255) :: text_line

      call random_seed()

      print *, "read mesh: ", configData%meshname
      open(FILE_ID, file="input/IceCascade/" // configData%meshname, status='old')

      do i = 1, 3
         read(FILE_ID, *) ! Skip first three header lines
      enddo

      ! The line looks like this:
      ! n=       49707, e=       98361, et=triangle, f=fepoint

      read(FILE_ID, "(A)") text_line

      pos1 = 3
      pos2 = index(text_line, ',') - 1

      read(text_line(pos1:pos2), *) num_of_nodes

      if (configData%nnode /= num_of_nodes) then
         print *, "Number of nodes in global config corrected (old: ", configData%nnode, ", new: ", num_of_nodes, ")"
         configData%nnode = num_of_nodes
      endif

      print *, "num_of_nodes: ", num_of_nodes

      ! VARIABLES =
      ! "1: x [km]"
      ! "2: y [km]"
      ! "3: z [km]"
      ! "4: node"
      ! "5: Precipitation [m/y]"
      ! "6: Fluvial Erosion Rate [m/y]"
      ! "7: Diffusion Erosion Rate [m/y]"
      ! "8: Landslide Erosion Rate [m/y]"
      ! "9: Total Erosion Rate [m/y]"
      ! "10: Catchment Color"
      ! "11: Catchment Number"
      ! "12: Glacial Erosion Rate [m/y]"
      ! "13: Ice Thickness [m]"
      ! "14: Mass Balance [1/y]"
      ! "15: Total Topography [m]"
      ! "16: Sliding Velocity [m/y]"
      ! "17: Gerode Term [m/y]"
      ! "18: Rock/Alluvium Contact [km]"
      ! "19: Isostatic deflection [m/y]"
      ! "20: Slope [m/km]"
      ! "21: Totalflexiso [m]"
      ! "22: Constriction"
      ! "23: Cumulative Erosion [m]"
      ! "24: Surface Area [km*km]"

      do i = 1, configData%nnode
         read(FILE_ID, *) x(i), &
            y(i), &
            htemp, &
            unused_data, & ! node id
            prec(i), &
            fluvial_erosion_rate, &
            diffusion_erosion_rate, &
            landslide_erosion_rate, &
            unused_data, & ! total erosion
            unused_data, & ! catchment color
            unused_data, & ! catchment number
            glacial_erosion_rate, &
            iceth(i), &
            gbalance(i), &
            tott(i), &
            slide(i), &
            gerode_term(i), &
            htemp0, &
            isodh(i), &
            slope(i), &
            memory(i,4), & ! total flexiso
            strict(i), &
            totalerosion(i), &
            memory(i,7) ! surface area

         !call random_number(delta_z)
         h(i) = (htemp * 1000.0_8) !+ (delta_z * 10.0)

         !call random_number(delta_z)
         h0(i) = (htemp0 * 1000.0_8) !+ (delta_z * 10.0)

         memory(i,2) = fluvial_erosion_rate * global_cascade_dt
         memory(i,3) = diffusion_erosion_rate * global_cascade_dt
         memory(i,8) = landslide_erosion_rate * global_cascade_dt
         dhg(i) = glacial_erosion_rate * global_cascade_dt

         memory(i, 5) = 0.0_8

         if (configData%addshelf) then
            shelf(i) = 1.0_8
         else
            shelf(i) = 0.0_8
         endif


         ! print *, "x: ", x(i), ", y: ", y(i), ", h: ", h(i)
      enddo

      close(FILE_ID)

      allocate(delta_p(configData%nnode))

      delta_p = huge(unused_data)

      do i = 1, configData%nnode
         do j = 1, configData%nnode
            if (i /= j) then
               dist = hypot(x(i) - x(j), y(i) - y(j))
               if (dist > 0.0 .and. dist < delta_p(i)) then
                  delta_p(i) = dist
               endif
            endif
         enddo
      enddo

      delta = sum(delta_p) / dble(configData%nnode)

      deallocate(delta_p)

      print *, "New delta (average node spacing): ", delta

      ! delta = configData%sidex / dble(configData%nx - 1)
      ! configData%surfscale = configData%sidex * configData%sidey / dble(configData%nnode)
      ! configData%surfmin = configData%surfscale / 4.0_8
   end subroutine read_nodal_geometry2
end module m_read_nodal_geometry2

