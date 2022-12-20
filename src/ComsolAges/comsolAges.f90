program comsolAges

! Written by Willi Kappler, 2013.05
! willi.kappler@uni-tuebingen.de
!
! Using age calculation from Pecube and Terra
!

  use m_file_io
  use m_ageCalculation

  implicit none

  integer(4), parameter :: fileID = 30

  integer(4) :: commandArgCount, numParticles, numTimeSteps
  integer(4) :: i, j, ioStatus, timeIndex, diffIndex

  integer(4), dimension(:), allocatable :: particleID

  character(1024) :: inputFileName
  character(1024) :: outputFileName

  logical, dimension(:), allocatable :: valid_surface_particle

  real(8), dimension(:), allocatable :: time, zmin, zmax, time_rev
  real(8), dimension(:,:), allocatable :: temperature, all_ages
  real(8), dimension(:,:), allocatable :: xValue, yValue, zValue

  real(8) :: apatiteHe, zirconHe
  real(8) :: zirconFT, apatiteFT
  real(8) :: muscovite_age

  real(8), parameter :: zTop = 99.990 ! [km]
  real(8) :: lastZValue, diffTime

  commandArgCount = command_argument_count()

  print *, "number of command line argumets given: ", commandArgCount

  if (commandArgCount /= 2) then
     print *, "No file name given! Please provide the input and output file name on the command line"
     stop
  end if

  call get_command_argument(1, inputFileName)

  call open_input_file(inputFileName, time, temperature, &
       particleID, xValue, yValue, zValue, &
       zmin, zmax, numParticles, numTimeSteps)

  call get_command_argument(2, outputFileName)


  ! open(fileID, file="matlab_out_before.dat", iostat=ioStatus)

  ! if (ioStatus /= 0) then
  !    print *, "error while opening file: 'matlab_out_before.dat'"
  !    stop
  ! end if

  ! do i = 1, numTimeSteps
  !    write(fileID, *) time(i), temperature(:, numTimeSteps - i + 1)
  ! end do

  ! close(fileID)


!  print *, time

  allocate(valid_surface_particle(numParticles))

  ! get rid of entries that we are not interested in
  ! ( a.k.a. z filter)

  do i = 1, numParticles
     timeIndex = 0
     lastZValue = 0.0
     valid_surface_particle(i) = .false.

     do j = 1, numTimeSteps
        if (zValue(i, j) >= zTop) then
           valid_surface_particle(i) = .true.
           if (lastZValue >= zTop) then
              timeIndex = j -1
              exit
           else
              lastZValue = zValue(i, j)
           end if
        end if
     end do

     if (timeIndex > 0) then
        diffIndex = numTimeSteps - timeIndex
        diffTime = time(numTimeSteps) - time(timeIndex)

!        print *, i, zValue(i, timeIndex + 1), zValue(i, timeIndex), &
!             zValue(i, timeIndex - 1), zValue(i, timeIndex - 2), diffIndex

        do j = numTimeSteps, diffIndex, -1
           temperature(i, j) = temperature(i, j - diffIndex + 1)
           xValue(i, j) = xValue(i, j - diffIndex + 1)
           yValue(i, j) = yValue(i, j - diffIndex + 1)
           zValue(i, j) = zValue(i, j - diffIndex + 1)
        end do

        do j=diffIndex, 1, -1
           temperature(i,j) = temperature(i, 1)
           xValue(i, j) = xValue(i, 1)
           yValue(i, j) = yValue(i, 1)
           zValue(i, j) = zValue(i, 1)
        end do
     end if
  end do

  ! open(fileID, file="matlab_out_after.dat", iostat=ioStatus)

  ! if (ioStatus /= 0) then
  !    print *, "error while opening file: 'matlab_out_after.dat'"
  !    stop
  ! end if

  ! do i = 1, numTimeSteps
  !    write(fileID, "(F32.16)", advance='no') time(i)
  !    do j = 1, numParticles
  !       if (valid_surface_particle(j)) then
  !          write(fileID, "(F32.16)", advance='no') temperature(j, numTimeSteps - i + 1)
  !       end if
  !    end do
  !    write(fileID, *) ""
  ! end do

  ! close(fileID)

  print *, "num of valid particles: ", count(valid_surface_particle)

  allocate(all_ages(numParticles, 5))

! reverse time for age calculation

  allocate(time_rev(numTimeSteps))

  time_rev = time(numTimeSteps:1:-1)

!  print *, time_rev

  apatiteHe = 0.0
  zirconHe = 0.0
  zirconFT = 0.0
  apatiteFT = 0.0
  muscovite_age = 0.0

  print *, "write tecplot output file for the last time step"

  print *, "opening output file: '", trim(outputFileName), "'"

  open(fileID, file=trim(outputFileName), iostat=ioStatus)

  if (ioStatus /= 0) then
     print *, "error while opening file: '", trim(outputFileName), "'"
     stop
  end if

  write(fileID, *) 'TITLE = "ComsolAges"'
  write(fileID, *) 'VARIABLES = "ID", "Time (Ma)", "Temp (C)", "x (km)", "y (km)", "z (km)", &
       & "ZMin (km)", "ZMax (km)", "AHe (Ma)", "AFT (Ma)", "ZHe (Ma)", "ZFT (Ma)", "MAR (Ma)"'
  write(fileID, *) 'ZONE T = "comsolAges", DATAPACKING = POINT,  I = ', count(valid_surface_particle)

! TODO: run loop multi-threaded with OpenMP.
! Store the result in a big array and write data after the loop is finished.

  do i = 1, numParticles
     if (valid_surface_particle(i)) then
        call Mad_He(time_rev, temperature(i, :), numTimeSteps, apatiteHe, zirconHe)
        call correctNumber(apatiteHe)
        call correctNumber(zirconHe)

        call ZFT(time_rev, temperature(i, :), numTimeSteps, zirconFT)
        call correctNumber(zirconFT)

!        call Mad_Trax_AFT(time_rev, temperature(i, :), numTimeSteps, apatiteFT)
        call get_aft_age(time, temperature(i, :), numTimeSteps, apatiteFT)
        call correctNumber(apatiteFT)

        call Muscovite(time_rev, temperature(i, :), numTimeSteps, muscovite_age)
        call correctNumber(muscovite_age)

        write(fileID, *) particleID(i), time(numTimeSteps), temperature(i, numTimeSteps), &
             xValue(i, numTimeSteps), yValue(i, numTimeSteps), zValue(i, numTimeSteps), & 
             zmin(i), zmax(i), apatiteHe, apatiteFT, zirconHe, zirconFT, muscovite_age

        all_ages(i, 1) = apatiteHe
        all_ages(i, 2) = apatiteFT
        all_ages(i, 3) = zirconHe
        all_ages(i, 4) = zirconFT
        all_ages(i, 5) = muscovite_age

     end if

     call printProgress(i, numParticles)

  end do

  close(fileID)

  ! open(fileID, file="matlab_out_ages.dat", iostat=ioStatus)

  ! if (ioStatus /= 0) then
  !    print *, "error while opening file: 'matlab_out_ages.dat'"
  !    stop
  ! end if

  ! do i = 1, 5
  !    do j = 1, numParticles
  !       if (valid_surface_particle(j)) then
  !          write(fileID, "(F32.16)", advance='no') all_ages(j, i)
  !       end if
  !    end do
  !    write(fileID, *) ""
  ! end do

  ! close(fileID)


  deallocate(all_ages)
  deallocate(time_rev)
  deallocate(time)
  deallocate(temperature)
  deallocate(particleID)
  deallocate(xValue)
  deallocate(yValue)
  deallocate(zValue)
  deallocate(zmin)
  deallocate(zmax)

end program comsolAges

subroutine correctNumber(x)
! "fix" NaN and Inf numbers so that they do not appear in the Tecplot
! output file. Tecplot does not like invalid numbers and will refuse to
! read in the data file

  implicit none

  real(8), intent(inout) :: x

  integer(8) :: IPinf, IMinf

  real(8) ::  Pinf, Minf

  data IPinf/B'01111111100000000000000000000000'/    ! +Infinity
  data IMinf/B'11111111100000000000000000000000'/    ! -Infinity

  Pinf = transfer(IPinf,Pinf)
  Minf = transfer(IMinf,Minf)

  if ((x /= x) .or. (x == Pinf) .or. (x == Minf) .or. (x < 0.0)) then
     x = 0.0
  end if
end subroutine correctNumber

subroutine printProgress(i, max_num)
  implicit none

  integer(4), intent(in) :: i, max_num
  integer(4) :: j, five_percent, current_percent
  integer(4) :: command_status

  five_percent = max_num / 20

  do j = 1, 20
     current_percent = (five_percent * j)
     if (i == current_percent) then
        print "(A, I3, A)", "progress: ", 100 * current_percent / max_num, "%"
        command_status = system("date")
     end if
  end do

end subroutine printProgress
