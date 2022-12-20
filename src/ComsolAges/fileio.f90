module m_file_io
contains

  subroutine open_input_file(inputFileName, &
                             time, &
                             temperature, &
                             particleID, &
                             xValue, &
                             yValue, &
                             zValue, &
                             zmin, &
                             zmax, &
                             numParticles, &
                             numTimeSteps)

! Reads in the text file exported from Comsol

    implicit none

! input parameter
    character(1024), intent(in) :: inputFileName

! output parameter
    real(8), dimension(:), allocatable, intent(out) :: time, zmin, zmax
    real(8), dimension(:,:), allocatable, intent(out) :: xValue
    real(8), dimension(:,:), allocatable, intent(out) :: yValue
    real(8), dimension(:,:), allocatable, intent(out) :: zValue
    real(8), dimension(:,:), allocatable, intent(out) :: temperature

    integer(4), intent(out) :: numParticles
    integer(4), intent(out) :: numTimeSteps

    integer(4), dimension(:), allocatable, intent(out) :: particleID

! local variables
    real(8), parameter :: sec_per_year = 60.0 * 60.0 * 24.0 * 365
    real(8), parameter :: kelvin_to_degree = 273.15
    real(8), parameter :: meter_to_km = 1000.0
    real(8), parameter :: mil_years = 1.0e6

    integer(4), parameter :: property_per_particle = 4
    integer(4), parameter :: fileID = 20

    integer(4) :: ioStatus, i, j, last_index
    integer(4) :: particle_id, numExpressions, counter

    character(1) :: col1
    character(20) :: col2
    character(40) :: col3

    character(10), dimension(:), allocatable :: time_values

    real(8) :: time_value, min_temperature, max_temperature
    real(8), dimension(:,:), allocatable :: particle_values

    print *, "open file '", trim(inputFileName), "'"

    open(fileID, file=inputFileName, iostat=ioStatus)

    if (ioStatus /= 0) then
       print *, "error while opening file: '", inputFileName, "'"
       stop
    end if


    ! parse the file header
    do i = 1, 8
       read (fileID, *, iostat=ioStatus) col1, col2, col3

       if (trim(col2) == "Nodes:") then
          read (col3, *) numParticles
          print *, "numParticles: ", numParticles
       else if (trim(col2) == "Expressions:") then
          read (col3, *) numExpressions
          numTimeSteps = numExpressions / property_per_particle
          print *, "numTimeSteps: ", numTimeSteps
       endif
    end do

    allocate(time_values((numExpressions * 4) + 1))
    allocate(time(numTimeSteps))
    allocate(particleID(numParticles))
    allocate(temperature(numParticles, numTimeSteps))
    allocate(xValue(numParticles, numTimeSteps))
    allocate(yValue(numParticles, numTimeSteps))
    allocate(zValue(numParticles, numTimeSteps))
    allocate(zmin(numParticles))
    allocate(zmax(numParticles))

    ! initialize variables
    time = 0.0
    temperature = 0.0
    xValue = 0.0
    yValue = 0.0
    zValue = 0.0
    zmin = 0.0
    zmax = 0.0

    print *, "read in time values"

    read (fileID, *) time_values

    counter = 1
    last_index = -10

    ! parse the first data row that contains the time information
    do i = 1, numExpressions * 4
       if (time_values(i)(1:2) == "t=" ) then
          if (i - last_index == 16) then
             read (time_values(i)(3:), *) time_value
             time(counter) = time_value / (sec_per_year * mil_years)
             counter = counter + 1
             last_index = i
          end if
       endif
    end do

    deallocate(time_values)

    allocate(particle_values(property_per_particle, numTimeSteps))

    min_temperature = 100000.0
    max_temperature = 0.0

    print *, "read in particle data"

    do i = 1, numParticles
       call read_particle_data(fileID, property_per_particle, &
       numTimeSteps, particle_id, particle_values)
       particleID(i) = particle_id
       ! layout of input file: temperature, x, y, z

       do j = 1, numTimeSteps
          temperature(i, j) = particle_values(1, j) - kelvin_to_degree
          xValue(i, j) = particle_values(2, j) / meter_to_km
          yValue(i, j) = particle_values(3, j) / meter_to_km
          zValue(i, j) = particle_values(4, j) / meter_to_km
          if (temperature(i, j) > max_temperature) then
             max_temperature = temperature(i, j)
          end if

          if (temperature(i, j) < min_temperature) then
             min_temperature = temperature(i, j)
          end if

       end do
       zmin(i) = minval(particle_values(4, :)) / meter_to_km
       zmax(i) = maxval(particle_values(4, :)) / meter_to_km
    end do

    close(fileID)

    deallocate(particle_values)

    print *, "min_temperature: ", min_temperature
    print *, "max_temperature: ", max_temperature

  end subroutine open_input_file

  subroutine read_particle_data(fileID, property_per_particle, numTimeSteps, &
    particle_id, particle_values)

    ! read in the actual particle data.
    ! we can't use the fortran built in read() function to read the whole line at once,
    ! since the number of characters in each line is way too high.

    implicit none


    integer(4), intent(in) :: fileID, property_per_particle, numTimeSteps
    integer(4), intent(out) :: particle_id

    real(8), intent(out), dimension(property_per_particle, numTimeSteps) :: particle_values

    character(1) :: singleChar
    character(64) :: singleValue

    integer(4) :: i, j, k, ioStatus, valueCounter, charCounter, columnCounter, intValue

    logical :: isInt

    read(fileID, "(I24)", iostat=ioStatus, advance='no') particle_id

    if (ioStatus /= 0) then
       print *, "error while reading input file: ", ioStatus
       stop
    end if

!    print *, "particle_id: ", particle_id

    valueCounter = 1
    columnCounter = 0

    ! read in each character one at a time
    do i = 1, numTimeSteps
       do j = 1, property_per_particle
          do k = 1, 64
             singleValue(k:k) = ' '
          end do

          charCounter = 1
          isInt = .true.
          do
             read(fileID, "(A1)", iostat=ioStatus, advance='no') singleChar

             if (singleChar == '\n') then ! end of line reached
                print *, "eol"
                exit
             end if

             columnCounter = columnCounter + 1

             if (ioStatus > 0) then
                print *, "error while reading input file: ", ioStatus
                print *, "particle id: ", particle_id
                print *, "valueCounter: ", valueCounter
                print *, "columnCounter: ", columnCounter
                print *, "charCounter: ", charCounter
                print *, "singleChar: '", singleChar, "'"
                stop
             end if


             if ((singleChar == ' ') .or. (singleChar == '\t')) then ! ignore whitespace
                if (charCounter > 1) then ! end of number, delimiter reached
                   exit
                else
                   continue ! no valid data found yet
                end if
             else
                if ((singleChar == '.') .or. (singleChar == 'N')) then
                   isInt = .false.
                end if
                singleValue(charCounter:charCounter) = singleChar ! store data
                charCounter = charCounter + 1 ! and go to next character
             end if
          end do

!          print *, "singleValue: '", trim(singleValue), "', charCounter: ", charCounter

          if (isInt) then
             read(singleValue, "(I24)") intValue
             particle_values(j, i) = intValue
          else
             read(singleValue, "(D24.8)") particle_values(j, i)
          end if

          if (particle_values(j, i) /= particle_values(j, i)) then
!             print *, "particle id: ", particle_id, particle_values(j, i)

             if (i > 1) then
                particle_values(j, i) = particle_values(j, i - 1)
!                print *, "replaced with: ", particle_values(j, i - 1)
             else
                ! set everything to zero since we can't replace the NaN with s.th useful
                particle_values = 0.0
                return
             end if

          endif

          valueCounter = valueCounter + 1

       end do
    end do

!    print *, "valueCounter: ", valueCounter

!    print *, "first values: ", particle_values(:, 1)
!    print *, "last values: ", particle_values(:, numTimeSteps)
!    print *, ""

  end subroutine read_particle_data
end module m_file_io
