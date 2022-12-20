c change_sea_level

      subroutine change_sea_level (time,sea_level)

c subroutine to change sea level through time

c INPUT: time      = current time (in yrs)
c        sea_level = current sea_level (in m)

c OUTPUT : sea_level  = new sea_level

c subroutines called:
c NONE

      common /vocal/ ivocal

      sea_level=0.

      return
      end
