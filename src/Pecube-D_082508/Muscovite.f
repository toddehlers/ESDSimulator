      subroutine Muscovite (temp,time,ntime,age)

      real temp(ntime),time(ntime)

      ! Kinetic parameters to use:
      ! 1=Unknown default
      ! 2=Hames and Bowring, 1994
      ! 3=Robbins, 1972; Hames and Bowring, 1994
      kmod=3

      if (kmod.eq.1) then
        energy_cal=41.8		! Activation energy [kcal/mol]
        diff_sec=0.963		! D0/a^2 [1/s]
      elseif (kmod.eq.2) then
        energy=183.		! Activation energy [kJ/mol]
        D0=0.033		! D0 [cm^2/s]
        a=750			! a [um]
      elseif (kmod.eq.3) then
        energy=180.		! Activation energy [kJ/mol]
        D0=4.e-4		! D0 [cm^2/s]
        a=750			! a [um]
      endif

      ! Geometry factor for plane sheet
      geom=8.65

      ! Conversion factors
      if (kmod.eq.1) then
        energy=energy_cal*4.184*1.e3
        diff=diff_sec*365.25*24.*3600.*1.e6
      else
        energy=energy*1.e3
        diff=((D0*0.01*0.01)/(a*1e-6)**2)*365.25*24.*3600.*1.e6
      endif

      cooling=0.
      closure=0.

      r=8.314
      age=time(1)
        do i=ntime,2,-1
          if (i.eq.1) then
          cooling=(temp(i+1)-temp(i))/(time(i+1)-time(i))
          elseif (i.eq.ntime) then
          cooling=(temp(i)-temp(i-1))/(time(i)-time(i-1))
          else
          cooling=(temp(i+1)-temp(i-1))/(time(i+1)-time(i-1))
          endif
c note that cooling cannot be nil and is therefore forced to
c be at least 1deg/10My
        cooling=max(cooling,1./10.)
        tau=r*(temp(i)+273.)**2/energy/cooling
        closure=energy/r/log(geom*tau*diff)-273.
        !print *, 'Closure: ',closure
          if (temp(i).gt.closure) then
          ratio=(closurep-tempp)/(closurep-tempp+temp(i)-closure)
          age=time(i)+(time(i-1)-time(i))*ratio
          return
          endif
        closurep=closure
        tempp=temp(i)
        enddo

      return
      end
