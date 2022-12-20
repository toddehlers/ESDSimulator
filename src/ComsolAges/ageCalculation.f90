module m_ageCalculation

! Code taken from Pecube and Terra

contains
  subroutine Mad_He(time, temperature, ntime, apatiteHe, zirconHe)
    use m_tridag

    implicit none

! input and output parameters
    integer(4), intent(in) :: ntime
    real(8), intent(in), dimension(ntime) :: time, temperature

    real(8), intent(out) :: apatiteHe, zirconHe

! local variables
    real(8), dimension(:,:), allocatable :: age, diag, sup, inf, f
    real(8), dimension(:), allocatable :: agei, f1, f2, Ft

    real(8) :: a, al, alpha, am, ard, as, beta, D0, D0a2, D0a2rd
    real(8) :: D0al2, D0am2, D0as2, D0z, D0za2, Da2now, Da2now_h
    real(8) :: Da2now_l, Da2now_m, Da2then, Da2then_h, Da2then_l
    real(8) :: Da2then_m, Dal2now, Dal2then, Dam2now, Dam2then
    real(8) :: Das2now, Das2then, dr, dt, dt0, Dza2now, Dza2then
    real(8) :: Ea, Eard, Eaz, Et, fact, fstep, He_h, He_l, He_m
    real(8) :: P_h, P_l, P_m, psi, temp, tempf, temps, Thh
    real(8) :: Thhc, Thl, Thlc, Thm, Thmc, Uh, Uhc, Ul, Ulc, Um
    real(8) :: Umc


    real(8), parameter :: pi = 3.141592
    real(8), parameter :: R = 8.314 ! Universal gas constant [J/(K*mol)]
    real(8), parameter :: S = 20.0 ! Alpha stopping distance [um]

    integer(4) :: i, istep, itime, k, n, nstep
    integer(4), dimension(8)  :: age_flags

    ! Default diffusion kinetic model
    ! 1 = Wolf et al., 1996
    ! 2 = Farley, 2000
    integer(4), parameter :: kmod = 2

    age_flags = 1

!    print *, "age_flags: ", age_flags

    ! Allocates alpha ejection correction array
    allocate (Ft(8))

    if (kmod == 1) then
       ! Kinetics from Wolf et al., 1996
       D0a2=10.**7.7*3600.*24.*365.25e6              ! [1/My]
       Ea=36.2e3*4.184                               ! Activation energy in undamaged apatite grain [J/mol]
    elseif (kmod == 2) then
       ! Kinetic parameters from Farley, 2000 (added by dwhipp 03/07)
       D0=50.                                        ! Diffusivity at infinite temp. [cm^2/s]
       a=42.5                                        ! Diffusion domain [um]
       Ea=33.0e3*4.184                               ! Activation energy in undamaged apatite grain [J/mol]
       Ft(1)=1.-(3.*S)/(4.*a)+(S**3.)/(16.*a**3.)    ! Alpha ejection correction factor []
       ! Convert diffusion kinetics to same units as Wolf et al., 1996 (above)
       D0=D0*0.01*0.01                               ! [m^2/s]
       a=a*1e-6                                      ! [m]
       D0a2=(D0/(a*a))*3600.*24.*365.25e6            ! D0 over a^2 [1/My]
    endif

    ! Zircon He diffusion kinetics (Reiners et al., 2004)
    D0z=0.46                  ! Diffusivity of zircon grain [cm^2/s]
    Eaz=169.*1000.            ! Activation energy in undamaged zircon grain [J/mol]
    
    ! Radiation damaged apatite diffusion kinetics (Schuster et al., 2006)
    Eard=120.e3               ! Activation energy of undamaged apatite [J/mol]
    Et=29.e3                  ! Activation energy in damaged zones of crystal [J/mol]
    psi=1.26e-4               ! Radiation damage constant [gm/nmol]
    D0a2rd=1.58e4*3600.*24.*365.25e6   ! D0/a^2 [1/My]
    ard=60.                   ! Radiation damaged grain size *****ONLY FOR FT CORRECTION***** [um]
    
    ! User-input parameters (this and conversions below added by dwhipp 03/07)
    as=20.            ! Diffusion domain (small grain) [um]
    am=40.            ! Diffusion domain (medium grain)[um]
    al=70.            ! Diffusion domain (large grain) [um]
    Ul=10.            ! Uranium concentration (low damage) [ppm]
    Um=50.            ! Uranium concentration (moderate damage) [ppm]
    Uh=200.           ! Uranium concentration (high damage) [ppm]
    Thl=10.           ! Thorium concentration (low damage) [ppm]
    Thm=50.           ! Thorium concentration (moderate damage) [ppm]
    Thh=200.          ! Thorium concentration (high damage) [ppm]
    

    ! Other parameter conversions
    D0z=D0z*0.01*0.01                                 ! [m^2/s]
    Ft(2)=1.-(3.*S)/(4.*as)+(S**3.)/(16.*as**3.)      ! Alpha ejection correction factor, small grain []
    Ft(3)=1.-(3.*S)/(4.*am)+(S**3.)/(16.*am**3.)      ! Alpha ejection correction factor, medium grain []
    Ft(4)=1.-(3.*S)/(4.*al)+(S**3.)/(16.*al**3.)      ! Alpha ejection correction factor, large grain []
    Ft(5)=1.-(3.*S)/(4.*ard)+(S**3.)/(16.*ard**3.)    ! Alpha ejection correction factor, low r.d. []
    Ft(6)=Ft(5)                                       ! Alpha ejection correction factor, medium r.d. []
    Ft(7)=Ft(5)                                       ! Alpha ejection correction factor, high r.d. []
    Ft(8)=1.                                          ! No alpha ejection correction for ZHe
    as=as*1e-6                                ! [m]
    am=am*1e-6                                ! [m]
    al=al*1e-6                                ! [m]
    D0as2=(D0/(as*as))*3600.*24.*365.25e6     ! D0 over a^2 [1/My]
    D0am2=(D0/(am*am))*3600.*24.*365.25e6     ! D0 over a^2 [1/My]
    D0al2=(D0/(al*al))*3600.*24.*365.25e6     ! D0 over a^2 [1/My]
    D0za2=(D0z/(a*a))*3600.*24.*365.25e6      ! D0 over a^2 [1/My]
    


    ! Convert U/Th concentrations [ppm]-->[nmol/g]
    Ulc=Ul*(1/238.03)*1e3
    Umc=Um*(1/238.03)*1e3
    Uhc=Uh*(1/238.03)*1e3
    Thlc=Thl*(1/232.038)*1e3
    Thmc=Thm*(1/232.038)*1e3
    Thhc=Thh*(1/232.038)*1e3


    ! He production rate (low, med, high damage) [nmol/(gm*My)]
    P_l=(8.*4.95e-18*Ulc+6.*1.56e-18*Thlc)*3600.*24.*365.25e6
    P_m=(8.*4.95e-18*Umc+6.*1.56e-18*Thmc)*3600.*24.*365.25e6
    P_h=(8.*4.95e-18*Uhc+6.*1.56e-18*Thhc)*3600.*24.*365.25e6
    
    ! Define number of iterations and time step
    dt0=0.1
    n=100

    ! Allocate arrays
    allocate (agei(8),age(8,n),diag(8,n),sup(8,n),inf(8,n),f(8,n))
    allocate (f1(8),f2(8))

    ! Initialize age and He concentrations (added by dwhipp 03/07)
    age=0.            ! []
    agei=0.           ! []
    He_l=0.           ! [nmol/gm]
    He_m=0.           ! [nmol/gm]
    He_h=0.           ! [nmol/gm]
    
    ! Loop through iteratively to calculate ages
    do itime=1,ntime-1

       ! Update time step
       dt=dt0
       nstep=max(1,int((time(itime)-time(itime+1)+tiny(dt))/dt))
       dt=(time(itime)-time(itime+1))/nstep
       alpha=0.5
       dr=1./(n-1)

!       print *, "ageCalulation.f90, itime: ", itime, ", nstep: ", nstep

       ! beta determines the geometry: 2 = spherical
       !                               1 = cylindrical
       beta=2.
       
       ! Grab current temperature info
       temps=temperature(itime)
       tempf=temperature(itime+1)
       
       ! Calculate new D0/a^2 (added by cspath 3/12/07)
       Da2now=D0a2*exp(-Ea/R/(temps+273.))                                               ! [1/My]
       Das2now=D0as2*exp(-Ea/R/(temps+273.))                                             ! [1/My]
       Dam2now=D0am2*exp(-Ea/R/(temps+273.))                                             ! [1/My]
       Dal2now=D0al2*exp(-Ea/R/(temps+273.))                                             ! [1/My]
       Da2now_l=(D0a2rd*exp(-Eard/R/(temps+273.)))/((psi*He_l*exp(Et/R/(temps+273.)))+1) ! [1/My]
       Da2now_m=(D0a2rd*exp(-Eard/R/(temps+273.)))/((psi*He_m*exp(Et/R/(temps+273.)))+1) ! [1/My]
       Da2now_h=(D0a2rd*exp(-Eard/R/(temps+273.)))/((psi*He_h*exp(Et/R/(temps+273.)))+1) ! [1/My]
       Dza2now=D0za2*exp(-Eaz/R/(temps+273.))                                            ! [1/My]
       
       ! Loop through all timesteps in cooling history to calculate cooling ages
       do istep=1,nstep
          
          ! Determine normalized time and temperature gradient
          fstep=float(istep)/(nstep)
          temp=temps+(tempf-temps)*fstep
          
          ! Calculate various diffusion parameters (added by cspath 03/07)
          ! Calculate new diffusion parameters (standard grain)
          if (age_flags(1).eq.1) then
             Da2then=Da2now
             Da2now=D0a2*exp(-Ea/R/(temp+273.))
             f1(1)=alpha*dt*Da2now/dr**2
             f2(1)=(1.-alpha)*dt*Da2then/dr**2
          endif
          
          ! Calculate new diffusion parameters (small grain)
          if (age_flags(2).eq.1) then
             Das2then=Das2now
             Das2now=D0as2*exp(-Ea/R/(temp+273.))
             f1(2)=alpha*dt*Das2now/dr**2
             f2(2)=(1.-alpha)*dt*Das2then/dr**2
          endif
          
          ! Calculate new diffusion parameters (medium grain)
          if (age_flags(3).eq.1) then
             Dam2then=Dam2now
             Dam2now=D0am2*exp(-Ea/R/(temp+273.))
             f1(3)=alpha*dt*Dam2now/dr**2
             f2(3)=(1.-alpha)*dt*Dam2then/dr**2
          endif
          
          ! Calculate new diffusion parameters (large grain)
          if (age_flags(4).eq.1) then
             Dal2then=Dal2now
             Dal2now=D0al2*exp(-Ea/R/(temp+273.))
             f1(4)=alpha*dt*Dal2now/dr**2
             f2(4)=(1.-alpha)*dt*Dal2then/dr**2
          endif
          
          ! Calculate new diffusion parameters (low damage)
          if (age_flags(5).eq.1) then
             Da2then_l=Da2now_l
             Da2now_l=(D0a2rd*exp(-Eard/R/(temp+273.)))/((psi*He_l*exp(Et/R/(temp+273.)))+1)
             f1(5)=alpha*dt*Da2now_l/dr**2
             f2(5)=(1.-alpha)*dt*Da2then_l/dr**2
          endif

          ! Calculate new diffusion parameters (moderate damage)
          if (age_flags(6).eq.1) then
             Da2then_m=Da2now_m
             Da2now_m=(D0a2rd*exp(-Eard/R/(temp+273.)))/((psi*He_m*exp(Et/R/(temp+273.)))+1)
             f1(6)=alpha*dt*Da2now_m/dr**2
             f2(6)=(1.-alpha)*dt*Da2then_m/dr**2
          endif
          
          ! Calculate new diffusion parameters (high damage)
          if (age_flags(7).eq.1) then
             Da2then_h=Da2now_h
             Da2now_h=(D0a2rd*exp(-Eard/R/(temp+273.)))/((psi*He_h*exp(Et/R/(temp+273.)))+1)
             f1(7)=alpha*dt*Da2now_h/dr**2
             f2(7)=(1.-alpha)*dt*Da2then_h/dr**2
          endif
          
          ! Calculate new diffusion parameters (zircon grain)
          if (age_flags(8).eq.1) then
             Dza2then=Dza2now
             Dza2now=D0za2*exp(-Eaz/R/(temp+273.))
             f1(8)=alpha*dt*Dza2now/dr**2
             f2(8)=(1.-alpha)*dt*Dza2then/dr**2
          endif
          
          ! Fill 2-D arrays to hold the seven age calculations (added by cspath 03/07)
          do k=1,8
             if (age_flags(k).eq.1) then
                do i=2,n-1
                   diag(k,i)=1.+2.*f1(k)
                   sup(k,i)=-f1(k)*(1.+beta/(i-1)/2.)
                   inf(k,i)=-f1(k)*(1.-beta/(i-1)/2.)
                   f(k,i)=age(k,i)+f2(k)* &
                        ((age(k,i+1)-2.*age(k,i)+age(k,i-1)) + &
                        beta*(age(k,i+1)-age(k,i-1))/(i-1)/2.)+dt
                enddo
             endif
          enddo
          
          ! Calculate the seven different ages (added by cspath 03/07)
          do k=1,8
             if (age_flags(k).eq.1) then
                diag(k,1)=1.
                sup(k,1)=-1.
                f(k,1)=0.
                diag(k,n)=1.
                inf(k,n)=0.
                f(k,n)=0.
                
                ! Invert tri-diagonal matrices
                call tridag(inf(k,:),diag(k,:),sup(k,:),f(k,:),age(k,:),n)
                
                ! Calculate apparent age for current time step
                agei(k)=0.
                do i=1,n
                   fact=1.
                   if (i.eq.1.or.i.eq.n) fact=0.5
                   agei(k)=agei(k)+age(k,i)*fact*dr**3*(i-1)*(i-1)
                enddo
                agei(k)=3.*agei(k)!*Ft(k)
             endif
          enddo
          
          ! Calculates radiation damage factor as function of total alpha decay
          He_l=He_l+dt*P_l
          He_m=He_m+dt*P_m
          He_h=He_h+dt*P_h
          
       enddo
    enddo
    
    ! Store ages array in new array
    apatiteHe = agei(1)
    zirconHe = agei(8)
    
    deallocate (age,diag,sup,inf,f)
    deallocate (agei,f1,f2)
    
  end subroutine Mad_He

  subroutine ZFT(time, temperature, ntime, age) ! WK: TODO

    implicit none

! input and output parameters
    integer(4), intent(in) :: ntime

    real(8), intent(in), dimension(ntime) :: time, temperature

    real(8), intent(out) :: age

! local variables

    real(8) :: a, B, closure, ratio, cooling, D0, diff, energy, energy_cal, geom
    real(8) :: r, tau, tempp, closurep

    integer(4) :: i, kmod
    
    ! Kinetic parameters to use:
    ! 1=Batt et al., 2001
    ! 2=
    ! 3=
    kmod = 1
    
    if (kmod == 1) then
       energy_cal = 49.77        ! Activation energy [kcal/mol]
       B=3.16e-22              ! 1/(D0/a^2) [My]
    elseif (kmod == 2) then
       energy = 0.0              ! Activation energy [kJ/mol]
       D0 = 0.0                  ! D0 [cm^2/s]
       a = 0.0                   ! a [um]
    elseif (kmod == 3) then
       energy = 0.0              ! Activation energy [kJ/mol]
       D0 = 0.0                  ! D0 [cm^2/s]
       a = 0.0                   ! a [um]
    endif
    
    ! Geometry factor for plane sheet
    geom = 55.0
    
    ! Conversion factors
    if (kmod == 1) then
       energy = energy_cal * 4.184 * 1.e3
       diff = 1.0 / B
    else
       energy = energy * 1.e3
       diff = ((D0 * 0.01 * 0.01) / (a * 1e-6)**2) * 365.25 * 24.0 * 3600.0 * 1.e6
    endif
    
    cooling = 0.0 ! unit: [deg / My]
    closure = 0.0 ! unit: deg cel
    closurep = 0.0
    tempp = 0.0
    
    r = 8.314
    age = time(1)

    do i=ntime,2,-1
       if (i.eq.1) then
          cooling = (temperature(i+1) - temperature(i)) / (time(i+1) - time(i))
       elseif (i.eq.ntime) then
          cooling = (temperature(i) - temperature(i-1)) / (time(i) - time(i-1))
       else
          cooling = (temperature(i+1) - temperature(i-1)) / (time(i+1) - time(i-1))
       endif
       ! note that cooling cannot be nil and is therefore forced to
       ! be at least 1deg/10My
       cooling = max(cooling, 1.0/10.0)
       tau = r * (temperature(i) + 273.0)**2 / energy / cooling
       closure = energy / r / log(geom * tau * diff) - 273.0
       if (temperature(i) > closure) then
          ratio = (closurep - tempp) / (closurep - tempp + temperature(i) - closure)
          age = time(i) + (time(i-1) - time(i)) * ratio
          ! print *, 'Closure: ', closure, ", age: ", age, ", temperature: ", temperature(i)
          ! print *, "min temp: ", minval(temperature), ", max temp: ", maxval(temperature)
          ! print *, temperature
          return
       endif
       closurep = closure
       tempp = temperature(i)
    enddo
  end subroutine ZFT

  subroutine Mad_Trax_AFT(time, temperature, ntime, age) ! TODO

! subroutine Mad_Trax to calculate fission track age, track length
! distribution and statistics from a given thermal history
!
! in input:
! real*4  time_i(n)  :   the time values (in Myr) in descending order
!                        at which the thermal history is given
!                        (ex: 100,50,20,10,0); the last value
!                        should always be 0; the first value
!                        should be smaller than 1000.
! real*4  temp_i(n)  :   the thermal history in degree Celsius
! integer n          :   the number of time-temperature pairs used
!                        to describe the temperature history
! integer out_flag   :   =0 only calculate fission track age
!                        =1 also calculate track length distribution
!                           and statistics
! integer param_flag :   =1 uses Laslett et al, 1987 parameters
!                        =2 uses Crowley et al., 1991 Durango parameters
!                        =3 uses Crowley et al., 1991 F-apatite parameters
!
! in output:
! real*4  fta        :   fission track age in Myr
! real*4  ftld(17)   :   normalised track length distribution where
!                        ftld(k) is the percentage of track with
!                        length between k-0.5 and k+0.5 microns
! real*4  ftldmean   :   mean track length in microns
! real*4  ftldsd     :   track length standard deviation in microns
!
! This subroutine is based on the subroutine "ftmod.pas" provided by
! Peter vanderBeek in December 1995. The algorithm is explained in
! Peter's PhD thesis and is based on the work by Lutz and Omar (1991)
!
! References:
!
! VanderBeek, P., 1995.Tectonic evolution of continental rifts, PhD Thesis,
!       Faculty of Earth Sicences, Free University, Amsterdam.
!
! Lutz, T.M. and Omar, G.I., 1991. An inverse method of modeling thermal
!       histories from apatite fission-track data. EPSL, 104, 181-195.
!
! Laslett, G.M., Green, P.F., Duddy, I.R. and Gleadow, A.J.W., 1987. Thermal
!       annealing of fission tracks in apatite 2. A quantitative analysis.
!       Chem. Geol. (Isot. Geosci. Sect.) 65, 1-13.
!
! Crowley, K.D., Cameron, M. and Schaefer, R.L., 1991. Experimental studies
!       of annealing of etched fission tracks in fluorapatite. Geochim.
!       Cosmochim. Acta, 55, 1449-1465.
!

    implicit none

! input and output parameters
    integer(4), intent(in) :: ntime

    real(8), intent(in), dimension(ntime) :: time, temperature

    real(8), intent(out) :: age

! local variables
    real(8), dimension(1000) :: r

    real(8) :: a, b, c0, c1, c2, c3, deltat, dj, dt, gr
    real(8) :: rp, sumdj, temp, tempm, tempp, teq
    real(8) :: time_interval, xfct, xind, time_local

    integer(4) :: i, j, nstep, param_flag


    param_flag = 3

    if (param_flag.eq.1) then

! from Laslett et al., 1987

       a=0.35
       b=2.7
       c0=-4.87
       c1=0.000168
       c2=0.00472416
       c3=0.
       
    elseif (param_flag.eq.2) then
       
! from Crowley et al., 1991 (Durango)

       a=0.49
       b=3.0
       c0=-3.202
       c1=0.0000937
       c2=0.001839
       c3=0.00042
       
    elseif (param_flag.eq.3) then
       
! from Crowley et al., 1991 (F-apatite)
! note: modified to use correct beta value for fluorapatite from
! Crowley et al., 1991 - dwhipp 02/07

       a=0.76
       b=4.30
       c0=-1.508
       c1=0.00002076
       c2=0.0002143
       c3=0.0009967
       
    else

! lc,mod fanning Arrhenius model from Ketcham et al., 1999 table 5e
! *** DO NOT USE: Not yet functional ***
! added by dwhipp (08/07)

       a=-0.05771
       b=-13.218
       c0=-9.0722
       c1=0.00029896
! Original c2 of Ketcham et al., 1999 modified for version of track length reduction
! equation used below
       c2=-15.846
       c3=0.00076370
       
    endif

! unannealed fission track length
    xind=16.

! mean length of spontaneous tracks in standards
    xfct=14.5

! calculate the number of time steps assuming 1My time step length
! if run time > 100 My, else take 100 time steps
    nstep=int(time(1))
!     if (nstep.gt.8000) stop 'Your final time is greater than '//
!    &                        'the age of the Universe...'
!     if (nstep.gt.4500) stop 'Your final time is greater than '//
!    &                        'the age of the Earth...'
!     if (nstep.gt.1000) stop 'Fission track does not work very well '//
!    &                        'for time spans greater than 1Byr...'
    if (nstep.gt.1000) then
       print *, "ageCalculation.f90: nstep > 1000"
       stop 
    end if
    time_interval=1.0
    if (nstep.lt.100) then
       nstep=100
       time_interval=time(1)/100.0
    endif
    deltat=time_interval*1.e6*365.24*24.0*3600.0
    
! calculate final temperature


    tempp = calc_temperature(temperature, time, 0.0_8, ntime) + 273.0
    rp = 0.5


! begining of time stepping

    do i=1,nstep
       
       time_local = float(i) * time_interval

       
! calculate temperature by linear interpolation

       temp = calc_temperature(temperature, time, time_local, ntime) + 273.0


! calculate mean temperature over the time step

       tempm = (temp + tempp) / 2.0


! calculate the "equivalent time", teq

       teq = exp((-c2 / c1) + ((mad_trax_g(rp, a, b) - c0) / c1) * (1.0 / tempm - c3))
       if (i == 1) teq = 0.0

! check if we are not getting too close to r=0
! in which case r remains 0 for all following time steps

       if (dlog(teq + deltat) > &
          (expos(1.0 / b, a) - a * c0 - 1.0) / a / c1 * (1.0 / tempm - c3) - c2 / c1) then

          do j = i, nstep
             r(j) = 0.0
          enddo
          nstep = i
       else
! otherwise calculate reduction in length, r, over the time step, dt

          dt = teq + deltat
          gr = c0 + ((c1 * dlog(dt) + c2) / ((1.0 / tempm) - c3))
! equation for decrease in normalized track length from Ketcham et al., 1999
! note: this is slightly different than above, and may require other code mods to implement
! also note: this is currently not functional
! dwhipp - (08/07)
!        gr=c0+c1*((dlog(dt)+c2)/((1./tempm)-c3))
          r(i) = xinv(gr, a, b)

! update variables for next time step

          tempp = temp
          rp = r(i)

!       print*,i,time,temp,r(i)
       endif
    enddo
    
! all reduction factors for all time steps have been calculated
! now estimate the fission track age by simple summation
! (here it helps to use 1Myr time steps)

    sumdj = 0.0
      
    do i=1,nstep
       
       if (r(i) < 0.35) then
          dj = 0.0
       elseif (r(i) < 0.66) then
          dj = 2.15 * r(i) - 0.76
       else
          dj = r(i)
       endif
       
       sumdj = sumdj + dj
       
    enddo


    age = (xind / xfct) * sumdj * time_interval
    
    end subroutine Mad_Trax_AFT

!---
    function mad_trax_g(r,a,b)
      implicit none

      real(8), intent(in) :: r, a, b
      real(8) :: mad_trax_g
      
      mad_trax_g = (expos((1. -expos(r, b)) / b, a) - 1.0) / a
      
      return
    end function mad_trax_g
    
!---
    function xinv(gr,a,b)
      implicit none
      
      real(8), intent(in) :: gr, a, b
      real(8) :: xinv

      xinv=expos(1. -b * expos(a * gr + 1.0, 1.0 / a), 1.0 / b)
      
      return
    end function xinv
    
!---
    function calc_temperature(temp,time,t,n)
      implicit none
      
      integer(4), intent(in) :: n
      real(8), intent(in), dimension(n) ::  temp, time
      real(8), intent(in) :: t
      real(8) :: calc_temperature, rat
      
      integer(4) :: i

      do i=1, n-1
         if ((t - time(i)) * (t - time(i + 1)) < 0.0) then
            rat = (t - time(i)) / (time(i + 1) - time(i))
            calc_temperature = temp(i) + rat * (temp(i + 1) - temp(i))
            return
         endif
      enddo
      
      calc_temperature = temp(n)
      
      return
    end function calc_temperature
    
!---
    function expos (x,a)
      implicit none

      real(8), intent(in) :: a, x
      real(8) :: expos
      
      expos = exp(a * dlog(x))
      
      return
    end function expos
    




    subroutine Muscovite (time, temp, ntime, age)
      implicit none

! input and output parameters
      integer(4), intent(in) :: ntime

      real(8), intent(in), dimension(ntime) :: temp, time
      real(8), intent(out) :: age

! local variables
      real(8) :: a, closure, closurep, D0, diff_sec, energy_cal
      real(8) :: geom, r, ratio, tau, tempp, cooling, diff, energy

      integer(4) :: i, kmod

      ! Kinetic parameters to use:
      ! 1 = Unknown default
      ! 2 = Hames and Bowring, 1994
      ! 3 = Robbins, 1972; Hames and Bowring, 1994
      kmod = 3
      
      if (kmod == 1) then
         energy_cal = 41.8_8       ! Activation energy [kcal/mol]
         diff_sec = 0.963_8        ! D0/a^2 [1/s]
      elseif (kmod == 2) then
         energy = 183.0_8          ! Activation energy [kJ/mol]
         D0 = 0.033_8              ! D0 [cm^2/s]
         a = 750.0_8               ! a [um]
      elseif (kmod == 3) then
         energy = 180.0_8          ! Activation energy [kJ/mol]
         D0 = 4.e-4_8              ! D0 [cm^2/s]
         a = 750.0_8               ! a [um]
      endif
      
      ! Geometry factor for plane sheet
      geom = 8.65_8
      
      ! Conversion factors
      if (kmod == 1) then
         energy = energy_cal * 4.184_8 * 1.e3_8
         diff = diff_sec * 365.25_8 * 24.0_8 * 3600.0_8 * 1.e6_8
      else
         energy = energy * 1.e3_8
         diff = ((D0 * 0.01_8 * 0.01_8) / (a * 1e-6_8)**2.0_8) * 365.25_8 * 24.0_8 * 3600.0_8 * 1.e6_8
      endif
      
      cooling = 0.0_8
      closure = 0.0_8
      closurep = 0.0_8
      tempp = 0.0_8
      
      r = 8.314_8
      
      age = time(1)
      do i = ntime,2,-1
         if (i == 1) then
            cooling = (temp(i+1) - temp(i)) / (time(i+1) - time(i))
         elseif (i == ntime) then
            cooling = (temp(i) - temp(i-1)) / (time(i) - time(i-1))
         else
            cooling = (temp(i+1) - temp(i-1)) / (time(i+1) - time(i-1))
         endif
         ! note that cooling cannot be nil and is therefore forced to
         ! be at least 1deg/10My
         cooling = max(cooling, 1.0_8 / 10.0_8)
         tau = r * (temp(i) + 273.0_8)**2 / energy / cooling
         closure = energy / r / log(geom * tau * diff) - 273.0_8
         if (temp(i) > closure) then
!            print *, 'Closure: ',closure
            ratio = (closurep - tempp) / (closurep - tempp + temp(i) - closure)
            age = time(i) + (time(i-1) - time(i)) * ratio
            return
         endif
         closurep=closure
         tempp=temp(i)
      enddo
    end subroutine Muscovite

end module m_ageCalculation
