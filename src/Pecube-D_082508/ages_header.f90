
subroutine ages_header (nsurf,xlonmin,xlonmax,xsurf,xmin,xmax,&
                        xlatmin,xlatmax,ysurf,ymin,ymax,&
                        zsurf_all,fteg,ftag,run,nstep,ftmg,&
                        header_info,ftzg,thermflag,age_flags,&
                        nrun,nx,ny,nelemsurf)

  ! Subroutine for Pecube that writes ages from every time step to
  ! Tecplot formatted files called Ages_tec000.dat, Ages_tec001.dat, etc

  real*8 xlonmin,xlonmax,xmin,xmax,xlatmin,xlatmax,ymin,ymax,xxx
  real*8 xsurf(nsurf),ysurf(nsurf),thermflag(nstep+1),yyy
  real*8,dimension(nstep,nsurf) :: zsurf_all
  real*4 fteg(nstep,nsurf,8),ftag(nstep,nsurf),ftmg(nstep,nsurf)
  real*4 ftzg(nstep,nsurf)
  real*8 header_info(6)
  character run*100,count3*3,as*5,am*5,al*5,eUl*5,eUm*5,eUh*5
  integer age_flags(11),age_count,nrun,nx,ny,counter,nelemsurf

  age_count = 0
  do m=1,11
    if (age_flags(m).eq.1) age_count = age_count + 1
  enddo

  do k=1,nstep
    if (thermflag(k).eq.1) then
      write (count3,'(i3)') k
      if (k.lt.100) count3(1:1)='0'
      if (k.lt.10) count3(1:2)='00'

      open (92,file=run(1:nrun)//'/Ages_tec'//count3//'.dat',status='unknown')

      ! Writes the various grain sizes and radiation levels to
      ! character arrays printed to ages output Tecplot file
      write (as,'(f5.1)') header_info(1)
      write (am,'(f5.1)') header_info(2)
      write (al,'(f5.1)') header_info(3)
      write (eUl,'(f5.1)') header_info(4)
      write (eUm,'(f5.1)') header_info(5)
      write (eUh,'(f5.1)') header_info(6)

      ! Writes Tecplot header in Ages_tec*.dat file
      write (92,*) 'TITLE = "Pecube Ages"'
      write (92,'(A)',ADVANCE="no") 'VARIABLES = "x (km)" "y (km)" "z (km)"'

      if (age_flags(1).eq.1) write (92,'(A)',ADVANCE="no") ' "AHe Age - Farley, 2000 (Ma)"'
      if (age_flags(2).eq.1) write (92,'(A)',ADVANCE="no") ' "AHe Age - a='//as//' um (Ma)"'
      if (age_flags(3).eq.1) write (92,'(A)',ADVANCE="no") ' "AHe Age - a='//am//' um (Ma)"'
      if (age_flags(4).eq.1) write (92,'(A)',ADVANCE="no") ' "AHe Age - a='//al//' um (Ma)"'
      if (age_flags(5).eq.1) write (92,'(A)',ADVANCE="no") ' "AHe Age - eU='//eUl//' ppm (Ma)"'
      if (age_flags(6).eq.1) write (92,'(A)',ADVANCE="no") ' "AHe Age - eU='//eUm//' ppm (Ma)"'
      if (age_flags(7).eq.1) write (92,'(A)',ADVANCE="no") ' "AHe Age - eU='//eUh//' ppm (Ma)"'
      if (age_flags(9).eq.1) write (92,'(A)',ADVANCE="no") ' "AFT Age (Ma)"'
      if (age_flags(8).eq.1) write (92,'(A)',ADVANCE="no") ' "ZHe Age (Ma)"'
      if (age_flags(10).eq.1) write (92,'(A)',ADVANCE="no") ' "ZFT Age (Ma)"'
      if (age_flags(11).eq.1) write (92,'(A)',ADVANCE="no") ' "MAr Age (Ma)"'

      write (92,*)
      write (92,*) 'ZONE T="Ages"'
      write(92,'(A2,i10)',advance="no") 'n=',nsurf
      write(92,'(A4,i10)',advance="no") ', e=',nelemsurf
      write(92,*) ', et=quadrilateral, f=fepoint'

!       write (92,*) 'I=',nsurf,', J=1, K=1, ZONETYPE=Ordered'
!       write (92,*) 'DATAPACKING=POINT'
!       write (92,'(A)',ADVANCE="no") ' DT=(DOUBLE DOUBLE DOUBLE'
!       do i=1,age_count-1
!         write (92,'(A)',ADVANCE="no") ' DOUBLE'
!       enddo
!       if (age_count.ne.0) then
!         write (92,'(A)') ' DOUBLE)'
!       else
!         write (92,'(A)') ')'
!       endif

      ! Loop to write latitude,longitude,elevation,AHe age, and
      ! AFT age. This loop is taken directly from Pecube.f90
      ! Note: The x, y, and z values may not be accurate with
      ! ages at every time step
      ! The arrays xdepth, ydepth, and zdepth may be accurate
      do i=1,nsurf
        xxx=xlonmin+(xlonmax-xlonmin)*(xsurf(i)-xmin)/(xmax-xmin)
        yyy=xlatmin+(xlatmax-xlatmin)*(ysurf(i)-ymin)/(ymax-ymin)
        write (92,'(3f12.4)',ADVANCE="no") xxx,yyy,zsurf_all(k,i)
        do j=1,7
          if (age_flags(j).eq.1) then
            write (92,'(f12.4)',ADVANCE="no") fteg(k,i,j)
          endif
        enddo
        if (age_flags(9).eq.1) write (92,'(f12.4)',ADVANCE="no") ftag(k,i)
        if (age_flags(8).eq.1) write (92,'(f12.4)',ADVANCE="no") fteg(k,i,8)
        if (age_flags(10).eq.1) write (92,'(f12.4)',ADVANCE="no") ftzg(k,i)
        if (age_flags(11).eq.1) write (92,'(f12.4)',ADVANCE="no") ftmg(k,i)
        write(92,*)
      enddo

      counter=0
      do i=1,nelemsurf
        counter=counter+1
        if (mod(counter,nx).eq.0) counter=counter+1
        write (92,*) counter+nx,counter+nx+1,counter+1,counter
      enddo

      close(92)

    endif
  enddo

  return
end
