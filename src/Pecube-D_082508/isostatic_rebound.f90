      subroutine isostatic_rebound (topoa,topob,xsurf,ysurf,nsurf,rsurf,rhoc,rhom, &
                                    nx,ny,nxiso,nyiso,xstep,ystep,young,poisson,thickness)

      implicit real*8 (a-h,o-z)

      real*8 topoa(nsurf),topob(nsurf),xsurf(nsurf),ysurf(nsurf),rsurf(nsurf)
      real*8,dimension(:,:),allocatable :: load
      real*8,dimension(:),allocatable :: workx,worky

      allocate (load(nxiso,nyiso))
      allocate (workx(nxiso),worky(nyiso))

      g=9.81

! flexural isostasy

1     continue

      fact=rhoc*g*xstep*ystep
      xl=(nxiso-1)*xstep
      yl=(nyiso-1)*ystep
      pi=3.141592654
      pixl=pi/xl
      piyl=pi/yl
      d=young/12./(1.-poisson**2)*thickness**3
      xk=(rhom-rhoc)*g

      nx1=(nxiso-nx)/2
      ny1=(nyiso-ny)/2
      load=0.d0

        do j=1,ny
          do i=1,nx
          ij=(j-1)*nx+i
          ii=nx1+i-1
          jj=ny1+j-1
          if ((ii-1)*(ii-nxiso).gt.0 .or. (jj-1)*(jj-nyiso).gt.0) then
          print*,ii,jj,i,j,nx1,ny1
          stop 'error in isostatic_rebound '
          endif
          load(ii,jj)=(topoa(ij)-topob(ij))*1.e3*fact
          enddo
        enddo

        do j=1,nyiso
        workx=load(:,j)
        call sinft (workx,nxiso)
        load(:,j)=workx
        enddo

        do i=1,nxiso
        worky=load(i,:)
        call sinft (worky,nyiso)
        load(i,:)=worky
        enddo

      load=load*4./xl/yl

        do j=1,nyiso
        fj=(j*piyl)**2
          do i=1,nxiso
          fi=(i*pixl)**2
          tij=d/xk*(fi**2+2.*fi*fj+fj**2)+1.
          load(i,j)=load(i,j)/xk/tij
          enddo
        enddo

        do i=1,nxiso
        worky=load(i,:)
        call sinft (worky,nyiso)
        load(i,:)=worky
        enddo

        do j=1,nyiso
        workx=load(:,j)
        call sinft (workx,nxiso)
        load(:,j)=workx
        enddo

        do j=1,ny
          do i=1,nx
          ij=(j-1)*nx+i
          ii=nx1+i-1
          jj=ny1+j-1
          if ((ii-1)*(ii-nxiso).gt.0 .or. (jj-1)*(jj-nyiso).gt.0) then
          print*,ii,jj,i,j,nx1,ny1
          stop 'error in isostatic_rebound '
          endif
          rsurf(ij)=load(ii,jj)/1.e3
          enddo
        enddo

      deallocate (load)
      deallocate (workx,worky)

      return
      end
