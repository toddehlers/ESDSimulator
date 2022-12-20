!******************************************************************
!   Inputs now include geoflag,theta,phi,mft_ratein,mct_ratein,   *
!   st_ratein,mbt_ratein                                          *
!                                                                 *
!   Calls (2) to geometry.f90 now include geoflag,theta,phi,      *
!   mft_ratein,mct_ratein,stf_ratein,mbt_ratein                   *
!   LAST MODIFIED ON: September 10, 2007                          *
!   MODIFIED BY: Dave Whipp                                       *
!******************************************************************

      subroutine make_matrix (mpe,ael,bel,icon,xg,yg,zg,xgp,ygp,zgp,kfix,&
                              diffusivity,heatproduction,Peclet,x1f,y1f,x2f,y2f,&
                              def,dif,alpha0,dt0,time,tempp,nnode,istatic,zl,&
                              xming,xmaxg,yming,ymaxg,rebound,geoflag,theta,phi,&
                              mft_ratein,mbt_ratein,mct_ratein,stf_ratein,tlapse,&
                              tmsl,efold,rhocrust,elenum,nelem,eheat,econd,isdiff,&
                              shdiff,lhdiff,ghdiff,thdiff,ishp,shhp,lhhp,ghhp,thhp,&
                              shearint,friction,Peclet2)

! this is where the FE matrices and RHS vectors are built
      implicit real*8 (a-h,o-z)
      real*8 xg(nnode),yg(nnode),zg(nnode)
      real*8 xgp(nnode),ygp(nnode),zgp(nnode)
      real*8 tempp(nnode),eheat(nelem),econd(nelem)
      real*8 ael(mpe,mpe),bel(mpe)
      real*8 jcb(3,3),jcbi(3,3),jcbp(3,3),jcbip(3,3)
      integer icon(mpe),kfix(nnode),elenum
      real*8,dimension(:),allocatable :: rr,ss,tt,ww
      real*8,dimension(:,:),allocatable :: aelp,b,bp,vv
      real*8,dimension(:),allocatable :: x,y,z,xp,yp,zp,h,ht,dhdr,dhds,dhdt,rb
      real*8 ishp,shhp,lhhp,ghhp,thhp,ishpin,shhpin,lhhpin,ghhpin,thhpin
      real*8 istc,shtc,lhtc,ghtc,thtc,isdiff,shdiff,lhdiff,ghdiff,thdiff
      real*8 vergin,mft_dipin,mbt_dipin,mct_dipin,stf_dipin,mht1_dipin
      real*8 mht2_dipin,mht3_dipin,mht4_dipin,mfts,mfti,mbts,mbti,mcts,mcti
      real*8 sfts,stfi,mht1s,mht1i,mht2s,mht2i,mht3s,mht3i,mht4s,mht4i
      real*8 kb0s,kb0i,kb1s,kb1i,kb2s,kb2i,kb3s,kb3i,kb4s,kb4i
      real*8 kb5s,kb5i,kb6s,kb6i,pie,halfpi,mftdip,mbtdip,mctdip
      real*8 stfdip,mht1dip,mht2dip,mht3dip,mht4dip,verge,tprev,heatnow,condnow
      real*8 kb7as,kb7ai,kb7bs,kb7bi,heatsum,condsum
      integer znum,shearint,shup,lhup
      real*8 friction,stf2as,stf2ai
      real*8 stf2bs,stf2bi

      allocate (aelp(mpe,mpe),b(mpe,3),bp(mpe,3),rb(mpe))
      allocate (x(mpe),y(mpe),z(mpe),xp(mpe),yp(mpe),zp(mpe))
      allocate (vv(mpe,3))
      allocate (h(mpe),ht(mpe),dhdr(mpe),dhds(mpe),dhdt(mpe))

      eps=tiny(eps)
      pie=atan(1.)*4.                                                           ! define pi

! Set nint = number of nodes per element
! Brick: 8; Triangular: 6
        if (mpe.eq.8) then
          nint=8
        else
          nint=6
        endif


! Define material constants
      rhoc=1.
      rho=2700.
      c=1000.
      fric=friction
 !      rhog=9.81*0.1 (original version - 02/08)
      rhog=9.81*tan(pie/6.)*0.05                                                ! rhog for version with shear heating
 !      rhog=0.                                                                  ! rhog for version without shear heating
      cond=diffusivity
      heat=heatproduction
      velo=Peclet
      dynamic=1.
      dt=dt0
      alpha=alpha0
        if (istatic.eq.1) then
        dynamic=0.
        dt=1.
        alpha=1.
        endif

! Define Gauss integration points for triangular domain
      allocate (rr(nint),ss(nint),tt(nint),ww(nint))

      if (nint.eq.6) then
      rr(1)=0.16667
      ss(1)=0.16667
      tt(1)=-.57735
      ww(1)=1.
      rr(2)=0.66667
      ss(2)=0.16667
      tt(2)=-.57735
      ww(2)=1.
      rr(3)=0.16667
      ss(3)=0.66667
      tt(3)=-.57735
      ww(3)=1.
      rr(4)=0.16667
      ss(4)=0.16667
      tt(4)=.57735
      ww(4)=1.
      rr(5)=0.66667
      ss(5)=0.16667
      tt(5)=.57735
      ww(5)=1.
      rr(6)=0.16667
      ss(6)=0.66667
      tt(6)=.57735
      ww(6)=1.

! Define Gauss integration points and weights for brick element
      elseif (nint.eq.8) then
      rr(1)=-.57735
      ss(1)=-.57735
      tt(1)=-.57735
      ww(1)=1.
      rr(2)=.57735
      ss(2)=-.57735
      tt(2)=-.57735
      ww(2)=1.
      rr(3)=.57735
      ss(3)=.57735
      tt(3)=-.57735
      ww(3)=1.
      rr(4)=-.57735
      ss(4)=.57735
      tt(4)=-.57735
      ww(4)=1.
      rr(5)=-.57735
      ss(5)=-.57735
      tt(5)=.57735
      ww(5)=1.
      rr(6)=.57735
      ss(6)=-.57735
      tt(6)=.57735
      ww(6)=1.
      rr(7)=.57735
      ss(7)=.57735
      tt(7)=.57735
      ww(7)=1.
      rr(8)=-.57735
      ss(8)=.57735
      tt(8)=.57735
      ww(8)=1.
      endif

! Define nodal coordinates for current element
! icon is the element connectivity matrix
        do k=1,mpe
        x(k)=xg(icon(k))
        y(k)=yg(icon(k))
        z(k)=zg(icon(k))
        xp(k)=xgp(icon(k))
        yp(k)=ygp(icon(k))
        zp(k)=zgp(icon(k))
        enddo

! Initialize a, ap and b element matricies
      ael=0.
      aelp=0.
      bel=0.

! Initialize material property sum variables
      heatsum=0.
      condsum=0.

      timep=time-dt

! Velocity values should change when istep increments in Pecube.f90
! Values should stay the same for every itime value within each istep
! Pass itime and istep, or new variable to flag when find_velo should or
! should not calculate new velocities

        do i=1,mpe
          call find_velo (x(i),y(i),z(i),dt,velo,x1f,y1f,x2f,y2f,def,dif,time,&
          xming,xmaxg,yming,ymaxg,zl,vv(i,1),vv(i,2),vv(i,3),geoflag,theta,phi,&
          mft_ratein,mbt_ratein,mct_ratein,stf_ratein,Peclet2)
        enddo

! Start loop through number of nodes per element
        do iint=1,nint

! Define Gauss points and weights for current node
        r=rr(iint)
        s=ss(iint)
        t=tt(iint)
        w=ww(iint)

! Define shape functions for brick and triangular elements
          if (mpe.eq.8) then
          h(1)=(1.-r)*(1.-s)*(1.-t)/8.
          h(2)=(1.+r)*(1.-s)*(1.-t)/8.
          h(3)=(1.+r)*(1.+s)*(1.-t)/8.
          h(4)=(1.-r)*(1.+s)*(1.-t)/8.
          h(5)=(1.-r)*(1.-s)*(1.+t)/8.
          h(6)=(1.+r)*(1.-s)*(1.+t)/8.
          h(7)=(1.+r)*(1.+s)*(1.+t)/8.
          h(8)=(1.-r)*(1.+s)*(1.+t)/8.
          else
          h(1)=(1.-r-s)*(1.-t)/2.
          h(2)=r*(1.-t)/2.
          h(3)=s*(1.-t)/2.
          h(4)=(1.-r-s)*(1.+t)/2.
          h(5)=r*(1.+t)/2.
          h(6)=s*(1.+t)/2.
          endif

! Not sure what is going on here, but looks like multiplying shape functions
! by nodal positions (loading?)
        xint=0.
        yint=0.
        zint=0.
        rint=0.
          do k=1,mpe
          xint=xint+h(k)*x(k)
          yint=yint+h(k)*y(k)
          zint=zint+h(k)*z(k)
          rint=rint+h(k)*rb(k)
          enddo

! Commented out because veloz is redefined below...not sure why this was here.
! dwhipp - 10/07
!        veloz=veloz+rint/dt

! Define shape function derivatives for brick and triangular elements
          if (mpe.eq.8) then
          dhdr(1)=-(1.-s)*(1.-t)/8.
          dhdr(2)=(1.-s)*(1.-t)/8.
          dhdr(3)=(1.+s)*(1.-t)/8.
          dhdr(4)=-(1.+s)*(1.-t)/8.
          dhdr(5)=-(1.-s)*(1.+t)/8.
          dhdr(6)=(1.-s)*(1.+t)/8.
          dhdr(7)=(1.+s)*(1.+t)/8.
          dhdr(8)=-(1.+s)*(1.+t)/8.
          dhds(1)=-(1.-r)*(1.-t)/8.
          dhds(2)=-(1.+r)*(1.-t)/8.
          dhds(3)=(1.+r)*(1.-t)/8.
          dhds(4)=(1.-r)*(1.-t)/8.
          dhds(5)=-(1.-r)*(1.+t)/8.
          dhds(6)=-(1.+r)*(1.+t)/8.
          dhds(7)=(1.+r)*(1.+t)/8.
          dhds(8)=(1.-r)*(1.+t)/8.
          dhdt(1)=-(1.-r)*(1.-s)/8.
          dhdt(2)=-(1.+r)*(1.-s)/8.
          dhdt(3)=-(1.+r)*(1.+s)/8.
          dhdt(4)=-(1.-r)*(1.+s)/8.
          dhdt(5)=(1.-r)*(1.-s)/8.
          dhdt(6)=(1.+r)*(1.-s)/8.
          dhdt(7)=(1.+r)*(1.+s)/8.
          dhdt(8)=(1.-r)*(1.+s)/8.
          else
          dhdr(1)=-(1.-t)/2.
          dhdr(2)=(1.-t)/2.
          dhdr(3)=0.
          dhdr(4)=-(1.+t)/2.
          dhdr(5)=(1.+t)/2.
          dhdr(6)=0.
          dhds(1)=-(1.-t)/2.
          dhds(2)=0.
          dhds(3)=(1.-t)/2.
          dhds(4)=-(1.+t)/2.
          dhds(5)=0.
          dhds(6)=(1.+t)/2.
          dhdt(1)=-(1.-r-s)/2.
          dhdt(2)=-r/2.
          dhdt(3)=-s/2.
          dhdt(4)=(1.-r-s)/2.
          dhdt(5)=r/2.
          dhdt(6)=s/2.
          endif

! Calculate Jacobian of shape function derivatives (?)
        jcb=0.
        jcbp=0.
          do k=1,mpe
          jcb(1,1)=jcb(1,1)+dhdr(k)*x(k)
          jcb(1,2)=jcb(1,2)+dhdr(k)*y(k)
          jcb(1,3)=jcb(1,3)+dhdr(k)*z(k)
          jcb(2,1)=jcb(2,1)+dhds(k)*x(k)
          jcb(2,2)=jcb(2,2)+dhds(k)*y(k)
          jcb(2,3)=jcb(2,3)+dhds(k)*z(k)
          jcb(3,1)=jcb(3,1)+dhdt(k)*x(k)
          jcb(3,2)=jcb(3,2)+dhdt(k)*y(k)
          jcb(3,3)=jcb(3,3)+dhdt(k)*z(k)
          jcbp(1,1)=jcbp(1,1)+dhdr(k)*xp(k)
          jcbp(1,2)=jcbp(1,2)+dhdr(k)*yp(k)
          jcbp(1,3)=jcbp(1,3)+dhdr(k)*zp(k)
          jcbp(2,1)=jcbp(2,1)+dhds(k)*xp(k)
          jcbp(2,2)=jcbp(2,2)+dhds(k)*yp(k)
          jcbp(2,3)=jcbp(2,3)+dhds(k)*zp(k)
          jcbp(3,1)=jcbp(3,1)+dhdt(k)*xp(k)
          jcbp(3,2)=jcbp(3,2)+dhdt(k)*yp(k)
          jcbp(3,3)=jcbp(3,3)+dhdt(k)*zp(k)
          enddo

! Calculate element volume for current and previous time step (?)
        volume=jcb(1,1)*jcb(2,2)*jcb(3,3)+jcb(1,2)*jcb(2,3)*jcb(3,1) &
              +jcb(2,1)*jcb(3,2)*jcb(1,3) &
              -jcb(1,3)*jcb(2,2)*jcb(3,1)-jcb(1,2)*jcb(2,1)*jcb(3,3) &
              -jcb(2,3)*jcb(3,2)*jcb(1,1)
        volumep=jcbp(1,1)*jcbp(2,2)*jcbp(3,3)+jcbp(1,2)*jcbp(2,3)*jcbp(3,1) &
               +jcbp(2,1)*jcbp(3,2)*jcbp(1,3) &
               -jcbp(1,3)*jcbp(2,2)*jcbp(3,1)-jcbp(1,2)*jcbp(2,1)*jcbp(3,3) &
               -jcbp(2,3)*jcbp(3,2)*jcbp(1,1)
	
          if (volume.le.eps) then
          print*,'volume= ',volume
          stop 'Element bow-tied or collapsed'
          endif

          if (volumep.le.eps) then
          print*,'volumep= ',volumep
          stop 'Element bow-tied or collapsed'
          endif

! Not sure what this does
        jcbi(1,1)=(jcb(2,2)*jcb(3,3)-jcb(2,3)*jcb(3,2))/volume
        jcbi(2,1)=(jcb(2,3)*jcb(3,1)-jcb(2,1)*jcb(3,3))/volume
        jcbi(3,1)=(jcb(2,1)*jcb(3,2)-jcb(2,2)*jcb(3,1))/volume
        jcbi(1,2)=(jcb(1,3)*jcb(3,2)-jcb(1,2)*jcb(3,3))/volume
        jcbi(2,2)=(jcb(1,1)*jcb(3,3)-jcb(1,3)*jcb(3,1))/volume
        jcbi(3,2)=(jcb(1,2)*jcb(3,1)-jcb(1,1)*jcb(3,2))/volume
        jcbi(1,3)=(jcb(1,2)*jcb(2,3)-jcb(1,3)*jcb(2,2))/volume
        jcbi(2,3)=(jcb(1,3)*jcb(2,1)-jcb(1,1)*jcb(2,3))/volume
        jcbi(3,3)=(jcb(1,1)*jcb(2,2)-jcb(1,2)*jcb(2,1))/volume
        jcbip(1,1)=(jcbp(2,2)*jcbp(3,3)-jcbp(2,3)*jcbp(3,2))/volumep
        jcbip(2,1)=(jcbp(2,3)*jcbp(3,1)-jcbp(2,1)*jcbp(3,3))/volumep
        jcbip(3,1)=(jcbp(2,1)*jcbp(3,2)-jcbp(2,2)*jcbp(3,1))/volumep
        jcbip(1,2)=(jcbp(1,3)*jcbp(3,2)-jcbp(1,2)*jcbp(3,3))/volumep
        jcbip(2,2)=(jcbp(1,1)*jcbp(3,3)-jcbp(1,3)*jcbp(3,1))/volumep
        jcbip(3,2)=(jcbp(1,2)*jcbp(3,1)-jcbp(1,1)*jcbp(3,2))/volumep
        jcbip(1,3)=(jcbp(1,2)*jcbp(2,3)-jcbp(1,3)*jcbp(2,2))/volumep
        jcbip(2,3)=(jcbp(1,3)*jcbp(2,1)-jcbp(1,1)*jcbp(2,3))/volumep
        jcbip(3,3)=(jcbp(1,1)*jcbp(2,2)-jcbp(1,2)*jcbp(2,1))/volumep

! Define b matricies (not sure what this does)
          do k=1,mpe
            do i=1,3
              b(k,i)=jcbi(i,1)*dhdr(k)+jcbi(i,2)*dhds(k)+jcbi(i,3)*dhdt(k)
              bp(k,i)=jcbip(i,1)*dhdr(k)+jcbip(i,2)*dhds(k)+jcbip(i,3)*dhdt(k)
            enddo
          enddo

 ! Define epsij array and shear heating contribution to heat production
           sheat=0.                                                              ! Initialize shear heating variable
           do j=1,3                                                              ! Shear heating is implemented using the same method as by F. Herman, but
             do i=1,3                                                            !    the empirical values for the calculation have been changed to reflect
               epsij=0.                                                          !    those observed by the moderate friction case for Granite of Hansen and
               do k=1,mpe                                                        !    Carter, 1982
                 epsij=epsij+b(k,i)*vv(k,j)                                      ! Strain rate tensor component of interest [1/My]
               enddo
!                heat=heat+abs(epsij)*rhog*zint

! 2011.07.25, WK: shearint is an integer and not a boolean!
               if (shearint == 1) then                                                ! Calculate shear heating if flag is set
 !                stress_bri=fric*rho*9.81*(zl-zint)*1.e3-fric*1.e3*9.81 &       ! Define brittle stress [Pa]
 !                           *(zl-zint)*1.e3
                 stress_bri=fric*rho*9.81*(zsint-zint)*1.e3-fric*1.e3*9.81 &     ! Define brittle stress [Pa]
                            *(zsint-zint)*1.e3                                   ! This version calculates depth from mean surface elevaton for given element
                                                                                 !   (see above)
                 xn=3.4                                                          ! Define xn (empirical exponent) [unitless]
                 A=2.512e-9*(1.e6)**(-xn)                                        ! Define A (Moderate friction case of Hansen and Carter, 1982) [Pa^n/s]
                 expo=exp(-139.e3/8.3144/(273.15+tempp(icon(j))))                ! Define expo [unitless]
                 epsijsec=abs(epsij)/(1.e6*3600.*24.*365.25)                     ! Define epsijsec [1/s]
                 stress_duc=(epsijsec/(A*expo))**(1./xn)                         ! Define ductile stress [Pa]
                 if (epsijsec.gt.epmax) epmax=epsijsec
                 stress_min=stress_bri                                           ! Define brittle stress as minimum stress [Pa]
                 if (stress_min.gt.stress_duc) stress_min=stress_duc             ! Redefine minimum stress if ductile is smaller than brittle [Pa]
                 if (stress_min.gt.50.e6) stress_min=50.e6                       ! Limit maximum stress to 50 MPa [Pa]
                 shearheat=stress_min*abs(epsij)/rho/c                           ! Calculate shear heating [C/My]
                 !heat=heat+abs(epsij)*rhog*zint                                 ! Original version (02/08)
                 sheat=sheat+shearheat                                           ! Sum up shear heating for each tensor component
               else
                 sheat=0.                                                        ! Set shear heating contribution to zero
               endif
            enddo
          enddo

! Velocity values should change when istep increments in Pecube.f90
! Values should stay the same for every itime value within each istep
! Pass itime and istep, or new variable to flag when find_velo should or
! should not calculate new velocities

        call find_velo (xint,yint,zint,dt,velo,x1f,y1f,x2f,y2f,def,dif,time,&
                         xming,xmaxg,yming,ymaxg,zl,vx,vy,vz,geoflag,theta,phi,&
                         mft_ratein,mbt_ratein,mct_ratein,stf_ratein,Peclet2)

        velox=velo*vx
        veloy=velo*vy
        veloz=velo*vz

! Calculate x, y, z ranges and max velocity (from components)
        uvnorm=velox**2+veloy**2+veloz**2
        xmin=minval(x)
        xmax=maxval(x)
        dx=xmax-xmin
        ymin=minval(y)
        ymax=maxval(y)
        dy=ymax-ymin
        zmin=minval(z)
        zmax=maxval(z)
        dz=zmax-zmin

! Calculate tau_a (timescale for advection)
! see Quantitative Thermochronology pg. 96, eqn. 5.49
        if (uvnorm.le.tiny(uvnorm)) then
          tau=0.
        else
          tau=(abs(velox)*dx+abs(veloy)*dy+abs(veloz)*dz)/sqrt(15.)/uvnorm
        endif

! Not certain what this does, uses shape functions, tau and velocities
        do k=1,mpe
          ht(k)=h(k)+tau*(velox*b(k,1)+veloy*b(k,2)+veloz*b(k,3))
        enddo

        if (geoflag.eq.4) then
          !***
          ! Geometry variables (don't mess with these unless you know what you're doing)
          !***
          vergin=180			! Convergence direction from GPS (Jouanne et al., 2004)
          mft_dipin=-40			! From Avouac cross-section (degrees)
          mbt_dipin=-40			! From Avouac cross-section (degrees)
          ! Original MCT dip commented out, see below
          !mct_dipin=-28                ! From Searle & Godin, 2003 - (22-28 deg is range they list)
          mct_dipin=-40                 ! From argument that high relief, physiographic transition zone should be south of MCT kink band plane
                                        !   Also, DeCelles et al., 2001 lists a range of 30-48 deg for MCT dip
          stf_dipin=-28		! Assume same dip as MCT - Searle & Godin, 2003 say the foliation dips less than 30 deg NNE
          mht1_dipin=-4			! From Avouac cross-section (degrees)
          mht2_dipin=-7			! From Avouac cross-section (degrees)
          mht3_dipin=-19		! From Avouac cross-section (degrees)
          mht4_dipin=-10		! From Avouac cross-section (degrees)
          stf_seg=51.3			! Distance from left edge of model to tear in STF (km)

          mfts=-0.70			! MFT slope - From Avouac cross-section (unitless)
          mfti=21.25			! MFT z-intercept - From Avouac cross-section (km)
          mbts=-0.70			! MBT slope - From Avouac cross-section (unitless)
          mbti=40.10			! MBT z-intercept - From Avouac cross-section (km)
          ! Original MCT slope/intercept commented out, see below
          !mcts=-0.49			! MCT slope - From Searle & Godin, 2003 (unitless)
          !mcti=59.34			! MCT z-intercept - Calculated from Searle & Godin, 2003 slope with
					! Lave & Avouac, 2000 cross-section (km)
          mcts=-0.70                    ! MCT slope - From argument that high relief, physiographic transition zone should be south of MCT kink band plane (unitless)
                                        !   Also, DeCelles et al., 2001 lists a range of 30-48 deg for MCT dip
          mcti=82.09                    ! MCT z-intercept - Calculated from above with Avouac cross-section (km)
          stfas=-0.49			! STF-W slope - From Searle & Godin, 2003 (unitless)
          stfai=70.42			! STF-W z-intercept - From Avouac cross-section (km)
          stfbs=-0.49			! STF-E slope - From Searle & Godin, 2003 (unitless)
          stfbi=88.99			! STF-E z-intercept - From Avouac cross-section (km)
          mht1s=-0.07			! MHT section 1 slope - From Avouac cross-section (unitless)
          mht1i=-0.67			! MHT section 1 z-intercept - From Avouac cross-section (km)
          mht2s=-0.12			! MHT section 2 slope - From Avouac cross-section (unitless)
          mht2i=2.72			! MHT section 2 z-intercept - From Avouac cross-section (km)
          mht3s=-0.33			! MHT section 3 slope - From Avouac cross-section (unitless)
          mht3i=24.46			! MHT section 3 z-intercept - From Avouac cross-section (km)
          mht4s=-0.17			! MHT section 4 slope - From Avouac cross-section (unitless)
          mht4i=0.94			! MHT section 4 z-intercept - From Avouac cross-section (km)

          ! Kink band plane slopes/intercepts
          kb0s=2.75			! Kink band plane 0 slope - From Avouac cross-section (unitless)
          kb0i=-84.62			! Kink band plane 0 z-intercept - From Avouac cross-section (km)
          kb1s=2.48			! Kink band plane 1 slope - From Avouac cross-section (unitless)
          kb1i=-89.46			! Kink band plane 1 z-intercept - From Avouac cross-section (km)
          kb2s=10.39			! Kink band plane 2 slope - From Avouac cross-section (unitless)
          kb2i=-679.11			! Kink band plane 2 z-intercept - From Avouac cross-section (km)
          kb3s=2.30			! Kink band plane 3 slope - From Avouac cross-section (unitless)
          kb3i=-154.44			! Kink band plane 3 z-intercept - From Avouac cross-section (km)
          kb4s=4.33			! Kink band plane 4 slope - From Avouac cross-section (unitless)
          kb4i=-459.99			! Kink band plane 4 z-intercept - From Avouac cross-section (km)
          kb5s=3.87			! Kink band plane 5 slope - From Avouac cross-section (unitless)
          kb5i=-588.57			! Kink band plane 5 z-intercept - From Avouac cross-section (km)
          ! Original kb6 slope and intercept commented out, see above for reason why (line ~158)
          !kb6s=2.90			! Kink band plane 6 slope - From Avouac cross-section (unitless)
          !kb6i=-571.37			! Kink band plane 6 z-intercept - From Avouac cross-section (km)
          kb6s=2.14                     ! Kink band plane 6 slope - From Avouac cross-section (unitless)
          kb6i=-358.20                  ! Kink band plane 6 z-intercept - From Avouac cross-section (km)
          kb7as=2.90			! Kink band plane 7a slope - From Avouac cross-section (unitless)
          kb7ai=-660.71			! Kink band plane 7a z-intercept - From Avouac cross-section (km)
          kb7bs=2.90			! Kink band plane 7b slope - From Avouac cross-section (unitless)
          kb7bi=-848.62			! Kink band plane 7b z-intercept - From Avouac cross-section (km)

          yrsec=3.15569259747d7		! Number of seconds in a year [s]

          !***
          !  Convert units of input variables / Additional math
          !***
          pie=3.14159265						! define pi
          halfpi=pie/2							! define pi over two
          mftdip=(pie/180.0)*mft_dipin					! convert from deg to rad
          mbtdip=(pie/180.0)*mbt_dipin					! convert from deg to rad
          mctdip=(pie/180.0)*mct_dipin					! convert from deg to rad
          stfdip=(pie/180.0)*stf_dipin					! convert from deg to rad
          mht1dip=(pie/180.0)*mht1_dipin				! convert from deg to rad
          mht2dip=(pie/180.0)*mht2_dipin				! convert from deg to rad
          mht3dip=(pie/180.0)*mht3_dipin				! convert from deg to rad
          mht4dip=(pie/180.0)*mht4_dipin				! convert from deg to rad
          verge = (pie/180.0)*vergin					! convert from deg to rad
        endif

          do j=1,mpe
            ! Modify heat production distribution to decrease exponentially if efold
            !   does not equal 0.
            ! Note: Decreases exponentially below msl, constant at input value above
            ! dwhipp - 10/07
            if (geoflag.ne.4) then
              if (efold.ne.0.) then
                if (z(j).le.zl) then
                  heatnow=heat*exp((z(j)-zl)/efold)
                else
                  heatnow=heat
                endif
              else
                 heatnow=heat
              endif

            ! Define heat production for Nepal model geometry
            ! dwhipp - 10/07
            else
              ! Shift z-values to be relative to msl
              xmod=x(j)
              ymod=y(j)
              zmod=z(j)-zl

              ! Define material properties for different lithostratigraphic units
              if (ymod.le.(zmod-kb0i)/kb0s) then    				! If true, point is in front of thrust belt
                heatnow=ishp
              else if (ymod.le.(zmod-kb1i)/kb1s) then			! Point is south of KB1
                if (zmod.gt.mfts*ymod+mfti) then    				! Point is in MFT hanging wall
                  heatnow=shhp
                else								! Point is in MFT footwall
                  heatnow=ishp
                endif
              else if (ymod.le.(zmod-kb2i)/kb2s .and. zmod.le.mht1s*ymod+mht1i) then	! Point is in MHT1 footwall
                heatnow=ishp
              else if (zmod.le.mbts*ymod+mbti .and. zmod.gt.mht1s*ymod+mht1i) then	! Point is in MHT1 (MFT) hanging wall
                heatnow=shhp
              else if (ymod.le.(zmod-kb3i)/kb3s .and. zmod.gt.mbts*ymod+mbti) then	! Point is in MBT hanging wall
                heatnow=lhhp
              else if (ymod.le.(zmod-kb4i)/kb4s) then			! Point is south of KB4
                if (zmod.gt.mht2s*ymod+mht2i) then				! Point is in MHT2 (MBT) hanging wall
                  heatnow=lhhp
                else								! Point is in MHT2 footwall
                  heatnow=ishp
                endif
              else if (xmod.le.stf_seg) then					! Point is west of STF tear
                if (ymod.le.(zmod-kb5i)/kb5s) then				! Point is south of KB5
                  if (zmod.gt.stfas*ymod+stfai) then				! Point is in STFa hanging wall
                    heatnow=thhp
                  else if (zmod.gt.mcts*ymod+mcti) then			! Point is in MCT hanging wall
                    heatnow=ghhp
                  else if (zmod.gt.mht3s*ymod+mht3i) then			! Point is in MHT3 hanging wall
                    heatnow=lhhp
                  else								! Point is in MHT3 footwall
                    heatnow=ishp
                  endif
                else if (ymod.le.(zmod-kb6i)/kb6s) then			! Point is south of KB6
                  if (zmod.gt.stfas*ymod+stfai) then				! Point is in STFa hanging wall
                    heatnow=thhp
                  else if (zmod.gt.mcts*ymod+mcti) then			! Point is in MCT hanging wall
                    heatnow=ghhp
                  else if (zmod.gt.mht4s*ymod+mht4i) then			! Point is in MHT4 hanging wall
                    heatnow=lhhp
                  else								! Point is in MHT4 footwall
                    heatnow=ishp
                  endif
                else if (ymod.le.(zmod-kb7ai)/kb7as) then			! Point is south of KB7a
                  if (zmod.gt.stfas*ymod+stfai) then				! Point is in STFa hanging wall
                    heatnow=thhp
                  else if (zmod.gt.mht4s*ymod+mht4i) then			! Point is in MHT4 hanging wall
                    heatnow=ghhp
                  else								! Point is in MHT4 footwall
                    heatnow=ishp
                  endif
                else								! Point is north of KB7a
                  if (zmod.gt.mht4s*ymod+mht4i) then				! Point is in MHT4 hanging wall
                    heatnow=thhp
                  else								! Point is in MHT4 footwall
                    heatnow=ishp
                  endif
                endif
              else								! Point is east of STF tear
                if (ymod.le.(zmod-kb5i)/kb5s) then				! Point is south of KB5
                  if (zmod.gt.stfbs*ymod+stfbi) then				! Point is in STFb hanging wall
                    heatnow=thhp
                  else if (zmod.gt.mcts*ymod+mcti) then			! Point is in MCT hanging wall
                    heatnow=ghhp
                  else if (zmod.gt.mht3s*ymod+mht3i) then			! Point is in MHT3 hanging wall
                    heatnow=lhhp
                  else								! Point is in MHT3 footwall
                    heatnow=ishp
                  endif
                else if (ymod.le.(zmod-kb6i)/kb6s) then			! Point is south of KB6
                  if (zmod.gt.stfbs*ymod+stfbi) then				! Point is in STFb hanging wall
                    heatnow=thhp
                  else if (zmod.gt.mcts*ymod+mcti) then			! Point is in MCT hanging wall
                    heatnow=ghhp
                  else if (zmod.gt.mht4s*ymod+mht4i) then			! Point is in MHT4 hanging wall
                    heatnow=lhhp
                  else								! Point is in MHT4 footwall
                    heatnow=ishp
                  endif
                else if (ymod.le.(zmod-kb7bi)/kb7bs) then			! Point is south of KB7b
                  if (zmod.gt.stfbs*ymod+stfbi) then				! Point is in STFb hanging wall
                    heatnow=thhp
                  else if (zmod.gt.mht4s*ymod+mht4i) then			! Point is in MHT4 hanging wall
                    heatnow=ghhp
                  else								! Point is in MHT4 footwall
                    heatnow=ishp
                  endif
                else								! Point is north of KB7b
                  if (zmod.gt.mht4s*ymod+mht4i) then				! Point is in MHT4 hanging wall
                    heatnow=thhp
                  else								! Point is in MHT4 footwall
                    heatnow=ishp
                  endif
                endif
              endif
            endif

            heatsum=heatsum+heatnow
            bel(j)=bel(j)+h(j)*dt*((1.-alpha)*heatnow*volumep+alpha*heatnow*volume)*w
            do i=1,mpe
              if (geoflag.ne.4) then
                condnow=cond
              else
                ! Shift z-values to be relative to msl
                xmod=x(i)
                ymod=y(i)
                zmod=z(i)-zl

                ! MFT active (perhaps w/ MBT and MCT)
                if (ymod.le.(zmod-kb0i)/kb0s) then    					! If true, point is in front of thrust belt
                  condnow=isdiff
                else if (ymod.le.(zmod-kb1i)/kb1s) then					! Point is south of KB1
                  if (zmod.gt.mfts*ymod+mfti) then    					! Point is in MFT hanging wall
                    condnow=shdiff
                  else										! Point is in MFT footwall
                    condnow=isdiff
                  endif
                else if (ymod.le.(zmod-kb2i)/kb2s .and. zmod.le.mht1s*ymod+mht1i) then	! Point is in MHT1 footwall
                  condnow=isdiff
                else if (zmod.le.mbts*ymod+mbti .and. zmod.gt.mht1s*ymod+mht1i) then	! Point is in MHT1 (MFT) hanging wall
                  condnow=shdiff
                else if (ymod.le.(zmod-kb3i)/kb3s .and. zmod.gt.mbts*ymod+mbti) then	! Point is in MBT hanging wall
                  condnow=lhdiff
                else if (ymod.le.(zmod-kb4i)/kb4s) then					! Point is south of KB4
                  if (zmod.gt.mht2s*ymod+mht2i) then						! Point is in MHT2 (MBT) hanging wall
                    condnow=lhdiff
                  else										! Point is in MHT2 footwall
                    condnow=isdiff
                  endif
                else if (xmod.le.stf_seg) then						! Point is west of STF tear
                  if (ymod.le.(zmod-kb5i)/kb5s) then						! Point is south of KB5
                    if (zmod.gt.stfas*ymod+stfai) then					! Point is in STFa hanging wall
                      condnow=thdiff
                    else if (zmod.gt.mcts*ymod+mcti) then					! Point is in MCT hanging wall
                      condnow=ghdiff
                    else if (zmod.gt.mht3s*ymod+mht3i) then					! Point is in MHT3 hanging wall
                      condnow=lhdiff
                    else									! Point is in MHT3 footwall
                      condnow=isdiff
                    endif
                  else if (ymod.le.(zmod-kb6i)/kb6s) then					! Point is south of KB6
                    if (zmod.gt.stfas*ymod+stfai) then					! Point is in STFa hanging wall
                      condnow=thdiff
                    else if (zmod.gt.mcts*ymod+mcti) then					! Point is in MCT hanging wall
                      condnow=ghdiff
                    else if (zmod.gt.mht4s*ymod+mht4i) then					! Point is in MHT4 hanging wall
                      condnow=lhdiff
                    else									! Point is in MHT4 footwall
                      condnow=isdiff
                    endif
                  else if (ymod.le.(zmod-kb7ai)/kb7as) then					! Point is south of KB7a
                    if (zmod.gt.stfas*ymod+stfai) then					! Point is in STFa hanging wall
                      condnow=thdiff
                    else if (zmod.gt.mht4s*ymod+mht4i) then					! Point is in MHT4 hanging wall
                      condnow=ghdiff
                    else									! Point is in MHT4 footwall
                      condnow=isdiff
                    endif
                  else										! Point is north of KB7a
                    if (zmod.gt.mht4s*ymod+mht4i) then					! Point is in MHT4 hanging wall
                      condnow=thdiff
                    else									! Point is in MHT4 footwall
                      condnow=isdiff
                    endif
                  endif
                else										! Point is east of STF tear
                  if (ymod.le.(zmod-kb5i)/kb5s) then						! Point is south of KB5
                    if (zmod.gt.stfbs*ymod+stfbi) then					! Point is in STFb hanging wall
                      condnow=thdiff
                    else if (zmod.gt.mcts*ymod+mcti) then					! Point is in MCT hanging wall
                      condnow=ghdiff
                    else if (zmod.gt.mht3s*ymod+mht3i) then					! Point is in MHT3 hanging wall
                      condnow=lhdiff
                    else									! Point is in MHT3 footwall
                      condnow=isdiff
                    endif
                  else if (ymod.le.(zmod-kb6i)/kb6s) then					! Point is south of KB6
                    if (zmod.gt.stfbs*ymod+stfbi) then					! Point is in STFb hanging wall
                      condnow=thdiff
                    else if (zmod.gt.mcts*ymod+mcti) then					! Point is in MCT hanging wall
                      condnow=ghdiff
                    else if (zmod.gt.mht4s*ymod+mht4i) then					! Point is in MHT4 hanging wall
                      condnow=lhdiff
                    else									! Point is in MHT4 footwall
                      condnow=isdiff
                    endif
                  else if (ymod.le.(zmod-kb7bi)/kb7bs) then					! Point is south of KB7b
                    if (zmod.gt.stfbs*ymod+stfbi) then					! Point is in STFb hanging wall
                      condnow=thdiff
                    else if (zmod.gt.mht4s*ymod+mht4i) then					! Point is in MHT4 hanging wall
                      condnow=ghdiff
                    else									! Point is in MHT4 footwall
                      condnow=isdiff
                    endif
                  else										! Point is north of KB7b
                    if (zmod.gt.mht4s*ymod+mht4i) then					! Point is in MHT4 hanging wall
                      condnow=thdiff
                    else									! Point is in MHT4 footwall
                      condnow=isdiff
                    endif
                  endif
                endif
              endif
              condsum=condsum+condnow
              xcond=condnow*(b(i,1)*b(j,1)+b(i,2)*b(j,2)+b(i,3)*b(j,3))*volume &
                  +rhoc*ht(i)*(velox*b(j,1)+veloy*b(j,2)+veloz*b(j,3))*volume
              xmass=h(i)*rhoc*h(j)*volume*dynamic
              xcondp=condnow*(bp(i,1)*bp(j,1)+bp(i,2)*bp(j,2)+bp(i,3)*bp(j,3))*volumep &
                  +rhoc*ht(i)*(velox*bp(j,1)+veloy*bp(j,2)+veloz*bp(j,3))*volumep
              xmassp=h(i)*rhoc*h(j)*volumep*dynamic
              ael(i,j)=ael(i,j)+(xmass+dt*alpha*xcond)*w
              aelp(i,j)=aelp(i,j)+(xmass-dt*(1.-alpha)*xcondp)*w
            enddo
          enddo

        enddo

! Store heat production and conductivity in material property arrays
        eheat(elenum)=heatsum/mpe**2
        econd(elenum)=condsum/mpe**3

        do i=1,mpe
          do j=1,mpe
            ! Modified to allow transient surface temperatures (original version below)
            ! dwhipp (09/07)
            !bel(i)=bel(i)+aelp(i,j)*tempp(icon(j))
            if (kfix(icon(j)).ne.2) tprev=tempp(icon(j))
            if (kfix(icon(j)).eq.2) tprev=(-zp(j)+zl)*tlapse+tmsl
            bel(i)=bel(i)+aelp(i,j)*tprev
          enddo
        enddo

        do i=1,mpe
          ! Modified to allow transient surface temperatures (original version below)
          ! dwhipp (09/07)
          !if (kfix(icon(i)).eq.1) then
          if (kfix(icon(i)).ne.0) then
            ! NEW - set surface temperature to new value from lapse rate (original below)
            ! dwhipp (09/07)
            if (kfix(icon(i)).eq.1) fixt=tempp(icon(i))
            if (kfix(icon(i)).eq.2) fixt=(-zp(i)+zl)*tlapse+tmsl
            do j=1,mpe
              bel(j)=bel(j)-ael(j,i)*fixt
              ael(i,j)=0.
              ael(j,i)=0.
            enddo
          ael(i,i)=1.
          bel(i)=fixt
          endif
        enddo

      deallocate (rr,ss,tt,ww)
      deallocate (aelp,b,bp,vv)
      deallocate (x,y,z,xp,yp,zp,rb)
      deallocate (h,ht,dhdr,dhds,dhdt)

      return
      end
