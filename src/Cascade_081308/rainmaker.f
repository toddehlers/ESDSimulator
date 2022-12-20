c-----------------------------------------------------------------
      subroutine rainmaker(nnode,x,y,h,water,surface,
     &           prec,iflag_oro,iflag_uni)
c-----------------------------------------------------------------
c
c NOTE SOME USER DEFINED INPUT IS REQUIRED IN THE TOP 64 LINES
c OF THIS SUBROUTINE 
c This routine was written by Gerard Roe and implemented by Todd
c Todd Ehlers, 9/2001, Univ. Washington, Quaternary Research Center
c
      implicit none
      integer nnode
      integer iflag_oro
      integer iflag_uni
      integer i,j

c domain sizes [km]
      real sidex, sidey
      parameter(sidex = 400.)
      parameter(sidey = 150.)

c grid size [km]
      real del_gr
      parameter(del_gr = 2)

c regular grid dimension
      integer nxs, nys
      parameter (nxs = sidex/del_gr + 1)
      parameter (nys = sidey/del_gr + 1)

c Precipitation Parameters
c     T0 = climate temperature (default is 0C)
c     a0 [m/yr] &  a1 [m/yr per (m/s)]
c     alf =  reciprocal vertical velocity variance [1/(m/s)].
c     wnd = wind speed [m/s], 
c     angle = clockwise angle rel. to y=0 line [degrees]
      real T0,a0,a1,alf,wnd,angle
      parameter(T0 = 0.0)
c Note no need to change a0 or a1.
c      change a0 is same as changing To (approx)
c      changing a1 is approx the same as changing wnd
      parameter(a0=1.0, a1=110.0)
      parameter(alf = 1./0.01)
      parameter(wnd=1.5)
      parameter(angle = 90.)

c smoothing scales [km].
c xwind_s = cross wind smoothing scale
c upwnd_s = upwind smoothing scale
c if less than grid scale, NO smoothing is applied in that direction
      real xwnd_s,upwnd_s
      parameter(xwnd_s = 50.0, upwnd_s = 50.0)

c subsampling to speed up precip
      integer subsam
      parameter(subsam=1)

      real x_gr(nxs), y_gr(nys), prec_gr(nxs,nys)
      real prec(nnode),water(nnode),surface(nnode)

      real x(nnode), y(nnode),h(nnode)

c Save prec_gr so it can be used when iflag_oro == 0
      save prec_gr

      call oro_prec(nnode,x,y,h,surface,prec,water,
     &     prec_gr,x_gr,y_gr, del_gr,
     &     xwnd_s,upwnd_s,subsam,T0,
     &     a0,a1,alf,wnd,angle,
     &     iflag_oro,iflag_uni,nxs,nys)

      end

c*****************************************************************
c NO USER DEFINED INPUT NEEDED BELOW HERE
c*****************************************************************
c-----------------------------------------------------------------
      subroutine oro_prec(nnode,x,y,h,surface,prec,water,
     &     prec_gr,x_gr,y_gr, del_gr,
     &     xwnd_s,upwnd_s,subsam,T0,
     &     a0,a1,alf,wnd, angle,
     &     iflag_oro,iflag_uni,nxs,nys)
c-----------------------------------------------------------------
c code to produce precip param that works for CASCADE...
      implicit none
c     inputs: nnode             no. of nodes
c              x,y,h             x,y, and h of nodes
c              surface           surface area assoc. with each node
c              del_gr            grid scale for prec. calc.
c              nxs,nys           no xs and ys 
c              xwnd_s,upwnd_s    cross winds and upwind smooth. scales
c              subsam            node subsampling interval to 
c                                speed up interpolation 
c              T0,a0,a1,alf         precip parameters.
c              wnd,angle         wnd strength (m/s) and direction.
c                                angle = clockwise angle 
c                                from x=0 line.
c              iflag_oro         if = 1, recalc. precip.
c                                if = 0, just interpolate onto nodes.
c              iflag_uni         if = 1, uniform precip.
c                                if = 0, full calculation.
c outputs:
c              prec              precipitation rate at each node
c              water             precip *area at each node
c              x_gr,y_gr,prec_gr x,y,precip on reg. grid


c grid size [km]
      real del_gr
      real xwnd_s,upwnd_s

c subsampling to speed up precip
      integer subsam

      integer nnode
      integer iflag_oro, iflag_uni

c Precip parameters
c wind speed [m/s]
      real T0,a0,a1,alf
      real wnd,angle

      real x(nnode), y(nnode),h(nnode)
      real prec(nnode),water(nnode),surface(nnode)

c declare variables & define parameters
      integer i,j,k
      real dum

c regular grid dimension
      integer nxs, nys

c field variables
      real x_gr(nxs), y_gr(nys), prec_gr(nxs,nys)
      real z(nxs,nys), divQ(nxs,nys)
      real dzdx(nxs,nys), dzdy(nxs,nys)

      real bckgr

c begin routine.

c get background precip rate for h=w=gradh=0
c but including prob. dist. of vert velo. given by alf param.
         call back_div(T0,a0,a1,alf,bckgr)

c if uniform precip requested, set it equal to the 
c background rate and (sorry) go to bottom of routine.

c      print *,'Number of nodes according to rainmaker.f',nnode

      if (iflag_uni.eq.1) then
         do i = 1,nnode
c            print *,'i = ',i
c            print *,'nnode = ',nnode
c            prec(i) = 1.
            prec(i) = bckgr
c            print *,'prec = ',prec(i)
         enddo
         goto 100
      endif


c if update of precip on regular grid required
      if (iflag_oro.eq.1) then

c define regular grid.
      do i  = 1,nxs
         x_gr(i) = (i-1)*del_gr
      enddo
      do j  = 1,nys
         y_gr(j) = (j-1)*del_gr
      enddo

c Interpolate from nodal points to regular grid.
c uses bivar.f in this case.
      call grid(nnode,x,y,h,x_gr,y_gr,nxs,nys,subsam,z)
      
c output topography
c      open(unit=21, file = 'bivar.dat', status = 'new')
c      do i = 1,nxs
c         write(21,1000) (z(i,j),j=1,nys)
c      enddo
c      close(21)

c calculate slopes perpendicular to prevailing wind direction.
c normal to front in this case.
c note have to shift horiz. units to m
c      call xderiv(dAdx,A,delx,imx,jmx)
      call xderiv(dzdx,z,1000*del_gr,nxs,nys)
      call yderiv(dzdy,z,1000*del_gr,nxs,nys)

c output moist. conv.
c      open(unit=21, file = 'grady.dat', status = 'new')
c      do i = 1,nxs
c         write(21,1000) (dzdy(i,j),j=1,nys)
c      enddo
c      close(21)
c      open(unit=21, file = 'gradx.dat', status = 'new')
c      do i = 1,nxs
c         write(21,1000) (dzdx(i,j),j=1,nys)
c      enddo
c      close(21)

c calculate moisture convergence.
c use standard routine form my thesis.
      call moist_div(z,dzdy,dzdx,nxs,nys,wnd,angle,T0,
     &     a0,a1,alf,divq)

c output moist conv.
c      open(unit=21, file = 'divq.dat', status = 'new')
c      do i = 1,nxs
c         write(21,1000) (divq(i,j),j=1,nys)
c      enddo
c      close(21)

c average divQ over defined smoothing scales.
c gaussian upwind weighting applied in this case.
      call smooth(divq,x_gr,y_gr,nxs,nys,xwnd_s,upwnd_s,
     &     del_gr,angle,bckgr,prec_gr)
c end conditional
      endif

c do linear interpolation back onto nodes
      call linear(nnode,prec_gr,x,y,nxs,nys,del_gr,prec)

c multiply precip rate by nodal area to get rate of water input.
 100   continue
       do i =1,nnode
         water(i) = surface(i)*prec(i)
c         water(i) = surface(i)*1.0
      enddo


c 1000 format(500(E15.6))
c er, stop, end
      return 
      end

c---------------------------------------------------------------
      subroutine grid(nnode,x,y,h,x_gr,y_gr,nxs,nys,subsam,z)
c---------------------------------------------------------------
c routine to call bivar and to return a sampled version of the input 
c topography on a regular grid.
      implicit none
c     inputs: 
c             nnode     number of nodes 
c             x,y        co-ords of nodes 
c             h          elevation of nodes 
c             x_gr, y_gr co-ords of reg grid 
c             nxs, nys   kinda obvious 
c             subsam   subsampling of input nodes.
c
c      outputs: 
c             z       elevation on regular grid 
c

c declare variables.
      integer i,j,k,kmax
      integer nxs,nys
      integer nnode
      integer subsam
      real x(nnode), y(nnode), h(nnode)
      real x1(nnode/subsam),y1(nnode/subsam), h1(nnode/subsam)
      real x_gr(nxs), y_gr(nys), z(nxs,nys)

c work arrays for bivar      
      integer IWK(31*nnode+nxs*nys)
      real     WK(6*nnode)

c subsample initial data (to speed up interpolation).
      kmax = int(nnode/subsam)
      do k = 1,kmax
         x1(k) = x(k*subsam)
         y1(k) = y(k*subsam)
         h1(k) = h(k*subsam)
      enddo
      
c call idsfft to get smoothed fields.
       call IDSFFT (1,kmax,x1,y1,h1,nxs,nys,nxs,x_gr,y_gr,z,IWK,WK)

c er, return, end
      return
      end
c----------------------------------------------------------------
      subroutine moist_div(z,dzdy,dzdx,nxs,nys,wnd,angle,Tbck,a0,a1,
     &     alf,divq)
c----------------------------------------------------------------
c code to calculate moisture divergence from the w and T fields passed 
c calculated in the main program.
      implicit none
c inputs:  h, dhdy/dhdx   height, gradient in x and y directions
c                            (calc. on reg. grid)
c          nxs,nys       number of grid points.
c          Tbck          background temperature 
c          a0,a1         background and slope term, respectively.
c          wnd           prescribed prevailing wind speed.
c          angle         prescribed wind direction
c          alf           prescribed variance in atmospheric vert. velocity.
c
c outputs: divq          moisture divergence at each point.
c
c other:   a,b,e0        parameters in Clasius-Clapeyron relation.
c          p,a1,a2,a3    parameters in error function approximation.
c          gamma         atmospheric lapse rate. 

c Other variables etc.
      integer i,j,nxs,nys
      real divq(nxs,nys),w,T,Tbck,wnd, angle
      real z(nxs,nys), dzdy(nxs,nys), dzdx(nxs,nys)

c Precipitation parameters taken from my thesis.
      real a0,a1
      real a,b,e0,satfac
      real p,c1,c2,c3

c Other stuff for precipitation.
      real alf,alf_p
      real x_0
      real erf,t0
      real I1,I2

      real pi
      parameter(pi=3.1415927)

c lapse rate
      real gamma
      parameter(gamma = -6.5E-03)
c Sat. vap. pres. params...  
c and error function parameters for esat distribution.
      parameter(a=17.67,b=243.5,e0=611.2)
      parameter(p=0.47047,c1=0.3480242,c2=-0.0958798,c3=0.7478556)

c loop over grid
      do i = 1,nxs
         do j = 1,nys

c Probability distribution of velocities around the time-mean
c upslope velocity integrated over time to give an effective upslope
c velocity. See my thesis...

c transform for precip.
           w = wnd*sin(3.14159*angle/180.)*dzdy(i,j)
     &        - wnd*cos(3.14159*angle/180.)*dzdx(i,j)
           alf_p = alf/a1
           x_0   = a1*w + a0

c analytic expression for precip. integral.
           t0=1/(1+p*alf_p*abs(x_0))
           erf=1-(c1*t0+c2*t0**2+c3*t0**3)*exp(-alf_p**2*x_0**2)
           I1=x_0/2*(1+sign(1.,x_0)*erf)
           I2=1./(2*alf_p*sqrt(pi))*exp(-alf_p**2*x_0**2)
c           write(6,*) alf_p,erf,x_0,I1,I2
          
c divq=esat*(I1+I2)
c change temperatue by changing 0.0 here.
           T = Tbck+gamma*z(i,j)
c sat. vap. factor
           satfac=exp(a*T/(b+T))
           divq(i,j)=satfac*(I1+I2)
        enddo
      enddo
c er, return, end
      return
      end

c-----------------------------------------
      subroutine back_div(Tbck,a0,a1,alf,bckgr)
c-----------------------------------------
c code to calculate background moisture divergence (ok, convergence) 
c at z=0 in absence of surface slopes. Uses same prob. distribution of
c vertical velocities as divq
       implicit none
c          a0,a1         background and slope term, respectively.
c          alf           prescribed variance in atmospheric vert. velocity.
c          Tbck             background temperature.
c
c outputs: bckgr          moisture divergence at each point.
c
c other:   a,b,e0        parameters in Clasius-Clapeyron relation.
c          p,c1,c2,c3    parameters in error function approximation.
c          gamma         atmospheric lapse rate. 

c Other variables etc.
      integer i,j,nxs,nys
      real bckgr,Tbck

c Precipitation parameters taken from my thesis.
      real a0,a1
      real a,b,e0,satfac
      real p,c1,c2,c3

c Other stuff for precipitation.
      real alf,alf_p
      real x_0
      real erf,t0
      real I1,I2

      real pi
      parameter(pi=3.1415927)

c lapse rate.
      real gamma
      parameter(gamma = -6.5E-03)

c Sat. vap. pres. params...  
c and error function parameters for esat distribution.
      parameter(a=17.67,b=243.5,e0=611.2)
      parameter(p=0.47047,c1=0.3480242,c2=-0.0958798,c3=0.7478556)

c transform for precip.
      alf_p = alf/a1
      x_0   = a0

c analytic expression for precip. integral.
      t0=1/(1+p*alf_p*abs(x_0))
      erf=1-(c1*t0+c2*t0**2+c3*t0**3)*exp(-alf_p**2*x_0**2)
      I1=x_0/2*(1+sign(1.,x_0)*erf)
      I2=1./(2*alf_p*sqrt(pi))*exp(-alf_p**2*x_0**2)
c divq=esat*(I1+I2)

c sat. vap. factor
      satfac=exp((a*Tbck)/(b+Tbck))
      bckgr=satfac*(I1+I2)
c er, return, end
      return
      end

c------------------------------------------------------
      subroutine smooth(divq,x_gr,y_gr,nxs,nys,xwnd_s,upwnd_s,
     &     del_gr,angle,bckgr,prec_gr)
c------------------------------------------------------
      implicit none
c subroutine to smooth moisture divergence field to get
c precipitation field
c
c     inputs:
c             divq           moisture conv. on reg. grid
c             x_gr,y_gr            reg. grid
c             nxs,nys        kinda obvious
c             xwnd_s, upwnd_s assumed smoothing scales across 
c                            and along the prevailing wind directions
c             del_gr         grid size
c             bckgr          background precip rate in [m/yr]
c             angle          prescribed angle of prevailing wind
c
c     outputs:
c             prec_gr        also kinda obvious

c others:
c             phi =          distance from gd. pt. perp. to 
c                            prevailing wind direction
c             lam =          distance from gd. pt. parallel. to 
c                            prevailing wind direction

c delare variables
      integer i,j,i1,j1, ixwnd_s, iupwnd_s
      integer nxs,nys
      integer iwk1,iwk2,jwk1,jwk2
      real xwnd_s, upwnd_s, del_gr, a0
      real divq(nxs,nys),prec_gr(nxs,nys)
      real x_gr(nxs,nys),y_gr(nxs,nys)
      real wt1,wt2, sum_wt
      real phi,lam
      real bckgr,angle

      real pi
      parameter(pi = 3.1415926)

c smoothing scales in terms of number grid points
      ixwnd_s = int(xwnd_s/del_gr)
      iupwnd_s = int(upwnd_s/del_gr)

c Determine boundaries of smoothing region.
c Ensure at least 2 smoothing widths upwind and
c 2 smoothing widths either side.
      jwk1 = 2*int(ixwnd_s*cos(angle*pi/180.) 
     &     + iupwnd_s*sin(angle*pi/180.))
      jwk2 = 2*int(ixwnd_s*cos(angle*pi/180.))
      if (angle.le.90) then
         iwk1 = 2*int(ixwnd_s*sin(angle*pi/180.))
         iwk2 = 2*int(ixwnd_s*sin(angle*pi/180.)+
     &                iupwnd_s*cos(angle*pi/180.))
      else
         iwk1 = 2*int(ixwnd_s*sin(angle*pi/180.)-
     &        iupwnd_s*cos(angle*pi/180.))
         iwk2 = 2*int(ixwnd_s*sin(angle*pi/180.))
      endif
c      write(6,*) iwk1,iwk2,jwk1,jwk2
    
c begin loop over variables.
      do i = 1,nxs
         do j = 1,nys
            prec_gr(i,j) = 0
            sum_wt = 0

            do j1 = j-jwk1,j+jwk2
               do i1 = i-iwk1, i+iwk2
c phi = cross wind distance
c lam  = upwind distance
                  phi=(i1-i)*sin(angle*pi/180.)
     &                 +(j1-j)*cos(angle*pi/180.)
                  lam=-(i1-i)*cos(angle*pi/180.) 
     &                 +(j1-j)*sin(angle*pi/180.)

c define weighting factors, allow for zero smoothing.
c if smoothing scale is less than xwnd or upwind then only
c count points within 1 grid pt of the wind direction.
                  if (ixwnd_s.eq.0) then
                     if(abs(phi).lt.1.0) then
                        wt1=1.
                     else
                        wt1=0.   
                     endif
                  else
                     wt1 = exp(-1.*(phi/ixwnd_s)**2) 
                  endif

                  if (iupwnd_s.eq.0) then
                     if(abs(lam).lt.1.0) then
                        wt2=1.
                     else
                        wt2=0.  
                     endif
                  else
                     wt2 = exp(-1.*(lam/iupwnd_s)**2) 
                  endif


                  if ((i1.ge.1).and.(i1.le.nxs)) then
                     if (j1.lt.1) then
                        prec_gr(i,j) = prec_gr(i,j) + bckgr*wt1*wt2
                        sum_wt = sum_wt+wt1*wt2          
                     elseif (j1.le.nys) then
                        prec_gr(i,j) = prec_gr(i,j) + 
     &                       divq(i1,j1)*wt1*wt2
                        sum_wt = sum_wt+wt1*wt2
                     endif  
                  endif
               enddo
            enddo            

c finally divide by weighting factor.
            prec_gr(i,j) = prec_gr(i,j)/sum_wt

         enddo
      enddo

c, er return, end.
      return
      end


c-----------------------------------------------------------------------
      subroutine linear(nnode,prec_gr,x,y,nxs,nys,del_gr,prec)
c-----------------------------------------------------------------------
      implicit none
c subroutine to do linear interpolation from regular grid back onto nodes
c     inputs: 
c            nnode             number of nodes
c            prec_gr           precipitation on regular grid
c            x,y               coords at nodes
c            del_gr            regular grid spacing
c     outputs:
c            prec              precipitation on nodes

      integer i,nnode,nxs,nys,i0,j0
      real del_gr
      real prec_gr(nxs,nys)
      real x(nnode),y(nnode),prec(nnode)
      real iwk,jwk,t,u

c begin loop over nodes
      do i = 1,nnode
         iwk = x(i)/del_gr + 1
         jwk = y(i)/del_gr + 1
         i0 = int(x(i)/del_gr)+1
         j0 = int(y(i)/del_gr)+1
         if ((iwk.lt.1) .or. (iwk.gt.nxs) .or. 
     &        (jwk.lt.1) .or. (jwk.gt.nys)) then
            write(6,*) 'oops','iwk,jwk,i',iwk,jwk,i,x(i),y(i)
            stop
         endif
         
c working numbers.
         t = (iwk-i0)
         u = (jwk-j0)
c linear interpolation.
         prec(i) = (1-t)*(1-u)*prec_gr(i0,j0) +
     &              t*(1-u)*   prec_gr(i0+1,j0) +
     &              t*u*       prec_gr(i0+1,j0+1) +
     &              (1-t)*u*   prec_gr(i0,j0+1)
      enddo

c er, return, end
      return
      end

c------------------------------------------------------------
      subroutine xderiv(dAdx,A,delx,imx,jmx)
c-------------------------------------------------------------
      implicit none
c subroutine to take a real 2-d matrix and take a deriv in the 
c x direction or the first array address. 2nd order accurate.
      integer i,j,imx,jmx
      real A(imx,jmx), dAdx(imx,jmx)
      real delx

      do j = 1,jmx
         dAdx(1,j) = (-3*A(1,j)+4*A(2,j)-A(3,j))/(2*delx)
         do i = 2,imx-1
            dAdx(i,j) = (A(i+1,j)-A(i-1,j))/(2*delx)
         enddo
         dAdx(imx,j) = (3*A(imx,j)-4*A(imx-1,j)+A(imx-2,j))/(2*delx)
      enddo

      return
      end

c------------------------------------------------------------
      subroutine yderiv(dAdy,A,dely,imx,jmx)
c-------------------------------------------------------------
      implicit none
c subroutine to take a real 2-d matrix and take a deriv in the 
c y direction or the second array address. 2nd order accurate.
      integer i,j,imx,jmx
      real A(imx,jmx), dAdy(imx,jmx)
      real dely
      
      do i = 1,imx
         dAdy(i,1) = (-3*A(i,1)+4*A(i,2)-A(i,3))/(2*dely)
         do j = 2,jmx-1
            dAdy(i,j) = (A(i,j+1)-A(i,j-1))/(2*dely)
         enddo
         dAdy(i,jmx) = (3*A(i,jmx)-4*A(i,jmx-1)+A(i,jmx-2))/(2*dely)
      enddo

      return
      end

C
C $Id: idsfft.f,v 1.5 2000/08/22 15:06:54 haley Exp $
C                                                                      
C                Copyright (C)  2000
C        University Corporation for Atmospheric Research
C                All Rights Reserved
C
C This file is free software; you can redistribute it and/or modify
C it under the terms of the GNU General Public License as published
C by the Free Software Foundation; either version 2 of the License, or
C (at your option) any later version.
C
C This software is distributed in the hope that it will be useful, but
C WITHOUT ANY WARRANTY; without even the implied warranty of
C MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
C General Public License for more details.
C
C You should have received a copy of the GNU General Public License
C along with this software; if not, write to the Free Software
C Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
C USA.
C
      SUBROUTINE IDSFFT (MD,NDP,XD,YD,ZD,NXI,NYI,NZI,XI,YI,ZI,IWK,WK)
C
C DIMENSION OF           XD(NDP), YD(NDP), ZD(NDP), XI(NXI),
C ARGUMENTS              YI(NYI), ZI(NZI,NYI), WK(6*NDP),
C                        IWK(31*NDP + NXI*NYI)
C
C PURPOSE                This subroutine performs smooth surface
C                        fitting when the projections of the data
C                        points in the X-Y plane are irregularly
C                        distributed in the plane.
C
C USAGE                  CALL IDSFFT (MD,NDP,XD,YD,ZD,NXI,NYI,NZI,
C                                     XI,YI,ZI,IWK,WK)
C
C ARGUMENTS
C
C ON INPUT               MD
C                          Mode of computation (must be 1, 2, or 3,
C                          else an error return will occur).
C                           = 1 if this is the first call to this
C                               subroutine, or if the value of NDP
C                               has been changed from the previous
C                               call, or if the contents of the XD
C                               or YD arrays have been changed from
C                               the previous call.
C                           = 2 if the values of NDP and the XD and
C                               YD arrays are unchanged from the
C                               previous call, but new values for
C                               XI and YI are being used.  If MD = 2
C                               and NDP has been changed since the
C                               previous call to IDSFFT, an error
C                               return occurs.
C                           = 3 if the values of NDP, NXI, NYI, XD,
C                               YD, XI, and YI are unchanged from the
C                               previous call, i.e., if the only change
C                               on input to IDSFFT is in the ZD array.
C                               If MD = 3 and NDP, NXI, or NYI has been
C                               changed since the previous call to
C                               IDSFFT, an error return occurs.
C
C                           Between the call with MD = 2 or MD = 3 and
C                           the preceding call, the IWK and WK work
C                           arrays should not be disturbed.
C
C                        NDP
C                          Number of data points (must be 4 or
C                          greater, else an error return will occur).
C
C                        XD
C                          Array of dimension NDP containing the X
C                          coordinates of the data points.
C
C                        YD
C                          Array of dimension NDP containing the Y
C                          coordinates of the data points.
C
C                        ZD
C                          Array of dimension NDP containing the Z
C                          coordinates of the data points.
C
C                        NXI
C                          Number of output grid points in the X
C                          direction (must be 1 or greater, else
C                          an error return will occur).
C
C                        NYI
C                          Number of output grid points in the Y
C                          direction (must be 1 or greater, else
C                          an error return will occur).
C
C                        NZI
C                          First dimension of ZI as declared in the
C                          calling program.  NZI must be greater than
C                          or equal to NXI, else an error return will
C                          occur.
C
C                        XI
C                          Array of dimension NXI containing the
C                          X coordinates of the output grid points.
C
C                        YI
C                         Array of dimension NYI containing the
C                         Y coordinates of the output grid points.
C
C                        IWK
C                          Integer work array of dimension at
C                          least 31*NDP + NXI*NYI.
C
C                        WK
C                          Real work array of dimension at least 6*NDP.
C
C ON OUTPUT              ZI
C                           Real, two-dimensional array of dimension
C                           (NZI,NYI), storing the interpolated Z
C                           values at the output grid points.
C
C SPECIAL CONDITIONS     Inadequate work space IWK and WK may cause
C                        incorrect results.
C
C                        The data points must be distinct and their
C                        projections in the X-Y plane must not be
C                        collinear, else an error return occurs.
C
C IDSFFT calls the subroutines IDGRID, IDPDRV, IDPTIP, and IDTANG.
C
C Declaration statements.
C
      DIMENSION XD(NDP),YD(NDP),ZD(NDP),XI(NXI),YI(NYI),ZI(NZI,NYI),
     +          IWK(31*NDP+NXI*NYI),WK(6*NDP)
C
      COMMON /IDPT/ ITPV,DMMY(27)
      SAVE   /IDPT/
C
      COMMON /IDCOMN/ INTY,ITTY,ALSP,BLSP,CLSP,XAVG,YAVG
      SAVE   /IDCOMN/
C
C Check for an uncleared prior error.
C
        IF (ICFELL('IDSFFT (BIVAR) - UNCLEARED PRIOR ERROR',1).NE.0)
     +                                                        RETURN
C
C Check for input errors.
C
      IF (MD.LT.1.OR.MD.GT.3) THEN
        CALL SETER ('IDSFFT (BIVAR) - INPUT VARIABLE MD IS OUT OF RANGE'
     +,2,1)
        RETURN
      END IF
C
      IF (NDP.LT.4) THEN
        CALL SETER ('IDSFFT (BIVAR) - INPUT VARIABLE NDP IS OUT OF RANGE
     +',3,1)
        RETURN
      END IF
C
      IF (NXI.LT.1.OR.NYI.LT.1) THEN
        CALL SETER ('IDSFFT (BIVAR) - INPUT VARIABLE NXI OR NYI IS OUT O
     +F RANGE',4,1)
        RETURN
      END IF
C
      IF (NXI.GT.NZI) THEN
        CALL SETER ('IDSFFT (BIVAR) - INPUT VARIABLE NZI IS LESS THAN NX
     +I',5,1)
        RETURN
      END IF
C
      IF (MD.LE.1) THEN
C
        IWK(1)=NDP
C
      ELSE
C
        NDPPV=IWK(1)
C
        IF (NDP.NE.NDPPV) THEN
          CALL SETER  ('IDSFFT (BIVAR) - MD = 2 OR 3 BUT NDP WAS CHANGED
     + SINCE LAST CALL',6,1)
          RETURN
        END IF
C
      END IF
C
      IF (MD.LE.2) THEN
C
        IWK(2)=INTY
        IWK(3)=NXI
        IWK(4)=NYI
C
      ELSE
C
        INTYPV=IWK(2)
C
        IF (INTY.NE.INTYPV) THEN
          CALL SETER ('IDSFFT (BIVAR) - MD = 3 BUT ITY WAS CHANGED SINCE
     + LAST CALL',7,1)
           RETURN
        END IF
C
        NXIPV=IWK(3)
C
        IF (NXI.NE.NXIPV) THEN
          CALL SETER ('IDSFFT (BIVAR) - MD = 3 BUT NXI WAS CHANGED SINCE
     + LAST CALL',8,1)
           RETURN
        END IF
C
        NYIPV=IWK(4)
C
        IF (NYI.NE.NYIPV) THEN
          CALL SETER ('IDSFFT (BIVAR) - MD = 3 BUT NYI WAS CHANGED SINCE
     + LAST CALL',9,1)
          RETURN
        END IF
C
      END IF
C
C Allocate storage areas in the array IWK.
C
      JWIPT=16
      JWIWL=6*NDP+1
      JWNGP0=JWIWL-1
      JWIPL=24*NDP+1
      JWIWP=30*NDP+1
      JWIGP0=31*NDP
      JWWPD=5*NDP+1
C
C If MD = 1, triangulate the X/Y plane.  (If MD = 2 or 3, this has
C already been done and need not been redone.)
C
      IF (MD.EQ.1) THEN
        CALL IDTANG (NDP,XD,YD,NT,IWK(JWIPT),NL,IWK(JWIPL),IWK(JWIWL),
     +                                                  IWK(JWIWP),WK)
        IF (ICFELL('IDSFFT',10).NE.0) RETURN
        IWK(5)=NT
        IWK(6)=NL
      ELSE
        NT=IWK(5)
        NL=IWK(6)
      END IF
C
      IF (NT.EQ.0) RETURN
C
C If linear interpolation is activated, compute the coefficients of the
C least-squares-fit plane and the mean values of X and Y.
C
      IF (INTY.NE.0) THEN
        CALL IDLSQF (XD,YD,ZD,NDP,ALSP,BLSP,CLSP,XAVG,YAVG)
      END IF
C
C If MD = 1 or 2, sort output grid points in ascending order by the
C number of the containing region (a triangle number or a combined
C pair of border line segment numbers).  (If MD = 3, this has already
C been done and need not be redone.)
C
      IF (MD.LE.2) THEN
        CALL IDGRID (XD,YD,NT,IWK(JWIPT),NL,IWK(JWIPL),NXI,
     +               NYI,XI,YI,IWK(JWNGP0+1),IWK(JWIGP0+1))
        IF (ICFELL('IDSFFT',11).NE.0) RETURN
      END IF
C
C If quintic interpolation is activated, estimate partial derivatives.
C
      IF (INTY.EQ.0) THEN
        CALL IDPDRV (NDP,XD,YD,ZD,NT,IWK(JWIPT),WK,WK(JWWPD))
      END IF
C
C Interpolate to get ZI values.
C
      ITPV=0
      JIG0MX=0
      JIG1MN=NXI*NYI+1
      NNGP=NT+2*NL
C
      DO 105 JNGP=1,NNGP
C
        ITI=JNGP
        IF (JNGP.LE.NT) GO TO 101
        IL1=(JNGP-NT+1)/2
        IL2=(JNGP-NT+2)/2
        IF (IL2.GT.NL) IL2=1
        ITI=IL1*(NT+NL)+IL2
C
  101   JWNGP=JWNGP0+JNGP
        NGP0=IWK(JWNGP)
        IF (NGP0.EQ.0) GO TO 103
        JIG0MN=JIG0MX+1
        JIG0MX=JIG0MX+NGP0
C
        DO 102 JIGP=JIG0MN,JIG0MX
          JWIGP=JWIGP0+JIGP
          IZI=IWK(JWIGP)
          IYI=(IZI-1)/NXI+1
          IXI=IZI-NXI*(IYI-1)
          CALL IDPTIP(XD,YD,ZD,NT,IWK(JWIPT),NL,IWK(JWIPL),WK,
     +                        ITI,XI(IXI),YI(IYI),ZI(IXI,IYI))
  102   CONTINUE
C
  103   JWNGP=JWNGP0+2*NNGP+1-JNGP
        NGP1=IWK(JWNGP)
        IF (NGP1.EQ.0) GO TO 105
        JIG1MX=JIG1MN-1
        JIG1MN=JIG1MN-NGP1
C
        DO 104 JIGP=JIG1MN,JIG1MX
          JWIGP=JWIGP0+JIGP
          IZI=IWK(JWIGP)
          IYI=(IZI-1)/NXI+1
          IXI=IZI-NXI*(IYI-1)
          CALL IDPTIP(XD,YD,ZD,NT,IWK(JWIPT),NL,IWK(JWIPL),WK,
     +                        ITI,XI(IXI),YI(IYI),ZI(IXI,IYI))
  104   CONTINUE
C
  105 CONTINUE
C
C Done.
C
        RETURN
C
      END
C
C $Id: idgrid.f,v 1.4 2000/08/22 15:06:53 haley Exp $
C                                                                      
C                Copyright (C)  2000
C        University Corporation for Atmospheric Research
C                All Rights Reserved
C
C This file is free software; you can redistribute it and/or modify
C it under the terms of the GNU General Public License as published
C by the Free Software Foundation; either version 2 of the License, or
C (at your option) any later version.
C
C This software is distributed in the hope that it will be useful, but
C WITHOUT ANY WARRANTY; without even the implied warranty of
C MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
C General Public License for more details.
C
C You should have received a copy of the GNU General Public License
C along with this software; if not, write to the Free Software
C Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
C USA.
C
      SUBROUTINE IDGRID (XD,YD,NT,IPT,NL,IPL,NXI,NYI,XI,YI,NGP,IGP)
C
C This subroutine organizes grid points for surface fitting by
C sorting them in ascending order of triangle number or border
C line segment number.
C
C The input arguments are as follows:
C
C   XD,YD  arrays of dimension NDP containing the X and Y
C          coordinates of the data points, where NDP is the
C          number of the data points.
C
C   NT     number of triangles.
C
C   IPT    integer array of dimension 3*NT containing the
C          point numbers of the vertices of the triangles.
C
C   NL     number of border line segments.
C
C   IPL    integer array of dimension 3*NL containing the
C          point numbers of the end points of the border
C          line segments and their respective triangle
C          numbers.
C
C   NXI    number of grid points in the X coordinate.
C
C   NYI    number of grid points in the Y coordinate.
C
C   XI,YI  arrays of dimension NXI and NYI containing
C          the X and Y coordinates of the grid points,
C          respectively.
C
C The output arguments are as follows:
C
C   NGP    integer array of dimension 2*(NT+2*NL) where the
C          number of grid points that belong to each of the
C          triangles or of the border line segments are to
C          be stored.
C
C   IGP    integer array of dimension NXI*NYI where the
C          grid point numbers are to be stored in ascending
C          order of the triangle number or border line
C          segment number.
C
C Declaration statements:
C
      DIMENSION XD(*),YD(*),IPT(*),IPL(*),XI(*),YI(*),NGP(*),IGP(*)
C
      COMMON /IDCOMN/ INTY,ITTY,ALSP,BLSP,CLSP,XAVG,YAVG
      SAVE   /IDCOMN/
C
C Statement functions.
C
      SPDT(U1,V1,U2,V2,U3,V3)=(U1-U2)*(U3-U2)+(V1-V2)*(V3-V2)
      VPDT(U1,V1,U2,V2,U3,V3)=(U1-U3)*(V2-V3)-(V1-V3)*(U2-U3)
C
C Preliminary processing.
C
      NXINYI=NXI*NYI
C
      XIMN=MIN(XI(1),XI(NXI))
      XIMX=MAX(XI(1),XI(NXI))
      YIMN=MIN(YI(1),YI(NYI))
      YIMX=MAX(YI(1),YI(NYI))
C
C Determine grid points inside the data area.
C
      JNGP0=0
      JNGP1=2*(NT+2*NL)+1
C
      JIGP0=0
      JIGP1=NXINYI+1
C
      DO 114 IT0=1,NT
C
        NGP0=0
        NGP1=0
C
        X1=XD(IPT(IT0*3-2))
        Y1=YD(IPT(IT0*3-2))
        X2=XD(IPT(IT0*3-1))
        Y2=YD(IPT(IT0*3-1))
        X3=XD(IPT(IT0*3))
        Y3=YD(IPT(IT0*3))
C
        XMN=MIN(X1,X2,X3)
        XMX=MAX(X1,X2,X3)
        YMN=MIN(Y1,Y2,Y3)
        YMX=MAX(Y1,Y2,Y3)
C
        INSD=0
C
        DO 102 IXI=1,NXI
          IF (XI(IXI).GE.XMN.AND.XI(IXI).LE.XMX) GO TO 101
          IF (INSD.EQ.0) GO TO 102
          IXIMX=IXI-1
          GO TO 103
  101     IF (INSD.EQ.1) GO TO 102
          INSD=1
          IXIMN=IXI
  102   CONTINUE
C
        IF (INSD.EQ.0) GO TO 113
C
        IXIMX=NXI
C
  103   DO 112 IYI=1,NYI
C
          YII=YI(IYI)
C
          IF (YII.LT.YMN.OR.YII.GT.YMX) GO TO 112
C
          DO 111  IXI=IXIMN,IXIMX
C
            XII=XI(IXI)
C
            L=0
            IF (VPDT(X1,Y1,X2,Y2,XII,YII)) 111,104,105
  104       L=1
  105       IF (VPDT(X2,Y2,X3,Y3,XII,YII)) 111,106,107
  106       L=1
  107       IF (VPDT(X3,Y3,X1,Y1,XII,YII)) 111,108,109
  108       L=1
  109       IZI=NXI*(IYI-1)+IXI
C
            IF (L.EQ.0) THEN
              NGP0=NGP0+1
              JIGP0=JIGP0+1
              IF (JIGP0.GE.JIGP1) THEN
                CALL SETER ('IDGRID (BIVAR) - INTERNAL ERROR - SEE CONSU
     +LTANT',1,1)
                RETURN
              END IF
              IGP(JIGP0)=IZI
            ELSE
              IF (JIGP1.LE.NXINYI) THEN
                DO 110 JIGP1I=JIGP1,NXINYI
                  IF (IZI.EQ.IGP(JIGP1I)) GO TO 111
  110           CONTINUE
              END IF
              NGP1=NGP1+1
              JIGP1=JIGP1-1
              IF (JIGP1.LE.JIGP0) THEN
                CALL SETER ('IDGRID (BIVAR) - INTERNAL ERROR - SEE CONSU
     +LTANT',2,1)
                RETURN
              END IF
              IGP(JIGP1)=IZI
            END IF
C
  111     CONTINUE
C
  112   CONTINUE
C
  113   JNGP0=JNGP0+1
        NGP(JNGP0)=NGP0
C
        JNGP1=JNGP1-1
        NGP(JNGP1)=NGP1
C
  114 CONTINUE
C
C If linear interpolation is selected, jump to do the semi-infinite
C quadrangular areas formed by rays emanating from (XAVG,YAVG).
C
      IF (INTY.NE.0) GO TO 300
C
C Determine grid points outside the data area in each of the
C semi-infinite rectangular and triangular areas.
C
      DO 263 IL0=1,NL
C
C Rectangular area ...
C
        NGP0=0
        NGP1=0
C
        X1=XD(IPL(IL0*3-2))
        Y1=YD(IPL(IL0*3-2))
        X2=XD(IPL(IL0*3-1))
        Y2=YD(IPL(IL0*3-1))
C
        XMN=XIMN
        XMX=XIMX
        YMN=YIMN
        YMX=YIMX
C
        IF (Y2.GE.Y1) XMN=MIN(X1,X2)
        IF (Y2.LE.Y1) XMX=MAX(X1,X2)
        IF (X2.LE.X1) YMN=MIN(Y1,Y2)
        IF (X2.GE.X1) YMX=MAX(Y1,Y2)
C
        INSD=0
C
        DO 202 IXI=1,NXI
          IF (XI(IXI).GE.XMN.AND.XI(IXI).LE.XMX) GO TO 201
          IF (INSD.EQ.0) GO TO 202
          IXIMX=IXI-1
          GO TO 203
  201     IF (INSD.EQ.1) GO TO 202
          INSD=1
          IXIMN=IXI
  202   CONTINUE
C
        IF (INSD.EQ.0) GO TO 213
C
        IXIMX=NXI
C
  203   DO 212 IYI=1,NYI
C
          YII=YI(IYI)
C
          IF (YII.LT.YMN.OR.YII.GT.YMX) GO TO 212
C
          DO 211 IXI=IXIMN,IXIMX
C
            XII=XI(IXI)
C
            L=0
            IF (VPDT(X1,Y1,X2,Y2,XII,YII)) 205,204,211
  204       L=1
  205       IF (SPDT(X2,Y2,X1,Y1,XII,YII)) 211,206,207
  206       L=1
  207       IF (SPDT(X1,Y1,X2,Y2,XII,YII)) 211,208,209
  208       L=1
  209       IZI=NXI*(IYI-1)+IXI
C
            IF (L.EQ.0) THEN
              NGP0=NGP0+1
              JIGP0=JIGP0+1
              IF (JIGP0.GE.JIGP1) THEN
                CALL SETER ('IDGRID (BIVAR) - INTERNAL ERROR - SEE CONSU
     +LTANT',3,1)
                RETURN
              END IF
              IGP(JIGP0)=IZI
            ELSE
              IF (JIGP1.LE.NXINYI) THEN
                DO 210 JIGP1I=JIGP1,NXINYI
                  IF (IZI.EQ.IGP(JIGP1I)) GO TO 211
  210           CONTINUE
              END IF
              NGP1=NGP1+1
              JIGP1=JIGP1-1
              IF (JIGP1.LE.JIGP0) THEN
                CALL SETER ('IDGRID (BIVAR) - INTERNAL ERROR - SEE CONSU
     +LTANT',4,1)
                RETURN
              END IF
              IGP(JIGP1)=IZI
            END IF
C
  211     CONTINUE
C
  212   CONTINUE
C
  213   JNGP0=JNGP0+1
        NGP(JNGP0)=NGP0
C
        JNGP1=JNGP1-1
        NGP(JNGP1)=NGP1
C
C Triangular area ...
C
        NGP0=0
        NGP1=0
C
        ILP1=MOD(IL0,NL)+1
C
        X3=XD(IPL(ILP1*3-1))
        Y3=YD(IPL(ILP1*3-1))
C
        XMN=XIMN
        XMX=XIMX
        YMN=YIMN
        YMX=YIMX
C
        IF (Y3.GE.Y2.AND.Y2.GE.Y1) XMN=X2
        IF (Y3.LE.Y2.AND.Y2.LE.Y1) XMX=X2
        IF (X3.LE.X2.AND.X2.LE.X1) YMN=Y2
        IF (X3.GE.X2.AND.X2.GE.X1) YMX=Y2
C
        INSD=0
C
        DO 252 IXI=1,NXI
          IF (XI(IXI).GE.XMN.AND.XI(IXI).LE.XMX) GO TO 251
          IF (INSD.EQ.0) GO TO 252
          IXIMX=IXI-1
          GO TO 253
  251     IF (INSD.EQ.1) GO TO 252
          INSD=1
          IXIMN=IXI
  252   CONTINUE
C
        IF (INSD.EQ.0) GO TO 261
C
        IXIMX=NXI
C
  253   DO 260 IYI=1,NYI
C
          YII=YI(IYI)
C
          IF (YII.LT.YMN.OR.YII.GT.YMX) GO TO 260
C
          DO 259 IXI=IXIMN,IXIMX
C
            XII=XI(IXI)
C
            L=0
            IF (SPDT(X1,Y1,X2,Y2,XII,YII)) 255,254,259
  254       L=1
  255       IF (SPDT(X3,Y3,X2,Y2,XII,YII)) 257,256,259
  256       L=1
  257       IZI=NXI*(IYI-1)+IXI
C
            IF (L.EQ.0) THEN
              NGP0=NGP0+1
              JIGP0=JIGP0+1
              IF (JIGP0.GE.JIGP1) THEN
                CALL SETER ('IDGRID (BIVAR) - INTERNAL ERROR - SEE CONSU
     +LTANT',5,1)
                RETURN
              END IF
              IGP(JIGP0)=IZI
            ELSE
              IF (JIGP1.LE.NXINYI) THEN
                DO 258 JIGP1I=JIGP1,NXINYI
                  IF (IZI.EQ.IGP(JIGP1I)) GO TO 259
  258           CONTINUE
              END IF
              NGP1=NGP1+1
              JIGP1=JIGP1-1
              IF (JIGP1.LE.JIGP0) THEN
                CALL SETER ('IDGRID (BIVAR) - INTERNAL ERROR - SEE CONSU
     +LTANT',6,1)
                RETURN
              END IF
              IGP(JIGP1)=IZI
            END IF
C
  259     CONTINUE
C
  260   CONTINUE
C
  261   JNGP0=JNGP0+1
        NGP(JNGP0)=NGP0
C
        JNGP1=JNGP1-1
        NGP(JNGP1)=NGP1
C
  263 CONTINUE
C
      RETURN
C
C Linear interpolation was selected.  Examine semi-infinite quadrangular
C areas formed by rays emanating from (XAVG,YAVG) and passing through
C the endpoints of a line segment on the edge of the data area.
C
  300 DO 314 IL0=1,NL
C
        NGP0=0
        NGP1=0
C
        X1=XD(IPL(IL0*3-2))
        Y1=YD(IPL(IL0*3-2))
        X2=XD(IPL(IL0*3-1))
        Y2=YD(IPL(IL0*3-1))
C
        XMN=XIMN
        XMX=XIMX
        YMN=YIMN
        YMX=YIMX
C
        IF (XAVG.LE.X1.AND.XAVG.LE.X2) XMN=MIN(X1,X2)
        IF (XAVG.GE.X1.AND.XAVG.GE.X2) XMX=MAX(X1,X2)
        IF (YAVG.LE.Y1.AND.YAVG.LE.Y2) YMN=MIN(Y1,Y2)
        IF (YAVG.GE.Y1.AND.YAVG.GE.Y2) YMX=MAX(Y1,Y2)
C
        INSD=0
C
        DO 302 IXI=1,NXI
          IF (XI(IXI).GE.XMN.AND.XI(IXI).LE.XMX) GO TO 301
          IF (INSD.EQ.0) GO TO 302
          IXIMX=IXI-1
          GO TO 303
  301     IF (INSD.EQ.1) GO TO 302
          INSD=1
          IXIMN=IXI
  302   CONTINUE
C
        IF (INSD.EQ.0) GO TO 313
C
        IXIMX=NXI
C
  303   DO 312 IYI=1,NYI
C
          YII=YI(IYI)
C
          IF (YII.LT.YMN.OR.YII.GT.YMX) GO TO 312
C
          DO 311 IXI=IXIMN,IXIMX
C
            XII=XI(IXI)
C
            L=0
            IF (VPDT(X1  ,Y1  ,X2,Y2,XII,YII)) 305,304,311
  304       L=1
  305       IF (VPDT(XAVG,YAVG,X1,Y1,XII,YII)) 311,306,307
  306       L=1
  307       IF (VPDT(XAVG,YAVG,X2,Y2,XII,YII)) 309,308,311
  308       L=1
  309       IZI=NXI*(IYI-1)+IXI
C
            IF (L.EQ.0) THEN
              NGP0=NGP0+1
              JIGP0=JIGP0+1
              IGP(JIGP0)=IZI
            ELSE
              IF (JIGP1.LE.NXINYI) THEN
                DO 310 JIGP1I=JIGP1,NXINYI
                  IF (IZI.EQ.IGP(JIGP1I)) GO TO 311
  310           CONTINUE
              END IF
              NGP1=NGP1+1
              JIGP1=JIGP1-1
              IGP(JIGP1)=IZI
            END IF
C
  311     CONTINUE
C
  312   CONTINUE
C
  313   JNGP0=JNGP0+1
        NGP(JNGP0)=NGP0
C
        JNGP1=JNGP1-1
        NGP(JNGP1)=NGP1
C
        JNGP0=JNGP0+1
        NGP(JNGP0)=0
C
        JNGP1=JNGP1-1
        NGP(JNGP1)=0
C
  314 CONTINUE
C
      RETURN
C
      END
C
C $Id: idpdrv.f,v 1.4 2000/08/22 15:06:54 haley Exp $
C                                                                      
C                Copyright (C)  2000
C        University Corporation for Atmospheric Research
C                All Rights Reserved
C
C This file is free software; you can redistribute it and/or modify
C it under the terms of the GNU General Public License as published
C by the Free Software Foundation; either version 2 of the License, or
C (at your option) any later version.
C
C This software is distributed in the hope that it will be useful, but
C WITHOUT ANY WARRANTY; without even the implied warranty of
C MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
C General Public License for more details.
C
C You should have received a copy of the GNU General Public License
C along with this software; if not, write to the Free Software
C Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
C USA.
C
      SUBROUTINE IDPDRV (NDP,XD,YD,ZD,NT,IPT,PD,WK)
C
C This subroutine estimates partial derivatives of the first and
C second order at the data points.
C
C The input arguments are as follows:
C
C   NDP       number of data points.
C
C   XD,YD,ZD  arrays of dimension NDP containing the X,
C             Y, and Z coordinates of the data points.
C
C   NT        number of triangles.
C
C   IPT       integer array of dimension 3*NT containing the
C             point numbers of the vertices of the triangles.
C
C The output argument is as follows:
C
C   PD        array of dimension 5*NDP, where the estimated
C             ZX, ZY, ZXX, ZXY, and ZYY values at the Ith
C             data point are to be stored as the (5*I-4)th,
C             (5*I-3)rd, (5*I-2)nd, (5*I-1)st, and (5*I)th
C             elements, respectively, where I = 1, 2, ...,
C             NDP.
C
C The other argument is as follows:
C
C   WK        array of dimension NDP used internally as a
C             work area.
C
C Declaration statements.
C
      DIMENSION XD(NDP),YD(NDP),ZD(NDP),IPT(3*NT),PD(5*NDP),WK(NDP)
      DIMENSION IPTI(3),XV(3),YV(3),ZV(3),ZXV(3),ZYV(3),W1(3),W2(3)
C
      DATA EPSLN / 1.E-6 /
C
C Clear the PD array.
C
      JPDMX=5*NDP
C
      DO 101 JPD=1,JPDMX
        PD(JPD)=0.0
  101 CONTINUE
C
      DO 102 IDP=1,NDP
        WK(IDP)=0.0
  102 CONTINUE
C
C Estimate ZX and ZY.
C
      DO 105 IT=1,NT
        JPT0=3*(IT-1)
        DO 103 IV=1,3
          JPT=JPT0+IV
          IDP=IPT(JPT)
          IPTI(IV)=IDP
          XV(IV)=XD(IDP)
          YV(IV)=YD(IDP)
          ZV(IV)=ZD(IDP)
  103   CONTINUE
        DX1=XV(2)-XV(1)
        DY1=YV(2)-YV(1)
        DZ1=ZV(2)-ZV(1)
        DX2=XV(3)-XV(1)
        DY2=YV(3)-YV(1)
        DZ2=ZV(3)-ZV(1)
        VPX=DY1*DZ2-DZ1*DY2
        VPY=DZ1*DX2-DX1*DZ2
        VPZ=DX1*DY2-DY1*DX2
        VPZMN=ABS(DX1*DX2+DY1*DY2)*EPSLN
        IF (ABS(VPZ).LE.VPZMN) GO TO 105
        D12=SQRT((XV(2)-XV(1))**2+(YV(2)-YV(1))**2)
        D23=SQRT((XV(3)-XV(2))**2+(YV(3)-YV(2))**2)
        D31=SQRT((XV(1)-XV(3))**2+(YV(1)-YV(3))**2)
        W1(1)=1.0/(D31*D12)
        W1(2)=1.0/(D12*D23)
        W1(3)=1.0/(D23*D31)
        W2(1)=VPZ*W1(1)
        W2(2)=VPZ*W1(2)
        W2(3)=VPZ*W1(3)
        DO 104 IV=1,3
          IDP=IPTI(IV)
          JPD0=5*(IDP-1)
          WI=(W1(IV)**2)*W2(IV)
          PD(JPD0+1)=PD(JPD0+1)+VPX*WI
          PD(JPD0+2)=PD(JPD0+2)+VPY*WI
          WK(IDP)=WK(IDP)+VPZ*WI
  104   CONTINUE
  105 CONTINUE
      DO 106 IDP=1,NDP
        IF (WK(IDP).NE.0.) THEN
          JPD0=5*(IDP-1)
          PD(JPD0+1)=-PD(JPD0+1)/WK(IDP)
          PD(JPD0+2)=-PD(JPD0+2)/WK(IDP)
        END IF
  106 CONTINUE
C
C Estimate ZXX, ZXY, and ZYY.
C
      DO 109 IT=1,NT
        JPT0=3*(IT-1)
        DO 107 IV=1,3
          JPT=JPT0+IV
          IDP=IPT(JPT)
          IPTI(IV)=IDP
          XV(IV)=XD(IDP)
          YV(IV)=YD(IDP)
          JPD0=5*(IDP-1)
          ZXV(IV)=PD(JPD0+1)
          ZYV(IV)=PD(JPD0+2)
  107   CONTINUE
        DX1=XV(2)-XV(1)
        DY1=YV(2)-YV(1)
        DZX1=ZXV(2)-ZXV(1)
        DZY1=ZYV(2)-ZYV(1)
        DX2=XV(3)-XV(1)
        DY2=YV(3)-YV(1)
        DZX2=ZXV(3)-ZXV(1)
        DZY2=ZYV(3)-ZYV(1)
        VPXX=DY1*DZX2-DZX1*DY2
        VPXY=DZX1*DX2-DX1*DZX2
        VPYX=DY1*DZY2-DZY1*DY2
        VPYY=DZY1*DX2-DX1*DZY2
        VPZ=DX1*DY2-DY1*DX2
        VPZMN=ABS(DX1*DX2+DY1*DY2)*EPSLN
        IF (ABS(VPZ).LE.VPZMN) GO TO 109
        D12=SQRT((XV(2)-XV(1))**2+(YV(2)-YV(1))**2)
        D23=SQRT((XV(3)-XV(2))**2+(YV(3)-YV(2))**2)
        D31=SQRT((XV(1)-XV(3))**2+(YV(1)-YV(3))**2)
        W1(1)=1.0/(D31*D12)
        W1(2)=1.0/(D12*D23)
        W1(3)=1.0/(D23*D31)
        W2(1)=VPZ*W1(1)
        W2(2)=VPZ*W1(2)
        W2(3)=VPZ*W1(3)
        DO 108 IV=1,3
          IDP=IPTI(IV)
          JPD0=5*(IDP-1)
          WI=(W1(IV)**2)*W2(IV)
          PD(JPD0+3)=PD(JPD0+3)+VPXX*WI
          PD(JPD0+4)=PD(JPD0+4)+(VPXY+VPYX)*WI
          PD(JPD0+5)=PD(JPD0+5)+VPYY*WI
  108   CONTINUE
  109 CONTINUE
      DO 110 IDP=1,NDP
        JPD0=5*(IDP-1)
        IF (WK(IDP).NE.0.) THEN
          PD(JPD0+3)=-PD(JPD0+3)/WK(IDP)
          PD(JPD0+4)=-PD(JPD0+4)/(2.0*WK(IDP))
          PD(JPD0+5)=-PD(JPD0+5)/WK(IDP)
        END IF
  110 CONTINUE
C
      RETURN
C
      END
C
C $Id: idptip.f,v 1.4 2000/08/22 15:06:54 haley Exp $
C                                                                      
C                Copyright (C)  2000
C        University Corporation for Atmospheric Research
C                All Rights Reserved
C
C This file is free software; you can redistribute it and/or modify
C it under the terms of the GNU General Public License as published
C by the Free Software Foundation; either version 2 of the License, or
C (at your option) any later version.
C
C This software is distributed in the hope that it will be useful, but
C WITHOUT ANY WARRANTY; without even the implied warranty of
C MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
C General Public License for more details.
C
C You should have received a copy of the GNU General Public License
C along with this software; if not, write to the Free Software
C Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
C USA.
C
      SUBROUTINE IDPTIP (XD,YD,ZD,NT,IPT,NL,IPL,PDD,ITI,XII,YII,ZII)
C
C This subroutine performs punctual interpolation or extrapolation,
C i.e., determines the Z value at a point.
C
C The input arguments are as follows:
C
C   XD,YD,ZD  Arrays of dimension NDP containing the X, Y, and
C             Z coordinates of the data points, where NDP is
C             the number of the data points.
C
C   NT        number of triangles.
C
C   IPT       integer array of dimension 3*NT containing the
C             point numbers of the vertices of the triangles.
C
C   NL        number of border line segments.
C
C   IPL       integer array of dimension 3*NL containing the
C             point numbers of the end points of the border
C             line segments and their respective triangle
C             numbers.
C
C   PDD       array of dimension 5*NDP containing the partial
C             derivatives at the data points.
C
C   ITI       triangle number of the triangle in which lies
C             the point for which interpolation is to be
C             performed.
C
C   XII,YII   X and Y coordinates of the point for which
C             interpolation is to be performed.
C
C The output argument is as follows:
C
C   ZII       interpolated Z value.
C
C Declaration statements.
C
      DIMENSION XD(*),YD(*),ZD(*),IPT(*),IPL(*),PDD(*)
C
      COMMON /IDCOMN/ INTY,ITTY,ALSP,BLSP,CLSP,XAVG,YAVG
      SAVE   /IDCOMN/
C
      COMMON /IDPT/ ITPV,X0,Y0,AP,BP,CP,DP,
     1              P00,P10,P20,P30,P40,P50,P01,P11,P21,P31,P41,
     2              P02,P12,P22,P32,P03,P13,P23,P04,P14,P05
      SAVE   /IDPT/
C
      DIMENSION   X(3),Y(3),Z(3),PD(15),
     1            ZU(3),ZV(3),ZUU(3),ZUV(3),ZVV(3)
C
      EQUIVALENCE (P5,P50)
C
C Statement function.
C
      VPDT(U1,V1,U2,V2,U3,V3)=(U1-U3)*(V2-V3)-(V1-V3)*(U2-U3)
C
C Preliminary processing.
C
      NTL=NT+NL
      IF (ITI.LE.NTL) GO TO 101
      IL1=ITI/NTL
      IF (INTY.NE.0) GO TO 115
      IL2=ITI-IL1*NTL
      IF (IL1.EQ.IL2) GO TO 107
      GO TO 112
C
C Calculation of ZII by interpolation.  Jump if the necessary
C coefficients have already been calculated.
C
  101 IF (ITI.EQ.ITPV) GO TO 105
C
C Load coordinate and partial-derivative values at the vertices.
C
      JIPT=3*(ITI-1)
      JPD=0
      DO 103 I=1,3
        JIPT=JIPT+1
        IDP=IPT(JIPT)
        X(I)=XD(IDP)
        Y(I)=YD(IDP)
        Z(I)=ZD(IDP)
        IF (INTY.NE.0) GO TO 103
        JPDD=5*(IDP-1)
        DO 102 KPD=1,5
          JPD=JPD+1
          JPDD=JPDD+1
          PD(JPD)=PDD(JPDD)
  102   CONTINUE
  103 CONTINUE
C
C Determine the coefficients for the coordinate system
C transformation from the X-Y system to the U-V system
C and vice-versa.
C
      X0=X(1)
      Y0=Y(1)
      A=X(2)-X0
      B=X(3)-X0
      C=Y(2)-Y0
      D=Y(3)-Y0
      AD=A*D
      BC=B*C
      DLT=AD-BC
      IF (INTY.NE.0) THEN
        AP=(D*(Z(2)-Z(1))-C*(Z(3)-Z(1)))/DLT
        BP=(B*(Z(1)-Z(2))-A*(Z(1)-Z(3)))/DLT
        CP=Z(1)-AP*X(1)-BP*Y(1)
        GO TO 106
      END IF
      AP= D/DLT
      BP=-B/DLT
      CP=-C/DLT
      DP= A/DLT
C
C Convert the partial derivatives at the vertices of the
C triangle for the U-V coordinate system.
C
      AA=A*A
      ACT2=2.0*A*C
      CC=C*C
      AB=A*B
      ADBC=AD+BC
      CD=C*D
      BB=B*B
      BDT2=2.0*B*D
      DD=D*D
      DO 104 I=1,3
        JPD=5*I
        ZU(I)=A*PD(JPD-4)+C*PD(JPD-3)
        ZV(I)=B*PD(JPD-4)+D*PD(JPD-3)
        ZUU(I)=AA*PD(JPD-2)+ACT2*PD(JPD-1)+CC*PD(JPD)
        ZUV(I)=AB*PD(JPD-2)+ADBC*PD(JPD-1)+CD*PD(JPD)
        ZVV(I)=BB*PD(JPD-2)+BDT2*PD(JPD-1)+DD*PD(JPD)
  104 CONTINUE
C
C Calculate the coefficients of the polynomial.
C
      P00=Z(1)
      P10=ZU(1)
      P01=ZV(1)
      P20=0.5*ZUU(1)
      P11=ZUV(1)
      P02=0.5*ZVV(1)
      H1=Z(2)-P00-P10-P20
      H2=ZU(2)-P10-ZUU(1)
      H3=ZUU(2)-ZUU(1)
      P30= 10.0*H1-4.0*H2+0.5*H3
      P40=-15.0*H1+7.0*H2    -H3
      P50=  6.0*H1-3.0*H2+0.5*H3
      H1=Z(3)-P00-P01-P02
      H2=ZV(3)-P01-ZVV(1)
      H3=ZVV(3)-ZVV(1)
      P03= 10.0*H1-4.0*H2+0.5*H3
      P04=-15.0*H1+7.0*H2    -H3
      P05=  6.0*H1-3.0*H2+0.5*H3
      RLU=SQRT(AA+CC)
      RLV=SQRT(BB+DD)
      THXU=ATAN2(C,A)
      THUV=ATAN2(D,B)-THXU
      CSUV=COS(THUV)
      P41=5.0*RLV*CSUV/RLU*P50
      P14=5.0*RLU*CSUV/RLV*P05
      H1=ZV(2)-P01-P11-P41
      H2=ZUV(2)-P11-4.0*P41
      P21= 3.0*H1-H2
      P31=-2.0*H1+H2
      H1=ZU(3)-P10-P11-P14
      H2=ZUV(3)-P11-4.0*P14
      P12= 3.0*H1-H2
      P13=-2.0*H1+H2
      THUS=ATAN2(D-C,B-A)-THXU
      THSV=THUV-THUS
      AA= SIN(THSV)/RLU
      BB=-COS(THSV)/RLU
      CC= SIN(THUS)/RLV
      DD= COS(THUS)/RLV
      AC=AA*CC
      AD=AA*DD
      BC=BB*CC
      G1=AA*AC*(3.0*BC+2.0*AD)
      G2=CC*AC*(3.0*AD+2.0*BC)
      H1=-AA*AA*AA*(5.0*AA*BB*P50+(4.0*BC+AD)*P41)
     1   -CC*CC*CC*(5.0*CC*DD*P05+(4.0*AD+BC)*P14)
      H2=0.5*ZVV(2)-P02-P12
      H3=0.5*ZUU(3)-P20-P21
      P22=(G1*H2+G2*H3-H1)/(G1+G2)
      P32=H2-P22
      P23=H3-P22
      ITPV=ITI
C
C Convert XII and YII to the U-V system.
C
  105 IF (INTY.NE.0) GO TO 106
      DX=XII-X0
      DY=YII-Y0
      U=AP*DX+BP*DY
      V=CP*DX+DP*DY
C
C Evaluate the polynomial.
C
      P0=P00+V*(P01+V*(P02+V*(P03+V*(P04+V*P05))))
      P1=P10+V*(P11+V*(P12+V*(P13+V*P14)))
      P2=P20+V*(P21+V*(P22+V*P23))
      P3=P30+V*(P31+V*P32)
      P4=P40+V*P41
      ZII=P0+U*(P1+U*(P2+U*(P3+U*(P4+U*P5))))
C
      RETURN
C
C Linear interpolation.
C
  106 ZII=AP*XII+BP*YII+CP
C
      RETURN
C
C Calculation of ZII by extrapolation in the rectangle.  First,
C check if the necessary coefficients have been calculated.
C
  107 IF (ITI.EQ.ITPV) GO TO 111
C
C Load coordinate and partial-derivative values at the end
C points of the border line segment.
C
      JIPL=3*(IL1-1)
      JPD=0
      DO 109 I=1,2
        JIPL=JIPL+1
        IDP=IPL(JIPL)
        X(I)=XD(IDP)
        Y(I)=YD(IDP)
        Z(I)=ZD(IDP)
        JPDD=5*(IDP-1)
        DO 108 KPD=1,5
          JPD=JPD+1
          JPDD=JPDD+1
          PD(JPD)=PDD(JPDD)
  108   CONTINUE
  109 CONTINUE
C
C Determine the coefficients for the coordinate system
C transformation from the X-Y system to the U-V system
C and vice-versa.
C
      X0=X(1)
      Y0=Y(1)
      A=Y(2)-Y(1)
      B=X(2)-X(1)
      C=-B
      D=A
      AD=A*D
      BC=B*C
      DLT=AD-BC
      AP= D/DLT
      BP=-B/DLT
      CP=-BP
      DP= AP
C
C Convert the partial derivatives at the end points of the
C border line segment for the U-V coordinate system.
C
      AA=A*A
      ACT2=2.0*A*C
      CC=C*C
      AB=A*B
      ADBC=AD+BC
      CD=C*D
      BB=B*B
      BDT2=2.0*B*D
      DD=D*D
      DO 110 I=1,2
        JPD=5*I
        ZU(I)=A*PD(JPD-4)+C*PD(JPD-3)
        ZV(I)=B*PD(JPD-4)+D*PD(JPD-3)
        ZUU(I)=AA*PD(JPD-2)+ACT2*PD(JPD-1)+CC*PD(JPD)
        ZUV(I)=AB*PD(JPD-2)+ADBC*PD(JPD-1)+CD*PD(JPD)
        ZVV(I)=BB*PD(JPD-2)+BDT2*PD(JPD-1)+DD*PD(JPD)
  110 CONTINUE
C
C Calculate the coefficients of the polynomial.
C
      P00=Z(1)
      P10=ZU(1)
      P01=ZV(1)
      P20=0.5*ZUU(1)
      P11=ZUV(1)
      P02=0.5*ZVV(1)
      H1=Z(2)-P00-P01-P02
      H2=ZV(2)-P01-ZVV(1)
      H3=ZVV(2)-ZVV(1)
      P03= 10.0*H1-4.0*H2+0.5*H3
      P04=-15.0*H1+7.0*H2    -H3
      P05=  6.0*H1-3.0*H2+0.5*H3
      H1=ZU(2)-P10-P11
      H2=ZUV(2)-P11
      P12= 3.0*H1-H2
      P13=-2.0*H1+H2
      P21=0.0
      P23=-ZUU(2)+ZUU(1)
      P22=-1.5*P23
      ITPV=ITI
C
C Convert XII and YII to the U-V system.
C
  111 DX=XII-X0
      DY=YII-Y0
      U=AP*DX+BP*DY
      V=CP*DX+DP*DY
C
C Evaluate the polynomial.
C
      P0=P00+V*(P01+V*(P02+V*(P03+V*(P04+V*P05))))
      P1=P10+V*(P11+V*(P12+V*P13))
      P2=P20+V*(P21+V*(P22+V*P23))
      ZII=P0+U*(P1+U*P2)
C
      RETURN
C
C Calculate ZII by extrapolation in the triangle.  First,
C check if the necessary coefficients have been calculated.
C
  112 IF (ITI.EQ.ITPV) GO TO 114
C
C Load coordinate and partial-derivative values at the vertex
C of the triangle.
C
      JIPL=3*IL2-2
      IDP=IPL(JIPL)
      X0=XD(IDP)
      Y0=YD(IDP)
      Z0=ZD(IDP)
      JPDD=5*(IDP-1)
      DO 113 KPD=1,5
        JPDD=JPDD+1
        PD(KPD)=PDD(JPDD)
  113 CONTINUE
C
C Calculate the coefficients of the polynomial.
C
      P00=Z0
      P10=PD(1)
      P01=PD(2)
      P20=0.5*PD(3)
      P11=PD(4)
      P02=0.5*PD(5)
      ITPV=ITI
C
C Convert XII and YII to the U-V system.
C
  114 U=XII-X0
      V=YII-Y0
C
C Evaluate the polynomial.
C
      P0=P00+V*(P01+V*P02)
      P1=P10+V*P11
      ZII=P0+U*(P1+U*P20)
C
      RETURN
C
C Linear extrapolation.
C
  115 DO 116 I=1,2
        X(I)=XD(IPL(3*(IL1-1)+I))
        Y(I)=YD(IPL(3*(IL1-1)+I))
        Z(I)=ZD(IPL(3*(IL1-1)+I))
  116 CONTINUE
C
C Use linear interpolation in one of the two triangles formed from the
C the quadrilateral defined by the ends of the edge segment and two
C points that are twice as far from (XAVG,YAVG) as those endpoints and
C have Z values on the least-squares-fit plane.
C
      IF ((X(1)-XAVG)*(X(1)-XAVG)+(Y(1)-YAVG)*(Y(1)-YAVG).LT.
     +    (X(2)-XAVG)*(X(2)-XAVG)+(Y(2)-YAVG)*(Y(2)-YAVG)) THEN
        X(3)=XAVG+2.*(X(1)-XAVG)
        Y(3)=YAVG+2.*(Y(1)-YAVG)
        Z(3)=ALSP*X(3)+BLSP*Y(3)+CLSP
        IF (VPDT(X(2),Y(2),X(3),Y(3),XII,YII).GT.0.) THEN
          X(1)=XAVG+2.*(X(2)-XAVG)
          Y(1)=YAVG+2.*(Y(2)-YAVG)
          Z(1)=ALSP*X(1)+BLSP*Y(1)+CLSP
          IF (VPDT(X(1),Y(1),X(3),Y(3),XII,YII).GT.0.) THEN
            ZII=ALSP*XII+BLSP*YII+CLSP
            RETURN
          END IF
        END IF
      ELSE
        X(3)=XAVG+2.*(X(2)-XAVG)
        Y(3)=YAVG+2.*(Y(2)-YAVG)
        Z(3)=ALSP*X(3)+BLSP*Y(3)+CLSP
        IF (VPDT(X(1),Y(1),X(3),Y(3),XII,YII).LT.0.) THEN
          X(2)=XAVG+2.*(X(1)-XAVG)
          Y(2)=YAVG+2.*(Y(1)-YAVG)
          Z(2)=ALSP*X(2)+BLSP*Y(2)+CLSP
          IF (VPDT(X(2),Y(2),X(3),Y(3),XII,YII).LT.0.) THEN
            ZII=ALSP*XII+BLSP*YII+CLSP
            RETURN
          END IF
        END IF
      END IF
C
      A=X(2)-X(1)
      B=X(3)-X(1)
      C=Y(2)-Y(1)
      D=Y(3)-Y(1)
      AD=A*D
      BC=B*C
      DLT=AD-BC
      AP=(D*(Z(2)-Z(1))-C*(Z(3)-Z(1)))/DLT
      BP=(B*(Z(1)-Z(2))-A*(Z(1)-Z(3)))/DLT
      CP=Z(1)-AP*X(1)-BP*Y(1)
      ZII=AP*XII+BP*YII+CP
C
      RETURN
C
      END
C
C $Id: idtang.f,v 1.4 2000/08/22 15:06:55 haley Exp $
C                                                                      
C                Copyright (C)  2000
C        University Corporation for Atmospheric Research
C                All Rights Reserved
C
C This file is free software; you can redistribute it and/or modify
C it under the terms of the GNU General Public License as published
C by the Free Software Foundation; either version 2 of the License, or
C (at your option) any later version.
C
C This software is distributed in the hope that it will be useful, but
C WITHOUT ANY WARRANTY; without even the implied warranty of
C MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
C General Public License for more details.
C
C You should have received a copy of the GNU General Public License
C along with this software; if not, write to the Free Software
C Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
C USA.
C
      SUBROUTINE IDTANG (NDP,XD,YD,NT,IPT,NL,IPL,IWL,IWP,WK)
C
C This subroutine performs triangulation.  It divides the X-Y
C plane into a number of triangles according to given data
C points in the plane, determines line segments that form the
C border of the data area, and determines the triangle numbers
C corresponding to the border line segments.
C
C At completion, point numbers of the vertices of each triangle
C are listed counter-clockwise.  Point numbers of the end points
C of each border line segment are listed counter-clockwise,
C listing order of the line segments being counter-clockwise.
C
C This subroutine calls the function IDXCHG.
C
C The input arguments are as follows:
C
C   NDP  number of data points.
C
C   XD   array of dimension NDP containing the
C        X coordinates of the data points.
C
C   YD   array of dimension NDP containing the
C        Y coordinates of the data points.
C
C The output arguments are as follows:
C
C   NT    number of triangles.
C
C   IPT  integer array of dimension 6*NDP-15, where the
C        point numbers of the vertices of the (IT)th
C        triangle are to be stored as the (3*IT-2)nd,
C        (3*IT-1)st, and (3*IT)th elements,
C        IT=1,2,...,NT.
C
C   NL   number of border line segments,
C
C   IPL  integer array of dimension 6*NDP, where the
C        point numbers of the end points of the (IL)th
C        border line segment and its respective triangle
C        number are to be stored as the (3*IL-2)nd,
C        (3*IL-1)st, and (3*IL)th elements,
C        IL=1,2,..., NL.
C
C The other arguments are as follows:
C
C   IWL  integer array of dimension 18*NDP used internally as a
C        work area.
C
C   IWP  integer array of dimension NDP used internally as a
C         work area.
C
C   WK   real array of dimension NDP used internally as a
C        work area.
C
C Declaration statements.
C
      DIMENSION XD(NDP),YD(NDP),IPT(6*NDP-15),IPL(6*NDP),IWL(18*NDP),
     +          IWP(NDP),WK(NDP),ITF(2)
C
      DATA EPSLN / 1.E-6 /, NREP / 200 /
C
C Statement functions.
C
      DSQF(U1,V1,U2,V2)=(U2-U1)**2+(V2-V1)**2
      SPDT(U1,V1,U2,V2,U3,V3)=(U2-U1)*(U3-U1)+(V2-V1)*(V3-V1)
      VPDT(U1,V1,U2,V2,U3,V3)=(V3-V1)*(U2-U1)-(U3-U1)*(V2-V1)
C
C Preliminary processing.
C
      IF (NDP.LT.4) THEN
        CALL SETER ('IDTANG (BIVAR) - INPUT PARAMETER NDP OUT OF RANGE',
     +1,1)
        RETURN
      END IF
C
      NDPM1=NDP-1
C
C Determine the closest pair of data points and the midpoint.
C
      DSQMN=DSQF(XD(1),YD(1),XD(2),YD(2))
      IPMN1=1
      IPMN2=2
      DO 102 IP1=1,NDPM1
        X1=XD(IP1)
        Y1=YD(IP1)
        IP1P1=IP1+1
        DO 101 IP2=IP1P1,NDP
          DSQI=DSQF(X1,Y1,XD(IP2),YD(IP2))
          IF (DSQI.EQ.0.) THEN
            CALL SETER ('IDTANG (BIVAR) - TWO OF THE INPUT DATA POINTS A
     +RE IDENTICAL',2,1)
            RETURN
          END IF
          IF (DSQI.GE.DSQMN) GO TO 101
          DSQMN=DSQI
          IPMN1=IP1
          IPMN2=IP2
  101   CONTINUE
  102 CONTINUE
      XDMP=(XD(IPMN1)+XD(IPMN2))/2.0
      YDMP=(YD(IPMN1)+YD(IPMN2))/2.0
C
C Sort the other (NDP-2) data points in ascending order of
C distance from the midpoint and store the sorted data point
C numbers in the IWP array.
C
      JP1=2
C
      DO 103 IP1=1,NDP
        IF (IP1.EQ.IPMN1.OR.IP1.EQ.IPMN2) GO TO 103
        JP1=JP1+1
        IWP(JP1)=IP1
        WK(JP1)=DSQF(XDMP,YDMP,XD(IP1),YD(IP1))
  103 CONTINUE
C
      DO 105 JP1=3,NDPM1
        DSQMN=WK(JP1)
        JPMN=JP1
        DO 104 JP2=JP1,NDP
          IF (WK(JP2).GE.DSQMN) GO TO 104
          DSQMN=WK(JP2)
          JPMN=JP2
  104   CONTINUE
        ITS=IWP(JP1)
        IWP(JP1)=IWP(JPMN)
        IWP(JPMN)=ITS
        WK(JPMN)=WK(JP1)
  105 CONTINUE
C
C If necessary, modify the ordering in such a way that the
C first three data points are not collinear.
C
      X1=XD(IPMN1)
      Y1=YD(IPMN1)
      X2=XD(IPMN2)
      Y2=YD(IPMN2)
C
      DO 106 JP=3,NDP
        IP=IWP(JP)
        SP=SPDT(XD(IP),YD(IP),X1,Y1,X2,Y2)
        VP=VPDT(XD(IP),YD(IP),X1,Y1,X2,Y2)
        IF (ABS(VP).GT.(ABS(SP)*EPSLN)) GO TO 107
  106 CONTINUE
C
      CALL SETER ('IDTANG (BIVAR) - ALL COLLINEAR DATA POINTS',3,1)
      RETURN
C
  107 IF (JP.EQ.3) GO TO 109
      JPMX=JP
C
      DO 108 JPC=4,JPMX
        JP=JPMX+4-JPC
        IWP(JP)=IWP(JP-1)
  108 CONTINUE
C
      IWP(3)=IP
C
C Form the first triangle.  Store point numbers of the vertices
C of the triangle in the IPT array, and store point numbers of
C the border line segments and the triangle number in the IPL
C array.
C
  109 IP1=IPMN1
      IP2=IPMN2
      IP3=IWP(3)
      IF (VPDT(XD(IP1),YD(IP1),XD(IP2),YD(IP2),XD(IP3),YD(IP3))
     +                                        .GE.0.) GO TO 110
      IP1=IPMN2
      IP2=IPMN1
  110 NT0=1
      NTT3=3
      IPT(1)=IP1
      IPT(2)=IP2
      IPT(3)=IP3
      NL0=3
      NLT3=9
      IPL(1)=IP1
      IPL(2)=IP2
      IPL(3)=1
      IPL(4)=IP2
      IPL(5)=IP3
      IPL(6)=1
      IPL(7)=IP3
      IPL(8)=IP1
      IPL(9)=1
C
C Add the remaining (NDP-3) data points, one by one.
C
      DO 130 JP1=4,NDP
C
        IP1=IWP(JP1)
        X1=XD(IP1)
        Y1=YD(IP1)
C
C Determine the first invisible and visible border line segments,
C ILIV and ILVS.
C
        DO 113 IL=1,NL0
          IP2=IPL(3*IL-2)
          IP3=IPL(3*IL-1)
          X2=XD(IP2)
          Y2=YD(IP2)
          X3=XD(IP3)
          Y3=YD(IP3)
          SP=SPDT(X1,Y1,X2,Y2,X3,Y3)
          VP=VPDT(X1,Y1,X2,Y2,X3,Y3)
          IF (IL.NE.1) GO TO 111
          IXVS=0
          IF (VP.LE.(ABS(SP)*(-EPSLN))) IXVS=1
          ILIV=1
          ILVS=1
          GO TO 113
  111     IXVSPV=IXVS
          IF (VP.GT.(ABS(SP)*(-EPSLN))) GO TO 112
          IXVS=1
          IF (IXVSPV.EQ.1) GO TO 113
          ILVS=IL
          IF (ILIV.NE.1) GO TO 114
          GO TO 113
  112     IXVS=0
          IF (IXVSPV.EQ.0) GO TO 113
          ILIV=IL
          IF (ILVS.NE.1) GO TO 114
  113   CONTINUE
C
        IF (ILIV.EQ.1.AND.ILVS.EQ.1) ILVS=NL0
  114   IF (ILVS.LT.ILIV) ILVS=ILVS+NL0
C
C Shifts (rotates) the IPL array to have the invisible border
C line segments contained in the first part of the IPL array.
C
        IF (ILIV.EQ.1) GO TO 117
        NLSH=ILIV-1
        NLSHT3=NLSH*3
C
        DO 115 JL1=1,NLSHT3
          JL2=JL1+NLT3
          IPL(JL2)=IPL(JL1)
  115   CONTINUE
C
        DO 116 JL1=1,NLT3
          JL2=JL1+NLSHT3
          IPL(JL1)=IPL(JL2)
  116   CONTINUE
C
        ILVS=ILVS-NLSH
C
C Add triangles to the IPT array, update border line segments
C in the IPL array, and set flags for the border lines segments
C to be re-examined in the IWL array.
C
  117   JWL=0
C
        DO 121 IL=ILVS,NL0
C
          ILT3=IL*3
          IPL1=IPL(ILT3-2)
          IPL2=IPL(ILT3-1)
          IT  =IPL(ILT3)
C
C Add a triangle to the IPT array.
C
          NT0=NT0+1
          NTT3=NTT3+3
          IPT(NTT3-2)=IPL2
          IPT(NTT3-1)=IPL1
          IPT(NTT3)  =IP1
C
C Update border line segments in the IPL array.
C
          IF (IL.NE.ILVS) GO TO 118
          IPL(ILT3-1)=IP1
          IPL(ILT3)  =NT0
  118     IF (IL.NE.NL0) GO TO 119
          NLN=ILVS+1
          NLNT3=NLN*3
          IPL(NLNT3-2)=IP1
          IPL(NLNT3-1)=IPL(1)
          IPL(NLNT3)  =NT0
C
C Determine the vertex that does not lie on the border
C line segments.
C
  119     ITT3=IT*3
          IPTI=IPT(ITT3-2)
          IF (IPTI.NE.IPL1.AND.IPTI.NE.IPL2) GO TO 120
          IPTI=IPT(ITT3-1)
          IF (IPTI.NE.IPL1.AND.IPTI.NE.IPL2) GO TO 120
          IPTI=IPT(ITT3)
C
C Check if an exchange is necessary.
C
  120     IF (IDXCHG(XD,YD,IP1,IPTI,IPL1,IPL2).EQ.0) GO TO 121
C
C Modify the IPT array when necessary.
C
          IPT(ITT3-2)=IPTI
          IPT(ITT3-1)=IPL1
          IPT(ITT3)  =IP1
          IPT(NTT3-1)=IPTI
          IF (IL.EQ.ILVS)  IPL(ILT3)=IT
          IF (IL.EQ.NL0.AND.IPL(3).EQ.IT) IPL(3)=NT0
C
C Set flags in the IWL array.
C
          JWL=JWL+4
          IWL(JWL-3)=IPL1
          IWL(JWL-2)=IPTI
          IWL(JWL-1)=IPTI
          IWL(JWL)  =IPL2
C
  121   CONTINUE
C
        NL0=NLN
        NLT3=NLNT3
        NLF=JWL/2
        IF (NLF.EQ.0) GO TO 130
C
C Improve triangulation.
C
        NTT3P3=NTT3+3
C
        DO 129 IREP=1,NREP
C
          DO 127 ILF=1,NLF
C
            IPL1=IWL(2*ILF-1)
            IPL2=IWL(2*ILF)
C
C Locate in the IPT array two triangles, one on each side of
C the flagged line segment.
C
            NTF=0
C
            DO 122 ITT3R=3,NTT3,3
              ITT3=NTT3P3-ITT3R
              IPT1=IPT(ITT3-2)
              IPT2=IPT(ITT3-1)
              IPT3=IPT(ITT3)
              IF (IPL1.NE.IPT1.AND.IPL1.NE.IPT2.AND.
     +            IPL1.NE.IPT3) GO TO 122
              IF (IPL2.NE.IPT1.AND.IPL2.NE.IPT2.AND.
     +            IPL2.NE.IPT3) GO TO 122
              NTF=NTF+1
              ITF(NTF)=ITT3/3
              IF (NTF.EQ.2) GO TO 123
  122       CONTINUE
C
            IF (NTF.LT.2) GO TO 127
C
C Determine the vertices of the triangles that do not lie
C on the line segment.
C
  123       IT1T3=ITF(1)*3
            IPTI1=IPT(IT1T3-2)
            IF (IPTI1.NE.IPL1.AND.IPTI1.NE.IPL2) GO TO 124
            IPTI1=IPT(IT1T3-1)
            IF (IPTI1.NE.IPL1.AND.IPTI1.NE.IPL2) GO TO 124
            IPTI1=IPT(IT1T3)
  124       IT2T3=ITF(2)*3
            IPTI2=IPT(IT2T3-2)
            IF (IPTI2.NE.IPL1.AND.IPTI2.NE.IPL2) GO TO 125
            IPTI2=IPT(IT2T3-1)
            IF (IPTI2.NE.IPL1.AND.IPTI2.NE.IPL2) GO TO 125
            IPTI2=IPT(IT2T3)
C
C Check if an exchange is necessary.
C
  125       IF (IDXCHG(XD,YD,IPTI1,IPTI2,IPL1,IPL2).EQ.0) GO TO 127
C
C Modify the IPT array when necessary.
C
            IPT(IT1T3-2)=IPTI1
            IPT(IT1T3-1)=IPTI2
            IPT(IT1T3)  =IPL1
            IPT(IT2T3-2)=IPTI2
            IPT(IT2T3-1)=IPTI1
            IPT(IT2T3)  =IPL2
C
C Set new flags.
C
            JWL=JWL+8
            IWL(JWL-7)=IPL1
            IWL(JWL-6)=IPTI1
            IWL(JWL-5)=IPTI1
            IWL(JWL-4)=IPL2
            IWL(JWL-3)=IPL2
            IWL(JWL-2)=IPTI2
            IWL(JWL-1)=IPTI2
            IWL(JWL)  =IPL1
C
            DO 126 JLT3=3,NLT3,3
              IPLJ1=IPL(JLT3-2)
              IPLJ2=IPL(JLT3-1)
              IF ((IPLJ1.EQ.IPL1.AND.IPLJ2.EQ.IPTI2).OR.
     +            (IPLJ2.EQ.IPL1.AND.IPLJ1.EQ.IPTI2))
     +                                  IPL(JLT3)=ITF(1)
              IF ((IPLJ1.EQ.IPL2.AND.IPLJ2.EQ.IPTI1).OR.
     +            (IPLJ2.EQ.IPL2.AND.IPLJ1.EQ.IPTI1))
     +                                  IPL(JLT3)=ITF(2)
  126       CONTINUE
C
  127     CONTINUE
C
          NLFC=NLF
          NLF=JWL/2
          IF (NLF.EQ.NLFC) GO TO 130
C
C Reset the IWL array for the next round.
C
          JWL1MN=2*NLFC+1
          NLFT2=NLF*2
C
          DO 128 JWL1=JWL1MN,NLFT2
            JWL=JWL1+1-JWL1MN
            IWL(JWL)=IWL(JWL1)
  128     CONTINUE
C
          NLF=JWL/2
C
  129   CONTINUE
C
  130 CONTINUE
C
C Re-arrange the IPT array so that the vertices of each triangle
C are listed counter-clockwise.
C
      DO 131 ITT3=3,NTT3,3
        IP1=IPT(ITT3-2)
        IP2=IPT(ITT3-1)
        IP3=IPT(ITT3)
        IF (VPDT(XD(IP1),YD(IP1),XD(IP2),YD(IP2),XD(IP3),YD(IP3))
     1                                         .GE.0.0) GO TO 131
        IPT(ITT3-2)=IP2
        IPT(ITT3-1)=IP1
  131 CONTINUE
C
      NT=NT0
      NL=NL0
C
      RETURN
C
      END
C
C $Id: icfell.f,v 1.5 2000/08/22 15:06:52 haley Exp $
C                                                                      
C                Copyright (C)  2000
C        University Corporation for Atmospheric Research
C                All Rights Reserved
C
C This file is free software; you can redistribute it and/or modify
C it under the terms of the GNU General Public License as published
C by the Free Software Foundation; either version 2 of the License, or
C (at your option) any later version.
C
C This software is distributed in the hope that it will be useful, but
C WITHOUT ANY WARRANTY; without even the implied warranty of
C MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
C General Public License for more details.
C
C You should have received a copy of the GNU General Public License
C along with this software; if not, write to the Free Software
C Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
C USA.
C
      FUNCTION ICFELL (MESSG,NERRF)
C
        CHARACTER*(*) MESSG
C
C ICFELL (which stands for "I Check For Errors on Lower Level") is used
C to check for the occurrence of a recoverable error in a lower-level
C routine and (perhaps) to update the current error flag and message.
C
C The old value of the current error flag is returned as the value of
C the function ICFELL.  If that value is zero, nothing else has been
C done.  If the value of ICFELL is non-zero, the following actions have
C been taken:
C
C   If MESSG is blank, the current error message has not been changed.
C
C   If MESSG is non-blank and its length is 6 or less, it should be the
C   name of the routine referencing ICFELL; the current error message
C   has been altered by prepending first a slash and then the value of
C   MESSG.  (Note that, if an error occurs several levels deep, the
C   effect of using ICFELL at each level is, effectively, to generate
C   traceback information in the error message.)
C
C   If MESSG is non-blank and its length is 7 or greater, its value
C   has become the new value of the current error message and the
C   previous error message has been printed.  This is used at the
C   beginning of an NCAR Graphics routine to check for an outstanding
C   error that the user has not recovered from and to ensure that the
C   message for the outstanding error gets printed.
C
C   If the expression NERRF has the value zero, the current error
C   flag has not been changed.
C
C   If the expression NERRF has a non-zero value, the value of the
C   current error flag has been made equal to that value.
C
C An example:  Assume that the routine "A" calls the routine "B" and
C that "B" detects an error and calls SETER with error number "32" and
C error message "B - ERROR HAS OCCURRED".  If recovery mode is not in
C effect, SETER prints the error message and STOPs, but, if recovery
C mode is in effect, control returns from "SETER" to "B" and thence to
C "A".  At that point, the statement
C
C   IF (ICFELL('A',13).NE.0) RETURN
C
C detects the fact that an error has occurred in "B" and results in a
C return from "A" to whatever routine called it.  It also changes the
C current error message to read "A/B - ERROR HAS OCCURRED" and changes
C the error number from "32" to "13".
C
C Another example:  Assume that the NCAR Graphics routine "A" is called
C when recovery mode is set and that it detects an error, calls SETER,
C and RETURNs to the user.  If the user neglects to check the error
C state and calls the routine "B" next, the statement
C
C   IF (ICFELL('B - UNCLEARED PRIOR ERROR',1).NE.0) RETURN
C
C ensures that the error message from routine "A" will be printed, that
C it will be replaced by "B - UNCLEARED PRIOR ERROR", and that "B" will
C not attempt to execute.
C
C
C The common blocks SECOMI and SECOMC are used to hold shared variables
C of types INTEGER and CHARACTER, respectively, for the routine SETER
C and associated routines.  For descriptions of these variables and for
C default values of them, see the block data routine SEBLDA.
C
        COMMON /SECOMI/ IERRU,IERRF,IRECF,LOMSG
        SAVE   /SECOMI/
C
        COMMON /SECOMC/ ERMSG
          CHARACTER*256 ERMSG
        SAVE   /SECOMC/
C
C Define a character temporary to use.
C
        CHARACTER*256 CTEMP
C
C Set the value of the function to the current value of the error flag.
C
        ICFELL=IERRF
C
C If the current error flag is non-zero, update its value and the value
C of the error message as directed by the values of the arguments of
C ICFELL.
C
        IF (IERRF.NE.0) THEN
C
          IF (MESSG.NE.' ') THEN
C
            LENMSG=ICLOEM(MESSG)
C
            IF (LENMSG.LE.6) THEN
              CTEMP=ERMSG
              ERMSG=MESSG(1:LENMSG)//'/'//CTEMP
            ELSE
              CALL EPRIN
              ERMSG=MESSG
            END IF
C
            LOMSG=ICLOEM(ERMSG)
C
          END IF
C
          IF (NERRF.NE.0) IERRF=NERRF
C
        END IF
C
C Done.
C
        RETURN
C
      END
C
C $Id: seter.f,v 1.6 2000/08/22 15:06:58 haley Exp $
C                                                                      
C                Copyright (C)  2000
C        University Corporation for Atmospheric Research
C                All Rights Reserved
C
C This file is free software; you can redistribute it and/or modify
C it under the terms of the GNU General Public License as published
C by the Free Software Foundation; either version 2 of the License, or
C (at your option) any later version.
C
C This software is distributed in the hope that it will be useful, but
C WITHOUT ANY WARRANTY; without even the implied warranty of
C MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
C General Public License for more details.
C
C You should have received a copy of the GNU General Public License
C along with this software; if not, write to the Free Software
C Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
C USA.
C
      SUBROUTINE SETER (MESSG,NERRF,IROPT)
C
        CHARACTER*(*) MESSG
C
C This routine is called when an error occurs.  The arguments specify
C an error message, an error number, and an option indicating whether
C the error is recoverable or fatal.  Exactly how the error is handled
C depends not only on the values of these arguments, but also on the
C values of internal variables specifying whether recovery mode is in
C effect and whether a previous error has occurred and not yet been
C recovered from and cleared.
C
C If no uncleared recoverable error has occurred and there are no
C errors in the arguments, then the following apply:
C
C   IROPT = 1, recovery mode active      -  just remember the error.
C   IROPT = 1, recovery mode not active  -  print and stop.
C   IROPT = 2                            -  print, dump, and stop.
C
C Input:
C
C   MESSG  - the error message (113 characters maximum).  The error
C            message should contain the name of the routine in which
C            the error occurred, followed by a blank, then a hyphen
C            (minus sign), then another blank, and then a short
C            description of the error.
C   NERRF  - the error number - must be non-zero.
C   IROPT  - the option - must have IROPT = 1 or 2.
C
C Error states:
C
C   1 - message length not positive.
C   2 - NERRF equal to 0.
C   3 - an unrecovered error followed by another error.
C   4 - bad value for IROPT (less than 1 or greater than 2).
C
C Force load of the BLOCK DATA subroutine.
C
        EXTERNAL SEBLDA
C
C The common blocks SECOMI and SECOMC are used to hold shared variables
C of types INTEGER and CHARACTER, respectively, for the routine SETER
C and associated routines.  For descriptions of these variables and for
C default values of them, see the block data routine SEBLDA.
C
        COMMON /SECOMI/ IERRU,IERRF,IRECF,LOMSG
        SAVE   /SECOMI/
C
        COMMON /SECOMC/ ERMSG
          CHARACTER*256 ERMSG
        SAVE   /SECOMC/
C
C The unit number for error messages is I1MACH(4).  Save that value,
C if it has not already been done.
C
        IF (IERRU.EQ.0) IERRU=I1MACH(4)
C
C Check for various error conditions.  The low-order bits of IERRC
C are used to keep track of which such errors have occurred.
C
        IERRC=0
C
C Check for a message of length zero or less.
C
        IF (LEN(MESSG).LE.0) THEN
          IERRC=IERRC+1
          WRITE (IERRU,1001)
        END IF
C
C Check for NERRF = 0.
C
        IF (NERRF.EQ.0) THEN
          IERRC=IERRC+2
          WRITE (IERRU,1002)
        END IF
C
C Check for a previous unrecovered error.
C
        IF (IERRF.NE.0) THEN
          IERRC=IERRC+4
          WRITE (IERRU,1003)
        END IF
C
C Check for an illegal value of the recovery flag.
C
        IF (IROPT.NE.1.AND.IROPT.NE.2) THEN
          IERRC=IERRC+8
          WRITE (IERRU,1004)
        END IF
C
C If one of the error conditions applies, print the appropriate
C information and quit.
C
        IF (IERRC.NE.0) THEN
          IF (MOD(IERRC/4,2).NE.0) THEN
            WRITE (IERRU,1005) IERRF,ERMSG(1:LOMSG)
          END IF
          IF (MOD(IERRC,2).EQ.0) THEN
            WRITE (IERRU,1006) NERRF,MESSG(1:ICLOEM(MESSG))
          END IF
          CALL FDUM
          STOP
        END IF
C
C Save the error message and error number.
C
        IERRF=NERRF
        ERMSG=MESSG
        LOMSG=ICLOEM(ERMSG)
C
C If recovery mode is activated and the error is recoverable, return
C to the caller for recovery action.
C
        IF (IRECF.EQ.1.AND.IROPT.EQ.1) RETURN
C
C Otherwise, print the error message.
C
        WRITE (IERRU,1007) IERRF,ERMSG(1:LOMSG)
C
C If the error is fatal, call the dump routine.
C
        IF (IROPT.EQ.2) CALL FDUM
C
C Quit.
C
        STOP
C
C Formats used above.
C
 1001 FORMAT (' ERROR    1 IN SETER - MESSAGE LENGTH IS NOT POSITIVE')
 1002 FORMAT (' ERROR    2 IN SETER - ILLEGAL VALUE FOR ERROR NUMBER')
 1003 FORMAT (' ERROR    3 IN SETER - AN UNCLEARED PRIOR ERROR EXISTS')
 1004 FORMAT (' ERROR    4 IN SETER - ILLEGAL VALUE FOR RECOVERY FLAG')
 1005 FORMAT (' ... MESSAGE FOR UNCLEARED PRIOR ERROR IS AS FOLLOWS:'/
     +        ' ... ERROR ',I4,' IN ',A)
 1006 FORMAT (' ... MESSAGE FOR CURRENT CALL TO SETER IS AS FOLLOWS:'/
     +        ' ... ERROR ',I4,' IN ',A)
 1007 FORMAT (' ERROR ',I4,' IN ',A)
C
      END
C
C $Id: idlsqf.f,v 1.3 2000/08/22 15:06:53 haley Exp $
C                                                                      
C                Copyright (C)  2000
C        University Corporation for Atmospheric Research
C                All Rights Reserved
C
C This file is free software; you can redistribute it and/or modify
C it under the terms of the GNU General Public License as published
C by the Free Software Foundation; either version 2 of the License, or
C (at your option) any later version.
C
C This software is distributed in the hope that it will be useful, but
C WITHOUT ANY WARRANTY; without even the implied warranty of
C MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
C General Public License for more details.
C
C You should have received a copy of the GNU General Public License
C along with this software; if not, write to the Free Software
C Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
C USA.
C
      SUBROUTINE IDLSQF (X,Y,Z,N,A,B,C,XAVG,YAVG)
C
        DIMENSION X(N),Y(N),Z(N)
C
C IDLSQF fits a plane to the data defined by points (X(I),Y(I),Z(I)),
C for I from 1 to N.  It returns the coefficients A, B, and C in the
C equation "Z=AX+BY+C".  It also returns the average X and average Y.
C
        SA=0.
        SB=0.
        SC=0.
        SD=0.
        SE=0.
        SF=0.
        SG=0.
        SH=0.
C
        DO 101 I=1,N
          SA=SA+X(I)
          SB=SB+Y(I)
          SC=SC+Z(I)
          SD=SD+X(I)*X(I)
          SE=SE+Y(I)*Y(I)
          SF=SF+X(I)*Y(I)
          SG=SG+Y(I)*Z(I)
          SH=SH+X(I)*Z(I)
  101   CONTINUE
C
        RN=REAL(N)
C
        A=((RN*SH-SA*SC)*(RN*SE-SB*SB)-(RN*SF-SA*SB)*(RN*SG-SB*SC))/
     +    ((RN*SD-SA*SA)*(RN*SE-SB*SB)-(RN*SF-SA*SB)*(RN*SF-SA*SB))
        B=((RN*SD-SA*SA)*(RN*SG-SB*SC)-(RN*SF-SA*SB)*(RN*SH-SA*SC))/
     +    ((RN*SD-SA*SA)*(RN*SE-SB*SB)-(RN*SF-SA*SB)*(RN*SF-SA*SB))
        C=(SC-SA*A-SB*B)/RN
C
        XAVG=SA/RN
        YAVG=SB/RN
C
C Done.
C
        RETURN
C
      END
C
C $Id: idxchg.f,v 1.4 2000/08/22 15:06:55 haley Exp $
C                                                                      
C                Copyright (C)  2000
C        University Corporation for Atmospheric Research
C                All Rights Reserved
C
C This file is free software; you can redistribute it and/or modify
C it under the terms of the GNU General Public License as published
C by the Free Software Foundation; either version 2 of the License, or
C (at your option) any later version.
C
C This software is distributed in the hope that it will be useful, but
C WITHOUT ANY WARRANTY; without even the implied warranty of
C MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
C General Public License for more details.
C
C You should have received a copy of the GNU General Public License
C along with this software; if not, write to the Free Software
C Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
C USA.
C
      FUNCTION IDXCHG (X,Y,I1,I2,I3,I4)
C
C This function determines whether or not two adjacent triangles in the
C triangulation should be "exchanged".  By default, it uses a criterion,
C due to C. L. Lawson, that maximizes the minimum angle occurring in the
C triangles, but if ITTY is set non-zero, the function IDXCHN is called
C instead.
C
C The input arguments are as follows:
C
C   X and Y are arrays containing the coordinates of the data points.
C
C   I1, I2, I3, and I4 are point numbers of four points P1, P2, P3,
C   and P4 forming a quadrilateral with P3 and P4 connected diagonally
C   to form two triangles.
C
C The function returns an integer 1 if P1 and P2 should be connected
C instead of P3 and P4 (an "exchange"); otherwise, it returns 0.
C
C Declaration statements.
C
      COMMON /IDCOMN/ INTY,ITTY,ALSP,BLSP,CLSP,XAVG,YAVG
      SAVE   /IDCOMN/
C
      DIMENSION X(*), Y(*)
C
      EQUIVALENCE (C2SQ,C1SQ),(A3SQ,B2SQ),(B3SQ,A1SQ),
     1            (A4SQ,B1SQ),(B4SQ,A2SQ),(C4SQ,C3SQ)
C
      DATA  EPSLN / 1.E-6 /
C
C See if the new triangulation is to be produced and, if so, do it.
C
      IF (ITTY.NE.0) THEN
        IDXCHG=IDXCHN(X,Y,I1,I2,I3,I4)
        RETURN
      END IF
C
C Preliminary processing.
C
      X1=X(I1)
      Y1=Y(I1)
      X2=X(I2)
      Y2=Y(I2)
      X3=X(I3)
      Y3=Y(I3)
      X4=X(I4)
      Y4=Y(I4)
C
C Calculation.
C
      IDX=0
C
      U3=(Y2-Y3)*(X1-X3)-(X2-X3)*(Y1-Y3)
      U4=(Y1-Y4)*(X2-X4)-(X1-X4)*(Y2-Y4)
C
      IF (U3*U4.GT.0.) THEN
        U1=(Y3-Y1)*(X4-X1)-(X3-X1)*(Y4-Y1)
        U2=(Y4-Y2)*(X3-X2)-(X4-X2)*(Y3-Y2)
        A1SQ=(X1-X3)**2+(Y1-Y3)**2
        B1SQ=(X4-X1)**2+(Y4-Y1)**2
        C1SQ=(X3-X4)**2+(Y3-Y4)**2
        A2SQ=(X2-X4)**2+(Y2-Y4)**2
        B2SQ=(X3-X2)**2+(Y3-Y2)**2
        C3SQ=(X2-X1)**2+(Y2-Y1)**2
        S1SQ=U1*U1/(C1SQ*AMAX1(A1SQ,B1SQ))
        S2SQ=U2*U2/(C2SQ*AMAX1(A2SQ,B2SQ))
        S3SQ=U3*U3/(C3SQ*AMAX1(A3SQ,B3SQ))
        S4SQ=U4*U4/(C4SQ*AMAX1(A4SQ,B4SQ))
        IF ((MIN(S3SQ,S4SQ)-MIN(S1SQ,S2SQ)).GT.EPSLN) IDX=1
      END IF
C
      IDXCHG=IDX
C
C Done.
C
      RETURN
C
      END
C
C $Id: icloem.f,v 1.4 2000/08/22 15:06:52 haley Exp $
C                                                                      
C                Copyright (C)  2000
C        University Corporation for Atmospheric Research
C                All Rights Reserved
C
C This file is free software; you can redistribute it and/or modify
C it under the terms of the GNU General Public License as published
C by the Free Software Foundation; either version 2 of the License, or
C (at your option) any later version.
C
C This software is distributed in the hope that it will be useful, but
C WITHOUT ANY WARRANTY; without even the implied warranty of
C MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
C General Public License for more details.
C
C You should have received a copy of the GNU General Public License
C along with this software; if not, write to the Free Software
C Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
C USA.
C
      FUNCTION ICLOEM (MESSG)
C
        CHARACTER*(*) MESSG
C
C ICLOEM(MESSG) is the index of the last non-blank character in MESSG.
C
        DO 101 I=LEN(MESSG),1,-1
         IF (MESSG(I:I).NE.' ') THEN
           ICLOEM=I
           RETURN
         END IF
  101   CONTINUE
C
        ICLOEM=1
C
        RETURN
C
      END
C
C $Id: eprin.f,v 1.5 2000/08/22 15:06:51 haley Exp $
C                                                                      
C                Copyright (C)  2000
C        University Corporation for Atmospheric Research
C                All Rights Reserved
C
C This file is free software; you can redistribute it and/or modify
C it under the terms of the GNU General Public License as published
C by the Free Software Foundation; either version 2 of the License, or
C (at your option) any later version.
C
C This software is distributed in the hope that it will be useful, but
C WITHOUT ANY WARRANTY; without even the implied warranty of
C MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
C General Public License for more details.
C
C You should have received a copy of the GNU General Public License
C along with this software; if not, write to the Free Software
C Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
C USA.
C
      SUBROUTINE EPRIN
C
C This routine just prints the current error message, if any.
C
C The common blocks SECOMI and SECOMC are used to hold shared variables
C of types INTEGER and CHARACTER, respectively, for the routine SETER
C and associated routines.  For descriptions of these variables and for
C default values of them, see the block data routine SEBLDA.
C
        COMMON /SECOMI/ IERRU,IERRF,IRECF,LOMSG
        SAVE   /SECOMI/
C
        COMMON /SECOMC/ ERMSG
          CHARACTER*256 ERMSG
        SAVE   /SECOMC/
C
C Do it.
C
        IF (IERRF.NE.0) WRITE (IERRU,'('' ERROR '',I4,'' IN '',A)')
     +                                              IERRF,ERMSG(1:LOMSG)
C
C Done.
C
        RETURN
C
      END
C
C $Id: seblda.f,v 1.4 2000/08/22 15:06:58 haley Exp $
C                                                                      
C                Copyright (C)  2000
C        University Corporation for Atmospheric Research
C                All Rights Reserved
C
C This file is free software; you can redistribute it and/or modify
C it under the terms of the GNU General Public License as published
C by the Free Software Foundation; either version 2 of the License, or
C (at your option) any later version.
C
C This software is distributed in the hope that it will be useful, but
C WITHOUT ANY WARRANTY; without even the implied warranty of
C MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
C General Public License for more details.
C
C You should have received a copy of the GNU General Public License
C along with this software; if not, write to the Free Software
C Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
C USA.
C
      BLOCK DATA SEBLDA
C
C The common blocks SECOMI and SECOMC are used to hold shared variables
C of types INTEGER and CHARACTER, respectively, for the routine SETER
C and associated routines.  For descriptions of these variables and for
C default values of them, see the block data routine SEBLDA.
C
        COMMON /SECOMI/ IERRU,IERRF,IRECF,LOMSG
        SAVE   /SECOMI/
C
        COMMON /SECOMC/ ERMSG
          CHARACTER*256 ERMSG
        SAVE   /SECOMC/
C
C IERRU is the logical unit for error messages.  Its default value is
C zero, which serves as a signal to SETER that the proper value should
C be gotten from I1MACH.
C
        DATA IERRU / 0 /
C
C IERRF is the error flag.  Initially, it is zero.  A non-zero value
C means that there is an uncleared prior error.
C
        DATA IERRF / 0 /
C
C IRECF is the recovery flag, which can have values of 1, which turns
C recovery mode on, or 2, which turns recovery mode off and causes all
C errors to be treated as fatal errors.  The latter is the default.
C
        DATA IRECF / 2 /
C
C LOMSG is the actual length of the error message in ERMSG.
C
        DATA LOMSG / 1 /
C
C ERMSG is the text of the current error message, only meaningful when
C IERRF is non-zero.
C
        DATA ERMSG / ' ' /
C
      END
C
C $Id: idxchn.f,v 1.3 2000/08/22 15:06:55 haley Exp $
C                                                                      
C                Copyright (C)  2000
C        University Corporation for Atmospheric Research
C                All Rights Reserved
C
C This file is free software; you can redistribute it and/or modify
C it under the terms of the GNU General Public License as published
C by the Free Software Foundation; either version 2 of the License, or
C (at your option) any later version.
C
C This software is distributed in the hope that it will be useful, but
C WITHOUT ANY WARRANTY; without even the implied warranty of
C MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
C General Public License for more details.
C
C You should have received a copy of the GNU General Public License
C along with this software; if not, write to the Free Software
C Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
C USA.
C
      FUNCTION IDXCHN (X,Y,I1,I2,I3,I4)
C
        DIMENSION X(*),Y(*)
C
C This function determines whether or not two adjacent triangles in
C the triangulation should be "exchanged".  It generates the Delaunay
C triangulation, used by the package "nngridr", by preferring triangles
C whose circumscribed circles have no data points in their interiors.
C
C The input arguments are the X and Y coordinate arrays of the data and
C the indices of four data points - P1, P2, P3, and P4 - forming a
C quadrilateral, with P3 and P4 connected diagonally.  The function
C returns a 1 when it decides that P1 and P2 should be connected
C instead (an "exchange") and a 0 otherwise.
C
        VPDT(U1,V1,U2,V2,U3,V3)=(U1-U3)*(V2-V3)-(V1-V3)*(U2-U3)
C
C Extract the required coordinates.
C
        X1=X(I1)
        Y1=Y(I1)
        X2=X(I2)
        Y2=Y(I2)
        X3=X(I3)
        Y3=Y(I3)
        X4=X(I4)
        Y4=Y(I4)
C
C Assume no exchange (value = 0) until we find out otherwise.
C
        IVAL=0
C
C If interior angles are less than 180 degrees, and either the point
C (X2,Y2) is inside the circle passing through the points (X1,Y1),
C (X3,Y3), and (X4,Y4) or the point (X1,Y1) is inside the circle
C passing through the points (X2,Y2), (X3,Y3), and (X4,Y4), return
C the value 1; this will cause the triangulation to be changed so
C as to use the diagonal from point 1 to point 2, rather than the
C diagonal from point 3 to point 4.  (Note that, mathematically, only
C one of the tests inside the IF-block should need to be made, because
C passing one of the tests implies passing the other, as well.)
C
        IF (VPDT(X1,Y1,X2,Y2,X3,Y3)*VPDT(X1,Y1,X2,Y2,X4,Y4).LT.0.) THEN
          CALL IDGTCP (X1,Y1,X3,Y3,X4,Y4,XC,YC)
          IF ((X2-XC)*(X2-XC)+(Y2-YC)*(Y2-YC).LT.
     +        (X1-XC)*(X1-XC)+(Y1-YC)*(Y1-YC)) IVAL=1
          CALL IDGTCP (X2,Y2,X3,Y3,X4,Y4,XC,YC)
          IF ((X1-XC)*(X1-XC)+(Y1-YC)*(Y1-YC).LT.
     +        (X2-XC)*(X2-XC)+(Y2-YC)*(Y2-YC)) IVAL=1
        END IF
C
C Set the value of the function for return to the caller.
C
        IDXCHN=IVAL
C
C Done.
C
        RETURN
C
      END










C
C $Id: idgtcp.f,v 1.3 2000/08/22 15:06:53 haley Exp $
C                                                                      
C                Copyright (C)  2000
C        University Corporation for Atmospheric Research
C                All Rights Reserved
C
C This file is free software; you can redistribute it and/or modify
C it under the terms of the GNU General Public License as published
C by the Free Software Foundation; either version 2 of the License, or
C (at your option) any later version.
C
C This software is distributed in the hope that it will be useful, but
C WITHOUT ANY WARRANTY; without even the implied warranty of
C MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
C General Public License for more details.
C
C You should have received a copy of the GNU General Public License
C along with this software; if not, write to the Free Software
C Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
C USA.
C
      SUBROUTINE IDGTCP (X1,Y1,X2,Y2,X3,Y3,XC,YC)
C
C This subroutine, given the coordinates of three points in the plane -
C (X1,Y1), (X2,Y2), and (X3,Y3) - returns the coordinates of the center
C of the circle passing through those three points: (XC,YC).
C
        A=2.*(X1-X2)
        B=2.*(Y1-Y2)
        C=X1*X1-X2*X2+Y1*Y1-Y2*Y2
        D=2.*(X1-X3)
        E=2.*(Y1-Y3)
        F=X1*X1-X3*X3+Y1*Y1-Y3*Y3
        XC=(C*E-B*F)/(A*E-B*D)
        YC=(C*D-A*F)/(B*D-A*E)
C
        RETURN
C
      END
C
C $Id: fdum.f,v 1.5 2000/08/22 15:06:52 haley Exp $
C                                                                      
C                Copyright (C)  2000
C        University Corporation for Atmospheric Research
C                All Rights Reserved
C
C This file is free software; you can redistribute it and/or modify
C it under the terms of the GNU General Public License as published
C by the Free Software Foundation; either version 2 of the License, or
C (at your option) any later version.
C
C This software is distributed in the hope that it will be useful, but
C WITHOUT ANY WARRANTY; without even the implied warranty of
C MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
C General Public License for more details.
C
C You should have received a copy of the GNU General Public License
C along with this software; if not, write to the Free Software
C Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
C USA.
C
      SUBROUTINE FDUM
C
C This routine is called when a fatal error occurs.  The default
C version does nothing.  A version may be supplied by the user to
C do something more useful.
C
        RETURN
C
      END
      INTEGER FUNCTION I1MACH (I)       
C***BEGIN PROLOGUE  I1MACH      
C***PURPOSE  Return integer machine dependent constants.
C***LIBRARY   SLATEC    
C***CATEGORY  R1
C***TYPE      INTEGER (I1MACH-I)
C***KEYWORDS  MACHINE CONSTANTS 
C***AUTHOR  Fox, P. A., (Bell Labs)     
C           Hall, A. D., (Bell Labs)    
C           Schryer, N. L., (Bell Labs) 
C***DESCRIPTION 
C       
C   I1MACH can be used to obtain machine-dependent parameters for the   
C   local machine environment.  It is a function subprogram with one    
C   (input) argument and can be referenced as follows:  
C       
C        K = I1MACH(I)  
C       
C   where I=1,...,16.  The (output) value of K above is determined by   
C   the (input) value of I.  The results for various values of I are    
C   discussed below.    
C       
C   I/O unit numbers:   
C     I1MACH( 1) = the standard input unit.     
C     I1MACH( 2) = the standard output unit.    
C     I1MACH( 3) = the standard punch unit.     
C     I1MACH( 4) = the standard error message unit.     
C       
C   Words:      
C     I1MACH( 5) = the number of bits per integer storage unit. 
C     I1MACH( 6) = the number of characters per integer storage unit.   
C       
C   Integers:   
C     assume integers are represented in the S-digit, base-A form       
C       
C                sign ( X(S-1)*A**(S-1) + ... + X(1)*A + X(0) ) 
C       
C                where 0 .LE. X(I) .LT. A for I=0,...,S-1.      
C     I1MACH( 7) = A, the base. 
C     I1MACH( 8) = S, the number of base-A digits.      
C     I1MACH( 9) = A**S - 1, the largest magnitude.     
C       
C   Floating-Point Numbers:     
C     Assume floating-point numbers are represented in the T-digit,     
C     base-B form       
C                sign (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )   
C       
C                where 0 .LE. X(I) .LT. B for I=1,...,T,
C                0 .LT. X(1), and EMIN .LE. E .LE. EMAX.
C     I1MACH(10) = B, the base. 
C       
C   Single-Precision:   
C     I1MACH(11) = T, the number of base-B digits.      
C     I1MACH(12) = EMIN, the smallest exponent E.       
C     I1MACH(13) = EMAX, the largest exponent E.
C       
C   Double-Precision:   
C     I1MACH(14) = T, the number of base-B digits.      
C     I1MACH(15) = EMIN, the smallest exponent E.       
C     I1MACH(16) = EMAX, the largest exponent E.
C       
C   To alter this function for a particular environment, the desired    
C   set of DATA statements should be activated by removing the C from   
C   column 1.  Also, the values of I1MACH(1) - I1MACH(4) should be      
C   checked for consistency with the local operating system.    
C       
C***REFERENCES  P. A. Fox, A. D. Hall and N. L. Schryer, Framework for  
C                 a portable library, ACM Transactions on Mathematical  
C                 Software 4, 2 (June 1978), pp. 177-188.       
C***ROUTINES CALLED  (NONE)     
C***REVISION HISTORY  (YYMMDD)  
C   750101  DATE WRITTEN
C   891012  Added VAX G-floating constants.  (WRB)      
C   891012  REVISION DATE from Version 3.2      
C   891214  Prologue converted to Version 4.0 format.  (BAB)    
C   900618  Added DEC RISC constants.  (WRB)    
C   900723  Added IBM RS 6000 constants.  (WRB) 
C   901009  Correct I1MACH(7) for IBM Mainframes. Should be 2 not 16.   
C           (RWC)       
C   910710  Added HP 730 constants.  (SMR)      
C   911114  Added Convex IEEE constants.  (WRB) 
C   920121  Added SUN -r8 compiler option constants.  (WRB)     
C   920229  Added Touchstone Delta i860 constants.  (WRB)       
C   920501  Reformatted the REFERENCES section.  (WRB)  
C   920625  Added CONVEX -p8 and -pd8 compiler option constants.
C           (BKS, WRB)  
C   930201  Added DEC Alpha and SGI-IRIS constants.  (RWC and WRB)      
c   010734  commented out some of the c preprocessing options at the end 
C***END PROLOGUE  I1MACH
C       
      INTEGER IMACH(16),OUTPUT  
      SAVE IMACH
      EQUIVALENCE (IMACH(4),OUTPUT)     
C       
C     MACHINE CONSTANTS FOR THE AMIGA   
C     ABSOFT COMPILER   
C       
C     DATA IMACH( 1) /          5 /     
C     DATA IMACH( 2) /          6 /     
C     DATA IMACH( 3) /          5 /     
C     DATA IMACH( 4) /          6 /     
C     DATA IMACH( 5) /         32 /     
C     DATA IMACH( 6) /          4 /     
C     DATA IMACH( 7) /          2 /     
C     DATA IMACH( 8) /         31 /     
C     DATA IMACH( 9) / 2147483647 /     
C     DATA IMACH(10) /          2 /     
C     DATA IMACH(11) /         24 /     
C     DATA IMACH(12) /       -126 /     
C     DATA IMACH(13) /        127 /     
C     DATA IMACH(14) /         53 /     
C     DATA IMACH(15) /      -1022 /     
C     DATA IMACH(16) /       1023 /     
C       
C     MACHINE CONSTANTS FOR THE APOLLO  
C       
C     DATA IMACH( 1) /          5 /     
C     DATA IMACH( 2) /          6 /     
C     DATA IMACH( 3) /          6 /     
C     DATA IMACH( 4) /          6 /     
C     DATA IMACH( 5) /         32 /     
C     DATA IMACH( 6) /          4 /     
C     DATA IMACH( 7) /          2 /     
C     DATA IMACH( 8) /         31 /     
C     DATA IMACH( 9) / 2147483647 /     
C     DATA IMACH(10) /          2 /     
C     DATA IMACH(11) /         24 /     
C     DATA IMACH(12) /       -125 /     
C     DATA IMACH(13) /        129 /     
C     DATA IMACH(14) /         53 /     
C     DATA IMACH(15) /      -1021 /     
C     DATA IMACH(16) /       1025 /     
C       
C     MACHINE CONSTANTS FOR THE BURROUGHS 1700 SYSTEM   
C       
C     DATA IMACH( 1) /          7 /     
C     DATA IMACH( 2) /          2 /     
C     DATA IMACH( 3) /          2 /     
C     DATA IMACH( 4) /          2 /     
C     DATA IMACH( 5) /         36 /     
C     DATA IMACH( 6) /          4 /     
C     DATA IMACH( 7) /          2 /     
C     DATA IMACH( 8) /         33 /     
C     DATA IMACH( 9) / Z1FFFFFFFF /     
C     DATA IMACH(10) /          2 /     
C     DATA IMACH(11) /         24 /     
C     DATA IMACH(12) /       -256 /     
C     DATA IMACH(13) /        255 /     
C     DATA IMACH(14) /         60 /     
C     DATA IMACH(15) /       -256 /     
C     DATA IMACH(16) /        255 /     
C       
C     MACHINE CONSTANTS FOR THE BURROUGHS 5700 SYSTEM   
C       
C     DATA IMACH( 1) /          5 /     
C     DATA IMACH( 2) /          6 /     
C     DATA IMACH( 3) /          7 /     
C     DATA IMACH( 4) /          6 /     
C     DATA IMACH( 5) /         48 /     
C     DATA IMACH( 6) /          6 /     
C     DATA IMACH( 7) /          2 /     
C     DATA IMACH( 8) /         39 /     
C     DATA IMACH( 9) / O0007777777777777 /      
C     DATA IMACH(10) /          8 /     
C     DATA IMACH(11) /         13 /     
C     DATA IMACH(12) /        -50 /     
C     DATA IMACH(13) /         76 /     
C     DATA IMACH(14) /         26 /     
C     DATA IMACH(15) /        -50 /     
C     DATA IMACH(16) /         76 /     
C       
C     MACHINE CONSTANTS FOR THE BURROUGHS 6700/7700 SYSTEMS     
C       
C     DATA IMACH( 1) /          5 /     
C     DATA IMACH( 2) /          6 /     
C     DATA IMACH( 3) /          7 /     
C     DATA IMACH( 4) /          6 /     
C     DATA IMACH( 5) /         48 /     
C     DATA IMACH( 6) /          6 /     
C     DATA IMACH( 7) /          2 /     
C     DATA IMACH( 8) /         39 /     
C     DATA IMACH( 9) / O0007777777777777 /      
C     DATA IMACH(10) /          8 /     
C     DATA IMACH(11) /         13 /     
C     DATA IMACH(12) /        -50 /     
C     DATA IMACH(13) /         76 /     
C     DATA IMACH(14) /         26 /     
C     DATA IMACH(15) /     -32754 /     
C     DATA IMACH(16) /      32780 /     
C       
C     MACHINE CONSTANTS FOR THE CDC 170/180 SERIES USING NOS/VE 
C       
C     DATA IMACH( 1) /          5 /     
C     DATA IMACH( 2) /          6 /     
C     DATA IMACH( 3) /          7 /     
C     DATA IMACH( 4) /          6 /     
C     DATA IMACH( 5) /         64 /     
C     DATA IMACH( 6) /          8 /     
C     DATA IMACH( 7) /          2 /     
C     DATA IMACH( 8) /         63 /     
C     DATA IMACH( 9) / 9223372036854775807 /    
C     DATA IMACH(10) /          2 /     
C     DATA IMACH(11) /         47 /     
C     DATA IMACH(12) /      -4095 /     
C     DATA IMACH(13) /       4094 /     
C     DATA IMACH(14) /         94 /     
C     DATA IMACH(15) /      -4095 /     
C     DATA IMACH(16) /       4094 /     
C       
C     MACHINE CONSTANTS FOR THE CDC 6000/7000 SERIES    
C       
C     DATA IMACH( 1) /          5 /     
C     DATA IMACH( 2) /          6 /     
C     DATA IMACH( 3) /          7 /     
C     DATA IMACH( 4) /    6LOUTPUT/     
C     DATA IMACH( 5) /         60 /     
C     DATA IMACH( 6) /         10 /     
C     DATA IMACH( 7) /          2 /     
C     DATA IMACH( 8) /         48 /     
C     DATA IMACH( 9) / 00007777777777777777B /  
C     DATA IMACH(10) /          2 /     
C     DATA IMACH(11) /         47 /     
C     DATA IMACH(12) /       -929 /     
C     DATA IMACH(13) /       1070 /     
C     DATA IMACH(14) /         94 /     
C     DATA IMACH(15) /       -929 /     
C     DATA IMACH(16) /       1069 /     
C       
C     MACHINE CONSTANTS FOR THE CELERITY C1260  
C       
C     DATA IMACH( 1) /          5 /     
C     DATA IMACH( 2) /          6 /     
C     DATA IMACH( 3) /          6 /     
C     DATA IMACH( 4) /          0 /     
C     DATA IMACH( 5) /         32 /     
C     DATA IMACH( 6) /          4 /     
C     DATA IMACH( 7) /          2 /     
C     DATA IMACH( 8) /         31 /     
C     DATA IMACH( 9) / Z'7FFFFFFF' /    
C     DATA IMACH(10) /          2 /     
C     DATA IMACH(11) /         24 /     
C     DATA IMACH(12) /       -126 /     
C     DATA IMACH(13) /        127 /     
C     DATA IMACH(14) /         53 /     
C     DATA IMACH(15) /      -1022 /     
C     DATA IMACH(16) /       1023 /     
C       
C     MACHINE CONSTANTS FOR THE CONVEX  
C     USING THE -fn COMPILER OPTION     
C       
C     DATA IMACH( 1) /          5 /     
C     DATA IMACH( 2) /          6 /     
C     DATA IMACH( 3) /          7 /     
C     DATA IMACH( 4) /          6 /     
C     DATA IMACH( 5) /         32 /     
C     DATA IMACH( 6) /          4 /     
C     DATA IMACH( 7) /          2 /     
C     DATA IMACH( 8) /         31 /     
C     DATA IMACH( 9) / 2147483647 /     
C     DATA IMACH(10) /          2 /     
C     DATA IMACH(11) /         24 /     
C     DATA IMACH(12) /       -127 /     
C     DATA IMACH(13) /        127 /     
C     DATA IMACH(14) /         53 /     
C     DATA IMACH(15) /      -1023 /     
C     DATA IMACH(16) /       1023 /     
C       
C     MACHINE CONSTANTS FOR THE CONVEX  
C     USING THE -fi COMPILER OPTION     
C       
C     DATA IMACH( 1) /          5 /     
C     DATA IMACH( 2) /          6 /     
C     DATA IMACH( 3) /          7 /     
C     DATA IMACH( 4) /          6 /     
C     DATA IMACH( 5) /         32 /     
C     DATA IMACH( 6) /          4 /     
C     DATA IMACH( 7) /          2 /     
C     DATA IMACH( 8) /         31 /     
C     DATA IMACH( 9) / 2147483647 /     
C     DATA IMACH(10) /          2 /     
C     DATA IMACH(11) /         24 /     
C     DATA IMACH(12) /       -125 /     
C     DATA IMACH(13) /        128 /     
C     DATA IMACH(14) /         53 /     
C     DATA IMACH(15) /      -1021 /     
C     DATA IMACH(16) /       1024 /     
C       
C     MACHINE CONSTANTS FOR THE CONVEX  
C     USING THE -p8 COMPILER OPTION     
C       
C     DATA IMACH( 1) /          5 /     
C     DATA IMACH( 2) /          6 /     
C     DATA IMACH( 3) /          7 /     
C     DATA IMACH( 4) /          6 /     
C     DATA IMACH( 5) /         32 /     
C     DATA IMACH( 6) /          4 /     
C     DATA IMACH( 7) /          2 /     
C     DATA IMACH( 8) /         63 /     
C     DATA IMACH( 9) / 9223372036854775807 /    
C     DATA IMACH(10) /          2 /     
C     DATA IMACH(11) /         53 /     
C     DATA IMACH(12) /      -1023 /     
C     DATA IMACH(13) /       1023 /     
C     DATA IMACH(14) /        113 /     
C     DATA IMACH(15) /     -16383 /     
C     DATA IMACH(16) /      16383 /     
C       
C     MACHINE CONSTANTS FOR THE CONVEX  
C     USING THE -pd8 COMPILER OPTION    
C       
C     DATA IMACH( 1) /          5 /     
C     DATA IMACH( 2) /          6 /     
C     DATA IMACH( 3) /          7 /     
C     DATA IMACH( 4) /          6 /     
C     DATA IMACH( 5) /         32 /     
C     DATA IMACH( 6) /          4 /     
C     DATA IMACH( 7) /          2 /     
C     DATA IMACH( 8) /         63 /     
C     DATA IMACH( 9) / 9223372036854775807 /    
C     DATA IMACH(10) /          2 /     
C     DATA IMACH(11) /         53 /     
C     DATA IMACH(12) /      -1023 /     
C     DATA IMACH(13) /       1023 /     
C     DATA IMACH(14) /         53 /     
C     DATA IMACH(15) /      -1023 /     
C     DATA IMACH(16) /       1023 /     
C       
C     MACHINE CONSTANTS FOR THE CRAY    
C     USING THE 46 BIT INTEGER COMPILER OPTION  
C       
C     DATA IMACH( 1) /        100 /     
C     DATA IMACH( 2) /        101 /     
C     DATA IMACH( 3) /        102 /     
C     DATA IMACH( 4) /        101 /     
C     DATA IMACH( 5) /         64 /     
C     DATA IMACH( 6) /          8 /     
C     DATA IMACH( 7) /          2 /     
C     DATA IMACH( 8) /         46 /     
C     DATA IMACH( 9) / 1777777777777777B /      
C     DATA IMACH(10) /          2 /     
C     DATA IMACH(11) /         47 /     
C     DATA IMACH(12) /      -8189 /     
C     DATA IMACH(13) /       8190 /     
C     DATA IMACH(14) /         94 /     
C     DATA IMACH(15) /      -8099 /     
C     DATA IMACH(16) /       8190 /     
C       
C     MACHINE CONSTANTS FOR THE CRAY    
C     USING THE 64 BIT INTEGER COMPILER OPTION  
C       
C     DATA IMACH( 1) /        100 /     
C     DATA IMACH( 2) /        101 /     
C     DATA IMACH( 3) /        102 /     
C     DATA IMACH( 4) /        101 /     
C     DATA IMACH( 5) /         64 /     
C     DATA IMACH( 6) /          8 /     
C     DATA IMACH( 7) /          2 /     
C     DATA IMACH( 8) /         63 /     
C     DATA IMACH( 9) / 777777777777777777777B / 
C     DATA IMACH(10) /          2 /     
C     DATA IMACH(11) /         47 /     
C     DATA IMACH(12) /      -8189 /     
C     DATA IMACH(13) /       8190 /     
C     DATA IMACH(14) /         94 /     
C     DATA IMACH(15) /      -8099 /     
C     DATA IMACH(16) /       8190 /     
C       
C     MACHINE CONSTANTS FOR THE DATA GENERAL ECLIPSE S/200      
C       
C     DATA IMACH( 1) /         11 /     
C     DATA IMACH( 2) /         12 /     
C     DATA IMACH( 3) /          8 /     
C     DATA IMACH( 4) /         10 /     
C     DATA IMACH( 5) /         16 /     
C     DATA IMACH( 6) /          2 /     
C     DATA IMACH( 7) /          2 /     
C     DATA IMACH( 8) /         15 /     
C     DATA IMACH( 9) /      32767 /     
C     DATA IMACH(10) /         16 /     
C     DATA IMACH(11) /          6 /     
C     DATA IMACH(12) /        -64 /     
C     DATA IMACH(13) /         63 /     
C     DATA IMACH(14) /         14 /     
C     DATA IMACH(15) /        -64 /     
C     DATA IMACH(16) /         63 /     
C       
C     MACHINE CONSTANTS FOR THE DEC ALPHA       
C     USING G_FLOAT     
C       
C     DATA IMACH( 1) /          5 /     
C     DATA IMACH( 2) /          6 /     
C     DATA IMACH( 3) /          5 /     
C     DATA IMACH( 4) /          6 /     
C     DATA IMACH( 5) /         32 /     
C     DATA IMACH( 6) /          4 /     
C     DATA IMACH( 7) /          2 /     
C     DATA IMACH( 8) /         31 /     
C     DATA IMACH( 9) / 2147483647 /     
C     DATA IMACH(10) /          2 /     
C     DATA IMACH(11) /         24 /     
C     DATA IMACH(12) /       -127 /     
C     DATA IMACH(13) /        127 /     
C     DATA IMACH(14) /         53 /     
C     DATA IMACH(15) /      -1023 /     
C     DATA IMACH(16) /       1023 /     
C       
C     MACHINE CONSTANTS FOR THE DEC ALPHA       
C     USING IEEE_FLOAT  
C       
C     DATA IMACH( 1) /          5 /     
C     DATA IMACH( 2) /          6 /     
C     DATA IMACH( 3) /          6 /     
C     DATA IMACH( 4) /          6 /     
C     DATA IMACH( 5) /         32 /     
C     DATA IMACH( 6) /          4 /     
C     DATA IMACH( 7) /          2 /     
C     DATA IMACH( 8) /         31 /     
C     DATA IMACH( 9) / 2147483647 /     
C     DATA IMACH(10) /          2 /     
C     DATA IMACH(11) /         24 /     
C     DATA IMACH(12) /       -125 /     
C     DATA IMACH(13) /        128 /     
C     DATA IMACH(14) /         53 /     
C     DATA IMACH(15) /      -1021 /     
C     DATA IMACH(16) /       1024 /     
C       
C     MACHINE CONSTANTS FOR THE DEC RISC
C       
C     DATA IMACH( 1) /          5 /     
C     DATA IMACH( 2) /          6 /     
C     DATA IMACH( 3) /          6 /     
C     DATA IMACH( 4) /          6 /     
C     DATA IMACH( 5) /         32 /     
C     DATA IMACH( 6) /          4 /     
C     DATA IMACH( 7) /          2 /     
C     DATA IMACH( 8) /         31 /     
C     DATA IMACH( 9) / 2147483647 /     
C     DATA IMACH(10) /          2 /     
C     DATA IMACH(11) /         24 /     
C     DATA IMACH(12) /       -125 /     
C     DATA IMACH(13) /        128 /     
C     DATA IMACH(14) /         53 /     
C     DATA IMACH(15) /      -1021 /     
C     DATA IMACH(16) /       1024 /     
C       
C     MACHINE CONSTANTS FOR THE DEC VAX 
C     USING D_FLOATING  
C       
C     DATA IMACH( 1) /          5 /     
C     DATA IMACH( 2) /          6 /     
C     DATA IMACH( 3) /          5 /     
C     DATA IMACH( 4) /          6 /     
C     DATA IMACH( 5) /         32 /     
C     DATA IMACH( 6) /          4 /     
C     DATA IMACH( 7) /          2 /     
C     DATA IMACH( 8) /         31 /     
C     DATA IMACH( 9) / 2147483647 /     
C     DATA IMACH(10) /          2 /     
C     DATA IMACH(11) /         24 /     
C     DATA IMACH(12) /       -127 /     
C     DATA IMACH(13) /        127 /     
C     DATA IMACH(14) /         56 /     
C     DATA IMACH(15) /       -127 /     
C     DATA IMACH(16) /        127 /     
C       
C     MACHINE CONSTANTS FOR THE DEC VAX 
C     USING G_FLOATING  
C       
C     DATA IMACH( 1) /          5 /     
C     DATA IMACH( 2) /          6 /     
C     DATA IMACH( 3) /          5 /     
C     DATA IMACH( 4) /          6 /     
C     DATA IMACH( 5) /         32 /     
C     DATA IMACH( 6) /          4 /     
C     DATA IMACH( 7) /          2 /     
C     DATA IMACH( 8) /         31 /     
C     DATA IMACH( 9) / 2147483647 /     
C     DATA IMACH(10) /          2 /     
C     DATA IMACH(11) /         24 /     
C     DATA IMACH(12) /       -127 /     
C     DATA IMACH(13) /        127 /     
C     DATA IMACH(14) /         53 /     
C     DATA IMACH(15) /      -1023 /     
C     DATA IMACH(16) /       1023 /     
C       
C     MACHINE CONSTANTS FOR THE ELXSI 6400      
C       
C     DATA IMACH( 1) /          5 /     
C     DATA IMACH( 2) /          6 /     
C     DATA IMACH( 3) /          6 /     
C     DATA IMACH( 4) /          6 /     
C     DATA IMACH( 5) /         32 /     
C     DATA IMACH( 6) /          4 /     
C     DATA IMACH( 7) /          2 /     
C     DATA IMACH( 8) /         32 /     
C     DATA IMACH( 9) / 2147483647 /     
C     DATA IMACH(10) /          2 /     
C     DATA IMACH(11) /         24 /     
C     DATA IMACH(12) /       -126 /     
C     DATA IMACH(13) /        127 /     
C     DATA IMACH(14) /         53 /     
C     DATA IMACH(15) /      -1022 /     
C     DATA IMACH(16) /       1023 /     
C       
C     MACHINE CONSTANTS FOR THE HARRIS 220      
C       
C     DATA IMACH( 1) /          5 /     
C     DATA IMACH( 2) /          6 /     
C     DATA IMACH( 3) /          0 /     
C     DATA IMACH( 4) /          6 /     
C     DATA IMACH( 5) /         24 /     
C     DATA IMACH( 6) /          3 /     
C     DATA IMACH( 7) /          2 /     
C     DATA IMACH( 8) /         23 /     
C     DATA IMACH( 9) /    8388607 /     
C     DATA IMACH(10) /          2 /     
C     DATA IMACH(11) /         23 /     
C     DATA IMACH(12) /       -127 /     
C     DATA IMACH(13) /        127 /     
C     DATA IMACH(14) /         38 /     
C     DATA IMACH(15) /       -127 /     
C     DATA IMACH(16) /        127 /     
C       
C     MACHINE CONSTANTS FOR THE HONEYWELL 600/6000 SERIES       
C       
C     DATA IMACH( 1) /          5 /     
C     DATA IMACH( 2) /          6 /     
C     DATA IMACH( 3) /         43 /     
C     DATA IMACH( 4) /          6 /     
C     DATA IMACH( 5) /         36 /     
C     DATA IMACH( 6) /          6 /     
C     DATA IMACH( 7) /          2 /     
C     DATA IMACH( 8) /         35 /     
C     DATA IMACH( 9) / O377777777777 /  
C     DATA IMACH(10) /          2 /     
C     DATA IMACH(11) /         27 /     
C     DATA IMACH(12) /       -127 /     
C     DATA IMACH(13) /        127 /     
C     DATA IMACH(14) /         63 /     
C     DATA IMACH(15) /       -127 /     
C     DATA IMACH(16) /        127 /     
C       
C     MACHINE CONSTANTS FOR THE HP 730  
C       
C     DATA IMACH( 1) /          5 /     
C     DATA IMACH( 2) /          6 /     
C     DATA IMACH( 3) /          6 /     
C     DATA IMACH( 4) /          6 /     
C     DATA IMACH( 5) /         32 /     
C     DATA IMACH( 6) /          4 /     
C     DATA IMACH( 7) /          2 /     
C     DATA IMACH( 8) /         31 /     
C     DATA IMACH( 9) / 2147483647 /     
C     DATA IMACH(10) /          2 /     
C     DATA IMACH(11) /         24 /     
C     DATA IMACH(12) /       -125 /     
C     DATA IMACH(13) /        128 /     
C     DATA IMACH(14) /         53 /     
C     DATA IMACH(15) /      -1021 /     
C     DATA IMACH(16) /       1024 /     
C       
C     MACHINE CONSTANTS FOR THE HP 2100 
C     3 WORD DOUBLE PRECISION OPTION WITH FTN4  
C       
C     DATA IMACH( 1) /          5 /     
C     DATA IMACH( 2) /          6 /     
C     DATA IMACH( 3) /          4 /     
C     DATA IMACH( 4) /          1 /     
C     DATA IMACH( 5) /         16 /     
C     DATA IMACH( 6) /          2 /     
C     DATA IMACH( 7) /          2 /     
C     DATA IMACH( 8) /         15 /     
C     DATA IMACH( 9) /      32767 /     
C     DATA IMACH(10) /          2 /     
C     DATA IMACH(11) /         23 /     
C     DATA IMACH(12) /       -128 /     
C     DATA IMACH(13) /        127 /     
C     DATA IMACH(14) /         39 /     
C     DATA IMACH(15) /       -128 /     
C     DATA IMACH(16) /        127 /     
C       
C     MACHINE CONSTANTS FOR THE HP 2100 
C     4 WORD DOUBLE PRECISION OPTION WITH FTN4  
C       
C     DATA IMACH( 1) /          5 /     
C     DATA IMACH( 2) /          6 /     
C     DATA IMACH( 3) /          4 /     
C     DATA IMACH( 4) /          1 /     
C     DATA IMACH( 5) /         16 /     
C     DATA IMACH( 6) /          2 /     
C     DATA IMACH( 7) /          2 /     
C     DATA IMACH( 8) /         15 /     
C     DATA IMACH( 9) /      32767 /     
C     DATA IMACH(10) /          2 /     
C     DATA IMACH(11) /         23 /     
C     DATA IMACH(12) /       -128 /     
C     DATA IMACH(13) /        127 /     
C     DATA IMACH(14) /         55 /     
C     DATA IMACH(15) /       -128 /     
C     DATA IMACH(16) /        127 /     
C       
C     MACHINE CONSTANTS FOR THE HP 9000 
C       
C     DATA IMACH( 1) /          5 /     
C     DATA IMACH( 2) /          6 /     
C     DATA IMACH( 3) /          6 /     
C     DATA IMACH( 4) /          7 /     
C     DATA IMACH( 5) /         32 /     
C     DATA IMACH( 6) /          4 /     
C     DATA IMACH( 7) /          2 /     
C     DATA IMACH( 8) /         31 /     
C     DATA IMACH( 9) / 2147483647 /     
C     DATA IMACH(10) /          2 /     
C     DATA IMACH(11) /         24 /     
C     DATA IMACH(12) /       -126 /     
C     DATA IMACH(13) /        127 /     
C     DATA IMACH(14) /         53 /     
C     DATA IMACH(15) /      -1015 /     
C     DATA IMACH(16) /       1017 /     
C       
C     MACHINE CONSTANTS FOR THE IBM 360/370 SERIES,     
C     THE XEROX SIGMA 5/7/9, THE SEL SYSTEMS 85/86, AND 
C     THE PERKIN ELMER (INTERDATA) 7/32.
C       
C     DATA IMACH( 1) /          5 /     
C     DATA IMACH( 2) /          6 /     
C     DATA IMACH( 3) /          7 /     
C     DATA IMACH( 4) /          6 /     
C     DATA IMACH( 5) /         32 /     
C     DATA IMACH( 6) /          4 /     
C     DATA IMACH( 7) /          2 /     
C     DATA IMACH( 8) /         31 /     
C     DATA IMACH( 9) /  Z7FFFFFFF /     
C     DATA IMACH(10) /         16 /     
C     DATA IMACH(11) /          6 /     
C     DATA IMACH(12) /        -64 /     
C     DATA IMACH(13) /         63 /     
C     DATA IMACH(14) /         14 /     
C     DATA IMACH(15) /        -64 /     
C     DATA IMACH(16) /         63 /     
C       
C     MACHINE CONSTANTS FOR THE IBM PC  
C       
C     DATA IMACH( 1) /          5 /     
C     DATA IMACH( 2) /          6 /     
C     DATA IMACH( 3) /          0 /     
C     DATA IMACH( 4) /          0 /     
C     DATA IMACH( 5) /         32 /     
C     DATA IMACH( 6) /          4 /     
C     DATA IMACH( 7) /          2 /     
C     DATA IMACH( 8) /         31 /     
C     DATA IMACH( 9) / 2147483647 /     
C     DATA IMACH(10) /          2 /     
C     DATA IMACH(11) /         24 /     
C     DATA IMACH(12) /       -125 /     
C     DATA IMACH(13) /        127 /     
C     DATA IMACH(14) /         53 /     
C     DATA IMACH(15) /      -1021 /     
C     DATA IMACH(16) /       1023 /     
C       
C     MACHINE CONSTANTS FOR THE IBM RS 6000     
C       
C     DATA IMACH( 1) /          5 /     
C     DATA IMACH( 2) /          6 /     
C     DATA IMACH( 3) /          6 /     
C     DATA IMACH( 4) /          0 /     
C     DATA IMACH( 5) /         32 /     
C     DATA IMACH( 6) /          4 /     
C     DATA IMACH( 7) /          2 /     
C     DATA IMACH( 8) /         31 /     
C     DATA IMACH( 9) / 2147483647 /     
C     DATA IMACH(10) /          2 /     
C     DATA IMACH(11) /         24 /     
C     DATA IMACH(12) /       -125 /     
C     DATA IMACH(13) /        128 /     
C     DATA IMACH(14) /         53 /     
C     DATA IMACH(15) /      -1021 /     
C     DATA IMACH(16) /       1024 /     
C       
C     MACHINE CONSTANTS FOR THE INTEL i860      
C       
C     DATA IMACH( 1) /          5 /     
C     DATA IMACH( 2) /          6 /     
C     DATA IMACH( 3) /          6 /     
C     DATA IMACH( 4) /          6 /     
C     DATA IMACH( 5) /         32 /     
C     DATA IMACH( 6) /          4 /     
C     DATA IMACH( 7) /          2 /     
C     DATA IMACH( 8) /         31 /     
C     DATA IMACH( 9) / 2147483647 /     
C     DATA IMACH(10) /          2 /     
C     DATA IMACH(11) /         24 /     
C     DATA IMACH(12) /       -125 /     
C     DATA IMACH(13) /        128 /     
C     DATA IMACH(14) /         53 /     
C     DATA IMACH(15) /      -1021 /     
C     DATA IMACH(16) /       1024 /     
C       
C     MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING   
C     32-BIT INTEGER ARITHMETIC.
C       
C     DATA IMACH( 1) /          5 /     
C     DATA IMACH( 2) /          6 /     
C     DATA IMACH( 3) /          5 /     
C     DATA IMACH( 4) /          6 /     
C     DATA IMACH( 5) /         32 /     
C     DATA IMACH( 6) /          4 /     
C     DATA IMACH( 7) /          2 /     
C     DATA IMACH( 8) /         31 /     
C     DATA IMACH( 9) / 2147483647 /     
C     DATA IMACH(10) /          2 /     
C     DATA IMACH(11) /         24 /     
C     DATA IMACH(12) /       -127 /     
C     DATA IMACH(13) /        127 /     
C     DATA IMACH(14) /         56 /     
C     DATA IMACH(15) /       -127 /     
C     DATA IMACH(16) /        127 /     
C       
C     MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING   
C     16-BIT INTEGER ARITHMETIC.
C       
C     DATA IMACH( 1) /          5 /     
C     DATA IMACH( 2) /          6 /     
C     DATA IMACH( 3) /          5 /     
C     DATA IMACH( 4) /          6 /     
C     DATA IMACH( 5) /         16 /     
C     DATA IMACH( 6) /          2 /     
C     DATA IMACH( 7) /          2 /     
C     DATA IMACH( 8) /         15 /     
C     DATA IMACH( 9) /      32767 /     
C     DATA IMACH(10) /          2 /     
C     DATA IMACH(11) /         24 /     
C     DATA IMACH(12) /       -127 /     
C     DATA IMACH(13) /        127 /     
C     DATA IMACH(14) /         56 /     
C     DATA IMACH(15) /       -127 /     
C     DATA IMACH(16) /        127 /     
C       
C     MACHINE CONSTANTS FOR THE SGI-IRIS
C       
C     DATA IMACH( 1) /          5 /     
C     DATA IMACH( 2) /          6 /     
C     DATA IMACH( 3) /          6 /     
C     DATA IMACH( 4) /          6 /     
C     DATA IMACH( 5) /         32 /     
C     DATA IMACH( 6) /          4 /     
C     DATA IMACH( 7) /          2 /     
C     DATA IMACH( 8) /         31 /     
C     DATA IMACH( 9) / 2147483647 /     
C     DATA IMACH(10) /          2 /     
C     DATA IMACH(11) /         24 /     
C     DATA IMACH(12) /       -125 /     
C     DATA IMACH(13) /        128 /     
C     DATA IMACH(14) /         53 /     
C     DATA IMACH(15) /      -1021 /     
C     DATA IMACH(16) /       1024 /     
C       
C     MACHINE CONSTANTS FOR THE SUN     
C       
C     DATA IMACH( 1) /          5 /     
C     DATA IMACH( 2) /          6 /     
C     DATA IMACH( 3) /          6 /     
C     DATA IMACH( 4) /          6 /     
C     DATA IMACH( 5) /         32 /     
C     DATA IMACH( 6) /          4 /     
C     DATA IMACH( 7) /          2 /     
C     DATA IMACH( 8) /         31 /     
C     DATA IMACH( 9) / 2147483647 /     
C     DATA IMACH(10) /          2 /     
C     DATA IMACH(11) /         24 /     
C     DATA IMACH(12) /       -125 /     
C     DATA IMACH(13) /        128 /     
C     DATA IMACH(14) /         53 /     
C     DATA IMACH(15) /      -1021 /     
C     DATA IMACH(16) /       1024 /     
C       
C     MACHINE CONSTANTS FOR THE SUN     
C     USING THE -r8 COMPILER OPTION     
C       
C     DATA IMACH( 1) /          5 /     
C     DATA IMACH( 2) /          6 /     
C     DATA IMACH( 3) /          6 /     
C     DATA IMACH( 4) /          6 /     
C     DATA IMACH( 5) /         32 /     
C     DATA IMACH( 6) /          4 /     
C     DATA IMACH( 7) /          2 /     
C     DATA IMACH( 8) /         31 /     
C     DATA IMACH( 9) / 2147483647 /     
C     DATA IMACH(10) /          2 /     
C     DATA IMACH(11) /         53 /     
C     DATA IMACH(12) /      -1021 /     
C     DATA IMACH(13) /       1024 /     
C     DATA IMACH(14) /        113 /     
C     DATA IMACH(15) /     -16381 /     
C     DATA IMACH(16) /      16384 /     
C       
C     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES FTN COMPILER 
C       
C     DATA IMACH( 1) /          5 /     
C     DATA IMACH( 2) /          6 /     
C     DATA IMACH( 3) /          1 /     
C     DATA IMACH( 4) /          6 /     
C     DATA IMACH( 5) /         36 /     
C     DATA IMACH( 6) /          4 /     
C     DATA IMACH( 7) /          2 /     
C     DATA IMACH( 8) /         35 /     
C     DATA IMACH( 9) / O377777777777 /  
C     DATA IMACH(10) /          2 /     
C     DATA IMACH(11) /         27 /     
C     DATA IMACH(12) /       -128 /     
C     DATA IMACH(13) /        127 /     
C     DATA IMACH(14) /         60 /     
C     DATA IMACH(15) /      -1024 /     
C     DATA IMACH(16) /       1023 /     
C       
C     MACHINE CONSTANTS FOR THE Z80 MICROPROCESSOR      
C       
C     DATA IMACH( 1) /          1 /     
C     DATA IMACH( 2) /          1 /     
C     DATA IMACH( 3) /          0 /     
C     DATA IMACH( 4) /          1 /     
C     DATA IMACH( 5) /         16 /     
C     DATA IMACH( 6) /          2 /     
C     DATA IMACH( 7) /          2 /     
C     DATA IMACH( 8) /         15 /     
C     DATA IMACH( 9) /      32767 /     
C     DATA IMACH(10) /          2 /     
C     DATA IMACH(11) /         24 /     
C     DATA IMACH(12) /       -127 /     
C     DATA IMACH(13) /        127 /     
C     DATA IMACH(14) /         56 /     
C     DATA IMACH(15) /       -127 /     
C     DATA IMACH(16) /        127 /     
C       
C***FIRST EXECUTABLE STATEMENT  I1MACH  
c ***start messing here***
c#if defined(CRAY)
c      DATA IMACH( 1) /   5 /                                          
c      DATA IMACH( 2) /   6 /                                          
c      DATA IMACH( 3) /   102 /                                          
c      DATA IMACH( 4) /   0 /                                          
c      DATA IMACH( 5) /    64 /                                          
c      DATA IMACH( 6) /     8 /                                          
c      DATA IMACH( 7) /     2 /                                          
c      DATA IMACH( 8) /    63 /                                          
c      DATA IMACH( 9) /  777777777777777777777B /                        
c      DATA IMACH(10) /     2 /                                          
c      DATA IMACH(11) /    47 /                                          
c      DATA IMACH(12) / -8189 /                                          
c      DATA IMACH(13) /  8190 /                                          
c      DATA IMACH(14) /    94 /                                          
c      DATA IMACH(15) / -8099 /                                          
c      DATA IMACH(16) /  8190 /                                          
c#else
c#if defined(__alpha) && defined(__osf__)
c      DATA IMACH( 1) /          5 /     
c      DATA IMACH( 2) /          6 /     
c      DATA IMACH( 3) /          6 /     
c      DATA IMACH( 4) /          0 /     
c      DATA IMACH( 5) /         32 /     
c      DATA IMACH( 6) /          4 /     
c      DATA IMACH( 7) /          2 /     
c      DATA IMACH( 8) /         31 /     
c      DATA IMACH( 9) / 2147483647 /     
c      DATA IMACH(10) /          2 /     
c      DATA IMACH(11) /         24 /     
c      DATA IMACH(12) /       -125 /     
c      DATA IMACH(13) /        128 /     
c      DATA IMACH(14) /         53 /     
c      DATA IMACH(15) /      -1021 /     
c      DATA IMACH(16) /       1024 /     
c#else
      DATA IMACH(1)  /   5/
      DATA IMACH(2)  /   6/
      DATA IMACH(3)  /   6/
      DATA IMACH(4)  /   0/
      DATA IMACH(5)  /  32/
      DATA IMACH(6)  /   4/
      DATA IMACH(7)  /   2/
      DATA IMACH(8)  /  31/
      DATA IMACH(9)  /2147483647/
      DATA IMACH(10) /   2/
      DATA IMACH(11) /  24/
      DATA IMACH(12) /-125/
      DATA IMACH(13) / 128/
      DATA IMACH(14) /  53/
      DATA IMACH(15) /-1021/
      DATA IMACH(16) / 1024/
c#endif       
c#endif       

      IF (I .LT. 1  .OR.  I .GT. 16) GO TO 10   
C       
      I1MACH = IMACH(I) 
      RETURN    
C       
   10 CONTINUE  
      WRITE (UNIT = OUTPUT, FMT = 9000) 
 9000 FORMAT ('1ERROR    1 IN I1MACH - I OUT OF BOUNDS')
C       
C     CALL FDUMP
C       
      STOP      
      END       
