      subroutine interpolate (t,tp,z,zp,n)

      implicit real*8 (a-h,o-z)

! interpolation routine

      real*8 t(n),tp(n),z(n),zp(n)

      do ip=2,n-1
      xp=zp(ip)
        do i=1,n-1
        x0=z(i)
        x1=z(i+1)
          if ((xp-x0)*(xp-x1).le.0.) then
          r=(xp-x0)/(x1-x0)
          tp(ip)=t(i)+(t(i+1)-t(i))*r
          goto 1
          endif
        enddo
      print*,'problem in interpolation'
      stop
1     continue
      enddo

      return
      end
