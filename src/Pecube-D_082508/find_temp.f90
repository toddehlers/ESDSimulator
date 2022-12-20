      subroutine find_temp (t,z,n,depth,tnow)

      implicit real*8 (a-h,o-z)

! interpolation routine

      real*8 t(n),z(n)

      tnow=t(n)
      if (depth.gt.z(n)) return

      tnow=t(1)
      if (depth.lt.z(1)) return

        do i=1,n-1
        x0=z(i)
        x1=z(i+1)
          if ((depth-x0)*(depth-x1).le.0.) then
          r=(depth-x0)/(x1-x0)
          tnow=t(i)+(t(i+1)-t(i))*r
          return
          endif
        enddo

      return
      end

