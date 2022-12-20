      subroutine interpol3d (x,y,z,f,xi,yi,zi,fi,ni,in)

      implicit real*8 (a-h,o-z)

      real*8 xi(ni),yi(ni),zi(ni),fi(ni)
      real*8,dimension(:),allocatable :: ri,si,ti
      logical in

      allocate (ri(ni),si(ni),ti(ni))

      if (ni.eq.8) then
      ri=(/-1.,1.,1.,-1.,-1.,1.,1.,-1./)
      si=(/-1.,-1.,1.,1.,-1.,-1.,1.,1./)
      ti=(/-1.,-1.,-1.,-1.,1.,1.,1.,1./)
      else
      ri=(/0.,1.,0.,0.,1.,0./)
      si=(/0.,0.,1.,0.,0.,1./)
      ti=(/-1.,-1.,-1.,1.,1.,1./)
      endif

      r=0.
      s=0.
      t=0.

    1 continue
      rp=r
      sp=s
      tp=t
      sum1=0.
      sum2=0.
        do i=1,ni
        sum=(1.+s*si(i))*(1.+t*ti(i))*xi(i)
        sum1=sum1+sum
        sum2=sum2+sum*ri(i)
        enddo
      r=(8.*x-sum1)/sum2
      sum1=0.
      sum2=0.
        do i=1,ni
        sum=(1.+t*ti(i))*(1.+r*ri(i))*yi(i)
        sum1=sum1+sum
        sum2=sum2+sum*si(i)
        enddo
      s=(8.*y-sum1)/sum2
      sum1=0.
      sum2=0.
        do i=1,ni
        sum=(1.+r*ri(i))*(1.+s*si(i))*zi(i)
        sum1=sum1+sum
        sum2=sum2+sum*ti(i)
        enddo
      t=(8.*z-sum1)/sum2
      drst=(r-rp)**2+(s-sp)**2+(t-tp)**2
      if (drst.gt.1.e-12) goto 1

      f=0.
        do i=1,ni
        f=f+fi(i)*(1.+ri(i)*r)*(1.+si(i)*s)*(1.+ti(i)*t)/8.
        enddo

      in=.false.
      if ((t+1.)*(t-1.).le.0.) in=.true.

      deallocate (ri,si,ti)

      return
      end
