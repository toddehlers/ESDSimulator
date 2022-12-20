      subroutine find_element (xdepth,ydepth,zdepth,tnow,x,y,z,t, &
                                 xsurf,ysurf,zsurf,ielsurf,neighbour, &
                                 iconsurf,icon,nelemsurf,nelem,nsurf,nz,nnode, &
                                 npe,mpe,interpolate,obflag)

      implicit real*8 (a-h,o-z)

      real*8 x(nnode),y(nnode),z(nnode),t(nnode)
      real*8 xsurf(nsurf),ysurf(nsurf),zsurf(nsurf)
      integer iconsurf(npe,nelemsurf),icon(mpe,nelem)
      integer neighbour(npe,nelemsurf)
      logical interpolate,in,obflag
      real*8,dimension(:),allocatable :: xi,yi,zi,fi

      allocate (xi(mpe),yi(mpe),zi(mpe),fi(mpe))

      ! Set logical interpolate to true
      ! This will be set to false below if the point of interest is not
      ! found within an element
      ! obflag will be set to true if point is out side or below base of model
      obflag=.FALSE.
      interpolate=.TRUE.

      ! Define eps and reset ielsurf to 1
      eps=tiny(eps)
      if (ielsurf.eq.0) ielsurf=1

1     continue
        ! Loop through all nodes in element surface (4 for brick)
        do k=1,npe
          ! Define i1 and i2 (connected nodes for surface)
          i1=iconsurf(k,ielsurf)
          i2=iconsurf(1+mod(k,npe),ielsurf)
          side=(xdepth-xsurf(i1))*(ysurf(i2)-ysurf(i1))-(ydepth-ysurf(i1))*(xsurf(i2)-xsurf(i1))
          if (side.gt.eps) then
            ielsurf=neighbour(k,ielsurf)
            if (ielsurf.eq.0) then
              obflag=.TRUE.
              interpolate=.FALSE.
              deallocate (xi,yi,zi,fi)
              return
            endif
            goto 1
          endif
        enddo

! we found the element

        do k=1,nz-1
        iel3D=(nz-1)*(ielsurf-1)+k
          do i=1,mpe
          ic=icon(i,iel3D)
          xi(i)=x(ic)
          yi(i)=y(ic)
          zi(i)=z(ic)
          fi(i)=t(ic)
          enddo
        call interpol3d (xdepth,ydepth,zdepth,tnow,xi,yi,zi,fi,mpe,in)
        if (in) then
          deallocate (xi,yi,zi,fi)
          return
        endif
        enddo

      if (zdepth.lt.0.) obflag=.TRUE.
      interpolate=.FALSE.
      deallocate (xi,yi,zi,fi)
      return

      end