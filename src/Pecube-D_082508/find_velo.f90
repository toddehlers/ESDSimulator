! Subroutine that calls function to calculate velocities in the x, y, and z directions 

      subroutine find_velo (x,y,z,dt,Peclet,x1f,y1f,x2f,y2f,def,dif,time,&
                              xmin,xmax,ymin,ymax,zl,vx,vy,vz,geoflag,theta,phi,&
                              mft_ratein,mbt_ratein,mct_ratein,stf_ratein,Peclet2)

! this extra layer between Pecube and the velocity function geometry ensures
! that the advection of the "lagrangian particles" for which we track the thermal
! history (to calculate the age) is second order accurate; it uses a mid-point
! algorithm.
      implicit real*8 (a-h,o-z)

      xx0=x
      yy0=y
      zz0=z
      xx=xx0
      yy=yy0
      zz=zz0
      niter=0


1     call geometry (xx,yy,zz,vx,vy,vz,time,xmin,xmax,ymin,ymax,zl,x1f,y1f,&
                     x2f,y2f,def,dif,geoflag,theta,phi,mft_ratein,&
                     mbt_ratein,mct_ratein,stf_ratein,Peclet,Peclet2)
      xxn=xx0+dt*Peclet*vx/2.
      yyn=yy0+dt*Peclet*vy/2.
      zzn=zz0+dt*Peclet*vz/2.


      xnorm=sqrt((xxn-xx)**2+(yyn-yy)**2+(zzn-zz)**2)
      xx=xxn
      yy=yyn
      zz=zzn
      niter=niter+1
! In some cases (low resolution mesh) this mid-point algorithm cannot converge
! we assume a simple explicit estimate for the velocity
        if (niter.gt.10) then
        call geometry (xx0,yy0,zz0,vx,vy,vz,time,xmin,xmax,ymin,ymax,zl,&
                       x1f,y1f,x2f,y2f,def,dif,geoflag,theta,phi,mft_ratein,&
                       mbt_ratein,mct_ratein,stf_ratein,Peclet,Peclet2)
        return
        endif
      if (xnorm.gt.zl*1.e-6) goto 1

      return
      end
