! Subroutine to output x, y, and z postions along with node number,
! temperature value, and x, y, and z velocities

      subroutine tec_mat_output (x,y,z,t,nnode,icon,nelem,mpe,run,&
                     c3,dt,Peclet,x1f,y1f,x2f,y2f,def,dif,time,&
                     xmin,xmax,ymin,ymax,zl,vx,vy,vz,geoflag,theta,phi,&
                     mft_ratein,mbt_ratein,mct_ratein,stf_ratein,xlonmin,&
                     xlatmin,xlonmax,xlatmax,nrun,Peclet2)

! This subroutine formats the output from Pecube so that it can easily
! be read into Tecplot or Matlab.

! Declare variables
      implicit real*8 (a-h,o-z)

      real*8 x(nnode),y(nnode),z(nnode),t(nnode),z_mod,counter,vx,vy,vz
!      real*8 eheat(nelem)
      integer istep,ie,k
      integer icon(mpe,nelem),nrun
      character run*100,c3*3,mpe_char*10

! Open files to create
      open(77,file=run(1:nrun)//'/Temps_tec'//c3//'.dat',status='unknown')
!      open(78,file=run//'/Temps_mat.dat',status='unknown')

! Write headers to files
      write (77,*) 'TITLE = "Pecube temperature output"'
      write (77,'(a100)') 'VARIABLES = "x (km)" "y (km)" "z (km)" "node" "temperature (C)" "U (mm/yr)" "V (mm/yr)" "W (mm/yr)"'
      !write (77,'(a200)') 'VARIABLES = "x (km)" "y (km)" "z (km)" "node" "temperature (C)" "U (mm/yr)" "V (mm/yr)" "W (mm/yr)" "heat production (<math>0</math>C/My)" "thermal diffusivity (km<sup>2</sup>/My)"'
      write (77,*) 'ZONE T = "Pecube"'
      write (77,*) 'n=',nnode,', e=',nelem,'et=brick,f=fepoint'

! Counter variable that replaced the variable 'i' when being printed
! to the output files
! Reason is that 'i' is an integer and to format the output file, all variables
! need to be floating point.  So, the counting variable was created
	counter=1.0

! Write temperatures, velocities, and positions to files
      do i=1,nnode
        call find_velo(x(i),y(i),z(i),dt,Peclet,x1f,y1f,x2f,y2f,def,dif,time,&
                       xmin,xmax,ymin,ymax,zl,vx,vy,vz,geoflag,theta,phi,&
                       mft_ratein,mbt_ratein,mct_ratein,stf_ratein,Peclet2)
        vx=Peclet*vx
        vy=Peclet*vy
        vz=Peclet*vz

! Calculates elevation from model thickness to node elevation
        z_mod=z(i)-zl
! Write thermal field data to file
        write (77,'(8f12.3)') x(i)+xlonmin,y(i)+xlatmin,z_mod,counter,t(i),vx,vy,vz
        counter=counter+1.0
      enddo

      write (mpe_char,'(i10)') mpe

! Write out node-element connectivities for tecplot version
      do ie=1,nelem
        write (77,'('//mpe_char//'i10)') (icon(k,ie),k=1,mpe)
      enddo

! Close files
      close(77)
      end
