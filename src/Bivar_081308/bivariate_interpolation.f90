   program bivar

   implicit none

   real ( kind = 8 ),dimension(:),allocatable :: xd,yd,zd,xi,yi,zi,wk
   integer,dimension(:),allocatable :: iwk
   integer cnt,md,ncp,ndp,nip,timeint,nx,ny,i,j,ios,noutdir,counter
   real ( kind = 8 ) junk1,junk2,dx,dy,junkt2,x_val,y_val,x_orig,y_orig,i_float,n_float
   character*4 time,junkt
   CHARACTER :: line*1024,outdir*100

   integer(4) :: command_status

   md=1
   ncp=4
   ndp=0
   nip=0
   dx=0.
   dy=0.
   x_val=0.
   y_val=0.
   nx=0
   ny=0
   x_orig=0.
   y_orig=0.

   write (6,*) 'Bivar execution started'
   write (6,*)

   open (14,file='input/Bivar/Bivar.in',status='unknown',iostat=ios)
   if (ios.ne.0) then
     open (14,file='Bivar.in',status='unknown')
   endif

   open (55,status='scratch')
 1 read (14,'(a1024)',end=2) line
   if (line(1:1).ne.'$'.and. line(1:1).ne.' ') write (55,'(a)') line
   goto 1
 2 close (14)
   rewind (55)

   do i=1,100
     outdir(i:i)=' '
   enddo
   read (55,'(a100)') outdir
   do i=1,100
     if (outdir(i:i).ne.' ') noutdir=i
   enddo

   open (10,file='output/IceCascade/topography',status='unknown',iostat=ios)
   if (ios.ne.0) then
     open (10,file='topography',status='unknown')
   endif

   open (13,file='output/IceCascade/geometry',status='unknown',iostat=ios)
   if (ios.ne.0) then
     open (13,file='geometry',status='unknown')
   endif

   write (6,*) 'topography and geometry file opend'

   ! Read in number of points in x and y directions
   read (55,*) nx,ny
   write (6,*) 'number of points: ', nx, ny

   ! Read in the spacing for x and y nodes
   read (55,*) dx,dy
   write (6,*) 'spacing: ', dx, dy

   ! Read in x and y origin points
   read (55,*) x_orig,y_orig
   write (6,*) 'origin: ', x_orig, y_orig

!    print *, 'Enter number of points wanted in x direction: '
!      read *, nx
!    print *, 'Enter number of points wanted in y direction: '
!      read *, ny
!    print *, 'Enter spacing of points wanted in x direction (utm): '
!      read *, dx
!    print *, 'Enter spacing of points wanted in y direction (utm): '
!      read *, dy
!    print *, 'Enter x origin wanted (utm): '
!      read *, x_orig
!    print *, 'Enter y origin wanted (utm): '
!      read *, y_orig

   nip=nx*ny

   allocate (xi(nip),yi(nip),zi(nip))

   zi=0.

   x_val = (nx-1)*dx
   y_val = (ny-1)*dy

   cnt=1
   n_float=y_orig
   do while (n_float .le. y_val)
     i_float=x_orig
     do while (i_float .le. x_val)
       xi(cnt)=i_float/1000.
       yi(cnt)=n_float/1000.
       cnt=cnt+1
       i_float=i_float+dx
     enddo
     n_float=n_float+dy
   enddo

   timeint=0
! loop for input files
   do
     command_status = system("date")

     read (10,*,end=500) junkt,junkt2,ndp
     read (13,*,end=500) junkt,junkt2,junk1

     print *, 'Processing Time ',junkt2
     print *, 'ndp ', ndp

     allocate (xd(ndp),yd(ndp),zd(ndp))
     allocate (iwk(max(31,27+ncp)*ndp+nip),wk(8*ndp))

   do i=1,ndp
     read (10,*) junk1,zd(i),junk2
     read (13,*) junk1,xd(i),yd(i),junk2
   enddo

   print *, "call to interpolate"

!   print *, "md: ", md, ", ncp: ", ncp, ", ndp: ", ndp, ", nip: ", nip
!   call idbvip (md,ncp,ndp,xd,yd,zd,nip,xi,yi,zi,iwk,wk)
   call interpolate(ndp, nip, xd, yd, zd, xi, yi, zi)

   print *, "timeint: ", timeint

   write (time,'(I4.4)') timeint
   timeint=timeint+1

   ! contains just z values:
   open (11,file=outdir(1:noutdir)//'/topo_pecube_'//time//'.dat',status='unknown')

   ! contains x, y, z values:
   open (12,file=outdir(1:noutdir)//'/topo_tec_interp_'//time//'.dat',status='unknown')

   write (12,*) 'TITLE = "Bivariate Interpolation of Cascade Topography"'
   write (12,*) 'VARIABLES = "x" "y" "z"'
   write (12,*) 'ZONE T = "Bivar"'
   write (12,'(A2,i10)',advance="no") 'n=',nip
   write (12,'(A4,i10)',advance="no") ', e=',(nx-1)*(ny-1)
   write (12,*) ', et=quadrilateral, f=fepoint'
!    write (12,*) 'I=',nip,', J=1, K=1, ZONETYPE=Ordered'
!    write (12,*) 'DATAPACKING=POINT'
!    write (12,*) 'DT=(DOUBLE DOUBLE DOUBLE)'

    do i=1,nip
     write (11,*) zi(i)
     write (12,*) xi(i),yi(i),zi(i)
    enddo

    counter=0
    do i=1,ny-1
      do j=1,nx-1
        counter=counter+1
        if (mod(counter,nx).eq.0) counter=counter+1
        write (12,*) counter+nx,counter+nx+1,counter+1,counter
      enddo
    enddo

   close(11)
   close(12)

   deallocate (xd,yd,zd)
   deallocate (iwk,wk)

! end loop for input files
   enddo

500   close(10)
      close(13)

   deallocate (xi,yi,zi)

   write (6,*)
   write (6,*) 'Bivar execution ended'

   end
