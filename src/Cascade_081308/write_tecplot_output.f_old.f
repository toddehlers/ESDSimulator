c write_tecplot_output

      subroutine write_tecplot_output
     &     (h,h0,x,y,nnode,vertices,nt,prec,time,timeint,
     &	    dhfluvial,dhdiffusion,dhls,dt,temp,catch)

c subroutine to write out files in Tecplot format

c INPUT: h         = current topography
c        h0        = height of the bedrock-alluvium interface
c        x,y       = x- and y-nodal coordinates
c        nnode     = number of nodes
c        vertices  = triangles connectivity
c        nt        = number of triangles
c        prec      = precip
c        time      = current time

      common /vocal/ ivocal

      real     h(nnode),h0(nnode),x(nnode),y(nnode)
      real     prec(nnode)
      real     dhfluvial(nnode),dhdiffusion(nnode),dhls(nnode),dt
      integer  vertices(3,nt),timeint
      character*4 timenice
      integer  temp(nnode),catch(nnode)

c round current time to something nice for filename
      timeint=timeint+1
c      print *,'timeint: ',timeint
      write (timenice,'(I4.4)') timeint
c      if (timeint.ge.1000) then
c        write (timenice,'I4') timeint
c      else if (timeint.ge.100) then
c        write (timenice,'I3') timeint
c      else if (timeint.ge.10) then
c        write (timenice,'I2') timeint
c      else if (timeint.ge.0) then
c        write (timenice,'I1') timeint
c      else
c        print *,'Unusual timesteps for tecplot output.  Exiting...'
c        stop
c      endif

c open files - hard-coded for now (11/06)
      open (47,file='RUN1/tec/topo_tec_'//timenice//'.dat',
     &      status='unknown')
c Write header(s)
      write(47,*) 'TITLE = "Cascade model topography"'
      write(47,*) 'VARIABLES = "x" "y" "z" "node" "Precipitation"
     & "Fluvial Erosion Rate" "Diffusion Erosion Rate"
     & "Landslide Erosion Rate" "Total Erosion Rate"
     & "Catchment Color" "Catchment Number"'
      write(47,*) 'ZONE T = "Cascade"'
      write(47,*) 'n=',nnode,', e=',nt,', et=triangle, f=fepoint'

c Write topography
      do i=1,nnode
      write(47,*) x(i),y(i),h(i)/1000,i,prec(i),
     &		  dhfluvial(i)/dt,dhdiffusion(i)/dt,
     &            dhls(i)/dt,
     &            (dhfluvial(i)+dhdiffusion(i)+dhls(i))/dt,
     &            temp(i),catch(i)
      enddo

c Write connectivities
      do k=1,nt
      write (47,*) (vertices(i,k),i=1,3)
      enddo

c close files
      close(47)
      end
