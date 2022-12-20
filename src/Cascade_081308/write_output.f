c write_output

      subroutine write_output
     &     (h,h0,x,y,fix,nnode,nnodemax,iadapt,ifirst,vertices,nt,
     &      ndon,param,nparam,water,dhfluvial,dhdiffusion,dhls,
     &      ncat,nlake,time,dt,prec,slope,dslope,timeint,tecflag)

c subroutine to write output to a series of files in ASCII

c INPUT: h         = current topography
c        h0        = height of the bedrock-alluvium interface
c        x,y       = x- and y-nodal coordinates
c        fix       = boundary conditions
c        nnode     = number of nodes
c        nnodemax  = maximum number of nodes
c        iadapt    = flag to allow dynamic remeshing (=0 no; =1 yes)
c        ifirst    = =1 means first zeroth time step (used to store
c                    initial geometry and parameters)
c        vertices  = triangles connectivity
c        nt        = number of triangles
c        ndon      = donor array
c        param     = erosion parameters
c        nparam    = number of erosion parameters
c        water     = discharge
c        slope     = slope [m/km]
c        dhfluvial = fluvial erosion increment over the time step
c        dhdiffusion= diffusion erosion increment over the time step
c        ncat      = catchment names
c        nlake     = lake flag
c        time      = current time
c        slope     = downstream slope found at beginning of time step
c        dslope    = downstream slope found after all erosion [m/m]

c subroutines called:
c NONE

      common /vocal/ ivocal

      real     h(nnode),h0(nnode),x(nnode),y(nnode)
      real     dhfluvial(nnode),dhdiffusion(nnode),dhls(nnode)
      real     param(nnodemax,nparam),water(nnode),prec(nnode)
      real     slope(nnode),dslope(nnode)
      integer  vertices(3,nt),ndon(nnode),ncat(nnode),nlake(nnode)
      real     fix(nnode)
      integer  timeint,tecflag
      integer  temp(nnode),catch(nnode)

c the 'temp' variable stores the catchment color value of each node for use
c in the tecplot output
c The 'catch' variable stores the catchment number for each node
  
c topography 
      write (7,*) 'TIME ',time,nnode     
        do i=1,nnode
        write (7,*) i,h(i),h(i)-h0(i)
        enddo

c connections
c write connections at all time steps DS 8/14/1
c        if (iadapt.eq.1 .or. ifirst.eq.1) then
      if (tecflag .eq. 0) then
        write (8,*) 'TIME ',time,nt
          do k=1,nt
          write (8,*) k,(vertices(i,k),i=1,3)
          enddo
c        endif
      endif
c donors
      if (tecflag .eq. 0) then
        if (ifirst.eq.0) then
        write (9,*) 'TIME ',time,nnode
          do i=1,nnode
          write (9,*) i,ndon(i)
          enddo
        endif
      endif

c geometry
c write geometry at all time steps DS 8/14/1
c        if (iadapt.eq.1 .or. ifirst.eq.1) then
        write (10,*) 'TIME ',time,nnode
          do i=1,nnode
          write (10,*) i,x(i),y(i),fix(i)
          enddo
c        endif

c parameters
      if (tecflag .eq. 0) then
      write (11,*) 'TIME ',time,nparam,nnode
        do i=1,nnode
        write (11,*) (param(i,k),k=1,nparam)
        enddo
      endif

c discharge
      if (tecflag .eq. 0) then
        if (ifirst.eq.0) then
        write (12,*) 'TIME ',time,nnode
          do i=1,nnode
          write (12,*) i,water(i)
          enddo
        endif
      endif

c precipitation - added by DW 10/06
      if (tecflag .eq. 0) then
        if (ifirst.eq.0) then
        write (77,*) 'TIME ',time,nnode
          do i=1,nnode
          write (77,*) i,prec(i)
          enddo
        endif
      endif

c erosion rate
      if (tecflag .eq. 0) then
        if (ifirst.eq.0) then
        write (13,*) 'TIME ',time,nnode
          do i=1,nnode
          write (13,*) i,dhfluvial(i)/dt,dhdiffusion(i)/dt,
     &                 dhls(i)/dt,
     &                 (dhfluvial(i)+dhdiffusion(i)+dhls(i))/dt
          enddo
        endif
      endif

c catchements - the 'temp' variable is used to store the catchment values
c of each node and then is sent to the write_tecplot_output.f function
c to be writin to the tecplot formatted output files
      if (tecflag .eq. 0) then
        if (ifirst.eq.0) then
        mcat=0
        write (14,*) 'TIME ',time,nnode
1111    lcat=0
          do i=1,nnode
            if (ncat(i).gt.0) then
            lcat=1
            ncat0=ncat(i)
            mcat=mcat+1
            write (14,*) 'Cat ',mcat
              do j=i,nnode
                if (ncat(j).eq.ncat0) then
                write (14,*) j
                ncat(j)=-ncat(j)
                catch(j)=mcat
                if (MOD(mcat,10) .eq. 0) then
                  temp(j)=10
                else
                  temp(j)=MOD(mcat,10)
                endif
                endif
              enddo
            endif
          if (lcat.eq.1) goto 1111
          enddo
          do i=1,nnode
          ncat(i)=-ncat(i)
          enddo
        endif
      endif

c lakes
      if (tecflag .eq. 0) then
        if (ifirst.eq.0) then
        write (15,*) 'TIME ',time,nnode
          do i=1,nnode
          if (nlake(i).eq.1) write (15,*) i
          enddo
        endif
      endif

c slope-area
      if (tecflag .eq. 0) then
        if (ifirst.eq.0) then
         write (17,*) 'TIME',time,nnode
         do i=1,nnode
          write (17,*) i,-dslope(i),water(i)
         enddo
        endif
      endif

c regolith
      
c      write (16,*) 'TIME ',time,nnode
c        do i=1,nnode
c        reg_thick(i) = h(i)-h0(i)
c       if (reg_thick(i).lt.1.E-8) reg_thick(i) = 1.E-8
c        write (16,*) i, h(i)-h0(i) 
c        enddo

c Make tecplot formatted output files (dwhipp 11/06)
      call write_tecplot_output (h,h0,x,y,nnode,vertices,nt,prec,time,
     &                           timeint,dhfluvial,dhdiffusion,dhls,dt,
     &                           temp,catch)

c Make topography files to be loaded into Pecube (dwhipp 02/07)
c      call write_pecube_topo (h,nnode,timeint,run_name,nrun_name)

      return
      end
