c fluvial_erosion

      subroutine fluvial_erosion (xkf,xlf_BR,xlf_AL,
     &                            x,y,h,h0,hi,ndon,nnode,
     &                            surface,slope,length,
     &                            water,sediment,
     &                            ibucket,dt,fix,
     &                            dhh,
     &                            nb,nn,nbmax,
     &                            iorder,itype_node,
     &                            dhminfluv,dhmaxfluv,
     &                            nlake,
     &                            sea_level,outflux,ideposition,
     &                            width_c,thresh)

c subroutine to calculate fluvial erosion

c INPUT: xkf         = fluvial erosion constant []
c        xlf_BR      = bedrock erosion length scale [m]
c        xlf_AL      = alluvials erosion length scale [m?]
c        x,y         = x- and y-nodal positions [km]
c        h           = present topography [m]
c        h0          = bedrock-alluvion interface
c        hi          = initial topography (at time 0)
c        ndon        = donor array
c        nnode       = number of nodes
c        surface     = surface associated with each node, [km2]
c        slope       = slope associated with each nodal link (stream) [m/km]
c        length      = length associated with each nodal link (stream) [km]
c        ibucket     = working array
c        dt          = time step length [yr]
c        fix         = boundary conditoin array
c        dhh         = amount of material eroded/deposited over the time step
c        nb          = number of natural neighbours per node
c        nn          = list of natural neighbours
c        nbmax       = maximum of natural neighbours
c        iorder      = working array containing the proper sequence in node
c                      number in which to perform the river erosion operations
c        itype_node  = type associated to each node (see later in code)
c        dhminfluv   = minimum amount removed by river erosion
c        dhmaxfluv   = maximum amount removed by river erosion
c        nlake       = determines whether a node is part of a lake or not
c        sea_level   = sea level
c        ideposition = to prevent deposition by rivers (=0)
c        width_c     = coef controling channel width as fnct of discharge [sqrt(yr/m)]
c        width_ch    = effective channel width [km]
c        thresh      = theshold for discharge for channel formation [m km2/yr]

c OUTPUT:  several arrays are updated:
c        h           = new current topography [m]

c          the contribution from river incision to outflux is calculated

c          the following arrays are filled:
c        water       = water discharge at each point
c        sediment    = sediment load at each point

c subroutines called:
c - debug
c - find_order

c changed elevation update statements to include fixed BC's DS 11/18/1

      common /vocal/ ivocal

      real*4     x(nnode),y(nnode),h(nnode),h0(nnode),hi(nnode)
      real*4     xkf(nnode),xlf_BR(nnode)
      real*4     surface(nnode),dhh(nnode)
      real*4     water(nnode),sediment(nnode)
      real*4     slope(nnode),length(nnode)
      integer*4  ndon(nnode)
      integer*4  ibucket(nnode)
      real       fix(nnode)

      integer*4  nb(*),nn(nbmax,*)
      integer*4  iorder(nnode),itype_node(nnode)

      integer*4  nlake(nnode)

c sets all ibucket to 0
c sets all water to 1 and all sediment to 0
c note that water(i) can be different for every node; this is where a more
c  "complex" orographic model should be included...


      do i=1,nnode
       ibucket(i)=0
       sediment(i)=0.
       dhh(i)=0.
      enddo

      dhminfluv=0.
      dhmaxfluv=0.

c finds proper ordering (cascade algorithm)

      if (ivocal.eq.1) call debug ('find_order$',0)
      call find_order (ibucket,ndon,fix,iorder,
     &                 nnode,norder)
      if (ivocal.eq.1) call debug ('fluvial_erosion$',1)

c itype_node determines the type of node:
c   -1:   local minima next to (or on) a fixed boundary
c    0:   diffusion only
c    1:   channel because one of its parents was a channel

        do i=1,nnode
        itype_node(i)=0
        enddo

        do jorder=1,norder
          i=iorder(jorder)

          slop=-slope(i)
c special treatment for self donors
          if (ndon(i).eq.i) then
          dh=0.
          itype_node(i)=-2
          outflux=outflux+sediment(i)*(1.-fix(i))
          sediment(i)=0.
          water(i)=0.

c special treatment for nodes below sea level
          elseif (h(i).lt.sea_level) then
          dh=sediment(i)/surface(i)
          h(i)=h(i)+dh*real(fix(i))
          itype_node(i)=-1
          sediment(i)=0.

          else

c normal nodes
         if (itype_node(i).ne.0) itype_node(i)=1

         sedeqb=xkf(i)*slop*water(i)*dt

         width_ch=width_c*sqrt(water(i))

         if(water(i).lt.thresh) width_ch=999999999.

c special treatment for lake nodes
         if (nlake(i).eq.1) then
          dh=sediment(i)/surface(i)
          h(i)=h(i)+dh*real(fix(i))
          sediment(i)=0.
         elseif (sediment(i).ge.sedeqb.and.ideposition.eq.1) then
          dh=(sediment(i)-sedeqb)/surface(i)
          if (water(i).ne.0.) then
           dhmax=sediment(i)*length(i)/xkf(i)/water(i)/dt+
     &             h(ndon(i))-h(i)
          else
           dhmax=dh
          endif
          if (dh.le.dhmax) then
           sediment(i)=sedeqb
           h(i)=h(i)+dh*real(fix(i))
          else
           dh=dhmax
           sediment(i)=sediment(i)-dhmax*surface(i)
           h(i)=h(i)+dh*real(fix(i))
          endif

c erosion
         else

c first case: eroding bedrock only
          if (h(i).le.h0(i)) then
           dsed_ch=(sedeqb)*(length(i)/xlf_BR(i))
           dh=-(dsed_ch)/(width_ch*length(i))
           dsediment=-dh*surface(i)

           h(i)=h(i)+dh*real(fix(i))

           sediment(i)=sediment(i)+dsediment

c second case: eroding alluvions only
          else
           dsed_ch=(sedeqb)*(length(i)/xlf_AL)
           dh=-(dsed_ch)/(width_ch*length(i))
           if (dh.ge.h0(i)-h(i)) then
            dsediment=-dh*surface(i)

            h(i)=h(i)+dh*real(fix(i))

            sediment(i)=sediment(i)+dsediment

c third case: eroding alluvions and bedrock
           else
            dh1=h(i)-h0(i)
            dsed_al_ch=dh1*length(i)*width_ch
            dsed_ch=sedeq*(length(i)/xlf_BR(i))
            dh=-(dsed_ch)/(width_ch*length(i))
            dsediment=-dh*surface(i)
            dsed_al_cell=dh1*surface(i)

            h(i)=h0(i)+dh*real(fix(i))

            sediment(i)=sediment(i)+dsediment+dsed_al_cell
           endif
          endif

         endif

         dhminfluv=amin1(dhminfluv,dh)
         dhmaxfluv=amax1(dhmaxfluv,dh)
         if (ndon(i).ne.i) then
          water(ndon(i))=water(ndon(i))+water(i)
          sediment(ndon(i))=sediment(ndon(i))+sediment(i)
          if (slop.ne.0.) itype_node(ndon(i))=1

         endif
         endif

         dhh(i)=dh*real(fix(i))

c re-setting the height of the border nodes to zero, in case they have
c been changed within this routine
c         if (fix(i).eq.0) h(i)=0.


       enddo

      return
      end
