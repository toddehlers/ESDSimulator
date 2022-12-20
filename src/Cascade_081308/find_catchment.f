c find_catchment

      subroutine find_catchment (ndon,nwork,ncat,nsill,nempty,nlake,
     &                           nn,nb,nbmax,h,nnode,nself,nselfu,fix,
     &                           x,y,slope,length)

c subroutine to find the distribution of catchment areas, their sills
c and the distribution of "lakes"

c a catchment area is characterized by a minimum in topography

c every node is part of a catchment area
c ncat(i) is the node number of the minimum in which node i
c is draining

c nsill(i) is 0 unless i is a minimum (giving its name to a catchment area)
c in which case nsill(i) gives the number of the node at the sill of
c the catchment of name i

c nlake(i) is 0 if node i is not in a lake
c             1 if node i is in a lake

c INPUT: ndon      = donor array
c        nwork     = working array
c        nn        = neighbour array
c        nb        = number of neighbour per node
c        nbmax     = maximum number of neighbour per node
c        h         = present topography
c        nnode     = number of nodes
c        fix       = boundary condition array
c        x,y       = x- and y-nodal corrdinates 
c        slope      = slope associated to each node/stream segment
c        length     = length associated to each node/stream segment

c OUTPUT:  ndon is modified
c        ncat      = catchement array
c        nsill     = sill array
c        nempty    = array containing where a sill empties
c        nlake     = lake array
c        nself     = number of nodes presenting a problem (no output was found)
c        nselfu    = number of local minima that could not be resolved

      common /vocal/ ivocal

c subroutines called:
c NONE

      integer    ndon(nnode),nwork(nnode),ncat(nnode),nsill(nnode)
      integer    nlake(nnode),nempty(nnode)
      integer    nb(nnode),nn(nbmax,nnode)
      real       fix(nnode)
      real       h(nnode),x(nnode),y(nnode),slope(nnode),length(nnode)

      hmax=h(1)
        do i=1,nnode
        nlake(i)=0
        hmax=max(hmax,h(i))
        enddo

c9999  continue

        do i=1,nnode
        ncat(i)=0
        nwork(i)=i
        enddo

c finds catchment names

11111 nnn=0
c loop over all nodes
        do i=1,nnode
c if one has not yet found the cathcment number for node i
         if (ncat(i).eq.0) then
          nl=nwork(i)
          nd=ndon(nl)
          if (nd.eq.nl) then
           ncat(i)=nl
          else
           nwork(i)=nd
           nnn=nnn+1
          endif
         endif
        enddo
      if (nnn.ne.0) goto 11111

c check if a catchment area is touching a fixed node
c in which case there is no need to find a sill

9999  continue

        do i=1,nnode
        nwork(i)=0
        enddo
        do i=1,nnode
          if (fix(i).lt.0.5) then
          nc=ncat(i)
          nwork(nc)=1
          endif
        enddo

c finds sills
        do i=1,nnode
        nc=ncat(i)
          if (nwork(nc).eq.0) then
c look for the lowest point on the catchment boundary that has a donor below it
          hmin=hmax
          jjj=0
            do j=1,nnode
              if (ncat(j).eq.nc) then
                do k=1,nb(j)
                nnn=nn(k,j)
                  if (ncat(nnn).ne.nc .and.h(nnn).lt.h(j)) then
                    if (h(j).lt.hmin) then
                    hmin=h(j)
                    jjj=j
                    endif
                  endif
                enddo
              endif
            enddo
c sometimes one has to look at a second ring of nodes
c ie neighbours of neighbours
          if (jjj.eq.0) then
            do j=1,nnode
              if (ncat(j).eq.nc) then
                do k=1,nb(j)
                nnn=nn(k,j)
                  if (ncat(nnn).ne.nc) then
                    do kk=1,nb(nnn)
                    mmm=nn(kk,nnn)
                      if (ncat(mmm).ne.nc .and. h(mmm).lt.h(nnn)) then
                        if (h(nnn).lt.hmin) then
                        hmin=h(nnn)
                        jjj=nnn
                        endif
                      endif
                    enddo
                  endif
                enddo
              endif
            enddo
          endif
c otherwise give up...
          if (jjj.eq.0) then
          stop 'cannot find a sill...'
          endif
c a sill has been found in node jjj
c we now find the one of its neighbours that is in another catchment
c and has the steepest slope
          xslopemin=0.
            do k=1,nb(jjj)
            nnn=nn(k,jjj)
              if (ncat(nnn).ne.nc .and. nnn.ne.jjj) then
              xlength=sqrt((x(nnn)-x(jjj))**2+(y(nnn)-y(jjj))**2)
              if (xlength.eq.0) print*,nnn,jjj
              xslope=(h(nnn)-h(jjj))/xlength
                if (xslope.lt.xslopemin) then
                xslopemin=xslope
                ndon(nc)=nnn
                length(nc)=xlength
                slope(nc)=xslope
                endif
              endif
            enddo
c find lake nodes in this catchment
            do k=1,nnode
            if (ncat(k).eq.nc .and. h(k).lt.h(jjj)) nlake(k)=1
            enddo
c reset catchment names
            do k=1,nnode
            if (ncat(k).eq.nc) ncat(k)=ncat(ndon(nc))
            enddo
          goto 9999
          endif
        enddo

      return
      end
