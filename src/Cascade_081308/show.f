c show

      subroutine show (x,y,z,surface,ndon,nb,nn,nbmax,ncat,nlake,nnode,
     &                 mode,zmin,zmax,istream,ilake,inum)

c subroutine to plot the a primary field with the option of
c rivers, lakes and catchments superimposed

c mode determines whether a primary field must be contoured
c 0 = none
c 1 = z is contoured using zmin, zmax
c 2 = z is contoured using dynamical range in z (zmin, zmax not used)
c 3 = catchments are shown

c istream = 1 streams are drawn
c ilake   = 1 lakes are drawn

c INPUT: x,y       = x- and y-nodal coordinates
c        z         = field to be contoured
c        surface   = surface attached to each node
c        ndon      = donor array
c        nb        = number of neighbours per node
c        nn        = neighbour list (per node)
c        nbmax     = maximum number of neighbours per node
c        ncat      = catchment array
c        nlake     = lake array
c        nnode     = number of nodes
c        mode      = see above
c        zmin,zmax = see above
c        istream   = see above
c        ilake     = see above
c        inum      = 1 the first time the subroutine is called and the
c                    plot window is created
c                  = 0 the last time the subroutine is called and the
c                    user is prompted for mouse action before the
c                    window is destroyed
c                  > 1 the window stays up but the use is not prompted for
c                    mouse action

c subroutines called:
c - xopen
c - xcmap
c - xscale
c - xplot
c - xcirclef
c - xclear
c - xsave
c - xchoice
c - xclose
c - xpen

      real    x(*),y(*),z(*),surface(*)
      integer ndon(*),ncat(*),nlake(*)
      integer nb(*),nn(nbmax,*)

      integer window(4),red(230),green(230),blue(230)
      integer red1(6),green1(6),blue1(6)

      data window /100,100,800,800/
      data red,green,blue /230*0,230*0,230*0/
      data red1 /0,0,0,255,255,255/
      data green1 /0,255,255,255,0,255/
      data blue1 /255,255,0,0,0,255/

      if (mode.eq.2) then
      zzmin=z(1)
      zzmax=z(1)
        do i=1,nnode
        zzmin=min(zzmin,z(i))
        zzmax=max(zzmax,z(i))
        enddo
      elseif (mode.eq.1) then
      zzmin=zmin
      zzmax=zmax
      endif

      xmin=x(1)
      xmax=x(1)
      ymin=y(1)
      ymax=y(1)
        do i=1,nnode
        xmin=min(xmin,x(i))
        xmax=max(xmax,x(i))
        ymin=min(ymin,y(i))
        ymax=max(ymax,y(i))
        enddo

      if (inum.eq.1) then
      call xopen (window,char(0))
        do k=1,5
          do j=1,35
          fact=float(j)/35.
          red((k-1)*35+j+30)=red1(k)+fact*(red1(k+1)-red1(k))
          green((k-1)*35+j+30)=green1(k)+fact*(green1(k+1)-green1(k))
          blue((k-1)*35+j+30)=blue1(k)+fact*(blue1(k+1)-blue1(k))
          enddo
        enddo
      call xcmap (red,green,blue)
      call xscale (0,window(3),0,window(4),xmin,xmax,ymax,ymin,1)
      else
      call xclear
      endif

      pi=3.141592654
      if (mode.eq.1 .or. mode.eq.2) then
        do i=1,nnode
        radius=sqrt(surface(i)/pi)
        ipen=30+175*(z(i)-zzmin)/(zzmax-zzmin)
        ipen=min(ipen,215)
        ipen=max(ipen,31)
        call xpen (ipen)
        call xcirclef (x(i),y(i),radius)
        enddo
      elseif (mode.eq.3) then
      mcat=0
1111  lcat=0
        do i=1,nnode
          if (ncat(i).gt.0) then
          lcat=1
          ncat0=ncat(i)
          mcat=mcat+1
          ipen=31+mod(mcat*37,175)
          call xpen (ipen)
            do j=i,nnode
              if (ncat(j).eq.ncat0) then
                if (surface(j).le.0.) then
                print*,surface(j)
                stop
                endif
              radius=sqrt(surface(j)/pi)
              call xcirclef (x(j),y(j),radius)
              ncat(j)=-ncat(j)
              endif
            enddo
          endif
        if (lcat.eq.1) goto 1111
        enddo
        do i=1,nnode
        ncat(i)=-ncat(i)
        enddo
      endif

      if (istream.eq.1) then
      call xpen (1)
        do i=1,nnode
        call xplot (x(i),y(i),3)
        nd=ndon(i)
        call xplot (x(nd),y(nd),2)
        enddo
      endif

      if (ilake.eq.1) then
      call xpen (1)
        do i=1,nnode
          if (nlake(i).eq.1) then
          radius=sqrt(surface(i)/pi)
          call xcirclef (x(i),y(i),radius)
          endif
        enddo
      endif

      call xsave

        if (inum.eq.0) then
        ichoice=xchoice('Continue            ',1)
        call xclose
        endif

      return
      end
