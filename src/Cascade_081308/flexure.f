c flexure

      subroutine flexure (x,y,h,h0,hi,hiso,nnode,flex,work,nflex,
     &                    surface,points,vertices,neighbour,eps,
     &                    hflex,rhocflex,thickflex,ym,pratio,rhoaflex,
     &                    hisomin,hisomax,ixflex,iyflex,fix)

c subroutine to calculate the flexural response due to loading/unloading
c by deposition/erosion on the landscape
c This routine uses a spectral method to solve the thin elastic plate
c equation. It therefore assumes that the plate has uniform elastic
c properties.
c Note that the "load" applied to the thin elastic plate is calculated
c from the eroded/deposited material since time t=0. The initial topography
c does not constitue a load in itself.

c INPUT: x,y        = x- and y- coordinates of the nodes
c        h          = present topography
c        h0         = initial topography
c        hi         = bedrock-alluvials interface
c        hiso       = amount of isistatic uplift
c        nnode      = number of nodes
c        flex       = working array
c        work       = working array
c        nflex      = size of FFT (must be a power of 2)
c        surface    = surface area attached to each node
c        points     = working array
c        vertices   = triangle list
c        neighbour  = neighbour list
c        eps        = precision
c        hflex      = length of the thin elastic plate (where it is pinned)
c        rhocflex   = density of crustal rocks (in kg/m3)
c        thickflex  = elastic thickness (in km)
c        ym         = young modulus (in Pa)
c        pratio     = poisson's ratio
c        rhoaflex   = asthenospheric density (in kg/m3)
c        hisomin    = minimum isostatic deflection
c        hisomax    = maximum isostatic deflection
c        ixflex     = flag to turn x-isostasy on
c        iyflex     = flag to turn y-isostasy on
c        fix        = boundary conditions

c OUTPUT: in output several arrays are updated:
c        h,hi,hiso 

c subroutines called
c - debug
c - Triloc_del
c - sinft

      common /vocal/ ivocal

      real*4    x(nnode),y(nnode),h(nnode)
      real*4    h0(nnode),hiso(nnode),hi(nnode)
      real*4    flex(nflex,nflex),work(nflex,nflex)
      real*4    surface(nnode)
      real*8    xfl,yfl
      integer   vertices(3,*)
      integer	neighbour(3,*)
      real      fix(*)
      real*8    points(2,*),eps
      logical   out

      hx=hflex*1.e3
      drho=rhocflex*9.81*1.e3
      elt=thickflex*1.e3
      dflex=ym/12./(1.-pratio**2)
      d=dflex*elt**3
      xk=rhoaflex*9.81

c iref=0 reference height is sea-level
c iref=1 reference height is h0

      iref=1
      jref=1-iref

      xmin=x(1)
      xmax=x(1)
      ymin=y(1)
      ymax=y(1)
        do i=1,nnode
        xmin=amin1(x(i),xmin)
        xmax=amax1(x(i),xmax)
        ymin=amin1(y(i),ymin)
        ymax=amax1(y(i),ymax)
        points(1,i)=x(i)
        points(2,i)=y(i)
        enddo

      xmean=(xmin+xmax)/2.
      xfmin=xmean*1.e3-hx/2.
      xfmax=xmean*1.e3+hx/2.
      dxf=xfmax-xfmin
      ddxf=dxf/(nflex-1)
      ymean=(ymin+ymax)/2.
      yfmin=ymean*1.e3-hx/2.
      yfmax=ymean*1.e3+hx/2.
      dyf=yfmax-yfmin
      ddyf=dyf/(nflex-1)


        do j=1,nflex
        do i=1,nflex
        flex(i,j)=0.
        enddo
        enddo

      iflexmin=1+(xmin*1.e3-xfmin)/ddxf
      iflexmax=1+(xmax*1.e3-xfmin)/ddxf
      jflexmin=1+(ymin*1.e3-yfmin)/ddyf
      jflexmax=1+(ymax*1.e3-yfmin)/ddyf

      loc=1
      if (ivocal.eq.1) call debug ('Triloc_del$',0)
        do j=jflexmin,jflexmax
        yloc=((j-1)*ddyf+yfmin)/1.e3
        yfl=yloc
        do i=iflexmin,iflexmax
        xloc=((i-1)*ddxf+xfmin)/1.e3
        xfl=xloc
        call Triloc_del (xfl,yfl,points,vertices,neighbour,
     &                   loc,eps,out,kface,iface)
        if (out) then
        flex(i,j)=0.
        else
        i1=vertices(1,loc)
        i2=vertices(2,loc)
        i3=vertices(3,loc)
        x1=x(i1)
        x2=x(i2)
        x3=x(i3)
        y1=y(i1)
        y2=y(i2)
        y3=y(i3)
        a1=x2*y3-x3*y2
        a2=x3*y1-x1*y3
        a3=x1*y2-x2*y1
        b1=y2-y3
        b2=y3-y1
        b3=y1-y2
        c1=x3-x2
        c2=x1-x3
        c3=x2-x1
        a=x1*y2+x2*y3+x3*y1-y1*x2-y2*x3-y3*x1
        h1=(a1+b1*xloc+c1*yloc)/a
        h2=(a2+b2*xloc+c2*yloc)/a
        h3=(a3+b3*xloc+c3*yloc)/a
        hh=h(i1)*h1+h(i2)*h2+h(i3)*h3
        hh0=(hi(i1)*iref-hiso(i1)*iref)*h1
     &     +(hi(i2)*iref-hiso(i2)*iref)*h2
     &     +(hi(i3)*iref-hiso(i3)*iref)*h3
        flex(i,j)=-(hh-hh0)*drho*ddxf*ddyf
        endif
        enddo
        enddo
      if (ivocal.eq.1) call debug ('flexure$',1)

      if (ixflex.eq.0) then
        do j=1,nflex
        fmean=0.
          do i=iflexmin,iflexmax
          fmean=fmean+flex(i,j)
          enddo
        fmean=fmean/(iflexmax-iflexmin+1)
          do i=1,nflex
          flex(i,j)=fmean
          enddo
        enddo
      endif

      if (iyflex.eq.0) then
        do i=1,nflex
        fmean=0.
          do j=jflexmin,jflexmax
          fmean=fmean+flex(i,j)
          enddo
        fmean=fmean/(jflexmax-jflexmin+1)
          do j=1,nflex
          flex(i,j)=fmean
          enddo
        enddo
      endif

      if (ivocal.eq.1) call debug ('sinft$',0)
        do j=1,nflex
        call sinft (flex(1,j),nflex)
        enddo
      if (ivocal.eq.1) call debug ('flexure$',1)

        do j=1,nflex
        do i=1,nflex
        work(j,i)=flex(i,j)
        enddo
        enddo

      if (ivocal.eq.1) call debug ('sinft$',0)
        do i=1,nflex
        call sinft (work(1,i),nflex)
        enddo
      if (ivocal.eq.1) call debug ('flexure$',1)

        do j=1,nflex
        do i=1,nflex
        work(j,i)=work(j,i)*4./hx/hx
        enddo
        enddo

      pi=3.141592654
      pihx=pi/hx
        do j=1,nflex
        fj=(j*pihx)**2*iyflex
          do i=1,nflex
          fi=(i*pihx)**2*ixflex
          tij=d/xk*(fi**2+2.*fi*fj+fj**2)+1.
          work(j,i)=work(j,i)/xk/tij
          enddo
        enddo

      if (ivocal.eq.1) call debug ('sinft$',0)
        do i=1,nflex
        call sinft (work(1,i),nflex)
        enddo
      if (ivocal.eq.1) call debug ('flexure$',1)

        do j=1,nflex
        do i=1,nflex
        flex(i,j)=work(j,i)
        enddo
        enddo

      if (ivocal.eq.1) call debug ('sinft$',0)
        do j=1,nflex
        call sinft (flex(1,j),nflex)
        enddo
      if (ivocal.eq.1) call debug ('flexure$',1)

      if (ixflex.eq.0) then
      iflexfirst=1+(nflex-1)*(xmin*1.e3-xfmin)/dxf
        do j=1,nflex
        flex(iflexfirst,j)=flex(iflexfirst+1,j)
        enddo
      iflexlast=1+(nflex-1)*(xmax*1.e3-xfmin)/dxf
        do j=1,nflex
        flex(iflexlast,j)=flex(iflexlast-1,j)
        enddo
      endif

      if (iyflex.eq.0) then
      jflexfirst=1+(nflex-1)*(ymin*1.e3-yfmin)/dyf
        do i=1,nflex
        flex(i,jflexfirst)=flex(i,jflexfirst+1)
        enddo
      jflexlast=1+(nflex-1)*(ymax*1.e3-yfmin)/dyf
        do i=1,nflex
        flex(i,jflexlast)=flex(i,jflexlast-1)
        enddo
      endif


      hisomax=0.
      hisomin=0.
        do i=1,nnode
        xp=x(i)*1.e3
        yp=y(i)*1.e3
        iflex=1+(nflex-1)*(xp-xfmin)/dxf
        jflex=1+(nflex-1)*(yp-yfmin)/dyf
        xflex=(iflex-1)*ddxf+xfmin
        rx=(xp-xflex)/dxf*2.-1.
        yflex=(jflex-1)*ddyf+yfmin
        ry=(yp-yflex)/dyf*2.-1.
        hisotot=(flex(iflex,jflex)*(1.-rx)*(1.-ry)/4.
     &          +flex(iflex+1,jflex)*(1.+rx)*(1.-ry)/4.
     &          +flex(iflex+1,jflex+1)*(1.+rx)*(1.+ry)/4.
     &          +flex(iflex,jflex+1)*(1.-rx)*(1.+ry)/4.)/1.e3
        dh=hisotot-hiso(i)
        hiso(i)=hisotot
        hisomin=amin1(hisomin,dh)
        hisomax=amax1(hisomax,dh)
        h(i)=h(i)+dh
        h0(i)=h0(i)+dh
        hi(i)=hi(i)+dh
        enddo

      return
      end 
