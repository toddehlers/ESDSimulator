c DS 11/17/1
c this is a simplified version of the landslide routine originally by 
c  van der Beek.  it uses a simple threshold slope criterion to determine
c  sliding.  it only fails one node above the initiation site so there
c  is only ever one node in the failure.  
c  this version should reduce the cost of the sliding routine and should 
c  produce the slope-area relationship shown in tucker and bras, 1999.
c

      subroutine landslide_simple(x,y,h,nn,nb,nn2,nb2,fix,surface,time,
     &                     nnode,nbmax,pmax,dhls,cell,sidex,sidey)
                           
       common /vocal/ ivocal

       real x(nnode),y(nnode),h(nnode)
       real fix(nnode),dhls(nnode)

       integer nn(nbmax,nnode),nb(nnode),nn2(nbmax,nnode),nb2(nnode)
       integer cell(nnode,nbmax,2)
       
       real surface(nnode)
       real pmax

       real pi

c pot  = array of potential initiation sites
c npot = no of potential initiation sites
       integer npot 
       real dij
       real havt,hapr,vol,vol2

c used for runout
       integer cont,nvpt,nombre

c check the dimensions of pts if you get seg faults
       integer pts(nnode*4,2),hi(nnode*4)
       real hp(nnode)

c record number of successful depositions
       integer nsuccess

       integer niter
       integer ideposit

c dh, minimum height change [m] required for slide to occur
       real diff,dh

c iterative variables
       integer i,j,n

c number of calls to find_unstable
       integer ctr

       pi = 3.1415926
       dh = 0.1
       ideposit = 0
c **********************************
c done with variable declaration
c **********************************

c fill array of potential initiation sites
       npot=0

c       print *,nnode,pmax
c initialize dhls and save previous height
       do i=1,nnode
        dhls(i) = 0.
        hp(i) = h(i)
       enddo

c find unstable sites
c       print *,npot
       call find_unstable(x,y,h,nnode,nbmax,nb,nn,pts,npot,dh,
     & pmax,1,sidex,sidey)
c       print *,npot
       ctr = 1
       
c loop thru all potential slide sites

       niter = 1
 1     nsuccess = 0
       do n=1,npot
       
c assign id of site
        i = pts(n,1)
        j = pts(n,2)

c check if slope has changed
        dij=(x(i)-x(j))**2.+(y(i)-y(j))**2.
        dij=1000.*sqrt(dij)
        s=atan((h(j)-h(i))/dij)
        diff = h(j) - (h(i)+dij*tan(pmax))

        if (s.gt.pmax) then

         vol=0.
        
         havt=h(j)
         h(j) = h(i) + dij*tan(pmax)
         hapr=h(j)
         dhls(j) = dhls(j) + havt-hapr

c deposit slide
         if (ideposit.eq.1) then

         call ccvol(j,nb,vol,havt,hapr,x,y,h,nbmax,nnode,cell)

         vol=abs(vol)

         if(vol.gt.0.) then
       
          cont=1

c initiation site
          nvpt=i

c number of deposition sites
          nombre=0

          do while(cont.eq.1.and.nombre.lt.20)
           nombre=nombre+1
           i0=nvpt
        
           call eetalt(i0,x,y,h,vol,vol2,nnode,nb2,nn2,nbmax,
     &                 cont,nvpt,dhls,cell,nb)

          enddo
          if(cont.eq.1.and.nombre.ge.20) then
c           write(16,*) 'impossible to distribute landslide volume'
c           write(16,*) ': last point ',nvpt
	      else
c	       print *,'successful in distributing landslide volume'
           nsuccess = nsuccess + 1
          endif
        
c end of vol.gt.0 block
         endif

c end of ideposit block
         endif
        
        endif 

       enddo           

c write out to check for slides
       if (npot.gt.0) then
c        write (16,998) time,npot,nsuccess
       endif

c check for unstable sites and go back to fail them

c       call find_unstable(x,y,h,nnode,nbmax,nb,nn,pts,npot,dh,
c     & pmax,0,sidex,sidey)

       npot = 0

c increment number of calls to find_unstable
       ctr = ctr + 1

       if (npot.gt.0) then
        niter = niter + 1
c        print *,time,': iteration ',niter,' in landslide_simple.f'

c stop if too many iterations
        if (niter.eq.100) then
         print *,'STOPPING in landslide_simple.f!',time,': '
         stop
        endif

        goto 1
       endif

c write out to file
c       if (ctr.gt.2) write (16,*) ctr,'calls to find_unstable at time',time
       if (ctr.gt.2) write (16,*),ctr,time

c reset boundary elevations to fix problems with diffusion
c       do i=1,nnode
c        if (fix(i).eq.0) then
c         h(i) = hp(i)
c        endif
c       enddo
        
       return             

997    format (E12.5,E12.5,E12.5)
998    format (E13.5,I10,I10)
999    format (E13.5,E13.5,E13.5,E13.5,E13.5,E14.5,E14.5,E14.5,E14.5)
       end
       
c******************************************************************************
c******************************************************************************

       subroutine find_unstable(x,y,h,nnode,nbmax,nb,nn,pts,npot,dh,
     & pmax,ifirst,sidex,sidey)

       integer nb(nnode),nn(nbmax,nnode),pts(nnode*4,2),npot
       real x(nnode),y(nnode),h(nnode),dh,pmax

       integer bdry
       integer i,j,k
       real s,diff,dij

       npot = 0

       do i=1,nnode

        do k=1,nb(i)

         bdry = 0

         j=nn(k,i)

c ignore points on bdry after first pass
         if (ifirst.eq.0) then
          if (x(j).eq.0.0.or.x(j).eq.sidex.or.
     &        y(j).eq.0.0.or.y(j).eq.sidey) bdry=1 
         endif

         if (j.ne.i.and.bdry.eq.0) then
          dij=(x(i)-x(j))**2+(y(i)-y(j))**2
          dij=1000.*sqrt(dij)
          s=atan((h(j)-h(i))/dij)
          diff = h(j) - (h(i)+dij*tan(pmax))
          if (s.ge.pmax.and.diff.gt.dh) then
           npot = npot + 1
           pts(npot,1) = i
           pts(npot,2) = j
          endif
         endif
        enddo
       enddo

       return
       end
       
c******************************************************************************

       subroutine ccvol(j,nb,vol,havt,hapr,x,y,h,nbmax,nnode,cell)

       integer nb(nnode)

       real x(nnode),y(nnode),h(nnode)
       integer cell(nnode,nbmax,2)
       
c this routine calculates the volume removed or added when a nodes
c elevation changes
       do k=1,nb(j)-1

        aid=(x(cell(j,k,2))-x(cell(j,k,1)))*(y(cell(j,k,1))-y(j))
        aid=aid+(y(cell(j,k,2))-y(cell(j,k,1)))*(x(j)-x(cell(j,k,1)))

c convert km2 to m2
        aid=abs(aid)*1.e6

c divide cross product by 2. to get area DS 8/13/1
        aid = aid/2. 

c quadrature formula??, is the 2 already factored into this equation??
        vol=vol+abs(hapr-havt)*(aid)/6.

       enddo 

        
       return
       end

c******************************************************************************

       subroutine eetalt(i0,x,y,h,vol,vol2,nnode,nb2,nn2,
     &                  nbmax,cont,nvpt,dhls,cell,nb)
       
       real x(nnode),y(nnode),h(nnode),dhls(nnode)
       integer nn2(nbmax,nnode),nb2(nnode),nb(nnode)
       integer cell(nnode,nbmax,2),cont
      
       pi=3.141592 
       pmax2=0.
       jmax2=i0
       do k=1,nb2(i0)
        j=nn2(k,i0)
        if (j.ne.i0) then
          distj=(x(i0)-x(j))**2+(y(j)-y(i0))**2

          distj=1000.*sqrt(distj)

          p=(h(i0)-h(j))/distj
        if(p.gt.pmax2)then
         pmax2=p
         jmax2=j
         distj2=distj
        endif
        endif
       enddo
       if(pmax2.gt.2.*pi/180.) then
        havt=h(jmax2)
        hapr=h(i0)-distj2*tan(2.*pi/180.)

        call ccvol(jmax2,nb,vol2,havt,hapr,x,y,h,nbmax,nnode,cell)

        if(vol2.gt.vol)then
         h(jmax2)=(hapr-havt)*vol/vol2+havt
         dhls(jmax2) = dhls(jmax2) + h(jmax2)-havt
         cont=0
        else
         h(jmax2)=hapr
         dhls(jmax2) = dhls(jmax2) + hapr-havt
         nvpt=jmax2
         vol=vol-vol2
        endif 
       else
        pmax2=0.
        jsup=i0
        do k=1,nb2(i0)
         j=nn2(k,i0)
         if (j.ne.i0) then
          distj=(x(i0)-x(j))**2+(y(j)-y(i0))**2

          distj=1000.*sqrt(distj)

          p=(h(j)-h(i0))/distj
          if(p.gt.pmax2)then
           pmax2=p
           jsup=j
           distj2=distj
          endif
         endif
        enddo
        nvpt=jsup
        cont=1
       endif 
        
       return
       end
