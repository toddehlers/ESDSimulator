c subroutine generating stochastic landslides 
c added by B.Champel-Duffait, january 2000

c DS 11/01
c This should be a working and correct version of the code developed by
c Peter van der Beek.  Many corrections have been made to the code as well
c as simple conversions from m to km.  Check the pvdb/ directory for slightly
c more explanation and archived codes.  Do not base future codes off of the
c archived codes as they are unstable and incorrect.
c

c added nb to end of parameter list DS 8/13/1
c added dtc to parameter list DS 8/15/1
c
c added lsmeth to parameter list DS 11/14/1
c  lsmeth = 1, probabilistic sliding
c         = 2, threshold slope of pmax
c         = 3, threshold slope of pmax but no deposition
c

      subroutine landslide(x,y,h,nn2,nb2,fix,surface,time,
     &                     nnode,nbmax,pmax,smax,dhls,
     &                     cohes,rho,grav,xk1,distmax,
     &                     dt,seed,xk0,dtc,tt,cell,nb,lsmeth)
                           
       common /vocal/ ivocal

       real x(nnode),y(nnode),h(nnode)
       real fix(nnode),dhls(nnode)

       integer nn2(nbmax,nnode),nb2(nnode),nb(nnode)
       integer cell(nnode,nbmax,2)
       
       real surface(nnode)
       real smax(nnode)
       
       real proba
       integer alea
       integer cont
    
c time since last landslide at this node :
       real tt(nnode)
       
       real pi

c smax = max slope at a node
c pot  = array of potential initiation sites
c npot = no of potential initiation sites
c cop  = array of i's lowest neighbour
c enl  = array of nodes that will be removed
c ptbaiss = array of lowered nodes
c ptaire = array of points for area calculation

       integer pot(nnode),cop(nnode),enl(nnode),pot2(nnode),enl2(nnode)
       integer npot,ienl 
       real dij,dist,distj
       real havt,hapr,vol,vol2
       integer ptbaiss(nnode),test
       integer kbaiss
       real dhls_ave

c record number of successful depositions
       integer nsuccess

c hc = maximum stable height at a node
       real hc      

       pi=3.14159
       
c fill array of potential initiation sites
       npot=1

       do i=1,nnode
        smax(i)=0.
        cop(i)=i
        
        do k=1,nb2(i)-1

         j=nn2(k,i)
         if (j.ne.i) then
           dij=(x(i)-x(j))**2+(y(i)-y(j))**2
           dij=1000.*sqrt(dij)
           s=atan((h(j)-h(i))/dij)
         if (s.ge.smax(i)) then
            smax(i)=s
            cop(i)=j
         endif
         endif
        enddo
        
        if(smax(i).gt.pmax) then 
         pot(npot)=i
         npot=npot+1
        endif
       enddo 

c fill array of nodes to be removed

       ienl=1
       do k=1,npot-1
        i=pot(k)

        do kk=1,nb2(i)

         j=nn2(kk,i)

         do kkk=1,npot-1
          if(j.eq.pot(kkk).and.h(j).lt.h(i)) then
           enl(ienl)=i
           ienl=ienl+1
          endif
		 enddo
        enddo
       enddo    
               
       
c order enl from highest to lowest

       do iess=1,ienl-2
        ja=enl(iess)
        jb=iess
        do jess=iess,ienl-1
         if(ja.lt.enl(jess)) then
          ja=enl(jess)
          jb=jess
         endif
        enddo
        enl(jb)=enl(iess)
        enl(iess)=ja
       enddo
              
c eliminate double points in enl

       i0=1
       enl2(1)=enl(1)
       do iess=2,ienl-1
        if(enl(iess).ne.enl2(i0)) then
         enl2(i0+1)=enl(iess)
         i0=i0+1
        endif
       enddo
       do i=1,i0-1
        enl(i)=enl2(i)
       enddo
       ienl=i0   
       
c eliminate points that are too high

       k=1
       do i=1,npot-1
        test=0
        do j=1,ienl-1
         if(pot(i).eq.enl(j)) test=1
        enddo
        if(test.eq.0) then
         pot2(k)=pot(i)
         k=k+1
        endif
       enddo    
       npot=k
       do i=1,npot-1
        pot(i)=pot2(i)
       enddo
c now pot only has correct points

       dist=0.
       
c loop thru all potential slide sites
       nsuccess = 0
       do i=1,npot-1
       
c assign id of site
        j=pot(i)

        alea = 0

c probabilistic failure 
        if (lsmeth.eq.1) then

         if (cos(smax(j)-pmax).lt.1.) then
          hc=(4.*cohes/(rho*grav))
          hc=hc*sin(smax(j))*cos(pmax)
          hc=hc/(1.-cos(smax(j)-pmax))
          k=cop(j)
          do while(smax(cop(k)).gt.pmax)
           k=cop(k)
          enddo
          htemp=h(k)-h(j)

          proba=htemp/hc+xk0*tt(j)/dtc
          proba=min(0.999,proba)

          proba=exp((xk1*dt)*log(1.-proba))
c          proba=exp((xk1*dt/dtc)*log(1.-proba))
          proba=1.-proba
          proba=min(1.,proba)

c       check if landsliding takes place
          call alea0(alea,proba,seed)
         endif

c slope threshold failure
        elseif (lsmeth.eq.2) then
         alea = 1

        endif
         
        if(alea.ne.0) then
        
         vol=0.
        
c dist is the maximum length of the slide [m]
         dist = distmax
        
c lower points and calculate volume

c kbaiss: # of nodes in slide
         kbaiss=0
         dhls_ave = 0.

c find nodes in slide
         do jj=1,nnode
          if (jj.ne.j) then

           distj=(x(jj)-x(j))**2.+(y(jj)-y(j))**2.
           distj=1000.*sqrt(distj)

c use different failure planes for different sliding methods
           if (lsmeth.eq.1) then
            cou=distj*tan((pmax+smax(j))/2.)
           elseif (lsmeth.eq.2) then
            cou = distj*tan(pmax)
           endif

           if((distj.le.dist).and.(h(jj).gt.(h(j)+cou))) then

            kbaiss=kbaiss+1
c            print *,kbaiss,j,jj,distj,distmax,h(jj),h(jj)+cou

            ptbaiss(kbaiss)=jj

            havt=h(jj)

c use different failure planes for different sliding methods
            if (lsmeth.eq.1) then
             h(jj)=h(j)+distj*tan((pmax+smax(j))/2.)
            elseif (lsmeth.eq.2) then
             h(jj) = h(j) + distj*tan(pmax)
            endif

            hapr=h(jj)
            dhls(jj)=hapr-havt
            dhls_ave = dhls_ave + dhls(jj)

            call cvol(jj,nb,vol,havt,hapr,x,y,h,nbmax,nnode,cell)

           endif
          endif
         enddo
         vol=abs(vol)

c get average change in elevation
         dhls_ave = dhls_ave / real(kbaiss)
        
c calculate area
c  get surface area from map area
c  where (pmax+smax)/2 is the angle of the failure plane

         call caire(cell,ptbaiss,kbaiss,x,y,nnode,aire,nbmax,nb)
         if (lsmeth.eq.1) then
          aire=aire/cos((pmax+smax(j))/2.)
         elseif (lsmeth.eq.2) then
          aire = aire/cos(pmax)
         endif


c continue programme
         if(vol.gt.0.) then
        
c write out landslide info to file for probabilistic sliding
          if (lsmeth.eq.1) then
           write (16,999) time,x(j),y(j),aire,vol,dhls_ave,
     &     htemp/hc,xk0*tt(j)/dtc,proba
          endif

c distribute volume below point i
          tt(j)=0.
          cont=1
          nvpt=j
          nombre=0
          do while(cont.eq.1.and.nombre.lt.20)
           nombre=nombre+1
           i0=nvpt
        
           call etalt(i0,x,y,h,vol,vol2,nnode,nb2,nn2,nbmax,
     &                cont,nvpt,dhls,cell,nb)

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
        
       
c end of alea.ne.0 block
        endif 

c end of i=1,npot-1 block
       enddo           

c write out to check for slides
       if (nsuccess.ne.(npot-1)) then
        write (16,998) time,npot-1,nsuccess
       endif
       
       return             

998    format (E13.5,I10,I10)
999    format (E13.5,E13.5,E13.5,E13.5,E13.5,E14.5,E14.5,E14.5,E14.5)
       end
       
c******************************************************************************
c******************************************************************************
       subroutine alea0(alea,proba,seed)      
       
       integer alea

       alea=0
        
c      seed=rand(seed)
       call random(seed,1)
       if(seed.lt.proba) then   
        alea=1
       endif
       
       return
       end  
         
c******************************************************************************

       subroutine cvol(j,nb,vol,havt,hapr,x,y,h,nbmax,nnode,cell)

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

       subroutine etalt(i0,x,y,h,vol,vol2,nnode,nb2,nn2,
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

        call cvol(jmax2,nb,vol2,havt,hapr,x,y,h,nbmax,nnode,cell)

        if(vol2.gt.vol)then
         h(jmax2)=(hapr-havt)*vol/vol2+havt
         dhls(jmax2)=h(jmax2)-havt
         cont=0
        else
         h(jmax2)=hapr
         dhls(jmax2)=hapr-havt
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
       
c************************************************************************

       subroutine caire(cell,ptbaiss,kbaiss,x,y,nnode,aire,nbmax,nb)     

       real x(nnode),y(nnode)
       integer ptbaiss(nnode)
       integer nb(nnode)

       integer cell(nnode,nbmax,2)
       integer tri(nnode,3)
       integer ktri
       
       ktri=0

c get triangles connected to this node
       do k=1,kbaiss
        j=ptbaiss(k)

        do knb=1,nb(j)-1

         ktri = ktri + 1

c  sort vertices by id in descending order
         tri(ktri,1)=max0(j,cell(j,knb,1),cell(j,knb,2))
         if(tri(ktri,1).eq.j) then
          tri(ktri,3)=min0(cell(j,knb,1),cell(j,knb,2))
          tri(ktri,2)=max0(cell(j,knb,1),cell(j,knb,2))
         endif 
         if(tri(ktri,1).eq.cell(j,knb,1)) then
          tri(ktri,3)=min0(j,cell(j,knb,2))
          tri(ktri,2)=max0(j,cell(j,knb,2))
         endif
         if(tri(ktri,1).eq.cell(j,knb,2)) then
          tri(ktri,3)=min0(j,cell(j,knb,1))
          tri(ktri,2)=max0(j,cell(j,knb,1))
         endif

        enddo
       enddo


c take out double points
       k=1
       do while(k.le.ktri-1)
        k0=k+1
        do while(k0.le.ktri)
         if((tri(k,1).eq.tri(k0,1)).and.(tri(k,2).eq.tri(k0,2))) then
          if(tri(k,3).eq.tri(k0,3)) then 
           do kfin=k0,ktri-1
            tri(kfin,1)=tri(kfin+1,1)
            tri(kfin,2)=tri(kfin+1,2)
            tri(kfin,3)=tri(kfin+1,3)
           enddo
           ktri=ktri-1
          else
           k0=k0+1
          endif
         else
          k0=k0+1
         endif 
        enddo
        k=k+1
       enddo
       
c calculate area
       aire=0.
       do k=1,ktri
        call caire2(tri(k,1),tri(k,2),tri(k,3),x,y,nnode,aire)
       enddo              
                
      return
      end
       
c * * * * * * * * * * * * * * * * * * * * * 

       subroutine caire2(i,j,k,x,y,nnode,aire)
       
       real x(nnode),y(nnode)
       real aide
               
       aide=(y(j)-y(i))*(x(k)-x(i))+(x(i)-x(j))*(y(k)-y(i))

c divide by 2 to get area
       aide = aide/2.

       aire=aire+abs(aide)*1.e6
      
       return 

       end
