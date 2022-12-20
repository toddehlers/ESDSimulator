c find greatest downstream slope for each node
c
      subroutine find_dslope(x,y,h,nbmax,nnode,nn,nb,dslope)
                           
       common /vocal/ ivocal

       real x(nnode),y(nnode),h(nnode)
       real dslope(nnode)

       integer nn(nbmax,nnode),nb(nnode)

       real s,dij

c iterative variables
       integer i,j,k

       pi = 3.1415926
c **********************************
c done with variable declaration
c **********************************
       
       do i=1,nnode

        dslope(i) = 0.

        do k=1,nb(i)

         j=nn(k,i)
         if (j.ne.i) then
          dij=(x(i)-x(j))**2.+(y(i)-y(j))**2.
          dij=1000.*sqrt(dij)
          s=(h(j)-h(i))/dij

          if (s.lt.dslope(i)) then
           dslope(i) = s
          endif

         endif
        enddo
       enddo

       return             

       end
