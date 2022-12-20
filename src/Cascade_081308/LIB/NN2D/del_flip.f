c
c------------------------------------------------------------------------
c
c	del_flip - flips a triangulation into a Delaunay triangulation in 2-D
c
c	Input:
c		points(2,np)		list of node co-ordinates	
c	        v(3,*)			list of triangle vertices	
c	        e(3,*)			list of neighbouring triangles.	
c					E(i,j) is the triangle which
c					shares the face containing 
c					node mod((i+1),3) and mod(i+1),3)+1
c					in triangle j, stored counterclockwise 
c					about j.
c					(This is the `opposite' definition) 
c		nt			number of triangles in current
c					triangulation.
c		mask			a logical work array which stores
c					all triangles that have been flipped
c					at sometime during the calculation.
c		mask_e			a logical work array which stores
c					the list of edges on the stack.
c					This is used to prevent the same
c					edge being put on the stack twice.
c
c	Output:
c
c               v			updated
c               e			updated
c		nswaps			the number of swaps required to
c					flip into a Delaunay. If nswaps = 0
c					then original triangulation was 
c					Delaunay.
c
c	Comments:
c
c	The method uses a combination of a standard flipping algorithm
c	using a LIFO stack and a mask to record all triangles touched
c	during the stack.
c
c	This uses an algorithm that guarantees a Delaunay after one pass.
c
c	Converted to use an adjacency array defined using the `opposite'
c	definition 31/1/96.
c
c					M. Sambridge, RSES, Aug. 1995.
c
c------------------------------------------------------------------------
c
	Subroutine del_flip (points,e,v,nt,mask,mask_e,nswaps) 

	real*8		points(2,*)
 	integer		v(3,*)
 	integer		e(3,*)
	integer		erl,era,erb,erc,edg,v1,v2,v3,a,b,c,d,l,r
	integer		p,t,pos
	logical		swap
	logical		mask(*)
	logical		mask_e(3,*)
        integer		c1(3),c2(3)
        data		c1/2,3,1/
        data		c2/3,1,2/

c					initialize masks
	nswaps = 0
	do 10 t=1,nt
           mask(t) = .false.
           mask_e(1,t) = .false.
           mask_e(2,t) = .false.
           mask_e(3,t) = .false.
 10     continue
c					main loop over initial triangles
	do 20 t=1,nt

	   if(mask(t))go to 20

cd         write(*,*)' '
cd         write(*,*)' Starting stack from triangle',t
cd         write(*,*)' '

           v1 = v(1,t)
           v2 = v(2,t)
           v3 = v(3,t)
	   a = e(c2(1),t)
	   b = e(c2(2),t)
	   c = e(c2(3),t)
c					loop over edges of base triangle
           if(a.ne.0)then
              era = edg(a,t,e)
              era = c1(era)
           end if

           if(b.ne.0)then
              erb = edg(b,t,e)
              erb = c1(erb)
           end if

           if(c.ne.0)then
              erc = edg(c,t,e)
              erc = c1(erc)
           end if

cd         write(*,*)' i = 1',' t = ',t
cd         write(*,*)' v ',v(1,t),v(2,t),v(3,t)
cd         write(*,*)' e ',e(c2(1),t),e(c2(2),t),e(c2(3),t)
cd         write(*,*)' a b c ',a,b,c

c					initialize stackpair
 	   call stackpairinit

c					add edges of base triangle to stack
c
	   if(a.ne.0)then
              call pushpair(a,era)
              mask_e(era,a) = .true.
cd	      write(*,*)' pushed pair ',a,era,
cd   &                  ' nodes: ',v(era,a),v(c1(era),a)
           end if

	   if(b.ne.0)then
              call pushpair(b,erb)
              mask_e(erb,b) = .true.
cd	      write(*,*)' pushed pair ',b,erb,
cd   &                  ' nodes: ',v(erb,b),v(c1(erb),b)
           end if

	   if(c.ne.0)then
              call pushpair(c,erc)
              mask_e(erc,c) = .true.
cd	      write(*,*)' pushed pair ',c,erc,
cd   &                  ' nodes: ',v(erc,c),v(c1(erc),c)
           end if
c					check stack is not empty

           if(a.eq.0.and.b.eq.0.and.c.eq.0)go to 100

c					loop while stack is not empty
 50        continue
      
	   call poppair(r,pos)

cd	   write(*,*)' popped pair ',r,pos,
cd   &               ' nodes: ',v(pos,r),v(c1(pos),r)

	   l = e(c2(pos),r)
c					reset edge mask for edge l-r
           mask_e(pos,r) = .false.
c					check if new point is in circumcircle
c
           era = c1(pos)
           erb = c2(pos)
	   v1 = v(pos,r)
	   v2 = v(era,r)
	   v3 = v(erb,r)
	   erl=edg(l,r,e)
           erl = c1(erl)
           p = v(c2(erl),l)
      
	   if(swap(points(1,v1),points(1,v2),
     &             points(1,v3),points(1,p)))then

c					update triangle mask
	      mask(r) = .true.
	      mask(l) = .true.

cd            write(*,*)' triangles flipped =',r,l

	      nswaps = nswaps + 1
c					new point is inside circumcircle
c					for triangle r so swap diagonal
              a=e(c2(era),r)
              b=e(c2(erb),r)
              erc=c1(erl)
              c=e(c2(erc),l)
              d=e(c2(c2(erl)),l)
c					update adjacency list for triangle l
	      v(erc,l) = v3
	      e(c2(erc),l) = r
	      e(c2(erl),l)  = a

c					update adjacency list for triangle r
	      v(1,r)=p
	      v(2,r)=v3
	      v(3,r)=v1
	      e(c2(1),r)=l
	      e(c2(2),r)=b
	      e(c2(3),r)=c
c					put edges l-a and r-b on stack
c					update adjacency list for 
c					triangle a 
	      if(a.ne.0)then
	         pos = edg(a,r,e)
	         pos = c1(pos)
	         e(c2(pos),a)=l
                 if(.not.mask_e(pos,a))then
                    call pushpair(a,pos)
cd	            write(*,*)' pushed pair ',a,pos,
cd   &                        ' nodes: ',v(pos,a),v(c1(pos),a)
		 end if
	      end if

	      if(b.ne.0)then
	         pos = edg(b,r,e)
	         pos = c1(pos)
                 if(.not.mask_e(pos,b))then
                    call pushpair(b,pos)
cd	            write(*,*)' pushed pair ',b,pos,
cd   &                        ' nodes: ',v(pos,b),v(c1(pos),b)
		 end if
              end if
c					put edges l-d and r-c on stack
c					update adjacency list for 
c					triangle c
	      if(c.ne.0) then
                 pos = edg(c,l,e)
	         pos = c1(pos)
                 e(c2(pos),c)=r
                 if(.not.mask_e(pos,c))then
                    call pushpair(c,pos)
cd	            write(*,*)' pushed pair ',c,pos,
cd   &                        ' nodes: ',v(pos,c),v(c1(pos),c)
		 end if
	      end if

	      if(d.ne.0) then
                 pos = edg(d,l,e)
	         pos = c1(pos)
                 if(.not.mask_e(pos,d))then
                    call pushpair(d,pos)
cd	            write(*,*)' pushed pair ',d,pos,
cd   &              ' nodes: ',v(pos,d),v(c1(pos),d)
		 end if
	      end if
	   end if

	   call stackpairempty(k)
	   if(k.ne.1)go to 50
 100       continue
           call stackpairflush()
           
 20     continue
   
	return
	end
c
