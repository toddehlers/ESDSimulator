c build_a

      subroutine build_a (x1,x2,x3,y1,y2,y3,xk,alpha,ae1,ae2)

c this subroutine builds the conductivity and mass matrices
c (ae1 and ae2) stored in compact form (11,12,22,13,23,33) for
c one three-noded element of geometry (x1,y1),(x2,y2),(x3,y3)
c conductivity*dt=xk and alpha is the time integration parameter
c (0<alpha<1)

c INPUT: x1, x2, x3 = x-coordinates of the nodes of the element
c        y1, y2, y3 = y-coordinates of ...
c        xk         = conductivity averaged over triangle
c        alpha      = time integration parameter (0 = explicit)
c                                                (1 = implicit)

c OUTPUT: ae1 = elemental conductivity matrix
c         ae2 = mass matrix

c subroutines called:
c NONE

      common /vocal/ ivocal

      real   ae1(6),ae2(6)

      xkalpha=xk*alpha
      xk1malpha=xk*(1.-alpha)
      surf2=x1*y2+x2*y3+x3*y1-y1*x2-y2*x3-y3*x1
      surf=surf2/2.
      a1=(x2*y3-x3*y2)/surf2
      b1=(y2-y3)/surf2
      c1=(x3-x2)/surf2
      a2=(x3*y1-x1*y3)/surf2
      b2=(y3-y1)/surf2
      c2=(x1-x3)/surf2
      a3=(x1*y2-x2*y1)/surf2
      b3=(y1-y2)/surf2
      c3=(x2-x1)/surf2
      xint=(x1+x2+x3)/3.
      yint=(y1+y2+y3)/3.
      xn1=a1+b1*xint+c1*yint
      xn2=a2+b2*xint+c2*yint
      xn3=a3+b3*xint+c3*yint
      xn11=xn1*xn1
      xn12=xn1*xn2
      xn22=xn2*xn2
      xn13=xn1*xn3
      xn23=xn2*xn3
      xn33=xn3*xn3
      bc11=b1*b1+c1*c1
      bc12=b1*b2+c1*c2
      bc22=b2*b2+c2*c2
      bc13=b1*b3+c1*c3
      bc23=b2*b3+c2*c3
      bc33=b3*b3+c3*c3
      ae1(1)=(xn11+xkalpha*bc11)*surf
      ae1(2)=(xn12+xkalpha*bc12)*surf
      ae1(3)=(xn22+xkalpha*bc22)*surf
      ae1(4)=(xn13+xkalpha*bc13)*surf
      ae1(5)=(xn23+xkalpha*bc23)*surf
      ae1(6)=(xn33+xkalpha*bc33)*surf
      ae2(1)=(xn11-xk1malpha*bc11)*surf
      ae2(2)=(xn12-xk1malpha*bc12)*surf
      ae2(3)=(xn22-xk1malpha*bc22)*surf
      ae2(4)=(xn13-xk1malpha*bc13)*surf
      ae2(5)=(xn23-xk1malpha*bc23)*surf
      ae2(6)=(xn33-xk1malpha*bc33)*surf

      return
      end

c compared equations with ones from fem book.
c  confirmed ...
c  surf2,surf, all the a,b,c's xint,yint, xn1,xn2,xn3, bc's

