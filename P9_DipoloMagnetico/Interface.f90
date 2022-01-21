
!**********************************************************************************************************************
SUBROUTINE circle(np,n,r,dr,pt,ptcon)
IMPLICIT NONE
!Variable declaration
INTEGER np,n,i,j,ptcon(np)
REAL*4 dr,th,dth,amax,amin,am,pi,r,thc
REAL*4 pt(np,2)

pi=ACOS(-1.0)

amax=0.4*dr	!max allowed separation between particles
amin=0.2*dr	!min allowed separation between particles
am=0.5*(amax+amin) !mean value of separation

n=INT(2.0*pi*r/am)	!Number of particles for the circle
WRITE(*,*)n
dth=2.0*pi/FLOAT(n)	!delta theta 

!Initial position of the particles
OPEN(1,FILE='pt.txt',STATUS='REPLACE')
DO i=1, n
  th=FLOAT(i-1)*dth
  pt(i,1)=r*COS(th)
  pt(i,2)=r*SIN(th)
  WRITE(1,'(2f12.8)')pt(i,1),pt(i,2)
END DO
WRITE(1,'(2f12.8)')pt(1,1),pt(1,2)
CLOSE(1)

!Set initial connectivity
DO i=1, n-1
  ptcon(i)=i+1
END DO
ptcon(n)=1
END SUBROUTINE

!***************************************************************************
SUBROUTINE move_particles(nx,ny,n,np,x0,y0,dx,dy,dt,ptcon,pt,x,y,xc,yc,u,v)
IMPLICIT NONE
INTEGER nx,ny,n,np,i,ix,iy,ptcon(np),m
REAL*4 x0,y0,dx,dy,dt,pt(np,2),xc(0:nx+1),yc(0:ny+1),u(0:nx,0:ny+1),v(0:nx+1,0:ny)
REAL*4 x(0:nx),y(0:ny)
REAL*4 up,vp,us,un,ve,vw
REAL*4 amax,d,x1,x2,y1,y2

!Identification of the control volume
DO i=1, n
  ix=INT((pt(i,1)-x0)/dx)+1
  iy=INT((pt(i,2)-y0)/dy)+1
  
  !Horizontal velocity interpolation
  IF ((ix .LT. 0) .OR. (ix .GT. nx)) THEN
    up=0.0
  ELSEIF (pt(i,2) .GT. yc(iy)) THEN
    un=u(ix-1,iy+1)+(pt(i,1)-x(ix-1))*(u(ix,iy+1)-u(ix-1,iy+1))/dx
    us=u(ix-1,iy)+(pt(i,1)-x(ix-1))*(u(ix,iy)-u(ix-1,iy))/dx
    up=us+(pt(i,2)-yc(iy))*(un-us)/dy
  ELSEIF (pt(i,2) .LE. yc(iy)) THEN
    un=u(ix-1,iy)+(pt(i,1)-x(ix-1))*(u(ix,iy)-u(ix-1,iy))/dx
    us=u(ix-1,iy-1)+(pt(i,1)-x(ix-1))*(u(ix,iy-1)-u(ix-1,iy-1))/dx
    up=us+(pt(i,2)-yc(iy-1))*(un-us)/dy
  END IF
  
  !Vertical velocity interpolation
  IF ((iy .LT. 0) .OR. (iy .GT. ny)) THEN
    vp=0.0
  ELSEIF (pt(i,1) .GT. xc(ix)) THEN
    ve=v(ix+1,iy-1)+(pt(i,2)-y(iy-1))*(v(ix+1,iy)-v(ix+1,iy-1))/dy
    vw=v(ix,iy-1)+(pt(i,2)-y(iy-1))*(v(ix,iy)-v(ix,iy-1))/dy
    vp=vw+(pt(i,1)-xc(ix))*(ve-vw)/dx
  ELSEIF (pt(i,1) .LE. xc(ix)) THEN
    ve=v(ix,iy-1)+(pt(i,2)-y(iy-1))*(v(ix,iy)-v(ix,iy-1))/dy
    vw=v(ix-1,iy-1)+(pt(i,2)-y(iy-1))*(v(ix-1,iy)-v(ix-1,iy-1))/dy
    vp=vw+(pt(i,1)-xc(ix-1))*(ve-vw)/dx
  END IF
 
  !Calculation of the new position of particles
  pt(i,1)=pt(i,1)+up*dt
  pt(i,2)=pt(i,2)+vp*dt
END DO

!To add new points to the interface
amax=0.5*dx
m=n
DO i=1, m
  x1=pt(i,1); x2=pt(ptcon(i),1)
  y1=pt(i,2); y2=pt(ptcon(i),2)
  !Distance between two connected points
  d=ABS(SQRT((x1-x2)**2+(y1-y2)**2))
  
  IF (d .GT. amax) THEN !If the points are too separated
    n=n+1	!Add a new point
    pt(n,1)=0.5*(pt(i,1)+pt(ptcon(i),1))	!x coordinate of the new point
    pt(n,2)=0.5*(pt(i,2)+pt(ptcon(i),2))	!y coordinate of the new point
    ptcon(n)=ptcon(i)	!Set new connection
    ptcon(i)=n		!Set new connection
  END IF
END DO
END SUBROUTINE
