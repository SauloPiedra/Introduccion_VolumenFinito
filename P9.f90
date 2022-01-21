PROGRAM Dipole
IMPLICIT NONE
!Declaración de variables
CHARACTER*100 itchar
INTEGER i,j,nx,ny,maxiter_S,itmax,ei,ej,it,maxiter,c, n, np
REAL*4 dx,dy,x0,xl,y0,yl,dv,Se,Sw,Sn,Ss,tolerance,residual,time,dt,ue,uw,vn,vs,Div,Re,gamma,dbkx,dbky
REAL*4 xi, yi, Q, r, ai,bi,ci
REAL*4, ALLOCATABLE :: aE(:,:),aW(:,:),aN(:,:),aS(:,:),aP(:,:),Sp(:,:),u(:,:),u1(:,:),v(:,:),v1(:,:)
REAL*4, ALLOCATABLE :: x(:),xc(:),y(:),yc(:),P(:,:),Pp(:,:),d_h(:,:),d_v(:,:),uc(:,:),vc(:,:)
REAL*4, ALLOCATABLE :: B(:,:), pt(:,:)
INTEGER, ALLOCATABLE :: ptcon(:)

!Inicialización de variables
nx=200; ny=200
x0=-10.0; xl=10.0
y0=-10.0; yl=10.0
xi=0.0; yi=0.0
Re=1.0
Q=10.0
n=0
r=2.
!Magnet dimension in [m]
ai=1.e-2  
np=20000
time=100.0
dt=0.001
maxiter_S=100
maxiter=1000
tolerance=1e-5

ALLOCATE(aE(nx,ny),aW(nx,ny),aN(nx,ny),aS(nx,ny),aP(nx,ny),Sp(nx,ny),P(0:nx+1,0:ny+1),Pp(0:nx+1,0:ny+1))
ALLOCATE(u(0:nx,0:ny+1),u1(0:nx,0:ny+1),v(0:nx+1,0:ny),v1(0:nx+1,0:ny),x(0:nx),xc(0:nx+1),y(0:ny),yc(0:ny+1))
ALLOCATE(uc(0:nx+1,0:ny+1),vc(0:nx+1,0:ny+1),d_h(0:nx+1,0:ny+1),d_v(0:nx+1,0:ny+1))
ALLOCATE (B(0:nx+1,0:ny+1),pt(np,2),ptcon(np))

gamma=1.0/Re
dx=(xl-x0)/FLOAT(nx)
dy=(yl-y0)/FLOAT(ny)
itmax=INT(time/dt)+1
dv=dx*dy
Se=dy; Sw=dy
Sn=dx; Ss=dx

CALL MESH_1D(nx,x0,xl,x,xc)
CALL MESH_1D(ny,y0,yl,y,yc)

!Condiciones de frontera e iniciales
u=0.0; v=0.0; P=0.0; Pp=0.0

d_h=0.0; d_v=0.0

! u(:,ny+1)=1.0		!Deslizamiento en cara NORTE

u1=u; v1=v		!u1 y v1 son al tiempo actual

CALL B0z_Field(nx,ny,ai,ai*xc,ai*yc,B)
B=B/MAXVAL(ABS(B))
CALL WriteScalarField2D('B',0,B,xc,yc,nx,ny)

CALL circle(np,n,r,dx,pt,ptcon)

!archivo de animación
OPEN(1,FILE='anim.gnp',STATUS='REPLACE')
WRITE(1,*)'set size square; set xrange[-10:10]; set yrange[-10:10]; unset key'
WRITE(1,*)'set style data filledcurves'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!INICIA CICLO TEMPORAL

DO it=1,itmax
  Div=1.0
  c=1
  
  !INICIA CICLO SIMPLEC!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  DO WHILE ((c .LT. maxiter_S) .AND. (Div .GT. tolerance))
  aE=0.0; aW=0.0; aN=0.0; aS=0.0; aP=0.0; Sp=0.0
  
  !ECUACIÓN COMPONENTE U #######################################################################################
  ei=UBOUND(u,1)-1
  ej=UBOUND(u,2)-1
  
  !Cálculo de coeficientes
  DO i=1,ei
    DO j=1,ej
      !Flujos en las caras
      ue=(u(i,j)+u(i+1,j))*0.5
      uw=(u(i,j)+u(i-1,j))*0.5
      vn=(v(i,j)+v(i+1,j))*0.5
      vs=(v(i,j-1)+v(i+1,j-1))*0.5
      
      aE(i,j)=gamma*Se/dx-0.5*ue*Se
      aW(i,j)=gamma*Sw/dx+0.5*uw*Sw
      aN(i,j)=gamma*Sn/dy-0.5*vn*Sn
      aS(i,j)=gamma*Ss/dy+0.5*vs*Ss
      aP(i,j)=aE(i,j)+aW(i,j)+aN(i,j)+aS(i,j)+dv/dt
      Sp(i,j)=u(i,j)*dv/dt-(P(i+1,j)-P(i,j))*dv/dx
    END DO
  END DO
  
  !Corrección por condiciones de frontera
  dbkx=0.0; dbky=1.0
  
  !ESTE
  aP(ei,1:ej)=aP(ei,1:ej)+dbkx*aE(ei,1:ej)
  Sp(ei,1:ej)=Sp(ei,1:ej)+(1+dbkx)*aE(ei,1:ej)*u(ei+1,1:ej)
  aE(ei,1:ej)=0.0
  
  !OESTE
  aP(1,1:ej)=aP(1,1:ej)+dbkx*aW(1,1:ej)
  Sp(1,1:ej)=Sp(1,1:ej)+(1+dbkx)*aw(1,1:ej)*u(0,1:ej)
  aw(1,1:ej)=0.0
  
  !NORTE
  aP(1:ei,ej)=aP(1:ei,ej)+dbky*aN(1:ei,ej)
  Sp(1:ei,ej)=Sp(1:ei,ej)+(1+dbky)*aN(1:ei,ej)*u(1:ei,ej+1)
  aN(1:ei,ej)=0.0
  
  !SUR
  aP(1:ei,1)=aP(1:ei,1)+dbky*aS(1:ei,1)
  Sp(1:ei,1)=Sp(1:ei,1)+(1+dbky)*aS(1:ei,1)*u(1:ei,0)
  aS(1:ei,1)=0.0
  
  !Cálculo de coeficientes para corrección de la presión
  d_h(1:ei,1:ej)=Se/(aP(1:ei,1:ej)-(aE(1:ei,1:ej)+aW(1:ei,1:ej)+aN(1:ei,1:ej)+aS(1:ei,1:ej)))
  
  !Solución del sistema de ecuaciones
  CALL Gauss_TDMA2D(u1,ei,ej,aP,aE,aW,aN,aS,sP,nx,ny,maxiter,tolerance,residual)
!###########################################################################################################
  aP=0.0; aE=0.0; aW=0.0; aN=0.0; aS=0.0; Sp=0.0 
  !ECUACIÓN COMPONENTE V

  ei=UBOUND(v,1)-1
  ej=UBOUND(v,2)-1
  
  !Cálculo de coeficientes
  DO i=1,ei
    DO j=1,ej
      !Flujos en las caras
      ue=(u(i,j)+u(i,j+1))*0.5
      uw=(u(i-1,j)+u(i-1,j+1))*0.5
      vn=(v(i,j)+v(i,j+1))*0.5
      vs=(v(i,j)+v(i,j-1))*0.5
      
      aE(i,j)=gamma*Se/dx-0.5*ue*Se
      aW(i,j)=gamma*Sw/dx+0.5*uw*Sw
      aN(i,j)=gamma*Sn/dy-0.5*vn*Sn
      aS(i,j)=gamma*Ss/dy+0.5*vs*Ss
      aP(i,j)=aE(i,j)+aW(i,j)+aN(i,j)+aS(i,j)+dv/dt
      Sp(i,j)=v(i,j)*dv/dt-(P(i,j+1)-P(i,j))*dv/dy-Q*0.5*(B(i,j)+B(i,j+1))*dv
    END DO
  END DO
  
  !Corrección por condiciones de frontera
  dbkx=1.0; dbky=0.0
  
  !ESTE
  aP(ei,1:ej)=aP(ei,1:ej)+dbkx*aE(ei,1:ej)
  Sp(ei,1:ej)=Sp(ei,1:ej)+(1+dbkx)*aE(ei,1:ej)*v(ei+1,1:ej)
  aE(ei,1:ej)=0.0
  
  !OESTE
  aP(1,1:ej)=aP(1,1:ej)+dbkx*aW(1,1:ej)
  Sp(1,1:ej)=Sp(1,1:ej)+(1+dbkx)*aW(1,1:ej)*v(0,1:ej)
  aw(1,1:ej)=0.0
  
  !NORTE
  aP(1:ei,ej)=aP(1:ei,ej)+dbky*aN(1:ei,ej)
  Sp(1:ei,ej)=Sp(1:ei,ej)+(1+dbky)*aN(1:ei,ej)*v(1:ei,ej+1)
  aN(1:ei,ej)=0.0
  
  !SUR
  aP(1:ei,1)=aP(1:ei,1)+dbky*aS(1:ei,1)
  Sp(1:ei,1)=Sp(1:ei,1)+(1+dbky)*aS(1:ei,1)*v(1:ei,0)
  aS(1:ei,1)=0.0
  
  !Cálculo de coeficientes para corrección de la presión
  d_v(1:ei,1:ej)=Sn/(aP(1:ei,1:ej)-(aE(1:ei,1:ej)+aW(1:ei,1:ej)+aN(1:ei,1:ej)+aS(1:ei,1:ej)))
  
  !Solución del sistema de ecuaciones
  CALL Gauss_TDMA2D(v1,ei,ej,aP,aE,aW,aN,aS,sP,nx,ny,maxiter,tolerance,residual)
  
!##########################################################################################################

  !CORRECIÓN DE LA PRESIÓN
aP=0.0; aE=0.0; aW=0.0; aN=0.0; aS=0.0; Sp=0.0; Pp=0.0
  ei=UBOUND(Pp,1)-1
  ej=UBOUND(Pp,2)-1
  
  !Cálculo de coeficientes
  DO i=1,ei
    DO j=1,ej
      !Flujos en las caras
      ue=u1(i,j)
      uw=u1(i-1,j)
      vn=v1(i,j)
      vs=v1(i,j-1)
      
      aE(i,j)=d_h(i,j)*Se
      aW(i,j)=d_h(i-1,j)*Sw
      aN(i,j)=d_v(i,j)*Sn
      aS(i,j)=d_v(i,j-1)*Ss
      aP(i,j)=aE(i,j)+aW(i,j)+aN(i,j)+aS(i,j)
      Sp(i,j)=-ue*Se+uw*Sw-vn*Sn+vs*Ss
    END DO
  END DO
  
  !Solución del sistema de ecuaciones
  CALL Gauss_TDMA2D(Pp,ei,ej,aP,aE,aW,aN,aS,sP,nx,ny,maxiter,tolerance,residual)
  
  !Correción de la presión
  P=P+Pp
  
  !CORRECIÓN DE VELOCIDADES
    !COMPONENTE U
    DO i=1, nx-1
      DO j=1,ny
	u1(i,j)=u1(i,j)+d_h(i,j)*(Pp(i,j)-Pp(i+1,j))
      END DO
    END DO
    
    !COMPONENTE V
    DO i=1,nx
      DO j=1,ny-1
	v1(i,j)=v1(i,j)+d_v(i,j)*(Pp(i,j)-Pp(i,j+1))
      END DO
    END DO

    !Para la divergencia y el criterio de paro
    Div=MAXVAL(ABS(Sp(1:ei,1:ej))) !Divergencia de u=0, ecuación de continuidad
    c=c+1 !Aumenta el contador de iteraciones SIMPLEC
    
  END DO !TERMINA CICLO SIMPLEC &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  
  u=u1; v=v1
  WRITE(*,*)it,c,Div,n
  
  CALL interpolateToNodesUs(uc,u,nx,ny)
  CALL interpolateToNodesVs(vc,v,nx,ny)
  
  CALL move_particles(nx,ny,n,np,x0,y0,dx,dy,dt,ptcon,pt,x,y,xc,yc,u,v)
  
  !Escritura de campo de velocidades 
  IF (MOD(it,200) .EQ. 0) THEN
  
    CALL WriteVectorField('vel',it,uc,vc,xc,yc,nx,ny)
  
  !para el archivo de animación
    WRITE(itchar,'(i6)')it
    itchar=ADJUSTL(itchar)
    itchar=itchar(1:LEN_TRIM(itchar))
!     WRITE(1,*)"p 'vel"//itchar(1:LEN_TRIM(itchar))//".txt' u 1:2:(1*$3):(1*$4) w vec"
!     WRITE(1,*)'pause 0.05'

    OPEN(50,FILE='pt'//itchar(1:LEN_TRIM(itchar))//'.txt',STATUS='REPLACE')
    j=1
    DO WHILE (ptcon(j) .NE. 1)
      WRITE(50,*)pt(j,1),pt(j,2)
      j=ptcon(j)
    END DO
    WRITE(50,*)pt(j,1),pt(j,2)
    WRITE(50,*)pt(1,1),pt(1,2)
    CLOSE(50)
  
    WRITE(1,*)"p 'pt"//itchar(1:LEN_TRIM(itchar))//".txt' u 1:2 lc 14, 'iman.txt' w filledcurves lc -1"
  END IF

END DO !TERMINA CICLO TEMPORAL

CLOSE(1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
END PROGRAM

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
SUBROUTINE MESH_1D(nx,x0,xl,x,xc)
IMPLICIT NONE
INTEGER i, nx
REAL*4 x0,xl,dx
REAL*4 x(0:nx),xc(0:nx+1)

dx=1.0/FLOAT(nx)

DO i=0,nx
  x(i)=x0+FLOAT(i)*(xl-x0)*dx
END DO

xc(0)=x(0); xc(nx+1)=x(nx)

DO i=1,nx
  xc(i)=0.5*(x(i)+x(i-1))
END DO

END SUBROUTINE

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

subroutine TDMA(x,a,b,c,d,n)

implicit none
 integer n,k
 real a(n),b(n),c(n),d(n),x(n),m

 do k=2,N
  m=a(k)/b(k-1)
  b(k)=b(k)-m*c(k-1)
  d(k)=d(k)-m*d(k-1)
 end do

 x(n)=d(n)/b(n)

 do k=n-1,1,-1
  x(k)=(d(k)-c(k)*x(k+1))/b(k)
 end do

end subroutine

!***************************************************************

subroutine lineX_2D(phi,ei,ej,aP,aE,aW,aN,aS,sP,nx,ny)

implicit none
integer bi,ei,bj,ej,i,j,nx,ny
real phi(0:ei+1,0:ej+1),aP(nx,ny),aE(nx,ny),aW(nx,ny),aN(nx,ny),aS(nx,ny),sP(nx,ny)
real a(ei),b(ei),c(ei),d(ei)
bi=1; bj=1
do j=bj,ej
do i=bi,ei
	a(i)=-aW(i,j)
	b(i)=aP(i,j)
	c(i)=-aE(i,j)
	d(i)=sp(i,j) + aN(i,j) * phi(i,j+1) + aS(i,j) * phi(i,j-1)
end do
 call TDMA(phi(bi:ei,j), a, b, c ,d ,ei)
end do

end subroutine

!****************************************************************

subroutine lineY_2D(phi,ei,ej,aP,aE,aW,aN,aS,sP,nx,ny)

implicit none
integer bi,ei,bj,ej,i,j,nx,ny
real phi(0:ei+1,0:ej+1),aP(nx,ny),aE(nx,ny),aW(nx,ny),aN(nx,ny),aS(nx,ny),sP(nx,ny)
real a(ej),b(ej),c(ej),d(ej)
bi=1; bj=1

do i=bi,ei
do j=bj,ej
	a(j)=-aS(i,j)
	b(j)=aP(i,j)
	c(j)=-aN(i,j)
	d(j)=sp(i,j) + aE(i,j) * phi(i+1,j) + aW(i,j) * phi(i-1,j) 
end do
 call TDMA(phi(i,bj:ej), a, b, c ,d ,ej)
end do

end subroutine

!****************************
!****************************

real*4 function calcResidual(phi,ei,ej,aP,aE,aW,aN,aS,sP,nx,ny)
implicit none
integer bi,ei,bj,ej,i,j,nx,ny
real*4 phi(0:ei+1,0:ej+1),aP(nx,ny),aE(nx,ny),aW(nx,ny),aN(nx,ny),aS(nx,ny),sP(nx,ny)
real*4 acum(ei,ej),residual,NINV
bi=1; bj=1
acum=0
NINV = 1.0 / dfloat(ei*ej)
do i=bi,ei
do j=bj,ej
acum(i,j) = aE(i,j) * phi(i+1,j) + aW(i,j) * phi(i-1,j) + aN(i,j) * phi(i,j+1) + &
aS(i,j) * phi(i,j-1) + sp(i,j) - aP(i,j) * phi(i,j)
end do
end do
residual = sqrt( NINV * sum(acum * acum) )
calcResidual=residual
end function

!************************************************************************************************************

subroutine Gauss_TDMA2D(phi,ei,ej,aP,aE,aW,aN,aS,sP,nx,ny,max_iter,tolerance,residual)

implicit none

integer bi,ei,bj,ej,i,j,nx,ny,count_iter,max_iter
real*4 phi(0:ei+1,0:ej+1),aP(nx,ny),aE(nx,ny),aW(nx,ny),aN(nx,ny),aS(nx,ny),sP(nx,ny)
real*4 residual,tolerance
	
	interface
		real*4 function calcResidual(phi,ei,ej,aP,aE,aW,aN,aS,sP,nx,ny)
		implicit none
		integer bi,ei,bj,ej,i,j,nx,ny
		real*4 phi(0:ei+1,0:ej+1),aP(nx,ny),aE(nx,ny),aW(nx,ny),aN(nx,ny),aS(nx,ny),sP(nx,ny)
		end function
	end interface

count_iter=0;  residual=1.0

do while((count_iter <= max_iter).and.(residual > tolerance)) 
call lineX_2D(phi,ei,ej,aP,aE,aW,aN,aS,sP,nx,ny)
call lineY_2D(phi,ei,ej,aP,aE,aW,aN,aS,sP,nx,ny)
residual = calcResidual(phi,ei,ej,aP,aE,aW,aN,aS,sP,nx,ny)
count_iter=count_iter+1
end do

end subroutine

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
subroutine interpolateToNodesUs(uc,us,nx,ny)

implicit none
integer nx,ny
real:: us(0:nx,0:ny+1),uc(0:nx+1,0:ny+1)
integer bi,ei,bj,ej,i,j
bi=1; bj=1;ej=ny; ei=nx

! Internal points
do i=bi,ei
do j=bj-1,ej+1
	uc(I,J) = ( us(I-1,J) + us(I,J) ) * 0.5
end do 
end do

uc(bi-1,:)=us(0,:)
uc(ei+1,:)=us(nx,:)

end subroutine

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

subroutine interpolateToNodesVs(vc,vs,nx,ny)
implicit none
integer nx,ny
real:: vs(0:nx+1,0:ny),vc(0:nx+1,0:ny+1)
integer bi,ei,bj,ej,i,j
bi=1; bj=1;ej=ny; ei=nx

do i=bi-1,ei+1
do j=bj,ej
	vc(I,J) = ( vs(I,J) + vs(I,J-1) ) * 0.5
end do 
end do

vc(:,bj-1)=vs(:,0)
vc(:,ej+1)=vs(:,ny)
end subroutine

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

Subroutine WriteVectorField(Name,kx,uc,vc,xc,yc,nx,ny)
integer i,j,nx,ny,kx
real*4 uc(0:nx+1,0:ny+1),vc(0:nx+1,0:ny+1),xc(0:nx+1),yc(0:ny+1)
character*(*)Name
character*50 txt,Filename
write(txt,'(i6)')kx
txt=ADJUSTL(txt)
Filename=name//txt(1:len_trim(txt))//".txt"

open(11,file=Filename(1:len_trim(Filename)))
	do i=0,nx+1
	do j=0,ny+1
	write(11,*)xc(i),yc(j),uc(i,j),vc(i,j)
	end do
	write(11,*)''
	end do
close(11)
End Subroutine

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::.
Subroutine WriteScalarField2D(Name,kx,T,xc,yc,nx,ny)
integer i,j,nx,ny,kx
real*4 T(0:nx+1,0:ny+1),xc(0:nx+1),yc(0:ny+1)
character*(*)Name
character*50 txt,Filename
write(txt,'(i6)')kx
txt=ADJUSTL(txt)
Filename=name//txt(1:len_trim(txt))//".txt"

open(10,file=Filename(1:len_trim(Filename)))
	do j=0,ny+1
	do i=0,nx+1
	write(10,*)xc(i),yc(j),T(i,j)
	end do
	write(10,*)''
	end do
close(10)
End Subroutine

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
  
  !Horizontal velocity
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
  
  !Vertical velocity
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


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! NON UNIFORM MAGNETIC FIELD
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!MAGNETIC FIELD FROM FURLANI
SUBROUTINE B0z_Field(nx,ny,dc,xc,yc,Bz)
IMPLICIT NONE
INTEGER i, j, nx, ny,in,im,ik
REAL*4 x,y,z,a,b,c,d,dc,E1,xn(2),ym(2),zk(2),g
REAL*4 B0(0:nx+1,0:ny+1),Bz(0:nx+1,0:ny+1),xc(0:nx+1),yc(0:ny+1)
Bz=0.0; B0=0.0
a=0.5*dc; b=0.5*dc; c=0.5*dc !magnet dimensions in [m]
xn(1)=-a; xn(2)=a; ym(1)=-b; ym(2)=b; zk(1)=-c; zk(2)=c
d=3.e-3		!Distance between the bottom of the cavity and the surface of the magnet

E1=2.9e-2		!Normalization parameter (magnetization of the magnet)

DO in=1,2
DO im=1,2
DO ik=1,2
  DO i=0, nx+1
    DO j=0, ny+1
      x=xc(i); y=yc(j)
	z=c+d
	g=1.0/SQRT((x-xn(in))**2+(y-ym(im))**2+(z-zk(ik))**2)
	B0(i,j)=((-1.0)**(in+im+ik))*ATAN(((x-xn(in))*(y-ym(im)))/(z-zk(ik))*g)
    END DO
  END DO
  Bz=Bz+B0
END DO
END DO
END DO

Bz=E1*Bz

END SUBROUTINE B0z_Field
