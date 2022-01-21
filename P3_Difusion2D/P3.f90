PROGRAM Poisson2D !Difusión 2D
!Declaración de variables
INTEGER i,j,nx,ny,max_iter
REAL*4 dx,dy,x0,xl,y0,yl,dv,Se,Sw,Sn,Ss,tolerance,residual
REAL*4, ALLOCATABLE :: aE(:,:),aW(:,:),aN(:,:),aS(:,:),aP(:,:),Sp(:,:),phi(:,:),x(:),xc(:),y(:),yc(:)
REAL*4, ALLOCATABLE :: phia(:,:)

max_iter=1000
tolerance=1e-6

nx=50; ny=50
x0=-1.0; xl=1.0
y0=-1.0; yl=1.0

dx=(xl-x0)/FLOAT(nx)
dy=(yl-y0)/FLOAT(ny)

Se=dy; Sw=dy
Sn=dx; Ss=dx

dv=dx*dy

!Alojamiento de memoria dinámica
ALLOCATE(phi(0:nx+1,0:ny+1),aE(nx,ny),aW(nx,ny),aN(nx,ny),aS(nx,ny),aP(nx,ny),Sp(nx,ny))
ALLOCATE(x(0:nx),xc(0:nx+1),y(0:ny),yc(0:ny+1),phia(0:nx+1,0:ny+1))

!creación de la malla
CALL MESH_1D(nx,x0,xl,x,xc)
CALL MESH_1D(ny,y0,yl,y,yc)

!Condiciones de frontera
!ESTE
phi(nx+1,:)=yc(:)+yc(:)*yc(:)+1.0

!OESTE
phi(0,:)=yc(:)-yc(:)*yc(:)+1.0

!NORTE
phi(:,ny+1)=xc(:)+xc(:)*xc(:)+1.0

!SUR
phi(:,0)=xc(:)-xc(:)*xc(:)+1.0

!Cálculo de coeficientes
DO i=1,nx
  DO j=1,ny
    aE(i,j)=Se/dx
    aW(i,j)=Sw/dx
    aN(i,j)=Sn/dy
    aS(i,j)=Ss/dy
    aP(i,j)=aE(i,j)+aW(i,j)+aN(i,j)+aS(i,j)
    Sp(i,j)=-(2.0*xc(i)+2.0*yc(j))*dv
  END DO
END DO

!Correción por condiciones de frontera
!ESTE
aP(nx,:)=aP(nx,:)+aE(nx,:)
Sp(nx,:)=Sp(nx,:)+2.0*aE(nx,:)*phi(nx+1,1:ny)
aE(nx,:)=0.0

!OESTE
aP(1,:)=aP(1,:)+aW(1,:)
Sp(1,:)=Sp(1,:)+2.0*aW(1,:)*phi(0,1:ny)
aW(1,:)=0.0

!NORTE
aP(:,ny)=aP(:,ny)+aN(:,ny)
Sp(:,ny)=Sp(:,ny)+2.0*aN(:,ny)*phi(1:nx,ny+1)
aN(:,ny)=0.0

!SUR
aP(:,1)=aP(:,1)+aS(:,1)
Sp(:,1)=Sp(:,1)+2.0*aS(:,1)*phi(1:nx,0)
aS(:,1)=0.0

!Solución del sistema de ecuaciones
CALL Gauss_TDMA2D(phi,nx,ny,aP,aE,aW,aN,aS,sP,nx,ny,max_iter,tolerance,residual)

WRITE(*,*)residual

CALL WriteScalarField2D('phi',0,phi,xc,yc,nx,ny)

!Solución analítica
DO i=0, nx+1
  DO j=0, ny+1
    phia(i,j)=xc(i)*xc(i)*yc(j)+xc(i)*yc(j)*yc(j)+1.0
  END DO
END DO

CALL WriteScalarField2D('phia',0,phia,xc,yc,nx,ny)

CALL WriteScalarField2D('phid',0,phia-phi,xc,yc,nx,ny)

END PROGRAM

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE MESH_1D(nx,x0,xl,x,xc)
IMPLICIT NONE
INTEGER nx,i
REAL*4 dx,x0,xl
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