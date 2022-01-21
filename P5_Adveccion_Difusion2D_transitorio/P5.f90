PROGRAM Poisson2D !Advección-Difusión 2D temporal
!Declaración de variables
CHARACTER*50 itchar,name
INTEGER i,j,nx,ny,max_iter,it,itmax,ix
REAL*4 dx,dy,x0,xl,y0,yl,dv,Se,Sw,Sn,Ss,tolerance,residual,pi,dt,time,ue,uw,vn,vs,gamma,Pe,Tpe,Tpw
REAL*4, ALLOCATABLE :: aE(:,:),aW(:,:),aN(:,:),aS(:,:),aP(:,:),Sp(:,:),phi(:,:),x(:),xc(:),y(:),yc(:)
REAL*4, ALLOCATABLE :: uc(:,:),vc(:,:)

!Definición de pi
pi=ACOS(-1.0)

max_iter=10000
tolerance=1e-6

nx=50; ny=50
x0=0.0; xl=1.0
y0=0.0; yl=1.0
time=10.0 !Time=15
dt=0.001 !dt=0.005

Tpe=0.0; tpw=0.0

Pe=100.0 !Pe=200
gamma=1.0/Pe

itmax=INT(time/dt)+1

dx=(xl-x0)/FLOAT(nx)
dy=(yl-y0)/FLOAT(ny)

Se=dy; Sw=dy
Sn=dx; Ss=dx

dv=dx*dy

!Alojamiento de memoria dinámica
ALLOCATE(phi(0:nx+1,0:ny+1),aE(nx,ny),aW(nx,ny),aN(nx,ny),aS(nx,ny),aP(nx,ny),Sp(nx,ny))
ALLOCATE(x(0:nx),xc(0:nx+1),y(0:ny),yc(0:ny+1),uc(0:nx+1,0:ny+1),vc(0:nx+1,0:ny+1))

!creación de la malla
CALL MESH_1D(nx,x0,xl,x,xc)
CALL MESH_1D(ny,y0,yl,y,yc)

!Campo de velocidades
DO i=0,nx+1
  DO j=0,ny+1
    uc(i,j)=-SIN(pi*xc(i))*COS(pi*yc(j))
    vc(i,j)=COS(pi*xc(i))*SIN(pi*yc(j))
  END DO
END DO

phi=0.0

!Condiciones de frontera
!NORTE
phi(:,ny+1)=1.0

!SUR
phi(:,0)=0.0

!Creación del archivo de animación
OPEN(2,FILE='anim.gnp',STATUS='REPLACE')
WRITE(2,*)'set xrange[0:1.0]; set yrange[0:1.0]; set view map; set size square'
WRITE(2,*)"set xlabel 'x'; set ylabel 'y' rotate by 0"


!Inicia ciclo temporal
DO it=1,itmax

  !Cálculo de coeficientes
  DO i=1,nx
    DO j=1,ny
      !Cálculo de los flujos por las caras
      ue=0.5*(uc(i,j)+uc(i+1,j))
      uw=0.5*(uc(i,j)+uc(i-1,j))
      vn=0.5*(vc(i,j)+vc(i,j+1))
      vs=0.5*(vc(i,j)+vc(i,j-1))
      
      !cálculo de coeficientes
      aE(i,j)=gamma*Se/dx-MIN(ue*Se,0.0)!0.5*ue*Se
      aW(i,j)=gamma*Sw/dx+MAX(uw*Sw,0.0)!0.5*uw*Sw
      aN(i,j)=gamma*Sn/dy-MIN(vn*Sn,0.0)!0.5*vn*Sn
      aS(i,j)=gamma*Ss/dy+MAX(vs*Ss,0.0)!0.5*vs*Ss
      aP(i,j)=aE(i,j)+aW(i,j)+aN(i,j)+aS(i,j)+dv/dt
      Sp(i,j)=phi(i,j)*dv/dt
    END DO
  END DO

  !Correción por condiciones de frontera
  !ESTE (Neumann)
  aP(nx,:)=aP(nx,:)-aE(nx,:)
  Sp(nx,:)=Sp(nx,:)+aE(nx,:)*dx*Tpe !Tpe derivada de T  en frontera este
  aE(nx,:)=0.0
  
  !OESTE (Neumann)
  aP(1,:)=aP(1,:)-aW(1,:)
  Sp(1,:)=Sp(1,:)-aW(1,:)*dx*Tpw !Tpw derivada de T en la frontera oeste
  aW(1,:)=0.0
  
  !NORTE (Dirichlet)
  aP(:,ny)=aP(:,ny)+aN(:,ny)
  Sp(:,ny)=Sp(:,ny)+2.0*aN(:,ny)*phi(1:nx,ny+1)
  aN(:,ny)=0.0
  
  !SUR (Dirichlet)
  aP(:,1)=aP(:,1)+aS(:,1)
  Sp(:,1)=Sp(:,1)+2.0*aS(:,1)*phi(1:nx,0)
  aS(:,1)=0.0

  !Solución del sistema de ecuaciones
    CALL Gauss_TDMA2D(phi,nx,ny,aP,aE,aW,aN,aS,sP,nx,ny,max_iter,tolerance,residual)
  
  !Cálculo de phi en fronteras este y oeste
  !ESTE
  phi(nx+1,:)=phi(nx,:)+0.5*dx*Tpe
  
  !OESTE
  phi(0,:)=phi(1,:)-0.5*dx*Tpw
  
  !Escribir ene el archivo de animación
  IF (MOD(it,100) .EQ. 0) THEN
    CALL WriteScalarField2D('temp',it,phi,xc,yc,nx,ny)
    WRITE(itchar,'(i6)')it
    itchar=ADJUSTL(itchar)
    name='temp'//itchar(1:LEN_TRIM(itchar))//'.txt'
    
    WRITE(2,*)'sp "'//name(1:LEN_TRIM(name))//'" w pm3d'
    WRITE(2,*)'pause 0.5'
  END IF
  
END DO !Termina ciclo temporal

CLOSE(2)

CALL WriteScalarField2D('temperature',1,phi,xc,yc,nx,ny)

!velocidades
! DO i=0, nx+1
!   DO j=0, ny+1
!     uc(i,j)= 0.015*uc(i,j)/(uc(i,j)**2+vc(i,j)**2)**0.5
!     vc(i,j)= 0.015*vc(i,j)/(uc(i,j)**2+vc(i,j)**2)**0.5
!   END DO
! END DO
  
OPEN(3,FILE='vel.txt',STATUS='replace') !Para abrir un archivo en FORTRAN
DO i=0, nx+1
  DO j=0, ny+1
    WRITE(3,*)xc(i), yc(j), uc(i,j), vc(i,j)
  END DO
END DO

CLOSE(3)
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

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
SUBROUTINE Jacobi(A, T_vec, S, n, max_iter, tolerance)
IMPLICIT NONE
REAL*4 :: T_vec(n), S(n), A(n,n), T_veck_1(n)
INTEGER it, k, i, j, max_iter, n
REAL*4 suma, error, tolerance
error=1.0
k=1
T_veck_1=T_vec
  DO WHILE ((k .LE. max_iter) .AND. (error .GT. tolerance))
    DO i=1, n
      suma=0.0
      DO j=1, n
        suma=suma+A(i,j)*T_veck_1(j) !no considera que j debe ser diferente de i
      END DO
      suma=suma-A(i,i)*T_veck_1(i) !dado que i diferente de j
      T_vec(i)=(S(i)-suma)/A(i,i)
    END DO
  error=MAXVAL(ABS(T_vec-T_veck_1))
  k=k+1
  T_veck_1=T_vec
END DO
! WRITE(*,*)k, error
END SUBROUTINE