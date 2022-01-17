PROGRAM Calor1D !Ecuación de Poisson para temperatura
IMPLICIT NONE
!Declaración de variables
INTEGER nx, i, j
REAL*4 dx, dv, k, L, Ta, Tb, Se, Sw, q
REAL*4, ALLOCATABLE :: aE(:), aW(:), aP(:), Sp(:), Tn(:),T(:),x(:)

!Inicialización de variables
nx=20; Se=1.0; Sw=1.0
Ta=100.0; Tb=200.0
q=1000E3
L=0.02
k=0.5

dx=L/FLOAT(nx)

dv=dx	!en 1D el volumen de los elementos es la longitud de los mismos

!Alojamiento de memoria dinámica
ALLOCATE(Tn(0:nx+1), aE(nx), aW(nx), aP(nx), Sp(nx), T(0:nx+1),x(0:nx+1))

Tn=0.0

!Condiciones de frontera
Tn(0)=Ta; Tn(nx+1)=Tb

!Algoritmo de solución
  !Cálculo de coeficientes
  
  DO i=1, nx
    aE(i)=k*Se/dx
    aW(i)=k*Sw/dx
    aP(i)=aE(i)+aW(i)
    Sp(i)=q*dv
  END DO
  
  !Correción por condiciones de frontera
  !OESTE
  aP(1)=aP(1)+aW(1)
  Sp(1)=Sp(1)+2.0*aW(1)*Tn(0)
  aW(1)=0.0
  
  !ESTE
  aP(nx)=aP(nx)+aE(nx)
  Sp(nx)=Sp(nx)+2.0*aE(nx)*Tn(nx+1)
  aE(nx)=0.0
  
  !Solución del sistema de ecuaciones por el método de Jacobi
  DO j=1, 1000
    DO i=1, nx
      Tn(i)=(aE(i)*Tn(i+1)+aW(i)*Tn(i-1)+Sp(i))/aP(i)
    END DO
  END DO
  
  OPEN(1,FILE='temperatura.txt',STATUS='REPLACE')
  
  
  !Solución analítica
  T(0)=Ta
  DO i=1, nx
    x(i)=(FLOAT(i)-0.5)*dx
    T(i)=Ta+(Tb-Ta)/L*x(i)-(q/(2.*k))*x(i)*(x(i)-L)
  END DO
  T(nx+1)=Tb
  
!   WRITE(1,*)2.0*0.0,Tn(0),T(0)
  WRITE(1,'(f8.6,2f10.4)')2.0*0.0,Tn(0),T(0)
  DO i=1,nx
!     WRITE(1,*)x(i),Tn(i),T(i)
    WRITE(1,'(f8.6,2f10.4)')x(i),Tn(i),T(i)
  END DO
  WRITE(1,'(f8.6,2f10.4)')L,Tn(nx+1),T(nx+1)
!   WRITE(1,*)L,Tn(nx+1),T(nx+1)
  
  CLOSE(1)
END PROGRAM