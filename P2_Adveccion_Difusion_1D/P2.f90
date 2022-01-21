PROGRAM adv_diff !div(u*phi)=div(gamma div(phi))
IMPLICIT NONE
!Declaración de variables
INTEGER nx, i, j
REAL*4 dx, dv, L, T0, Tl, Se, Sw, Pe, gamma, U
REAL*4, ALLOCATABLE :: aE(:), aW(:), aP(:), Sp(:), T(:)

!Inicialización de variables
nx=50; Se=1.0; Sw=1.0
T0=1.0; Tl=0.0
L=1.0
U=1.0
Pe=100.0

dx=L/FLOAT(nx)

gamma=1.0/Pe

dv=dx	!en 1D el volumen de los elementos es la longitud de los mismos

!Alojamiento de memoria dinámica
ALLOCATE(T(0:nx+1), aE(nx), aW(nx), aP(nx), Sp(nx))

!Condiciones de frontera
T(0)=T0; T(nx+1)=Tl

T(1:nx)=0.5*(T(0)+T(nx+1))

!Algoritmo de solución
  !Cálculo de coeficientes
  
  DO i=1, nx
    !esquema central
!     aE(i)=gamma*Se/dx-0.5*U*Se
!     aW(i)=gamma*Sw/dx+0.5*U*Sw
    !esquema upwind
    aE(i)=gamma*Se/dx-MIN(U*Se,0.0)!0.5*U*Se
    aW(i)=gamma*Sw/dx+MAX(U*Sw,0.0)!0.5*U*Sw
    aP(i)=aE(i)+aW(i)
    Sp(i)=0.0
  END DO
  
  !Correción por condiciones de frontera
  !OESTE
  aP(1)=aP(1)+aW(1)
  Sp(1)=Sp(1)+2.0*aW(1)*T(0)
  aW(1)=0.0
  
  !ESTE
  aP(nx)=aP(nx)+aE(nx)
  Sp(nx)=Sp(nx)+2.0*aE(nx)*T(nx+1)
  aE(nx)=0.0
  
  !Solución del sistema de ecuaciones por el método de Gauss-Seidel
  DO j=1, 10000
    DO i=1, nx
      T(i)=(aE(i)*T(i+1)+aW(i)*T(i-1)+Sp(i))/aP(i)
    END DO
  END DO
  
  OPEN(1,FILE='temperatura.txt',STATUS='REPLACE')
  
  WRITE(1,*)2.0*0.0,T(0)
  
  DO i=1,nx
    WRITE(1,*)(FLOAT(i)-0.5)*dx,T(i)
  END DO
  
  WRITE(1,*)L,T(nx+1)
  
  CLOSE(1)
END PROGRAM