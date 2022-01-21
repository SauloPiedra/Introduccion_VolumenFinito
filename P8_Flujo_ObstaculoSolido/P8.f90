PROGRAM Flujo_canal_obs_sol
IMPLICIT NONE

!Declaración de variables
INTEGER nx, ny, i, j, itmax, max_iter, maxit_SIMPLEC, c, ei, ej, it, np, n
REAL*4 dx, dy, dv, dt, time, tolerance, Se, Sw, Sn, Ss, ue, uw, vn, vs, x0, xl, y0, yl, Re, gama, tolerance_s, Div, residual
REAL*4 dbkx, dbky, pvor
CHARACTER*50 Filename, itchar,x0c,y0c,xlc,ylc
REAL*4, ALLOCATABLE :: aE(:,:), aW(:,:), aN(:,:), aS(:,:), aP(:,:), Sp(:,:), x(:), xc(:), y(:), yc(:)
REAL*4, ALLOCATABLE :: u(:,:), v(:,:), u1(:,:), v1(:,:), P(:,:), Pp(:,:), d_h(:,:), d_v(:,:), uc(:,:), vc(:,:)
REAL*4 x0_obs,y0_obs,xl_obs,yl_obs, v_obs,x_obs !variables para definir el obstáculo sólido
INTEGER, ALLOCATABLE :: mark_cells(:,:) !arreglo para marcar celdas
REAL*4, ALLOCATABLE :: pt(:,:) !arreglo para el seguimiento Lagrangiano
!Pp = P prima; u1 y v1 son al tiempo actual

!Inicialización de variables

nx=200; ny=80
np=80000
pvor=7.0
x0=0.0; xl=20.0
y0=-2.5; yl=2.5
time=50.0
dt=0.001
it=0
tolerance=1e-5
tolerance_s=1e-5
max_iter=1000
maxit_SIMPLEC=100
Re=100.0

ALLOCATE(aE(nx,ny), aW(nx,ny), aN(nx,ny), aS(nx,ny), aP(nx,ny), Sp(nx,ny))
ALLOCATE(u(0:nx,0:ny+1), u1(0:nx,0:ny+1), v(0:nx+1,0:ny), v1(0:nx+1,0:ny), P(0:nx+1,0:ny+1), Pp(0:nx+1,0:ny+1))
ALLOCATE(d_h(0:nx+1,0:ny+1), d_v(0:nx+1,0:ny+1), x(0:nx), xc(0:nx+1), y(0:ny), yc(0:ny+1), uc(0:nx+1,0:ny+1), vc(0:nx+1,0:ny+1))
ALLOCATE(mark_cells(0:nx+1,0:ny+1), pt(np,2))

n=0
dx=(xl-x0)/FLOAT(nx)
dy=(yl-y0)/FLOAT(ny)
dv=dx*dy
Se=dy; Sw=dy
Sn=dx; Ss=dx
gama=1.0/Re
itmax=INT(time/dt)+1

CALL MESH_1D(x0, xl, x, xc, nx)
CALL MESH_1D(y0, yl, y, yc, ny)


!Condiciones de frontera
  !no deslizamiento
v=0.0  !no deslizamiento
P=0.0
u=1.0
u(:,1:ny)=1.0 !velocidad de entrada
u(:,ny+1)=0.0; u(:,0)=0.0
u1=u; v1=v
d_h=0.0
d_v=0.0

x0_obs=4.0; xl_obs=5.; y0_obs=-0.5; yl_obs=0.5 !localización del obstáculo
mark_cells=0

OPEN(3,FILE='obstaculo.txt',STATUS='replace')

!cilco para marcar el obstaculo y anular velocidades dentro de él
DO i=1, nx
  DO j=1, ny
    IF ((xc(i) .GT. x0_obs) .AND. (xc(i) .LT. xl_obs) .AND. (yc(j) .GT. y0_obs) .AND. (yc(j) .LT. yl_obs)) THEN !Obstáculo cuadrado
      mark_cells(i,j)=1
      u(i,j)=0.0; v(i,j)=0.0
      WRITE(3,*)xc(i),yc(j)
    END IF
  END DO
END DO
CLOSE(3)

!velocidad en las fronteras del obstáculo sólido
DO i=1, nx
  DO j=1,ny
    IF ((mark_cells(i,j) .EQ. 0) .AND. (mark_cells(i+1,j) .EQ. 1)) THEN
      u(i,j)=0.0
    END IF
    IF ((mark_cells(i,j) .EQ. 0) .AND. (mark_cells(i,j+1) .EQ. 1))  THEN
      v(i,j)=0.0
    END IF
  END DO
END DO
u1=u; v1=v

CALL particle_inyection(pt,yl_obs,y0_obs,n,np,dx,dy)
!Escritura de la posición inicial de las partículas
WRITE(itchar,'(i6)')it
  itchar=ADJUSTL(itchar)
  itchar='lag'//itchar(1:LEN_TRIM(itchar))//'.txt'
OPEN(2,FILE=itchar(1:LEN_TRIM(itchar)),STATUS='replace')

DO i=1,n
  WRITE(2,*)pt(i,1), pt(i,2)
END DO
CLOSE(2)

OPEN(1,FILE='anim.gnp',STATUS='replace')
  !Para el archivo gnp de la animación
    WRITE(1,*)'set xrange[0:20]; set yrange[-2.5:2.5]; unset key; set size ratio 0.25'
    WRITE(1,*)"p '"//itchar(1:LEN_TRIM(itchar))//"' w d ls 2"
!     WRITE(1,*)'pause 0.05'
    
OPEN(3, FILE='anim2.gnp', STATUS='REPLACE')
    WRITE(3,*)'set xrange[0:20]; set yrange[-2.5:2.5]; unset key; set size ratio 0.25'
    
DO it=1, itmax !Inicia ciclo temporal
  Div=1.0
  c=1
  DO WHILE ((DIV .GT. tolerance_s) .AND. (c .LT. maxit_SIMPLEC)) !inicia ciclo SIMPLEC
    aP=0.0; aE=0.0; aW=0.0; aN=0.0; aS=0.0; Sp=0.0 
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  !Ecuación componente u
    ei=UBOUND(u,1)-1
    ej=UBOUND(u,2)-1
    !Cálculo de coeficientes
    DO j=1, ej
      DO i=1, ei
        ue=(u(i,j)+u(i+1,j))*0.5
        uw=(u(i,j)+u(i-1,j))*0.5
        vn=(v(i,j)+v(i+1,j))*0.5
        vs=(v(i,j-1)+v(i+1,j-1))*0.5
        aE(i,j)=gama*Se/dx-MIN(ue*Se,0.0)!ue*Se*0.5
        aW(i,j)=gama*Sw/dx+MAX(uw*Sw,0.0)!uw*Sw*0.5
        aN(i,j)=gama*Sn/dy-MIN(vn*Sn,0.0)!vn*Sn*0.5
        aS(i,j)=gama*Ss/dy+MAX(vs*Ss,0.0)!vs*Ss*0.5
        aP(i,j)=aE(i,j)+aW(i,j)+aN(i,j)+aS(i,j)+dv/dt
        Sp(i,j)=u(i,j)*dv/dt-(P(i+1,j)-P(i,j))*dv/dx
        
        !Correción de frontera inmersa
        !correción de frontera ESTE
        IF ((mark_cells(i,j) .EQ. 0) .AND. (mark_cells(i+2,j) .EQ. 1)) THEN
          Sp(i,j)=Sp(i,j)+aE(i,j)*u(i+1,j)
          aE(i,j)=0.0
        END IF
        !correción frontera OESTE
        IF ((mark_cells(i,j).EQ. 0) .AND. (mark_cells(i-1,j) .EQ. 1)) THEN
          Sp(i,j)=Sp(i,j)+aW(i,j)*u(i-1,j)
          aW(i,j)=0.0
        END IF
        !correción frontera NORTE
        IF ((mark_cells(i,j) .EQ. 0) .AND. ((mark_cells(i+1,j+1) .EQ. 1) .OR. (mark_cells(i,j+1) .EQ. 1))) THEN
          aP(i,j)=aP(i,j)+aN(i,j)
          Sp(i,j)=Sp(i,j)+2*aN(i,j)*u(i,j+1)
          aN(i,j)=0.0
        END IF
        !Corrección frontera SUR
        IF ((mark_cells(i,j) .EQ. 0) .AND. ((mark_cells(i+1,j-1) .EQ. 1) .OR. (mark_cells(i,j-1) .EQ. 1))) THEN
          aP(i,j)=aP(i,j)+aS(i,j)
          Sp(i,j)=Sp(i,j)+2*aS(i,j)*u(i,j-1)
          aS(i,j)=0.0
        END IF
        !dentro del obstáculo
        IF ((mark_cells(i,j) .EQ. 1) .OR. ((mark_cells(i,j) .EQ. 0) .AND. (mark_cells(i+1,j) .EQ. 1))) THEN
          aE(i,j)=0.0; aW(i,j)=0.0; aN(i,j)=0.0; aS(i,j)=0.0; u(i,j)=0.0
          aP(i,j)=dv/dt
          Sp(i,j)=dv/dt*u(i,j)
        END IF  
      END DO
    END DO
    !Corección de condiciones de frontera
    dbkx=0.0; dbky=1.0
    !ESTE
    aP(ei,1:ej)=aP(ei,1:ej)-aE(ei,1:ej)
    !Sp(ei,1:ej)=Sp(ei,1:ej)+(1+dbkx)*aE(ei,1:ej)*u(ei+1,1:ej)  Condicion tipo Neumann
    aE(ei,1:ej)=0.0
    !OESTE
    aP(1,1:ej)=aP(1,1:ej)+dbkx*aW(1,1:ej)
    Sp(1,1:ej)=Sp(1,1:ej)+(1+dbkx)*aW(1,1:ej)*u(0,1:ej)
    aW(1,1:ej)=0.0
    !NORTE
    aP(1:ei,ej)=aP(1:ei,ej)+dbky*aN(1:ei,ej)
    Sp(1:ei,ej)=Sp(1:ei,ej)+(1+dbky)*aN(1:ei,ej)*u(1:ei,ej+1)
    aN(1:ei,ej)=0.0
    !SUR
    aP(1:ei,1)=aP(1:ei,1)+dbky*aS(1:ei,1)
    Sp(1:ei,1)=Sp(1:ei,1)+(1+dbky)*aS(1:ei,1)*u(1:ei,0)
    aS(1:ei,1)=0.0
    !cálculo de coeficientes para la corrección de la presión
    d_h(1:ei,1:ej)=Se/(aP(1:ei,1:ej)-aE(1:ei,1:ej)-aW(1:ei,1:ej)-aN(1:ei,1:ej)-aS(1:ei,1:ej))
    !Solución del sistema de ecuaciones
    CALL Gauss_TDMA2D(u1,ei,ej,aP,aE,aW,aN,aS,Sp,nx,ny,max_iter,tolerance,residual) !Tolerance es del sistema de ecuaciones
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  aP=0.0; aE=0.0; aW=0.0; aN=0.0; aS=0.0; Sp=0.0 
  !Ecuación componente v
    ei=UBOUND(v,1)-1
    ej=UBOUND(v,2)-1
    !Cálculo de coeficientes
    DO j=1, ej
      DO i=1, ei
        ue=(u(i,j)+u(i,j+1))*0.5
        uw=(u(i-1,j)+u(i-1,j+1))*0.5
        vn=(v(i,j)+v(i,j+1))*0.5
        vs=(v(i,j)+v(i,j-1))*0.5
        aE(i,j)=gama*Se/dx-MIN(ue*Se,0.0)!ue*Se*0.5
        aW(i,j)=gama*Sw/dx+MAX(uw*Sw,0.0)!uw*Sw*0.5
        aN(i,j)=gama*Sn/dy-MIN(vn*Sn,0.0)!vn*Sn*0.5
        aS(i,j)=gama*Ss/dy+MAX(vs*Ss,0.0)!vs*Ss*0.5
        aP(i,j)=aE(i,j)+aW(i,j)+aN(i,j)+aS(i,j)+dv/dt
        Sp(i,j)=v(i,j)*dv/dt-(P(i,j+1)-P(i,j))*dv/dy

        !Correción de frontera inmersa
        !correción de frontera ESTE
        IF ((mark_cells(i,j) .EQ. 0) .AND. ((mark_cells(i+1,j+1) .EQ. 1) .OR. (mark_cells(i+1,j) .EQ. 1))) THEN
          aP(i,j)=aP(i,j)+aE(i,j)
          Sp(i,j)=Sp(i,j)+2*aE(i,j)*v(i+1,j)
          aE(i,j)=0.0
        END IF
        !correción frontera OESTE
        IF ((mark_cells(i,j) .EQ. 0) .AND. ((mark_cells(i-1,j+1) .EQ. 1) .OR. (mark_cells(i-1,j) .EQ. 1))) THEN
          aP(i,j)=aP(i,j)+aW(i,j)
          Sp(i,j)=Sp(i,j)+2*aW(i,j)*v(i-1,j)
          aW(i,j)=0.0
        END IF
        !correción frontera NORTE
        IF ((mark_cells(i,j) .EQ. 0) .AND. (mark_cells(i,j+2) .EQ.1)) THEN
          !aP(i,j)=aP(i,j)!+aN(i,j)
          Sp(i,j)=Sp(i,j)+aN(i,j)*v(i,j+1)
          aN(i,j)=0.0
        END IF
        !Corrección frontera SUR
        IF ((mark_cells(i,j) .EQ. 0) .AND. (mark_cells(i,j-1) .EQ. 1)) THEN
          !aP(i,j)=aP(i,j)!+aS(i,j)
          Sp(i,j)=Sp(i,j)+aS(i,j)*v(i,j-1)
          aS(i,j)=0.0
        END IF
        !dentro del obstáculo
        IF ((mark_cells(i,j) .EQ. 1) .OR. ((mark_cells(i,j) .EQ. 0) .AND. (mark_cells(i,j+1) .EQ. 1))) THEN
          aE(i,j)=0.0; aW(i,j)=0.0; aN(i,j)=0.0; aS(i,j)=0.0; v(i,j)=0.0
          aP(i,j)=dv/dt
          Sp(i,j)=dv/dt*v(i,j)
        END IF
      END DO
    END DO
    !Corección de condiciones de frontera
    dbkx=1.0; dbky=0.0
    !ESTE
    aP(ei,1:ej)=aP(ei,1:ej)-dbkx*aE(ei,1:ej)
    !Sp(ei,1:ej)=Sp(ei,1:ej)+(1+dbkx)*aE(ei,1:ej)*v(ei+1,1:ej) Condicion tipo Neumann
    aE(ei,1:ej)=0.0
    !OESTE
    aP(1,1:ej)=aP(1,1:ej)+dbkx*aW(1,1:ej)
    Sp(1,1:ej)=Sp(1,1:ej)+(1+dbkx)*aW(1,1:ej)*v(0,1:ej)
    aW(1,1:ej)=0.0
    !NORTE
    aP(1:ei,ej)=aP(1:ei,ej)+dbky*aN(1:ei,ej)
    Sp(1:ei,ej)=Sp(1:ei,ej)+(1+dbky)*aN(1:ei,ej)*v(1:ei,ej+1)
    aN(1:ei,ej)=0.0
    !SUR
    aP(1:ei,1)=aP(1:ei,1)+dbky*aS(1:ei,1)
    Sp(1:ei,1)=Sp(1:ei,1)+(1+dbky)*aS(1:ei,1)*v(1:ei,0)
    aS(1:ei,1)=0.0
    !cálculo de coeficientes para la corrección de la presión
    d_v(1:ei,1:ej)=Sn/(aP(1:ei,1:ej)-aE(1:ei,1:ej)-aW(1:ei,1:ej)-aN(1:ei,1:ej)-aS(1:ei,1:ej))
    !Solución del sistema de ecuaciones
    CALL Gauss_TDMA2D(v1,ei,ej,aP,aE,aW,aN,aS,Sp,nx,ny,max_iter,tolerance,residual) !Tolerance es del sistema de ecuaciones

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!Correción de la presión.
aP=0.0; aE=0.0; aW=0.0; aN=0.0; aS=0.0; Sp=0.0; Pp=0.0
    ei=UBOUND(Pp,1)-1
    ej=UBOUND(Pp,2)-1
    !Cálculo de flujos en las caras
    DO i=1, ei
      DO j=1,ej
        ue=u1(i,j)
        uw=u1(i-1,j)
        vn=v1(i,j)
        vs=v1(i,j-1)
        !cálculo de coeficientes
        aE(i,j)=d_h(i,j)*Se
        aW(i,j)=d_h(i-1,j)*Sw
        aN(i,j)=d_v(i,j)*Sn
        aS(i,j)=d_v(i,j-1)*Ss
        aP(i,j)=aE(i,j)+aW(i,j)+aN(i,j)+aS(i,j)
        Sp(i,j)=uw*Sw-ue*Se+vs*Ss-vn*Sn   !Divergencia de U*
      END DO
    END DO
    CALL Gauss_TDMA2D(Pp,ei,ej,aP,aE,aW,aN,aS,Sp,nx,ny,max_iter,tolerance,residual) !Cálculo de los coeficientes
    P=P+Pp  !Suma de la correción de la presión y la presión supuesta
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    !Corrección de las componentes de la velocidad
    !Componente u
    DO i=1,nx-1
      DO j=1,ny
        IF ((mark_cells(i,j) .EQ. 0) .AND. (mark_cells(i+1,j) .NE. 1)) THEN
          u1(i,j)=u1(i,j)+(Pp(i,j)-Pp(i+1,j))*d_h(i,j)
        END IF
      END DO
    END DO
    !Componente v
    DO i=1, nx
      DO j=1, ny-1
        IF ((mark_cells(i,j) .EQ. 0) .AND. (mark_cells(i,j+1) .NE. 1)) THEN
          v1(i,j)=v1(i,j)+(Pp(i,j)-Pp(i,j+1))*d_v(i,j)
        END IF  
      END DO
    END DO
      u1(nx,:)=u1(nx-1,:)
      v1(nx+1,:)=v1(nx,:)
    Div=MAXVAL(ABS(Sp(1:ei,1:ej))) !Divergencia de u=0, ecuación de continuidad
    c=c+1
  END DO !Termina ciclo SIMPLEC
  u=u1; v=v1

  WRITE(*,*)it, c, Div

  CALL interpolateToNodesUs(uc,u,nx,ny)
  CALL interpolateToNodesVs(vc,v,nx,ny)

  CALL move_particles(u,v,pt,x,y,dx,dy,dt,x0,y0,nx,ny,n,np)

  IF (MOD(it,200) .EQ. 0) THEN
  !Para escribirlas posiciones de las partículas en un archivo de texto
    WRITE(itchar,'(i6)')it
    itchar=ADJUSTL(itchar)
    itchar='lagrange'//itchar(1:LEN_TRIM(itchar))//'.txt'
    OPEN(2,FILE=itchar(1:LEN_TRIM(itchar)),STATUS='replace')
    DO i=1,n
      WRITE(2,*)pt(i,1), pt(i,2)
    END DO
    CLOSE(2)
    CALL particle_inyection(pt,yl_obs,y0_obs,n,np,dx,dy)
    !WRITE(*,*)n
  !Para el archivo gnp de la animación
    WRITE(1,*)"p '"//itchar(1:LEN_TRIM(itchar))//"' w d ls 4, 'obstaculo.txt' w p ls 20"
    WRITE(1,*)'pause 0.01'
    
    CALL WriteVectorField('vel',it,uc,vc,xc,yc,nx,ny,Filename)
    WRITE(itchar,'(i6)')it
    itchar=ADJUSTL(itchar)
    itchar='vel'//itchar(1:LEN_TRIM(itchar))//'.txt'
    
    WRITE(3,*)"p '"//itchar(1:LEN_TRIM(itchar))//"' u 1:2:(0.3*$3):(0.3*$4) w vec"
  END IF
END DO !Termina ciclo temporal

CLOSE(1)
CLOSE(3)

CALL WriteVectorField('vel',1,uc,vc,xc,yc,nx,ny,Filename)

END PROGRAM

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
SUBROUTINE MESH_1D(x0, xl, x, xc, nx)
IMPLICIT NONE
INTEGER m, nx
REAL*4 x(0:nx), xc(0:nx+1), x0, xl, dx

dx=1.0/FLOAT(nx)

DO m=0, nx
  x(m)=x0+FLOAT(m)*(xl-x0)*dx
END DO

xc(0)=x(0); xc(nx+1)=x(nx)

DO m=1, nx
  xc(m)=(x(m)+x(m-1))*0.5
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

Subroutine WriteVectorField(Name,kx,uc,vc,xc,yc,nx,ny,Filename)
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
	end do
close(11)
End Subroutine

!**************************************************

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
	end do
close(10)
End Subroutine

!******************************************************************************************************************
SUBROUTINE particle_inyection(pt,yl_obs,y0_obs,n,np,dx,dy)
IMPLICIT NONE
REAL*4 pt(np,2)
REAL*4 xp0,ym, dy, dx, yl_obs,y0_obs
INTEGER n,np,ny_obs, i

ym=(yl_obs+y0_obs)/2
ny_obs=INT((1.2*yl_obs-ym)/dy)
xp0=0.1

DO i=1,10*ny_obs
  n=n+1
  pt(n,1)=xp0
  pt(n,2)=ym-dy*FLOAT(i-1)/10
END DO
DO i=1,10*ny_obs
  n=n+1
  pt(n,1)=xp0
  pt(n,2)=ym+dy*FLOAT(i)/10
END DO
END SUBROUTINE

!***************************************************************************************************
SUBROUTINE move_particles(u,v,pt,x,y,dx,dy,dt,x0,y0,nx,ny,n,np)
IMPLICIT NONE
REAL*4 u(0:nx,0:ny+1),v(0:nx+1,0:ny),pt(np,2), x(0:nx),y(0:ny)
REAL*4 dx,dy,dt,up,vp,x0,y0
INTEGER nx,ny,n,np,i,j,ix,iy

DO i=1,n
  ix=INT((pt(i,1)-x0)/dx)+1
  iy=INT((pt(i,2)-y0)/dy)+1
  
  !Interpolación de la velocidad horizontal de la partícula
  up=u(ix-1,iy)+(pt(i,1)-x(ix-1))*(u(ix,iy)-u(ix-1,iy))/dx
  
  !Interpolación de la velocidad vertical de la partícula
  vp=v(ix,iy-1)+(pt(i,2)-y(iy-1))*(v(ix,iy)-v(ix,iy-1))/dy
  
  !Cálculo de la nueva posición de la partícula
  IF ((ix .GE. nx+1) .OR. (iy .LE. 0) .OR. (iy .GE. ny+1)) THEN
    up=0.0
    vp=0.0
  END IF
  
    pt(i,1)=pt(i,1)+up*dt
    pt(i,2)=pt(i,2)+vp*dt
    
END DO
END SUBROUTINE