%Ejemplo pagina 130 libro G. R. Liu - Mesh Free Methods_ Moving Beyond the Finite Element Method-CRC Press (2003)

clear all
xa(1)=0.5;
xa(2)=0;
xa(3)=0;
xa(4)=1;
xa(5)=1;
xa(6)=2;
xa(7)=2;
ya(1)=0.5;
ya(2)=1;
ya(3)=0;
ya(4)=1;
ya(5)=0;
ya(6)=1+0.1;
ya(7)=0;




   %number of points 
    M=length(xa);
    N=M-1;
    u=ya(2:M).*xa(2:M);
    
    % coordenates
    x0=xa(1) ; y0=ya(1); 
    
    %normalizamos distancias
    dmax=max(max(xa(2:M)-x0),max(ya(2:M)-y0));
    
  %local coordinates
  x1r=(xa(2:M)-x0)/dmax
  y1r=(ya(2:M)-y0)/dmax
  x2r=x1r.^2 ; y2r=y1r.^2;
  x1ry1r=y1r.*x1r;
 %Normalize Moment matrix
 PQ=zeros(N,N);
 PQ(:,1)=1;
 PQ(:,2)=x1r';
 PQ(:,3)=y1r';
 PQ(:,4)=x1ry1r';
 PQ(:,5)=x2r';
 PQ(:,6)=y2r';
 %inv(PQ)
 [L, U] = lu(PQ);

A=inv(PQ)
  %para mi punto central
 pt=[1,x1r(1),y1r(1),x1r(1)*y1r(1),x1r(1).^2,y1r(1)^2];

  phi=pt*A;
  up=phi*u'
  
   
 phix=[0,1,0,x1r(1),2*x1r(1),0]*A/dmax;
 phiy=[0,0,1,y1r(1),0,2*y1r(1)]*A/dmax;
 phiyx=[0,0,0,1,0,0]*A/dmax^2;
 phixx=[0,0,0,0,2,0]*A/dmax^2;
 phiyy=[0,0,0,0,0,2]*A/dmax^2;
 upx=phixx*u'
 upy=phiyy*u'
 upyx=phiyx*u'
  