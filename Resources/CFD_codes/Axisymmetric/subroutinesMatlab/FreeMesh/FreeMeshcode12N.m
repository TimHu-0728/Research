%Ejemplo pagina 130 libro G. R. Liu - Mesh Free Methods_ Moving Beyond the Finite Element Method-CRC Press (2003)

clear all
clear all
eps=0
xa(1)=0+0.4*eps;
xa(2)=0+eps;
xa(3)=0-eps;
xa(4)=0+2*eps;
xa(5)=2;
xa(6)=2+eps;
xa(7)=2;
xa(8)=2;
xa(9)=3-eps;
xa(10)=3;
xa(11)=3+eps;
xa(12)=3;


ya(1)=3;
ya(2)=2-eps;
ya(3)=1;
ya(4)=0+eps;
ya(5)=3-eps;
ya(6)=2;
ya(7)=1+eps;
ya(8)=0;
ya(9)=3+eps;

ya(10)=2+eps;
ya(11)=1-eps;
ya(12)=1;
%ya(7)=0;




   %number of points 
    M=length(xa);
    N=M;


plot(xa,ya,'o')
hold on
 plot(xa(1),ya(1),'x')
u=1+xa.^2;
   %number of points 
%    idcoli
M=12
    % coordenates
     x0=xa(1) ; y0=ya(1); 
%     
%     %normalizamos distancias
     dmax=max(max(xa(1:M)-x0),max(ya(1:M)-y0));
%     
  %local coordinates
  x1r=(xa-x0)/dmax;
  y1r=(ya-y0)/dmax;
  x2r=x1r.^2 ; y2r=y1r.^2;
  x3r=x1r.^3 ; y3r=y1r.^3;
  x1ry1r=y1r.*x1r;
 %Normalize Moment matrix
 PQ=zeros(N,N);
 PQ(:,1)=1;
 PQ(:,2)=x1r';
 PQ(:,3)=y1r';
 PQ(:,4)=x1ry1r';
 PQ(:,5)=x2r';
 PQ(:,6)=y2r';
 PQ(:,7)=x2r'.*y1r';
 PQ(:,8)=x1r'.*y2r';
 PQ(:,9)=x2r'.*y2r';
 PQ(:,10)=x3r';
 PQ(:,11)=y3r';
 PQ(:,12)=x3r'.*y1r';
%  [R,idx] = rref(PQ);
%   [R,idy] = rref(PQ')
 
 %a=PQ(:,jb)
  [Xsub,idy]=licols(PQ,1e-6)
  [Xsub1,idx]=licols(PQ',1e-6);
  PQ1=PQ(idx,idy);
  [L, U] = lu(PQ);
%L
A=inv(PQ1)


 %para mi punto central (1,1)
  x=x1r(1);
  y=y1r(1);
 pt=[1,x,y,x*y,x^2,y^2,x^2*y,y^2*x, x^2*y^2,x^3,y^3,x^3*y];
 ptx=[ 0, 1, 0, y, 2*x, 0, 2*x*y, y^2, 2*x*y^2, 3*x^2, 0, 3*x^2*y];
 pty=[ 0, 0, 1, x, 0, 2*y, x^2, 2*x*y, 2*x^2*y, 0, 3*y^2, x^3];
 ptxx=[ 0, 0, 0, 0, 2, 0, 2*y, 0, 2*y^2, 6*x, 0, 6*x*y];
 ptyy=[ 0, 0, 0, 0, 0, 2, 0, 2*x, 2*x^2, 0, 6*y, 0];
 %listv={'pt','ptx','pty',
 
 
  pt=pt(idy);
 ptx=ptx(idy);
 pty=pty(idy);
 ptxx=ptxx(idy);
 ptyy=ptyy(idy);

 
 
 up=pt*A*u(idx)'
 upx=ptx*A*u(idx)'/dmax
 upxx=ptxx*A*u(idx)'/dmax^2
 upy=pty*A*u(idx)'/dmax
 upyy=ptyy*A*u(idx)'/dmax^2
% upyx=ptyx*A*u(idx)'/dmax^2
 
 
 
 

 
 
 
 stop

    