function [dx,dy,dxx,dyy,C]=collocationmatrixFreeMesh7(la,f1,g1)

%la: los marcadores a las 6 columnas distintas de cero la(n,6);
% xa posición radial de los 6 puntos xa(n,6);  x0=xa(:,1)  
% ya posición axil de los 6 puntos xa(n,6);    y0=ya(:,) punto central
n=length(la(:,1));
M=7;
  dx=zeros(n,n);
  dy=zeros(n,n);
  dxx=zeros(n,n);
  dyy=zeros(n,n);
  C=zeros(n,1);
  tic
 for l=1:n

 xa=f1(la(l,:));
 ya=g1(la(l,:));
% coordenates
    x0=xa(1) ; y0=ya(1); 
    
    %normalizamos distancias
    dmax=max(max(xa(1:M)-x0),max(ya(1:M)-y0));
 
 
  %local coordinates
  %local coordinates
  x1r=(xa(1:M)-x0)/dmax;
  y1r=(ya(1:M)-y0)/dmax;
  x2r=x1r.^2 ; y2r=y1r.^2;
  x1ry1r=y1r.*x1r;
 %Normalize Moment matrix
 N=M;
 PQ=zeros(N,N);
 PQ(:,1)=1;
 PQ(:,2)=x1r';
 PQ(:,3)=y1r';
 PQ(:,4)=x1ry1r';
 PQ(:,5)=x2r';
 PQ(:,6)=y2r';
 PQ(:,7)=x2r'.*y1r';
 %a=PQ(:,jb)
 % idx=1:7;
% idy=1:7;
 PQ1=PQ;
   [Xsub,idy]=licols(PQ,1e-12);
   [Xsub,idx]=licols(PQ',1e-12);
   PQ1=PQ(idx,idy);
l

 A=inv(PQ1);
  %para mi punto central (1,1)
  x=x1r(1);
  y=y1r(1);
 pt=[1,x,y,x*y,x^2,y^2,x^2*y];
 ptx=[ 0, 1, 0, y, 2*x, 0, 2*x*y];
 pty=[ 0, 0, 1, x, 0, 2*y, x^2];
 ptxx=[ 0, 0, 0, 0, 2, 0, 2*y];
 ptyy=[ 0, 0, 0, 0, 0, 2,0];
 %listv={'pt','ptx','pty',
 %pt=pt(idy);
 ptx=ptx(idy);
 pty=pty(idy);
 ptxx=ptxx(idy);
 ptyy=ptyy(idy); 
 dx(l,la(l,idx))=ptx*A/dmax;
 dy(l,la(l,idx))=pty*A/dmax;
 dxx(l,la(l,idx))=ptxx*A/dmax^2;
 dyy(l,la(l,idx))=ptyy*A/dmax^2;
 end

dx=sparse(dx);
dy=sparse(dy);
dxx=sparse(dxx);
dyy=sparse(dyy);


 
 

    
  


    
    
        


