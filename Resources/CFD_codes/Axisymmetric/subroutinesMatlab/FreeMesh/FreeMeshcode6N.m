%Ejemplo pagina 130 libro G. R. Liu - Mesh Free Methods_ Moving Beyond the Finite Element Method-CRC Press (2003)

clear all
clear all
xa(1)=0.5;
xa(2)=0;
xa(3)=0;
xa(4)=1;
xa(5)=1;
xa(6)=2;
%xa(7)=2;
ya(1)=0.5;
ya(2)=1;
ya(3)=0;
ya(4)=1;
ya(5)=0;
ya(6)=1;
%ya(7)=0;




   %number of points 
    M=length(xa);
    N=M;


plot(xa,ya,'o')
u=xa;
   %number of points 
    
    % coordenates
%     x0=xa(1) ; y0=ya(1); 
%     
%     %normalizamos distancias
%     dmax=max(max(xa(2:M)-x0),max(ya(2:M)-y0));
%     
  %local coordinates
  x1r=xa(1:M)
  y1r=ya(1:M);
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
 A=inv(PQ)
 

 


 %para mi punto central (1,1)
 pt=[1,xa(1),ya(1),xa(1)*ya(1),xa(1).^2,ya(1)^2];
 ptx=[0,1,0,xa(1),2*xa(1),ya(1)];
 pty=[0,1,0,ya(1),0,2*ya(1)];
 ptyx=[0,0,0,1,0,0];
 ptxx=[0,0,0,0,2,0];
 ptyy=[0,0,0,0,0,2];

 up=pt*A*u'
 upx=ptx*A*u'
 upxx=ptxx*A*u'
 upy=pty*A*u'
 upyy=ptyy*A*u'
 upyx=ptyx*A*u'
 
 
 
 
hold on
 plot(xa(1),ya(1),'x')
 
 
 
 stop

    