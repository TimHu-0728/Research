function [dx,dy,dxx,dyy]=collocationmatrix(la,xa,ya)

%la: los marcadores a las 6 columnas distintas de cero la(n,6);
% xa posición radial de los 6 puntos xa(n,6);  x0=xa(:,1)  
% ya posición axil de los 6 puntos xa(n,6);    y0=ya(:,) punto central
n=length(la(:,1));

    
    % coordenates
    x0=xa(:,1) ; y0=ya(:,1);  
x1r=zeros(n,5);
y1r=zeros(n,5);
    for k=1:5
    x1r(:,k)=xa(:,k+1)-x0;
    y1r(:,k)=ya(:,k+1)-y0;
    end

     %x1r=xa(:,2:6)-x0;
    %y1r=ya(:,2:6)-y0;
   
  

%     x1r=x1r' ; y1r=y1r'; va=va';
    x2r=x1r.*x1r;
    y2r=y1r.*y1r;
    x1ry1r=x1r.*y1r;
    
 
% calculamos las 6 componetes no ceros
 [duxva,duyva,duxxva,duyyva]=collocation(x1r,y1r,x2r,y2r,x1ry1r);

 %motamos matrices
 duxva=reshape(duxva,n*6,1);
 duyva=reshape(duyva,n*6,1);
 duxxva=reshape(duxxva,n*6,1);
 duyyva=reshape(duyyva,n*6,1);
 ja1=reshape(la,n*6,1);
% ja1=[la(:,1);la(:,2);la(:,3);la(:,4);la(:,5);la(:,6)];
 ia=[1:n,1:n,1:n,1:n,1:n,1:n];
 dx=sparse(ia',ja1,duxva);
 dy=sparse(ia',ja1,duyva);
 dxx=sparse(ia',ja1,duxxva);
 dyy=sparse(ia',ja1,duyyva);
%  dx=speye(n,n);
%  dy=speye(n,n);
%  dxx=speye(n,n);
%  dyy=speye(n,n);
% %  for k=1:6
% 
% 
%  for l=1:n
%  dx(l,la(l,1:6))=duxva(l,1:6);
%  dy(l,la(l,1:6))=duyva(l,1:6);
%  dxx(l,la(l,1:6))=duxxva(l,1:6);
%  dyy(l,la(l,1:6))=duyyva(l,1:6);
%  end

 
 

    
  


    
    
        


