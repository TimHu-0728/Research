function [dx,dy,dxx,dyy]=collocationmatrixN8(la,f1,g1)

%la: los marcadores a las 6 columnas distintas de cero la(n,6);
% xa posición radial de los 6 puntos xa(n,6);  x0=xa(:,1)  
% ya posición axil de los 6 puntos xa(n,6);    y0=ya(:,) punto central
n=length(la(:,1));

    
 
  duxva=zeros(n,8);
  duyva=zeros(n,8);
  duyyva=zeros(n,8);
  duxxva=zeros(n,8);
  tic
 for l=1:n

 rbb=f1(la(l,:));
 zbb=g1(la(l,:));

 [duxva(l,:),duyva(l,:),duxxva(l,:),duyyva(l,:)]=collocationmatrixN8a(rbb,zbb);
 end

 %motamos matrices
 duxva=reshape(duxva,n*8,1);
 duyva=reshape(duyva,n*8,1);
 duxxva=reshape(duxxva,n*8,1);
 duyyva=reshape(duyyva,n*8,1);
 ja1=reshape(la,n*8,1);
% ja1=[la(:,1);la(:,2);la(:,3);la(:,4);la(:,5);la(:,6)];
 ia=[1:n,1:n,1:n,1:n,1:n,1:n,1:n,1:n];
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

 
 

    
  


    
    
        


