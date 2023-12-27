function [dx,dy,dxx,dyy]=collocationmatrixtaylorN5(la,f1,g1)

%la: los marcadores a las 5 columnas distintas de cero la(n,5);
% xa posición radial de los 5 puntos xa(n,5);  x0=xa(:,1)  
% ya posición axil de los 5 puntos xa(n,5);    y0=ya(:,) punto central
n=length(la(:,1));

    
 
  duxva=zeros(n,5);
  duyva=zeros(n,5);
  duyyva=zeros(n,5);
  duxxva=zeros(n,5);
  tic
 for l=1:n

 rbb=f1(la(l,:));
 zbb=g1(la(l,:));

 [duxva(l,:),duyva(l,:),duxxva(l,:),duyyva(l,:)]=collocationmatrixtaylorN5a(rbb,zbb);
 end

 %motamos matrices
 duxva=reshape(duxva,n*5,1);
 duyva=reshape(duyva,n*5,1);
 duxxva=reshape(duxxva,n*5,1);
 duyyva=reshape(duyyva,n*5,1);
 ja1=reshape(la,n*5,1);
% ja1=[la(:,1);la(:,2);la(:,3);la(:,4);la(:,5);la(:,5)];
 ia=[1:n,1:n,1:n,1:n,1:n];
 dx=sparse(ia',ja1,duxva);
 dy=sparse(ia',ja1,duyva);
 dxx=sparse(ia',ja1,duxxva);
 dyy=sparse(ia',ja1,duyyva);
%  dx=speye(n,n);
%  dy=speye(n,n);
%  dxx=speye(n,n);
%  dyy=speye(n,n);
% %  for k=1:5
% 
% 
%  for l=1:n
%  dx(l,la(l,1:5))=duxva(l,1:5);
%  dy(l,la(l,1:5))=duyva(l,1:5);
%  dxx(l,la(l,1:5))=duxxva(l,1:5);
%  dyy(l,la(l,1:5))=duyyva(l,1:5);
%  end

 
 

    
  


    
    
        


