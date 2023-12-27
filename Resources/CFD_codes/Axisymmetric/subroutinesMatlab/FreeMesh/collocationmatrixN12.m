function [dx,dy,dxx,dyy,C]=collocationmatrixN12(la,f1,g1)

%la: los marcadores a las 6 columnas distintas de cero la(n,6);
% xa posición radial de los 6 puntos xa(n,6);  x0=xa(:,1)  
% ya posición axil de los 6 puntos xa(n,6);    y0=ya(:,) punto central
n=length(la(:,1));

    
 
  duxva=zeros(n,12);
  duyva=zeros(n,12);
  duyyva=zeros(n,12);
  duxxva=zeros(n,12);
  C=zeros(n,1);
  tic
 for l=1:n

 rbb=f1(la(l,:));
 zbb=g1(la(l,:));

 [duxva(l,:),duyva(l,:),duxxva(l,:),duyyva(l,:), C(l)]=collocationmatrixN12a(rbb,zbb);
 end

 %motamos matrices
 duxva=reshape(duxva,n*12,1);
 duyva=reshape(duyva,n*12,1);
 duxxva=reshape(duxxva,n*12,1);
 duyyva=reshape(duyyva,n*12,1);
 ja1=reshape(la,n*12,1);
% ja1=[la(:,1);la(:,2);la(:,3);la(:,4);la(:,5);la(:,6)];
 


ia=[1:n,1:n,1:n,1:n,1:n,1:n,1:n,1:n,1:n,1:n,1:n,1:n];
 dx=sparse(ia',ja1,duxva);
 dy=sparse(ia',ja1,duyva);
 dxx=sparse(ia',ja1,duxxva);
 dyy=sparse(ia',ja1,duyyva);
nx=length(dx(:,1));
ny=length(dx(1,:)); 
if(abs(nx-ny)>0)
dx(n,n)=0;
dy(n,n)=0;
dxx(n,n)=0;
dyy(n,n)=0;
end
% 
%  for l=1:n
%  dx(l,la(l,1:12))=duxva(l,1:12);
%  dy(l,la(l,1:12))=duyva(l,1:12);
%  dxx(l,la(l,1:12))=duxxva(l,1:12);
%  dyy(l,la(l,1:12))=duyyva(l,1:12);
%  end

 
 

    
  


    
    
        


