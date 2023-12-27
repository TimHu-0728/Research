function[z,d,d2]=finitas4ordentotalperiodic(nz,H)
z=zeros(1,nz);
dz=H/(nz-1);
for j=1:nz
z(j)=(j-1)*dz;
end
%use c,x,dx,dx2,dp,u,v,du
dx=0.5/dz;
dx1=1/(12*dz);
dx2=1/(12*dz^2);
% Matrix of derivatives

d=zeros(nz,nz);
d2=zeros(nz,nz);




for j=3:nz-2; %forth order

d(j,j-2)=dx1;
d(j,j-1)=-8*dx1;
d(j,j)=0;
d(j,j+1)=8*dx1;
d(j,j+2)=-dx1;


d2(j,j-2)=-dx2;
d2(j,j-1)=16*dx2;
d2(j,j)=-30*dx2;
d2(j,j+1)=16*dx2;
d2(j,j+2)=-dx2;

end
j=1
d(j,nz-2)=dx1;
d(j,nz-1)=-8*dx1;
d(j,j)=0;
d(j,j+1)=8*dx1;
d(j,j+2)=-dx1;
d2(j,nz-2)=-dx2;
d2(j,nz-1)=16*dx2;
d2(j,j)=-30*dx2;
d2(j,j+1)=16*dx2;
d2(j,j+2)=-dx2;

j=2
d(j,nz-1)=dx1;
d(j,1)=-8*dx1;
d(j,j)=0;
d(j,j+1)=8*dx1;
d(j,j+2)=-dx1;

d2(j,nz-1)=-dx2;
d2(j,1)=16*dx2;
d2(j,j)=-30*dx2;
d2(j,j+1)=16*dx2;
d2(j,j+2)=-dx2;



j=nz
d(j,j-2)=dx1;
d(j,j-1)=-8*dx1;
d(j,j)=0;
d(j,2)=8*dx1;
d(j,3)=-dx1;


d2(j,j-2)=-dx2;
d2(j,j-1)=16*dx2;
d2(j,j)=-30*dx2;
d2(j,2)=16*dx2;
d2(j,3)=-dx2;


j=nz-1
d(j,j-2)=dx1;
d(j,j-1)=-8*dx1;
d(j,j)=0;
d(j,j+1)=8*dx1;
d(j,2)=-dx1;


d2(j,j-2)=-dx2;
d2(j,j-1)=16*dx2;
d2(j,j)=-30*dx2;
d2(j,j+1)=16*dx2;
d2(j,2)=-dx2;


d2=sparse(d2);
d=sparse(d);

