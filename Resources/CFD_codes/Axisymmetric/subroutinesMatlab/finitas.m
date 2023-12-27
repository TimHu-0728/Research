function[z,d,d2]=finitas(nz,H)
z=zeros(1,nz);
dz=H/(nz-1);
for j=1:nz
z(j)=(j-1)*dz;
end
%use c,x,dx,dx2,dp,u,v,du
dx=0.5/dz;
dx1=1/(12*dz);
% Matrix of derivatives

d=zeros(nz,nz);


d(1,1)=-3*dx; %second order
d(1,2)=+4*dx;
d(1,3)=-dx;

d(2,1)=-dx;
d(2,2)=0;
d(2,3)=dx;
for j=3:nz-2; %forth order
d(j,j-2)=dx1;
d(j,j-1)=-8*dx1;
d(j,j)=0;
d(j,j+1)=8*dx1;
d(j,j+2)=-dx1;
end
d(nz-1,nz-2)=-dx;
d(nz-1,nz-1)=0;
d(nz-1,nz)=dx;
d(nz,nz-2)=dx;
d(nz,nz-1)=-4*dx;
d(nz,nz)=3*dx;

dx=1/(dz*dz);
dx1=1/(12*dz*dz);
d2=zeros(nz,nz);
d2(1,1)=0;
d2(1,2)=0;
d2(1,3)=0;

d2(2,1)=dx;
d2(2,2)=-2*dx;
d2(2,3)=dx;
for j=3:nz-2; %forth order
d2(j,j-2)=-dx1;
d2(j,j-1)=16*dx1;
d2(j,j)=0-30*dx1;
d2(j,j+1)=16*dx1;
d2(j,j+2)=-dx1;
end
d2(nz,1)=0;
d2(nz,2)=0;
d2(nz,3)=0;

d2(nz-1,nz-2)=dx;
d2(nz-1,nz-1)=-2*dx;
d2(nz-1,nz)=dx;




