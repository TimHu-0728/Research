function[z,d,d2]=finitasvariable2thfullsymetric(nz,z)
%Fnite differences 4 order
%second matrix
ip(1,1)=1;
ip(1,2)=2;
ip(1,3)=3;

for j=2:nz-1
ip(j,1)=j;
ip(j,2)=j-1;
ip(j,3)=j+1;
end

ip(nz,1)=nz;
ip(nz,2)=nz-1;
ip(nz,3)=nz-2;


for j=1:nz

      
f0a(j,1)=1;
f1a(j,1)=0;
f2a(j,1)=0;



f0a(j,2)=0;
f1a(j,2)=1;
f2a(j,2)=0;



f0a(j,3)=0;
f1a(j,3)=0;
f2a(j,3)=1;





end



for j=1:nz
dz=z(ip(j,2))-z(ip(j,1));
dz1=z(ip(j,3))-z(ip(j,1));
for k=1:3    
f0=f0a(j,k);
f1=f1a(j,k);
f2=f2a(j,k);
dz2a(j,ip(j,k))=(2*(dz*f0 - dz1*f0 - dz*f2 + dz1*f1))/(dz*dz1*(dz - dz1));
dza(j,ip(j,k))=-(dz^2*f0 - dz1^2*f0 - dz^2*f2 + dz1^2*f1)/(dz*dz1*(dz - dz1));
 
end
end
dza(1,:)=0;
dz2a(1,:)=0;
dz=(z(2)-z(1));
dx=1/(dz*dz);
dz2a(1,1)=-2*dx;
dz2a(1,2)=2*dx;

dza(nz,:)=0;
dz2a(nz,:)=0;
dz=(z(nz)-z(nz-1));
dx=1/(dz*dz);
dz2a(nz,nz)=-2*dx;
dz2a(nz,nz-1)=2*dx;


d=sparse(dza);
d2=sparse(dz2a);
%d2=d*d;


