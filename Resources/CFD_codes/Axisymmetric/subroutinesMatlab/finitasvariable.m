function[z,d,d2]=finitasvariable(nz,z)
%Fnite differences 4 order
%second matrix
ip(1,1)=1;
ip(1,2)=2;
ip(1,3)=3;
ip(1,4)=4;
ip(1,5)=5;
ip(2,1)=2;
ip(2,2)=1;
ip(2,3)=3;
ip(2,4)=4;
ip(2,5)=5;
for j=3:nz-2
ip(j,1)=j;
ip(j,2)=j-2;
ip(j,3)=j-1;
ip(j,4)=j+1;
ip(j,5)=j+2;
end
ip(nz-1,1)=nz-1;
ip(nz-1,2)=nz;
ip(nz-1,3)=nz-2;
ip(nz-1,4)=nz-3;
ip(nz-1,5)=nz-4;
ip(nz,1)=nz;
ip(nz,2)=nz-1;
ip(nz,3)=nz-2;
ip(nz,4)=nz-3;
ip(nz,5)=nz-4;
f0a=zeros(nz,5);
f1a=zeros(nz,5);
f2a=zeros(nz,5);
f3a=zeros(nz,5);
f4a=zeros(nz,5);

dz2a=speye(nz,nz);
dza=speye(nz,nz);
for j=1:nz

      
f0a(j,1)=1;
f1a(j,1)=0;
f2a(j,1)=0;
f3a(j,1)=0;
f4a(j,1)=0;

f0a(j,2)=0;
f1a(j,2)=1;
f2a(j,2)=0;
f3a(j,2)=0;
f4a(j,2)=0;

f0a(j,3)=0;
f1a(j,3)=0;
f2a(j,3)=1;
f3a(j,3)=0;
f4a(j,3)=0;

f0a(j,4)=0;
f1a(j,4)=0;
f2a(j,4)=0;
f3a(j,4)=1;
f4a(j,4)=0;

f0a(j,5)=0;
f1a(j,5)=0;
f2a(j,5)=0;
f3a(j,5)=0;
f4a(j,5)=1;
end



for j=1:nz
dz=z(ip(j,2))-z(ip(j,1));
dz1=z(ip(j,3))-z(ip(j,1));
dz2=z(ip(j,4))-z(ip(j,1));
dz3=z(ip(j,5))-z(ip(j,1));


dz1=dz1/dz;
dz2=dz2/dz;
dz3=dz3/dz;

for k=1:5    
f0=f0a(j,k);
f1=f1a(j,k);
f2=f2a(j,k);
f3=f3a(j,k);
f4=f4a(j,k);
a1=((2*dz1^3*dz2*f0 - 2*dz1*dz2^3*f0 + 2*dz1*dz3^3*f0 - 2*dz1^3*dz3*f0 - 2*dz2*dz3^3*f0 + 2*dz2^3*dz3*f0 + 2*dz1*dz2^3*f4 - 2*dz1*dz3^3*f3 + 2*dz2*dz3^3*f2 - 2*dz1^3*dz2*f4 + 2*dz1^3*dz3*f3 - 2*dz2^3*dz3*f2)*dz^4 + (2*dz1*dz2^4*f0 - 2*dz1^4*dz2*f0 - 2*dz1*dz3^4*f0 + 2*dz1^4*dz3*f0 + 2*dz2*dz3^4*f0 - 2*dz2^4*dz3*f0 - 2*dz1*dz2^4*f4 + 2*dz1*dz3^4*f3 - 2*dz2*dz3^4*f2 + 2*dz1^4*dz2*f4 - 2*dz1^4*dz3*f3 + 2*dz2^4*dz3*f2)*dz^3 + (2*dz1^4*dz2^3*f0 - 2*dz1^3*dz2^4*f0 + 2*dz1^3*dz3^4*f0 - 2*dz1^4*dz3^3*f0 - 2*dz2^3*dz3^4*f0 + 2*dz2^4*dz3^3*f0 + 2*dz1^3*dz2^4*f4 - 2*dz1^3*dz3^4*f3 - 2*dz1^4*dz2^3*f4 + 2*dz1^4*dz3^3*f3 + 2*dz2^3*dz3^4*f2 - 2*dz2^4*dz3^3*f2)*dz + 2*dz1*dz2^3*dz3^4*f0 - 2*dz1*dz2^4*dz3^3*f0 - 2*dz1^3*dz2*dz3^4*f0 + 2*dz1^3*dz2^4*dz3*f0 + 2*dz1^4*dz2*dz3^3*f0 - 2*dz1^4*dz2^3*dz3*f0 - 2*dz1*dz2^3*dz3^4*f1 + 2*dz1*dz2^4*dz3^3*f1 + 2*dz1^3*dz2*dz3^4*f1 - 2*dz1^3*dz2^4*dz3*f1 - 2*dz1^4*dz2*dz3^3*f1 + 2*dz1^4*dz2^3*dz3*f1)/((dz1^3*dz2^2*dz3 - dz1^3*dz2*dz3^2 - dz1^2*dz2^3*dz3 + dz1^2*dz2*dz3^3 + dz1*dz2^3*dz3^2 - dz1*dz2^2*dz3^3)*dz^4 + (- dz1^4*dz2^2*dz3 + dz1^4*dz2*dz3^2 + dz1^2*dz2^4*dz3 - dz1^2*dz2*dz3^4 - dz1*dz2^4*dz3^2 + dz1*dz2^2*dz3^4)*dz^3 + (dz1^4*dz2^3*dz3 - dz1^4*dz2*dz3^3 - dz1^3*dz2^4*dz3 + dz1^3*dz2*dz3^4 + dz1*dz2^4*dz3^3 - dz1*dz2^3*dz3^4)*dz^2 + (- dz1^4*dz2^3*dz3^2 + dz1^4*dz2^2*dz3^3 + dz1^3*dz2^4*dz3^2 - dz1^3*dz2^2*dz3^4 - dz1^2*dz2^4*dz3^3 + dz1^2*dz2^3*dz3^4)*dz);
a2=((dz1^2*dz2^3*f0 - dz1^3*dz2^2*f0 - dz1^2*dz3^3*f0 + dz1^3*dz3^2*f0 + dz2^2*dz3^3*f0 - dz2^3*dz3^2*f0 - dz1^2*dz2^3*f4 + dz1^2*dz3^3*f3 + dz1^3*dz2^2*f4 - dz1^3*dz3^2*f3 - dz2^2*dz3^3*f2 + dz2^3*dz3^2*f2)*dz^4 + (dz1^4*dz2^2*f0 - dz1^2*dz2^4*f0 + dz1^2*dz3^4*f0 - dz1^4*dz3^2*f0 - dz2^2*dz3^4*f0 + dz2^4*dz3^2*f0 + dz1^2*dz2^4*f4 - dz1^2*dz3^4*f3 + dz2^2*dz3^4*f2 - dz1^4*dz2^2*f4 + dz1^4*dz3^2*f3 - dz2^4*dz3^2*f2)*dz^3 + (dz1^3*dz2^4*f0 - dz1^4*dz2^3*f0 - dz1^3*dz3^4*f0 + dz1^4*dz3^3*f0 + dz2^3*dz3^4*f0 - dz2^4*dz3^3*f0 - dz1^3*dz2^4*f4 + dz1^3*dz3^4*f3 + dz1^4*dz2^3*f4 - dz1^4*dz3^3*f3 - dz2^3*dz3^4*f2 + dz2^4*dz3^3*f2)*dz^2 + dz1^2*dz2^4*dz3^3*f0 - dz1^2*dz2^3*dz3^4*f0 + dz1^3*dz2^2*dz3^4*f0 - dz1^3*dz2^4*dz3^2*f0 - dz1^4*dz2^2*dz3^3*f0 + dz1^4*dz2^3*dz3^2*f0 + dz1^2*dz2^3*dz3^4*f1 - dz1^2*dz2^4*dz3^3*f1 - dz1^3*dz2^2*dz3^4*f1 + dz1^3*dz2^4*dz3^2*f1 + dz1^4*dz2^2*dz3^3*f1 - dz1^4*dz2^3*dz3^2*f1)/((dz1^3*dz2^2*dz3 - dz1^3*dz2*dz3^2 - dz1^2*dz2^3*dz3 + dz1^2*dz2*dz3^3 + dz1*dz2^3*dz3^2 - dz1*dz2^2*dz3^3)*dz^4 + (- dz1^4*dz2^2*dz3 + dz1^4*dz2*dz3^2 + dz1^2*dz2^4*dz3 - dz1^2*dz2*dz3^4 - dz1*dz2^4*dz3^2 + dz1*dz2^2*dz3^4)*dz^3 + (dz1^4*dz2^3*dz3 - dz1^4*dz2*dz3^3 - dz1^3*dz2^4*dz3 + dz1^3*dz2*dz3^4 + dz1*dz2^4*dz3^3 - dz1*dz2^3*dz3^4)*dz^2 + (- dz1^4*dz2^3*dz3^2 + dz1^4*dz2^2*dz3^3 + dz1^3*dz2^4*dz3^2 - dz1^3*dz2^2*dz3^4 - dz1^2*dz2^4*dz3^3 + dz1^2*dz2^3*dz3^4)*dz);

a1=-(dz1^2*dz2^3*f0 - dz1^3*dz2^2*f0 - dz1^2*dz2^4*f0 - dz1^2*dz3^3*f0 + dz1^3*dz3^2*f0 + dz1^4*dz2^2*f0 + dz1^2*dz3^4*f0 + dz1^3*dz2^4*f0 + dz2^2*dz3^3*f0 - dz1^4*dz2^3*f0 - dz1^4*dz3^2*f0 - dz2^3*dz3^2*f0 - dz1^3*dz3^4*f0 - dz2^2*dz3^4*f0 + dz1^4*dz3^3*f0 + dz2^4*dz3^2*f0 - dz1^2*dz2^3*f4 + dz1^2*dz3^3*f3 + dz1^3*dz2^2*f4 - dz1^3*dz3^2*f3 - dz2^2*dz3^3*f2 + dz2^3*dz3^2*f2 + dz2^3*dz3^4*f0 - dz2^4*dz3^3*f0 + dz1^2*dz2^4*f4 - dz1^2*dz3^4*f3 + dz2^2*dz3^4*f2 - dz1^4*dz2^2*f4 + dz1^4*dz3^2*f3 - dz2^4*dz3^2*f2 - dz1^3*dz2^4*f4 + dz1^3*dz3^4*f3 + dz1^4*dz2^3*f4 - dz1^4*dz3^3*f3 - dz2^3*dz3^4*f2 + dz2^4*dz3^3*f2 - dz1^2*dz2^3*dz3^4*f0 + dz1^2*dz2^4*dz3^3*f0 + dz1^3*dz2^2*dz3^4*f0 - dz1^3*dz2^4*dz3^2*f0 - dz1^4*dz2^2*dz3^3*f0 + dz1^4*dz2^3*dz3^2*f0 + dz1^2*dz2^3*dz3^4*f1 - dz1^2*dz2^4*dz3^3*f1 - dz1^3*dz2^2*dz3^4*f1 + dz1^3*dz2^4*dz3^2*f1 + dz1^4*dz2^2*dz3^3*f1 - dz1^4*dz2^3*dz3^2*f1)/(dz1*dz2*dz3*(dz1 - dz2)*(dz1 - dz3)*(dz2 - dz3)*(dz1 - 1)*(dz2 - 1)*(dz3 - 1));
a2=(2*(dz1^3*dz2^4*f0 - dz1^4*dz2^3*f0 - dz1^3*dz3^4*f0 + dz1^4*dz3^3*f0 + dz2^3*dz3^4*f0 - dz2^4*dz3^3*f0 - dz1^3*dz2^4*f4 + dz1^3*dz3^4*f3 + dz1^4*dz2^3*f4 - dz1^4*dz3^3*f3 - dz2^3*dz3^4*f2 + dz2^4*dz3^3*f2 + dz1*dz2^3*f0 - dz1^3*dz2*f0 - dz1*dz2^4*f0 - dz1*dz3^3*f0 + dz1^3*dz3*f0 + dz1^4*dz2*f0 + dz1*dz3^4*f0 + dz2*dz3^3*f0 - dz1^4*dz3*f0 - dz2^3*dz3*f0 - dz2*dz3^4*f0 + dz2^4*dz3*f0 - dz1*dz2^3*f4 + dz1*dz3^3*f3 - dz2*dz3^3*f2 + dz1^3*dz2*f4 - dz1^3*dz3*f3 + dz2^3*dz3*f2 + dz1*dz2^4*f4 - dz1*dz3^4*f3 + dz2*dz3^4*f2 - dz1^4*dz2*f4 + dz1^4*dz3*f3 - dz2^4*dz3*f2 - dz1*dz2^3*dz3^4*f0 + dz1*dz2^4*dz3^3*f0 + dz1^3*dz2*dz3^4*f0 - dz1^3*dz2^4*dz3*f0 - dz1^4*dz2*dz3^3*f0 + dz1^4*dz2^3*dz3*f0 + dz1*dz2^3*dz3^4*f1 - dz1*dz2^4*dz3^3*f1 - dz1^3*dz2*dz3^4*f1 + dz1^3*dz2^4*dz3*f1 + dz1^4*dz2*dz3^3*f1 - dz1^4*dz2^3*dz3*f1))/(dz1*dz2*dz3*(dz1 - dz2)*(dz1 - dz3)*(dz2 - dz3)*(dz1 - 1)*(dz2 - 1)*(dz3 - 1));


 dz2a(j,ip(j,k))=a2/dz^2;
 dza(j,ip(j,k))=a1/dz;
end
end
d=dza;
d2=dz2a;
d=sparse(d);
d2=sparse(d2);
%d2=d*d;


