function[zz,dz,dz2]=Chevidoble(nz1,H)

[zz,dz,dz2]=Chevi(nz1-1,2,-1); %[-1 1]  

zz1=sin(pi*zz/2.01);
dza=zz1*dz';
dz1=diag(1./dza)*dz;
dz2=dz1*dz1;
%final
zz=H*(zz1-zz1(1))/(zz1(nz1)-zz1(1));
dz=dz1*(zz1(nz1)-zz1(1))/H;
dz2=dz*dz;