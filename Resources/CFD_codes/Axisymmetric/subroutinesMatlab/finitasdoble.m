function[zz,dz,dz2]=finitasdoble(nz,H)
[zz1,dz,dz2]=finitas4ordentotal(nz,1); %z  liquid
a=0.1;

zz=H*(zz1-a*sin(2*pi*zz1));
[zz,dz,dz2]=finitasvariable(nz,zz);
dz=sparse(dz);
dz2=sparse(dz2);







