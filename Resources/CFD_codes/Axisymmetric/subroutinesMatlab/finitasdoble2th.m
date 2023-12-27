function[zz,dz,dz2]=finitasdoble2th(nz,H)
%[zz,dz,dz2]=Chevi(nz-1,H,0.0);  
[zz1,dz,dz2]=finitas4ordentotal(nz,1); %z  liquid
a=0.1;

zz=H*(zz1-a*sin(2*pi*zz1));
[zz,dz,dz2]=finitasvariable2th(nz,zz);
dz=sparse(dz);
dz2=sparse(dz2);







