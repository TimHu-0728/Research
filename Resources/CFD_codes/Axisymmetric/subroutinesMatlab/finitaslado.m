function[zz,dz,dz2]=finitaslado(nz,H,a)
[zz1,dz,dz2]=finitas4ordentotal(nz,1); %z  liquid
%a=0.3;
%a=0.4;
zz=H*(zz1-a*sin(pi*zz1));
[zz,dz,dz2]=finitasvariable(nz,zz);

dz=sparse(dz);
dz2=sparse(dz2);







