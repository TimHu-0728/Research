function[zz,dz,dz2]=finitasdoblegood(nz,H)

%[zz1,dza,dz2a]=finitas2(nz,1)
[zz1,dza,dz2a]=finitas4ordentotal(nz,1); %z  liquid
a=0.1;

zz=H*(zz1-a*sin(2*pi*zz1));
df1=-H*(2*a*pi*cos(2*pi*zz1) - 1);
df2=4*H*a*pi^2*sin(2*pi*zz1);
% df1=zz*dza';
% df2=zz*dz2a';



dz=dza;
dz2=dz2a;
for j=1:nz
    dz(j,:)=dza(j,:)/df1(j);
    dz2(j,:)=dz2a(j,:)/df1(j)^2-dza(j,:)*df2(j)/df1(j)^3;
end
%[zz,dz,dz2]=finitasvariable(nz,zz);
dz=sparse(dz);
dz2=sparse(dz2);







