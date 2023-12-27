

 function Fshape1=FdomainD(x,nz,zz,dz,dz2,nr,rr,dr,dr2,fl,gl,fcoil,gcoil,fb,ft,H2, Rout,zout,CoilD,nx)%  Ecuación de Laplace-Young + conservación del volumen
 %FdomainD(x,nzD,z0D,dz0D,dz0D,nrD,r0D,dr0D,drr0D,fl,gl, fcoil,gcoil,H2,Rout,zout),x0,options);
nt=nz*nr;
 %x vector N+1, con N primeros la forma y x(N+1), C
 f=x(1:nt);

% h=x(nt+1:2*nt);
 f=reshape(f,nr,nz);

 
 fzz=f*dz2';

 frr=dr2*f;


  
 fz=f*dz';

 fr=dr*f;


 
Fshapef=fzz+frr;




Fshapef(:,nz)=f(:,nz)-Rout;


Fshapef(1,:)=fr(1,:);

%Fshapef(1,nx:2*nx)=f(1,nx:2*nx)-fb(1,nx:2*nx);



Fshapef(nr,:)=fr(nr,:);
%Fshapef(nr,nx:2*nx)=f(nr,nx:2*nx)-ft(nr,nx:2*nx);


Fshapef(:,1)=f(:,1)-fl(:,1);

% %corner
 Fshapef(1,nz)=f(1,nz)-Rout;
% 
% 
 Fshapef(nr,nz)=f(nr,nz)-Rout;




 Fshapef=reshape(Fshapef,1,nz*nr);
 Fshapef(CoilD)=f(CoilD)-fcoil(CoilD);
% Fshapeh=reshape(Fshapeh,1,nz*nr);
% Fshapeh(CoilD)=h(CoilD)-gcoil(CoilD);



Fshape1=[Fshapef];
 

 

