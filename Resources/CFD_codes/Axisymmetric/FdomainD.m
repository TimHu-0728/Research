

 function Fshape1=FdomainD(x,nz,zz,dz,dz2,nr,rr,dr,dr2,fl,gl,fcoil,gcoil,H2, Rout,zout,CoilD)%  Ecuación de Laplace-Young + conservación del volumen
 %FdomainD(x,nzD,z0D,dz0D,dz0D,nrD,r0D,dr0D,drr0D,fl,gl, fcoil,gcoil,H2,Rout,zout),x0,options);
nt=nz*nr;
 %x vector N+1, con N primeros la forma y x(N+1), C
 f=x(1:nt);
 h=x(nt+1:2*nt);
% h=x(nt+1:2*nt);
 f=reshape(f,nr,nz);
 h=reshape(h,nr,nz);
 
 fzz=f*dz2';
 hzz=h*dz2';
 frr=dr2*f;
 hrr=dr2*h;

  
 fz=f*dz';
 hz=h*dz';
 fr=dr*f;
 hr=dr*h;

 
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




 

Fshapeh=hrr+hzz;
 Fshapeh(:,1)=h(:,1)-gl(:,1);
 %Fshapeh(:,1)=h(:,nz)-gl(:,1);
 Fshapeh(:,nz)=hz(:,nz);


Fshapeh(1,:)=h(1,:)-zout;
Fshapeh(nr,:)=h(nr,:)-H2;


%corner

% Fshapeh(1,nz)=h(1,nz)-zout;
% 
% Fshapeh(nr,nz)=h(nr,nz)-H2;



% Fshapef=reshape(Fshapef,1,nz*nr);
% Fshapef(CoilD)=f(CoilD)-fcoil(CoilD);
 Fshapeh=reshape(Fshapeh,1,nz*nr);
 Fshapeh(CoilD)=h(CoilD)-gcoil(CoilD);

Fshape1=[Fshapeh,Fshapef];
 

 

