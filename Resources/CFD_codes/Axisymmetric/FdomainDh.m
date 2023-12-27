

 function Fshape1=FdomainDh(x,nz,zz,dz,dz2,nr,rr,dr,dr2,fl,gl,fcoil,gcoil,H2, Rout,zout,CoilD)%  Ecuación de Laplace-Young + conservación del volumen
 %FdomainD(x,nzD,z0D,dz0D,dz0D,nrD,r0D,dr0D,drr0D,fl,gl, fcoil,gcoil,H2,Rout,zout),x0,options);
nt=nz*nr;
 %x vector N+1, con N primeros la forma y x(N+1), C

 h=x;
 h=reshape(h,nr,nz);
%  h(:,1)=gl(:,1);
%  h(:,nz)=gl(:,1);
 hzz=h*dz2';
 hrr=dr2*h;

  

 hz=h*dz';
 hr=dr*h;

 

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



Fshape1=[Fshapeh];
 

 

