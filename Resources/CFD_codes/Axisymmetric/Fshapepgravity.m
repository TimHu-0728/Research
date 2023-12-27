function Fshape=Fshapepgravity(x,nz,zz,dz,dz2,R,VF,theta,grav,rho)%  Ecuación de Laplace-Young + conservación del volumen
 %global k0
%rho is the density
%grav is the gravity
 %x vector N+1, con N primeros la forma y x(N+1), C
 f=x(1:nz);
 h=x(nz+1:2*nz);

 pref=x(2*nz+1);
 
 z=h;
 
 fz=f*dz';
 fzz=f*dz2';
 hz=h*dz';
 hzz=h*dz2';

 
 
  
  %Interface curvature
  
n=hz.^2+fz.^2;
k=-((hzz.*(fz)-hz.*fzz)./n.^1.5+hz./f./n.^0.5);
%k0=k;
 
Fshape0=pref-k+rho*grav*z;
  
 % Laplace-Young
 Fshape0(1)=f(1); 
 Fshape0(nz)=f(nz)-R;

 
 
 Fshape1=fz.*fzz+hz.*hzz;

%  dnf=-2*alpha-2*alpha^2*(1-z);
%  Fshape1=fz.*fzz+hz.*hzz-dnf;
%  
 Fshape1(1)=hz(1);
 Fshape1(nz)=hz(nz)+tan(pi/2-theta)*fz(nz);
 
 
  %VF
 VF0=0;
  for j=2:nz
     ap=h(j)*f(j)*fz(j);
     am=h(j-1)*f(j-1)*fz(j-1);
     
 VF0=VF0+0.5*(zz(j)-zz(j-1))*(ap+am);    
 end

 Fshape=[Fshape0,Fshape1,VF-VF0];
  
 
% %volumen
%  V0=0.;
%  for j=2:nz
% V0=V0+0.5*(h(j)-h(j-1))*(f(j)^2+f(j-1)^2); 
%  end
% V0=V0/V1;
%  
%  
%  Fshape(2*nz+1)=V-V0;


