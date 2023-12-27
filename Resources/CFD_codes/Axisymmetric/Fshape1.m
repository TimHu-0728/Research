

 function Fshape1=Fshape1(x,nz,zz,dz,dz2,Oh)%  Ecuación de Laplace-Young + conservación del volumen

 %x vector N+1, con N primeros la forma y x(N+1), C
 h=x(1:nz);
 %h=x(nz+1:2*nz);
 
 

 hz=h*dz';
 hzz=h*dz2';

  %Interface curvature
  
n=hz.^2+fz.^2;
k=((hzz.*(fz)-hz.*fzz)./n.^1.5+hz./f./n.^0.5);

 
Fshape0=k+2/R;
  
 % Laplace-Young
 Fshape0(1)=f(1)-R; 
 Fshape0(nz)=f(nz);
 
 Fshape1=fz.*fzz+hz.*hzz;

 
 Fshape1(1)=h(1);
 Fshape1(nz)=hz(nz);
 
 Fshape1=[Fshape0,Fshape1];
 

 
 

