function[z,dz,dz2]=Chevitanh4th(nr1,Rin, Rout,a)

[x,dx,dx2]=finitas4ordentotal(nr1,1); %[0 1]  
%
% parameter a 
%
% Transforming function z = Rin +(Rout-Rin)* Tanh(a*x)/Tanh(a); 0 < x < 1;
% Accumulates points at Rout 
%a = 3.;
z=Rin + (Rout-Rin)*tanh(a*x)/tanh(a);
dzdx=2.*a*(Rout-Rin)/tanh(a)./(1+cosh(2*a*x));
dz=diag(1./dzdx)*dx;
dz2=dz*dz;