function[z,dz,dz2]=Chevitanh_inverse(nr1,Rin, Rout,alphaf)

[x,dx,dx2]=Chevi(nr1,1,0.); %[0 1]  
%
% parameter a 
%
% Transforming function z = Rin +(Rout-Rin)* Tanh(a*(1-x)/Tanh(a); 0 < x < 1;
% Accumulates points at Rout 
a = alphaf;
z=Rout + (Rin-Rout)*tanh(a*(1-x))/tanh(a);
dzdx=-2.*a*(Rin-Rout)/tanh(a)./(1+cosh(2*a*(1-x)));
dz=diag(1./dzdx)*dx;
dz2=dz*dz;