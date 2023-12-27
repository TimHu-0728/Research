function[z,dz,dz2]=Chevitanh4th_i(nr1,Rin, Rout,a)

[x,dx,dx2]=finitas4ordentotal(nr1,1); %[0 1]  
%
% parameter a 
%
% Transforming function z = Rin +(Rout-Rin)* Tanh(a*x)/Tanh(a); 0 < x < 1;
% Accumulates points at Rout 
%a = 3.;
z    = Rin + (Rout-Rin)*(1-tanh(a*(1-x))/tanh(a));
dzdx = (a*(Rin - Rout)*(tanh(a*(x - 1)).^2 - 1))/tanh(a);
dz   = diag(1./dzdx)*dx;
dz2  = dz*dz;