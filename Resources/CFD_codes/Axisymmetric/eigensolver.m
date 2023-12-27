function[V,d]=eigensolver(A,B,omega,n)

 B=-B*1i;
 A=-A*1i;

options.disp=0;



 [V,d] = eigs(A,B,n,omega,options);
% A=full(A);
% B=full (B);
% [V,d] = eig(A,B);

