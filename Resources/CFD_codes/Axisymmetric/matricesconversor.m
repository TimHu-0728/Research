%% Matrix Algebra for Block A
% 


%%
% _r0A_ and _z0A_ must be extended to the entire domain
%

% zA=repmat(z0A, [nrA 1]);
% rA=repmat(r0A', [1 nzA]);

%%
% Testing purposes

% wA=zA.^2;
% uA=rA.^2;
% pA = rA.*zA;
% fA = 1.+0.2*sin(zA);
% phiA = Volt*zA.^2;
% sA = (1.+0.33*cos(zA/2));

%%
% This matrix operation allows to use an array of unknowns 1 x
% (nrA*NzA*NVAR_A) (It is much faster)
%
function [ddr0A,ddrr0A,ddz0A,ddzz0A,deltaA,ddrz0A]= matricesconversor(nrA,nzA,dr0A,drr0A,dz0A,dzz0A)
%converting matrice 
ddr0A=kron(speye(nzA,nzA),dr0A);
ddz0A=kron(dz0A,speye(nrA,nrA));
deltaA=kron(speye(nrA,nrA),speye(nzA,nzA));
ddrr0A=kron(speye(nzA,nzA),drr0A);
ddzz0A=kron(dzz0A,speye(nrA,nrA));
ddrz0A=kron(dz0A,dr0A); 

